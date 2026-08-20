#!/usr/bin/env python3
"""
Convert VTK field snapshots to POD coefficients.

Typical usage:

U only:
    python python/scripts/field_pod_from_vtk.py \
        --vtk-glob "openfoam/cases/cylinderLES3D/VTK/cylinderLES3D_*.vtk" \
        --fields U \
        --out data/processed/pod_U_volume_t140_270.npz \
        --time-start 140 \
        --time-step 0.25 \
        --energy-threshold 0.9999 \
        --n-modes-max 200 \
        --pod-weighting volume \
        --cell-stride 1

U and p:
    python python/scripts/field_pod_from_vtk.py \
        --vtk-glob "openfoam/cases/cylinderLES3D/VTK/cylinderLES3D_*.vtk" \
        --fields U p \
        --out data/processed/pod_U_p_t140_270.npz \
        --time-start 140 \
        --time-step 0.25 \
        --energy-threshold 0.9999 \
        --n-modes-max 30 \
        --cell-stride 1
"""

from __future__ import annotations

import argparse
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PYTHON_SRC = PROJECT_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from rom.pod import build_component_weights


def natural_key(path: str) -> list:
    """Sort paths with decimal time suffixes in physical time order."""
    name = Path(path).name
    match = re.search(r"_([0-9]+(?:\.[0-9]+)?)\.vtk$", name)
    if match:
        return [name[: match.start(1)], float(match.group(1))]

    parts = re.split(r"(\d+(?:\.\d+)?)", name)
    key = []
    for part in parts:
        if re.fullmatch(r"\d+(?:\.\d+)?", part):
            key.append(float(part))
        else:
            key.append(part)
    return key


def load_vtk_snapshot(
    path: Path,
    fields: list[str],
    cell_stride: int,
    p_scale: float,
    dtype: np.dtype,
) -> np.ndarray:
    try:
        import pyvista as pv
    except ImportError as exc:
        raise ImportError(
            "pyvista is required. Install with: pip install pyvista vtk"
        ) from exc

    mesh = pv.read(str(path))

    arrays = mesh.cell_data

    available = list(arrays.keys())
    missing = [f for f in fields if f not in available]

    if missing:
        raise KeyError(
            f"Missing fields {missing} in {path}. "
            f"Available cell arrays: {available}"
        )

    parts = []

    for field in fields:
        arr = np.asarray(arrays[field])

        if arr.ndim == 1:
            arr = arr[:, None]

        arr = arr[::cell_stride, :].astype(dtype, copy=False)

        if field == "p":
            arr = arr * dtype.type(p_scale)

        parts.append(arr.reshape(-1))

    return np.concatenate(parts).astype(dtype, copy=False)


def inspect_vtk_layout(
    path: Path,
    fields: list[str],
    cell_stride: int,
) -> tuple[np.ndarray, dict]:
    """Return sampled cell centers and the flattened feature layout."""
    try:
        import pyvista as pv
    except ImportError as exc:
        raise ImportError(
            "pyvista is required. Install with: pip install pyvista vtk"
        ) from exc

    mesh = pv.read(str(path))
    coordinates = np.asarray(mesh.cell_centers().points)[::cell_stride]
    field_entries = []
    feature_start = 0

    for field in fields:
        if field not in mesh.cell_data:
            raise KeyError(
                f"Missing field {field!r} in {path}. "
                f"Available cell arrays: {list(mesh.cell_data.keys())}"
            )

        array = np.asarray(mesh.cell_data[field])
        n_components = 1 if array.ndim == 1 else int(array.shape[1])
        feature_count = len(coordinates) * n_components
        field_entries.append(
            {
                "name": field,
                "n_components": n_components,
                "feature_start": feature_start,
                "feature_stop": feature_start + feature_count,
            }
        )
        feature_start += feature_count

    layout = {
        "n_cells": int(len(coordinates)),
        "coordinate_dimension": int(coordinates.shape[1]),
        "fields": field_entries,
        "n_features": int(feature_start),
        "flattening": "field-major, then cell-major, then component",
    }
    return coordinates, layout


def load_cell_volumes(path: Path, cell_stride: int) -> np.ndarray:
    """Read finite-volume cell volumes from a VTK/OpenFOAM mesh."""
    try:
        import pyvista as pv
    except ImportError as exc:
        raise ImportError(
            "pyvista is required. Install with: pip install pyvista vtk"
        ) from exc

    mesh = pv.read(str(path))
    sized = mesh.compute_cell_sizes(length=False, area=False, volume=True)
    if "Volume" not in sized.cell_data:
        raise ValueError(
            f"Could not obtain cell volumes from {path}. "
            "Volume-weighted POD requires finite-volume cell volumes."
        )
    volumes = np.asarray(sized.cell_data["Volume"], dtype=float)[::cell_stride]
    if volumes.ndim != 1 or volumes.size == 0:
        raise ValueError("Cell volumes are empty.")
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0):
        raise ValueError("Cell volumes must be finite and positive.")
    return volumes


def compute_times(
    n: int,
    time_start: float | None,
    time_step: float | None,
) -> np.ndarray:
    if time_start is not None and time_step is not None:
        return time_start + time_step * np.arange(n, dtype=float)

    return np.arange(n, dtype=float)


def method_of_snapshots_pod(
    matrix_path: Path,
    n_samples: int,
    n_features: int,
    dtype: np.dtype,
    chunk_features: int,
    energy_threshold: float,
    n_modes_max: int | None,
    feature_weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, np.ndarray]:
    """
    Returns:
        coefficients, singular_values, energy, cumulative_energy, n_modes,
        spatial_modes
    """

    X = np.load(matrix_path, mmap_mode="r")

    if feature_weights is None:
        weights = np.ones(n_features, dtype=np.float64)
    else:
        weights = np.asarray(feature_weights, dtype=np.float64)
    if weights.shape != (n_features,):
        raise ValueError(f"feature_weights must have shape ({n_features},), got {weights.shape}")
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ValueError("feature_weights must be finite and positive")

    print("")
    print("==> Computing mean field")

    mean = np.zeros(n_features, dtype=np.float64)

    for j0 in range(0, n_features, chunk_features):
        j1 = min(j0 + chunk_features, n_features)
        block = np.asarray(X[:, j0:j1], dtype=np.float64)
        mean[j0:j1] = block.mean(axis=0)
        print(f"    mean block {j0}:{j1}")

    print("")
    print("==> Computing snapshot covariance matrix")

    C = np.zeros((n_samples, n_samples), dtype=np.float64)

    for j0 in range(0, n_features, chunk_features):
        j1 = min(j0 + chunk_features, n_features)
        block = np.asarray(X[:, j0:j1], dtype=np.float64)
        block -= mean[j0:j1]
        C += (block * weights[j0:j1][None, :]) @ block.T
        print(f"    covariance block {j0}:{j1}")

    if n_samples > 1:
        C /= float(n_samples - 1)

    print("")
    print("==> Eigen-decomposition")

    evals, evecs = np.linalg.eigh(C)

    order = np.argsort(evals)[::-1]
    evals = np.maximum(evals[order], 0.0)
    evecs = evecs[:, order]

    total = evals.sum()

    if total <= 0:
        raise ValueError("Total POD energy is zero.")

    energy = evals / total
    cumulative_energy = np.cumsum(energy)

    n_modes = int(np.searchsorted(cumulative_energy, energy_threshold) + 1)

    if n_modes_max is not None:
        n_modes = min(n_modes, int(n_modes_max))

    singular_values = np.sqrt(evals * max(n_samples - 1, 1))
    singular_tolerance = (
        np.finfo(np.float64).eps
        * max(n_samples, n_features)
        * singular_values[0]
    )
    numerical_rank = int(np.count_nonzero(singular_values > singular_tolerance))
    if numerical_rank == 0:
        raise ValueError("Snapshot matrix has zero numerical POD rank.")

    n_modes = max(1, min(n_modes, numerical_rank))

    # Temporal POD coefficients: a_i(t) = sqrt(lambda_i*(N-1)) * temporal_mode_i(t)
    coefficients = evecs[:, :n_modes] * singular_values[:n_modes]

    print("")
    print("==> Reconstructing spatial POD modes")

    spatial_modes = np.empty((n_features, n_modes), dtype=dtype)
    retained_evecs = evecs[:, :n_modes]
    retained_singular_values = singular_values[:n_modes]

    for j0 in range(0, n_features, chunk_features):
        j1 = min(j0 + chunk_features, n_features)
        block = np.asarray(X[:, j0:j1], dtype=np.float64)
        block -= mean[j0:j1]
        # Physical modes are W^{-1/2} times the right singular vectors of
        # the weighted snapshot matrix X W^{1/2}, which simplifies to
        # X^T U / S for centered physical snapshots.
        mode_block = block.T @ retained_evecs
        mode_block /= retained_singular_values[None, :]
        spatial_modes[j0:j1, :] = mode_block.astype(dtype, copy=False)
        print(f"    spatial mode block {j0}:{j1}")

    return (
        coefficients,
        singular_values,
        energy,
        cumulative_energy,
        n_modes,
        spatial_modes,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--vtk-glob", required=True, help="Glob for internal VTK files.")
    parser.add_argument("--fields", nargs="+", default=["U"], help="Fields to flatten, e.g. U or U p.")
    parser.add_argument("--out", required=True, type=Path, help="Output .npz file.")

    parser.add_argument("--time-start", type=float, default=None)
    parser.add_argument("--time-step", type=float, default=None)

    parser.add_argument("--energy-threshold", type=float, default=0.9999)
    parser.add_argument("--n-modes-max", type=int, default=30)

    parser.add_argument(
        "--cell-stride",
        type=int,
        default=1,
        help="Use every N-th cell. 1 = full field, 2 = half spatial sampling, etc.",
    )

    parser.add_argument(
        "--p-scale",
        type=float,
        default=1.0,
        help="Scale pressure before concatenating with velocity.",
    )

    parser.add_argument(
        "--pod-weighting",
        choices=["volume", "uniform"],
        default="volume",
        help="POD inner-product weighting. Use volume for OpenFOAM finite-volume fields.",
    )

    parser.add_argument(
        "--dtype",
        choices=["float32", "float64"],
        default="float32",
        help="Storage dtype for snapshot matrix.",
    )

    parser.add_argument(
        "--chunk-features",
        type=int,
        default=200000,
        help="Feature chunk size for covariance computation.",
    )

    parser.add_argument(
        "--max-snapshots",
        type=int,
        default=None,
        help="Optional limit for testing.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.cell_stride < 1:
        raise ValueError("--cell-stride must be >= 1")

    dtype = np.dtype(args.dtype)

    files = sorted(glob.glob(args.vtk_glob), key=natural_key)

    if args.max_snapshots is not None:
        files = files[: args.max_snapshots]

    if not files:
        raise FileNotFoundError(f"No VTK files found for glob: {args.vtk_glob}")

    files = [Path(f) for f in files]

    print("==> VTK files")
    print(f"    count: {len(files)}")
    print(f"    first: {files[0]}")
    print(f"    last : {files[-1]}")
    print(f"    fields: {args.fields}")
    print(f"    cell_stride: {args.cell_stride}")
    print(f"    pod_weighting: {args.pod_weighting}")

    first_vec = load_vtk_snapshot(
        path=files[0],
        fields=args.fields,
        cell_stride=args.cell_stride,
        p_scale=args.p_scale,
        dtype=dtype,
    )

    n_samples = len(files)
    n_features = len(first_vec)
    coordinates, field_layout = inspect_vtk_layout(
        path=files[0],
        fields=args.fields,
        cell_stride=args.cell_stride,
    )
    cell_volumes = None
    feature_weights = None
    if args.pod_weighting == "volume":
        cell_volumes = load_cell_volumes(files[0], cell_stride=args.cell_stride)
        if len(cell_volumes) != field_layout["n_cells"]:
            raise ValueError(
                "Cell-volume count does not match sampled cells: "
                f"{len(cell_volumes)} != {field_layout['n_cells']}"
            )
        feature_weights = build_component_weights(cell_volumes, field_layout)
    else:
        feature_weights = np.ones(n_features, dtype=np.float64)

    if field_layout["n_features"] != n_features:
        raise ValueError(
            "Field layout does not match flattened snapshot size: "
            f"{field_layout['n_features']} != {n_features}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)

    matrix_path = args.out.with_suffix(".snapshot_matrix.npy")

    print("")
    print("==> Snapshot matrix")
    print(f"    shape: ({n_samples}, {n_features})")
    print(f"    dtype: {dtype}")
    print(f"    path : {matrix_path}")

    X = np.lib.format.open_memmap(
        matrix_path,
        mode="w+",
        dtype=dtype,
        shape=(n_samples, n_features),
    )

    X[0, :] = first_vec

    for i, path in enumerate(files[1:], start=1):
        vec = load_vtk_snapshot(
            path=path,
            fields=args.fields,
            cell_stride=args.cell_stride,
            p_scale=args.p_scale,
            dtype=dtype,
        )

        if len(vec) != n_features:
            raise ValueError(
                f"Feature length mismatch at {path}: "
                f"{len(vec)} != {n_features}"
            )

        X[i, :] = vec

        if i % 10 == 0 or i == n_samples - 1:
            print(f"    loaded {i + 1}/{n_samples}")

    del X

    time = compute_times(
        n=n_samples,
        time_start=args.time_start,
        time_step=args.time_step,
    )

    (
        coefficients,
        singular_values,
        energy,
        cumulative_energy,
        n_modes,
        spatial_modes,
    ) = method_of_snapshots_pod(
        matrix_path=matrix_path,
        n_samples=n_samples,
        n_features=n_features,
        dtype=dtype,
        chunk_features=args.chunk_features,
        energy_threshold=args.energy_threshold,
        n_modes_max=args.n_modes_max,
        feature_weights=feature_weights,
    )

    print("")
    print("==> POD result")
    print(f"    modes kept : {n_modes}")
    print(f"    energy kept: {cumulative_energy[n_modes - 1]:.8f}")
    print(f"    coeff shape: {coefficients.shape}")

    energy_kept = float(cumulative_energy[n_modes - 1])
    if energy_kept + 1e-12 < args.energy_threshold:
        print("")
        print(
            "WARNING: requested energy threshold was not reached. "
            f"Kept {energy_kept:.8f}, requested {args.energy_threshold:.8f}. "
            "Increase or remove --n-modes-max."
        )

    meta = {
        "schema_version": 2,
        "vtk_glob": args.vtk_glob,
        "n_files": len(files),
        "first_file": str(files[0]),
        "last_file": str(files[-1]),
        "fields": args.fields,
        "cell_stride": args.cell_stride,
        "p_scale": args.p_scale,
        "dtype": str(dtype),
        "n_samples": int(n_samples),
        "n_features": int(n_features),
        "energy_threshold": float(args.energy_threshold),
        "n_modes_max": args.n_modes_max,
        "n_modes": int(n_modes),
        "pod_weighting": args.pod_weighting,
        "inner_product": "sum_c V_c u_c dot v_c" if args.pod_weighting == "volume" else "uniform Euclidean feature inner product",
        "snapshot_matrix": str(matrix_path),
        "field_layout": field_layout,
    }

    np.savez_compressed(
        args.out,
        time=time,
        X=coefficients,
        pod_coefficients=coefficients,
        singular_values=singular_values,
        energy=energy,
        cumulative_energy=cumulative_energy,
        n_modes=np.array([n_modes]),
        spatial_modes=spatial_modes,
        coordinates=coordinates,
        cell_volumes=np.array([] if cell_volumes is None else cell_volumes),
        feature_weights=feature_weights,
        field_shape=np.array(json.dumps(field_layout)),
        fields=np.array(args.fields),
        vtk_files=np.array([str(f) for f in files]),
        meta_json=json.dumps(meta, indent=2),
    )

    print("")
    print("==> Wrote:")
    print(f"    {args.out}")
    print("")
    print("This .npz contains:")
    print("    time")
    print("    X                  # POD coefficients, shape (n_times, n_modes)")
    print("    pod_coefficients")
    print("    energy")
    print("    cumulative_energy")
    print("    singular_values")
    print("    spatial_modes")
    print("    coordinates")
    print("    field_shape")


if __name__ == "__main__":
    main()
