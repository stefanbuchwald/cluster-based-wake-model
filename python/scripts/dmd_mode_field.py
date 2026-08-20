#!/usr/bin/env python3
"""Project a DMD mode from POD space back to a physical field slice."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def load_spatial_pod_data(npz_path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Load schema-v2 spatial POD metadata."""
    if not npz_path.exists():
        raise FileNotFoundError(npz_path)

    with np.load(npz_path) as data:
        required = {"spatial_modes", "coordinates", "field_shape"}
        missing = sorted(required.difference(data.files))
        if missing:
            raise ValueError(
                f"{npz_path} uses the old POD schema and is missing {missing}. "
                "Re-run field_pod_from_vtk.py to create a schema-v2 NPZ."
            )

        spatial_modes = np.asarray(data["spatial_modes"])
        coordinates = np.asarray(data["coordinates"], dtype=float)
        field_layout = json.loads(str(data["field_shape"]))

    if spatial_modes.ndim != 2:
        raise ValueError("spatial_modes must have shape (n_features, n_modes)")
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("coordinates must have shape (n_cells, 3)")
    if field_layout.get("n_cells") != len(coordinates):
        raise ValueError("field_shape n_cells does not match coordinates")
    if field_layout.get("n_features") != spatial_modes.shape[0]:
        raise ValueError("field_shape n_features does not match spatial_modes")

    return spatial_modes, coordinates, field_layout


def load_dmd_data(
    dmd_dir: Path,
    summary_path: Path,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Load DMD modes, frequencies, and summary metadata."""
    modes_path = dmd_dir / "dmd_modes.npy"
    frequencies_path = dmd_dir / "dmd_frequencies_St.npy"

    if not modes_path.exists() or not frequencies_path.exists():
        raise FileNotFoundError(
            f"Expected {modes_path.name} and {frequencies_path.name} in {dmd_dir}"
        )
    if not summary_path.exists():
        raise FileNotFoundError(summary_path)

    modes = np.load(modes_path)
    frequencies = np.asarray(np.load(frequencies_path), dtype=float)
    with summary_path.open("r", encoding="utf-8") as handle:
        summary = json.load(handle)

    if modes.ndim != 2:
        raise ValueError("dmd_modes.npy must be a 2D matrix")
    if frequencies.shape != (modes.shape[1],):
        raise ValueError("DMD frequency count does not match DMD mode count")

    return modes, frequencies, summary


def select_mode(frequencies: np.ndarray, target_St: float) -> int:
    """Select the positive-frequency DMD mode nearest to ``target_St``."""
    candidates = np.flatnonzero(frequencies > 0.0)
    if candidates.size == 0:
        candidates = np.arange(len(frequencies))
    return int(candidates[np.argmin(np.abs(frequencies[candidates] - target_St))])


def select_field_layout(field_layout: dict, field_name: str) -> dict:
    for entry in field_layout.get("fields", []):
        if entry.get("name") == field_name:
            return entry
    available = [entry.get("name") for entry in field_layout.get("fields", [])]
    raise ValueError(f"Field {field_name!r} not found. Available fields: {available}")


def component_values(
    field_values: np.ndarray,
    component: str,
) -> tuple[np.ndarray, str, bool]:
    """Extract a real component or magnitude from a complex field."""
    n_components = field_values.shape[1]
    labels = {"x": 0, "y": 1, "z": 2}

    if component == "magnitude":
        values = np.linalg.norm(field_values.real, axis=1)
        return values, "magnitude of Re(phi)", False

    if component in labels:
        component_index = labels[component]
    else:
        try:
            component_index = int(component)
        except ValueError as exc:
            raise ValueError(
                "--component must be x, y, z, magnitude, or an integer"
            ) from exc

    if component_index < 0 or component_index >= n_components:
        raise ValueError(
            f"Component {component!r} is invalid for a field with "
            f"{n_components} component(s)."
        )

    return (
        field_values[:, component_index].real,
        f"Re(phi), component {component}",
        True,
    )


def extract_z_slice(
    coordinates: np.ndarray,
    values: np.ndarray,
    requested_z: float | None,
    tolerance: float | None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Extract the nearest available cell-center plane in the z direction."""
    z_values = coordinates[:, 2]
    unique_z = np.unique(z_values)
    target = (
        0.5 * (float(unique_z.min()) + float(unique_z.max()))
        if requested_z is None
        else float(requested_z)
    )
    selected_z = float(unique_z[np.argmin(np.abs(unique_z - target))])

    if tolerance is None:
        if len(unique_z) > 1:
            spacing = float(np.min(np.diff(unique_z)))
            tolerance = max(0.1 * spacing, 1e-10)
        else:
            tolerance = 1e-10

    mask = np.abs(z_values - selected_z) <= tolerance
    if np.count_nonzero(mask) < 3:
        raise ValueError(
            f"Only {np.count_nonzero(mask)} cells found near z={selected_z}. "
            "Increase --z-tolerance."
        )

    return coordinates[mask], values[mask], selected_z


def save_mode_field_plot(
    coordinates: np.ndarray,
    values: np.ndarray,
    actual_St: float,
    field_name: str,
    component_label: str,
    diverging: bool,
    z_slice: float,
    out: Path,
    cylinder_center: tuple[float, float],
    cylinder_diameter: float,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    triangulation = mtri.Triangulation(coordinates[:, 0], coordinates[:, 1])

    if diverging:
        bound = float(np.max(np.abs(values)))
        if bound <= 0:
            bound = 1.0
        levels = np.linspace(-bound, bound, 81)
        cmap = "RdBu_r"
    else:
        value_min = float(np.min(values))
        value_max = float(np.max(values))
        if np.isclose(value_min, value_max):
            value_max = value_min + 1.0
        levels = np.linspace(value_min, value_max, 81)
        cmap = "viridis"

    fig, axis = plt.subplots(figsize=(12, 4.8))
    contour = axis.tricontourf(
        triangulation,
        values,
        levels=levels,
        cmap=cmap,
        extend="both",
    )
    cylinder = plt.Circle(
        cylinder_center,
        radius=0.5 * cylinder_diameter,
        facecolor="white",
        edgecolor="black",
        linewidth=1.0,
        zorder=10,
    )
    axis.add_patch(cylinder)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("x / D")
    axis.set_ylabel("y / D")
    axis.set_title(
        f"DMD mode, St = {actual_St:.4f}, {component_label}, "
        f"{field_name}, z = {z_slice:.4g}"
    )
    fig.colorbar(contour, ax=axis, label="DMD mode amplitude")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a spatial DMD mode on one spanwise slice"
    )
    parser.add_argument("--npz", required=True, type=Path)
    parser.add_argument("--dmd-summary", required=True, type=Path)
    parser.add_argument("--dmd-dir", required=True, type=Path)
    parser.add_argument("--mode-st", required=True, type=float)
    parser.add_argument("--field", default="U")
    parser.add_argument("--component", default="x")
    parser.add_argument("--z-slice", type=float, default=None)
    parser.add_argument("--z-tolerance", type=float, default=None)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--cylinder-center",
        nargs=2,
        type=float,
        metavar=("X", "Y"),
        default=(0.0, 0.0),
    )
    parser.add_argument("--cylinder-diameter", type=float, default=1.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    spatial_modes, coordinates, field_layout = load_spatial_pod_data(args.npz)
    dmd_modes, frequencies, summary = load_dmd_data(
        dmd_dir=args.dmd_dir,
        summary_path=args.dmd_summary,
    )

    if spatial_modes.shape[1] != dmd_modes.shape[0]:
        raise ValueError(
            "POD/DMD mode mismatch: spatial_modes has "
            f"{spatial_modes.shape[1]} POD modes, while dmd_modes expects "
            f"{dmd_modes.shape[0]}. Re-run DMD with the same NPZ file."
        )

    mode_index = select_mode(frequencies, args.mode_st)
    actual_St = float(frequencies[mode_index])
    physical_mode = spatial_modes @ dmd_modes[:, mode_index]

    field_entry = select_field_layout(field_layout, args.field)
    start = int(field_entry["feature_start"])
    stop = int(field_entry["feature_stop"])
    n_components = int(field_entry["n_components"])
    field_values = physical_mode[start:stop].reshape(
        int(field_layout["n_cells"]),
        n_components,
    )

    values, component_label, diverging = component_values(
        field_values,
        args.component,
    )
    slice_coordinates, slice_values, selected_z = extract_z_slice(
        coordinates=coordinates,
        values=values,
        requested_z=args.z_slice,
        tolerance=args.z_tolerance,
    )

    save_mode_field_plot(
        coordinates=slice_coordinates,
        values=slice_values,
        actual_St=actual_St,
        field_name=args.field,
        component_label=component_label,
        diverging=diverging,
        z_slice=selected_z,
        out=args.out,
        cylinder_center=tuple(args.cylinder_center),
        cylinder_diameter=args.cylinder_diameter,
    )

    print("Spatial DMD mode field:")
    print(f"  requested St : {args.mode_st:.6f}")
    print(f"  matched St   : {actual_St:.6f}")
    print(f"  mode index   : {mode_index}")
    print(f"  DMD rank     : {summary.get('rank', 'unknown')}")
    print(f"  field        : {args.field}, component={args.component}")
    print(f"  z slice      : {selected_z:.6g}")
    print(f"  cells        : {len(slice_coordinates)}")
    print(f"  output       : {args.out}")


if __name__ == "__main__":
    main()
