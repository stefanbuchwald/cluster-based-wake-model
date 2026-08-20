#!/usr/bin/env python3
"""Run exact DMD on a time-ordered matrix stored in an NPZ file.

Example:
    python3 python/scripts/run_dmd_analysis.py \
        --input data/processed/pod_U_t140_270.npz \
        --out data/figures/dmd_U_t140_270 \
        --dt 0.25
"""

from __future__ import annotations

import argparse
import json
import sys
from itertools import permutations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PYTHON_SRC = PROJECT_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from rom.dmd import exact_dmd, reconstruct_dmd


CNM_FREQUENCIES = {
    "main_cycle": 0.203,
    "sub_harmonic": 0.136,
}


def read_npz_matrix(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Read the shared ``time``/``X`` NPZ schema used by the CNM scripts."""
    if not path.exists():
        raise FileNotFoundError(path)

    with np.load(path) as data:
        if "X" not in data:
            raise ValueError(f"{path} must contain key 'X'")

        X = np.asarray(data["X"], dtype=float)
        time = (
            np.asarray(data["time"], dtype=float)
            if "time" in data
            else np.arange(X.shape[0], dtype=float)
        )

        metadata = {}
        if "meta_json" in data:
            try:
                metadata = json.loads(str(data["meta_json"]))
            except (TypeError, ValueError, json.JSONDecodeError):
                metadata = {}

    if X.ndim != 2:
        raise ValueError(f"Expected X to be 2D, got {X.shape}")
    if time.ndim != 1 or len(time) != X.shape[0]:
        raise ValueError(f"time must have length {X.shape[0]}, got {time.shape}")
    if len(time) < 2:
        raise ValueError("DMD requires at least two time samples.")
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(X)):
        raise ValueError("Input contains non-finite values.")
    if np.any(np.diff(time) <= 0):
        raise ValueError("time must be strictly increasing.")

    return time, X, metadata


def resolve_dt(time: np.ndarray, requested_dt: float | None) -> float:
    """Validate uniform sampling and choose the DMD time step."""
    differences = np.diff(time)
    inferred = float(np.median(differences))

    if not np.allclose(differences, inferred, rtol=1e-6, atol=1e-12):
        raise ValueError(
            "DMD requires uniformly sampled data. "
            "Resample the input time series before running this script."
        )

    if requested_dt is None:
        return inferred
    if requested_dt <= 0 or not np.isfinite(requested_dt):
        raise ValueError("--dt must be a finite positive number.")
    if not np.isclose(requested_dt, inferred, rtol=1e-6, atol=1e-12):
        print(
            "WARNING: --dt differs from the spacing stored in the input: "
            f"{requested_dt:g} vs {inferred:g}."
        )

    return float(requested_dt)


def amplitude_order(result: dict) -> np.ndarray:
    return np.argsort(np.abs(result["amplitudes"]))[::-1]


def leading_positive_modes(
    result: dict,
    count: int,
    frequency_tolerance: float,
) -> list[int]:
    """Return strongest distinct positive-frequency modes by ``|b|``."""
    frequencies = np.asarray(result["frequencies_St"], dtype=float)
    selected: list[int] = []

    for index in amplitude_order(result):
        frequency = float(frequencies[index])
        if frequency <= frequency_tolerance:
            continue
        if any(
            abs(frequency - frequencies[old]) <= frequency_tolerance
            for old in selected
        ):
            continue
        selected.append(int(index))
        if len(selected) == count:
            break

    return selected


def growth_interpretation(growth_rate: float, tolerance: float) -> str:
    if growth_rate > tolerance:
        return (
            "growing - consistent with increasing oscillation amplitude "
            "during part of the analysis window"
        )
    if growth_rate < -tolerance:
        return "decaying"
    return "neutral / limit-cycle-like"


def compare_with_cnm(
    result: dict,
    growth_tolerance: float,
    frequency_tolerance: float,
) -> list[dict]:
    """Match each CNM cycle to the nearest distinct positive DMD frequency."""
    indices = leading_positive_modes(
        result=result,
        count=int(result["r"]),
        frequency_tolerance=frequency_tolerance,
    )

    if not indices:
        return []

    references = list(CNM_FREQUENCIES.items())
    if len(indices) == 1:
        reference = min(
            references,
            key=lambda item: abs(
                float(result["frequencies_St"][indices[0]]) - item[1]
            ),
        )
        assignments = [(reference, indices[0])]
    else:
        candidate_pairs = permutations(indices, len(references))
        best_indices = min(
            candidate_pairs,
            key=lambda assignment: sum(
                abs(float(result["frequencies_St"][idx]) - reference[1])
                for reference, idx in zip(references, assignment)
            ),
        )
        assignments = list(zip(references, best_indices))

    comparison = []
    for (reference_name, reference_frequency), index in assignments:
        frequency = float(result["frequencies_St"][index])
        growth_rate = float(result["growth_rates"][index])
        comparison.append(
            {
                "mode_index": int(index),
                "dmd_St": frequency,
                "cnm_cycle": reference_name,
                "cnm_St": float(reference_frequency),
                "absolute_difference": abs(frequency - reference_frequency),
                "growth_rate": growth_rate,
                "growth_interpretation": growth_interpretation(
                    growth_rate,
                    growth_tolerance,
                ),
            }
        )

    return comparison


def save_eigenvalue_plot(result: dict, out: Path) -> None:
    eigenvalues = np.asarray(result["eigenvalues"])
    amplitudes = np.abs(result["amplitudes"])
    amplitude_scale = np.log10(amplitudes + np.finfo(float).tiny)
    marker_sizes = 30.0 + 170.0 * (
        amplitudes / max(float(amplitudes.max()), np.finfo(float).tiny)
    )

    theta = np.linspace(0.0, 2.0 * np.pi, 400)
    plt.figure(figsize=(7, 6))
    plt.plot(np.cos(theta), np.sin(theta), "--", color="0.55", label="unit circle")
    scatter = plt.scatter(
        eigenvalues.real,
        eigenvalues.imag,
        c=amplitude_scale,
        s=marker_sizes,
        cmap="viridis",
        edgecolors="black",
        linewidths=0.5,
        alpha=0.85,
    )
    plt.axhline(0.0, color="0.75", linewidth=0.8)
    plt.axvline(0.0, color="0.75", linewidth=0.8)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel(r"$\mathrm{Re}(\mu)$")
    plt.ylabel(r"$\mathrm{Im}(\mu)$")
    plt.title("DMD eigenvalues")
    plt.grid(True, alpha=0.25)
    plt.legend(loc="best")
    plt.colorbar(scatter, label=r"$\log_{10}|b|$")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_spectrum_plot(result: dict, out: Path) -> None:
    frequencies = np.asarray(result["frequencies_St"], dtype=float)
    amplitudes = np.abs(result["amplitudes"])

    # Exclude the near-zero-frequency mode: this is the mean/steady field,
    # whose amplitude is orders of magnitude larger than any oscillatory
    # mode and would otherwise dominate the linear y-axis, hiding everything
    # else. It is not part of the oscillatory spectrum we want to compare
    # against the CNM cycle frequencies.
    oscillatory = frequencies > 1e-6
    mask = oscillatory & (frequencies >= -1e-12)

    order = np.argsort(frequencies[mask])
    plotted_frequencies = frequencies[mask][order]
    plotted_amplitudes = amplitudes[mask][order]

    plt.figure(figsize=(9, 4.8))
    markerline, stemlines, baseline = plt.stem(
        plotted_frequencies,
        plotted_amplitudes,
        basefmt=" ",
    )
    plt.setp(markerline, markersize=5, color="tab:blue")
    plt.setp(stemlines, linewidth=1.2, color="tab:blue")
    baseline.set_visible(False)

    plt.axvline(
        CNM_FREQUENCIES["main_cycle"],
        color="tab:red",
        linestyle="--",
        label="CNM main cycle: St = 0.203",
    )
    plt.axvline(
        CNM_FREQUENCIES["sub_harmonic"],
        color="tab:orange",
        linestyle="--",
        label="CNM sub-harmonic: St = 0.136",
    )
    plt.yscale("log")
    plt.xlabel("Strouhal number St")
    plt.ylabel(r"Amplitude $|b|$ (log scale)")
    plt.title("DMD amplitude spectrum (steady mode excluded)")
    plt.grid(True, alpha=0.3, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_reconstruction_plot(
    time: np.ndarray,
    X: np.ndarray,
    reconstruction: np.ndarray,
    out: Path,
) -> None:
    n_coefficients = min(4, X.shape[1])
    plt.figure(figsize=(12, 2.7 * n_coefficients))

    for index in range(n_coefficients):
        axis = plt.subplot(n_coefficients, 1, index + 1)
        axis.plot(time, X[:, index], linewidth=1.1, label=f"Data a{index + 1}")
        axis.plot(
            time,
            reconstruction[:, index],
            linewidth=1.1,
            linestyle="--",
            label=f"DMD a{index + 1}",
        )
        axis.set_ylabel(f"a{index + 1}")
        axis.grid(True, alpha=0.3)
        axis.legend(loc="best")

    plt.xlabel("Time")
    plt.suptitle("Leading coefficient reconstruction", y=0.995)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def mode_record(result: dict, index: int) -> dict:
    eigenvalue = complex(result["eigenvalues"][index])
    return {
        "mode_index": int(index),
        "eigenvalue_real": float(eigenvalue.real),
        "eigenvalue_imag": float(eigenvalue.imag),
        "St": float(result["frequencies_St"][index]),
        "growth_rate": float(result["growth_rates"][index]),
        "amplitude": float(abs(result["amplitudes"][index])),
    }


def parse_args() -> argparse.Namespace:
    example = """example:
  python3 python/scripts/run_dmd_analysis.py \\
      --input data/processed/pod_U_t140_270.npz \\
      --out data/figures/dmd_U_t140_270 \\
      --dt 0.25
"""
    parser = argparse.ArgumentParser(
        description="Exact DMD analysis for a time-ordered NPZ matrix",
        epilog=example,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--n-modes", type=int, default=None)
    parser.add_argument("--energy-threshold", type=float, default=0.9999)
    parser.add_argument("--dt", type=float, default=None)
    parser.add_argument("--reconstruction-modes", type=int, default=10)
    parser.add_argument("--diameter", type=float, default=1.0)
    parser.add_argument("--reference-velocity", type=float, default=1.0)
    parser.add_argument("--growth-tolerance", type=float, default=1e-3)
    parser.add_argument("--frequency-tolerance", type=float, default=1e-6)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    time, X, input_metadata = read_npz_matrix(args.input)
    dt = resolve_dt(time, args.dt)

    args.out.mkdir(parents=True, exist_ok=True)

    print("Loaded input:")
    print(f"  input file : {args.input}")
    print(f"  samples    : {X.shape[0]}")
    print(f"  features   : {X.shape[1]}")
    print(f"  time       : {time.min()} -> {time.max()}")
    print(f"  dt         : {dt}")
    if input_metadata.get("fields"):
        print(f"  fields     : {input_metadata['fields']}")
    if "snapshot_matrix" in input_metadata:
        print("  X meaning  : retained temporal POD coefficients")

    result = exact_dmd(
        X=X,
        dt=dt,
        n_modes=args.n_modes,
        energy_threshold=args.energy_threshold,
        diameter=args.diameter,
        reference_velocity=args.reference_velocity,
    )

    reconstruction_count = min(max(args.reconstruction_modes, 1), result["r"])
    reconstruction_indices = amplitude_order(result)[:reconstruction_count]
    reconstruction = reconstruct_dmd(
        result=result,
        n_samples=len(time),
        mode_indices=reconstruction_indices,
    )

    comparison = compare_with_cnm(
        result=result,
        growth_tolerance=args.growth_tolerance,
        frequency_tolerance=args.frequency_tolerance,
    )

    print("")
    print("DMD:")
    print(f"  rank       : {result['r']}")
    print(
        "  SVD energy : "
        f"{result['cumulative_energy'][result['r'] - 1]:.8f}"
    )
    print(f"  recon modes: {reconstruction_count}")

    print("")
    print("DMD vs CNM frequency comparison (nearest DMD modes):")
    if comparison:
        for number, item in enumerate(comparison, start=1):
            label = item["cnm_cycle"].replace("_", " ")
            print(
                f"  DMD mode {number}: St = {item['dmd_St']:.6f} "
                f"(CNM {label}: St = {item['cnm_St']:.3f}, "
                f"diff = {item['absolute_difference']:.6f})"
            )
            print(
                f"    growth rate = {item['growth_rate']:.6e}: "
                f"{item['growth_interpretation']}"
            )
    else:
        print("  No non-trivial positive-frequency DMD mode found.")

    np.save(args.out / "dmd_eigenvalues.npy", result["eigenvalues"])
    np.save(
        args.out / "dmd_continuous_eigenvalues.npy",
        result["continuous_eigenvalues"],
    )
    np.save(args.out / "dmd_frequencies_St.npy", result["frequencies_St"])
    np.save(args.out / "dmd_growth_rates.npy", result["growth_rates"])
    np.save(args.out / "dmd_amplitudes.npy", result["amplitudes"])
    np.save(args.out / "dmd_modes.npy", result["modes"])
    np.save(args.out / "dmd_reconstruction.npy", reconstruction)

    save_eigenvalue_plot(result, args.out / "09_dmd_eigenvalues.png")
    save_spectrum_plot(result, args.out / "10_dmd_spectrum.png")
    save_reconstruction_plot(
        time=time,
        X=X,
        reconstruction=reconstruction,
        out=args.out / "11_dmd_mode_reconstruction.png",
    )

    top_indices = amplitude_order(result)[: min(10, result["r"])]
    summary = {
        "input": str(args.input),
        "input_matrix": (
            "POD coefficients" if "snapshot_matrix" in input_metadata else "X"
        ),
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "time_min": float(time.min()),
        "time_max": float(time.max()),
        "dt": float(dt),
        "rank": int(result["r"]),
        "energy_threshold": float(args.energy_threshold),
        "svd_cumulative_energy_kept": float(
            result["cumulative_energy"][result["r"] - 1]
        ),
        "reconstruction_modes": reconstruction_indices.tolist(),
        "top_modes_by_amplitude": [
            mode_record(result, int(index)) for index in top_indices
        ],
        "cnm_frequency_comparison": comparison,
        "cnm_reference_frequencies": CNM_FREQUENCIES,
        "outputs": {
            "eigenvalues_plot": "09_dmd_eigenvalues.png",
            "spectrum_plot": "10_dmd_spectrum.png",
            "reconstruction_plot": "11_dmd_mode_reconstruction.png",
        },
    }

    with (args.out / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print("")
    print("DMD analysis complete.")
    print(f"  output: {args.out}")


if __name__ == "__main__":
    main()
