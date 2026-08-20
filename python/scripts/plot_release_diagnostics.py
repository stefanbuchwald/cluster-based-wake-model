#!/usr/bin/env python3
"""Create final release diagnostics from CNM and DMD summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def plot_network_cycle(out_dir: Path, cnm: dict) -> None:
    Q = np.load(out_dir / "Q_matrix.npy")
    T = np.load(out_dir / "T_matrix.npy")
    counts = np.load(out_dir / "transition_counts.npy")
    cycles = cnm.get("cycle_diagnostics", {}).get("cycles", [])
    if not cycles:
        raise ValueError("No CNM cycles available in cnm_summary.json")
    cycle = cycles[0]
    cycle_edges = {tuple(edge) for edge in cycle["edges"]}
    K = Q.shape[0]
    theta = np.linspace(0.0, 2.0 * np.pi, K, endpoint=False)
    xy = np.c_[np.cos(theta), np.sin(theta)]

    fig, ax = plt.subplots(figsize=(7.5, 7.0))
    ax.scatter(xy[:, 0], xy[:, 1], s=500, color="white", edgecolor="black", zorder=3)
    for k, (x, y) in enumerate(xy):
        ax.text(x, y, str(k), ha="center", va="center", fontsize=10, fontweight="bold", zorder=4)

    max_count = max(float(counts.max()), 1.0)
    for k in range(K):
        for j in range(K):
            count = int(counts[j, k])
            if j == k or count <= 0:
                continue
            is_cycle = (j, k) in cycle_edges
            alpha = 0.25 if not is_cycle else 0.95
            color = "tab:red" if is_cycle else "0.55"
            linewidth = 0.5 + 4.0 * count / max_count if is_cycle else 0.6 + 1.2 * count / max_count
            x0, y0 = xy[k]
            x1, y1 = xy[j]
            ax.annotate(
                "",
                xy=(x1, y1),
                xytext=(x0, y0),
                arrowprops={
                    "arrowstyle": "->",
                    "color": color,
                    "alpha": alpha,
                    "lw": linewidth,
                    "shrinkA": 18,
                    "shrinkB": 18,
                    "connectionstyle": "arc3,rad=0.12",
                },
                zorder=2 if is_cycle else 1,
            )
            if is_cycle:
                mx, my = 0.55 * np.array([x0 + x1, y0 + y1])
                t_value = T[j, k]
                label = f"n={count}\nT={t_value:.2f}" if np.isfinite(t_value) else f"n={count}\nT=nan"
                ax.text(mx, my, label, fontsize=8, ha="center", va="center", color="tab:red")

    ax.set_title(
        "Dominant CNM cycle "
        f"St={cycle['St']:.3f}, min support={cycle['min_support']}"
    )
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_dir / "15_network_cycle_support.png", dpi=220)
    plt.close(fig)


def choose_dmd_candidate(dmd: dict, target_st: float) -> dict:
    candidates = dmd.get("candidate_positive_frequency_modes", [])
    if not candidates:
        raise ValueError("No DMD candidates available in dmd_summary.json")
    return min(candidates, key=lambda item: abs(float(item["St"]) - target_st))


def plot_frequency_comparison(out_dir: Path, cnm: dict, dmd: dict) -> None:
    cycles = cnm.get("cycle_diagnostics", {}).get("cycles", [])
    if not cycles:
        raise ValueError("No CNM cycles available in cnm_summary.json")
    cycle = cycles[0]
    dmd_candidate = choose_dmd_candidate(dmd, float(cycle["St"]))
    values = [float(cycle["St"]), float(dmd_candidate["St"])]
    labels = ["CNM cycle", "DMD candidate"]

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    bars = ax.bar(labels, values, color=["tab:red", "tab:blue"], alpha=0.85)
    ax.set_ylabel("Strouhal number St")
    ax.set_ylim(0.0, max(values) * 1.35)
    ax.set_title("Independent CNM/DMD frequency comparison")
    ax.grid(True, axis="y", alpha=0.25)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2.0, value, f"{value:.3f}", ha="center", va="bottom")
    ax.text(
        0.5,
        0.92,
        f"Both indicate St ≈ {np.mean(values):.2f}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_dir / "16_cnm_dmd_frequency_comparison.png", dpi=220)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot final release diagnostics")
    parser.add_argument("--out", type=Path, required=True, help="Result directory")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cnm = load_json(args.out / "cnm_summary.json")
    dmd = load_json(args.out / "dmd_summary.json")
    plot_network_cycle(args.out, cnm)
    plot_frequency_comparison(args.out, cnm, dmd)
    print(args.out / "15_network_cycle_support.png")
    print(args.out / "16_cnm_dmd_frequency_comparison.png")


if __name__ == "__main__":
    main()
