#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def read_force_coeffs(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Could not find forceCoeffs file: {path}")

    header = None
    rows = []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()

            if not s:
                continue

            if s.startswith("#"):
                candidate = s.lstrip("#").strip().split()
                if candidate and candidate[0].lower() == "time":
                    header = candidate
                continue

            parts = s.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue

    if not rows:
        raise ValueError(f"No numeric data found in {path}")

    n_cols = len(rows[0])

    if header is None or len(header) != n_cols:
        fallback = ["Time", "Cm", "Cd", "Cl", "Cl(f)", "Cl(r)"]
        if n_cols <= len(fallback):
            header = fallback[:n_cols]
        else:
            header = ["Time"] + [f"col_{i}" for i in range(1, n_cols)]

    df = pd.DataFrame(rows, columns=header)

    if "Time" not in df.columns:
        raise ValueError(f"Expected Time column. Found: {list(df.columns)}")

    df = df.drop_duplicates(subset=["Time"], keep="last")
    df = df.sort_values("Time").reset_index(drop=True)

    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--t-min", type=float, default=None)
    parser.add_argument("--t-max", type=float, default=None)
    parser.add_argument("--no-auto-ylim", action="store_true")
    args = parser.parse_args()

    df = read_force_coeffs(args.file)

    if args.t_min is not None:
        df = df[df["Time"] >= args.t_min]

    if args.t_max is not None:
        df = df[df["Time"] <= args.t_max]

    if df.empty:
        raise ValueError("No data left after time filtering.")

    out = args.out
    out.parent.mkdir(parents=True, exist_ok=True)

    print("Loaded forceCoeffs:")
    print(f"  rows   : {len(df)}")
    print(f"  columns: {list(df.columns)}")
    print(f"  time   : {df['Time'].min()} -> {df['Time'].max()}")

    cols = [c for c in ["Cd", "Cl", "Cm"] if c in df.columns]

    if not cols:
        cols = [c for c in df.columns if c != "Time"]

    for col in cols:
        print(
            f"  {col}: min={df[col].min():.6g}, "
            f"max={df[col].max():.6g}, "
            f"last={df[col].iloc[-1]:.6g}"
        )

    plt.figure(figsize=(10, 5))

    for col in cols:
        plt.plot(df["Time"], df[col], linewidth=1.3, label=col)

    plt.xlabel("Time")
    plt.ylabel("Coefficient")
    plt.title("OpenFOAM force coefficients")
    plt.grid(True, alpha=0.3)
    plt.legend()

    if not args.no_auto_ylim:
        values = df[cols].to_numpy().ravel()
        q_low = pd.Series(values).quantile(0.01)
        q_high = pd.Series(values).quantile(0.99)
        margin = 0.15 * max(abs(q_high - q_low), 1e-12)
        plt.ylim(q_low - margin, q_high + margin)

    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()

    print(f"Wrote plot: {out}")


if __name__ == "__main__":
    main()
