#!/usr/bin/env python3
"""Generate a small periodic trajectory for the CNM/DMD quickstart.

This synthetic example exercises the analysis pipeline; it is not CFD data.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=Path("data/processed/demo_periodic.npz")
    )
    args = parser.parse_args()

    dt = 0.25
    frequency = 0.2
    time = dt * np.arange(241, dtype=float)
    phase = 2.0 * np.pi * frequency * time
    states = np.column_stack(
        (np.cos(phase), np.sin(phase), 0.3 * np.cos(2 * phase), 0.3 * np.sin(2 * phase))
    )
    metadata = {
        "source": "synthetic periodic trajectory; not CFD data",
        "frequency": frequency,
        "harmonic_frequency": 2 * frequency,
        "diameter": 1.0,
        "reference_velocity": 1.0,
        "sample_dt": dt,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, time=time, X=states, meta_json=json.dumps(metadata))
    print(f"Wrote {args.out}: {len(time)} samples, {states.shape[1]} features, dt={dt}")
    print("Synthetic frequencies: St=0.2 and St=0.4 for D=U_inf=1.")


if __name__ == "__main__":
    main()
