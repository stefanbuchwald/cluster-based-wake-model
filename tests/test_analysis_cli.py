import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_cnm_and_dmd_analyze_synthetic_oscillation(tmp_path):
    input_path = tmp_path / "pod.npz"
    output_path = tmp_path / "results"
    environment = {
        **os.environ,
        "MPLBACKEND": "Agg",
        "MPLCONFIGDIR": str(tmp_path / "matplotlib"),
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
    }
    for script, args in [
        ("scripts/generate_demo.py", ["--out", str(input_path)]),
        ("python/scripts/run_cnm_analysis.py", [
            "--input", str(input_path), "--out", str(output_path),
            "--input-type", "npz", "--K", "8", "--n-runs", "3",
        ]),
        ("python/scripts/run_dmd_analysis.py", [
            "--input", str(input_path), "--out", str(output_path),
        ]),
        ("python/scripts/plot_release_diagnostics.py", ["--out", str(output_path)]),
    ]:
        result = subprocess.run(
            [
                sys.executable,
                str(PROJECT_ROOT / script),
                *args,
            ],
            cwd=PROJECT_ROOT,
            env=environment,
            capture_output=True,
            text=True,
            timeout=90,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    cnm = json.loads((output_path / "cnm_summary.json").read_text())
    assert cnm["pod_retained_modes"] == 4
    assert cnm["sample_dt"] == 0.25
    assert cnm["cycle_diagnostics"]["cycles"]
    assert np.isfinite(cnm["cycle_diagnostics"]["cycles"][0]["St"])
    for name in ["04_Q_matrix.png", "15_network_cycle_support.png", "16_cnm_dmd_frequency_comparison.png"]:
        assert (output_path / name).is_file()
    np.testing.assert_allclose(
        np.load(output_path / "P_matrix.npy").sum(axis=0), 1.0
    )
    dmd = json.loads((output_path / "dmd_summary.json").read_text())
    assert dmd["full_retained_relative_frobenius_error"] < 1e-8
    for expected_St in [0.2, 0.4]:
        assert any(
            abs(candidate["St"] - expected_St) < 1e-8
            for candidate in dmd["candidate_positive_frequency_modes"]
        )
