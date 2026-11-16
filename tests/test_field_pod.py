import numpy as np
import pytest

from field_pod_from_vtk import method_of_snapshots_pod


@pytest.mark.parametrize("storage_dtype", ["float32", "float64"])
def test_snapshot_pod_reconstructs_readonly_input(tmp_path, storage_dtype):
    snapshots = np.random.default_rng(4).normal(size=(7, 4)).astype(storage_dtype)
    matrix_path = tmp_path / "snapshots.npy"
    np.save(matrix_path, snapshots)
    weights = np.array([1.0, 4.0, 0.5, 2.0])

    coefficients, _, _, _, n_modes, modes = method_of_snapshots_pod(
        matrix_path=matrix_path,
        n_samples=7,
        n_features=4,
        dtype=np.dtype(storage_dtype),
        chunk_features=2,
        energy_threshold=1.0,
        n_modes_max=4,
        feature_weights=weights,
    )

    assert n_modes == 4
    tolerance = 2e-6 if storage_dtype == "float32" else 1e-12
    centered = snapshots.astype(float) - snapshots.mean(axis=0, dtype=float)
    np.testing.assert_allclose(coefficients @ modes.T, centered, atol=tolerance)
    np.testing.assert_allclose(modes.T @ (weights[:, None] * modes), np.eye(4), atol=tolerance)
    np.testing.assert_array_equal(np.load(matrix_path), snapshots)
