"""POD helpers for finite-volume snapshot data."""

from __future__ import annotations

import numpy as np


def build_component_weights(cell_volumes: np.ndarray, field_layout: dict) -> np.ndarray:
    """Return one scalar weight per flattened feature.

    The repository flattens fields in field-major order. For a vector field such
    as velocity, the same finite-volume cell volume is repeated for each vector
    component belonging to that cell.
    """
    volumes = np.asarray(cell_volumes, dtype=float)
    if volumes.ndim != 1 or volumes.size == 0:
        raise ValueError("cell_volumes must be a non-empty 1D array")
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0):
        raise ValueError("cell_volumes must be finite and positive")

    parts = []
    for field in field_layout.get("fields", []):
        n_components = int(field["n_components"])
        if n_components < 1:
            raise ValueError("field n_components must be positive")
        parts.append(np.repeat(volumes, n_components))

    if not parts:
        raise ValueError("field_layout must contain at least one field")

    weights = np.concatenate(parts)
    expected = int(field_layout.get("n_features", weights.size))
    if weights.size != expected:
        raise ValueError(
            f"weight length does not match feature layout: {weights.size} != {expected}"
        )
    return weights


def weighted_snapshot_pod(
    X: np.ndarray,
    weights: np.ndarray | None = None,
    energy_threshold: float = 0.9999,
    n_modes_max: int | None = None,
) -> dict:
    """Compute POD with the weighted inner product ``sum_j w_j u_j v_j``.

    ``X`` stores one snapshot per row. The returned ``spatial_modes`` are in
    physical coordinates, so reconstructions use ``coefficients @ modes.T``.
    """
    snapshots = np.asarray(X, dtype=float)
    if snapshots.ndim != 2:
        raise ValueError(f"X must be 2D, got {snapshots.shape}")
    if snapshots.shape[0] < 2:
        raise ValueError("POD requires at least two snapshots")
    if not np.all(np.isfinite(snapshots)):
        raise ValueError("X contains non-finite values")

    n_samples, n_features = snapshots.shape
    if weights is None:
        weights_array = np.ones(n_features, dtype=float)
        weighting = "uniform"
    else:
        weights_array = np.asarray(weights, dtype=float)
        weighting = "custom"
    if weights_array.shape != (n_features,):
        raise ValueError(f"weights must have shape ({n_features},), got {weights_array.shape}")
    if not np.all(np.isfinite(weights_array)) or np.any(weights_array <= 0):
        raise ValueError("weights must be finite and positive")
    if not 0.0 < energy_threshold <= 1.0:
        raise ValueError("energy_threshold must be in (0, 1]")

    mean = snapshots.mean(axis=0)
    centered = snapshots - mean
    weighted_centered = centered * np.sqrt(weights_array)[None, :]
    U, singular_values, Vt_weighted = np.linalg.svd(weighted_centered, full_matrices=False)

    if np.allclose(singular_values, 0.0):
        raise ValueError("All singular values are zero; no fluctuation energy found")

    eigenvalues = singular_values**2 / float(n_samples - 1)
    energy = singular_values**2 / np.sum(singular_values**2)
    cumulative_energy = np.cumsum(energy)
    n_modes = int(np.searchsorted(cumulative_energy, energy_threshold) + 1)
    if n_modes_max is not None:
        n_modes = min(n_modes, int(n_modes_max))

    tolerance = np.finfo(float).eps * max(snapshots.shape) * singular_values[0]
    numerical_rank = int(np.count_nonzero(singular_values > tolerance))
    n_modes = max(1, min(n_modes, numerical_rank))

    coefficients = U[:, :n_modes] * singular_values[:n_modes]
    spatial_modes = Vt_weighted[:n_modes, :].T / np.sqrt(weights_array)[:, None]

    return {
        "mean": mean,
        "spatial_modes": spatial_modes,
        "coefficients": coefficients,
        "singular_values": singular_values,
        "eigenvalues": eigenvalues,
        "energy": energy,
        "cumulative_energy": cumulative_energy,
        "n_modes": int(n_modes),
        "weights": weights_array,
        "weighting": weighting,
    }
