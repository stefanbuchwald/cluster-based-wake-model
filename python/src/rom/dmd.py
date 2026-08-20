"""Exact Dynamic Mode Decomposition.

The routines in this module operate on time-ordered snapshots stored as rows,
while the standard DMD derivation uses snapshots as columns.
"""

from __future__ import annotations

import numpy as np


def _select_rank(
    singular_values: np.ndarray,
    n_modes: int | None,
    energy_threshold: float,
) -> tuple[int, np.ndarray, np.ndarray]:
    """Select a numerically valid truncation rank."""
    if singular_values.ndim != 1 or singular_values.size == 0:
        raise ValueError("singular_values must be a non-empty 1D array")

    if not 0.0 < energy_threshold <= 1.0:
        raise ValueError("energy_threshold must be in (0, 1]")

    tolerance = (
        np.finfo(float).eps
        * max(singular_values.size, 1)
        * float(singular_values[0])
    )
    numerical_rank = int(np.count_nonzero(singular_values > tolerance))

    if numerical_rank == 0:
        raise ValueError("Snapshot matrix has numerical rank zero.")

    energy = singular_values**2
    energy /= energy.sum()
    cumulative_energy = np.cumsum(energy)

    if n_modes is None:
        rank = int(np.searchsorted(cumulative_energy, energy_threshold) + 1)
    else:
        if n_modes < 1:
            raise ValueError("n_modes must be at least 1")
        rank = int(n_modes)

    rank = min(rank, numerical_rank)
    return rank, energy, cumulative_energy


def exact_dmd(
    X: np.ndarray,
    dt: float,
    n_modes: int | None = None,
    energy_threshold: float = 0.9999,
    diameter: float = 1.0,
    reference_velocity: float = 1.0,
) -> dict:
    """Compute exact DMD following Tu et al. (2014).

    Parameters
    ----------
    X:
        Time-ordered snapshot matrix with shape ``(n_samples, n_features)``.
    dt:
        Constant time difference between snapshots.
    n_modes:
        Fixed truncation rank. If omitted, the rank is selected from
        ``energy_threshold``.
    energy_threshold:
        Cumulative singular-value energy used when ``n_modes`` is omitted.
    diameter, reference_velocity:
        Reference quantities used for ``St = f * diameter / reference_velocity``.

    Notes
    -----
    Mode amplitudes are obtained from the standard least-squares fit
    ``Phi @ b ~= x_0`` to the first snapshot. This is intentionally the simple
    standard approach, not optimized DMD amplitudes over the full trajectory.
    """
    snapshots = np.asarray(X)

    if snapshots.ndim != 2:
        raise ValueError(f"X must be 2D, got shape {snapshots.shape}")
    if snapshots.shape[0] < 2:
        raise ValueError("DMD requires at least two snapshots.")
    if snapshots.shape[1] < 1:
        raise ValueError("DMD requires at least one feature.")
    if not np.all(np.isfinite(snapshots)):
        raise ValueError("X contains non-finite values.")
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be a finite positive number.")
    if not np.isfinite(diameter) or diameter <= 0:
        raise ValueError("diameter must be a finite positive number.")
    if not np.isfinite(reference_velocity) or reference_velocity <= 0:
        raise ValueError("reference_velocity must be a finite positive number.")

    X1 = snapshots[:-1].T
    X2 = snapshots[1:].T

    U, singular_values, Vh = np.linalg.svd(X1, full_matrices=False)
    rank, energy, cumulative_energy = _select_rank(
        singular_values=singular_values,
        n_modes=n_modes,
        energy_threshold=energy_threshold,
    )

    U_r = U[:, :rank]
    S_r = singular_values[:rank]
    V_r = Vh.conj().T[:, :rank]
    S_inv = np.diag(1.0 / S_r)

    A_tilde = U_r.conj().T @ X2 @ V_r @ S_inv
    eigenvalues, W = np.linalg.eig(A_tilde)
    modes = X2 @ V_r @ S_inv @ W

    amplitudes, *_ = np.linalg.lstsq(modes, X1[:, 0], rcond=None)

    safe_eigenvalues = eigenvalues.astype(complex, copy=True)
    tiny = np.finfo(float).tiny
    near_zero = np.abs(safe_eigenvalues) < tiny
    safe_eigenvalues[near_zero] = tiny + 0.0j

    omega = np.log(safe_eigenvalues) / float(dt)
    growth_rates = omega.real
    frequencies = omega.imag / (2.0 * np.pi)
    frequencies_St = frequencies * float(diameter) / float(reference_velocity)

    return {
        "eigenvalues": eigenvalues,
        "continuous_eigenvalues": omega,
        "frequencies_St": frequencies_St,
        "growth_rates": growth_rates,
        "amplitudes": amplitudes,
        "modes": modes,
        "dt": float(dt),
        "r": int(rank),
        "singular_values": singular_values,
        "energy": energy,
        "cumulative_energy": cumulative_energy,
    }


def reconstruct_dmd(
    result: dict,
    n_samples: int,
    mode_indices: np.ndarray | list[int] | None = None,
) -> np.ndarray:
    """Reconstruct snapshots from selected DMD modes.

    The returned array has shape ``(n_samples, n_features)``. For real-valued
    input data, small imaginary round-off components are removed.
    """
    if n_samples < 1:
        raise ValueError("n_samples must be at least 1")

    modes = np.asarray(result["modes"])
    amplitudes = np.asarray(result["amplitudes"])
    omega = np.asarray(result["continuous_eigenvalues"])
    dt = float(result["dt"])

    if mode_indices is None:
        indices = np.arange(modes.shape[1], dtype=int)
    else:
        indices = np.asarray(mode_indices, dtype=int)

    if indices.ndim != 1:
        raise ValueError("mode_indices must be one-dimensional")
    if indices.size == 0:
        return np.zeros((n_samples, modes.shape[0]), dtype=float)
    if np.any(indices < 0) or np.any(indices >= modes.shape[1]):
        raise IndexError("mode_indices contains an invalid DMD mode index")

    time = dt * np.arange(n_samples, dtype=float)
    dynamics = amplitudes[indices, None] * np.exp(
        omega[indices, None] * time[None, :]
    )
    reconstruction = modes[:, indices] @ dynamics
    return np.real_if_close(reconstruction.T, tol=1000).real
