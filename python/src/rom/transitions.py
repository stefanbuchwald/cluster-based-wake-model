"""Transition matrices and lightweight network dynamics for CMM/CNM."""

from __future__ import annotations

import math

import numpy as np


def transition_matrix_P(labels: np.ndarray, K: int) -> tuple[np.ndarray, np.ndarray]:
    labels = np.asarray(labels, dtype=int)
    counts = np.zeros((K, K), dtype=int)
    for old, new in zip(labels[:-1], labels[1:]):
        counts[int(new), int(old)] += 1
    col_sum = counts.sum(axis=0, keepdims=True)
    P = np.divide(counts, col_sum, out=np.zeros((K, K), dtype=float), where=col_sum > 0)
    return P, counts


def run_length_encoding(labels: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    labels = np.asarray(labels, dtype=int)
    if labels.size == 0:
        raise ValueError("labels must not be empty")
    starts = [0]
    values = [int(labels[0])]
    for i in range(1, labels.size):
        if labels[i] != labels[i - 1]:
            starts.append(i)
            values.append(int(labels[i]))
    starts_array = np.asarray(starts, dtype=int)
    ends_array = np.r_[starts_array[1:] - 1, labels.size - 1].astype(int)
    return np.asarray(values, dtype=int), starts_array, ends_array


def cnm_matrices(labels: np.ndarray, time: np.ndarray, K: int, exclude_censored: bool = True) -> dict:
    labels = np.asarray(labels, dtype=int)
    time = np.asarray(time, dtype=float)
    if labels.shape != time.shape:
        raise ValueError("labels and time must have the same 1D shape")
    if labels.ndim != 1 or labels.size < 2:
        raise ValueError("at least two labels are required")
    if np.any(labels < 0) or np.any(labels >= K):
        raise ValueError("labels contain cluster ids outside [0, K)")
    if np.any(np.diff(time) <= 0):
        raise ValueError("time must be strictly increasing")

    P, P_counts = transition_matrix_P(labels, K)
    run_labels, run_start, run_end = run_length_encoding(labels)
    dt = float(np.median(np.diff(time)))
    durations = np.maximum(time[run_end] - time[run_start] + dt, dt)
    complete = np.ones(run_labels.size, dtype=bool)
    if exclude_censored and complete.size:
        complete[0] = False
        complete[-1] = False

    Q_counts = np.zeros((K, K), dtype=int)
    T_sums = np.zeros((K, K), dtype=float)
    T_counts = np.zeros((K, K), dtype=int)

    for r in range(run_labels.size - 1):
        k = int(run_labels[r])
        j = int(run_labels[r + 1])
        if j == k:
            continue
        Q_counts[j, k] += 1
        if (not exclude_censored) or (complete[r] and complete[r + 1]):
            T_sums[j, k] += 0.5 * float(durations[r]) + 0.5 * float(durations[r + 1])
            T_counts[j, k] += 1

    col_sum = Q_counts.sum(axis=0, keepdims=True)
    Q = np.divide(Q_counts, col_sum, out=np.zeros((K, K), dtype=float), where=col_sum > 0)
    T = np.divide(T_sums, T_counts, out=np.full((K, K), np.nan), where=T_counts > 0)

    return {
        "P": P,
        "P_counts": P_counts,
        "Q": Q,
        "Q_counts": Q_counts,
        "T": T,
        "T_counts": T_counts,
        "run_labels": run_labels,
        "run_start_idx": run_start,
        "run_end_idx": run_end,
        "run_durations": durations,
        "run_complete": complete,
        "n_runs": int(run_labels.size),
        "n_complete_runs": int(np.count_nonzero(complete)),
        "left_censored": bool(exclude_censored and run_labels.size > 0),
        "right_censored": bool(exclude_censored and run_labels.size > 1),
    }


def _normalize_column(column: np.ndarray) -> np.ndarray:
    probs = np.asarray(column, dtype=float).copy()
    probs[~np.isfinite(probs)] = 0.0
    probs = np.clip(probs, 0.0, None)
    total = float(probs.sum())
    if total <= 0.0:
        return probs
    return probs / total


def deterministic_cnm(Q: np.ndarray, T: np.ndarray, centers: np.ndarray, start_cluster: int, n_transitions: int) -> dict:
    Q = np.asarray(Q, dtype=float)
    T = np.asarray(T, dtype=float)
    centers = np.asarray(centers, dtype=float)
    labels = [int(start_cluster)]
    event_time = [0.0]
    current = int(start_cluster)
    for _ in range(n_transitions):
        probs = _normalize_column(Q[:, current])
        probs[current] = 0.0
        if probs.sum() <= 0.0:
            break
        nxt = int(np.argmax(probs))
        duration = float(T[nxt, current])
        if not np.isfinite(duration) or duration <= 0.0:
            break
        labels.append(nxt)
        event_time.append(event_time[-1] + duration)
        current = nxt
    labels_array = np.asarray(labels, dtype=int)
    return {"labels": labels_array, "event_time": np.asarray(event_time), "states": centers[labels_array]}


def probabilistic_cmm(P: np.ndarray, centers: np.ndarray, start_cluster: int, n_steps: int, dt: float, rng: np.random.Generator) -> dict:
    P = np.asarray(P, dtype=float)
    centers = np.asarray(centers, dtype=float)
    K = P.shape[0]
    if P.shape != (K, K):
        raise ValueError("P must be square")
    if dt <= 0 or not np.isfinite(dt):
        raise ValueError("dt must be positive")
    labels = np.empty(n_steps + 1, dtype=int)
    labels[0] = int(start_cluster)
    for i in range(n_steps):
        current = int(labels[i])
        probs = _normalize_column(P[:, current])
        labels[i + 1] = current if probs.sum() <= 0.0 else int(rng.choice(K, p=probs))
    return {
        "labels": labels,
        "time": dt * np.arange(n_steps + 1, dtype=float),
        "states": centers[labels],
    }


def probabilistic_cnm(Q: np.ndarray, T: np.ndarray, centers: np.ndarray, start_cluster: int, n_transitions: int, rng: np.random.Generator) -> dict:
    Q = np.asarray(Q, dtype=float)
    T = np.asarray(T, dtype=float)
    centers = np.asarray(centers, dtype=float)
    K = Q.shape[0]
    labels = [int(start_cluster)]
    event_time = [0.0]
    current = int(start_cluster)
    for _ in range(n_transitions):
        probs = _normalize_column(Q[:, current])
        probs[current] = 0.0
        if probs.sum() <= 0.0:
            break
        probs = probs / probs.sum()
        nxt = int(rng.choice(K, p=probs))
        duration = float(T[nxt, current])
        if not np.isfinite(duration) or duration <= 0.0:
            break
        labels.append(nxt)
        event_time.append(event_time[-1] + duration)
        current = nxt
    labels_array = np.asarray(labels, dtype=int)
    return {"labels": labels_array, "event_time": np.asarray(event_time), "states": centers[labels_array]}


def _canonical_cycle(cycle: list[int]) -> tuple[int, ...]:
    rotations = [tuple(cycle[i:] + cycle[:i]) for i in range(len(cycle))]
    return min(rotations)


def enumerate_cycles(
    Q: np.ndarray,
    T: np.ndarray,
    counts: np.ndarray,
    min_count: int = 1,
    max_length: int | None = None,
    diameter: float = 1.0,
    reference_velocity: float = 1.0,
) -> list[dict]:
    Q = np.asarray(Q, dtype=float)
    T = np.asarray(T, dtype=float)
    counts = np.asarray(counts, dtype=int)
    K = Q.shape[0]
    if max_length is None:
        max_length = K
    adjacency = {
        k: [j for j in range(K) if j != k and counts[j, k] >= min_count and Q[j, k] > 0]
        for k in range(K)
    }
    seen: set[tuple[int, ...]] = set()
    cycles: list[dict] = []

    def dfs(start: int, current: int, path: list[int]) -> None:
        if len(path) > max_length:
            return
        for nxt in adjacency[current]:
            if nxt == start and len(path) >= 2:
                canonical = _canonical_cycle(path.copy())
                if canonical in seen:
                    continue
                seen.add(canonical)
                edges = [(path[(i + 1) % len(path)], path[i]) for i in range(len(path))]
                edge_times = [float(T[j, k]) for j, k in edges]
                if not all(np.isfinite(edge_times)) or any(t <= 0 for t in edge_times):
                    period = math.nan
                    frequency = math.nan
                    St = math.nan
                else:
                    period = float(sum(edge_times))
                    frequency = 1.0 / period
                    St = frequency * float(diameter) / float(reference_velocity)
                edge_support = [int(counts[j, k]) for j, k in edges]
                edge_prob = [float(Q[j, k]) for j, k in edges]
                cycles.append({
                    "nodes": list(canonical),
                    "edges": [[int(j), int(k)] for j, k in edges],
                    "period": period,
                    "frequency": frequency,
                    "St": St,
                    "min_support": int(min(edge_support)),
                    "total_support": int(sum(edge_support)),
                    "edge_support": edge_support,
                    "strength": float(np.prod(edge_prob)) if edge_prob else 0.0,
                })
            elif nxt not in path and nxt >= start:
                dfs(start, nxt, path + [nxt])

    for start in range(K):
        dfs(start, start, [start])

    cycles.sort(
        key=lambda c: (
            -c["min_support"],
            -c["strength"],
            c["period"] if np.isfinite(c["period"]) else np.inf,
        )
    )
    return cycles
