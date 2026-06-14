#!/usr/bin/env python3
"""
Standalone POD + k-means++ + CMM/CNM analysis.

Works directly with OpenFOAM forceCoeffs files, e.g.:

python3 python/scripts/run_cnm_analysis.py \
    --input-type force \
    --input openfoam/cases/cylinderLES3D/postProcessing/forceCoeffs/0/forceCoeffs.dat \
    --out data/figures/cnm_force_t20_from_t2 \
    --K 5 \
    --t-min 2

Also supports generic .npz files with:
    time: shape (n_samples,)
    X:    shape (n_samples, n_features)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.manifold import MDS
from sklearn.preprocessing import StandardScaler


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_force_coeffs(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Could not find forceCoeffs file: {path}")

    header = None
    data_rows = []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith("#"):
                candidate = stripped.lstrip("#").strip().split()
                if candidate and candidate[0].lower() == "time":
                    header = candidate
                continue

            parts = stripped.split()
            try:
                data_rows.append([float(x) for x in parts])
            except ValueError:
                continue

    if not data_rows:
        raise ValueError(f"No numeric data found in {path}")

    n_cols = len(data_rows[0])

    if header is None or len(header) != n_cols:
        fallback = ["Time", "Cm", "Cd", "Cl", "Cl(f)", "Cl(r)"]
        if n_cols <= len(fallback):
            header = fallback[:n_cols]
        else:
            header = ["Time"] + [f"col_{i}" for i in range(1, n_cols)]

    df = pd.DataFrame(data_rows, columns=header)

    if "Time" not in df.columns:
        raise ValueError(f"Expected a Time column. Found: {list(df.columns)}")

    df = df.drop_duplicates(subset=["Time"], keep="last")
    df = df.sort_values("Time").reset_index(drop=True)

    return df


def force_coeffs_to_matrix(
    df: pd.DataFrame,
    columns: list[str],
    t_min: float | None,
    t_max: float | None,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    mask = np.ones(len(df), dtype=bool)

    if t_min is not None:
        mask &= df["Time"].to_numpy() >= t_min

    if t_max is not None:
        mask &= df["Time"].to_numpy() <= t_max

    df = df.loc[mask].copy()

    used_columns = [c for c in columns if c in df.columns]

    if not used_columns:
        raise ValueError(
            f"None of requested force columns found: {columns}. "
            f"Available columns: {list(df.columns)}"
        )

    time = df["Time"].to_numpy(dtype=float)
    X = df[used_columns].to_numpy(dtype=float)

    valid = np.isfinite(time) & np.all(np.isfinite(X), axis=1)
    time = time[valid]
    X = X[valid]

    return time, X, used_columns


def read_npz_matrix(path: Path) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(path)

    data = np.load(path)

    if "X" not in data:
        raise ValueError(f"{path} must contain key 'X'")

    X = np.asarray(data["X"], dtype=float)

    if "time" in data:
        time = np.asarray(data["time"], dtype=float)
    else:
        time = np.arange(X.shape[0], dtype=float)

    if X.ndim != 2:
        raise ValueError(f"Expected X to be 2D, got {X.shape}")

    if time.ndim != 1 or len(time) != X.shape[0]:
        raise ValueError(f"time must have length {X.shape[0]}, got {time.shape}")

    return time, X


# ---------------------------------------------------------------------------
# POD / SVD
# ---------------------------------------------------------------------------

def compute_pod(
    X: np.ndarray,
    energy_threshold: float = 0.9999,
    n_modes_max: int | None = None,
    standardize: bool = False,
) -> dict:
    X = np.asarray(X, dtype=float)

    if X.ndim != 2:
        raise ValueError(f"X must be 2D, got {X.shape}")

    scaler = None

    if standardize:
        scaler = StandardScaler()
        X_work = scaler.fit_transform(X)
    else:
        X_work = X.copy()

    mean = X_work.mean(axis=0)
    X_fluc = X_work - mean

    U, S, Vt = np.linalg.svd(X_fluc, full_matrices=False)

    if np.allclose(S, 0):
        raise ValueError("All singular values are zero; no fluctuation energy found.")

    energy = S**2 / np.sum(S**2)
    cumulative_energy = np.cumsum(energy)

    n_modes = int(np.searchsorted(cumulative_energy, energy_threshold) + 1)

    if n_modes_max is not None:
        n_modes = min(n_modes, int(n_modes_max))

    n_modes = max(1, min(n_modes, len(S)))

    coefficients = U[:, :n_modes] * S[:n_modes]
    modes = Vt[:n_modes, :]

    return {
        "mean": mean,
        "modes": modes,
        "coefficients": coefficients,
        "singular_values": S,
        "energy": energy,
        "cumulative_energy": cumulative_energy,
        "n_modes": n_modes,
        "scaler": scaler,
    }


# ---------------------------------------------------------------------------
# Clustering / CNM
# ---------------------------------------------------------------------------

def reorder_by_population_and_transitions(
    labels: np.ndarray,
    centers: np.ndarray,
    K: int,
) -> tuple[np.ndarray, np.ndarray]:
    labels = labels.astype(int)
    counts = np.bincount(labels, minlength=K)

    trans_counts = np.zeros((K, K), dtype=int)

    for old, new in zip(labels[:-1], labels[1:]):
        trans_counts[new, old] += 1

    first = int(np.argmax(counts))

    order = [first]
    unused = set(range(K))
    unused.remove(first)

    current = first

    while unused:
        candidates = sorted(unused)
        weights = np.array([trans_counts[c, current] for c in candidates])

        if weights.max() > 0:
            nxt = candidates[int(np.argmax(weights))]
        else:
            nxt = max(candidates, key=lambda c: counts[c])

        order.append(nxt)
        unused.remove(nxt)
        current = nxt

    old_to_new = {old: new for new, old in enumerate(order)}
    new_labels = np.array([old_to_new[x] for x in labels], dtype=int)
    new_centers = centers[order, :]

    return new_labels, new_centers


def cluster_kmeans_pp(
    A: np.ndarray,
    K: int,
    random_state: int,
    n_init: int = 50,
    reorder: bool = True,
) -> dict:
    if K < 2:
        raise ValueError("K must be at least 2.")

    if K > A.shape[0]:
        raise ValueError(f"K={K} larger than number of samples {A.shape[0]}.")

    km = KMeans(
        n_clusters=K,
        init="k-means++",
        n_init=n_init,
        random_state=random_state,
        algorithm="lloyd",
    )

    labels = km.fit_predict(A)
    centers = km.cluster_centers_

    if reorder:
        labels, centers = reorder_by_population_and_transitions(labels, centers, K)

    populations = np.bincount(labels, minlength=K) / len(labels)

    return {
        "labels": labels,
        "centers": centers,
        "populations": populations,
    }


def transition_matrix_P(labels: np.ndarray, K: int) -> np.ndarray:
    P_counts = np.zeros((K, K), dtype=float)

    for k_old, j_new in zip(labels[:-1], labels[1:]):
        P_counts[j_new, k_old] += 1.0

    col_sum = P_counts.sum(axis=0, keepdims=True)

    P = np.divide(
        P_counts,
        col_sum,
        out=np.zeros_like(P_counts),
        where=col_sum > 0,
    )

    return P


def run_length_encoding(labels: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    labels = labels.astype(int)

    if len(labels) == 0:
        raise ValueError("labels must not be empty")

    values = [labels[0]]
    starts = [0]

    for i in range(1, len(labels)):
        if labels[i] != labels[i - 1]:
            values.append(labels[i])
            starts.append(i)

    values = np.asarray(values, dtype=int)
    starts = np.asarray(starts, dtype=int)
    ends = np.r_[starts[1:] - 1, len(labels) - 1].astype(int)

    return values, starts, ends


def cnm_matrices(labels: np.ndarray, time: np.ndarray, K: int) -> dict:
    labels = labels.astype(int)
    time = np.asarray(time, dtype=float)

    if len(labels) != len(time):
        raise ValueError(f"labels/time mismatch: {len(labels)} vs {len(time)}")

    P = transition_matrix_P(labels, K)

    run_labels, run_start, run_end = run_length_encoding(labels)

    dt_default = float(np.median(np.diff(time))) if len(time) >= 2 else 1.0

    residence_by_run = []
    residence_times = {k: [] for k in range(K)}

    for lab, i0, i1 in zip(run_labels, run_start, run_end):
        duration = float(time[i1] - time[i0] + dt_default)
        duration = max(duration, dt_default)

        residence_by_run.append(duration)
        residence_times[int(lab)].append(duration)

    Q_counts = np.zeros((K, K), dtype=float)
    T_sums = np.zeros((K, K), dtype=float)

    for r in range(len(run_labels) - 1):
        k = int(run_labels[r])
        j = int(run_labels[r + 1])

        if j == k:
            continue

        individual_T = 0.5 * residence_by_run[r] + 0.5 * residence_by_run[r + 1]

        Q_counts[j, k] += 1.0
        T_sums[j, k] += individual_T

    col_sum = Q_counts.sum(axis=0, keepdims=True)

    Q = np.divide(
        Q_counts,
        col_sum,
        out=np.zeros_like(Q_counts),
        where=col_sum > 0,
    )

    T = np.divide(
        T_sums,
        Q_counts,
        out=np.full_like(T_sums, np.nan),
        where=Q_counts > 0,
    )

    return {
        "P": P,
        "Q": Q,
        "T": T,
        "residence_times": residence_times,
        "run_labels": run_labels,
        "run_start_idx": run_start,
        "run_end_idx": run_end,
    }


def deterministic_cnm_prediction(
    time: np.ndarray,
    centers: np.ndarray,
    Q: np.ndarray,
    T: np.ndarray,
    start_cluster: int,
) -> tuple[np.ndarray, np.ndarray]:
    time = np.asarray(time, dtype=float)
    centers = np.asarray(centers, dtype=float)

    K, n_features = centers.shape

    pred_labels = np.zeros(len(time), dtype=int)
    pred_A = np.zeros((len(time), n_features), dtype=float)

    if len(time) == 0:
        return pred_labels, pred_A

    dt_default = float(np.median(np.diff(time))) if len(time) >= 2 else 1.0

    positive_T = T[np.isfinite(T) & (T > 0)]
    fallback_T = float(np.median(positive_T)) if positive_T.size else dt_default

    current = int(start_cluster)
    t_depart = float(time[0])
    idx = 0

    while idx < len(time):
        probs = Q[:, current].copy()
        probs[current] = 0.0

        if probs.sum() <= 0:
            pred_labels[idx:] = current
            pred_A[idx:, :] = centers[current]
            break

        nxt = int(np.argmax(probs))

        duration = T[nxt, current]

        if not np.isfinite(duration) or duration <= 0:
            duration = fallback_T

        t_arrive = t_depart + float(duration)

        while idx < len(time) and time[idx] <= t_arrive:
            alpha = (time[idx] - t_depart) / max(float(duration), 1e-12)
            alpha = float(np.clip(alpha, 0.0, 1.0))

            pred_labels[idx] = current if alpha < 0.5 else nxt
            pred_A[idx, :] = (1.0 - alpha) * centers[current] + alpha * centers[nxt]

            idx += 1

        current = nxt
        t_depart = t_arrive

    return pred_labels, pred_A


def mds_embedding(A: np.ndarray, centers: np.ndarray, random_state: int) -> tuple[np.ndarray, np.ndarray]:
    combined = np.vstack([A, centers])

    n_components = min(3, combined.shape[0] - 1, combined.shape[1])
    n_components = max(2, n_components)

    mds = MDS(
        n_components=n_components,
        random_state=random_state,
        metric_mds=True,
        n_init=4,
        init="random",
        max_iter=300,
    )

    emb = mds.fit_transform(combined)

    emb_A = emb[: len(A)]
    emb_C = emb[len(A):]

    if emb_A.shape[1] == 2:
        emb_A = np.c_[emb_A, np.zeros(len(emb_A))]
        emb_C = np.c_[emb_C, np.zeros(len(emb_C))]

    return emb_A, emb_C


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def save_energy_plot(energy: np.ndarray, cumulative: np.ndarray, n_modes: int, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)

    x = np.arange(1, len(energy) + 1)

    plt.figure(figsize=(8, 4.5))
    plt.semilogy(x, energy, "o-", label="Individual energy")
    plt.plot(x, cumulative, "s-", label="Cumulative energy")
    plt.axvline(n_modes, linestyle="--", linewidth=1.0, label=f"Kept modes = {n_modes}")
    plt.xlabel("Mode index")
    plt.ylabel("Energy fraction")
    plt.title("POD/SVD energy spectrum")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_matrix_plot(M: np.ndarray, title: str, out: Path, value_label: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7, 6))
    im = plt.imshow(M, origin="lower", aspect="equal")
    plt.colorbar(im, label=value_label)
    plt.xlabel("From cluster k")
    plt.ylabel("To cluster j")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_cluster_affiliation_plot(
    time: np.ndarray,
    labels: np.ndarray,
    pred_labels: np.ndarray | None,
    out: Path,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(11, 4))
    plt.step(time, labels, where="post", linewidth=1.2, label="Data")

    if pred_labels is not None:
        plt.step(time, pred_labels, where="post", linewidth=1.2, label="CNM deterministic")

    plt.xlabel("Time")
    plt.ylabel("Cluster")
    plt.title("Temporal evolution of cluster affiliation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_coefficients_plot(
    time: np.ndarray,
    A: np.ndarray,
    pred_A: np.ndarray | None,
    out: Path,
    n: int = 4,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)

    n = min(n, A.shape[1])

    plt.figure(figsize=(12, 2.7 * n))

    for i in range(n):
        ax = plt.subplot(n, 1, i + 1)
        ax.plot(time, A[:, i], linewidth=1.1, label=f"Data a{i + 1}")

        if pred_A is not None:
            ax.plot(time, pred_A[:, i], linewidth=1.1, linestyle="--", label=f"CNM a{i + 1}")

        ax.set_ylabel(f"a{i + 1}")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best")

    plt.xlabel("Time")
    plt.suptitle("Temporal evolution of leading coefficients", y=0.995)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_population_plot(
    labels: np.ndarray,
    pred_labels: np.ndarray | None,
    K: int,
    out: Path,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)

    data_pop = np.bincount(labels, minlength=K) / len(labels)
    x = np.arange(K)

    plt.figure(figsize=(9, 4.5))

    if pred_labels is not None:
        pred_pop = np.bincount(pred_labels, minlength=K) / len(pred_labels)
        width = 0.38
        plt.bar(x - width / 2, data_pop, width=width, label="Data")
        plt.bar(x + width / 2, pred_pop, width=width, label="CNM")
    else:
        plt.bar(x, data_pop, width=0.6, label="Data")

    plt.xlabel("Cluster")
    plt.ylabel("Probability")
    plt.title("Cluster probability distribution")
    plt.grid(True, axis="y", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_mds_plot(
    emb_A: np.ndarray,
    emb_C: np.ndarray,
    labels: np.ndarray,
    out: Path,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(
        emb_A[:, 0],
        emb_A[:, 1],
        emb_A[:, 2],
        c=labels,
        s=8,
        alpha=0.75,
    )

    ax.scatter(
        emb_C[:, 0],
        emb_C[:, 1],
        emb_C[:, 2],
        s=120,
        marker="X",
        edgecolors="black",
        linewidths=0.8,
    )

    for k, c in enumerate(emb_C):
        ax.text(c[0], c[1], c[2], str(k), fontsize=9)

    fig.colorbar(sc, ax=ax, label="Cluster")

    ax.set_title("MDS cluster map")
    ax.set_xlabel("x1")
    ax.set_ylabel("x2")
    ax.set_zlabel("x3")

    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_force_coefficients_plot(
    time: np.ndarray,
    X: np.ndarray,
    feature_names: list[str],
    out: Path,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 5))

    for i, name in enumerate(feature_names):
        plt.plot(time, X[:, i], linewidth=1.2, label=name)

    plt.xlabel("Time")
    plt.ylabel("Coefficient")
    plt.title("OpenFOAM force coefficients")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Robust y-limits to avoid initial startup spikes destroying the figure.
    values = X.ravel()
    q_low = pd.Series(values).quantile(0.01)
    q_high = pd.Series(values).quantile(0.99)
    margin = 0.15 * max(abs(q_high - q_low), 1e-12)
    plt.ylim(q_low - margin, q_high + margin)

    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="POD + k-means++ + CMM/CNM analysis")

    parser.add_argument("--input-type", choices=["force", "npz"], required=True)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)

    parser.add_argument("--K", type=int, default=10)
    parser.add_argument("--energy-threshold", type=float, default=0.9999)
    parser.add_argument("--n-modes-max", type=int, default=None)
    parser.add_argument("--random-state", type=int, default=1)

    parser.add_argument("--t-min", type=float, default=None)
    parser.add_argument("--t-max", type=float, default=None)

    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Use only every stride-th sample after time filtering.",
    )

    parser.add_argument(
        "--force-columns",
        nargs="+",
        default=["Cd", "Cl", "Cm"],
        help="forceCoeffs columns to use as features.",
    )

    parser.add_argument(
        "--standardize",
        action="store_true",
        help="Standardize features before POD for npz input.",
    )

    parser.add_argument(
        "--no-cnm-prediction",
        action="store_true",
        help="Skip deterministic CNM reconstruction.",
    )

    return parser.parse_args()


def load_input(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, list[str], bool]:
    if args.input_type == "force":
        df = read_force_coeffs(args.input)

        time, X, feature_names = force_coeffs_to_matrix(
            df=df,
            columns=args.force_columns,
            t_min=args.t_min,
            t_max=args.t_max,
        )

        # Force coefficients have different scales, so standardization is useful.
        standardize = True

        return time, X, feature_names, standardize

    if args.input_type == "npz":
        time, X = read_npz_matrix(args.input)

        if args.t_min is not None or args.t_max is not None:
            mask = np.ones(len(time), dtype=bool)

            if args.t_min is not None:
                mask &= time >= args.t_min

            if args.t_max is not None:
                mask &= time <= args.t_max

            time = time[mask]
            X = X[mask]

        feature_names = [f"x{i + 1}" for i in range(X.shape[1])]

        standardize = args.standardize

        return time, X, feature_names, standardize

    raise ValueError(f"Unsupported input type: {args.input_type}")


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    time, X, feature_names, standardize = load_input(args)

    if args.stride < 1:
        raise ValueError("--stride must be >= 1")

    if args.stride > 1:
        time = time[::args.stride]
        X = X[::args.stride]

    if len(time) < args.K:
        raise ValueError(f"Need at least K={args.K} samples, got {len(time)}.")

    print("Loaded input:")
    print(f"  input type : {args.input_type}")
    print(f"  input file : {args.input}")
    print(f"  samples    : {X.shape[0]}")
    print(f"  features   : {X.shape[1]} -> {feature_names}")
    print(f"  time       : {time.min()} -> {time.max()}")
    print(f"  standardize: {standardize}")

    if args.input_type == "force":
        save_force_coefficients_plot(
            time=time,
            X=X,
            feature_names=feature_names,
            out=args.out / "00_force_coefficients.png",
        )

    pod = compute_pod(
        X=X,
        energy_threshold=args.energy_threshold,
        n_modes_max=args.n_modes_max,
        standardize=standardize,
    )

    A = pod["coefficients"]

    print("")
    print("POD/SVD:")
    print(f"  modes kept : {pod['n_modes']}")
    print(f"  energy kept: {pod['cumulative_energy'][pod['n_modes'] - 1]:.8f}")

    clusters = cluster_kmeans_pp(
        A=A,
        K=args.K,
        random_state=args.random_state,
        n_init=50,
        reorder=True,
    )

    labels = clusters["labels"]
    centers = clusters["centers"]
    populations = clusters["populations"]

    print("")
    print("Clustering:")
    print(f"  K: {args.K}")
    print("  populations:")

    for k, p in enumerate(populations):
        print(f"    C{k}: {p:.5f}")

    cnm = cnm_matrices(labels=labels, time=time, K=args.K)

    if args.no_cnm_prediction:
        pred_labels = None
        pred_A = None
    else:
        pred_labels, pred_A = deterministic_cnm_prediction(
            time=time,
            centers=centers,
            Q=cnm["Q"],
            T=cnm["T"],
            start_cluster=int(labels[0]),
        )

    emb_A, emb_C = mds_embedding(
        A=A,
        centers=centers,
        random_state=args.random_state,
    )

    # Save numeric outputs
    np.save(args.out / "time.npy", time)
    np.save(args.out / "X.npy", X)
    np.save(args.out / "pod_coefficients.npy", A)
    np.save(args.out / "cluster_labels.npy", labels)
    np.save(args.out / "cluster_centers.npy", centers)
    np.save(args.out / "P_matrix.npy", cnm["P"])
    np.save(args.out / "Q_matrix.npy", cnm["Q"])
    np.save(args.out / "T_matrix.npy", cnm["T"])

    if pred_labels is not None:
        np.save(args.out / "cnm_pred_labels.npy", pred_labels)

    if pred_A is not None:
        np.save(args.out / "cnm_pred_coefficients.npy", pred_A)

    # Save figures
    save_energy_plot(
        energy=pod["energy"],
        cumulative=pod["cumulative_energy"],
        n_modes=pod["n_modes"],
        out=args.out / "01_pod_energy.png",
    )

    save_mds_plot(
        emb_A=emb_A,
        emb_C=emb_C,
        labels=labels,
        out=args.out / "02_mds_clusters.png",
    )

    save_matrix_plot(
        M=cnm["P"],
        title="CMM transition probability matrix P",
        out=args.out / "03_P_matrix.png",
        value_label="probability",
    )

    save_matrix_plot(
        M=cnm["Q"],
        title="CNM direct transition probability matrix Q",
        out=args.out / "04_Q_matrix.png",
        value_label="probability",
    )

    save_matrix_plot(
        M=cnm["T"],
        title="CNM transition time matrix T",
        out=args.out / "05_T_matrix.png",
        value_label="time",
    )

    save_cluster_affiliation_plot(
        time=time,
        labels=labels,
        pred_labels=pred_labels,
        out=args.out / "06_cluster_affiliation.png",
    )

    save_coefficients_plot(
        time=time,
        A=A,
        pred_A=pred_A,
        out=args.out / "07_first_coefficients.png",
        n=4,
    )

    save_population_plot(
        labels=labels,
        pred_labels=pred_labels,
        K=args.K,
        out=args.out / "08_cluster_probabilities.png",
    )

    summary = {
        "input": str(args.input),
        "input_type": args.input_type,
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "features": feature_names,
        "time_min": float(time.min()),
        "time_max": float(time.max()),
        "stride": int(args.stride),
        "K": int(args.K),
        "pod_n_modes": int(pod["n_modes"]),
        "pod_energy_threshold": float(args.energy_threshold),
        "pod_cumulative_energy_kept": float(
            pod["cumulative_energy"][pod["n_modes"] - 1]
        ),
        "cluster_populations": populations.tolist(),
        "outputs": {
            "force_coefficients": "00_force_coefficients.png",
            "energy": "01_pod_energy.png",
            "mds": "02_mds_clusters.png",
            "P": "03_P_matrix.png",
            "Q": "04_Q_matrix.png",
            "T": "05_T_matrix.png",
            "affiliation": "06_cluster_affiliation.png",
            "coefficients": "07_first_coefficients.png",
            "probabilities": "08_cluster_probabilities.png",
        },
    }

    with (args.out / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("")
    print("CNM analysis complete.")
    print(f"  output: {args.out}")
    print("")
    print("Key plots:")

    for name in summary["outputs"].values():
        print(f"  - {args.out / name}")


if __name__ == "__main__":
    main()
