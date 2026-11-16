#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# sweep_cnm_windows.sh
#
# Compare several time windows of forceCoeffs with POD + k-means++ + CNM.
#
# Usage:
#   ./scripts/sweep_cnm_windows.sh
#
# Optional:
#   K_LIST="6 8 10 12" STRIDE=10 RESAMPLE_DT=0.01 ./scripts/sweep_cnm_windows.sh
# =============================================================================

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo "Usage: K_LIST='6 8 10' STRIDE=10 RESAMPLE_DT=0.01 ./scripts/sweep_cnm_windows.sh"
    echo "RESAMPLE_DT is optional: interpolate force times to this fixed interval before stride."
    echo "Without RESAMPLE_DT, force data must already be uniformly sampled."
    exit 0
fi

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

CASE_DIR="$PROJECT_ROOT/openfoam/cases/cylinderLES3D"
PROCESSED_DIR="$PROJECT_ROOT/data/processed"
FIG_DIR="$PROJECT_ROOT/data/figures/window_sweep_$(date '+%Y%m%d_%H%M%S')"

K_LIST="${K_LIST:-6 8 10}"
STRIDE="${STRIDE:-10}"
RESAMPLE_DT="${RESAMPLE_DT:-}"

mkdir -p "$PROCESSED_DIR" "$FIG_DIR"

# -----------------------------------------------------------------------------
# Python selection
# -----------------------------------------------------------------------------

if [[ -x "$PROJECT_ROOT/.venv/bin/python" ]]; then
    PYTHON="$PROJECT_ROOT/.venv/bin/python"
elif [[ -x "$PROJECT_ROOT/.venv/bin/python3" ]]; then
    PYTHON="$PROJECT_ROOT/.venv/bin/python3"
else
    PYTHON="$(command -v python3)"
fi

echo "==> Project root : $PROJECT_ROOT"
echo "==> Case dir     : $CASE_DIR"
echo "==> Python       : $PYTHON"
echo "==> Figure dir   : $FIG_DIR"
echo "==> K list       : $K_LIST"
echo "==> Stride       : $STRIDE"
echo "==> Resample dt  : ${RESAMPLE_DT:-disabled; input must be uniform}"
echo ""

"$PYTHON" - "$RESAMPLE_DT" <<'PY'
import math
import sys
if sys.argv[1]:
    try:
        dt = float(sys.argv[1])
        if not math.isfinite(dt) or dt <= 0:
            raise ValueError
    except ValueError:
        raise SystemExit("ERROR: RESAMPLE_DT must be a finite positive number")
missing = []
for mod in ["numpy", "pandas", "matplotlib", "sklearn", "scipy"]:
    try:
        __import__(mod)
    except Exception:
        missing.append(mod)

if missing:
    print("ERROR: Missing Python modules:", ", ".join(missing))
    print("Run:")
    print("  source .venv/bin/activate")
    print("  pip install -r python/requirements.txt")
    sys.exit(1)

print("Python environment OK")
PY

echo ""

# -----------------------------------------------------------------------------
# Merge forceCoeffs files
# -----------------------------------------------------------------------------

echo "==> Searching forceCoeffs files"

mapfile -t FORCE_FILES < <(
    find "$CASE_DIR/postProcessing/forceCoeffs" \
        \( -name "forceCoeffs.dat" -o -name "coefficient.dat" \) \
        -type f \
        ! -path "*/combined/*" \
        ! -path "*old*" \
        2>/dev/null | sort -V
)

if [[ "${#FORCE_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: no forceCoeffs files found."
    exit 1
fi

printf '    %s\n' "${FORCE_FILES[@]}"
echo ""

MERGED_FORCE="$PROCESSED_DIR/forceCoeffs_merged_window_sweep_$(date '+%Y%m%d_%H%M%S').dat"

echo "==> Merging into:"
echo "    $MERGED_FORCE"

"$PYTHON" - "$MERGED_FORCE" "${FORCE_FILES[@]}" <<'PY'
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

out = Path(sys.argv[1])
files = [Path(p) for p in sys.argv[2:]]

header = None
rows = []

for path in files:
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
    raise SystemExit("No numeric forceCoeffs data found.")

n_cols = len(rows[0])

if header is None or len(header) != n_cols:
    fallback = ["Time", "Cm", "Cd", "Cl", "Cl(f)", "Cl(r)"]
    if n_cols <= len(fallback):
        header = fallback[:n_cols]
    else:
        header = ["Time"] + [f"col_{i}" for i in range(1, n_cols)]

df = pd.DataFrame(rows, columns=header)

if "Time" not in df.columns:
    raise SystemExit(f"No Time column found. Columns: {list(df.columns)}")

df = df.drop_duplicates(subset=["Time"], keep="last")
df = df.sort_values("Time").reset_index(drop=True)

out.parent.mkdir(parents=True, exist_ok=True)

with out.open("w", encoding="utf-8") as f:
    f.write("# Force coefficients merged by scripts/sweep_cnm_windows.sh\n")
    f.write("# " + " ".join(df.columns) + "\n")

    for _, row in df.iterrows():
        f.write(" ".join(f"{float(v):.12g}" for v in row.values) + "\n")

print(f"rows      : {len(df)}")
print(f"time range: {df['Time'].min()} -> {df['Time'].max()}")

for col in ["Cd", "Cl", "Cm"]:
    if col in df.columns:
        print(f"{col}: min={df[col].min():.6g}, max={df[col].max():.6g}, last={df[col].iloc[-1]:.6g}")
PY

MAX_TIME="$("$PYTHON" - "$MERGED_FORCE" <<'PY'
import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], comment="#", sep=r"\s+", header=None)
print(float(df.iloc[:, 0].max()))
PY
)"

echo ""
echo "==> Max available time: $MAX_TIME"
echo ""

# -----------------------------------------------------------------------------
# Window definitions
# -----------------------------------------------------------------------------

declare -a WINDOWS=(
  "0 20 startup"
  "20 40 wake_build_up"
  "40 60 shedding_onset"
  "60 100 developed_early"
  "100 140 developed_late"
  "60 140 final_window"
)

# -----------------------------------------------------------------------------
# Run window sweep
# -----------------------------------------------------------------------------

for item in "${WINDOWS[@]}"; do
    read -r TMIN TMAX TAG <<< "$item"

    TMAX_EFF="$("$PYTHON" - "$TMAX" "$MAX_TIME" <<'PY'
import sys
target = float(sys.argv[1])
max_time = float(sys.argv[2])
print(min(target, max_time))
PY
)"

    COUNT="$("$PYTHON" - "$MERGED_FORCE" "$TMIN" "$TMAX_EFF" "$STRIDE" "$RESAMPLE_DT" <<'PY'
import sys
import pandas as pd

path = sys.argv[1]
tmin = float(sys.argv[2])
tmax = float(sys.argv[3])
stride = int(sys.argv[4])
resample_dt = float(sys.argv[5]) if sys.argv[5] else None

df = pd.read_csv(path, comment="#", sep=r"\s+", header=None)
t = df.iloc[:, 0]
mask = (t >= tmin) & (t <= tmax)
selected_time = t.loc[mask].to_numpy()
n_samples = len(selected_time)
if resample_dt is not None and n_samples >= 2:
    import math
    n_samples = int(math.floor((selected_time[-1] - selected_time[0]) / resample_dt + 1e-9)) + 1
n = (n_samples + stride - 1) // stride
print(n)
PY
)"

    if [[ "$COUNT" -lt 2 ]]; then
        echo "==> Skipping window $TAG: not enough samples after stride."
        echo "    requested: $TMIN -> $TMAX, available clipped: $TMIN -> $TMAX_EFF, samples: $COUNT"
        echo ""
        continue
    fi

    echo "==> Window: $TAG"
    echo "    requested: $TMIN -> $TMAX"
    echo "    used     : $TMIN -> $TMAX_EFF"
    echo "    samples  : $COUNT"
    echo ""

    "$PYTHON" python/scripts/plot_force_coeffs.py \
        --file "$MERGED_FORCE" \
        --out "$FIG_DIR/window_${TAG}_force_coeffs.png" \
        --t-min "$TMIN" \
        --t-max "$TMAX_EFF"

    for K in $K_LIST; do
        if [[ "$COUNT" -lt "$K" ]]; then
            echo "    Skipping K=$K because samples=$COUNT < K=$K"
            continue
        fi

        echo "    Running CNM: K=$K"

        RESAMPLE_ARGS=()
        if [[ -n "$RESAMPLE_DT" ]]; then
            RESAMPLE_ARGS=(--resample-dt "$RESAMPLE_DT")
        fi
        "$PYTHON" python/scripts/run_cnm_analysis.py \
            --input-type force \
            --input "$MERGED_FORCE" \
            --out "$FIG_DIR/cnm_window_${TAG}_K${K}_stride${STRIDE}" \
            --K "$K" \
            --t-min "$TMIN" \
            --t-max "$TMAX_EFF" \
            --stride "$STRIDE" \
            "${RESAMPLE_ARGS[@]}"
    done

    echo ""
done

echo "==> Done."
echo ""
echo "Merged force file:"
echo "  $MERGED_FORCE"
echo ""
echo "Output directory:"
echo "  $FIG_DIR"
echo ""
echo "Example open commands:"
echo "  xdg-open \"$FIG_DIR/window_shedding_onset_force_coeffs.png\""
echo "  xdg-open \"$FIG_DIR/cnm_window_final_window_K10_stride${STRIDE}/04_Q_matrix.png\""
