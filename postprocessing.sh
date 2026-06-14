#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# postprocessing.sh – Post-processing for cfd_rom_regime_clustering_cylinder
#
# Usage examples:
#   ./postprocessing.sh --K 8 --stride 10 --t-min 20 --tag t60
#   ./postprocessing.sh --K 8 --stride 10 --t-min 20 --tag t60 --no-vtk
#   ./postprocessing.sh --K 10 --stride 20 --t-min 30 --tag t140 --vtk-latest
#
# Steps:
#   1. choose correct project-root venv
#   2. merge all forceCoeffs.dat files from continuation runs
#   3. plot force coefficients
#   4. run CNM analysis
#   5. optionally reconstruct and export VTK
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}"
OF_CASE="${PROJECT_ROOT}/openfoam/cases/cylinderLES3D"
PYTHON_DIR="${PROJECT_ROOT}/python"
OUT_BASE="${PROJECT_ROOT}/data/figures"
PROCESSED_DIR="${PROJECT_ROOT}/data/processed"

K=10
STRIDE=10
T_MIN=""
T_MAX=""
TAG=""
DO_RECONSTRUCT=1
DO_VTK=1
VTK_LATEST=0
VTK_TIME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --K)             K="$2"; shift 2 ;;
    --stride)        STRIDE="$2"; shift 2 ;;
    --t-min)         T_MIN="$2"; shift 2 ;;
    --t-max)         T_MAX="$2"; shift 2 ;;
    --tag)           TAG="$2"; shift 2 ;;
    --no-reconstruct) DO_RECONSTRUCT=0; shift ;;
    --no-vtk)        DO_VTK=0; shift ;;
    --vtk-latest)    VTK_LATEST=1; shift ;;
    --vtk-time)      VTK_TIME="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

mkdir -p "${OUT_BASE}" "${PROCESSED_DIR}"

LABEL="$(date '+%Y%m%d_%H%M%S')_K${K}_stride${STRIDE}${TAG:+_$TAG}"
OUT="${OUT_BASE}/cnm_${LABEL}"
mkdir -p "${OUT}"

# -----------------------------------------------------------------------------
# Python / venv
# -----------------------------------------------------------------------------

if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
  PYTHON="${PROJECT_ROOT}/.venv/bin/python"
elif [[ -x "${PROJECT_ROOT}/.venv/bin/python3" ]]; then
  PYTHON="${PROJECT_ROOT}/.venv/bin/python3"
elif [[ -x "${PYTHON_DIR}/.venv/bin/python" ]]; then
  PYTHON="${PYTHON_DIR}/.venv/bin/python"
elif [[ -x "${PYTHON_DIR}/.venv/bin/python3" ]]; then
  PYTHON="${PYTHON_DIR}/.venv/bin/python3"
else
  PYTHON="$(command -v python3)"
fi

echo "==> Project root : ${PROJECT_ROOT}"
echo "==> OpenFOAM case: ${OF_CASE}"
echo "==> Python       : ${PYTHON}"
echo "==> Output       : ${OUT}"
echo ""

"${PYTHON}" - <<'PY'
import sys
missing = []
for mod in ["numpy", "pandas", "matplotlib", "sklearn", "scipy"]:
    try:
        __import__(mod)
    except Exception:
        missing.append(mod)
if missing:
    print("ERROR: Missing Python modules:", ", ".join(missing))
    print("Fix with:")
    print("  python3 -m venv .venv")
    print("  source .venv/bin/activate")
    print("  pip install -r python/requirements.txt")
    sys.exit(1)
print("Python environment OK")
PY

# -----------------------------------------------------------------------------
# Check case
# -----------------------------------------------------------------------------

if [[ ! -d "${OF_CASE}" ]]; then
  echo "ERROR: OpenFOAM case not found: ${OF_CASE}"
  exit 1
fi

if [[ ! -f "${PYTHON_DIR}/scripts/plot_force_coeffs.py" ]]; then
  echo "ERROR: Missing ${PYTHON_DIR}/scripts/plot_force_coeffs.py"
  exit 1
fi

if [[ ! -f "${PYTHON_DIR}/scripts/run_cnm_analysis.py" ]]; then
  echo "ERROR: Missing ${PYTHON_DIR}/scripts/run_cnm_analysis.py"
  exit 1
fi

# -----------------------------------------------------------------------------
# 1. reconstructPar
# -----------------------------------------------------------------------------

if [[ "${DO_RECONSTRUCT}" == "1" ]]; then
  echo "==> reconstructPar"
  (
    cd "${OF_CASE}"
    mkdir -p logs
    reconstructPar 2>&1 | tee "${OUT}/log.reconstructPar"
  )
else
  echo "==> Skipping reconstructPar"
fi

echo ""

# -----------------------------------------------------------------------------
# 2. Optional VTK export
# -----------------------------------------------------------------------------

if [[ "${DO_VTK}" == "1" ]]; then
  echo "==> foamToVTK"

  VTK_ARGS=(-fields "(U p)" -noPointValues)

  if [[ "${VTK_LATEST}" == "1" ]]; then
    VTK_ARGS=(-latestTime "${VTK_ARGS[@]}")
  fi

  if [[ -n "${VTK_TIME}" ]]; then
    VTK_ARGS=(-time "${VTK_TIME}" "${VTK_ARGS[@]}")
  fi

  (
    cd "${OF_CASE}"
    foamToVTK "${VTK_ARGS[@]}" 2>&1 | tee "${OUT}/log.foamToVTK"
  )
else
  echo "==> Skipping foamToVTK"
fi

echo ""

# -----------------------------------------------------------------------------
# 3. Find and merge forceCoeffs files
# -----------------------------------------------------------------------------

echo "==> Searching forceCoeffs files"

mapfile -t FORCE_FILES < <(
	  find "${OF_CASE}/postProcessing/forceCoeffs" \
	    \( -name "forceCoeffs.dat" -o -name "coefficient.dat" \) \
	    -type f \
	    ! -path "*/combined/*" \
	    ! -path "*old*" \
	    2>/dev/null | sort -V
	)

if [[ "${#FORCE_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: no forceCoeffs files found under ${OF_CASE}/postProcessing/forceCoeffs"
  exit 1
fi

printf '    %s\n' "${FORCE_FILES[@]}"

MERGED_FORCE="${PROCESSED_DIR}/forceCoeffs_merged_${LABEL}.dat"

echo ""
echo "==> Merging forceCoeffs:"
echo "    ${MERGED_FORCE}"

"${PYTHON}" - "${MERGED_FORCE}" "${FORCE_FILES[@]}" <<'PY'
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
    f.write("# Force coefficients merged by postprocessing.sh\n")
    f.write("# " + " ".join(df.columns) + "\n")
    for _, row in df.iterrows():
        f.write(" ".join(f"{float(v):.12g}" for v in row.values) + "\n")

print(f"rows      : {len(df)}")
print(f"time range: {df['Time'].min()} -> {df['Time'].max()}")

for col in ["Cd", "Cl", "Cm"]:
    if col in df.columns:
        print(f"{col}: min={df[col].min():.6g}, max={df[col].max():.6g}, last={df[col].iloc[-1]:.6g}")
PY

echo ""

# -----------------------------------------------------------------------------
# 4. Auto T_MIN if not given
# -----------------------------------------------------------------------------

if [[ -z "${T_MIN}" ]]; then
  T_MIN=$("${PYTHON}" - "${MERGED_FORCE}" <<'PY'
import sys
import pandas as pd

path = sys.argv[1]
df = pd.read_csv(path, comment="#", sep=r"\s+", header=None)

t_end = float(df.iloc[:, 0].max())
print(0.2 * t_end)
PY
)
  echo "==> t_min auto: ${T_MIN}"
fi

# -----------------------------------------------------------------------------
# 5. Force coefficient plot
# -----------------------------------------------------------------------------

echo "==> plot_force_coeffs"

PLOT_ARGS=(
  "${PYTHON_DIR}/scripts/plot_force_coeffs.py"
  --file "${MERGED_FORCE}"
  --out "${OUT}/00_force_coefficients.png"
  --t-min "${T_MIN}"
)

if [[ -n "${T_MAX}" ]]; then
  PLOT_ARGS+=(--t-max "${T_MAX}")
fi

"${PYTHON}" "${PLOT_ARGS[@]}"

echo ""

# -----------------------------------------------------------------------------
# 6. CNM analysis
# -----------------------------------------------------------------------------

echo "==> run_cnm_analysis"
echo "    K=${K}, t_min=${T_MIN}, stride=${STRIDE}"

CNM_ARGS=(
  "${PYTHON_DIR}/scripts/run_cnm_analysis.py"
  --input-type force
  --input "${MERGED_FORCE}"
  --out "${OUT}"
  --K "${K}"
  --energy-threshold 0.9999
  --t-min "${T_MIN}"
  --stride "${STRIDE}"
  --force-columns Cd Cl Cm
)

if [[ -n "${T_MAX}" ]]; then
  CNM_ARGS+=(--t-max "${T_MAX}")
fi

"${PYTHON}" "${CNM_ARGS[@]}"

echo ""
echo "==> done"
echo "Output:"
echo "  ${OUT}"
echo ""
echo "Main files:"
echo "  ${OUT}/00_force_coefficients.png"
echo "  ${OUT}/03_P_matrix.png"
echo "  ${OUT}/04_Q_matrix.png"
echo "  ${OUT}/05_T_matrix.png"
echo "  ${OUT}/06_cluster_affiliation.png"
echo ""
echo "Open with:"
echo "  xdg-open \"${OUT}/00_force_coefficients.png\""
