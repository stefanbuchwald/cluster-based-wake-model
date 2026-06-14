#!/usr/bin/env bash
# =============================================================================
#  analyse.sh  –  Post-processing pipeline for cfd_rom_regime_clustering_cylinder
#
#  Usage:
#    ./analyse.sh                    # full run, auto-detect t_min
#    ./analyse.sh --t-min 10         # skip transient manually
#    ./analyse.sh --K 12             # use 12 clusters
#    ./analyse.sh --skip-reconstruct # skip reconstructPar (data already merged)
#    ./analyse.sh --help
#
#  What this script does:
#    1.  reconstructPar  (merge parallel time directories)
#    2.  foamToVTK       (export velocity/vorticity snapshots)
#    3.  plot_force_coeffs.py   (quick Cd/Cl overview plot)
#    4.  run_cnm_analysis.py    (POD + k-means++ + CMM/CNM + all plots)
#    5.  Print summary of outputs
#
#  Requires:  OpenFOAM 13, Python 3 + packages in python/requirements.txt
# =============================================================================

set -euo pipefail

# ── Resolve project root (script lives in project root or openfoam/cases/...) ─
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}"

# ── Default parameters ────────────────────────────────────────────────────────
OF_CASE="${PROJECT_ROOT}/openfoam/cases/cylinderLES3D"
PYTHON_DIR="${PROJECT_ROOT}/python"
DATA_DIR="${PROJECT_ROOT}/data"
FIGURE_DIR="${DATA_DIR}/figures"

K=10
T_MIN=""            # empty = auto-detect (skip first 20% of simulation time)
ENERGY_THRESH=0.9999
FORCE_COLUMNS="Cd Cl Cm"
SKIP_RECONSTRUCT=0
SKIP_VTK=0
RUN_TAG=""          # optional tag appended to output folder name

# ── Parse arguments ───────────────────────────────────────────────────────────
usage() {
  echo "Usage: $0 [options]"
  echo ""
  echo "Options:"
  echo "  --K N                  Number of clusters (default: ${K})"
  echo "  --t-min T              Skip snapshots before time T (default: auto 20%)"
  echo "  --energy-thresh E      POD energy threshold (default: ${ENERGY_THRESH})"
  echo "  --tag TAG              Append TAG to output folder name"
  echo "  --skip-reconstruct     Skip reconstructPar step"
  echo "  --skip-vtk             Skip foamToVTK export"
  echo "  --of-case PATH         Path to OpenFOAM case (default: ${OF_CASE})"
  echo "  --help                 Show this help"
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --K)               K="$2";             shift 2 ;;
    --t-min)           T_MIN="$2";         shift 2 ;;
    --energy-thresh)   ENERGY_THRESH="$2"; shift 2 ;;
    --tag)             RUN_TAG="$2";       shift 2 ;;
    --skip-reconstruct) SKIP_RECONSTRUCT=1; shift  ;;
    --skip-vtk)        SKIP_VTK=1;         shift  ;;
    --of-case)         OF_CASE="$2";       shift 2 ;;
    --help)            usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# ── Helpers ───────────────────────────────────────────────────────────────────
log()  { echo "[$(date '+%H:%M:%S')]  $*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "'$1' not found. Load OpenFOAM? Activate venv?"; }

# ── Sanity checks ─────────────────────────────────────────────────────────────
[[ -d "${OF_CASE}" ]]    || die "OpenFOAM case not found: ${OF_CASE}"
[[ -d "${PYTHON_DIR}" ]] || die "Python directory not found: ${PYTHON_DIR}"

FORCE_FILE="${OF_CASE}/postProcessing/forceCoeffs/0/coefficient.dat"

CNM_SCRIPT="${PYTHON_DIR}/scripts/run_cnm_analysis.py"
FORCE_PLOT_SCRIPT="${PYTHON_DIR}/scripts/plot_force_coeffs.py"

[[ -f "${CNM_SCRIPT}" ]]        || die "CNM script not found: ${CNM_SCRIPT}"
[[ -f "${FORCE_PLOT_SCRIPT}" ]] || die "Force plot script not found: ${FORCE_PLOT_SCRIPT}"

# ── Python environment ────────────────────────────────────────────────────────
VENV="${PYTHON_DIR}/.venv"
if [[ -f "${VENV}/bin/python3" ]]; then
  PYTHON="${VENV}/bin/python3"
else
  PYTHON="python3"
  log "No venv found at ${VENV} – using system python3"
fi

# ── Output folder ─────────────────────────────────────────────────────────────
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
OUT_TAG="${TIMESTAMP}_K${K}"
[[ -n "${RUN_TAG}" ]] && OUT_TAG="${OUT_TAG}_${RUN_TAG}"
CNM_OUT="${FIGURE_DIR}/cnm_${OUT_TAG}"
mkdir -p "${CNM_OUT}"
mkdir -p "${FIGURE_DIR}"

LOG_FILE="${CNM_OUT}/analyse.log"
exec > >(tee -a "${LOG_FILE}") 2>&1   # tee all output to log file

log "=================================================="
log "  CFD ROM – Post-processing pipeline"
log "  Project root : ${PROJECT_ROOT}"
log "  OF case      : ${OF_CASE}"
log "  Output       : ${CNM_OUT}"
log "  K            : ${K}"
log "  Energy thresh: ${ENERGY_THRESH}"
log "=================================================="

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 – reconstructPar
# ─────────────────────────────────────────────────────────────────────────────
if [[ ${SKIP_RECONSTRUCT} -eq 0 ]]; then
  log "STEP 1/4 – reconstructPar ..."
  need reconstructPar
  cd "${OF_CASE}"

  # Only reconstruct if processor directories exist
  if ls -d processor0 2>/dev/null | grep -q processor0; then
    reconstructPar 2>&1 | tee "${CNM_OUT}/log.reconstructPar"
    log "  reconstructPar done."
  else
    log "  No processor* dirs found – skipping reconstructPar."
  fi
  cd "${PROJECT_ROOT}"
else
  log "STEP 1/4 – reconstructPar SKIPPED (--skip-reconstruct)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 – foamToVTK
# ─────────────────────────────────────────────────────────────────────────────
if [[ ${SKIP_VTK} -eq 0 ]]; then
  log "STEP 2/4 – foamToVTK ..."
  need foamToVTK
  cd "${OF_CASE}"

  VTK_DIR="${OF_CASE}/VTK"
  mkdir -p "${VTK_DIR}"

  foamToVTK \
    -fields '(U p vorticity)' \
    -noPointValues \
    2>&1 | tee "${CNM_OUT}/log.foamToVTK"

  N_VTK=$(find "${VTK_DIR}" -name "*.vtk" 2>/dev/null | wc -l)
  log "  foamToVTK done: ${N_VTK} VTK files in ${VTK_DIR}"
  cd "${PROJECT_ROOT}"
else
  log "STEP 2/4 – foamToVTK SKIPPED (--skip-vtk)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 – Force coefficient plot
# ─────────────────────────────────────────────────────────────────────────────
log "STEP 3/4 – Force coefficient plot ..."

if [[ ! -f "${FORCE_FILE}" ]]; then
  # Try alternative filename used by some OF versions
  ALT="${OF_CASE}/postProcessing/forceCoeffs/0/forceCoeffs.dat"
  if [[ -f "${ALT}" ]]; then
    FORCE_FILE="${ALT}"
  else
    log "  WARNING: forceCoeffs file not found – skipping force plot."
    log "    Tried: ${FORCE_FILE}"
    log "    Tried: ${ALT}"
    FORCE_FILE=""
  fi
fi

if [[ -n "${FORCE_FILE}" ]]; then
  FORCE_PLOT_OUT="${CNM_OUT}/00_force_overview.png"
  "${PYTHON}" "${FORCE_PLOT_SCRIPT}" \
    --file "${FORCE_FILE}" \
    --out  "${FORCE_PLOT_OUT}"
  log "  Force plot written: ${FORCE_PLOT_OUT}"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 – CNM analysis (POD + k-means++ + CMM/CNM)
# ─────────────────────────────────────────────────────────────────────────────
log "STEP 4/4 – CNM analysis (POD + k-means++ + CMM/CNM) ..."

[[ -n "${FORCE_FILE}" ]] || die "Cannot run CNM analysis: no forceCoeffs file found."

# Auto-detect t_min: skip first 20% of simulation time
if [[ -z "${T_MIN}" ]]; then
  T_MIN=$(awk '
    /^[[:space:]]*[0-9]/ {
      if (NR == 2) t_start = $1
      t_end = $1
    }
    END {
      span = t_end - t_start
      printf "%.4f", t_start + 0.20 * span
    }
  ' "${FORCE_FILE}")
  log "  Auto t_min = ${T_MIN}  (skipping first 20% of simulation time)"
else
  log "  t_min = ${T_MIN}  (user-specified)"
fi

# Build force-columns arg
FC_ARGS=""
for col in ${FORCE_COLUMNS}; do
  FC_ARGS="${FC_ARGS} ${col}"
done

"${PYTHON}" "${CNM_SCRIPT}" \
  --input-type   force \
  --input        "${FORCE_FILE}" \
  --out          "${CNM_OUT}" \
  --K            "${K}" \
  --energy-threshold "${ENERGY_THRESH}" \
  --t-min        "${T_MIN}" \
  --force-columns ${FC_ARGS}

log "  CNM analysis done."

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
log ""
log "=================================================="
log "  DONE – Post-processing complete"
log "  Output folder : ${CNM_OUT}"
log "  Log file      : ${LOG_FILE}"
log ""
log "  Generated files:"
find "${CNM_OUT}" -maxdepth 1 -type f | sort | while read -r f; do
  SIZE=$(du -h "$f" | cut -f1)
  log "    ${SIZE}  $(basename "$f")"
done
log ""
log "  To open figures:"
log "    xdg-open ${CNM_OUT}/"
log ""
log "  To re-run with different K:"
log "    $0 --K 15 --t-min ${T_MIN} --skip-reconstruct --skip-vtk"
log "=================================================="
