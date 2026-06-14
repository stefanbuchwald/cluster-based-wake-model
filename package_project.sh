#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# package_project.sh
#
# Create a source-only ZIP of the project.
#
# Keeps:
#   - README/docs
#   - Python scripts and requirements
#   - OpenFOAM case setup: 0/, constant/, system/, Allrun/Allclean/etc.
#   - helper scripts
#
# Excludes:
#   - generated files below data/
#   - VTK/
#   - processor*/
#   - postProcessing/
#   - reconstructed time directories
#   - logs
#   - mesh output constant/polyMesh
#   - venv
#   - generated images, arrays, VTK files, videos, archives, and logs
# =============================================================================

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_NAME="$(basename "$PROJECT_ROOT")"
PARENT_DIR="$(dirname "$PROJECT_ROOT")"

STAMP="$(date +%Y%m%d_%H%M%S)"
DIST_DIR="$PROJECT_ROOT/dist"
OUT_ZIP="$DIST_DIR/${PROJECT_NAME}_code_only_${STAMP}.zip"

mkdir -p "$DIST_DIR"

if ! command -v zip >/dev/null 2>&1; then
    echo "ERROR: zip is not installed."
    echo "Install with:"
    echo "  sudo apt install zip"
    exit 1
fi

echo "==> Project root : $PROJECT_ROOT"
echo "==> Output zip   : $OUT_ZIP"
echo ""

# Optional tree snapshot of what should remain in the source package.
if command -v tree >/dev/null 2>&1; then
    echo "==> Writing PROJECT_TREE_PACKAGE.txt"

    tree -a \
        -I ".git|.venv|__pycache__|.pytest_cache|.ipynb_checkpoints|dist|raw|processed|figures|results|VTK|processor*|postProcessing|logs|polyMesh|*.vtk|*.vtu|*.pvd|*.foam|*.msh|*.png|*.jpg|*.jpeg|*.npy|*.npz|*.zip|*.avi|*.mp4|*.pvsm|TEMP|[1-9]*|[0-9]*.[0-9]*" \
        "$PROJECT_ROOT" > "$PROJECT_ROOT/PROJECT_TREE_PACKAGE.txt"
else
    echo "==> tree not installed; skipping PROJECT_TREE_PACKAGE.txt"
fi

cd "$PARENT_DIR"

echo "==> Creating code-only package..."

zip -r "$OUT_ZIP" "$PROJECT_NAME" \
    -x "$PROJECT_NAME/.git/*" \
       "$PROJECT_NAME/.venv/*" \
       "$PROJECT_NAME/.agents/*" \
       "$PROJECT_NAME/.codex/*" \
       "$PROJECT_NAME/dist/*" \
       "$PROJECT_NAME/PROJECT_TREE.txt" \
       "$PROJECT_NAME/PROJECT_TREE_PACKAGE.txt" \
       "$PROJECT_NAME/**/__pycache__/*" \
       "$PROJECT_NAME/**/.pytest_cache/*" \
       "$PROJECT_NAME/**/.ipynb_checkpoints/*" \
       "$PROJECT_NAME/data/raw/*" \
       "$PROJECT_NAME/data/processed/*" \
       "$PROJECT_NAME/data/figures/*" \
       "$PROJECT_NAME/data/results/*" \
       "$PROJECT_NAME/**/VTK/*" \
       "$PROJECT_NAME/**/processor*/*" \
       "$PROJECT_NAME/**/postProcessing/*" \
       "$PROJECT_NAME/**/logs/*" \
       "$PROJECT_NAME/**/log.*" \
       "$PROJECT_NAME/**/[1-9]*/*" \
       "$PROJECT_NAME/**/[0-9]*.[0-9]*/*" \
       "$PROJECT_NAME/**/constant/polyMesh/*" \
       "$PROJECT_NAME/**/0/Lambda2" \
       "$PROJECT_NAME/**/0/Q" \
       "$PROJECT_NAME/**/0/vorticity" \
       "$PROJECT_NAME/**/*.vtk" \
       "$PROJECT_NAME/**/*.vtu" \
       "$PROJECT_NAME/**/*.pvd" \
       "$PROJECT_NAME/**/*.foam" \
       "$PROJECT_NAME/**/*.msh" \
       "$PROJECT_NAME/**/*.png" \
       "$PROJECT_NAME/**/*.jpg" \
       "$PROJECT_NAME/**/*.jpeg" \
       "$PROJECT_NAME/**/*.npy" \
       "$PROJECT_NAME/**/*.npz" \
       "$PROJECT_NAME/**/*.avi" \
       "$PROJECT_NAME/**/*.mp4" \
       "$PROJECT_NAME/**/*.pvsm" \
       "$PROJECT_NAME/**/*.zip" \
       "$PROJECT_NAME/**/TEMP/*" \
       "$PROJECT_NAME/**/TEMP" \
       "$PROJECT_NAME/*.zip"

echo ""
echo "==> Done."
echo "Created:"
echo "    $OUT_ZIP"
echo ""
echo "Size:"
du -h "$OUT_ZIP"

echo ""
echo "Inspect with:"
echo "    unzip -l \"$OUT_ZIP\" | less"
echo ""
echo "Checking archive for forbidden generated files..."

FORBIDDEN_PATTERN='/(VTK|processor[^/]*|postProcessing|polyMesh|raw|processed|figures|results)/|\.(vtk|vtu|pvd|npy|npz|png|jpe?g|avi|mp4|msh|foam|pvsm)$'
ARCHIVE_LIST="$(unzip -l "$OUT_ZIP")"

if grep -Eiq "$FORBIDDEN_PATTERN" <<< "$ARCHIVE_LIST"; then
    echo "ERROR: generated data found in source archive:"
    grep -Ei "$FORBIDDEN_PATTERN" <<< "$ARCHIVE_LIST"
    exit 1
fi

echo "OK: no generated data found."
