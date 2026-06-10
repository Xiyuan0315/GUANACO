#!/usr/bin/env bash
#
# Build GUANACO as a standalone desktop app with Nuitka.
#
#   macOS  -> dist/GUANACO.app  (double-clickable bundle)
#   Linux  -> dist/desktop_entry.dist/  (folder with the executable)
#
# Run it from the repo root inside your environment, e.g.:
#   pixi run bash scripts/build_desktop.sh
#
# Notes
# -----
# * Data is loaded at import time, so the app needs its bundled web assets
#   (Dash/Plotly JS) AND the package data (assets/, cvd_color.json).
# * datashader/numba/llvmlite are excluded: they are only used by the optional
#   `embedding_render_backend: datashader` mode and are very hard to compile.
#   Build with the default `scattergl` backend. Set NUITKA_WITH_DATASHADER=1 to
#   attempt including them.
# * The first build is slow (10-30 min). If the app launches but errors with a
#   "No module named X", add `--include-module=X` (or `--include-package=X`)
#   below and rebuild.

set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$(pwd)"
ENTRY="scripts/desktop_entry.py"
ICON="src/guanaco/assets/logo.png"

# --- ensure Nuitka is available --------------------------------------------
python -c "import nuitka" 2>/dev/null || pip install "nuitka>=2.4" ordered-set
python -c "import webview" 2>/dev/null || pip install "pywebview>=5.0"

# --- packages whose bundled data files must ship (JS/CSS/templates/data) ----
DATA_PACKAGES=(
  guanaco                       # assets/, utils/cvd_color.json
  dash dash_bootstrap_components dash_bio dash_ag_grid dash_draggable
  dash_cytoscape plotly
  botocore certifi              # AWS endpoint data + CA bundle (S3 configs)
  matplotlib                    # mpl-data (used by logomaker)
)

# --- packages Nuitka may under-analyse (lazy / dynamic imports) -------------
FORCE_PACKAGES=(
  guanaco anndata mudata muon h5py scipy numpy pandas statsmodels
  logomaker pyfaidx pyjaspar Bio
  flask werkzeug jinja2 click
  boto3 botocore aiobotocore s3transfer
  webview
)

# --- things we never want to follow into ------------------------------------
EXCLUDE=(
  pytest IPython ipykernel jupyterlab notebook
)
if [ "${NUITKA_WITH_DATASHADER:-0}" != "1" ]; then
  EXCLUDE+=(datashader numba llvmlite)
fi

ARGS=(
  --standalone
  --assume-yes-for-downloads
  --output-dir=dist
  --enable-plugin=tk-inter        # config wizard (tkinter)
  --include-package-data=guanaco
)

for p in "${DATA_PACKAGES[@]}";  do ARGS+=("--include-package-data=$p"); done
for p in "${FORCE_PACKAGES[@]}"; do ARGS+=("--include-package=$p");      done
for p in "${EXCLUDE[@]}";        do ARGS+=("--nofollow-import-to=$p");   done

# --- platform specifics -----------------------------------------------------
case "$(uname -s)" in
  Darwin)
    ARGS+=(
      --macos-create-app-bundle
      --macos-app-name=GUANACO
      --macos-app-icon="$ICON"
      --macos-signed-app-name=com.guanaco.viz
    )
    ;;
esac

echo "Building GUANACO desktop app with Nuitka..."
echo "  entry : $ENTRY"
echo "  root  : $ROOT"
python -m nuitka "${ARGS[@]}" "$ENTRY"

echo
echo "Done. Output is under: $ROOT/dist/"
echo "  macOS: open dist/GUANACO.app   (or:  open dist/*.app)"
echo "  Linux: ./dist/desktop_entry.dist/desktop_entry"
