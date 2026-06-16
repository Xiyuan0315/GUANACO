# Installing GUANACO

GUANACO is a pure-Python package distributed as `guanaco-viz`. Nearly every
dependency ships as a prebuilt **binary wheel** on the supported platforms, so a
normal install needs **no C/C++/LLVM compiler** and typically finishes in
**1–3 minutes** (≈ 236 MB download for the core). One small pure-Python
dependency (`dash-cytoscape`) builds from its source archive instantly — still no
compiler required.

---

## 1. Requirements

### Python version

| Python | Supported |
|--------|-----------|
| 3.11   | ✅ (recommended) |
| 3.12   | ✅ |
| 3.13   | ✅ |
| 3.14   | ✅ (verified: stable wheels for the full stack, app runs) |
| 3.15+  | ❌ not yet tested |
| ≤ 3.10 | ❌ — `anndata ≥ 0.12` (and its `read_lazy` streaming) requires Python ≥ 3.11 |

### Operating system / architecture

Installs from prebuilt wheels (no build tools needed) on:

| Platform | Status |
|----------|--------|
| Linux x86-64 (most servers, WSL2) | ✅ |
| Windows 64-bit | ✅ |
| macOS — Apple Silicon (M1/M2/M3) | ✅ |
| macOS — Intel (x86-64) | ✅ |
| Linux ARM64 / aarch64 (e.g. AWS Graviton) | ⚠️ usually works; if a wheel is missing pip builds from source (needs a compiler) |

> GUANACO is a web app — it runs a local server you open in a browser. It works
> the same on a laptop, a lab workstation, an HPC login node, or a cloud VM.

---

## 2. Install (recommended: an isolated environment)

Always install into a **fresh virtual environment** so GUANACO's dependencies
don't collide with other projects.

### Option A — `venv` + `pip` (no extra tools)

```bash
# create and activate an isolated environment
python3.11 -m venv guanaco-env
source guanaco-env/bin/activate        # Windows: guanaco-env\Scripts\activate

# install
pip install --upgrade pip
pip install guanaco-viz
```

### Option B — conda / mamba

```bash
conda create -n guanaco python=3.11 -y
conda activate guanaco
pip install guanaco-viz
```

> Use `conda` only to create the Python environment; install GUANACO itself with
> `pip` (all dependencies are on PyPI).

### Option C — `uv` (fastest)

```bash
uv venv --python 3.11 guanaco-env
source guanaco-env/bin/activate
uv pip install guanaco-viz
```

### Option D — from source (development)

```bash
git clone https://github.com/Systems-Immunometabolism-Lab/guanaco-viz.git
cd guanaco-viz
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e .          # editable install of the core
```

---

## 3. Optional features (extras)

The core install runs the standard single-cell plots from local `.h5ad`/`.h5mu`
files. Heavier or specialized features are **opt-in** so the base install stays
small. Add them with `pip install "guanaco-viz[extra]"`:

| Extra | Adds | Install when you need… |
|-------|------|------------------------|
| `cloud` | fsspec, xarray, dask, s3fs, gcsfs, aiohttp | Loading **any** data from a remote URL — streaming remote `.zarr` (backed mode) **and** fetching a remote `.h5ad`/`.h5mu` over S3 / GCS / HTTPS |
| `tracks` | dash-bio, logomaker, pyjaspar, boto3 | The genome browser (bigWig tracks, motif logos, JASPAR) |
| `datashader` | datashader | The `datashader` render backend for very large embeddings |
| `desktop` | pywebview | Running as a desktop window (`guanaco-desktop`) |
| `notebook` | ipywidgets, ipycytoscape, marimo | Embedding GUANACO in Jupyter / marimo |
| `share` | pycloudflared, flask-compress, brotli | Public sharing via a Cloudflare quick tunnel (`--share`) |
| `server` | gunicorn | Production WSGI deployment (e.g. Plotly Cloud) |
| `all` | everything above | — |

Examples:

```bash
pip install "guanaco-viz[cloud]"             # remote .zarr streaming
pip install "guanaco-viz[cloud,tracks]"      # streaming + genome browser
pip install "guanaco-viz[all]"               # all features
```

> The local app server uses Flask/Werkzeug — you do **not** need `gunicorn`
> (`[server]`) to run GUANACO on your own machine.

---

## 4. Verify the installation

```bash
guanaco --help
python -c "import guanaco; print('GUANACO OK')"
```

---

## 5. Run it

GUANACO is driven by a small JSON config that points at your data. Minimal example
— save as `guanaco.json`:

```json
{
  "MyDataset": {
    "sc_data": "/absolute/path/to/your_data.h5ad",
    "optional_plot_components": ["heatmap", "violin", "dotplot", "stacked-bar"],
    "markers": ["CD3D", "CD8A", "MS4A1"]
  },
  "title": "My single-cell dataset"
}
```

Then launch:

```bash
guanaco -c guanaco.json
# open http://localhost:4399 in your browser
```

Prefer a GUI to build the config? Run `guanaco --config-builder` (needs `tkinter`;
on Debian/Ubuntu: `sudo apt install python3-tk`).

---

## 6. Troubleshooting

**`ERROR: Could not find a version that satisfies the requirement ...` / resolver
conflicts.** You're installing into an environment that already has conflicting
packages. Use a **fresh** virtual environment (Section 2).

**No matching distribution / install fails on Python ≤ 3.10 or ≥ 3.15.**
Use Python 3.11–3.14. GUANACO's floor is 3.11 (`anndata ≥ 0.12` requires it);
3.15 isn't tested yet.

**Loading from a URL fails with a `cloud` extra message.** Remote `sc_data` (an
`http(s)://`, `s3://`, or `gs://` link to a `.zarr`/`.h5ad`/`.h5mu`) needs the
cloud backends: `pip install "guanaco-viz[cloud]"`. Local file paths work with the
core install.

**Install tries to *build* `llvmlite`/`numba` from source and fails.** This only
happens on platforms without a prebuilt wheel (commonly Linux ARM64, or a too-new
Python). Either switch to Python 3.11/3.12, or install build tools
(`llvm`, a C/C++ compiler).

**`ModuleNotFoundError` for a feature** (e.g. `boto3`, `datashader`, `dash_bio`).
That feature lives in an extra — install it, e.g. `pip install "guanaco-viz[tracks]"`.
The error message names the right extra.

**`guanaco --config-builder` fails with a `tkinter` error.** `tkinter` ships with
most Python builds but not all Linux ones — install `python3-tk` (Debian/Ubuntu)
or `python3-tkinter` (Fedora).
