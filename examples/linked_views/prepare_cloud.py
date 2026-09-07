"""Build an allowlisted upload folder containing the current GUANACO wheel.

Run from this checkout with a Python environment containing pip and setuptools.
Only the bundled public spatial dataset is copied, not other research data,
notebook outputs, credentials or virtual environments.
"""

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def prepare():
    examples = Path(__file__).resolve().parent
    repository = examples.parents[1]
    spatial = examples / "data" / "visium_hne_spatial.h5ad"
    if not spatial.is_file():
        raise FileNotFoundError("Run prepare_spatial.py to prepare the original spatial example first.")
    destination = Path(tempfile.mkdtemp(prefix="guanaco-showcase-"))
    subprocess.run([
        sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation",
        str(repository), "--wheel-dir", str(destination),
    ], check=True)
    wheel = next(destination.glob("guanaco_viz-*.whl"))
    for name in ("app.py", "showcase.py", "cases.py", "demo_data.py", "synthetic_genes.gtf", "README.md", ".gitignore"):
        shutil.copy2(examples / name, destination / name)
    shutil.copytree(examples / "showcase_assets", destination / "showcase_assets")
    (destination / "data").mkdir()
    shutil.copy2(spatial, destination / "data" / spatial.name)
    (destination / "requirements.txt").write_text(f"./{wheel.name}\n", encoding="utf-8")
    print(f"\nUpload this folder to Plotly Cloud: {destination}")
    print("Main file: app.py | Python: 3.12 | No environment variables required")
    return destination


if __name__ == "__main__":
    prepare()
