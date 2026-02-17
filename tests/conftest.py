from pathlib import Path
import sys

# Ensure src-layout package is importable during test collection.
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
