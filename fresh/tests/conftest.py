import sys
from pathlib import Path


FRESH_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = FRESH_ROOT / "src"

if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))
