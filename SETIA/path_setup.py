from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent

sys.path.insert(0, str(ROOT / "helper"))
sys.path.insert(0, str(ROOT / "dynamics_module"))