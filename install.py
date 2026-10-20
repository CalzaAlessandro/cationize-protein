#!/usr/bin/env python3
"""Bootstrap this checkout without changing the Python used to start it."""
import sys
from pathlib import Path

sys.dont_write_bytecode = True
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from cationize_protein.dependencies import installer_main

if __name__ == "__main__":
    raise SystemExit(installer_main(root=PROJECT_ROOT))
