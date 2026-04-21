"""
Shared pytest configuration for bin/ script unit tests.

Adds the repository's bin/ directory to sys.path so test modules can import
pipeline scripts directly (e.g. ``import cellmeasurement as cm``).
"""

import sys
from pathlib import Path

# Insert bin/ at the front of sys.path once, for the whole test session.
_bin_dir = Path(__file__).resolve().parents[2] / "bin"
if str(_bin_dir) not in sys.path:
    sys.path.insert(0, str(_bin_dir))
