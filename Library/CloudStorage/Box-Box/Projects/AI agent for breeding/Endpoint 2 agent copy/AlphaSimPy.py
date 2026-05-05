"""
Backward-compatible shim for legacy imports.

Prefer `import alphasimpy` moving forward.
"""

import os
import sys

_SRC_DIR = os.path.join(os.path.dirname(__file__), "src")
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from alphasimpy import *  # noqa: F401,F403
