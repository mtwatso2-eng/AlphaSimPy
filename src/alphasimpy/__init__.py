"""
Canonical AlphaSimPy package namespace.

Re-exports the existing API from the legacy module implementation.
"""

import AlphaSimPy as _legacy

# Re-export only attributes that actually exist to avoid import-time
# failures if legacy __all__ contains stale names.
_legacy_all = getattr(_legacy, "__all__", [])
if _legacy_all:
    __all__ = [name for name in _legacy_all if hasattr(_legacy, name)]
else:
    __all__ = [name for name in dir(_legacy) if not name.startswith("_")]

globals().update({name: getattr(_legacy, name) for name in __all__})
