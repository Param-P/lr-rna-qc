"""Core bridge loading C++ fast QC engine with transparent Python fallback."""

import logging
from typing import Tuple, Any

logger = logging.getLogger(__name__)

BACKEND = "Python"
_engine: Any = None

try:
    from lr_rna_qc import _fast_qc  # type: ignore
    _engine = _fast_qc
    BACKEND = "C++"
    logger.debug("Successfully loaded C++ acceleration engine (_fast_qc).")
except ImportError:
    from lr_rna_qc import fallback
    _engine = fallback
    BACKEND = "Python"
    logger.info("C++ engine not available; using pure-Python fallback.")

def get_backend() -> str:
    """Return the active engine backend ('C++' or 'Python')."""
    return BACKEND

def get_engine(backend_override: str | None = None):
    """Return the requested or active computation engine."""
    if backend_override == "python":
        from lr_rna_qc import fallback
        return fallback
    elif backend_override == "cpp":
        if BACKEND != "C++":
            raise RuntimeError("C++ backend requested, but lr_rna_qc._fast_qc is not compiled or available.")
        return _engine
    return _engine
