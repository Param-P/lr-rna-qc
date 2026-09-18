"""Unit tests for core engine loading and fallback behavior."""

from lr_rna_qc.core import get_backend, get_engine
from lr_rna_qc import fallback


def test_backend_availability():
    backend = get_backend()
    assert backend in ("C++", "Python")


def test_fallback_engine_execution():
    engine = get_engine("python")
    assert engine is fallback
    res = engine.process_quality_string("IIII")
    assert res.length == 4
    assert res.mean_q > 0.0
