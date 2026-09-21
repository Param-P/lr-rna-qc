import math
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


def test_cpp_engine_if_available():
    if get_backend() != "C++":
        return

    cpp_engine = get_engine("cpp")

    # Empty quality string
    empty_qc = cpp_engine.process_quality_string("")
    assert empty_qc.length == 0
    assert empty_qc.mean_q == 0.0

    # Exact Phred calculation: ASCII 43 ('+') = Q10 (P=0.1), ASCII 63 ('?') = Q30 (P=0.001)
    # Expected mean P = (0.1 + 0.001) / 2 = 0.0505 -> Q = -10 * log10(0.0505) ~= 12.967
    qc = cpp_engine.process_quality_string("+?")
    expected_q = -10.0 * math.log10(0.0505)
    assert qc.length == 2
    assert abs(qc.mean_q - expected_q) < 1e-4

    # Chunk processing
    headers = ["@read1", "@read2"]
    seqs = ["ACGT", "AC"]
    quals = ["++++", "+?"]
    chunk_res = cpp_engine.process_fastq_chunk(headers, seqs, quals)
    assert len(chunk_res) == 2
    assert chunk_res[0].length == 4
    assert abs(chunk_res[0].mean_q - 10.0) < 1e-4
    assert chunk_res[1].length == 2

    # Numerical concordance between C++ and pure-Python fallback
    py_res = fallback.process_quality_string("+?")
    assert abs(qc.mean_q - py_res.mean_q) < 1e-6
