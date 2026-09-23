import math
from pathlib import Path
from lr_rna_qc.core import get_backend, get_engine, stream_fastq, scan_fastq, is_gzipped
from lr_rna_qc import fallback
from tests.make_synthetic_fastq import generate_synthetic_fastq


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

    empty_qc = cpp_engine.process_quality_string("")
    assert empty_qc.length == 0
    assert empty_qc.mean_q == 0.0

    # Exact Phred calculation: ASCII 43 ('+') = Q10 (P=0.1), ASCII 63 ('?') = Q30 (P=0.001)
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


def test_streaming_and_scan_pipeline(tmp_path: Path):
    fq_path = generate_synthetic_fastq(tmp_path / "test.fastq.gz", num_reads=35)
    assert is_gzipped(fq_path) is True

    # Stream in chunks of 10 reads (3 chunks of 10, 1 chunk of 5)
    chunks = list(stream_fastq(fq_path, chunk_size=10))
    assert len(chunks) == 4
    assert len(chunks[0][0]) == 10
    assert len(chunks[3][0]) == 5

    # Run scan with filtering and output file
    out_filtered = tmp_path / "filtered.fastq.gz"
    summary_py, lens_py, qs_py = scan_fastq(
        fq_path,
        min_q=10.0,
        min_len=200,
        out_filtered=out_filtered,
        backend_override="python",
        chunk_size=10,
    )
    assert summary_py.total_reads == 35
    assert summary_py.pass_reads + summary_py.fail_reads == 35
    assert summary_py.n50 > 0
    assert out_filtered.exists()

    # Compare with C++ backend if available
    if get_backend() == "C++":
        summary_cpp, lens_cpp, qs_cpp = scan_fastq(
            fq_path,
            min_q=10.0,
            min_len=200,
            backend_override="cpp",
            chunk_size=10,
        )
        assert summary_cpp.total_reads == summary_py.total_reads
        assert summary_cpp.pass_reads == summary_py.pass_reads
        assert summary_cpp.n50 == summary_py.n50
        assert abs(summary_cpp.mean_q - summary_py.mean_q) < 1e-2
