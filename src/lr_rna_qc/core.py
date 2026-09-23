import gzip
import logging
from pathlib import Path
from typing import Any, Iterator, List, Tuple
from lr_rna_qc.metrics import RunSummary, compute_run_summary

logger = logging.getLogger(__name__)

BACKEND = "Python"
_engine: Any = None

try:
    from lr_rna_qc import _fast_qc  # type: ignore
    _engine = _fast_qc
    BACKEND = "C++"
    logger.debug("Loaded C++ acceleration engine (_fast_qc).")
except ImportError:
    from lr_rna_qc import fallback
    _engine = fallback
    BACKEND = "Python"
    logger.info("C++ engine not available; using pure-Python fallback.")

def get_backend() -> str:
    """Return the active engine backend ('C++' or 'Python')."""
    return BACKEND


def get_engine(backend_override: str | None = None) -> Any:
    """Return computation engine according to selection or default."""
    if backend_override == "python":
        from lr_rna_qc import fallback
        return fallback
    elif backend_override == "cpp":
        if BACKEND != "C++":
            raise RuntimeError("C++ backend requested, but _fast_qc is not compiled or available.")
        return _engine
    return _engine


def is_gzipped(filepath: str | Path) -> bool:
    with open(filepath, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b" # Gzip magic bytes


def stream_fastq(
    filepath: str | Path,
    chunk_size: int = 10000,
) -> Iterator[Tuple[List[str], List[str], List[str]]]:
    """Stream FASTQ records in fixed-size chunks with O(1) memory footprint."""
    path = Path(filepath)
    opener = gzip.open if is_gzipped(path) else open

    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        headers: List[str] = []
        seqs: List[str] = []
        quals: List[str] = []

        for line in fh:
            header = line.rstrip("\r\n")
            if not header:
                continue
            try:
                seq = next(fh).rstrip("\r\n")
                _ = next(fh)  # '+' separator line
                qual = next(fh).rstrip("\r\n")
            except StopIteration:
                break

            headers.append(header)
            seqs.append(seq)
            quals.append(qual)

            if len(headers) >= chunk_size:
                yield headers, seqs, quals
                headers = []
                seqs = []
                quals = []

        if headers:
            yield headers, seqs, quals


def scan_fastq(
    filepath: str | Path,
    min_q: float = 0.0,
    min_len: int = 0,
    out_filtered: str | Path | None = None,
    backend_override: str | None = None,
    chunk_size: int = 10000,
) -> Tuple[RunSummary, List[int], List[float]]:
    """Scan FASTQ file through active engine, applying quality/length filters."""
    engine = get_engine(backend_override)

    all_lengths: List[int] = []
    all_mean_qs: List[float] = []
    pass_reads = 0
    fail_reads = 0

    out_fh = None
    if out_filtered is not None:
        out_path = Path(out_filtered)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if str(out_path).endswith(".gz"):
            out_fh = gzip.open(out_path, "wt", encoding="utf-8")
        else:
            out_fh = open(out_path, "wt", encoding="utf-8")

    try:
        for headers, seqs, quals in stream_fastq(filepath, chunk_size=chunk_size):
            qc_results = engine.process_fastq_chunk(headers, seqs, quals)

            for i in range(len(qc_results)):
                qc = qc_results[i]
                length = qc.length
                mean_q = qc.mean_q

                all_lengths.append(length)
                all_mean_qs.append(mean_q)

                if mean_q >= min_q and length >= min_len:
                    pass_reads += 1
                    if out_fh is not None:
                        out_fh.write(f"{headers[i]}\n{seqs[i]}\n+\n{quals[i]}\n")
                else:
                    fail_reads += 1
    finally:
        if out_fh is not None:
            out_fh.close()

    summary = compute_run_summary(all_lengths, all_mean_qs, pass_reads, fail_reads)
    return summary, all_lengths, all_mean_qs
