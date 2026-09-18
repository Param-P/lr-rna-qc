"""Generator for deterministic synthetic FASTQ and FASTQ.GZ test fixtures."""

import gzip
import random
from pathlib import Path


def generate_synthetic_fastq(
    output_path: str | Path,
    num_reads: int = 100,
    seed: int = 42,
    compress: bool = True,
) -> Path:
    """Generate deterministic synthetic Nanopore RNA reads for testing."""
    random.seed(seed)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    opener = gzip.open if compress or str(path).endswith(".gz") else open

    # Length distribution simulating Nanopore direct RNA / cDNA:
    # Most reads between 500 - 3000bp, some short debris (50bp), some ultra-long (15-20kb)
    with opener(path, "wt") as f:
        for i in range(num_reads):
            read_id = f"@synthetic_read_{i+1}"
            if i % 20 == 0:
                # Short debris
                length = random.randint(50, 150)
            elif i % 25 == 0:
                # Ultra-long transcript
                length = random.randint(10000, 18000)
            else:
                length = int(random.gauss(1200, 400))
                length = max(100, length)

            seq = "".join(random.choices("ACGT", k=length))

            # Quality scores: Nanopore reads typical Q range 8-15 (ASCII 41 - 48)
            # Add occasional low quality burst (ASCII 35 '#', Q=2)
            if i % 15 == 0:
                qual_chars = [chr(35)] * length
            else:
                qual_chars = [chr(random.randint(41, 48)) for _ in range(length)]
            qual = "".join(qual_chars)

            f.write(f"{read_id}\n{seq}\n+\n{qual}\n")

    return path


if __name__ == "__main__":
    out = generate_synthetic_fastq("tests/fixtures/sample.fastq.gz", num_reads=50)
    print(f"Generated synthetic fixture at {out}")
