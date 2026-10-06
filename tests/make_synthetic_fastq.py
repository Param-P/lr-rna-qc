import gzip
import random
from pathlib import Path


def generate_synthetic_fastq(
    output_path: str | Path,
    num_reads: int = 100,
    seed: int = 42,
    compress: bool = True,
    condition: str = "standard",
) -> Path:
    """Generate deterministic synthetic Nanopore RNA reads for testing edge cases."""
    random.seed(seed)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    is_gz = compress or str(path).endswith(".gz")
    opener = gzip.open if is_gz else open

    with opener(path, "wt", encoding="utf-8") as fh:
        for i in range(num_reads):
            read_id = f"@nanopore_rna_read_{i+1} sample={condition}"

            # Simulate transcriptome length distributions
            if condition == "mutant":
                # Simulated splicing mutant (e.g. truncated transcripts / nonsense decay)
                length = int(random.gauss(600, 200))
                length = max(80, length)
            elif i % 25 == 0:
                # Ultra-long full-length mRNA (10-18kb)
                length = random.randint(10000, 18000)
            elif i % 20 == 0:
                # Degraded debris / short primer artifact (50-120bp)
                length = random.randint(50, 120)
            else:
                # Standard transcript distribution (~1200-2500bp)
                length = int(random.gauss(1400, 450))
                length = max(150, length)

            seq = "".join(random.choices("ACGT", k=length))

            # Quality distribution: typical Nanopore Q8-Q15 (ASCII 41..48)
            # Add edge cases: localized low-Q bursts (Q=2 '#')
            if condition == "mutant":
                qual_chars = [chr(random.randint(37, 44)) for _ in range(length)]
            elif i % 12 == 0:
                # Severe burst of Q=2 across entire read
                qual_chars = [chr(35)] * length
            elif i % 5 == 0:
                # Front half degraded, back half high quality
                half = length // 2
                qual_chars = [chr(35)] * half + [chr(random.randint(43, 48))] * (length - half)
            else:
                qual_chars = [chr(random.randint(41, 48)) for _ in range(length)]

            qual = "".join(qual_chars)
            fh.write(f"{read_id}\n{seq}\n+\n{qual}\n")

    return path


if __name__ == "__main__":
    out_file = generate_synthetic_fastq("tests/fixtures/sample.fastq.gz", num_reads=100)
    print(f"Generated synthetic fixture at {out_file}")
