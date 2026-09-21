#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace lr_rna_qc {

struct ReadQC {
    double mean_q{0.0};
    size_t length{0};

    bool operator==(const ReadQC& other) const {
        return length == other.length && (mean_q == other.mean_q);
    }
};

// Computes exact mean Phred quality score and length from ASCII-33 quality string
ReadQC process_quality_string(std::string_view qual);

// Batch processes a chunk of FASTQ records
std::vector<ReadQC> process_fastq_chunk(
    const std::vector<std::string>& headers,
    const std::vector<std::string>& seqs,
    const std::vector<std::string>& quals
);

} // namespace lr_rna_qc
