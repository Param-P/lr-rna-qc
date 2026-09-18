#pragma once

#include <string>
#include <string_view>
#include <vector>
#include <cstddef>

namespace lr_rna_qc {

struct ReadQC {
    double mean_q;
    size_t length;
};

// Evaluates a single Phred quality string using exact error probability LUT
ReadQC process_quality_string(std::string_view qual);

// Evaluates a batch of sequences and quality strings
std::vector<ReadQC> process_fastq_chunk(
    const std::vector<std::string>& headers,
    const std::vector<std::string>& seqs,
    const std::vector<std::string>& quals
);

} // namespace lr_rna_qc
