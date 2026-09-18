#include "fast_qc.hpp"
#include <cmath>
#include <array>
#include <algorithm>

namespace lr_rna_qc {

// Precalculated Look-Up Table (LUT) mapping ASCII-33 Phred characters to error probabilities
// P = 10^(-Q / 10), where Q = ASCII - 33
static const std::array<float, 256> ASCII_TO_PROB = []() {
    std::array<float, 256> table{};
    for (size_t i = 0; i < 256; ++i) {
        if (i >= 33 && i <= 126) {
            double q = static_cast<double>(i - 33);
            table[i] = static_cast<float>(std::pow(10.0, -q / 10.0));
        } else {
            table[i] = 1.0f; // Maximum error probability for invalid / unassigned
        }
    }
    return table;
}();

ReadQC process_quality_string(std::string_view qual) {
    if (qual.empty()) {
        return {0.0, 0};
    }
    double sum_p = 0.0;
    for (unsigned char c : qual) {
        sum_p += ASCII_TO_PROB[c];
    }
    double mean_p = sum_p / static_cast<double>(qual.size());
    // Avoid log10(0) if mean_p is zero
    double mean_q = (mean_p > 0.0) ? -10.0 * std::log10(mean_p) : 93.0;
    return {mean_q, qual.size()};
}

std::vector<ReadQC> process_fastq_chunk(
    const std::vector<std::string>& /*headers*/,
    const std::vector<std::string>& /*seqs*/,
    const std::vector<std::string>& quals
) {
    std::vector<ReadQC> results;
    results.reserve(quals.size());
    for (const auto& q : quals) {
        results.push_back(process_quality_string(q));
    }
    return results;
}

} // namespace lr_rna_qc
