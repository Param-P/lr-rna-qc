#include "fast_qc.hpp"
#include <array>
#include <cmath>

namespace lr_rna_qc {

// Precalculated Look-Up Table (LUT) mapping ASCII-33 Phred characters to error probabilities.
// P = 10^(-Q / 10), where Q = ASCII - 33 
static const std::array<float, 256> ASCII_TO_PROB = []() {
    std::array<float, 256> table{};
    for (size_t i = 0; i < 256; i++) {
        if (i >= 33 && i <= 126) {
            double q = static_cast<double>(i - 33);
            table[i] = static_cast<float>(std::pow(10.0, -q / 10.0));
        } else {
            table[i] = 1.0f; // Default for non-phred chars
        }
    }
    return table;
}();

ReadQC process_quality_string(std::string_view qual) {
    const size_t len = qual.size();
    if (len == 0) {
        return {0.0, 0};
    }

    double sum_prob = 0.0;
    const char* ptr = qual.data();
    for (size_t i = 0; i < len; ++i) {
        sum_prob += ASCII_TO_PROB[static_cast<unsigned char>(ptr[i])];
    }

    const double mean_prob = sum_prob / static_cast<double>(len);
    // Q = -10 * log10(P_error); clamp to Phred 93.0 on numerical zero
    const double mean_q = (mean_prob > 0.0) ? -10.0 * std::log10(mean_prob) : 93.0;

    return {mean_q, len};
}

std::vector<ReadQC> process_fastq_chunk(
    const std::vector<std::string>& /*headers*/,
    const std::vector<std::string>& /*seqs*/,
    const std::vector<std::string>& quals
) {
    std::vector<ReadQC> results;
    results.reserve(quals.size());

    for (const std::string& qual : quals) {
        results.push_back(process_quality_string(qual));
    }

    return results;
}

} // namespace lr_rna_qc
