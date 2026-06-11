#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace gapstats {

// ─────────────────────────────────────────────────────────────────────────────
// BAM CIGAR constants (mirror BAM spec — no htslib dependency)
// ─────────────────────────────────────────────────────────────────────────────
static constexpr int CIGAR_MATCH    = 0;  // M
static constexpr int CIGAR_INS      = 1;  // I
static constexpr int CIGAR_DEL      = 2;  // D
static constexpr int CIGAR_REF_SKIP = 3;  // N
static constexpr int CIGAR_SOFT     = 4;  // S
static constexpr int CIGAR_HARD     = 5;  // H
static constexpr int CIGAR_PAD      = 6;  // P
static constexpr int CIGAR_EQUAL    = 7;  // =
static constexpr int CIGAR_DIFF     = 8;  // X

/// Decode CIGAR operation type from a raw BAM uint32_t element.
inline int cigar_op(uint32_t c) { return static_cast<int>(c & 0xfu); }

/// Decode CIGAR operation length from a raw BAM uint32_t element.
inline int cigar_len(uint32_t c) { return static_cast<int>(c >> 4u); }

/// Encode a CIGAR operation into a raw BAM uint32_t element.
inline uint32_t cigar_encode(int len, int op)
{
    return (static_cast<uint32_t>(len) << 4u) | static_cast<uint32_t>(op);
}

// ─────────────────────────────────────────────────────────────────────────────
// Pure utility functions
// ─────────────────────────────────────────────────────────────────────────────

/// Strip directory path and all extensions from a file path.
/// "path/to/HG002.minimap2.20x.cram" → "HG002"
inline std::string stem(const std::string& path)
{
    namespace fs = std::filesystem;
    std::string s = fs::path(path).filename().string();
    auto dot = s.find('.');
    return (dot == std::string::npos) ? s : s.substr(0, dot);
}

/**
 * Count CIGAR gap operations (INS or DEL) with length >= min_gap_size.
 *
 * @param cigar       Pointer to the raw BAM CIGAR array (uint32_t per op).
 * @param n_cigar     Number of CIGAR operations.
 * @param min_gap_size Minimum op length to count (inclusive).
 * @return Number of qualifying insertions + deletions.
 */
inline int count_gaps(const uint32_t* cigar, uint32_t n_cigar, int min_gap_size)
{
    int gaps = 0;
    for (uint32_t i = 0; i < n_cigar; ++i) {
        const int op  = cigar_op(cigar[i]);
        const int len = cigar_len(cigar[i]);
        if ((op == CIGAR_INS || op == CIGAR_DEL) && len >= min_gap_size)
            ++gaps;
    }
    return gaps;
}

/**
 * Compute decantiles P10, P20, …, P100 (10 values) using the nearest-rank
 * method.  The vector is sorted in-place.
 *
 * @return Array of 10 doubles [P10, P20, …, P100].  All zeros if v is empty.
 */
template <typename T>
std::array<double, 10> decantiles(std::vector<T>& v)
{
    std::array<double, 10> result{};
    if (v.empty()) return result;

    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();

    for (int k = 1; k <= 10; ++k) {
        std::size_t idx = static_cast<std::size_t>(
            std::ceil(k * 0.1 * static_cast<double>(n))) - 1;
        if (idx >= n) idx = n - 1;
        result[k - 1] = static_cast<double>(v[idx]);
    }
    return result;
}

} // namespace gapstats
