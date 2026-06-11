#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <numeric>

#include "core.hpp"

using namespace gapstats;
using Catch::Matchers::WithinRel;

// ─────────────────────────────────────────────────────────────────────────────
// stem()
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("stem strips directory and all extensions", "[stem]")
{
    CHECK(stem("/path/to/HG002.minimap2.20x.subset.cram") == "HG002");
    CHECK(stem("sample.bam") == "sample");
    CHECK(stem("/a/b/foo") == "foo");
    CHECK(stem("no_ext") == "no_ext");
    CHECK(stem("/deep/path/multi.dots.here.tsv") == "multi");
}

// ─────────────────────────────────────────────────────────────────────────────
// CIGAR encoding/decoding
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("cigar_op and cigar_len decode raw BAM uint32_t", "[cigar]")
{
    // encode: len=100, op=DEL(2) → (100 << 4) | 2 = 1602
    const uint32_t c = cigar_encode(100, CIGAR_DEL);
    CHECK(cigar_op(c) == CIGAR_DEL);
    CHECK(cigar_len(c) == 100);
}

TEST_CASE("cigar_encode round-trips correctly", "[cigar]")
{
    for (int op = 0; op <= 8; ++op) {
        for (int len : {1, 5, 50, 1000, 65535}) {
            const uint32_t c = cigar_encode(len, op);
            CHECK(cigar_op(c) == op);
            CHECK(cigar_len(c) == len);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// count_gaps()
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("count_gaps empty cigar returns 0", "[count_gaps]")
{
    CHECK(count_gaps(nullptr, 0, 1) == 0);
}

TEST_CASE("count_gaps counts only INS/DEL ops >= threshold", "[count_gaps]")
{
    // M(50)  D(10)  I(3)  D(5)  I(7)
    // threshold=5: D(10)✓  I(3)✗  D(5)✓  I(7)✓  → 3
    const std::vector<uint32_t> cigar = {
        cigar_encode(50, CIGAR_MATCH),
        cigar_encode(10, CIGAR_DEL),
        cigar_encode( 3, CIGAR_INS),
        cigar_encode( 5, CIGAR_DEL),
        cigar_encode( 7, CIGAR_INS),
    };
    CHECK(count_gaps(cigar.data(), (uint32_t)cigar.size(), 5) == 3);
    CHECK(count_gaps(cigar.data(), (uint32_t)cigar.size(), 6) == 2);  // D(10), I(7)
    CHECK(count_gaps(cigar.data(), (uint32_t)cigar.size(), 8) == 1);  // D(10)
    CHECK(count_gaps(cigar.data(), (uint32_t)cigar.size(), 11) == 0);
}

TEST_CASE("count_gaps ignores non-gap ops (M, S, H, N, =, X)", "[count_gaps]")
{
    const std::vector<uint32_t> cigar = {
        cigar_encode(100, CIGAR_MATCH),
        cigar_encode( 10, CIGAR_SOFT),
        cigar_encode(  5, CIGAR_HARD),
        cigar_encode( 20, CIGAR_REF_SKIP),
        cigar_encode( 30, CIGAR_EQUAL),
        cigar_encode(  8, CIGAR_DIFF),
    };
    CHECK(count_gaps(cigar.data(), (uint32_t)cigar.size(), 1) == 0);
}

TEST_CASE("count_gaps threshold=1 counts all INS and DEL", "[count_gaps]")
{
    const std::vector<uint32_t> cigar = {
        cigar_encode(1, CIGAR_INS),
        cigar_encode(1, CIGAR_DEL),
        cigar_encode(100, CIGAR_MATCH),
        cigar_encode(2, CIGAR_INS),
    };
    CHECK(count_gaps(cigar.data(), (uint32_t)cigar.size(), 1) == 3);
}

// ─────────────────────────────────────────────────────────────────────────────
// decantiles()
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("decantiles on empty vector returns all zeros", "[decantiles]")
{
    std::vector<int> v;
    const auto d = decantiles(v);
    for (auto x : d) CHECK(x == 0.0);
}

TEST_CASE("decantiles on single element returns that value for all", "[decantiles]")
{
    std::vector<int> v = {42};
    const auto d = decantiles(v);
    for (auto x : d) CHECK(x == 42.0);
}

TEST_CASE("decantiles on 10 sorted values [1..10]", "[decantiles]")
{
    // Nearest-rank: P10 = ceil(0.1*10)-1 = 0 → v[0]=1, P100 = v[9]=10
    std::vector<int> v = {5, 3, 1, 9, 2, 7, 4, 6, 8, 10};
    const auto d = decantiles(v);
    CHECK(d[0] == 1.0);   // P10
    CHECK(d[4] == 5.0);   // P50
    CHECK(d[9] == 10.0);  // P100
}

TEST_CASE("decantiles values are monotonically non-decreasing", "[decantiles]")
{
    std::vector<int> v = {8, 2, 5, 1, 9, 3, 7, 4, 6, 10, 11, 0};
    const auto d = decantiles(v);
    for (int i = 1; i < 10; ++i)
        CHECK(d[i] >= d[i - 1]);
}

TEST_CASE("decantiles on uniform distribution P50 equals median", "[decantiles]")
{
    // 100 values 1..100 → P50 = 50
    std::vector<int> v(100);
    std::iota(v.begin(), v.end(), 1);
    const auto d = decantiles(v);
    CHECK(d[4] == 50.0);   // P50
    CHECK(d[9] == 100.0);  // P100
    CHECK(d[0] == 10.0);   // P10
}
