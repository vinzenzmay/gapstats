#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "core.hpp"

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "CLI/CLI.hpp"

#ifndef GAPSTATS_VERSION
#define GAPSTATS_VERSION "dev"
#endif

// ─────────────────────────────────────────────────────────────────────────────
// Types
// ─────────────────────────────────────────────────────────────────────────────

struct FileStats {
    std::string          sample;
    std::vector<double>  gaps_per_mb;    // one entry per mapped read
    std::vector<int64_t> frag_lengths;   // reference span per mapped read
};

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
// Per-file processing
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Process a single BAM/CRAM file.
 *
 * Strategy: single sequential reader with htslib's built-in BGZF thread pool
 * (hts_tpool) for parallel decompression. CIGAR arithmetic is O(µs/read) and
 * is done inline in the main thread — the bottleneck is I/O and decompression,
 * not CIGAR parsing.
 *
 * Why not region-parallel N-file-handles?  For large (50 GB) CRAM files on a
 * spinning / USB drive, N concurrent readers each seeking to different byte
 * offsets creates random I/O that destroys throughput.  A single sequential
 * pass with N decompressor threads in the background is optimal for
 * I/O-bound workloads and gives the same CPU utilisation.
 */
static FileStats process_file(const std::string& path,
                               const std::string& reference,
                               int                min_gap_size,
                               int                n_threads)
{
    FileStats stats;
    stats.sample = gapstats::stem(path);

    htsFile* fp = hts_open(path.c_str(), "r");
    if (!fp) {
        std::cerr << "[gapstats] ERROR: cannot open " << path << "\n";
        return stats;
    }

    // Set reference for CRAM decoding
    if (!reference.empty()) {
        if (hts_set_fai_filename(fp, reference.c_str()) != 0) {
            std::cerr << "[gapstats] ERROR: cannot set reference for " << path << "\n";
            hts_close(fp);
            return stats;
        }
    }

    // ── Parallel BGZF decompression via htslib thread pool ──────────────────
    // This runs n_threads decompressor goroutines in the background while the
    // main thread calls sam_read1().  Delivers near-linear speedup up to the
    // I/O bandwidth limit.
    htsThreadPool pool = {nullptr, 0};
    if (n_threads > 1) {
        pool.pool = hts_tpool_init(n_threads);
        if (pool.pool)
            hts_set_thread_pool(fp, &pool);
    }

    sam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) {
        std::cerr << "[gapstats] ERROR: cannot read header from " << path << "\n";
        if (pool.pool) hts_tpool_destroy(pool.pool);
        hts_close(fp);
        return stats;
    }

    bam1_t* b = bam_init1();

    // ── Sequential scan ──────────────────────────────────────────────────────
    int ret;
    while ((ret = sam_read1(fp, hdr, b)) >= 0) {
        // Skip unmapped reads
        if (b->core.flag & BAM_FUNMAP) continue;

        const uint32_t  n_cigar = b->core.n_cigar;
        const uint32_t* cigar   = bam_get_cigar(b);

        // Reference span = genomic footprint of this alignment
        const int64_t ref_span =
            static_cast<int64_t>(bam_cigar2rlen(static_cast<int>(n_cigar), cigar));
        if (ref_span <= 0) continue;

        const int gap_count = gapstats::count_gaps(cigar, n_cigar, min_gap_size);

        stats.gaps_per_mb.push_back(
            static_cast<double>(gap_count) / (static_cast<double>(ref_span) / 1e6));
        stats.frag_lengths.push_back(ref_span);
    }

    if (ret < -1) {
        std::cerr << "[gapstats] WARNING: error reading records from " << path << "\n";
    }

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    if (pool.pool) hts_tpool_destroy(pool.pool);
    hts_close(fp);

    return stats;
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char** argv)
{
    CLI::App app{"gapstats — alignment gap & fragment-length QC statistics"};
    app.set_version_flag("-V,--version", GAPSTATS_VERSION);

    std::vector<std::string> inputs;
    std::string              output;
    std::string              reference;
    int                      min_gap_size = 6;
    int                      threads      = static_cast<int>(
                                                std::max(1u,
                                                    std::thread::hardware_concurrency()));

    app.add_option("-i,--input",        inputs,       "BAM/CRAM input file(s)")
        ->required()->check(CLI::ExistingFile);
    app.add_option("-o,--output",       output,       "Output TSV file")->required();
    app.add_option("-r,--reference",    reference,    "Reference FASTA (required for CRAM)")
        ->check(CLI::ExistingFile);
    app.add_option("-g,--min-gap-size", min_gap_size, "Min insertion/deletion size to count as gap [6]");
    app.add_option("-t,--threads",      threads,      "Worker threads [hardware_concurrency]");

    CLI11_PARSE(app, argc, argv);

    // Validate: CRAM inputs require a reference
    for (const auto& f : inputs) {
        const auto ext = std::filesystem::path(f).extension().string();
        if ((ext == ".cram") && reference.empty()) {
            std::cerr << "[gapstats] ERROR: --reference is required for CRAM input: "
                      << f << "\n";
            return 1;
        }
    }

    // ── Open output TSV ───────────────────────────────────────────────────────
    std::ofstream tsv(output);
    if (!tsv) {
        std::cerr << "[gapstats] ERROR: cannot open output file: " << output << "\n";
        return 1;
    }

    // Header
    tsv << "sample";
    for (int p = 10; p <= 100; p += 10) tsv << "\tgap_p" << p;
    for (int p = 10; p <= 100; p += 10) tsv << "\tflen_p" << p;
    tsv << "\n";

    // ── Process files sequentially; within each file use N threads ────────────
    for (const auto& f : inputs) {
        std::cerr << "[gapstats] Processing " << f << " ...\n";

        FileStats stats = process_file(f, reference, min_gap_size, threads);

        if (stats.gaps_per_mb.empty() || stats.frag_lengths.empty()) {
            std::cerr << "[gapstats] WARNING: no mapped reads found in " << f << "\n";
        }

        auto gap_dec  = gapstats::decantiles(stats.gaps_per_mb);
        auto flen_dec = gapstats::decantiles(stats.frag_lengths);

        tsv << stats.sample;
        for (auto v : gap_dec)  tsv << "\t" << v;
        for (auto v : flen_dec) tsv << "\t" << static_cast<int64_t>(v);
        tsv << "\n";
        tsv.flush();  // persist row immediately in case of interruption

        std::cerr << "[gapstats]   reads processed: " << stats.frag_lengths.size() << "\n";
    }

    std::cerr << "[gapstats] Done. Output written to " << output << "\n";
    return 0;
}
