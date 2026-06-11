# gapstats

Alignment gap & fragment-length QC statistics for BAM/CRAM long-read data.

[![CI](https://github.com/vinzenzmay/gapstats/actions/workflows/ci.yml/badge.svg)](https://github.com/vinzenzmay/gapstats/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/vinzenzmay/gapstats/branch/main/graph/badge.svg)](https://codecov.io/gh/vinzenzmay/gapstats)

## Install (pixi)

```sh
git clone https://github.com/vinzenzmay/gapstats
cd gapstats
pixi run install   # builds and installs gapstats into the pixi environment
```

Activate the environment and run:

```sh
pixi shell
gapstats --help
```

## Usage

```
gapstats --input sample.cram [sample2.cram ...] \
         --reference GRCh38.fa \
         --output stats.tsv \
         [--min-gap-size 6] [--threads 12]
```

## Output

TSV with 21 columns per input file:

| sample | gap_p10 … gap_p100 | flen_p10 … flen_p100 |
|--------|---------------------|----------------------|
| HG002  | gaps/Mb decantiles  | ref-span decantiles  |

## Build from source (without pixi)

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build
```

Requires: cmake ≥ 3.16, GCC ≥ 11, zlib, bz2, lzma, libcurl, libdeflate, autoconf/automake (for bundled htslib build).
