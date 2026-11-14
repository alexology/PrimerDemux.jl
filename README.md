# PrimerDemux.jl

[![CI](https://github.com/alexology/PrimerDemux.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/alexology/PrimerDemux.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://alexology.github.io/PrimerDemux.jl/dev/)

Primer and adapter removal and demultiplexing of FASTQ/FASTA reads in native
Julia, in the style of [cutadapt](https://github.com/marcelm/cutadapt) by
Marcel Martin. The alignment algorithm, its k-mer pre-filter, the rules that
choose between matches and the order of modifications and filters are
translated from cutadapt 5.2; the interface is a Julia API with keyword
arguments instead of a command line.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

With equivalent options the output is **identical to cutadapt 5.2** on the
41 scenarios of the test suite and on a real Oxford Nanopore sample (see
[Validation](#validation)), and it does not depend on the number of threads.

Julia ≥ 1.10 (tested with 1.12). Dependencies: FASTX.jl (FASTQ/FASTA
parsing), CodecZlibNG.jl (gzip input, zlib-ng), libdeflate_jll (gzip output),
TranscodingStreams.jl, StringViews.jl.

## Installation

The package is not registered. Install it from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/alexology/PrimerDemux.jl")
```

## Quick start

```julia
using PrimerDemux

# Linked, anchored primers on both strands; reads without them are dropped
# (cutadapt --revcomp -g ^FWD...REV_RC --discard-untrimmed)
primers = front_adapter("^AAGCTWMTGGGTTCATACCC...GTGCAATGGATTTAAGCTCCA")
stats = remove_primers("sample.fastq.gz", "sample.trimmed.fastq.gz";
                       adapters = primers, revcomp = true, discard_untrimmed = true)

# Demultiplexing by anchored 5' barcodes from a FASTA file
# (cutadapt -g ^file:barcodes.fasta -o "{name}.fastq.gz")
demultiplex("pool.fastq.gz", "demux/{name}.fastq.gz";
            adapters = adapters_from_fasta("barcodes.fasta"; anchored = true))

# Paired-end: primers in R1 and R2, pairs shorter than 100 bp dropped
remove_primers(("R1.fastq.gz", "R2.fastq.gz"), ("out_R1.fastq.gz", "out_R2.fastq.gz");
               adapters = front_adapter("^GGWACWGGWTGAACWGTWTAYCCYCC"),
               adapters2 = front_adapter("^TAAACTTCAGGGTGACCAAAAAATCA"),
               discard_untrimmed = true, min_length = 100)
```

Both functions return a `TrimStats` summary (printed like cutadapt's report).

## Describing adapters and primers

| Julia | cutadapt | Meaning |
|---|---|---|
| `front_adapter("ACGT")` | `-g ACGT` | 5' adapter anywhere, also partial at the 5' end; it and everything before it are removed |
| `front_adapter("^ACGT")` | `-g ^ACGT` | anchored 5': complete, at the start of the read |
| `front_adapter("XACGT")` | `-g XACGT` | non-internal 5': may be partial, must start the read |
| `back_adapter("ACGT")` | `-a ACGT` | 3' adapter anywhere, also partial at the 3' end; it and everything after it are removed |
| `back_adapter("ACGT\$")` | `-a 'ACGT$'` | anchored 3' |
| `back_adapter("ACGTX")` | `-a ACGTX` | non-internal 3' |
| `front_adapter("FWD...REV")` | `-g FWD...REV` | linked: 5' part, then 3' part in the rest; both required |
| `back_adapter("FWD...REV")` | `-a FWD...REV` | linked; only anchored parts (`^FWD`, `REV$`) required |
| `adapters_from_fasta(path)` | `-g file:path` | one adapter per record, named after the record; `anchored = true` for `^file:`, `type = :back` for `-a` |

The specification syntax is cutadapt's: `"name=SEQUENCE"` sets the name used
in reports and for demultiplexing (default: the specification itself),
`";e=0.15;o=5;noindels"` overrides the search parameters of one adapter,
`";required"`/`";optional"` apply to the parts of linked adapters, `"A{10}"`
repeats a base. IUPAC codes in adapters are wildcards.

Keywords of `front_adapter`, `back_adapter` and `adapters_from_fasta`:

| Keyword | Default | cutadapt |
|---|---|---|
| `max_errors` | `0.1` | `-e`: error rate; from 1 upwards, a number of errors |
| `min_overlap` | `3` | `-O` |
| `indels` | `true` | `--no-indels` |
| `adapter_wildcards` | `true` | `-N` / `--no-match-adapter-wildcards` |
| `read_wildcards` | `false` | `--match-read-wildcards` |
| `name` | | `name=` |

## Processing files

`remove_primers(input, output; ...)` writes one output (two for pairs);
`demultiplex(input, "dir/{name}.fastq.gz"; ...)` writes one file per adapter
name, plus `unknown` for reads without a match (or `untrimmed_output`, or
nothing with `discard_untrimmed = true`). For pairs, `{name}` follows R1, while
`{name1}` and `{name2}` in both paths route by the combination of R1 and R2
matches.

Inputs may be FASTQ or FASTA, plain or gzip. Outputs are gzip when the name
ends in `.gz`, FASTA when it ends in `.fasta`, `.fa`, `.fna` or `.fas`.
Outputs are written to temporary files in the same directory and renamed into
place at the end; an error leaves no partial output, and replaced files are
restored if putting the outputs in place fails (a crash during those final
renames is not covered; if the file system prevents a restore, the previous
content is kept in a temporary file whose name is logged). Existing files are
kept unless `overwrite = true`, also when they appear while the run is in
progress. Outputs that would be the same file (also through links, or names differing
only in case on Windows and macOS) or an input are refused. As in cutadapt,
paired inputs must have matching read IDs, up to a final 1, 2 or 3.

| Keyword | Default | cutadapt |
|---|---|---|
| `adapters`, `adapters2` | | `-a`/`-g`, `-A`/`-G` |
| `times` | `1` | `-n` |
| `action` | `:trim` | `--action` (`:trim`, `:retain`, `:mask`, `:lowercase`, `:none`) |
| `revcomp` | `false` | `--revcomp` (pairs: R1 and R2 swapped) |
| `rc_suffix` | `" rc"` | suffix that `--revcomp` adds to the names of reversed reads |
| `pair_adapters` | `false` | `--pair-adapters` |
| `cut`, `cut2` | `0` | `-u`, `-U` (an integer or up to two in a vector) |
| `quality_cutoff`, `quality_cutoff2` | `nothing` | `-q`, `-Q` (`20` or `(5, 20)`; without `quality_cutoff2`, R2 uses `quality_cutoff`, as cutadapt) |
| `quality_base` | `33` | `--quality-base` |
| `min_length`, `max_length` | `nothing` | `-m`, `-M` (pairs: `(R1, R2)`, `nothing` for no limit) |
| `max_n` | `nothing` | `--max-n` (count, or fraction below 1) |
| `max_expected_errors` | `nothing` | `--max-ee` |
| `discard_trimmed` | `false` | `--discard-trimmed` |
| `discard_untrimmed` | `false` | `--discard-untrimmed` |
| `untrimmed_output` | `nothing` | `--untrimmed-output` (`--untrimmed-paired-output`) |
| `pair_filter` | `:any` | `--pair-filter` (`:any`, `:both`, `:first`) |
| `threads` | `Threads.nthreads()` | `-j` |
| `batch_size` | `4096` | reads per batch |
| `compression_level` | `1` | `--compression-level` (libdeflate levels 0–12) |
| `overwrite` | `false` | |

Order of the steps, as in cutadapt: unconditional cuts, quality trimming,
adapters (with the reverse complement if requested), then the filters too
short, too long, too many N, too many expected errors, and finally trimmed or
untrimmed reads.

For single sequences, `find_adapter(adapters, seq)` returns the best
`AdapterMatch` and `trim_sequence(adapters, seq; times, action, revcomp)`
the trimmed sequence.

## Validation

- **Aligner and k-mer filter:** 3,000 random cases each (all placement flags,
  wildcards, indel costs, error rates) give the same results as cutadapt's
  `Aligner.locate` and `KmerFinder` (`test/reference/*.tsv`, generated by
  `validation/aligner_reference.py`).
- **Whole runs:** 41 scenarios on simulated reads (`test/cases.jl`: every
  adapter type, linked adapters, `--revcomp`, error rates, wildcards, actions,
  quality trimming and filters, demultiplexing single-end, paired-end and
  combinatorial) produce byte-identical files to cutadapt 5.2
  (`test/reference/`, generated by `validation/make_fixtures.jl`). During
  development the same comparison passed on three further seeds with 3,000
  reads per file.
- **Real data:** an Oxford Nanopore COI sample (502,718 reads, 1.27 Gbp),
  `revcomp`, linked anchored primers, `discard_untrimmed`: the 375,943 reads
  written (178,814 reverse-complemented) are identical to cutadapt's, line by
  line.

Speed is a goal, not a guarantee. On that sample (single measurements on a
12-thread Windows laptop): PrimerDemux about 25 s with 4 threads, 52 s with
1, 17 s with an uncompressed output; cutadapt 5.2 (WSL, `-j 4`) about 25 s.
Gzip inputs are decompressed (zlib-ng) by a task of their own while FASTX
parses, output blocks are compressed with libdeflate in the worker tasks,
and read buffers are reused between batches. With 3 or more threads one is
left to the reading tasks. Output size 0.98 GB against cutadapt's 1.02 GB.

## Differences from cutadapt

- Not available: `-b` (anywhere) and `;rightmost` adapters, `--action=crop`,
  NextSeq/poly-A trimming, `--trim-n`, `--length`, renaming templates, info
  and JSON report files, `--times` with `pair_adapters`, colorspace, bzip2/xz.
- The index cutadapt builds for several anchored adapters (an optimisation
  with its own rules for ambiguous sequences) is not used: results equal
  cutadapt with `--no-index`. With many similar barcodes cutadapt's default
  may assign some reads differently.
- Adapter names default to the specification, not to "1", "2", ...
- When a 5' search window of the k-mer pre-filter extends past the end of a
  short read, cutadapt reads past the string (an undefined result);
  PrimerDemux then always runs the aligner. No difference was observed.
- Gzip outputs consist of one gzip member per batch (like parallel gzip
  tools); the decompressed content is the same.

## Citation and license

Please cite cutadapt: Martin M. (2011) Cutadapt removes adapter sequences from
high-throughput sequencing reads. *EMBnet.journal* 17(1): 10–12.
doi:10.14806/ej.17.1.200. PrimerDemux.jl is released under the MIT license
([LICENSE](LICENSE)). cutadapt is also MIT-licensed; its copyright notice is
in [NOTICE](NOTICE).
