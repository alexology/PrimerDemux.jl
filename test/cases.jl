# Scenarios compared with cutadapt 5.2. Shared by the tests and by
# validation/make_fixtures.jl, which runs cutadapt to create test/reference/.
#
# kind: :se (one output), :se_demux ({name}), :pe (two outputs),
# :pe_demux ({name}), :pe_combi ({name1}-{name2}).
# `run(input, output, data)` calls PrimerDemux; `data` is the directory with
# barcodes.fasta. In `args`, BARCODES stands for that FASTA file.

const FWD = "AAGCTWMTGGGTTCATACCC"        # IUPAC primer
const REV_RC = "GTGCAATGGATTTAAGCTCCA"
const ADAPT = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
const REV = PrimerDemux.reverse_complement(REV_RC)

bc(data; kw...) = adapters_from_fasta(joinpath(data, "barcodes.fasta"); kw...)

const CASES = [
    (name = "front", kind = :se, args = "-g $FWD",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter(FWD))),
    (name = "back", kind = :se, args = "-a $REV_RC",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter(REV_RC))),
    (name = "prefix", kind = :se, args = "-g ^$FWD",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter("^" * FWD))),
    (name = "suffix", kind = :se, args = "-a $(REV_RC)\$",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter(REV_RC * "\$"))),
    (name = "prefix_noindels", kind = :se, args = "-g ^$FWD --no-indels",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter("^" * FWD; indels = false))),
    (name = "suffix_noindels", kind = :se, args = "-a $(REV_RC)\$ --no-indels -e 0.15",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter(REV_RC * "\$"; indels = false, max_errors = 0.15))),
    (name = "noninternal", kind = :se, args = "-g X$FWD -a $(REV_RC)X",
     run = (i, o, d) -> remove_primers(i, o; adapters = [front_adapter("X" * FWD), back_adapter(REV_RC * "X")])),
    (name = "linked_g", kind = :se, args = "-g ^$FWD...$(REV_RC)\$",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter("^$FWD...$(REV_RC)\$"))),
    (name = "linked_a", kind = :se, args = "-a ^$FWD...$REV_RC",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter("^$FWD...$REV_RC"))),
    (name = "linked_free", kind = :se, args = "-g $FWD...$REV_RC",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter("$FWD...$REV_RC"))),
    (name = "linked_a_free", kind = :se, args = "-a $FWD...$REV_RC",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter("$FWD...$REV_RC"))),
    (name = "linked_revcomp", kind = :se, args = "--revcomp -e 0.1 -g ^$FWD...$REV_RC --discard-untrimmed",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter("^$FWD...$REV_RC"), revcomp = true, discard_untrimmed = true)),
    (name = "revcomp_back", kind = :se, args = "--revcomp -a $ADAPT -O 5",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter(ADAPT; min_overlap = 5), revcomp = true)),
    (name = "e2_o5", kind = :se, args = "-e 0.2 -O 5 -g $FWD -a $ADAPT",
     run = (i, o, d) -> remove_primers(i, o; adapters = [front_adapter(FWD; max_errors = 0.2, min_overlap = 5), back_adapter(ADAPT; max_errors = 0.2, min_overlap = 5)])),
    (name = "abs_errors", kind = :se, args = "-e 2 -a $ADAPT",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter(ADAPT; max_errors = 2))),
    (name = "noindels", kind = :se, args = "--no-indels -a $ADAPT -g $FWD",
     run = (i, o, d) -> remove_primers(i, o; adapters = [back_adapter(ADAPT; indels = false), front_adapter(FWD; indels = false)])),
    (name = "read_wildcards", kind = :se, args = "--match-read-wildcards -a $ADAPT -g $FWD",
     run = (i, o, d) -> remove_primers(i, o; adapters = [back_adapter(ADAPT; read_wildcards = true), front_adapter(FWD; read_wildcards = true)])),
    (name = "no_adapter_wildcards", kind = :se, args = "-N -g $FWD",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter(FWD; adapter_wildcards = false))),
    (name = "times2", kind = :se, args = "-n 2 -a $ADAPT -a $REV_RC",
     run = (i, o, d) -> remove_primers(i, o; adapters = [back_adapter(ADAPT), back_adapter(REV_RC)], times = 2)),
    (name = "mask", kind = :se, args = "--action=mask -n 2 -g $FWD -a $ADAPT",
     run = (i, o, d) -> remove_primers(i, o; adapters = [front_adapter(FWD), back_adapter(ADAPT)], times = 2, action = :mask)),
    (name = "lowercase", kind = :se, args = "--action=lowercase -g $FWD...$REV_RC",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter("$FWD...$REV_RC"), action = :lowercase)),
    (name = "retain", kind = :se, args = "--action=retain -g $FWD -a $REV_RC",
     run = (i, o, d) -> remove_primers(i, o; adapters = [front_adapter(FWD), back_adapter(REV_RC)], action = :retain)),
    (name = "retain_linked", kind = :se, args = "--action=retain -a ^$FWD...$REV_RC",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter("^$FWD...$REV_RC"), action = :retain)),
    (name = "none_discard_trimmed", kind = :se, args = "--action=none -g $FWD --discard-trimmed",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter(FWD), action = :none, discard_trimmed = true)),
    (name = "filters", kind = :se, args = "-u 3 -u -2 -q 15,20 -a $ADAPT -m 60 -M 250 --max-n 2 --max-ee 3",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter(ADAPT), cut = [3, -2], quality_cutoff = (15, 20),
                                       min_length = 60, max_length = 250, max_n = 2, max_expected_errors = 3)),
    (name = "max_n_fraction", kind = :se, args = "-q 20 --max-n 0.01 -g $FWD",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter(FWD), quality_cutoff = 20, max_n = 0.01)),
    (name = "spec_parameters", kind = :se, args = "-a \"A{3}GATCGGAAGAGC;e=0.15;o=6\" -g \"$FWD;noindels\"",
     run = (i, o, d) -> remove_primers(i, o; adapters = [back_adapter("A{3}GATCGGAAGAGC;e=0.15;o=6"), front_adapter("$FWD;noindels")])),
    (name = "linked_optional", kind = :se, args = "-g \"$FWD;optional...$REV_RC;e=0.2\"",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter("$FWD;optional...$REV_RC;e=0.2"))),
    (name = "demux_barcodes", kind = :se_demux, args = "--no-index -g ^file:BARCODES -e 0.15",
     run = (i, o, d) -> demultiplex(i, o; adapters = bc(d; anchored = true, max_errors = 0.15))),
    (name = "demux_free_discard", kind = :se_demux, args = "-g file:BARCODES --discard-untrimmed",
     run = (i, o, d) -> demultiplex(i, o; adapters = bc(d), discard_untrimmed = true)),
    (name = "demux_primers_revcomp", kind = :se_demux, args = "--revcomp -g fwd=^$FWD -g rev=^$REV",
     run = (i, o, d) -> demultiplex(i, o; adapters = [front_adapter("fwd=^$FWD"), front_adapter("rev=^$REV")], revcomp = true)),
    (name = "pe_basic", kind = :pe, args = "-g $FWD -G $FWD -m 30",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter(FWD), adapters2 = front_adapter(FWD), min_length = 30)),
    (name = "pe_filter_both", kind = :pe, args = "-a $ADAPT -A $ADAPT -m 100:80 --pair-filter=both",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter(ADAPT), adapters2 = back_adapter(ADAPT), min_length = (100, 80), pair_filter = :both)),
    (name = "pe_untrimmed_r1_only", kind = :pe, args = "-g ^$FWD --discard-untrimmed",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter("^" * FWD), discard_untrimmed = true)),
    (name = "pe_untrimmed_any", kind = :pe, args = "-g $FWD -G $FWD --discard-untrimmed",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter(FWD), adapters2 = front_adapter(FWD), discard_untrimmed = true)),
    (name = "pe_revcomp", kind = :pe, args = "--revcomp -g ^$FWD -G ^$REV",
     run = (i, o, d) -> remove_primers(i, o; adapters = front_adapter("^" * FWD), adapters2 = front_adapter("^" * REV), revcomp = true)),
    (name = "pe_pair_adapters", kind = :pe, args = "--pair-adapters -g a=$FWD -g b=$ADAPT -G a=$FWD -G b=$REV_RC",
     run = (i, o, d) -> remove_primers(i, o; adapters = [front_adapter("a=$FWD"), front_adapter("b=$ADAPT")],
                                       adapters2 = [front_adapter("a=$FWD"), front_adapter("b=$REV_RC")], pair_adapters = true)),
    (name = "pe_quality", kind = :pe, args = "-q 20 -Q 10 -u 2 -U -3 -a $ADAPT -A $ADAPT --max-n 1 --pair-filter=first",
     run = (i, o, d) -> remove_primers(i, o; adapters = back_adapter(ADAPT), adapters2 = back_adapter(ADAPT), quality_cutoff = 20,
                                       quality_cutoff2 = 10, cut = 2, cut2 = -3, max_n = 1, pair_filter = :first)),
    (name = "pe_demux", kind = :pe_demux, args = "--no-index -g ^file:BARCODES -G $FWD",
     run = (i, o, d) -> demultiplex(i, o; adapters = bc(d; anchored = true), adapters2 = front_adapter(FWD))),
    (name = "pe_combinatorial", kind = :pe_combi, args = "--no-index -g ^file:BARCODES -G ^file:BARCODES",
     run = (i, o, d) -> demultiplex(i, o; adapters = bc(d; anchored = true), adapters2 = bc(d; anchored = true))),
    (name = "pe_combinatorial_discard", kind = :pe_combi, args = "--no-index -g file:BARCODES -G ^file:BARCODES --discard-untrimmed",
     run = (i, o, d) -> demultiplex(i, o; adapters = bc(d), adapters2 = bc(d; anchored = true), discard_untrimmed = true)),
]

# Output names used for cutadapt and for PrimerDemux.
function case_outputs(kind::Symbol)
    kind == :se && return ["out.fastq"]
    kind == :se_demux && return ["{name}.fastq"]
    kind == :pe && return ["out.1.fastq", "out.2.fastq"]
    kind == :pe_demux && return ["{name}.1.fastq", "{name}.2.fastq"]
    return ["{name1}-{name2}.1.fastq", "{name1}-{name2}.2.fastq"]
end
case_inputs(kind::Symbol) = kind in (:se, :se_demux) ? ["reads.1.fastq"] : ["reads.1.fastq", "reads.2.fastq"]
