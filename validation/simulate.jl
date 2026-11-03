# Simulated amplicon reads for comparisons with cutadapt: optional barcode,
# IUPAC forward primer, insert, reverse primer and Illumina adapter, each
# possibly mutated (substitutions, insertions, deletions) or truncated; random
# N bases, lowercase stretches and reverse-complemented reads.

const BARCODES = ["ACGTACGT", "TTGCAAGC", "GGATCCAA", "CATGCATG", "AACCGGTT"]
const SIM_FWD = "AAGCTWMTGGGTTCATACCC"
const SIM_REV_RC = "GTGCAATGGATTTAAGCTCCA"
const SIM_ADAPT = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

randseq(rng, n) = String(rand(rng, collect("ACGT"), n))

function concretize(rng, p)
    codes = Dict('W' => "AT", 'M' => "AC", 'R' => "AG", 'Y' => "CT", 'N' => "ACGT", 'S' => "CG", 'K' => "GT")
    return String([haskey(codes, c) ? rand(rng, codes[c]) : c for c in p])
end

function mutate(rng, s; rate = 0.06)
    out = Char[]
    for c in s
        r = rand(rng)
        if r < rate / 3
            push!(out, rand(rng, "ACGT"))
        elseif r < 2rate / 3
            continue
        elseif r < rate
            push!(out, c, rand(rng, "ACGT"))
        else
            push!(out, c)
        end
    end
    return String(out)
end

function truncate_randomly(rng, s)
    isempty(s) && return s
    r = rand(rng)
    r < 0.15 && return s[rand(rng, 1:length(s)):end]
    r < 0.3 && return s[1:rand(rng, 1:length(s))]
    return s
end

function simulate_read(rng; insert = 20:200)
    parts = String[]
    rand(rng) < 0.2 && push!(parts, randseq(rng, rand(rng, 1:15)))
    if rand(rng) < 0.7
        rand(rng) < 0.5 && push!(parts, mutate(rng, rand(rng, BARCODES)))
        push!(parts, truncate_randomly(rng, mutate(rng, concretize(rng, SIM_FWD))))
    end
    push!(parts, randseq(rng, rand(rng, insert)))
    rand(rng) < 0.6 && push!(parts, truncate_randomly(rng, mutate(rng, SIM_REV_RC)))
    rand(rng) < 0.4 && push!(parts, truncate_randomly(rng, mutate(rng, SIM_ADAPT)))
    rand(rng) < 0.2 && push!(parts, randseq(rng, rand(rng, 1:20)))
    s = collect(join(parts))
    for i in eachindex(s)
        rand(rng) < 0.005 && (s[i] = 'N')
    end
    rand(rng) < 0.5 && (s = collect(PrimerDemux.reverse_complement(String(s))))
    if rand(rng) < 0.05
        a = rand(rng, eachindex(s))
        b = min(length(s), a + 10)
        s[a:b] .= lowercase.(s[a:b])
    end
    return String(s)
end

function write_fastq(path, reads, tag, rng)
    open(path, "w") do io
        for (i, s) in enumerate(reads)
            q = String([Char(33 + (rand(rng) < 0.1 ? rand(rng, 2:15) : rand(rng, 20:40))) for _ in 1:length(s)])
            print(io, "@read", i, "_", tag, " extra=", i, "\n", s, "\n+\n", q, "\n")
        end
    end
end
