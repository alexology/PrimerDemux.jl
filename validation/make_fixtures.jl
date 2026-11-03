# Regenerates test/data/ (simulated reads) and test/reference/ (cutadapt 5.2
# outputs for every scenario of test/cases.jl).
#
# Usage: julia --project=. validation/make_fixtures.jl CUTADAPT [nreads] [seed]
# CUTADAPT is a command line that runs cutadapt, for instance on Windows with
# cutadapt in WSL: "wsl -e /home/user/miniforge3/envs/cutadapt-env/bin/cutadapt".
# Paths given to cutadapt are converted to /mnt/<drive>/... when the command
# starts with "wsl".
using PrimerDemux, Random, CodecZlibNG

const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "test", "cases.jl"))
include(joinpath(@__DIR__, "simulate.jl"))

cutadapt = split(ARGS[1])
nreads = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 150
seed = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 2026
use_wsl = first(cutadapt) == "wsl"
native(p) = use_wsl ? "/mnt/" * lowercase(p[1]) * replace(p[3:end], '\\' => '/') : p

data = joinpath(ROOT, "test", "data")
reference = joinpath(ROOT, "test", "reference")
mkpath(data)
mkpath(reference)
# Keep aligner.tsv/kmers.tsv (validation/aligner_reference.py); replace the scenarios.
for f in readdir(reference; join = true)
    isdir(f) && rm(f; recursive = true)
end

rng = MersenneTwister(seed)
write_fastq(joinpath(data, "reads.1.fastq"), [simulate_read(rng) for _ in 1:nreads], "1", rng)
write_fastq(joinpath(data, "reads.2.fastq"), [simulate_read(rng) for _ in 1:nreads], "2", rng)
open(joinpath(data, "barcodes.fasta"), "w") do io
    for (i, b) in enumerate(BARCODES)
        println(io, ">bc", i, " barcode ", i, "\n", b)
    end
end

# cutadapt's command-line arguments: a simple split that keeps "..." together.
splitargs(s) = [replace(m.match, "\"" => "") for m in eachmatch(r"\"[^\"]*\"|\S+", s)]

for case in CASES
    work = mktempdir()
    args = splitargs(replace(case.args, "BARCODES" => native(joinpath(data, "barcodes.fasta"))))
    outs = [native(joinpath(work, o)) for o in case_outputs(case.kind)]
    ins = [native(joinpath(data, i)) for i in case_inputs(case.kind)]
    outargs = length(outs) == 1 ? ["-o", outs[1]] : ["-o", outs[1], "-p", outs[2]]
    run(`$cutadapt -j 1 --quiet $args $outargs $ins`)
    dir = joinpath(reference, case.name)
    mkpath(dir)
    for f in readdir(work)
        open(GzipCompressorStream, joinpath(dir, f * ".gz"), "w") do io
            write(io, read(joinpath(work, f)))
        end
    end
    rm(work; recursive = true)
    println(case.name, ": ", length(readdir(dir)), " file(s)")
end
open(joinpath(reference, "VERSION"), "w") do io
    println(io, "cutadapt ", readchomp(`$cutadapt --version`))
    println(io, "reads: $nreads, seed: $seed")
end
