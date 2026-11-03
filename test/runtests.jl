using PrimerDemux
using Test
using CodecZlibNG: GzipDecompressorStream

const DATA = joinpath(@__DIR__, "data")
const REFERENCE = joinpath(@__DIR__, "reference")

include("cases.jl")

read_gz(path) = open(io -> read(io, String), GzipDecompressorStream, path)

@testset "PrimerDemux" begin
    include("test_units.jl")

    @testset "cutadapt 5.2 outputs: $(case.name)" for case in CASES
        work = mktempdir()
        inputs = joinpath.(DATA, case_inputs(case.kind))
        outputs = joinpath.(work, case_outputs(case.kind))
        case.run(length(inputs) == 1 ? inputs[1] : Tuple(inputs),
                 length(outputs) == 1 ? outputs[1] : Tuple(outputs), DATA)
        expected = sort(readdir(joinpath(REFERENCE, case.name)))
        @test sort(readdir(work)) == [chop(f; tail = 3) for f in expected]
        for f in expected
            got = joinpath(work, chop(f; tail = 3))
            isfile(got) || continue
            @test read(got, String) == read_gz(joinpath(REFERENCE, case.name, f))
        end
        rm(work; recursive = true)
    end

    include("test_io.jl")
    include("test_review.jl")
end
