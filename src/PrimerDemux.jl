"""
    PrimerDemux

Primer and adapter removal and demultiplexing of FASTQ/FASTA reads in the
style of [cutadapt](https://github.com/marcelm/cutadapt), with a Julia API.
The alignment algorithm, its k-mer pre-filter and the rules that choose
between matches are translated from cutadapt 5.2.

Main entry points: [`remove_primers`](@ref) and [`demultiplex`](@ref) for
files; [`front_adapter`](@ref), [`back_adapter`](@ref) and
[`adapters_from_fasta`](@ref) to describe what to search; [`find_adapter`](@ref)
and [`trim_sequence`](@ref) for single sequences.
"""
module PrimerDemux

using FASTX: FASTX, FASTQ, FASTA, description, identifier, sequence, quality
using StringViews: StringView
using TranscodingStreams: TranscodingStreams, NoopStream
using CodecZlibNG: GzipDecompressorStream
import libdeflate_jll

export front_adapter, back_adapter, adapters_from_fasta,
       Adapter, SingleAdapter, LinkedAdapter, Placement,
       FRONT, BACK, PREFIX, SUFFIX, FRONT_NOT_INTERNAL, BACK_NOT_INTERNAL,
       AdapterMatch, find_adapter, trim_sequence,
       remove_primers, demultiplex, TrimStats,
       reverse_complement, quality_trim_index

include("encoding.jl")
include("align.jl")
include("kmers.jl")
include("adapters.jl")
include("io.jl")
include("trimming.jl")
include("pipeline.jl")
include("api.jl")

end # module
