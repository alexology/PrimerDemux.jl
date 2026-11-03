# Functions for single sequences, mainly for exploration and testing.

_seqbytes(seq::AbstractString) = Vector{UInt8}(codeunits(String(seq)))

"""
    find_adapter(adapters, seq) -> Union{Nothing,AdapterMatch}

Best match of `adapters` (one adapter or a vector) in the sequence `seq`, with
the same rules as the file functions (highest score, then fewest errors, then
first adapter). `m.adapter` is the index of the matching adapter;
`PrimerDemux.score(m)` and `PrimerDemux.errors(m)` give score and errors.
"""
find_adapter(adapters, seq::AbstractString) =
    best_match(_adapters(adapters), _seqbytes(seq), AlignWorkspace())

"""
    trim_sequence(adapters, seq; times=1, action=:trim, revcomp=false)
        -> (sequence, matches, reverse_complemented)

Adapter removal on one sequence, as done for each read by
[`remove_primers`](@ref). With `revcomp=true` the returned sequence is the
reverse complement when that orientation matched better.
"""
function trim_sequence(adapters, seq::AbstractString; times::Integer = 1,
                       action::Symbol = :trim, revcomp::Bool = false)
    check_action(action, times)
    cutter = AdapterCutter(_adapters(adapters), times, action)
    ws = AlignWorkspace()
    read = Read(UInt8[], _seqbytes(seq), UInt8[])
    fwd, fm = match_and_trim(cutter, read, false, ws)
    if revcomp
        rev, rm = match_and_trim(cutter, reverse_complement(read, false), false, ws)
        if sum(score, rm; init = 0) > sum(score, fm; init = 0)
            return String(rev.seq), rm, true
        end
    end
    return String(fwd.seq), fm, false
end

function Base.show(io::IO, m::AdapterMatch)
    print(io, "AdapterMatch(adapter=", m.adapter)
    m.front === nothing || print(io, ", 5' read[", m.front.rstart, ":", m.front.rstop, "]")
    m.back === nothing || print(io, ", 3' read[", m.back.rstart, ":", m.back.rstop, "]")
    print(io, ", score=", score(m), ", errors=", errors(m), ")")
end
