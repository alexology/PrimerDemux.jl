# Read modifications applied before filtering, in cutadapt's order:
# unconditional cuts, quality trimming, adapter/primer removal.

const ACTIONS = (:trim, :retain, :mask, :lowercase, :none)

struct AdapterCutter
    adapters::Vector{Adapter}
    times::Int
    action::Symbol
end

# Checks shared by the file functions and trim_sequence.
function check_action(action::Symbol, times::Integer)
    action in ACTIONS || throw(ArgumentError("action must be one of $(ACTIONS)"))
    times >= 1 || throw(ArgumentError("times must be at least 1"))
    action == :retain && times > 1 && throw(ArgumentError("action :retain cannot be combined with times > 1"))
    return nothing
end

# Python's str.upper/lower on ASCII.
_upper(seq::Vector{UInt8}) = [UPPER_TABLE[b+1] for b in seq]
_lower(b::UInt8) = 0x41 <= b <= 0x5a ? b + 0x20 : b

"""
Search the best adapter up to `times` times, removing it each time; then apply
the action to the original read. Returns (read, matches).
"""
function match_and_trim(c::AdapterCutter, read::Read, has_qual::Bool, ws::AlignWorkspace)
    matches = AdapterMatch[]
    if c.action == :lowercase
        # cutadapt uppercases the read in this mode, matched or not
        read = Read(read.name, _upper(read.seq), read.qual)
    end
    trimmed = read
    for _ in 1:c.times
        m = best_match(c.adapters, trimmed.seq, ws)
        m === nothing && break
        push!(matches, m)
        start, stop = remainder_interval(m, length(trimmed))
        trimmed = slice(trimmed, start, stop, has_qual)
    end
    isempty(matches) && return trimmed, matches
    c.action == :trim && return trimmed, matches
    return apply_action(read, matches, c.action, has_qual), matches
end

# Applies `action` to the original read given the matches of all rounds.
function apply_action(read::Read, matches::Vector{AdapterMatch}, action::Symbol, has_qual::Bool)
    action == :lowercase && (read = Read(read.name, _upper(read.seq), read.qual))
    action == :none && return read
    if action == :retain
        start, stop = retained_adapter_interval(matches[end], length(read))
        return slice(read, start, stop, has_qual)
    end
    start, stop = remainder(matches)
    action == :trim && return slice(read, start, stop, has_qual)
    # :mask and :lowercase act outside the part that trimming would keep.
    seq = copy(read.seq)
    for i in eachindex(seq)
        if i <= start || i > stop
            seq[i] = action == :mask ? UInt8('N') : _lower(seq[i])
        end
    end
    return Read(read.name, seq, read.qual)
end

# Interval of the untrimmed part relative to the original read.
function remainder(matches::Vector{AdapterMatch})
    start = 0
    match_start = match_stop = 0
    for m in matches
        # length of the sequence the match was searched in
        seqlen = m.front !== nothing ? m.front.seqlen : m.back.seqlen
        match_start, match_stop = remainder_interval(m, seqlen)
        start += match_start
    end
    return start, start + (match_stop - match_start)
end

"""
    quality_trim_index(qualities, cutoff_front, cutoff_back; base=33) -> (start, stop)

BWA-style quality trimming as in cutadapt `-q front,back`: 0-based half-open
interval of the part to keep. `base` can also be given as a fourth positional
argument.
"""
quality_trim_index(qual::AbstractVector{UInt8}, cutoff_front::Integer, cutoff_back::Integer; base::Integer = 33) =
    quality_trim_index(qual, cutoff_front, cutoff_back, base)

function quality_trim_index(qual::AbstractVector{UInt8}, cutoff_front::Integer, cutoff_back::Integer, base::Integer)
    n = length(qual)
    start = 0
    stop = n
    s = 0
    max_qual = 0
    @inbounds for i in 1:n
        s += cutoff_front - (Int(qual[i]) - base)
        s < 0 && break
        if s > max_qual
            max_qual = s
            start = i
        end
    end
    max_qual = 0
    s = 0
    @inbounds for i in n:-1:1
        s += cutoff_back - (Int(qual[i]) - base)
        s < 0 && break
        if s > max_qual
            max_qual = s
            stop = i - 1
        end
    end
    start >= stop && return 0, 0
    return start, stop
end

function expected_errors(qual::AbstractVector{UInt8}, base::Int = 33)
    e = 0.0
    @inbounds for q in qual
        e += exp10(-(Int(q) - base) / 10)
    end
    return e
end

# Unconditional cuts and quality trimming of one mate.
struct MateSteps
    cuts::Vector{Int}
    quality::Union{Nothing,Tuple{Int,Int}}
    quality_base::Int
end

function apply_steps(s::MateSteps, read::Read, has_qual::Bool)
    for c in s.cuts
        if c > 0
            read = slice(read, c, length(read), has_qual)
        elseif c < 0
            read = slice(read, 0, max(length(read) + c, 0), has_qual)
        end
    end
    if s.quality !== nothing && has_qual
        start, stop = quality_trim_index(read.qual, s.quality[1], s.quality[2], s.quality_base)
        read = slice(read, start, stop, has_qual)
    end
    return read
end
