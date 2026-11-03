# Semiglobal alignment of an adapter (reference) against a read (query).
#
# This is a line-by-line translation of `Aligner.locate` and of
# `PrefixComparer`/`SuffixComparer` from cutadapt 5.2 (`src/cutadapt/_align.pyx`,
# MIT license, see NOTICE). Costs are edit distances; among the alignments
# whose error rate is acceptable, the one with the highest score (match +1,
# mismatch -1, indel -2) wins. Tie-breaking, Ukkonen's cut-off and the
# acceptance rules are kept exactly, including the use of the `origin` of the
# last computed cell when the final column is scanned.

const MATCH_SCORE = 1
const MISMATCH_SCORE = -1
const INSERTION_SCORE = -2
const DELETION_SCORE = -2

# Where an alignment may start or stop without cost (cutadapt `EndSkip`).
const REFERENCE_START = 1
const QUERY_START = 2
const REFERENCE_END = 4
const QUERY_STOP = 8

# Result of `locate`: 0-based half-open intervals, as in cutadapt.
struct Alignment
    ref_start::Int
    ref_stop::Int
    query_start::Int
    query_stop::Int
    score::Int
    errors::Int
end

abstract type AbstractLocator end

struct Aligner <: AbstractLocator
    ref::Vector{UInt8}          # translated reference
    m::Int
    max_error_rate::Float64
    start_in_reference::Bool
    start_in_query::Bool
    stop_in_reference::Bool
    stop_in_query::Bool
    wildcard_ref::Bool
    query_table::Vector{UInt8}
    compare_ascii::Bool
    indel_cost::Int
    min_overlap::Int
    n_counts::Vector{Int}       # n_counts[i+1] = number of N in reference[1:i]
    effective_length::Int
end

function Aligner(reference::AbstractString, max_error_rate::Real; flags::Integer = 15,
                 wildcard_ref::Bool = false, wildcard_query::Bool = false,
                 indel_cost::Integer = 1, min_overlap::Integer = 1)
    indel_cost >= 1 || throw(ArgumentError("indel_cost must be at least 1"))
    bytes = codeunits(String(reference))
    m = length(bytes)
    n_counts = zeros(Int, m + 1)
    nc = 0
    for i in 1:m
        n_counts[i] = nc
        (bytes[i] == UInt8('N') || bytes[i] == UInt8('n')) && (nc += 1)
    end
    n_counts[m+1] = nc
    effective_length = m
    if wildcard_ref
        effective_length = m - nc
        effective_length == 0 && throw(ArgumentError("Cannot have only N wildcards in the sequence"))
    end
    ref_table, query_table, compare_ascii = match_tables(wildcard_ref, wildcard_query)
    # Without wildcards cutadapt keeps the reference as given (it is already
    # uppercase for adapters) and uppercases only the query.
    ref = compare_ascii ? collect(bytes) : translate(bytes, ref_table)
    return Aligner(ref, m, Float64(max_error_rate),
                   (flags & REFERENCE_START) != 0, (flags & QUERY_START) != 0,
                   (flags & REFERENCE_END) != 0, (flags & QUERY_STOP) != 0,
                   wildcard_ref, query_table, compare_ascii, indel_cost, min_overlap,
                   n_counts, effective_length)
end

# Reusable buffers, one per task.
mutable struct AlignWorkspace
    query::Vector{UInt8}
    cost::Vector{Int}
    score::Vector{Int}
    origin::Vector{Int}
    rc_seq::Vector{UInt8}       # reverse complement of the current read
    rc_qual::Vector{UInt8}
end
AlignWorkspace() = AlignWorkspace(UInt8[], Int[], Int[], Int[], UInt8[], UInt8[])

function _column!(ws::AlignWorkspace, m::Int)
    if length(ws.cost) < m + 1
        resize!(ws.cost, m + 1)
        resize!(ws.score, m + 1)
        resize!(ws.origin, m + 1)
    end
    return ws.cost, ws.score, ws.origin
end

const NO_ALIGNMENT = Alignment(-1, -1, -1, -1, 0, -1)
found(a::Alignment) = a.errors >= 0

"""
    locate(aligner, query::AbstractVector{UInt8}, ws) -> Alignment

Best alignment of the adapter in `query`, or `NO_ALIGNMENT` (check with
`found`). Coordinates are 0-based and half-open, as in cutadapt.
"""
function locate(al::Aligner, query::AbstractVector{UInt8}, ws::AlignWorkspace)
    m = al.m
    n = length(query)
    s1 = al.ref
    s2 = translate!(ws.query, query, al.query_table)
    cost_col, score_col, origin_col = _column!(ws, m)
    max_error_rate = al.max_error_rate
    stop_in_query = al.stop_in_query
    compare_ascii = al.compare_ascii
    deletion_cost = al.indel_cost
    insertion_cost = al.indel_cost

    # maximum no. of errors (C cast of a non-negative double: truncation)
    k = unsafe_trunc(Int, max_error_rate * m)

    max_n = n
    min_n = 0
    if !al.start_in_query
        max_n = min(n, m + k)
    end
    if !stop_in_query
        min_n = max(0, n - m - k)
    end

    # Column min_n. Index i of cutadapt is i + 1 here.
    @inbounds if !al.start_in_reference && !al.start_in_query
        for i in 0:m
            score_col[i+1] = i * DELETION_SCORE
            cost_col[i+1] = max(i, min_n) * deletion_cost
            origin_col[i+1] = 0
        end
    elseif al.start_in_reference && !al.start_in_query
        for i in 0:m
            score_col[i+1] = 0
            cost_col[i+1] = min_n * deletion_cost
            origin_col[i+1] = min(0, min_n - i)
        end
    elseif !al.start_in_reference && al.start_in_query
        for i in 0:m
            score_col[i+1] = i * DELETION_SCORE
            cost_col[i+1] = i * deletion_cost
            origin_col[i+1] = max(0, min_n - i)
        end
    else
        for i in 0:m
            score_col[i+1] = 0
            cost_col[i+1] = min(i, min_n) * deletion_cost
            origin_col[i+1] = min_n - i
        end
    end

    best_ref_stop = m
    best_query_stop = n
    best_cost = m + n + 1
    best_origin = 0
    best_score = 0
    no_match = m + n + 1

    # Ukkonen's trick: index of the last cell that is at most k
    last = min(m, k + 1)
    al.start_in_reference && (last = m)

    origin_increment = al.start_in_query ? 1 : 0
    insertion_cost_increment = al.start_in_query ? 0 : insertion_cost
    insertion_score_increment = al.start_in_query ? 0 : INSERTION_SCORE
    last_filled_i = 0
    # `origin` survives the loops: cutadapt reuses it in the final column scan.
    origin = 0
    cost = 0
    score = 0

    @inbounds for j in (min_n+1):max_n
        diag_cost = cost_col[1]
        diag_score = score_col[1]
        diag_origin = origin_col[1]
        origin_col[1] += origin_increment
        cost_col[1] += insertion_cost_increment
        score_col[1] += insertion_score_increment
        qj = s2[j]
        for i in 1:last
            characters_equal = compare_ascii ? (s1[i] == qj) : ((s1[i] & qj) != 0x00)
            if characters_equal
                cost = diag_cost
                origin = diag_origin
                score = diag_score + MATCH_SCORE
            else
                cur_cost = cost_col[i+1]
                prev_cost = cost_col[i]
                cost_diag = diag_cost + 1
                cost_insertion = cur_cost + insertion_cost
                cost_deletion = prev_cost + deletion_cost
                if cost_diag <= cost_deletion && cost_diag <= cost_insertion
                    cost = cost_diag
                    origin = diag_origin
                    score = diag_score + MISMATCH_SCORE
                elseif cost_deletion <= cost_insertion
                    cost = cost_deletion
                    origin = origin_col[i]
                    score = score_col[i] + DELETION_SCORE
                else
                    cost = cost_insertion
                    origin = origin_col[i+1]
                    score = score_col[i+1] + INSERTION_SCORE
                end
            end
            diag_cost = cost_col[i+1]
            diag_score = score_col[i+1]
            diag_origin = origin_col[i+1]
            cost_col[i+1] = cost
            origin_col[i+1] = origin
            score_col[i+1] = score
        end
        last_filled_i = last
        while last >= 0 && cost_col[last+1] > k
            last -= 1
        end
        if last < m
            last += 1
        elseif stop_in_query
            cost = cost_col[m+1]
            score = score_col[m+1]
            origin = origin_col[m+1]
            len = m + min(origin, 0)
            cur_effective_length = len
            if al.wildcard_ref
                if len < m
                    cur_effective_length = len - (al.n_counts[m+1] - al.n_counts[m-len+1])
                else
                    cur_effective_length = al.effective_length
                end
            end
            is_acceptable = len >= al.min_overlap && cost <= cur_effective_length * max_error_rate
            best_length = m + min(best_origin, 0)
            if is_acceptable && (best_cost == no_match ||
                                 (origin <= best_origin + m ÷ 2 && score > best_score) ||
                                 (len > best_length && score > best_score))
                best_score = score
                best_cost = cost
                best_origin = origin
                best_ref_stop = m
                best_query_stop = j
                if cost == 0 && origin >= 0
                    break  # exact match, stop early
                end
            end
        end
    end

    @inbounds if max_n == n
        first_i = al.stop_in_reference ? 0 : m
        for i in last_filled_i:-1:first_i
            len = i + min(origin_col[i+1], 0)
            c = cost_col[i+1]
            s = score_col[i+1]
            if al.wildcard_ref
                if len < m
                    ref_start = -min(origin_col[i+1], 0)
                    cur_effective_length = len - (al.n_counts[i+1] - al.n_counts[ref_start+1])
                else
                    cur_effective_length = al.effective_length
                end
            else
                cur_effective_length = len
            end
            is_acceptable = len >= al.min_overlap && c <= cur_effective_length * max_error_rate
            best_length = best_ref_stop + min(best_origin, 0)
            # `origin` is deliberately the stale value from the loops above.
            if is_acceptable && (best_cost == no_match ||
                                 (origin <= best_origin + m ÷ 2 && s > best_score) ||
                                 (len > best_length && s > best_score))
                best_score = s
                best_cost = c
                best_origin = origin_col[i+1]
                best_ref_stop = i
                best_query_stop = n
            end
        end
    end
    best_cost == no_match && return NO_ALIGNMENT

    if best_origin >= 0
        ref_start, query_start = 0, best_origin
    else
        ref_start, query_start = -best_origin, 0
    end
    return Alignment(ref_start, best_ref_stop, query_start, best_query_stop, best_score, best_cost)
end

# Anchored adapters without indels: count mismatches of the aligned prefix
# (or suffix) only.
struct AffixComparer <: AbstractLocator
    ref::Vector{UInt8}          # translated; reversed for suffixes
    m::Int
    max_k::Int
    min_overlap::Int
    query_table::Vector{UInt8}
    compare_ascii::Bool
    suffix::Bool
end

function AffixComparer(reference::AbstractString, max_error_rate::Real; suffix::Bool,
                       wildcard_ref::Bool = false, wildcard_query::Bool = false,
                       min_overlap::Integer = 1)
    ref = String(reference)
    suffix && (ref = reverse(ref))
    m = ncodeunits(ref)
    effective_length = m
    if wildcard_ref
        effective_length -= count(==('N'), ref) - count(==('n'), ref)
        effective_length == 0 && throw(ArgumentError("Cannot have only N wildcards in the sequence"))
    end
    0 <= max_error_rate <= 1 || throw(ArgumentError("max_error_rate must be between 0 and 1"))
    min_overlap >= 1 || throw(ArgumentError("min_overlap must be at least 1"))
    ref_table, query_table, compare_ascii = match_tables(wildcard_ref, wildcard_query)
    return AffixComparer(translate(codeunits(ref), ref_table), m,
                         unsafe_trunc(Int, max_error_rate * effective_length),
                         min_overlap, query_table, compare_ascii, suffix)
end

function locate(c::AffixComparer, query::AbstractVector{UInt8}, ws::AlignWorkspace)
    n = length(query)
    len = min(c.m, n)
    q = translate!(ws.query, query, c.query_table)
    errors = 0
    @inbounds for i in 1:len
        qi = c.suffix ? q[n-i+1] : q[i]
        if c.compare_ascii
            c.ref[i] != qi && (errors += 1)
        else
            (c.ref[i] & qi) == 0x00 && (errors += 1)
        end
    end
    (errors > c.max_k || len < c.min_overlap) && return NO_ALIGNMENT
    score = (len - errors) * MATCH_SCORE + errors * MISMATCH_SCORE
    if c.suffix
        return Alignment(c.m - len, c.m, n - len, n, score, errors)
    end
    return Alignment(0, len, 0, len, score, errors)
end
