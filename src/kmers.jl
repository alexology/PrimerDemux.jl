# K-mer pre-filter that decides whether a read can contain an acceptable
# alignment of an adapter, translated from cutadapt 5.2 (`kmer_heuristic.py`
# and `_kmer_finder.pyx`, MIT license, see NOTICE). The aligner only runs when
# at least one k-mer occurs in its search window, exactly as in cutadapt, so
# results also follow cutadapt where the heuristic is not lossless.

# (start, stop, kmers): stop == nothing means "to the end of the read";
# negative values count from the end of the read.
const SearchSet = Tuple{Int,Union{Int,Nothing},Vector{String}}

function kmer_chunks(sequence::String, chunks::Int)
    len = ncodeunits(sequence)
    chunk_size, remainder = divrem(len, chunks)
    sizes = vcat(fill(chunk_size + 1, remainder), fill(chunk_size, chunks - remainder))
    result = String[]
    offset = 0
    for size in sizes
        kmer = sequence[(offset+1):(offset+size)]
        kmer in result || push!(result, kmer)
        offset += size
    end
    return result
end

function create_back_overlap_searchsets(adapter::String, min_overlap::Int, error_rate::Float64)
    adapter_length = ncodeunits(adapter)
    error_lengths = Tuple{Int,Int}[]
    max_error = 0
    for i in 0:adapter_length
        if unsafe_trunc(Int, i * error_rate) > max_error
            push!(error_lengths, (max_error, i - 1))
            max_error += 1
        end
    end
    push!(error_lengths, (max_error, adapter_length))

    search_sets = SearchSet[]
    minimum_length = min_overlap
    for (max_errors, len) in error_lengths
        minimum_length > len && continue
        if max_errors == 0
            # Directly matching 1- to 4-mers at the very end of the read.
            min_overlap_kmer_length = 5
            if minimum_length < min_overlap_kmer_length
                for i in minimum_length:(min_overlap_kmer_length-1)
                    push!(search_sets, (-i, nothing, [adapter[1:min(i, adapter_length)]]))
                end
                minimum_length = min_overlap_kmer_length
            end
        end
        kmer_sets = kmer_chunks(adapter[1:min(minimum_length, adapter_length)], max_errors + 1)
        push!(search_sets, (-len, nothing, kmer_sets))
        minimum_length = len + 1
    end
    return search_sets
end

function create_positions_and_kmers(adapter::String, min_overlap::Int, error_rate::Float64,
                                    back_adapter::Bool, front_adapter::Bool, internal::Bool)
    max_errors = unsafe_trunc(Int, ncodeunits(adapter) * error_rate)
    search_sets = SearchSet[]
    back_adapter && append!(search_sets, create_back_overlap_searchsets(adapter, min_overlap, error_rate))
    if front_adapter
        # Reverse the adapter, build the back searches, reverse them again.
        for (start, _, kmer_set) in create_back_overlap_searchsets(reverse(adapter), min_overlap, error_rate)
            push!(search_sets, (0, -start, unique!([reverse(k) for k in kmer_set])))
        end
    end
    internal && push!(search_sets, (0, nothing, kmer_chunks(adapter, max_errors + 1)))
    return remove_redundant_kmers(search_sets)
end

function minimize_kmer_search_list(kmer_search_list::Vector{Tuple{String,Int,Union{Int,Nothing}}})
    kmers = String[]
    positions = Dict{String,Vector{Tuple{Int,Union{Int,Nothing}}}}()
    for (kmer, start, stop) in kmer_search_list
        haskey(positions, kmer) || push!(kmers, kmer)
        push!(get!(positions, kmer, Tuple{Int,Union{Int,Nothing}}[]), (start, stop))
    end
    result = Tuple{String,Int,Union{Int,Nothing}}[]
    for kmer in kmers
        pos = positions[kmer]
        if length(pos) == 1
            push!(result, (kmer, pos[1][1], pos[1][2]))
            continue
        end
        if any(p -> p[1] == 0 && p[2] === nothing, pos)
            push!(result, (kmer, 0, nothing))
            continue
        end
        front = [p for p in pos if p[1] == 0]
        back = [p for p in pos if p[2] === nothing]
        any(p -> p[1] != 0 && p[2] !== nothing, pos) &&
            error("k-mer searches starting in the middle are not supported")
        isempty(front) || push!(result, (kmer, 0, maximum(p[2]::Int for p in front)))
        isempty(back) || push!(result, (kmer, minimum(p[1] for p in back), nothing))
    end
    return result
end

function remove_redundant_kmers(search_sets::Vector{SearchSet})
    flat = Tuple{String,Int,Union{Int,Nothing}}[]
    for (start, stop, kmer_set) in search_sets, kmer in kmer_set
        push!(flat, (kmer, start, stop))
    end
    keys_ = Tuple{Int,Union{Int,Nothing}}[]
    groups = Dict{Tuple{Int,Union{Int,Nothing}},Vector{String}}()
    for (kmer, start, stop) in minimize_kmer_search_list(flat)
        key = (start, stop)
        haskey(groups, key) || push!(keys_, key)
        push!(get!(groups, key, String[]), kmer)
    end
    return SearchSet[(k[1], k[2], groups[k]) for k in keys_]
end

# One shift-and search: up to 64 bits of concatenated words.
struct KmerSearch
    start::Int
    stop::Int                   # 0 = to the end of the read
    init_mask::UInt64
    found_mask::UInt64
    masks::Vector{UInt64}       # 256 entries, indexed by read byte + 1
end

struct KmerFinder
    searches::Vector{KmerSearch}
    always::Bool                # true when a search can never rule out a match
end

const MAX_WORD_SIZE = 64

function KmerFinder(positions_and_kmers::Vector{SearchSet}, ref_wildcards::Bool, query_wildcards::Bool)
    ref_table, query_table, compare_ascii = match_tables(ref_wildcards, query_wildcards)
    matches(r::UInt8, q::UInt8) = q != 0x00 && q < 0x80 &&
        (compare_ascii ? ref_table[r+1] == query_table[q+1] : (ref_table[r+1] & query_table[q+1]) != 0x00)
    searches = KmerSearch[]
    for (start, stop, kmers) in positions_and_kmers
        # An empty word is trivially present (cutadapt would shift by -1 bits).
        any(isempty, kmers) && return KmerFinder(KmerSearch[], true)
        index = 1
        while index <= length(kmers)
            offset = 0
            init_mask = UInt64(0)
            found_mask = UInt64(0)
            masks = zeros(UInt64, 256)
            while index <= length(kmers)
                kmer = codeunits(kmers[index])
                len = length(kmer)
                # Words longer than a machine word make cutadapt skip the filter.
                len > MAX_WORD_SIZE && return KmerFinder(KmerSearch[], true)
                offset + len > MAX_WORD_SIZE && break
                init_mask |= UInt64(1) << offset
                found_mask |= UInt64(1) << (offset + len - 1)
                for (p, r) in enumerate(kmer), q in 0x00:0x7f
                    matches(r, q) && (masks[q+1] |= UInt64(1) << (offset + p - 1))
                end
                offset += len
                index += 1
            end
            push!(searches, KmerSearch(start, stop === nothing ? 0 : stop, init_mask, found_mask, masks))
        end
    end
    return KmerFinder(searches, false)
end

function kmers_present(kf::KmerFinder, seq::AbstractVector{UInt8})
    kf.always && return true
    seq_length = length(seq)
    @inbounds for e in kf.searches
        start = e.start
        stop = e.stop
        if start < 0
            start = max(seq_length + start, 0)
        elseif start > seq_length
            continue
        end
        if stop < 0
            stop = seq_length + stop
            stop <= 0 && continue
        elseif stop == 0
            stop = seq_length
        elseif stop > seq_length
            # cutadapt reads past the end of the string here, with an outcome
            # that depends on memory contents; let the aligner decide instead.
            return true
        end
        stop - start <= 0 && continue
        R = UInt64(0)
        masks = e.masks
        for i in (start+1):stop
            R = ((R << 1) | e.init_mask) & masks[seq[i]+1]
            (R & e.found_mask) != 0 && return true
        end
    end
    return false
end
