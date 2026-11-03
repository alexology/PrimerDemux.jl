# Whole-file processing: modifications, filters, routing to outputs,
# statistics and multithreading.

"""
    TrimStats

Summary returned by [`remove_primers`](@ref) and [`demultiplex`](@ref). Counts
refer to reads (single-end) or pairs (paired-end); `bp_in`/`bp_out` sum both
mates. `with_adapters[i]` and `adapter_matches[i]` refer to mate `i`.
`written` gives the reads or pairs written to each output (by name for
demultiplexing, `"untrimmed"` for `untrimmed_output`).
"""
mutable struct TrimStats
    paired::Bool
    reads_in::Int
    bp_in::Int
    with_adapters::Vector{Int}
    reverse_complemented::Int
    adapter_names::Vector{Vector{String}}
    adapter_matches::Vector{Vector{Int}}
    too_short::Int
    too_long::Int
    too_many_n::Int
    too_many_expected_errors::Int
    discard_trimmed::Int
    discard_untrimmed::Int
    reads_out::Int
    bp_out::Int
    written::Vector{Pair{String,Int}}
end

function TrimStats(paired::Bool, adapters1, adapters2)
    names = [String[adapter_name(a) for a in adapters1], String[adapter_name(a) for a in adapters2]]
    paired || (names = names[1:1])
    return TrimStats(paired, 0, 0, zeros(Int, length(names)), 0, names,
                     [zeros(Int, length(n)) for n in names], 0, 0, 0, 0, 0, 0, 0, 0,
                     Pair{String,Int}[])
end

function merge!(a::TrimStats, b::TrimStats)
    a.reads_in += b.reads_in
    a.bp_in += b.bp_in
    a.with_adapters .+= b.with_adapters
    a.reverse_complemented += b.reverse_complemented
    for i in eachindex(a.adapter_matches)
        a.adapter_matches[i] .+= b.adapter_matches[i]
    end
    for f in (:too_short, :too_long, :too_many_n, :too_many_expected_errors,
              :discard_trimmed, :discard_untrimmed)
        setfield!(a, f, getfield(a, f) + getfield(b, f))
    end
    return a
end

_pct(x, total) = total == 0 ? "0.0" : string(round(100 * x / total; digits = 1))

function Base.show(io::IO, ::MIME"text/plain", s::TrimStats)
    unit = s.paired ? "pairs" : "reads"
    println(io, "PrimerDemux summary")
    println(io, "  Total $unit processed: ", s.reads_in)
    for (i, n) in enumerate(s.with_adapters)
        label = s.paired ? "  Read $i with adapter: " : "  Reads with adapters: "
        println(io, label, n, " (", _pct(n, s.reads_in), "%)")
    end
    s.reverse_complemented > 0 && println(io, "  Reverse-complemented: ", s.reverse_complemented)
    for (label, v) in (("too short", s.too_short), ("too long", s.too_long),
                       ("with too many N", s.too_many_n),
                       ("with too many expected errors", s.too_many_expected_errors),
                       ("discarded as trimmed", s.discard_trimmed),
                       ("discarded as untrimmed", s.discard_untrimmed))
        v > 0 && println(io, "  ", uppercasefirst(unit), " ", label, ": ", v, " (", _pct(v, s.reads_in), "%)")
    end
    println(io, "  ", uppercasefirst(unit), " written: ", s.reads_out, " (", _pct(s.reads_out, s.reads_in), "%)")
    println(io, "  Basepairs: ", s.bp_in, " processed, ", s.bp_out, " written (", _pct(s.bp_out, s.bp_in), "%)")
    for (mate, names) in enumerate(s.adapter_names), (j, name) in enumerate(names)
        prefix = s.paired ? "R$mate " : ""
        println(io, "  ", prefix, "adapter ", repr(name), ": ", s.adapter_matches[mate][j], " matches")
    end
    if !isempty(s.written)
        println(io, "  Written per output:")
        for (name, n) in s.written
            println(io, "    ", name, ": ", n)
        end
    end
end
Base.show(io::IO, s::TrimStats) = print(io, "TrimStats(", s.reads_in, " in, ", s.reads_out, " out)")

# --- configuration -----------------------------------------------------------

const PAIR_FILTERS = (:any, :both, :first)

struct Config
    paired::Bool
    has_qual::Bool
    steps1::MateSteps
    steps2::MateSteps
    cutter1::Union{Nothing,AdapterCutter}
    cutter2::Union{Nothing,AdapterCutter}
    pair_adapters::Bool
    revcomp::Bool
    rc_suffix::Vector{UInt8}
    min_length::Tuple{Union{Nothing,Int},Union{Nothing,Int}}
    max_length::Tuple{Union{Nothing,Int},Union{Nothing,Int}}
    max_n::Union{Nothing,Float64}
    max_ee::Union{Nothing,Float64}
    pair_filter::Symbol
    untrimmed_pair_filter::Symbol
    discard_trimmed::Bool
    discard_untrimmed::Bool
    mode::Symbol                # :single, :demux, :combinatorial
    route::Dict{Any,Int}        # demultiplexing: name(s) => output index
    untrimmed_dest::Int         # 0 = discarded
end

# Result of processing one read or pair.
struct Outcome
    dest::Int                   # 0 = filtered out
    r1::Read
    r2::Read
end

const EMPTY_READ = Read()

function _cutter_matches!(stats::TrimStats, mate::Int, matches::Vector{AdapterMatch})
    isempty(matches) && return
    stats.with_adapters[mate] += 1
    for m in matches
        stats.adapter_matches[mate][m.adapter] += 1
    end
end

function trim_single(cfg::Config, read::Read, ws::AlignWorkspace, stats::TrimStats)
    read = apply_steps(cfg.steps1, read, cfg.has_qual)
    cfg.cutter1 === nothing && return read, AdapterMatch[]
    if !cfg.revcomp
        read, matches = match_and_trim(cfg.cutter1, read, cfg.has_qual, ws)
        _cutter_matches!(stats, 1, matches)
        return read, matches
    end
    # The reverse complement lives in the workspace until the next read.
    rc = reverse_complement!(Read(UInt8[], ws.rc_seq, ws.rc_qual), read, cfg.has_qual)
    fwd, fwd_matches = match_and_trim(cfg.cutter1, read, cfg.has_qual, ws)
    rev, rev_matches = match_and_trim(cfg.cutter1, rc, cfg.has_qual, ws)
    if sum(score, rev_matches; init = 0) > sum(score, fwd_matches; init = 0)
        stats.reverse_complemented += 1
        rev = with_suffix(rev, cfg.rc_suffix)
        _cutter_matches!(stats, 1, rev_matches)
        return rev, rev_matches
    end
    _cutter_matches!(stats, 1, fwd_matches)
    return fwd, fwd_matches
end

function _trim_mate(c::Union{Nothing,AdapterCutter}, read::Read, has_qual::Bool, ws)
    c === nothing && return read, AdapterMatch[]
    return match_and_trim(c, read, has_qual, ws)
end

function trim_pair(cfg::Config, r1::Read, r2::Read, ws::AlignWorkspace, stats::TrimStats)
    r1 = apply_steps(cfg.steps1, r1, cfg.has_qual)
    r2 = apply_steps(cfg.steps2, r2, cfg.has_qual)
    if cfg.pair_adapters
        return trim_adapter_pairs(cfg, r1, r2, ws, stats)
    end
    t1, m1 = _trim_mate(cfg.cutter1, r1, cfg.has_qual, ws)
    t2, m2 = _trim_mate(cfg.cutter2, r2, cfg.has_qual, ws)
    if cfg.revcomp
        unswapped = sum(score, m1; init = 0) + sum(score, m2; init = 0)
        s1, n1 = _trim_mate(cfg.cutter1, r2, cfg.has_qual, ws)
        s2, n2 = _trim_mate(cfg.cutter2, r1, cfg.has_qual, ws)
        if sum(score, n1; init = 0) + sum(score, n2; init = 0) > unswapped
            stats.reverse_complemented += 1
            t1 = with_suffix(s1, cfg.rc_suffix)
            t2 = with_suffix(s2, cfg.rc_suffix)
            m1, m2 = n1, n2
        end
    end
    _cutter_matches!(stats, 1, m1)
    _cutter_matches!(stats, 2, m2)
    return t1, t2, m1, m2
end

# cutadapt --pair-adapters: adapters1[i] must match R1 and adapters2[i] R2.
function trim_adapter_pairs(cfg::Config, r1::Read, r2::Read, ws::AlignWorkspace, stats::TrimStats)
    a1 = cfg.cutter1.adapters
    a2 = cfg.cutter2.adapters
    best = nothing
    best_score = best_errors = 0
    for i in eachindex(a1)
        m1 = match_to(a1[i], i, r1.seq, ws)
        m1 === nothing && continue
        m2 = match_to(a2[i], i, r2.seq, ws)
        m2 === nothing && continue
        total_score = score(m1) + score(m2)
        total_errors = errors(m1) + errors(m2)
        if best === nothing || total_score > best_score ||
           (total_score == best_score && total_errors < best_errors)
            best = (m1, m2)
            best_score, best_errors = total_score, total_errors
        end
    end
    best === nothing && return r1, r2, AdapterMatch[], AdapterMatch[]
    m1, m2 = [best[1]], [best[2]]
    t1 = apply_action(r1, m1, cfg.cutter1.action, cfg.has_qual)
    t2 = apply_action(r2, m2, cfg.cutter1.action, cfg.has_qual)
    _cutter_matches!(stats, 1, m1)
    _cutter_matches!(stats, 2, m2)
    return t1, t2, m1, m2
end

# --- filters -----------------------------------------------------------------

function _pair_filtered(p1::Bool, p2::Bool, has1::Bool, has2::Bool, mode::Symbol)
    !has2 && return p1
    !has1 && return p2
    mode == :any && return p1 || p2
    mode == :both && return p1 && p2
    return p1
end

too_many_n(r::Read, cutoff::Float64) = begin
    n = count(b -> b == UInt8('N') || b == UInt8('n'), r.seq)
    cutoff < 1 ? (length(r) == 0 ? false : n / length(r) > cutoff) : n > cutoff
end

function route(cfg::Config, r1::Read, r2::Read, m1, m2, stats::TrimStats)
    paired = cfg.paired
    lmin1, lmin2 = cfg.min_length
    if lmin1 !== nothing || (paired && lmin2 !== nothing)
        f1 = lmin1 !== nothing && length(r1) < lmin1
        f2 = paired && lmin2 !== nothing && length(r2) < lmin2
        if _pair_filtered(f1, f2, lmin1 !== nothing, paired && lmin2 !== nothing, cfg.pair_filter)
            stats.too_short += 1
            return 0
        end
    end
    lmax1, lmax2 = cfg.max_length
    if lmax1 !== nothing || (paired && lmax2 !== nothing)
        f1 = lmax1 !== nothing && length(r1) > lmax1
        f2 = paired && lmax2 !== nothing && length(r2) > lmax2
        if _pair_filtered(f1, f2, lmax1 !== nothing, paired && lmax2 !== nothing, cfg.pair_filter)
            stats.too_long += 1
            return 0
        end
    end
    if cfg.max_n !== nothing
        if _pair_filtered(too_many_n(r1, cfg.max_n), paired && too_many_n(r2, cfg.max_n), true, paired, cfg.pair_filter)
            stats.too_many_n += 1
            return 0
        end
    end
    if cfg.max_ee !== nothing && cfg.has_qual
        e1 = expected_errors(r1.qual, cfg.steps1.quality_base) > cfg.max_ee
        e2 = paired && expected_errors(r2.qual, cfg.steps2.quality_base) > cfg.max_ee
        if _pair_filtered(e1, e2, true, paired, cfg.pair_filter)
            stats.too_many_expected_errors += 1
            return 0
        end
    end
    untrimmed1, untrimmed2 = isempty(m1), isempty(m2)
    if cfg.mode == :demux
        untrimmed1 || return cfg.route[adapter_name_of(cfg, 1, m1[end])]
        cfg.untrimmed_dest == 0 && (stats.discard_untrimmed += 1)
        return cfg.untrimmed_dest
    elseif cfg.mode == :combinatorial
        key = (untrimmed1 ? nothing : adapter_name_of(cfg, 1, m1[end]),
               untrimmed2 ? nothing : adapter_name_of(cfg, 2, m2[end]))
        dest = get(cfg.route, key, 0)
        dest == 0 && (stats.discard_untrimmed += 1)
        return dest
    end
    if cfg.discard_trimmed
        if _pair_filtered(!untrimmed1, paired && !untrimmed2, true, paired, cfg.pair_filter)
            stats.discard_trimmed += 1
            return 0
        end
    elseif cfg.discard_untrimmed || cfg.untrimmed_dest != 0
        if _pair_filtered(untrimmed1, untrimmed2, true, paired, cfg.untrimmed_pair_filter)
            cfg.untrimmed_dest == 0 && (stats.discard_untrimmed += 1)
            return cfg.untrimmed_dest
        end
    end
    return 1
end

adapter_name_of(cfg::Config, mate::Int, m::AdapterMatch) =
    adapter_name((mate == 1 ? cfg.cutter1 : cfg.cutter2).adapters[m.adapter])

function process(cfg::Config, r1::Read, r2::Read, ws::AlignWorkspace, stats::TrimStats)
    stats.reads_in += 1
    stats.bp_in += length(r1) + length(r2)
    if cfg.paired
        t1, t2, m1, m2 = trim_pair(cfg, r1, r2, ws, stats)
    else
        t1, m1 = trim_single(cfg, r1, ws, stats)
        t2, m2 = EMPTY_READ, AdapterMatch[]
    end
    dest = route(cfg, t1, t2, m1, m2, stats)
    return Outcome(dest, t1, t2)
end

# --- driver ------------------------------------------------------------------

_mates(x::Tuple) = (length(x) == 2 || throw(ArgumentError("expected two paths for paired-end data, got $(length(x))")); (String(x[1]), String(x[2])))
_mates(x::AbstractString) = (String(x),)

_lengths(x::Nothing, paired) = (nothing, nothing)
_lengths(x::Integer, paired) = (Int(x), paired ? Int(x) : nothing)
function _lengths(x::Tuple, paired)
    paired || throw(ArgumentError("two lengths are only allowed for paired-end data"))
    length(x) == 2 || throw(ArgumentError("expected two lengths"))
    return (x[1] === nothing ? nothing : Int(x[1]), x[2] === nothing ? nothing : Int(x[2]))
end

_cuts(x::Integer) = x == 0 ? Int[] : [Int(x)]
_cuts(x::AbstractVector) = Int[c for c in x if c != 0]
_cuts(x::Tuple) = _cuts(collect(x))

_quality(x::Nothing) = nothing
_quality(x::Integer) = (0, Int(x))
_quality(x::Tuple) = (Int(x[1]), Int(x[2]))

_adapters(x::Adapter) = Adapter[x]
_adapters(x::AbstractVector) = Adapter[a for a in x]

function run_pipeline(input, output; demux::Bool, adapters = Adapter[], adapters2 = Adapter[],
                      times::Integer = 1, action::Symbol = :trim, revcomp::Bool = false,
                      rc_suffix::AbstractString = " rc", pair_adapters::Bool = false,
                      cut = 0, cut2 = 0, quality_cutoff = nothing, quality_cutoff2 = nothing,
                      quality_base::Integer = 33, min_length = nothing, max_length = nothing,
                      max_n::Union{Nothing,Real} = nothing,
                      max_expected_errors::Union{Nothing,Real} = nothing,
                      discard_trimmed::Bool = false, discard_untrimmed::Bool = false,
                      untrimmed_output = nothing, pair_filter::Union{Nothing,Symbol} = nothing,
                      threads::Integer = Threads.nthreads(), batch_size::Integer = 4096,
                      compression_level::Integer = 1, overwrite::Bool = false)
    inputs = _mates(input)
    outputs = _mates(output)
    paired = length(inputs) == 2
    length(outputs) == length(inputs) ||
        throw(ArgumentError("give as many outputs as inputs (one for single-end, two for paired-end)"))
    adapters1 = _adapters(adapters)
    adapters2 = _adapters(adapters2)
    paired || isempty(adapters2) || throw(ArgumentError("adapters2 requires paired-end input"))
    check_action(action, times)
    threads >= 1 || throw(ArgumentError("threads must be at least 1"))
    batch_size >= 1 || throw(ArgumentError("batch_size must be at least 1"))
    discard_trimmed + discard_untrimmed + (untrimmed_output !== nothing) > 1 && throw(ArgumentError(
        "Only one of discard_trimmed, discard_untrimmed and untrimmed_output can be used at the same time."))
    pair_filter === nothing || pair_filter in PAIR_FILTERS ||
        throw(ArgumentError("pair_filter must be one of $(PAIR_FILTERS)"))
    if pair_adapters
        paired || throw(ArgumentError("pair_adapters requires paired-end input"))
        length(adapters1) == length(adapters2) && !isempty(adapters1) || throw(ArgumentError(
            "pair_adapters needs the same, non-zero number of adapters for R1 and R2"))
        times == 1 || throw(ArgumentError("times cannot be used with pair_adapters"))
        revcomp && throw(ArgumentError("revcomp cannot be used with pair_adapters"))
    end
    (max_n === nothing || max_n >= 0) || throw(ArgumentError("max_n must not be negative"))

    # Output layout
    mode = :single
    if demux
        if paired && all(o -> occursin("{name1}", o) && occursin("{name2}", o), outputs)
            mode = :combinatorial
        elseif all(o -> occursin("{name}", o), outputs)
            mode = :demux
        else
            throw(ArgumentError("demultiplex output must contain {name} (or {name1} and {name2} for combinatorial paired-end demultiplexing)"))
        end
        isempty(adapters1) && throw(ArgumentError("demultiplexing needs adapters"))
        discard_trimmed && throw(ArgumentError("discard_trimmed cannot be used for demultiplexing"))
        if mode == :combinatorial
            isempty(adapters2) && throw(ArgumentError("combinatorial demultiplexing needs adapters2"))
            pair_adapters && throw(ArgumentError("use {name} (not {name1} and {name2}) with pair_adapters"))
            untrimmed_output === nothing || throw(ArgumentError("untrimmed_output cannot be used with combinatorial demultiplexing"))
        end
    else
        any(o -> occursin(r"\{name[12]?\}", o), outputs) &&
            throw(ArgumentError("use demultiplex for outputs with {name}"))
    end
    untrimmed_paths = untrimmed_output === nothing ? nothing : _mates(untrimmed_output)
    untrimmed_paths === nothing || length(untrimmed_paths) == length(inputs) ||
        throw(ArgumentError("untrimmed_output needs one path per input"))

    # Destinations: names and paths of each output, in order.
    dests = Tuple{String,Vector{String}}[]
    routes = Dict{Any,Int}()
    untrimmed_dest = 0
    if mode == :single
        push!(dests, ("output", collect(outputs)))
        if untrimmed_paths !== nothing
            push!(dests, ("untrimmed", collect(untrimmed_paths)))
            untrimmed_dest = 2
        end
    elseif mode == :demux
        for a in adapters1
            name = adapter_name(a)
            haskey(routes, name) && continue
            push!(dests, (name, [replace(o, "{name}" => name) for o in outputs]))
            routes[name] = length(dests)
        end
        if !discard_untrimmed
            paths = untrimmed_paths === nothing ? [replace(o, "{name}" => "unknown") for o in outputs] :
                    collect(untrimmed_paths)
            push!(dests, ("unknown", paths))
            untrimmed_dest = length(dests)
        end
    else
        names1 = unique(adapter_name.(adapters1))
        names2 = unique(adapter_name.(adapters2))
        keys_ = Any[(a, b) for a in names1 for b in names2]
        if !discard_untrimmed
            push!(keys_, (nothing, nothing))
            append!(keys_, [(nothing, b) for b in names2])
            append!(keys_, [(a, nothing) for a in names1])
        end
        for key in keys_
            f1 = something(key[1], "unknown")
            f2 = something(key[2], "unknown")
            push!(dests, ("$f1-$f2", [replace(o, "{name1}" => f1, "{name2}" => f2) for o in outputs]))
            routes[key] = length(dests)
        end
    end
    all_paths = reduce(vcat, (d[2] for d in dests); init = String[])
    for p in all_paths
        isdir(p) && throw(ArgumentError("output $(repr(p)) is a directory"))
        !overwrite && ispath(p) && throw(ArgumentError("output $(repr(p)) exists; pass overwrite=true to replace it"))
    end
    seen = Dict{Any,String}()
    for p in all_paths, k in path_keys(p)
        haskey(seen, k) && throw(ArgumentError(
            "outputs $(repr(seen[k])) and $(repr(p)) would be the same file; use distinct names"))
        seen[k] = p
    end
    for p in inputs, k in path_keys(p)
        haskey(seen, k) && throw(ArgumentError("output $(repr(seen[k])) is also the input $(repr(p))"))
    end

    sources = InputSource[]
    try
        for p in inputs
            push!(sources, open_input(p; threaded = threads > 1))
        end
        is_fastq = [s.fastq for s in sources]
        if paired && is_fastq[1] != is_fastq[2]
            throw(ArgumentError("paired inputs must both be FASTQ or both FASTA"))
        end
        has_qual = is_fastq[1]
        untrimmed_pair_filter = something(pair_filter, :any)
        if paired && (isempty(adapters1) || isempty(adapters2)) &&
           (discard_untrimmed || untrimmed_output !== nothing)
            untrimmed_pair_filter = :both
        end
        cfg = Config(paired, has_qual,
                     MateSteps(_cuts(cut), _quality(quality_cutoff), quality_base),
                     MateSteps(_cuts(cut2), _quality(quality_cutoff2), quality_base),
                     isempty(adapters1) ? nothing : AdapterCutter(adapters1, times, action),
                     isempty(adapters2) ? nothing : AdapterCutter(adapters2, times, action),
                     pair_adapters, revcomp, collect(codeunits(String(rc_suffix))),
                     _lengths(min_length, paired), _lengths(max_length, paired),
                     max_n === nothing ? nothing : Float64(max_n),
                     max_expected_errors === nothing ? nothing : Float64(max_expected_errors),
                     something(pair_filter, :any), untrimmed_pair_filter,
                     discard_trimmed, discard_untrimmed, mode, routes, untrimmed_dest)
        files = Vector{Vector{OutputFile}}()
        try
            for (_, paths) in dests
                fs = OutputFile[]
                push!(files, fs)
                for p in paths
                    # registered at once, so that a later failure discards it
                    push!(fs, OutputFile(p, has_qual, Int(compression_level)))
                end
            end
            stats = TrimStats(paired, adapters1, adapters2)
            written = zeros(Int, length(dests))
            run_batches!(cfg, sources, files, stats, written, Int(threads), Int(batch_size))
            stats.written = [dests[i][1] => written[i] for i in eachindex(dests)]
            publish_all(reduce(vcat, files; init = OutputFile[]))
            return stats
        catch
            foreach(fs -> foreach(discard, fs), files)
            rethrow()
        end
    finally
        foreach(close_input, sources)
    end
end

# Stats with the adapter lists of `cfg`, for one batch.
batch_stats(cfg::Config) =
    TrimStats(cfg.paired, cfg.cutter1 === nothing ? Adapter[] : cfg.cutter1.adapters,
              cfg.cutter2 === nothing ? Adapter[] : cfg.cutter2.adapters)

# Reads of one batch; their vectors are reused by later batches.
mutable struct Batch
    r1::Vector{Read}
    r2::Vector{Read}
    n::Int
end
Batch() = Batch(Read[], Read[], 0)

# Per-task buffers: alignment scratch and formatted output per destination/mate.
struct Workspace
    align::AlignWorkspace
    bufs::Vector{Vector{Vector{UInt8}}}
end
Workspace(ndest::Int, nmates::Int) =
    Workspace(AlignWorkspace(), [[UInt8[] for _ in 1:nmates] for _ in 1:ndest])

struct BatchResult
    blocks::Vector{Vector{Vector{UInt8}}}   # [destination][mate], ready to append
    counts::Vector{Int}                     # reads or pairs per destination
    stats::TrimStats
end

# Processes one batch: trimming, filtering, formatting and compression. Each
# read is formatted before the next one is processed, so trimmed reads may
# share vectors with the batch and the workspace.
function process_batch(cfg::Config, files::Vector{Vector{OutputFile}}, batch::Batch, w::Workspace)
    st = batch_stats(cfg)
    counts = zeros(Int, length(files))
    for i in 1:batch.n
        o = process(cfg, batch.r1[i], cfg.paired ? batch.r2[i] : EMPTY_READ, w.align, st)
        o.dest == 0 && continue
        format_read!(w.bufs[o.dest][1], o.r1, files[o.dest][1].fastq)
        cfg.paired && format_read!(w.bufs[o.dest][2], o.r2, files[o.dest][2].fastq)
        counts[o.dest] += 1
        st.reads_out += 1
        st.bp_out += length(o.r1) + length(o.r2)
    end
    blocks = [[encode_block!(files[d][m], w.bufs[d][m]) for m in eachindex(w.bufs[d])] for d in eachindex(files)]
    return BatchResult(blocks, counts, st)
end

function next_batch!(batch::Batch, cfg::Config, sources, batch_size::Int)
    batch.n = read_batch!(batch.r1, sources[1], batch_size)
    if cfg.paired
        n2 = read_batch!(batch.r2, sources[2], batch_size)
        batch.n == n2 || throw(ArgumentError("paired input files have different numbers of reads"))
        for i in 1:batch.n
            n1, n2 = batch.r1[i].name, batch.r2[i].name
            record_names_match(n1, n2) || throw(ArgumentError(
                "Reads are improperly paired. Read name $(repr(String(copy(n1)))) in file 1 " *
                "does not match $(repr(String(copy(n2)))) in file 2."))
        end
    end
    return batch.n
end

function write_result!(files, stats::TrimStats, written::Vector{Int}, r::BatchResult)
    for d in eachindex(files), m in eachindex(r.blocks[d])
        write_block(files[d][m], r.blocks[d][m])
    end
    written .+= r.counts
    merge!(stats, r.stats)
    stats.reads_out += r.stats.reads_out
    stats.bp_out += r.stats.bp_out
    return nothing
end

function run_serial!(cfg::Config, sources, files, stats::TrimStats, written::Vector{Int}, batch_size::Int)
    w = Workspace(length(files), cfg.paired ? 2 : 1)
    batch = Batch()
    while next_batch!(batch, cfg, sources, batch_size) > 0
        write_result!(files, stats, written, process_batch(cfg, files, batch, w))
    end
    return nothing
end

# Reading, processing and writing overlap: a producer task parses batches
# (gzip inputs are decompressed by a further task), worker tasks process them,
# results are written in input order by the calling task. Julia tasks are not
# preemptive, so with three or more threads one is left to the reading tasks,
# which would otherwise wait for a worker to finish its batch.
function run_batches!(cfg::Config, sources, files, stats::TrimStats, written::Vector{Int},
                      threads::Int, batch_size::Int)
    threads == 1 && return run_serial!(cfg, sources, files, stats, written, batch_size)
    workers = threads >= 3 ? threads - 1 : threads
    nmates = cfg.paired ? 2 : 1
    # Batches being parsed, queued or processed; recycled afterwards.
    pool = Channel{Batch}(2 * workers + 2)
    foreach(_ -> put!(pool, Batch()), 1:(2*workers+2))
    workspaces = Channel{Workspace}(workers)
    foreach(_ -> put!(workspaces, Workspace(length(files), nmates)), 1:workers)
    producer = Ref{Task}()
    batches = Channel{Batch}(workers; spawn = true, taskref = producer) do ch
        while true
            batch = take!(pool)
            next_batch!(batch, cfg, sources, batch_size) == 0 && break
            put!(ch, batch)
        end
    end
    pending = Task[]
    try
        for batch in batches
            push!(pending, Threads.@spawn begin
                local w = take!(workspaces)
                try
                    process_batch(cfg, files, batch, w)
                finally
                    put!(workspaces, w)
                    put!(pool, batch)
                end
            end)
            if length(pending) >= workers
                write_result!(files, stats, written, fetch(popfirst!(pending)))
            end
        end
        while !isempty(pending)
            write_result!(files, stats, written, fetch(popfirst!(pending)))
        end
    catch e
        # report the error of the failed task itself (e.g. a malformed input)
        while e isa TaskFailedException
            e = e.task.exception
        end
        throw(e)
    finally
        # On failure, stop the producer and wait for every task before
        # returning, so that no task still uses the streams or the outputs.
        close(batches)
        close(pool)
        for t in pending
            try
                wait(t)
            catch
            end
        end
        try
            wait(producer[])
        catch
        end
    end
    return nothing
end

"""
    remove_primers(input, output; adapters, kwargs...) -> TrimStats

Remove adapters or primers from the reads of `input` (FASTQ or FASTA, plain
or gzip) and write them to `output` (gzip when the name ends in `.gz`; FASTA
when it ends in `.fasta`, `.fa`, `.fna` or `.fas`). For paired-end data give
tuples: `remove_primers(("R1.fq.gz", "R2.fq.gz"), ("o1.fq.gz", "o2.fq.gz"); ...)`.
Outputs are written to temporary files and renamed into place at the end;
existing files are kept unless `overwrite=true`.

Adapters (see [`front_adapter`](@ref), [`back_adapter`](@ref),
[`adapters_from_fasta`](@ref)):

- `adapters`: searched in R1 (cutadapt `-a`/`-g`); `adapters2`: in R2 (`-A`/`-G`).
- `times=1`: removal rounds (`-n`); `action=:trim` (`:retain`, `:mask`,
  `:lowercase`, `:none`) as `--action`.
- `revcomp=false`: also try the reverse complement and keep the better
  orientation (`--revcomp`); `rc_suffix=" rc"` is appended to the names of
  reversed reads. For pairs, R1 and R2 are swapped instead.
- `pair_adapters=false`: `adapters[i]` and `adapters2[i]` must be found together (`--pair-adapters`).

Other modifications, applied before the adapters: `cut`/`cut2` (`-u`/`-U`, an
integer or a vector of up to two), `quality_cutoff`/`quality_cutoff2`
(`-q`/`-Q`, `20` or `(front, back)`), `quality_base=33`.

Filters, applied after trimming in this order: `min_length` (`-m`), `max_length`
(`-M`) — integers, or `(R1, R2)` tuples with `nothing` for no limit —,
`max_n` (`--max-n`, a count or a fraction below 1), `max_expected_errors`
(`--max-ee`), then one of `discard_trimmed`, `discard_untrimmed` or
`untrimmed_output` (path, or tuple for pairs). `pair_filter` (`:any`, default,
`:both`, `:first`) decides when a pair is filtered; for untrimmed pairs it is
`:both` when adapters are searched in only one mate, as in cutadapt.

`threads` (default `Threads.nthreads()`), `batch_size=4096` reads per batch,
`compression_level=1` for gzip outputs (libdeflate, 0–12).
"""
remove_primers(input, output; kwargs...) = run_pipeline(input, output; demux = false, kwargs...)

"""
    demultiplex(input, output; adapters, kwargs...) -> TrimStats

Like [`remove_primers`](@ref), but every read goes to the file of the adapter
that matched it (in the last removal round): `output` must contain `{name}`,
replaced by the adapter name. Reads without a match go to `{name}` =
`"unknown"`, to `untrimmed_output` if given, or are dropped with
`discard_untrimmed=true`. One file is created for every adapter name, also
when it receives no reads.

For paired-end data `output` is a tuple. With `{name}` the R1 match decides;
with `{name1}` and `{name2}` in both paths, pairs are routed by the
combination of the R1 (`adapters`) and R2 (`adapters2`) matches, `"unknown"`
standing for a missing match (combinatorial demultiplexing).
"""
demultiplex(input, output; kwargs...) = run_pipeline(input, output; demux = true, kwargs...)
