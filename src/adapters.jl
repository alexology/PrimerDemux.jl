# Adapter and primer types, matching rules and the cutadapt specification
# syntax ("name=^ACGT...TTGC$;e=0.2;o=5").

"""
    Placement

Where a single adapter is searched (cutadapt option and syntax in brackets):

- `FRONT`: 5' adapter anywhere in the read, also partial at the 5' end (`-g ADAPTER`);
  the adapter and everything before it are removed.
- `BACK`: 3' adapter anywhere in the read, also partial at the 3' end (`-a ADAPTER`);
  the adapter and everything after it are removed.
- `PREFIX`: anchored 5' adapter, must be at the very start and complete (`-g ^ADAPTER`).
- `SUFFIX`: anchored 3' adapter, must be at the very end and complete (`-a ADAPTER\$`).
- `FRONT_NOT_INTERNAL`: 5' adapter that may be partial but must start the read (`-g XADAPTER`).
- `BACK_NOT_INTERNAL`: 3' adapter that may be partial but must end the read (`-a ADAPTERX`).
"""
@enum Placement FRONT BACK PREFIX SUFFIX FRONT_NOT_INTERNAL BACK_NOT_INTERNAL

removes_before(p::Placement) = p in (FRONT, PREFIX, FRONT_NOT_INTERNAL)

abstract type Adapter end

"""
    SingleAdapter

One adapter or primer sequence with its search parameters. Build it with
[`front_adapter`](@ref) or [`back_adapter`](@ref).
"""
struct SingleAdapter <: Adapter
    name::String
    sequence::String
    placement::Placement
    max_error_rate::Float64
    min_overlap::Int
    read_wildcards::Bool
    adapter_wildcards::Bool
    indels::Bool
    locator::Union{Aligner,AffixComparer}
    kmer_finder::KmerFinder
end

"""
    LinkedAdapter

A 5' adapter followed by a 3' adapter (`FRONT...BACK`): the 5' part is searched
first, the 3' part in what remains after it. Both parts are removed. Build it
with [`front_adapter`](@ref) or [`back_adapter`](@ref) and a `...` in the
specification.
"""
struct LinkedAdapter <: Adapter
    name::String
    front::SingleAdapter
    back::SingleAdapter
    front_required::Bool
    back_required::Bool
end

adapter_name(a::Adapter) = a.name

const IUPAC_CHARS = Set("ABCDGHKMNRSTUVWXY")

function SingleAdapter(sequence::AbstractString, placement::Placement;
                       name::AbstractString = "", max_errors::Real = 0.1,
                       min_overlap::Integer = 3, read_wildcards::Bool = false,
                       adapter_wildcards::Bool = true, indels::Bool = true)
    seq = replace(uppercase(String(sequence)), 'U' => 'T', 'I' => 'N')
    isempty(seq) && throw(ArgumentError("Adapter sequence is empty"))
    isascii(seq) || throw(ArgumentError("Adapter sequence must be ASCII: $(repr(seq))"))
    max_errors >= 0 || throw(ArgumentError("max_errors must not be negative"))
    min_overlap >= 1 || throw(ArgumentError("min_overlap must be at least 1"))
    n_count = count(==('N'), seq)
    rate = Float64(max_errors)
    if max_errors >= 1 && n_count != length(seq)
        rate = max_errors / (length(seq) - n_count)
    end
    overlap = min(Int(min_overlap), length(seq))
    if adapter_wildcards
        for c in seq
            c in IUPAC_CHARS || throw(ArgumentError(
                "Character '$c' in adapter sequence '$seq' is not a valid IUPAC code. " *
                "Use only characters 'ABCDGHIKMNRSTUVWXY'."))
        end
    end
    # Plain matching is used when the adapter contains only ACGT.
    wildcards = adapter_wildcards && !all(in("ACGT"), seq)
    anchored = placement in (PREFIX, SUFFIX)
    anchored && (overlap = length(seq))

    if anchored && !indels
        locator = AffixComparer(seq, rate; suffix = placement == SUFFIX, wildcard_ref = wildcards,
                                wildcard_query = read_wildcards, min_overlap = overlap)
        # No pre-filter: the comparison does not build a DP matrix.
        kmer_finder = KmerFinder(KmerSearch[], true)
    else
        flags = placement == FRONT ? QUERY_START | QUERY_STOP | REFERENCE_START :
                placement == BACK ? QUERY_START | QUERY_STOP | REFERENCE_END :
                placement == PREFIX ? QUERY_STOP :
                placement == SUFFIX ? QUERY_START :
                placement == FRONT_NOT_INTERNAL ? REFERENCE_START | QUERY_STOP :
                QUERY_START | REFERENCE_END
        locator = Aligner(seq, rate; flags, wildcard_ref = wildcards, wildcard_query = read_wildcards,
                          indel_cost = indels ? 1 : 100_000, min_overlap = overlap)
        back_kmers = placement in (BACK, SUFFIX, BACK_NOT_INTERNAL)
        internal = placement in (FRONT, BACK)
        positions = create_positions_and_kmers(seq, overlap, rate, back_kmers, !back_kmers, internal)
        kmer_finder = KmerFinder(positions, wildcards, read_wildcards)
    end
    return SingleAdapter(String(name), seq, placement, rate, overlap, read_wildcards,
                         wildcards, indels, locator, kmer_finder)
end

Base.length(a::SingleAdapter) = length(a.sequence)

function spec(a::SingleAdapter)
    a.placement == FRONT && return a.sequence * "..."
    a.placement == BACK && return a.sequence
    a.placement == PREFIX && return "^" * a.sequence * "..."
    a.placement == SUFFIX && return a.sequence * "\$"
    a.placement == FRONT_NOT_INTERNAL && return "X" * a.sequence * "..."
    return a.sequence * "X"
end
spec(a::LinkedAdapter) = spec(a.front) * "..." * spec(a.back)

Base.show(io::IO, a::SingleAdapter) =
    print(io, "SingleAdapter(", repr(a.name), ", ", repr(spec(a)), ", ", a.placement,
          ", max_error_rate=", a.max_error_rate, ", min_overlap=", a.min_overlap, ")")
Base.show(io::IO, a::LinkedAdapter) =
    print(io, "LinkedAdapter(", repr(a.name), ", ", repr(spec(a)), ", front_required=",
          a.front_required, ", back_required=", a.back_required, ")")

# ---------------------------------------------------------------------------
# Matches

# Alignment of one adapter part on the sequence it was searched in.
struct PartMatch
    astart::Int
    astop::Int
    rstart::Int
    rstop::Int
    score::Int
    errors::Int
    seqlen::Int                 # length of the searched sequence
end

"""
    AdapterMatch

Result of matching an adapter to a read. `adapter` is the index of the adapter
in the searched list. `front`/`back` are the parts that removed sequence before
or after the match (`nothing` when absent). Coordinates are 0-based and
half-open on the sequence the part was searched in (for the 3' part of a
linked adapter, the read after removal of the 5' part).
"""
struct AdapterMatch
    adapter::Int
    front::Union{Nothing,PartMatch}
    back::Union{Nothing,PartMatch}
end

score(m::AdapterMatch) = (m.front === nothing ? 0 : m.front.score) + (m.back === nothing ? 0 : m.back.score)
errors(m::AdapterMatch) = (m.front === nothing ? 0 : m.front.errors) + (m.back === nothing ? 0 : m.back.errors)

# Part of the matched sequence that remains after trimming (0-based, half-open).
function remainder_interval(m::AdapterMatch, seqlen::Int)
    start = m.front === nothing ? 0 : m.front.rstop
    stop = m.back === nothing ? seqlen : start + m.back.rstart
    return start, stop
end

# Part kept by action :retain (the adapter stays in the read).
function retained_adapter_interval(m::AdapterMatch, seqlen::Int)
    if m.front !== nothing && m.back === nothing
        return m.front.rstart, seqlen
    elseif m.front === nothing
        return 0, m.back.rstop
    end
    return m.front.rstart, m.back.rstop + m.front.rstop
end

function match_part(a::SingleAdapter, seq::AbstractVector{UInt8}, ws::AlignWorkspace)
    kmers_present(a.kmer_finder, seq) || return nothing
    al = locate(a.locator, seq, ws)
    found(al) || return nothing
    return PartMatch(al.ref_start, al.ref_stop, al.query_start, al.query_stop, al.score,
                     al.errors, length(seq))
end

function match_to(a::SingleAdapter, index::Int, seq::AbstractVector{UInt8}, ws::AlignWorkspace)
    part = match_part(a, seq, ws)
    part === nothing && return nothing
    return removes_before(a.placement) ? AdapterMatch(index, part, nothing) :
           AdapterMatch(index, nothing, part)
end

function match_to(a::LinkedAdapter, index::Int, seq::AbstractVector{UInt8}, ws::AlignWorkspace)
    front = match_part(a.front, seq, ws)
    a.front_required && front === nothing && return nothing
    rest = front === nothing ? seq : view(seq, (front.rstop+1):length(seq))
    back = match_part(a.back, rest, ws)
    back === nothing && (a.back_required || front === nothing) && return nothing
    return AdapterMatch(index, front, back)
end

# Best match among several adapters: highest score, then fewest errors, then
# the first adapter in the list.
function best_match(adapters::Vector{<:Adapter}, seq::AbstractVector{UInt8}, ws::AlignWorkspace)
    best = nothing
    for (i, a) in enumerate(adapters)
        m = match_to(a, i, seq, ws)
        m === nothing && continue
        if best === nothing || score(m) > score(best) ||
           (score(m) == score(best) && errors(m) < errors(best))
            best = m
        end
    end
    return best
end

# ---------------------------------------------------------------------------
# Specification syntax

const PARAMETER_ALIASES = Dict(
    "e" => "max_errors", "error_rate" => "max_errors", "max_error_rate" => "max_errors",
    "max_errors" => "max_errors", "o" => "min_overlap", "min_overlap" => "min_overlap",
    "required" => "required", "optional" => "optional", "indels" => "indels",
    "noindels" => "noindels")

function parse_search_parameters(spec::AbstractString)
    result = Dict{String,Any}()
    for field in split(spec, ';')
        field = strip(field)
        isempty(field) && continue
        key, eq, value = partition(field, "=")
        key = strip(key)
        haskey(PARAMETER_ALIASES, key) || throw(ArgumentError("Unknown parameter '$key'"))
        eq == "=" && isempty(strip(value)) && throw(ArgumentError("No value given for key '$key'"))
        key = PARAMETER_ALIASES[key]
        value = strip(value)
        parsed = isempty(value) ? true : something(tryparse(Int, value), tryparse(Float64, value), value)
        parsed isa AbstractString && throw(ArgumentError("Invalid value '$value' for '$key'"))
        haskey(result, key) && throw(ArgumentError("Key '$key' specified twice"))
        result[key] = parsed
    end
    haskey(result, "optional") && haskey(result, "required") &&
        throw(ArgumentError("'optional' and 'required' cannot be specified at the same time"))
    haskey(result, "indels") && haskey(result, "noindels") &&
        throw(ArgumentError("'indels' and 'noindels' cannot be specified at the same time"))
    haskey(result, "optional") && (delete!(result, "optional"); result["required"] = false)
    haskey(result, "noindels") && (delete!(result, "noindels"); result["indels"] = false)
    return result
end

function partition(s::AbstractString, sep::AbstractString)
    r = findfirst(sep, s)
    r === nothing && return String(s), "", ""
    return String(s[1:prevind(s, first(r))]), String(sep), String(s[nextind(s, last(r)):end])
end

# "TGA{5}CT" -> "TGAAAAACT"
function expand_braces(sequence::AbstractString)
    out = Char[]
    chars = collect(sequence)
    i = 1
    repeatable = false          # a plain character precedes the current position
    while i <= length(chars)
        c = chars[i]
        if c == '{'
            repeatable || throw(ArgumentError("\"{\" must be used after a character"))
            j = findnext(==('}'), chars, i + 1)
            j === nothing && throw(ArgumentError("Unterminated expression"))
            n = tryparse(Int, String(chars[(i+1):(j-1)]))
            (n === nothing || !(0 <= n <= 10000)) &&
                throw(ArgumentError("Invalid repeat count in $(repr(sequence))"))
            prev = pop!(out)
            append!(out, fill(prev, n))
            repeatable = false
            i = j + 1
        elseif c == '}'
            throw(ArgumentError("\"}\" cannot be used here"))
        else
            push!(out, c)
            repeatable = true
            i += 1
        end
    end
    return String(out)
end

struct ParsedSpec
    name::Union{Nothing,String}
    restriction::Symbol         # :none, :anchored, :noninternal
    sequence::String
    parameters::Dict{String,Any}
end

function parse_single_spec(spec::AbstractString, type::Symbol)
    spec, _, parameters_spec = partition(spec, ";")
    fields = split(spec, '='; limit = 2)
    name = nothing
    if length(fields) > 1
        name = String(strip(fields[1]))
        spec = fields[2]
    end
    spec = expand_braces(strip(spec))
    parameters = parse_search_parameters(parameters_spec)
    front_restriction = :none
    if startswith(spec, "^")
        front_restriction = :anchored
        spec = spec[2:end]
    end
    if startswith(uppercase(spec), "X")
        front_restriction == :none || throw(ArgumentError("two front restrictions in $(repr(spec))"))
        front_restriction = :noninternal
        spec = lstrip(c -> c in "xX", spec)
    end
    back_restriction = :none
    if endswith(spec, "\$")
        back_restriction = :anchored
        spec = spec[1:end-1]
    end
    if endswith(uppercase(spec), "X")
        back_restriction == :none || throw(ArgumentError("two back restrictions in $(repr(spec))"))
        back_restriction = :noninternal
        spec = rstrip(c -> c in "xX", spec)
    end
    front_restriction != :none && back_restriction != :none && throw(ArgumentError(
        "You cannot use multiple placement restrictions for an adapter at the same time. " *
        "Choose one of ^ADAPTER, ADAPTER\$, XADAPTER or ADAPTERX"))
    type == :front && back_restriction != :none &&
        throw(ArgumentError("Allowed placement restrictions for a 5' adapter are XADAPTER and ^ADAPTER"))
    type == :back && front_restriction != :none &&
        throw(ArgumentError("Allowed placement restrictions for a 3' adapter are ADAPTERX and ADAPTER\$"))
    restriction = front_restriction != :none ? front_restriction : back_restriction
    restriction == :anchored && haskey(parameters, "min_overlap") && throw(ArgumentError(
        "Setting 'min_overlap=' (or 'o=') for anchored adapters is not possible because " *
        "anchored adapters always need to match in full."))
    if get(parameters, "min_overlap", 0) > length(spec)
        parameters["min_overlap"] = length(spec)
    end
    return ParsedSpec(name, restriction, String(spec), parameters)
end

function placement_of(type::Symbol, restriction::Symbol)
    if type == :front
        return restriction == :none ? FRONT : restriction == :anchored ? PREFIX : FRONT_NOT_INTERNAL
    end
    return restriction == :none ? BACK : restriction == :anchored ? SUFFIX : BACK_NOT_INTERNAL
end

function build_single(p::ParsedSpec, type::Symbol, defaults::Dict{String,Any})
    params = merge(defaults, p.parameters)
    haskey(params, "required") && throw(ArgumentError("'required' and 'optional' can only be used within linked adapters"))
    return SingleAdapter(p.sequence, placement_of(type, p.restriction);
                         name = something(p.name, p.sequence),
                         max_errors = params["max_errors"], min_overlap = params["min_overlap"],
                         read_wildcards = params["read_wildcards"],
                         adapter_wildcards = params["adapter_wildcards"], indels = params["indels"])
end

function make_adapter(spec::AbstractString, type::Symbol, defaults::Dict{String,Any},
                      name::Union{Nothing,AbstractString})
    type in (:front, :back) || throw(ArgumentError("type must be :front or :back"))
    spec1, middle, spec2 = partition(spec, "...")
    if middle == "..." && !isempty(spec1) && !isempty(spec2)
        front = parse_single_spec(spec1, :front)
        back = parse_single_spec(spec2, :back)
        linked_name = something(name, front.name, spec)
        fparams = merge(defaults, front.parameters)
        bparams = merge(defaults, back.parameters)
        # -g requires both parts; -a only the anchored ones.
        front_required = type == :front ? true : front.restriction != :none
        back_required = type == :front ? true : back.restriction != :none
        front_required = pop!(fparams, "required", front_required)
        back_required = pop!(bparams, "required", back_required)
        fa = build_single(ParsedSpec("linked_front", front.restriction, front.sequence, fparams), :front, Dict{String,Any}())
        ba = build_single(ParsedSpec("linked_back", back.restriction, back.sequence, bparams), :back, Dict{String,Any}())
        return LinkedAdapter(String(linked_name), fa, ba, front_required, back_required)
    end
    if middle == "..."
        if isempty(spec1)
            type == :front && throw(ArgumentError("Invalid adapter specification"))
            spec = spec2         # -a ...ADAPTER: a plain 3' adapter
        else
            spec = spec1         # -a ADAPTER... is a 5' adapter
            type = :front
        end
    end
    p = parse_single_spec(spec, type)
    name === nothing || (p = ParsedSpec(String(name), p.restriction, p.sequence, p.parameters))
    return build_single(p, type, defaults)
end

function search_defaults(; max_errors::Real = 0.1, min_overlap::Integer = 3,
                         read_wildcards::Bool = false, adapter_wildcards::Bool = true,
                         indels::Bool = true)
    return Dict{String,Any}("max_errors" => max_errors, "min_overlap" => min_overlap,
                            "read_wildcards" => read_wildcards,
                            "adapter_wildcards" => adapter_wildcards, "indels" => indels)
end

"""
    front_adapter(spec; name, max_errors=0.1, min_overlap=3, indels=true,
                  adapter_wildcards=true, read_wildcards=false) -> Adapter

5' adapter or primer, as cutadapt `-g spec`. The specification follows the
cutadapt syntax:

- `"ACGT"` regular 5' adapter, `"^ACGT"` anchored, `"XACGT"` non-internal;
- `"name=ACGT"` gives the name used in reports and for demultiplexing
  (default: the specification itself);
- `"ACGT...TTGC"` linked adapter (5' part then 3' part, both required with
  `front_adapter`), e.g. `"^FWDPRIMER...REVPRIMER_RC\$"`;
- `";e=0.2;o=5;noindels"` and, in linked parts, `";required"`/`";optional"`
  override the keyword defaults for one adapter;
- `"A{5}"` expands to `"AAAAA"`; IUPAC codes are wildcards in the adapter.

`max_errors` is an error rate below 1 and an absolute number of errors from 1
upwards, as in cutadapt `-e`.
"""
front_adapter(spec::AbstractString; name::Union{Nothing,AbstractString} = nothing, kwargs...) =
    make_adapter(spec, :front, search_defaults(; kwargs...), name)

"""
    back_adapter(spec; kwargs...) -> Adapter

3' adapter or primer, as cutadapt `-a spec`: `"ACGT"` regular, `"ACGT\$"`
anchored, `"ACGTX"` non-internal. In a linked adapter (`"ACGT...TTGC"`) only
anchored parts are required by default. See [`front_adapter`](@ref) for the
syntax and the keywords.
"""
back_adapter(spec::AbstractString; name::Union{Nothing,AbstractString} = nothing, kwargs...) =
    make_adapter(spec, :back, search_defaults(; kwargs...), name)

"""
    adapters_from_fasta(path; type=:front, anchored=false, kwargs...) -> Vector{Adapter}

One adapter per FASTA record, named after the record identifier, as cutadapt
`-g file:path` (`type=:front`) or `-a file:path` (`type=:back`).
`anchored=true` corresponds to `-g ^file:path` and `-a file\$:path`. Records may
use the full specification syntax (for example linked `FWD...REV`). Typical
use: a FASTA of sample barcodes or primer pairs for [`demultiplex`](@ref).
"""
function adapters_from_fasta(path::AbstractString; type::Symbol = :front, anchored::Bool = false, kwargs...)
    type in (:front, :back) || throw(ArgumentError("type must be :front or :back"))
    defaults = search_defaults(; kwargs...)
    result = Adapter[]
    open(path, "r") do raw
        # gzip detected from the magic bytes, as for the reads
        gzipped = read(raw, 2) == UInt8[0x1f, 0x8b]
        seekstart(raw)
        stream = gzipped ? GzipDecompressorStream(raw) : raw
        reader = FASTA.Reader(stream)
        try
            for record in reader
                id = String(identifier(record))
                seq = String(sequence(record))
                spec = anchored ? (type == :front ? "^" * seq : seq * "\$") : seq
                push!(result, make_adapter(spec, type, defaults, isempty(id) ? nothing : id))
            end
        finally
            close(reader)
        end
    end
    return result
end
