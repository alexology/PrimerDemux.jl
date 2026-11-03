# Reading FASTQ/FASTA (plain or gzip) with FASTX and writing outputs through
# temporary files that are renamed into place only when all reads are done.

struct Read
    name::Vector{UInt8}         # full header line without '@' or '>'
    seq::Vector{UInt8}
    qual::Vector{UInt8}         # raw ASCII qualities; empty for FASTA
end

Read() = Read(UInt8[], UInt8[], UInt8[])

Base.length(r::Read) = length(r.seq)

# 0-based half-open slice, as Python's read[start:stop].
function slice(r::Read, start::Int, stop::Int, has_qual::Bool)
    start = clamp(start, 0, length(r.seq))
    stop = clamp(stop, start, length(r.seq))
    (start == 0 && stop == length(r.seq)) && return r
    return Read(r.name, r.seq[(start+1):stop], has_qual ? r.qual[(start+1):stop] : r.qual)
end

with_suffix(r::Read, suffix::Vector{UInt8}) = Read(vcat(r.name, suffix), r.seq, r.qual)

# Reverse complement into `dest` (its vectors are overwritten and reused).
function reverse_complement!(dest::Read, r::Read, has_qual::Bool)
    n = length(r.seq)
    resize!(dest.seq, n)
    @inbounds for i in 1:n
        dest.seq[n-i+1] = COMPLEMENT_TABLE[r.seq[i]+1]
    end
    if has_qual
        resize!(dest.qual, n)
        @inbounds for i in 1:n
            dest.qual[n-i+1] = r.qual[i]
        end
    end
    return Read(r.name, dest.seq, has_qual ? dest.qual : r.qual)
end

reverse_complement(r::Read, has_qual::Bool) = reverse_complement!(Read(), r, has_qual)

is_gzip_path(path::AbstractString) = endswith(lowercase(path), ".gz")

function is_fasta_path(path::AbstractString)
    p = lowercase(path)
    endswith(p, ".gz") && (p = p[1:end-3])
    return any(ext -> endswith(p, ext), (".fasta", ".fa", ".fna", ".fas"))
end

# --- inputs ------------------------------------------------------------------

# An opened input. For gzip files a separate task decompresses into `pipe`
# while FASTX parses from it, so decompression and parsing overlap.
mutable struct InputSource
    raw::IOStream
    stream::IO                  # what FASTX reads
    pipe::Union{Nothing,Base.BufferStream}
    inflater::Union{Nothing,Task}
    reader::Any                 # FASTQ.Reader, FASTA.Reader or nothing (empty file)
    fastq::Bool
end

# Decompressed bytes waiting for the parser before the inflater pauses.
const PIPE_LIMIT = 32 << 20

function inflate_into(pipe::Base.BufferStream, raw::IOStream)
    stream = GzipDecompressorStream(raw)
    buf = Vector{UInt8}(undef, 1 << 20)
    try
        while !eof(stream)
            n = readbytes!(stream, buf)
            while bytesavailable(pipe) > PIPE_LIMIT && isopen(pipe)
                sleep(0.001)
            end
            isopen(pipe) || return nothing      # the reader stopped
            write(pipe, view(buf, 1:n))
        end
    catch
        isopen(pipe) && rethrow()               # a failure, not a cancellation
    finally
        close(pipe)                             # end of data for the parser
        TranscodingStreams.finalize(stream.codec)
    end
    return nothing
end

# Opens `path`, detecting gzip from the magic bytes and FASTA/FASTQ from the
# first character.
function open_input(path::AbstractString; threaded::Bool = true)
    raw = open(path, "r")
    src = nothing
    try
        magic = read(raw, 2)
        seekstart(raw)
        if magic == UInt8[0x1f, 0x8b] && threaded
            pipe = Base.BufferStream()
            task = Threads.@spawn inflate_into(pipe, raw)
            src = InputSource(raw, NoopStream(pipe), pipe, task, nothing, true)
        else
            stream = magic == UInt8[0x1f, 0x8b] ? GzipDecompressorStream(raw) : NoopStream(raw)
            src = InputSource(raw, stream, nothing, nothing, nothing, true)
        end
        stream = src.stream
        if eof(stream)
            check_input(src)
            return src
        end
        first_byte = peek(stream, UInt8)
        if first_byte == UInt8('@')
            src.reader = FASTQ.Reader(stream; copy = false)
        elseif first_byte == UInt8('>')
            src.reader = FASTA.Reader(stream; copy = false)
            src.fastq = false
        else
            throw(ArgumentError("$(repr(path)) is neither FASTQ ('@') nor FASTA ('>')"))
        end
        return src
    catch
        src === nothing ? close(raw) : close_input(src)
        rethrow()
    end
end

# At the end of the parsed data: reports a decompression error, if any.
function check_input(src::InputSource)
    src.inflater === nothing && return nothing
    try
        wait(src.inflater)
    catch e
        e isa TaskFailedException ? throw(e.task.exception) : rethrow()
    end
    return nothing
end

function close_input(src::InputSource)
    if src.pipe !== nothing
        close(src.pipe)                         # stops the inflater
        try
            wait(src.inflater)
        catch
        end
    else
        try
            close(src.stream)
        catch
        end
    end
    close(src.raw)
    return nothing
end

function copy_bytes!(dest::Vector{UInt8}, src)
    resize!(dest, length(src))
    copyto!(dest, src)
    return dest
end

function fill_read!(r::Read, rec::FASTQ.Record)
    copy_bytes!(r.name, codeunits(description(rec)))
    copy_bytes!(r.seq, codeunits(sequence(StringView, rec)))
    copy_bytes!(r.qual, codeunits(quality(StringView, rec)))
    return r
end

function fill_read!(r::Read, rec::FASTA.Record)
    copy_bytes!(r.name, codeunits(description(rec)))
    copy_bytes!(r.seq, codeunits(sequence(StringView, rec)))
    empty!(r.qual)
    return r
end

# Reads up to `n` records into `reads`, reusing the vectors of earlier
# batches; returns the number read.
function read_batch!(reads::Vector{Read}, src::InputSource, n::Int)
    src.reader === nothing && return 0
    k = 0
    while k < n
        # FASTX readers keep their position internally.
        it = iterate(src.reader)
        if it === nothing
            check_input(src)
            break
        end
        k += 1
        k > length(reads) && push!(reads, Read())
        fill_read!(reads[k], it[1])
    end
    return k
end

# Whether two mates belong to the same pair, as dnaio's record_names_match
# (used by cutadapt): the IDs (up to the first space or tab) must be equal,
# ignoring a last character 1, 2 or 3 when both IDs end with one.
function record_names_match(h1::Vector{UInt8}, h2::Vector{UInt8})
    n = something(findfirst(b -> b == UInt8(' ') || b == UInt8('\t'), h2), length(h2) + 1) - 1
    length(h1) < n && return false
    length(h1) > n && h1[n+1] != UInt8(' ') && h1[n+1] != UInt8('\t') && return false
    n == 0 && return true
    UInt8('1') <= h1[n] <= UInt8('3') && UInt8('1') <= h2[n] <= UInt8('3') && (n -= 1)
    return view(h1, 1:n) == view(h2, 1:n)
end

# --- outputs -----------------------------------------------------------------

# Records are formatted and, for .gz names, compressed per batch in the worker
# tasks; each block is a complete gzip member (a multi-member file, like the
# output of parallel gzip tools), so the writer only appends bytes.
mutable struct OutputFile
    path::String
    temp::String
    io::IOStream
    fastq::Bool
    gzip::Bool
    level::Int
    blocks::Int                 # non-empty blocks written
end

# Keys identifying the file at `path`, to refuse outputs that would overwrite
# each other or an input. The name, with the real path of an existing parent,
# is compared ignoring case on Windows and macOS, whose usual file systems are
# case-insensitive (on a case-sensitive one a false match only gives an error);
# existing files are also identified by device and inode (links, short names).
function path_keys(path::AbstractString)
    full = normpath(abspath(path))
    parent, name = dirname(full), basename(full)
    isdir(parent) && (full = joinpath(realpath(parent), name))
    key = Sys.iswindows() || Sys.isapple() ? lowercase(full) : full
    st = stat(path)
    return ispath(st) ? Any[key, (st.device, st.inode)] : Any[key]
end

function new_temp(dir::AbstractString)
    flags = Base.Filesystem.JL_O_WRONLY | Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL
    while true
        temp = tempname(dir; cleanup = false)
        try
            close(Base.Filesystem.open(temp, flags, 0o666))
            return temp
        catch e
            e isa Base.IOError && e.code == Base.UV_EEXIST || rethrow()
        end
    end
end

function OutputFile(path::AbstractString, input_fastq::Bool, compression_level::Int)
    0 <= compression_level <= 12 || throw(ArgumentError("compression_level must be within 0..12"))
    dir = dirname(abspath(path))
    isdir(dir) || mkpath(dir)
    temp = new_temp(dir)
    io = try
        open(temp, "w")
    catch
        rm(temp; force = true)
        rethrow()
    end
    # A FASTA name turns FASTQ input into FASTA output, as in cutadapt.
    return OutputFile(String(path), temp, io, input_fastq && !is_fasta_path(path),
                      is_gzip_path(path), compression_level, 0)
end

function format_read!(buf::Vector{UInt8}, r::Read, fastq::Bool)
    push!(buf, fastq ? UInt8('@') : UInt8('>'))
    append!(buf, r.name)
    push!(buf, UInt8('\n'))
    append!(buf, r.seq)
    if fastq
        append!(buf, (UInt8('\n'), UInt8('+'), UInt8('\n')))
        append!(buf, r.qual)
    end
    push!(buf, UInt8('\n'))
    return nothing
end

# Whole-block gzip compression with libdeflate (faster than zlib at similar
# ratios; blocks are complete in memory, so no streaming is needed).
function gzip_member(data::Vector{UInt8}, level::Int)
    c = ccall((:libdeflate_alloc_compressor, libdeflate_jll.libdeflate), Ptr{Cvoid}, (Cint,), level)
    c == C_NULL && throw(OutOfMemoryError())
    try
        bound = ccall((:libdeflate_gzip_compress_bound, libdeflate_jll.libdeflate), Csize_t, (Ptr{Cvoid}, Csize_t), c, length(data))
        out = Vector{UInt8}(undef, bound)
        n = GC.@preserve data out ccall((:libdeflate_gzip_compress, libdeflate_jll.libdeflate), Csize_t,
                                        (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t),
                                        c, data, length(data), out, length(out))
        # cannot happen with the bound above
        n == 0 && error("libdeflate could not compress a block of $(length(data)) bytes")
        return resize!(out, n)
    finally
        ccall((:libdeflate_free_compressor, libdeflate_jll.libdeflate), Cvoid, (Ptr{Cvoid},), c)
    end
end

# Bytes to append for a block of formatted records (empty: nothing to write);
# `buf` is emptied and keeps its capacity for the next batch.
function encode_block!(out::OutputFile, buf::Vector{UInt8})
    block = isempty(buf) ? UInt8[] : out.gzip ? gzip_member(buf, out.level) : copy(buf)
    empty!(buf)
    return block
end

function write_block(out::OutputFile, block::Vector{UInt8})
    isempty(block) && return nothing
    write(out.io, block)
    out.blocks += 1
    return nothing
end

# A single rename replaces `dest` atomically or fails leaving it unchanged.
function rename_file(src::AbstractString, dest::AbstractString)
    err = ccall(:jl_fs_rename, Int32, (Cstring, Cstring), src, dest)
    err < 0 && Base.uv_error("rename($(repr(src)), $(repr(dest)))", err)
    return nothing
end

# Puts all outputs in place, or none of them: existing destinations are first
# moved aside, then the temporary files renamed; if a rename fails, the new
# files are removed and the old ones restored. A crash in between is not
# covered (several renames cannot be one atomic step).
function publish_all(outs::Vector{OutputFile})
    for out in outs
        # an empty gzip file still needs one (empty) member
        out.gzip && out.blocks == 0 && write(out.io, gzip_member(UInt8[], out.level))
        close(out.io)
    end
    backups = Pair{String,String}[]             # destination => backup
    placed = String[]
    try
        for out in outs
            ispath(out.path) || continue
            backup = new_temp(dirname(abspath(out.path)))
            try
                rename_file(out.path, backup)
            catch
                rm(backup; force = true)
                rethrow()
            end
            push!(backups, out.path => backup)
        end
        for out in outs
            rename_file(out.temp, out.path)
            push!(placed, out.path)
        end
    catch
        foreach(p -> rm(p; force = true), placed)
        for (path, backup) in backups
            try
                rename_file(backup, path)
            catch e
                @error "Could not restore $(repr(path)); its previous content is in $(repr(backup))" exception = e
            end
        end
        rethrow()
    end
    for (_, backup) in backups
        try
            rm(backup; force = true)
        catch e
            @warn "Could not remove the previous content of an output, kept in $(repr(backup))" exception = e
        end
    end
    return nothing
end

function discard(out::OutputFile)
    try
        close(out.io)
    catch
    end
    try
        rm(out.temp; force = true)
    catch e
        e isa Base.IOError || rethrow()
    end
    return nothing
end
