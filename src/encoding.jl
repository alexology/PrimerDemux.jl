# Byte translation tables used for matching, as in cutadapt's _match_tables.py.
#
# - UPPER_TABLE: ASCII uppercase; used when neither sequence allows wildcards.
# - ACGT_TABLE: A, C, G, T/U to one bit each, everything else to 0x80.
# - IUPAC_TABLE: IUPAC codes to the union of their bits; N also gets 0x80 so
#   that it matches the non-ACGT characters encoded by ACGT_TABLE.
# Two encoded characters x and y match when `x & y != 0`.

const UPPER_TABLE = let t = Vector{UInt8}(undef, 256)
    for b in 0x00:0xff
        t[b+1] = 0x61 <= b <= 0x7a ? b - 0x20 : b
    end
    t
end

const ACGT_TABLE = let t = fill(0x80, 256)
    for (c, v) in ('A' => 0x01, 'C' => 0x02, 'G' => 0x04, 'T' => 0x08, 'U' => 0x08)
        t[UInt8(c)+1] = v
        t[UInt8(lowercase(c))+1] = v
    end
    t
end

const IUPAC_TABLE = let t = zeros(UInt8, 256)
    A, C, G, T = 0x01, 0x02, 0x04, 0x08
    codes = ('X' => 0x00, 'A' => A, 'C' => C, 'G' => G, 'T' => T, 'U' => T,
             'R' => A | G, 'Y' => C | T, 'S' => G | C, 'W' => A | T,
             'K' => G | T, 'M' => A | C, 'B' => C | G | T, 'D' => A | G | T,
             'H' => A | C | T, 'V' => A | C | G, 'N' => (A | C | G | T) + 0x80)
    for (c, v) in codes
        t[UInt8(c)+1] = v
        t[UInt8(lowercase(c))+1] = v
    end
    t
end

# Tables for (reference, query) given the wildcard flags, and whether the
# comparison is plain equality of uppercased bytes.
function match_tables(wildcard_ref::Bool, wildcard_query::Bool)
    if wildcard_ref && wildcard_query
        return IUPAC_TABLE, IUPAC_TABLE, false
    elseif wildcard_ref
        return IUPAC_TABLE, ACGT_TABLE, false
    elseif wildcard_query
        return ACGT_TABLE, IUPAC_TABLE, false
    else
        return UPPER_TABLE, UPPER_TABLE, true
    end
end

@inline function translate!(dest::Vector{UInt8}, src::AbstractVector{UInt8}, table::Vector{UInt8})
    n = length(src)
    length(dest) < n && resize!(dest, n)
    @inbounds for j in 1:n
        dest[j] = table[src[j]+1]
    end
    return dest
end

translate(src::AbstractVector{UInt8}, table::Vector{UInt8}) = translate!(Vector{UInt8}(undef, length(src)), src, table)

# Complement of IUPAC nucleotide codes, preserving case; other bytes unchanged.
const COMPLEMENT_TABLE = let t = collect(0x00:0xff)
    pairs = ('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'U' => 'A',
             'R' => 'Y', 'Y' => 'R', 'S' => 'S', 'W' => 'W', 'K' => 'M',
             'M' => 'K', 'B' => 'V', 'V' => 'B', 'D' => 'H', 'H' => 'D', 'N' => 'N')
    for (a, b) in pairs
        t[UInt8(a)+1] = UInt8(b)
        t[UInt8(lowercase(a))+1] = UInt8(lowercase(b))
    end
    t
end

"""
    reverse_complement(seq::AbstractString) -> String

Reverse complement of a nucleotide sequence; IUPAC codes are complemented,
case is preserved and other characters are kept as they are.
"""
function reverse_complement(seq::AbstractString)
    s = codeunits(String(seq))
    n = length(s)
    out = Vector{UInt8}(undef, n)
    @inbounds for i in 1:n
        out[n-i+1] = COMPLEMENT_TABLE[s[i]+1]
    end
    return String(out)
end
