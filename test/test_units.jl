function read_tsv(path)
    rows = Vector{Vector{String}}()
    for line in eachline(path)
        startswith(line, "#") && continue
        push!(rows, split(line, '\t'; keepempty = true))
    end
    return rows[1], rows[2:end]
end

@testset "Aligner.locate against cutadapt" begin
    header, rows = read_tsv(joinpath(REFERENCE, "aligner.tsv"))
    @test header[end] == "result"
    ws = PrimerDemux.AlignWorkspace()
    mismatches = 0
    for r in rows
        reference, query = r[1], r[2]
        al = PrimerDemux.Aligner(reference, parse(Float64, r[7]); flags = parse(Int, r[3]),
                                 wildcard_ref = r[4] == "1", wildcard_query = r[5] == "1",
                                 indel_cost = parse(Int, r[6]), min_overlap = parse(Int, r[8]))
        a = PrimerDemux.locate(al, Vector{UInt8}(codeunits(query)), ws)
        got = PrimerDemux.found(a) ?
              join((a.ref_start, a.ref_stop, a.query_start, a.query_stop, a.score, a.errors), ",") : "none"
        got == r[9] || (mismatches += 1)
    end
    @test length(rows) >= 1000
    @test mismatches == 0
end

@testset "k-mer pre-filter against cutadapt" begin
    _, rows = read_tsv(joinpath(REFERENCE, "kmers.tsv"))
    mismatches = 0
    overflows = 0
    for r in rows
        pk = PrimerDemux.create_positions_and_kmers(String(r[1]), parse(Int, r[2]), parse(Float64, r[3]),
                                                    r[4] == "1", r[5] == "1", r[6] == "1")
        kf = PrimerDemux.KmerFinder(pk, r[7] == "1", r[8] == "1")
        seq = Vector{UInt8}(codeunits(r[9]))
        present = PrimerDemux.kmers_present(kf, seq)
        # A 5' window longer than the read makes cutadapt read past the string
        # (undefined result); PrimerDemux then always runs the aligner.
        if any(s -> s.stop > length(seq), kf.searches)
            overflows += 1
            present || (mismatches += 1)
        else
            present == (r[10] == "1") || (mismatches += 1)
        end
    end
    @test length(rows) - overflows >= 1000
    @test mismatches == 0
end

@testset "Adapter specifications" begin
    a = front_adapter("fwd=^ACGTACGTAC;e=0.2")
    @test a isa SingleAdapter
    @test a.name == "fwd"
    @test a.placement == PREFIX
    @test a.max_error_rate == 0.2
    @test a.min_overlap == 10                       # anchored: full length
    @test front_adapter("ACGT").placement == FRONT
    @test front_adapter("XACGT").placement == FRONT_NOT_INTERNAL
    @test back_adapter("ACGT").placement == BACK
    @test back_adapter("ACGT\$").placement == SUFFIX
    @test back_adapter("ACGTX").placement == BACK_NOT_INTERNAL
    @test back_adapter("ACGT...").placement == FRONT    # -a ADAPTER... is a 5' adapter
    @test back_adapter("...ACGT").placement == BACK
    @test front_adapter("ACGTACGTAC").name == "ACGTACGTAC"
    @test front_adapter("ACGT"; name = "x").name == "x"
    @test back_adapter("acgu").sequence == "ACGT"
    @test front_adapter("ACGTN"; max_errors = 2).max_error_rate == 0.5   # 2 errors / 4 non-N bases
    @test front_adapter("ACGT"; min_overlap = 10).min_overlap == 4
    @test back_adapter("A{3}C{2}G").sequence == "AAACCG"
    @test front_adapter("ACGTACGT;o=5;noindels").indels == false

    l = front_adapter("p=^AAAACCCC...GGGGTTTT\$")
    @test l isa LinkedAdapter
    @test l.name == "p"
    @test l.front.placement == PREFIX && l.back.placement == SUFFIX
    @test l.front_required && l.back_required
    l = back_adapter("AAAACCCC...GGGGTTTT")
    @test !l.front_required && !l.back_required          # -a: only anchored parts
    l = back_adapter("^AAAACCCC...GGGGTTTT")
    @test l.front_required && !l.back_required
    l = front_adapter("AAAACCCC;optional...GGGGTTTT")
    @test !l.front_required && l.back_required

    @test_throws ArgumentError front_adapter("")
    @test_throws ArgumentError front_adapter("ACGT\$")
    @test_throws ArgumentError back_adapter("^ACGT")
    @test_throws ArgumentError front_adapter("^XACGT")
    @test_throws ArgumentError front_adapter("ACGTZ")
    @test_throws ArgumentError front_adapter("ACGT;foo=1")
    @test_throws ArgumentError front_adapter("^ACGT;o=3")
    @test_throws ArgumentError front_adapter("ACGT;required")
    @test_throws ArgumentError front_adapter("ACGT;e=0.1;e=0.2")
    @test_throws ArgumentError front_adapter("{3}ACGT")
    @test_throws ArgumentError front_adapter("...ACGT")
    @test_throws ArgumentError front_adapter("NNNN")    # only wildcards

    fasta = adapters_from_fasta(joinpath(DATA, "barcodes.fasta"); anchored = true)
    @test length(fasta) == 5
    @test fasta[1].name == "bc1"
    @test fasta[1].placement == PREFIX
    @test all(a -> a.placement == SUFFIX, adapters_from_fasta(joinpath(DATA, "barcodes.fasta"); type = :back, anchored = true))
end

@testset "Single sequences" begin
    primer = front_adapter("ACGTTGCAAGGC")
    m = find_adapter(primer, "TTTACGTTGCAAGGCGATTACA")
    @test m isa AdapterMatch
    @test m.front.rstart == 3 && m.front.rstop == 15
    @test PrimerDemux.score(m) == 12 && PrimerDemux.errors(m) == 0
    @test find_adapter(primer, "GATTACAGATTACA") === nothing
    seq, matches, rc = trim_sequence(primer, "TTTACGTTGCAAGGCGATTACA")
    @test seq == "GATTACA" && length(matches) == 1 && !rc
    # one substitution (allowed: 0.1 * 12 = 1.2 errors)
    @test trim_sequence(primer, "ACGTAGCAAGGCGATTACA")[1] == "GATTACA"
    # two substitutions are too many
    @test trim_sequence(primer, "ACCTAGCAAGGCGATTACA")[1] == "ACCTAGCAAGGCGATTACA"
    # reverse complement
    seq, _, rc = trim_sequence(primer, reverse_complement("TTTACGTTGCAAGGCGATTACA"); revcomp = true)
    @test seq == "GATTACA" && rc
    # partial 3' adapter at the end of the read
    @test trim_sequence(back_adapter("AGATCGGAAGAGC"), "GATTACAGATTACAAGATC")[1] == "GATTACAGATTACA"
    # IUPAC primer
    @test trim_sequence(front_adapter("^AAGCTWMTGG"), "AAGCTAATGGCCCC")[1] == "CCCC"
    # linked primers, both required with front_adapter
    linked = front_adapter("^AAGGTTCCAA...GGCCAATTGG")
    @test trim_sequence(linked, "AAGGTTCCAAGATTACAGGCCAATTGGTT")[1] == "GATTACA"
    @test trim_sequence(linked, "AAGGTTCCAAGATTACA")[1] == "AAGGTTCCAAGATTACA"
    # best adapter: higher score wins over list order
    ads = [front_adapter("a=ACGTACGTAC"), front_adapter("b=TTTTGGGGCCCCAAAA")]
    @test find_adapter(ads, "TTTTGGGGCCCCAAAAGATTACA").adapter == 2
    # actions
    @test trim_sequence(back_adapter("GGCCAATTGG"), "GATTACAGGCCAATTGGTT"; action = :mask)[1] == "GATTACANNNNNNNNNNNN"
    @test trim_sequence(back_adapter("GGCCAATTGG"), "GATTACAGGCCAATTGGTT"; action = :lowercase)[1] == "GATTACAggccaattggtt"
    @test trim_sequence(back_adapter("GGCCAATTGG"), "GATTACAGGCCAATTGGTT"; action = :retain)[1] == "GATTACAGGCCAATTGG"
    @test trim_sequence(back_adapter("GGCCAATTGG"), "GATTACAGGCCAATTGGTT"; action = :none)[1] == "GATTACAGGCCAATTGGTT"
    @test_throws ArgumentError trim_sequence(primer, "ACGT"; action = :crop)
end

@testset "Quality trimming and helpers" begin
    q(s) = Vector{UInt8}(codeunits(s))
    # cutadapt documentation example: qualities 42 40 26 27 8 7 11 4 2 3, cutoff 10
    quals = UInt8.([42, 40, 26, 27, 8, 7, 11, 4, 2, 3] .+ 33)
    @test quality_trim_index(quals, 0, 10) == (0, 4)
    @test quality_trim_index(q("!!!!IIII"), 10, 0) == (4, 8)
    @test quality_trim_index(q("!!!!"), 10, 10) == (0, 0)
    @test quality_trim_index(UInt8[], 10, 10) == (0, 0)
    @test reverse_complement("ACGTRYKMNacgu-") == "-acgtNKMRYACGT"
    @test PrimerDemux.expand_braces("TGA{5}CT") == "TGAAAAACT"
    @test PrimerDemux.expected_errors(q("++")) ≈ 0.2       # Q10 twice
end
