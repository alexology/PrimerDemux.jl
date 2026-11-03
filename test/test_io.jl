using CodecZlibNG: GzipCompressorStream

@testset "Files and options" begin
    input = joinpath(DATA, "reads.1.fastq")
    input2 = joinpath(DATA, "reads.2.fastq")
    linked = front_adapter("^$FWD...$REV_RC")
    expected = read_gz(joinpath(REFERENCE, "linked_revcomp", "out.fastq.gz"))
    run(path; kw...) = remove_primers(input, path; adapters = linked, revcomp = true,
                                      discard_untrimmed = true, kw...)

    mktempdir() do dir
        # thread counts and batch sizes do not change the output
        for (threads, batch_size) in ((1, 4096), (4, 7), (3, 1))
            out = joinpath(dir, "t$(threads)_$(batch_size).fastq")
            stats = run(out; threads, batch_size)
            @test read(out, String) == expected
            @test stats.reads_in == 150
            @test stats.reads_out == count(==('\n'), expected) ÷ 4
            @test stats.with_adapters == [stats.reads_out]
            @test stats.written == ["output" => stats.reads_out]
        end

        # gzip input and output; FASTA output from FASTQ input
        gz_in = joinpath(dir, "in.fastq.gz")
        open(GzipCompressorStream, gz_in, "w") do io
            write(io, read(input))
        end
        for threads in (1, 4)
            gz_out = joinpath(dir, "out$threads.fastq.gz")
            remove_primers(gz_in, gz_out; adapters = linked, revcomp = true, discard_untrimmed = true, threads)
            @test read_gz(gz_out) == expected
        end

        # a gzip output (one member per batch) read back as input
        multi = joinpath(dir, "multi.fastq.gz")
        remove_primers(input, multi; batch_size = 7)
        @test read_gz(multi) == read(input, String)
        for threads in (1, 4)
            again = joinpath(dir, "multi_again$threads.fastq")
            remove_primers(multi, again; adapters = linked, revcomp = true, discard_untrimmed = true,
                           threads, batch_size = 11)
            @test read(again, String) == expected
        end

        # a malformed record early in a large gzip input: the decompressing
        # task, paused on a full buffer, is stopped and the error reported
        big = joinpath(dir, "big.fastq.gz")
        open(GzipCompressorStream, big, "w") do io
            write(io, "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nII\n")
            record = "@r\n" * "ACGT"^50 * "\n+\n" * "I"^200 * "\n"
            for _ in 1:120_000
                write(io, record)
            end
        end
        sub = mkpath(joinpath(dir, "out_big"))
        @test_throws Exception remove_primers(big, joinpath(sub, "o.fastq"); adapters = linked, threads = 4)
        @test isempty(readdir(sub))

        # a truncated or corrupt gzip input is an error, also when decompressed
        # by a separate task, and leaves no output
        bytes = read(gz_in)
        for (name, data) in (("truncated", bytes[1:end÷2]),
                             ("corrupt", vcat(bytes[1:200], fill(0xff, 300), bytes[501:end])))
            path = joinpath(dir, "$name.fastq.gz")
            write(path, data)
            sub = mkpath(joinpath(dir, "out_$name"))
            for threads in (1, 4)
                @test_throws Exception remove_primers(path, joinpath(sub, "o.fastq"); adapters = linked, threads)
                @test isempty(readdir(sub))
            end
        end
        fa_out = joinpath(dir, "out.fasta")
        run(fa_out)
        lines = split(expected, '\n'; keepempty = false)
        fasta = join([">" * lines[i][2:end] * "\n" * lines[i+1] * "\n" for i in 1:4:length(lines)])
        @test read(fa_out, String) == fasta

        # FASTA input gives FASTA output
        fa_out2 = joinpath(dir, "again.fastq")
        remove_primers(fa_out, fa_out2; adapters = back_adapter("GATTACA"))
        @test startswith(read(fa_out2, String), ">")

        # existing outputs are kept unless overwrite = true
        @test_throws ArgumentError run(fa_out)
        @test read(fa_out, String) == fasta
        run(fa_out; overwrite = true)
        @test read(fa_out, String) == fasta

        # a failure leaves neither output nor temporary files
        bad = joinpath(dir, "bad.fastq")
        write(bad, "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nII\n")
        sub = mkpath(joinpath(dir, "failure"))
        @test_throws Exception remove_primers(bad, joinpath(sub, "o.fastq"); adapters = linked)
        @test isempty(readdir(sub))

        # empty input
        empty = joinpath(dir, "empty.fastq")
        touch(empty)
        stats = remove_primers(empty, joinpath(dir, "empty_out.fastq"); adapters = linked)
        @test stats.reads_in == 0
        @test read(joinpath(dir, "empty_out.fastq"), String) == ""

        # paired files must have the same number of reads
        short = joinpath(dir, "short.fastq")
        write(short, join(collect(eachline(input2))[1:8], "\n") * "\n")
        @test_throws ArgumentError remove_primers((input, short), (joinpath(dir, "p1.fq"), joinpath(dir, "p2.fq"));
                                                  adapters = linked)
    end
end

@testset "Argument checks" begin
    input = joinpath(DATA, "reads.1.fastq")
    a = front_adapter("ACGTACGT")
    mktempdir() do dir
        o = joinpath(dir, "o.fastq")
        @test_throws ArgumentError remove_primers(input, joinpath(dir, "{name}.fastq"); adapters = a)
        @test_throws ArgumentError demultiplex(input, o; adapters = a)
        @test_throws ArgumentError demultiplex(input, joinpath(dir, "{name}.fastq"); adapters = Adapter[])
        @test_throws ArgumentError remove_primers(input, o; adapters = a, discard_untrimmed = true, discard_trimmed = true)
        @test_throws ArgumentError remove_primers(input, o; adapters = a, action = :retain, times = 2)
        @test_throws ArgumentError remove_primers(input, o; adapters = a, action = :crop)
        @test_throws ArgumentError remove_primers(input, o; adapters = a, adapters2 = a)
        @test_throws ArgumentError remove_primers(input, o; adapters = a, min_length = (10, 20))
        @test_throws ArgumentError remove_primers(input, (o, o * "2"); adapters = a)
        @test_throws ArgumentError remove_primers((input, input), (o, o); adapters = a)
        @test_throws ArgumentError remove_primers((input, input), (o, o * "2"); adapters = a, pair_filter = :all)
        @test_throws ArgumentError remove_primers((input, input), (o, o * "2"); adapters = a, adapters2 = [a, a],
                                                  pair_adapters = true)
        @test_throws ArgumentError remove_primers(input, input; adapters = a, overwrite = true)
        # two adapters with the same name share one demultiplexed file
        stats = demultiplex(input, joinpath(dir, "{name}.fastq");
                            adapters = [front_adapter("x=ACGTACGT"), front_adapter("x=TTGCAAGC")])
        @test sort(readdir(dir)) == ["unknown.fastq", "x.fastq"]
        @test [p.first for p in stats.written] == ["x", "unknown"]
        @test sum(p.second for p in stats.written) == 150
        @test sprint(show, MIME("text/plain"), stats) isa String
    end
end
