# Regressions for the review findings PD1-PD6.

fastq(names...) = join(["@$n\nACGTACGT\n+\nIIIIIIII\n" for n in names])

@testset "Review regressions PD1-PD6" begin
    a = front_adapter("ACGT")

    # PD1: mate IDs are checked as by dnaio's record_names_match (results of
    # dnaio 1.2.4, used by cutadapt 5.2)
    match = PrimerDemux.record_names_match
    for (h1, h2, expected) in (("A/1", "A/2", true), ("A", "B", false), ("X.1 a", "X.3 b", true),
                               ("X1", "X2", true), ("Xa", "Xb", false), ("X/1", "X/2 c", true),
                               ("X1", "Y2", false), ("A1 x", "A2", true), ("AB", "AB\tz", true),
                               ("X3", "X1", true), ("X4", "X5", false), ("X4", "X4", true),
                               ("X1", "X", false), ("X", "X1", false), ("A B", "A C", true),
                               ("AB", "A B", false), ("A/1 1:N", "A/2 2:N", true), ("X0", "X1", false))
        @test match(Vector{UInt8}(h1), Vector{UInt8}(h2)) == expected
    end
    mktempdir() do dir
        r1, r2 = joinpath(dir, "r1.fq"), joinpath(dir, "r2.fq")
        sub = mkpath(joinpath(dir, "out"))
        outs = (joinpath(sub, "o1.fq"), joinpath(sub, "o2.fq"))
        # Illumina headers, the pair number after the space
        write(r1, fastq("M:1:A:1:1:10:20 1:N:0:1", "M:1:A:1:1:10:30 1:N:0:1"))
        write(r2, fastq("M:1:A:1:1:10:20 2:N:0:1", "M:1:A:1:1:10:30 2:N:0:1"))
        @test remove_primers((r1, r2), outs; adapters = a).reads_out == 2
        rm.(outs)
        for names2 in (("B/2",), ("M:1:A:1:1:10:30 2:N:0:1", "M:1:A:1:1:10:20 2:N:0:1"))
            names1 = length(names2) == 1 ? ("A/1",) : ("M:1:A:1:1:10:20 1:N:0:1", "M:1:A:1:1:10:30 1:N:0:1")
            write(r1, fastq(names1...))
            write(r2, fastq(names2...))
            for threads in (1, 4)
                err = try
                    remove_primers((r1, r2), outs; adapters = a, threads)
                catch e
                    e
                end
                @test err isa ArgumentError && occursin("improperly paired", err.msg)
                @test isempty(readdir(sub))
            end
        end
    end

    # An error of the reading task is reported also when it comes before the
    # first batch or after the last one (on Julia 1.10, iterating the closed
    # channel used to end the loop silently and publish the outputs).
    mktempdir() do dir
        r1, r2 = joinpath(dir, "r1.fq"), joinpath(dir, "r2.fq")
        write(r1, fastq("A/1"))
        write(r2, fastq("B/2"))
        bad = joinpath(dir, "bad.fq")
        write(bad, fastq(string.("r", 1:50)...) * "@x\nACGT\n+\nII\n")
        sub = mkpath(joinpath(dir, "out"))
        outs = (joinpath(sub, "o1.fq"), joinpath(sub, "o2.fq"))
        failures = 0
        for _ in 1:20, threads in (2, 3, 4)
            failures += try
                remove_primers((r1, r2), outs; adapters = a, threads)
                0
            catch e
                e isa ArgumentError ? 1 : 0
            end
            failures += try
                remove_primers(bad, outs[1]; adapters = a, threads, batch_size = 7)
                0
            catch
                1
            end
            isempty(readdir(sub)) || break
        end
        @test failures == 120
        @test isempty(readdir(sub))
    end

    # PD2: outputs that are the same file under another name are refused
    mktempdir() do dir
        input = joinpath(DATA, "reads.1.fastq")
        o = joinpath(dir, "o.fastq")
        sub = mkpath(joinpath(dir, "s"))
        aliases = [joinpath(sub, "..", "o.fastq")]
        (Sys.iswindows() || Sys.isapple()) && push!(aliases, joinpath(dir, "O.FASTQ"))
        for alias in aliases
            @test_throws ArgumentError remove_primers((input, input), (o, alias); adapters = a)
            @test readdir(dir) == ["s"]
        end
        # a hard link to the input, or two links to one existing file
        copy_in = joinpath(dir, "in.fastq")
        cp(input, copy_in)
        hardlink(copy_in, joinpath(dir, "link.fastq"))
        @test_throws ArgumentError remove_primers(copy_in, joinpath(dir, "link.fastq"); adapters = a, overwrite = true)
        write(o, "old")
        hardlink(o, joinpath(dir, "o2.fastq"))
        @test_throws ArgumentError remove_primers((input, input), (o, joinpath(dir, "o2.fastq"));
                                                  adapters = a, overwrite = true)
        @test read(o, String) == "old"
        # a directory is not an output
        @test_throws ArgumentError remove_primers(input, sub; adapters = a, overwrite = true)
        @test sort(readdir(dir)) == ["in.fastq", "link.fastq", "o.fastq", "o2.fastq", "s"]
    end

    # PD3: if a rename fails, new files are removed and replaced ones restored
    mktempdir() do dir
        paths = [joinpath(dir, "a.fastq"), joinpath(dir, "b.fastq"), joinpath(dir, "c.fastq")]
        write(paths[1], "old a")
        write(paths[2], "old b")
        outs = [PrimerDemux.OutputFile(p, true, 1) for p in paths]
        foreach(out -> write(out.io, "new"), outs)
        temp3 = outs[3].temp
        outs[3].temp = joinpath(dir, "missing")         # its rename fails
        @test_throws Base.IOError PrimerDemux.publish_all(outs)
        @test read(paths[1], String) == "old a"
        @test read(paths[2], String) == "old b"
        @test !ispath(paths[3])
        @test sort(readdir(dir)) == sort(["a.fastq", "b.fastq", basename(temp3)])
        rm(temp3)
        # the same outputs put in place
        outs = [PrimerDemux.OutputFile(p, true, 1) for p in paths]
        foreach(out -> write(out.io, "new"), outs)
        PrimerDemux.publish_all(outs)
        @test all(p -> read(p, String) == "new", paths)
        @test sort(readdir(dir)) == ["a.fastq", "b.fastq", "c.fastq"]
    end

    # PD4: a failure opening a later output discards the earlier ones
    mktempdir() do dir
        input = joinpath(DATA, "reads.1.fastq")
        write(joinpath(dir, "file"), "")
        @test_throws Exception remove_primers((input, input), (joinpath(dir, "o1.fastq"),
                                              joinpath(dir, "file", "o2.fastq")); adapters = a)
        @test readdir(dir) == ["file"]
    end

    # PD5: trim_sequence checks action and times as the file functions do
    @test_throws ArgumentError trim_sequence(a, "ACGTCCCC"; times = 0)
    @test_throws ArgumentError trim_sequence(front_adapter("AAA"), "AAACCCAAAGGG"; action = :retain, times = 2)
    @test_throws ArgumentError trim_sequence(a, "ACGTCCCC"; action = :crop)
    @test trim_sequence(a, "ACGTCCCC")[1] == "CCCC"

    # PD6: base as documented (keyword) and as fourth positional argument
    @test quality_trim_index(UInt8[73, 73], 0, 20; base = 33) == (0, 2)
    @test quality_trim_index(UInt8[73, 73], 0, 20, 33) == (0, 2)
    @test quality_trim_index(UInt8[73, 73], 0, 20; base = 64) == (0, 0)
end
