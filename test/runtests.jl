using BinBencher
using BinBencherBackend
using Test

BBB_FILE_DIR::String = let
    bb = pathof(BinBencherBackend)
    joinpath(dirname(dirname(bb)), "files")
end

BB_FILE_DIR::String = let
    bb = pathof(BinBencher)
    joinpath(dirname(dirname(bb)), "files")
end

include(joinpath(dirname(BBB_FILE_DIR), "test", "sameref.jl"))

refdir = tempname()

@test BinBencher.main([
    "makeref",
    refdir,
    "--seq-mapping", joinpath(BB_FILE_DIR, "mapping.tsv"),
    "--seq-fasta", joinpath(BB_FILE_DIR, "seq.fasta"),
    "--tax", joinpath(BB_FILE_DIR, "tax.tsv"),
    "--tax-ncbi", joinpath(BB_FILE_DIR, "ncbi.tsv"),
    "--genome-directories",
    "organism=$(joinpath(BB_FILE_DIR, "organisms")),virus=$(joinpath(BB_FILE_DIR, "virus"))",
    "--quiet",
]) == 0

refpath = joinpath(refdir, "reference.json")
ref = Reference(refpath)
ref_old = Reference(joinpath(BBB_FILE_DIR, "ref.json"))

@test test_is_same(ref, ref_old) === nothing

outdir = tempname()

@test BinBencher.main(["bench", outdir, refpath, joinpath(BBB_FILE_DIR, "clusters.tsv"), "--quiet"]) == 0

for expected_file in ["log.txt", "stats.json", "recovery.json", "bins.json.gz"]
    # Note: This is zero on nonexisting files
    @test filesize(joinpath(outdir, expected_file)) > 0
end

# Multiple named binnings, with options
outdir = tempname()
clusters = joinpath(BBB_FILE_DIR, "clusters.tsv")
@test BinBencher.main([
    "bench", outdir, refpath, "a=$(clusters)", "b=$(clusters)",
    "--recalls", "0.5,0.9", "--precisions", "0.7", "--remove-flags", "virus", "--quiet",
]) == 0
for name in ["a", "b"]
    @test filesize(joinpath(outdir, name, "recovery.json")) > 0
end

@testset "CLI help and errors" begin
    redirect_stdio(; stdout = devnull, stderr = devnull) do
        @test BinBencher.main(String[]) == 2
        @test BinBencher.main(["--help"]) == 0
        @test BinBencher.main(["bench", "--help"]) == 0
        @test BinBencher.main(["help", "makeref"]) == 0
        @test BinBencher.main(["--version"]) == 0
        @test BinBencher.main(["notacommand"]) == 2
        @test BinBencher.main(["bench", tempname(), refpath]) == 2 # no bins
        @test BinBencher.main(["bench", tempname(), refpath, clusters, "--minsize", "0"]) == 2
    end
end

@testset "Invalid input" begin
    dir = mktempdir()
    bad_tsv = joinpath(dir, "bad.tsv")
    write(bad_tsv, "not a header\n")
    makeref_args(outdir, mapping, tax) = [
        "makeref", outdir,
        "--seq-mapping", mapping,
        "--seq-fasta", joinpath(BB_FILE_DIR, "seq.fasta"),
        "--tax", tax,
        "--genome-directories",
        "organism=$(joinpath(BB_FILE_DIR, "organisms")),virus=$(joinpath(BB_FILE_DIR, "virus"))",
        "--quiet",
    ]
    redirect_stdio(; stdout = devnull, stderr = devnull) do
        # Invalid flag
        @test BinBencher.main(["bench", tempname(), refpath, clusters, "--keep-flags", "bogus"]) == 1
        # Existing output directory
        @test BinBencher.main(["bench", outdir, refpath, clusters]) == 1
        # Missing bins file
        @test BinBencher.main(["bench", tempname(), refpath, joinpath(dir, "missing.tsv")]) == 1
        # Invalid reference file (BinBencherBackend.ReferenceError)
        @test BinBencher.main(["bench", tempname(), clusters, clusters]) == 1
        # Invalid binning file (BinBencherBackend.BinningError)
        @test BinBencher.main(["bench", tempname(), refpath, bad_tsv]) == 1
        # Same binning name twice, where the error is thrown before the threaded loop
        @test BinBencher.main(["bench", tempname(), refpath, "a=$(clusters)", "a=$(clusters)"]) == 1
        # Malformed sequence mapping file, which is parsed in a spawned task
        @test BinBencher.main(makeref_args(tempname(), bad_tsv, joinpath(BB_FILE_DIR, "tax.tsv"))) == 1
        # Malformed taxonomy file
        @test BinBencher.main(makeref_args(tempname(), joinpath(BB_FILE_DIR, "mapping.tsv"), bad_tsv)) == 1
    end
    # When called as a function, invalid input throws a `CLIError`
    @test_throws BinBencher.CLIError BinBencher.bench(tempname(), refpath, [joinpath(dir, "missing.tsv")]; quiet = true)
end
