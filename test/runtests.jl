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

BinBencher.makeref(
    refdir;
    seq_mapping = joinpath(BB_FILE_DIR, "mapping.tsv"),
    seq_fasta = joinpath(BB_FILE_DIR, "seq.fasta"),
    tax = joinpath(BB_FILE_DIR, "tax.tsv"),
    tax_ncbi = joinpath(BB_FILE_DIR, "ncbi.tsv"),
    genome_directories = "organism=$(joinpath(BB_FILE_DIR, "organisms")),virus=$(joinpath(BB_FILE_DIR, "virus"))",
    quiet = true,
)

refpath = joinpath(refdir, "reference.json")
ref = Reference(refpath)
ref_old = Reference(joinpath(BBB_FILE_DIR, "ref.json"))

@test test_is_same(ref, ref_old) === nothing

outdir = tempname()

BinBencher.bench(outdir, refpath, joinpath(BBB_FILE_DIR, "clusters.tsv"); quiet = true)

for expected_file in ["log.txt", "stats.json", "recovery.json", "bins.json.gz"]
    # Note: This is zero on nonexisting files
    @test filesize(joinpath(outdir, expected_file)) > 0
end
