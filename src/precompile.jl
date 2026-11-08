# Precompile workload: Run the CLI on the small example files in `files/`,
# so the code paths of the `binbench` executable are compiled ahead of time.

@setup_workload begin
    files = joinpath(dirname(@__DIR__), "files")
    genome_directories = "organism=$(joinpath(files, "organisms")),virus=$(joinpath(files, "virus"))"
    clusters = joinpath(files, "clusters.tsv")
    @compile_workload begin
        # `main` sets the global logger, so we restore it afterwards
        logger = global_logger()
        try
            mktempdir() do dir
                refdir = joinpath(dir, "ref")
                main([
                    "makeref", refdir,
                    "--seq-mapping", joinpath(files, "mapping.tsv"),
                    "--seq-fasta", joinpath(files, "seq.fasta"),
                    "--tax", joinpath(files, "tax.tsv"),
                    "--tax-ncbi", joinpath(files, "ncbi.tsv"),
                    "--genome-directories", genome_directories,
                    "--quiet",
                ])
                ref = joinpath(refdir, "reference.json")
                main(["bench", joinpath(dir, "single"), ref, clusters, "--quiet"])
                main([
                    "bench", joinpath(dir, "multiple"), ref, "a=$(clusters)", "b=$(clusters)",
                    "--keep-flags", "organism", "--recalls", "0.5,0.9", "--precisions", "0.9",
                    "--quiet",
                ])
            end
        finally
            global_logger(logger)
        end
        print_cli_help(devnull, String[])
        print_cli_help(devnull, ["bench"])
        print_cli_help(devnull, ["makeref"])
    end
end
