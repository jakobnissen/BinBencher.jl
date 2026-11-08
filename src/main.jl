# Command-line interface: Parsing of arguments, and the `main` entry point

const PROGNAME = "binbench"

@parser struct BenchArgs
    @description "Benchmark one or more binnings against a reference."
    "Binsplit separator; delimiter to split bins"
    sep = optional(option("-s", "--sep", str("SEP")))
    "Remove bins with shorter breadth"
    minsize = default(option("--minsize", integer("N"; min = 1)), 1)
    "Remove bins with fewer sequences"
    minseqs = default(option("--minseqs", integer("N"; min = 1)), 1)
    "Only keep genomes with all these flags set. Comma-separated, e.g. \"virus,plasmid\""
    keep_flags = optional(option("--keep-flags", str("FLAGS")))
    "Remove genomes with any of these flags set. Comma-separated, e.g. \"virus,plasmid\""
    remove_flags = optional(option("--remove-flags", str("FLAGS")))
    # Help is attached with `help`, since `@parser` only accepts literal strings as help
    recalls = optional(option("--recalls", str("RECALLS"))) |>
        help("Comma-separated list of recall thresholds [default: $(join(DEFAULT_RECALLS, ','))]")
    # Help is attached with `help`, since `@parser` only accepts literal strings as help
    precisions = optional(option("--precisions", str("PRECISIONS"))) |>
        help("Comma-separated list of precision thresholds [default: $(join(DEFAULT_PRECISIONS, ','))]")
    "Allow the same sequence to be in multiple bins"
    intersect = flag("--intersect")
    "Disable non-error logging to stderr"
    quiet = flag("--quiet")
    "Path of output directory to create (must not exist)"
    outdir = arg(str("OUTDIR"))
    "Path to reference JSON file (see the `makeref` subcommand)"
    ref = arg(str("REF"))
    "Path to a TSV file of bins, or multiple as name=path, e.g. a=a.tsv b=b.tsv"
    bins = many1(arg(str("BINS")))
end

@parser struct MakerefArgs
    @description "Create a new reference JSON file. See more details in the documentation."
    "Mapping of sequences to genomes in 4-column TSV file"
    seq_mapping = optional(option("--seq-mapping", str("PATH")))
    "FASTA file of sequences"
    seq_fasta = optional(option("--seq-fasta", str("PATH")))
    "JSON file with sequences and their mappings. Incompatible with --seq-mapping and --seq-fasta"
    seq_json = optional(option("--seq-json", str("PATH")))
    "TSV file with taxonomy for the genomes"
    tax = optional(option("--tax", str("PATH")))
    "Precomputed NCBI dump file"
    tax_ncbi = optional(option("--tax-ncbi", str("PATH")))
    "JSON file with full taxonomy. Incompatible with --tax and --tax-ncbi"
    tax_json = optional(option("--tax-json", str("PATH")))
    "Comma-separated flags=path pairs of genome directories, e.g. organism=orgs,virus+plasmid=phasmids"
    genome_directories = optional(option("--genome-directories", str("DIRS")))
    "JSON file with all genome information. Incompatible with --genome-directories"
    genome_json = optional(option("--genome-json", str("PATH")))
    "Do not error if output directory already exists"
    overwrite = flag("--overwrite")
    "Disable non-error logging to stderr"
    quiet = flag("--quiet")
    "Path to output directory"
    outdir = arg(str("OUTDIR"))
end

const CLI = or(
    command("bench", lift(BenchArgs)) |> help("Benchmark binnings against a reference"),
    command("makeref", lift(MakerefArgs)) |> help("Create a new reference JSON file"),
)

function run_command(x::BenchArgs)
    return bench(
        x.outdir, x.ref, x.bins;
        sep = x.sep,
        minsize = x.minsize,
        minseqs = x.minseqs,
        keep_flags = isnothing(x.keep_flags) ? FlagSet() : parse_flagset(x.keep_flags),
        remove_flags = isnothing(x.remove_flags) ? FlagSet() : parse_flagset(x.remove_flags),
        recalls = isnothing(x.recalls) ? nothing : parse_thresholds(x.recalls),
        precisions = isnothing(x.precisions) ? nothing : parse_thresholds(x.precisions),
        intersect = x.intersect,
        quiet = x.quiet,
    )
end

function run_command(x::MakerefArgs)
    return makeref(
        x.outdir;
        seq_mapping = x.seq_mapping,
        seq_fasta = x.seq_fasta,
        seq_json = x.seq_json,
        tax = x.tax,
        tax_ncbi = x.tax_ncbi,
        tax_json = x.tax_json,
        genome_directories = x.genome_directories,
        genome_json = x.genome_json,
        overwrite = x.overwrite,
        quiet = x.quiet,
    )
end

# The help text is generated as a `String`, and then printed. That way, the precompile
# workload can compile the generation of help, which is the expensive part, without
# printing to `Core.stdout`.
function print_cli_help(io::IO, argv::Vector{String})
    print(io, generate_help(CLI, argv; progname = PROGNAME))
    return nothing
end

"""
    main(args::Vector{String})::Cint

Entry point of the `binbench` executable. Return the exit code: 0 on success,
1 on invalid input (a `CLIError`), and 2 on invalid command-line arguments.
"""
function (@main)(args::Vector{String})::Cint
    if isempty(args)
        print_cli_help(Core.stderr, String[])
        return 2
    end
    # Help and version are handled here, and not by `OptParse.runparse`, since
    # that prints help to stderr, and does not distinguish help from parse errors.
    if args[1] in ("-V", "--version")
        println(Core.stdout, PROGNAME, " ", VERSION)
        return 0
    elseif args[1] == "help"
        print_cli_help(Core.stdout, args[2:end])
        return 0
    end
    before_terminator = something(findfirst(==("--"), args), lastindex(args) + 1) - 1
    if any(in(("-h", "--help")), view(args, 1:before_terminator))
        subcommand = filter(!in(("-h", "--help")), args[1:before_terminator])
        print_cli_help(Core.stdout, subcommand)
        return 0
    end
    parsed = optparse(CLI, args; progname = PROGNAME)
    parsed === nothing && return 2
    # Set the logger before running the command, so that errors thrown before
    # the command sets it (e.g. when parsing flags) are logged in the same format.
    set_global_logger!(; quiet = parsed.quiet)
    return try
        run_command(parsed)
        0
    catch e
        e isa CLIError || rethrow()
        @error e.msg
        1
    end
end
