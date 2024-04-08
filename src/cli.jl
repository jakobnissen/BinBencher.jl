# TODO
# Less whitespace in help string - https://github.com/comonicon/Comonicon.jl/issues/256
# In help string, improve the default value printing
# Help for the main function
# Filter for clades?

const DEFAULT_RECALLS = [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
const DEFAULT_PRECISIONS = [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

# TODO: This Base piracy is shady as fuck
Base.tryparse(::Type{Union{Nothing, String}}, x::String) = x
function Base.tryparse(::Type{Vector{Float64}}, x::String)
    map(eachsplit(x, ',')) do i
        @something tryparse(Float64, i) exitwith("Cannot parse as float: \"$(i)\"")
    end
end

function Base.tryparse(::Type{FlagSet}, x::String)
    Iterators.map(eachsplit(x, ',')) do str
        @something tryparse(Flag, str) parse_flag_error(str)
    end |> FlagSet
end

function parse_flag_error(s::AbstractString)
    exitwith(
        "Cannot parse as flag: \"$(s)\". " *
        "Available flags are: $(join(instances(Flag), ", "))",
    )
end

exitwith(s::AbstractString) = (println(stderr, s); exit(1))

function check_file(s::AbstractString, name::AbstractString)
    isfile(s) ? s : exitwith("$(name) does not point to an existing file: \"$(s)\"")
end

parentdir(dir::AbstractString) = dirname(rstrip(dir, ['\\', '/']))

function mkdir_checked(path::String, exists_ok::Bool)
    if exists_ok
        isfile(path) && exitwith("Output directory \"$(path)\" is an existing file")
    else
        ispath(path) && exitwith("Output directory \"$(path)\" already exists")
    end
    parent = parentdir(path)
    !isempty(parent) &&
        !isdir(parent) &&
        exitwith("Outdir's parent \"$(parent)\" is not an existing directory")
    mkdir(path)
end

"""
# Intro
Benchmark a set of bins agains a reference

# Args
- `ref`: Path to reference JSON file (see the `makeref subcommand`)
- `outdir`: Path of output directory to create (must not exist)
- `bins`: Path to one or more TSV files of bins.

# Options
- `-s, --sep`: Binsplit separator; delimiter to split bins
- `--minsize`: Remove bins with shorter breadth
- `--minseqs`: Remove bins fewer sequences
- `--keep-flags`: Only keep genomes with all these flags set.
  Takes a comma-separated list of flags, e.g. "virus,plasmid"
- `--remove-flags`: Remove genomes with any of these flags set. See above.
- `--recalls`: Comma-separated list of recall thresholds to use
- `--precisions`: Comma-separated list of precision thresholds to use

# Flags
- `--intersect`: Allow the same sequence to be in multiple bins
"""
Comonicon.@cast function bench(
    ref::String,
    outdir::String,
    bins...; # name=path strings, or a single path 
    sep::Union{String, Nothing}=nothing,
    minsize::Int=1,
    minseqs::Int=1,
    keep_flags::FlagSet=FlagSet(),
    remove_flags::FlagSet=FlagSet(),
    recalls::Union{Nothing, Vector{Float64}}=nothing,
    precisions::Union{Nothing, Vector{Float64}}=nothing,
    intersect::Bool=false,
)
    # Check validity of inputs
    pairs::Union{String, Vector{Pair{String, String}}} = if length(bins) == 1
        check_file(only(bins), "Binning file")
    else
        map(collect(bins)) do i
            p = split(i, '='; limit=2)
            length(p) == 1 && exitwith(
                "If multiple bins files is passed, pass them as e.g. bin1=path_to_bin1.tsv bin2=path_to_bin2.tsv",
            )
            # TODO: Handle spaces - perhaps strip quotes?
            String(p[1]) => String(p[2])
            #String(p[1]) => String(check_file(p[2], p[1]))
        end
    end
    println(pairs)
    exit(1)
    check_file(ref, "Reference")
    mkdir_checked(outdir, false)
    settings = OutputSettings(recalls, precisions)
    reference = Reference(ref)
    bin_paths::Vector{String} = pairs isa String ? [pairs] : last.(pairs)
    binnings = Vector{Binning}(undef, length(bin_paths))
    Threads.@threads for (i, path) in collect(enumerate(bin_paths))
        binnings[i] = Binning(
            path,
            reference;
            binsplit_separator=sep,
            min_size=minsize,
            min_seqs=minseqs,
            disjoint=!intersect,
            recalls=settings.recalls,
            precisions=settings.precisions,
            filter_genomes=g ->
                isdisjoint(flags(g), remove_flags) && issubset(keep_flags, flags(g)),
        )
    end
    if pairs isa String
        populate_output(outdir, only(binnings), settings)
    else
        for ((name, _), binning) in zip(pairs, binnings)
            subdir = joinpath(outdir, name)
            mkdir_checked(subdir, false)
            populate_output(subdir, binning, settings)
        end
    end
end

struct SeqArgs
    paths::Union{
        @NamedTuple{json::String},
        @NamedTuple{default_blast::String, fasta::String},
        @NamedTuple{custom_blast::String, fasta::String},
    }

    function SeqArgs(;
        json::Union{String, Nothing},
        default_blast::Union{String, Nothing},
        custom_blast::Union{String, Nothing},
        fasta::Union{String, Nothing},
    )
        if !isnothing(json)
            (isnothing(default_blast) && isnothing(custom_blast) && isnothing(fasta)) ||
                exitwith("If seq-json is set, neither seq-blast nor seq-fasta can be")
            json = expanduser(json)
            check_file(json, "seq-json")
            new((; json))
        elseif isnothing(fasta)
            exitwith("If seq-json is not set seq-fasta must be")
        else
            if !(isnothing(default_blast) ⊻ isnothing(custom_blast))
                exitwith(
                    "If seq-json is not set, exactly one of seq-custom-blast and seq-default-blast must be",
                )
            end
            fasta = expanduser(fasta)
            check_file(fasta, "seq-fasta")
            if !isnothing(default_blast)
                default_blast = expanduser(default_blast)
                check_file(default_blast, "seq-default-blast")
                new((; default_blast, fasta))
            else
                default_blast = expanduser(custom_blast::String)
                check_file(custom_blast::String, "seq-custom-blast")
                new((; custom_blast, fasta))
            end
        end
        error("unreachable!")
    end
end

# virus+plasmid=/path/to/these,organism=/path/to/those
function parse_genomes_dir(s::String)::Vector{Pair{FlagSet, String}}
    map(split(strip(s), ',')) do segment
        p = findfirst(==(UInt8('=')), codeunits(segment))
        isnothing(p) &&
            exitwith("No \"=\" symbol found in --genome-directories \"$(segment)]\"")
        left = segment[1:prevind(segment, p)]
        right = segment[(p + 1):end]
        flags =
            map(split(left, '+')) do flag_string
                @something tryparse(Flag, strip(flag_string)) parse_flag_error(
                    flag_string,
                )
            end |> FlagSet
        flags => expanduser(String(strip(right)))
    end
end

struct GenomeArgs
    x::Union{@NamedTuple{json::String}, Vector{Pair{FlagSet, String}}}

    function GenomeArgs(; json::Union{String, Nothing}, directories::Union{String, Nothing})
        if !isnothing(json)
            isnothing(directories) ||
                exitwith("If genome-json is set, genome-directories must not be set")
            json = expanduser(json)
            check_file(json, "genome-json")
            new((; json))
        else
            isnothing(directories) &&
                exitwith("If genome-json is not set, genome-directories must be set")
            v = parse_genomes_dir(directories)
            seen_flagsets = Set{FlagSet}()
            for (flagset, path) in v
                flagset ∈ seen_flagsets &&
                    exitwith("Flagset passed twice: $(join(flagset, '+'))")
                push!(seen_flagsets, flagset)
                isdir(path) || exitwith(
                    "Genome directory does not point to an existing directory: \"$(path)\"",
                )
            end
            new(v)
        end
    end
end

struct TaxArgs
    paths::Union{Nothing, @NamedTuple{json::String}, @NamedTuple{tax::String, ncbi::String}}

    function TaxArgs(;
        json::Union{String, Nothing},
        tax::Union{String, Nothing},
        ncbi::Union{String, Nothing},
    )
        if !isnothing(json)
            (!isnothing(ncbi) || !isnothing(tax)) &&
                exitwith("If tax-json is set, tax-ncbi and tax must not be set")
            json = expanduser(json)
            check_file(json, "tax-json")
            new((; json))
        elseif !(isnothing(tax) && isnothing(ncbi))
            (isnothing(tax) || isnothing(ncbi)) &&
                exitwith("If tax or tax-ncbi is set, both must be set")
            tax = expanduser(tax)
            check_file(tax, "tax")
            ncbi = expanduser(ncbi)
            check_file(ncbi, "tax-ncbi")
            new((; tax, ncbi))
        else
            isnothing(ncbi) || exitwith("If tax is not set, neither must tax-ncbi be")
            new(nothing)
        end
    end
end

"""
# Intro
Create a new reference JSON file. See more details in the documentation.

# Args
- `outdir`: Path to output directory

# Options
- `--seq-custom-blast`: Path to BLAST file of sequences to genomes (four-column format)
- `--seq-default-blast`: Path to BLAST file of sequences to genomes (default tabular)
- `--seq-fasta`: Path to FASTA file of sequences
- `--seq-json`: JSON file with sequences and their mappings.
  Incompatible with `seq-blast` and `seq-fasta`
- `--tax`: Path to TSV file with taxonomy for the genomes
- `--tax-ncbi`: Optional path to precomputed NCBI dump file
- `--tax-json`: JSON file with full taxonomy.
  Incompatible with `tax` and `tax-ncbi`
- `--genome_directories`: Comma-separated string of `genome+types=path`
  with the path to a directory of the genomes. Example:
  `organism=path/to/orgs,virus+plasmid=path/to/phasmids`
- `--genome-json`: JSON file with all genome information.
  Incompatible with `--genome-directories`

# Flags
- `--overwrite`: Do not error if output directory already exists
"""
Comonicon.@cast function makeref(
    outdir::String;
    # BLAST of sequences to genomes (four-column format, see docs)
    seq_custom_blast::Union{String, Nothing}=nothing,
    # BLAST of sequences to genomes (default tabular BLAST output)
    seq_default_blast::Union{String, Nothing}=nothing,
    # This is for unmapped sequences, so they still appear in reference
    seq_fasta::Union{String, Nothing}=nothing,
    # ... or they can pass in the JSON directly
    seq_json::Union{String, Nothing}=nothing,

    # The special tax file
    tax::Union{String, Nothing}=nothing,
    # Path to dump of NCBI - used to parse the `id=412` elements in `tax`
    tax_ncbi::Union{String, Nothing}=nothing,
    # ... or they can pass in the JSON directly
    tax_json::Union{String, Nothing}=nothing,

    # A string of the format
    # virus+plasmid=/path/to/phasmids,organisms=/path/to/organisms
    # I.e. (flagset, path) pairs comma-sep, with = to delimit the two, and +
    # to delimit flags in the flagset.
    genome_directories::Union{String, Nothing}=nothing,
    # ... or they can pass in the JSON directly
    genome_json::Union{String, Nothing}=nothing,

    # Do not error if output directory already exists
    overwrite::Bool=false,
)
    outdir = abspath(expanduser(outdir))
    mkdir_checked(outdir, overwrite)

    seq_args = SeqArgs(;
        json=seq_json,
        fasta=seq_fasta,
        default_blast=seq_default_blast,
        custom_blast=seq_custom_blast,
    )
    genome_args = GenomeArgs(; json=genome_json, directories=genome_directories)
    tax_args = TaxArgs(; json=tax_json, ncbi=tax_ncbi, tax)

    gns = genomes(genome_args)
    isnothing(genome_json) &&
        open(io -> JSON3.write(io, gns), joinpath(outdir, "genomes.json"), "w")

    taxes = taxonomy(tax_args, gns)
    isnothing(tax) ||
        open(io -> JSON3.write(io, taxes), joinpath(outdir, "taxonomy.json"), "w")

    seqs = sequences(seq_args, gns)
    isnothing(seq_json) &&
        open(io -> JSON3.write(io, seqs), joinpath(outdir, "seqs.json"), "w")
end
