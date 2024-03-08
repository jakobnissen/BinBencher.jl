# TODO
# Less whitespace in help string - https://github.com/comonicon/Comonicon.jl/issues/256
# In help string, improve the default value printing
# Help for the main function
# Filter for clades?

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

function checkfile(s::AbstractString, name::AbstractString)
    isfile(s) || exitwith("$(name) does not point to an existing file: \"$(s)\"")
end

parentdir(dir::AbstractString) = dirname(rstrip(dir, ['\\', '/']))

"""
# Intro
Benchmark a set of bins agains a reference

# Args
- `ref`: Path to reference JSON file (see the `makeref subcommand`)
- `bins`: Path to TSV file of bins.

# Options
- `-r, --rank`: Taxonomic rank to benchmark
- `-s, --sep`: Binsplit separator; delimiter to split bins
- `--minsize`: Remove bins with shorter breadth
- `--minseqs`: Remove bins fewer sequences
- `--keep-flags`: Only keep genomes with all these flags set.
  Takes a comma-separated list of flags, e.g. "virus,plasmid"
- `--remove-flags`: Remove genomes with any of these flags set. See above.
- `--recalls`: Comma-separated list of recall thresholds to use
- `--precisions`: Comma-separated list of precision thresholds to use

# Flags
- `--assembly`: Count number of recovered assemblies, not genomes
- `--intersect`: Allow the same sequence to be in multiple bins
"""
Comonicon.@cast function bench(
    ref::String,
    bins...; # TODO: Explicit type here, but Comonicon issue #261 prevents this
    rank::Int=0,
    sep::Union{String, Nothing}=nothing,
    minsize::Int=1,
    minseqs::Int=1,
    keep_flags::FlagSet=FlagSet(),
    remove_flags::FlagSet=FlagSet(),
    recalls::Vector{Float64}=[0.6, 0.7, 0.8, 0.9, 0.95, 0.99],
    precisions::Vector{Float64}=[0.6, 0.7, 0.8, 0.9, 0.95, 0.99],
    assembly::Bool=false,
    intersect::Bool=false,
)
    reference = Reference(ref)
    binnings = Vector{Binning}(undef, length(bins))
    Threads.@threads for (i, path) in collect(enumerate(bins))
        binnings[i] = Binning(
            path,
            reference;
            binsplit_separator=sep,
            min_size=minsize,
            min_seqs=minseqs,
            disjoint=!intersect,
            recalls,
            precisions,
            filter_genomes=g ->
                isdisjoint(flags(g), remove_flags) && issubset(keep_flags, flags(g)),
        )
    end
    for (path, binning) in zip(bins, binnings)
        println(path)
        BBB.print_matrix(binning; level=rank, assembly=assembly)
    end
end

struct SeqArgs
    paths::Union{@NamedTuple{json::String}, @NamedTuple{blast::String, fasta::String}}

    function SeqArgs(;
        json::Union{String, Nothing},
        blast::Union{String, Nothing},
        fasta::Union{String, Nothing},
    )
        if !isnothing(json)
            (isnothing(blast) && isnothing(fasta)) ||
                exitwith("If seq-json is set, neither seq-blast nor seq-fasta can be")
            json = expanduser(json)
            checkfile(json, "seq-json")
            new((; json))
        else
            (isnothing(blast) || isnothing(fasta)) &&
                exitwith("If seq-json is not set, both seq-fasta and seq-blast must be")
            blast = expanduser(blast)
            fasta = expanduser(fasta)
            checkfile(blast, "seq-blast")
            checkfile(fasta, "seq-fasta")
            new((; blast, fasta))
        end
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
            checkfile(json, "genome-json")
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
            checkfile(json, "tax-json")
            new((; json))
        elseif !(isnothing(tax) && isnothing(ncbi))
            (isnothing(tax) || isnothing(ncbi)) &&
                exitwith("If tax or tax-ncbi is set, both must be set")
            tax = expanduser(tax)
            checkfile(tax, "tax")
            ncbi = expanduser(ncbi)
            checkfile(ncbi, "tax-ncbi")
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
- `--seq-blast`: Path to BLAST file of sequences to genomes
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
    # BLAST of sequences to genomes
    seq_blast::Union{String, Nothing}=nothing,
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
    isfile(outdir) && exitwith("Output directory \"$(outdir)\" is an existing file")
    isdir(parentdir(outdir)) ||
        exitwith("Outdir's parent \"$(parentdir(outdir))\" is not an existing directory")
    if isdir(outdir)
        overwrite || exitwith("Output directory \"$(outdir)\" is an existing directory")
    else
        mkdir(outdir)
    end

    seq_args = SeqArgs(; json=seq_json, fasta=seq_fasta, blast=seq_blast)
    genome_args = GenomeArgs(; json=genome_json, directories=genome_directories)
    tax_args = TaxArgs(; json=tax_json, ncbi=tax_ncbi, tax)

    gns = genomes(genome_args)
    isnothing(genome_json) &&
        open(io -> JSON3.write(io, gns), joinpath(outdir, "genomes.json"), "w")
    println(gns)

    taxes = taxonomy(tax_args, gns)
    isnothing(tax) ||
        open(io -> JSON3.write(io, taxes), joinpath(outdir, "taxonomy.json"), "w")
    println(taxes)

    seqs = sequences(seq_args, gns)
    isnothing(seq_json) &&
        open(io -> JSON3.write(io, seqs), joinpath(outdir, "seqs.json"), "w")
    println(seqs)
end
