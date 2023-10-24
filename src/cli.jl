# Less whitespace in help string - https://github.com/comonicon/Comonicon.jl/issues/256
# In help string, improve the default value printing
# Help for the main function
# Filter for clades?

# TODO: This Base piracy is shady as fuck
Base.tryparse(::Type{Union{Nothing, String}}, x::String) = x
function Base.tryparse(::Type{Vector{Float64}}, x::String)
    [parse(Float64, i) for i in split(x, ',')]
end
function Base.tryparse(::Type{FlagSet}, x::String)
    FlagSet([tryparse(Flag, i) for i in split(x, ',')])
end

exitwith(s::AbstractString) = (println(stderr, s); exit(1))

function checkfile(s::AbstractString, name::AbstractString)
    isfile(s) || exitwith("$(name) does not point to an existing file: \"$(s)\"")
end

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
    bins::String;
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
    binning = Binning(
        bins,
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
    print_matrix(binning; level=rank, assembly=assembly)
end

struct SeqArgs
    paths::Union{
        @NamedTuple{json::String},
        @NamedTuple{blast::String, fasta::String}
    }

    function SeqArgs(;
        json::Union{String, Nothing},
        blast::Union{String, Nothing},
        fasta::Union{String, Nothing},
    )
        if !isnothing(json)
            (isnothing(blast) && isnothing(fasta)) || exitwith("If seq-json-path is set, neither seq-blast-path nor seq-fasta-path can be")
            json = expanduser(json)
            checkfile(json, "seq-json-path")
            new((;json))
        else
            isnothing(blast) || isnothing(fasta) && exitwith("If seq-json-path is not set, both seq-fasta-path and seq-blast-path must be")
            blast = expanduser(blast)
            fasta = expanduser(fasta)
            checkfile(blast, "seq-blast-path")
            checkfile(fasta, "seq-fasta-path")
            new((;blast, fasta))
        end
    end
end

struct TaxArgs
    paths::Union{
        Nothing,
        @NamedTuple{json::String},
        @NamedTuple{tax::String, ncbi::Union{String, Nothing}}
    }
end

# virus+plasmid=/path/to/these,organism=/path/to/those
function parse_genomes_dir(s::String)::Vector{Pair{FlagSet, String}}
    map(split(strip(s), ',')) do segment
        p = findfirst(==(UInt8('=')), codeunits(segment))
        isnothing(p) && exitwith("No =") # TODO
        left = segment[1:prevind(segment, p)]
        right = segment[p+1:end]
        flags = map(split(left, '+')) do flag_string
            @something tryparse(Flag, strip(flag_string)) exitwith("Bad flag $(flag_string)") #TODO
        end |> FlagSet
        flags => expanduser(String(strip(right)))
    end
end

struct GenomeArgs
    x::Union{
        @NamedTuple{json::String},
        Vector{Pair{FlagSet, String}}
    }

    function GenomeArgs(json::Union{String, Nothing}, genome_directories::Union{String, Nothing})
        if !isnothing(json)
            isnothing(genome_directories) || exitwith("If genome-json-path is set, genome-directories must not be set")
            json = expanduser(json)
            checkfile(json, "genome-json-path")
            new((;json))
        else
            isnothing(genome_directories) && exitwith("If genome-json-path is not set, genome-directories must be set")
            v = parse_genomes_dir(genome_directories)
            seen_flagsets = Set{FlagSet}()
            for (flagset, path) in v
                flagset ∈ seen_flagsets && exitwith("Foo!") # TODO
                push!(seen_flagsets, flagset)
                isdir(path) || exitwith("Genome directory does not point to an existing directory: \"$(path)\"")
            end
            new(v)
        end
    end
end

# BLAST + (contigs file)
# Taxmaps file, of what precise format?
# Genomes: 
#   - [(flagset, dir) ... ]
#   - Or a way to JSON cast FlagSet to int and have it passed as JSON?
Comonicon.@cast function makeref(;
    # Sequences: BLAST path contains mappings from seqs to refs
    seq_blast_path::Union{String, Nothing}=nothing,
    # This is for unmapped sequences, so they still appear
    seq_fasta_path::Union{String, Nothing}=nothing,
    seq_json_path::Union{String, Nothing}=nothing,

    # Nothing (flat tax from genomes), or
    # The special tax path, possibly with the NCBI
    tax_path::Union{String, Nothing}=nothing,
    tax_ncbi_path::Union{String, Nothing}=nothing,
    tax_json_path::Union{String, Nothing}=nothing,

    genome_directories::Union{String, Nothing}=nothing,
    genome_json_path::Union{String, Nothing}=nothing,
)
    #=
    seq_args = SeqArgs(json=seq_json_path, fasta=seq_fasta_path, blast=seq_blast_path)
    seqs = sequences(seq_args)
    println(seqs)
    =#

    genome_args = GenomeArgs(genome_json_path, genome_directories)
    gns = genomes(genome_args)
    println(gns)
end


