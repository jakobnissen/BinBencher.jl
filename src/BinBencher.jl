# Less whitespace in help string - https://github.com/comonicon/Comonicon.jl/issues/256
# In help string, improve the default value printing
# Help for the main function
# Filter for clades?

module BinBencher

using BinBencherBackend
using Comonicon: Comonicon

# TODO: This Base piracy is shady as fuck
Base.tryparse(::Type{Union{Nothing, String}}, x::String) = x
function Base.tryparse(::Type{Vector{Float64}}, x::String)
    [parse(Float64, i) for i in split(x, ',')]
end
function Base.tryparse(::Type{FlagSet}, x::String)
    FlagSet([tryparse(Flag, i) for i in split(x, ',')])
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

Comonicon.@cast function makeref(ref::String)
    reference = Reference(ref)
    println("It works!")
end

Comonicon.@main

end # module BinBencher
