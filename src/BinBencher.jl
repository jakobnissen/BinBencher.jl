module BinBencher

using BinBencherBackend:
    BinBencherBackend,
    FlagSet,
    Reference,
    Binning,
    Flags,
    Flag,
    flags,
    open_perhaps_gzipped,
    TAXMAPS_JSON_T,
    GENOMES_JSON_T,
    SEQUENCES_JSON_T,
    RANK_BY_NAME,
    RANKS

using Comonicon: Comonicon
using FASTX: FASTA
using JSON3: JSON3
using CodecZlib: GzipDecompressorStream

const BBB = BinBencherBackend

include("cli.jl")
include("parse.jl")

Comonicon.@main

end # module BinBencher
