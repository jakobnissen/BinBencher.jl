module BinBencher

using BinBencherBackend
using BinBencherBackend:
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

include("cli.jl")
include("parse.jl")

Comonicon.@main

end # module BinBencher
