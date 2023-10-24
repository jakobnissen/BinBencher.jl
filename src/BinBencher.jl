module BinBencher

using BinBencherBackend
using BinBencherBackend: open_perhaps_gzipped, TAXMAPS_JSON_T, GENOMES_JSON_T, SEQUENCES_JSON_T 
using Comonicon: Comonicon
using FASTX: FASTA
using JSON3: JSON3

include("cli.jl")
include("parse.jl")

Comonicon.@main

end # module BinBencher