# child_id => (rank, parent_id, name)
const NCBI_DICT_T = Dict{Int32, Tuple{Int32, Int32, String}}

# child_id  rank    parent_id   name
function parse_ncbi_taxes(io::IO)::NCBI_DICT_T
    result = NCBI_DICT_T()
    for (child_str, rank_str, parent_str, name) in Iterators.map(i -> split(i, '\t'), eachline(io))
        result[parse(Int32, child_str)] = (parse(Int32, rank_str), parse(Int32, parent_str), name)
    end
    result
end

# strain   child   parent
# class    child   id=3143
function parse_taxmaps(io::IO, ncbi::NCBI_DICT_T)::TAXMAPS_JSON_T
    result = Dict(i => Dict{String, String}() for i in 0:7)
    for (line_number, line) in enumerate(eachline(io))
        fields = split(rstrip(line), '\t')
        length(fields) == 3 || error(
            lazy"In taxmaps file, on line $(line_number), expected 3 tab-separated fields, got $(length(fields))"
        )
        (rank_str, child, parent) = fields
        rank = get(RANK_BY_NAME, lowercase(rank_str), 8)
        rank == 8 && error(
            lazy"Invalid rank name: \"$(rank_str)\". Valid names are:\n$(join(RANKS[1:end-1], \", \"))"
        )
        if !startswith(parent, "id=")
            existing = get!(result[rank], child, parent)
            existing == parent || error(lazy"At rank \"$(rank_str)\", child \"$(child)\" is has two distinct parents")
        else
            parent_id = parse(Int, parent[4:end])
            for rank in rank:6
                lookup = get(ncbi, parent_id, nothing)
                isnothing(lookup) && error(
                    lazy"NCBI taxid $(parent_id) is not a known NCBI taxid. Please note that we only accept taxids with a rank " *
                    lazy"of $(join(RANKS[1:end-1], \", \")). Verify the NCBI taxid, or manually include the taxonomic levels for child $(child)"
                ) # we have only the normal ranks
                (parent_rank, grandparent_id, parent) = lookup
                parent_rank == rank + 1 || error(
                    lazy"At rank $(rank_str), NCBI taxid $(parent_id) is listed as a parent. " *
                    lazy"However, the parent of a $(rank_str) must be a $(RANKS[rank+1]), whereas NCBI tax $(parent_id) is a $(RANKS[parent_rank])"
                )
                existing = get!(result[rank], child, parent)
                existing == parent || error(
                    lazy"At rank $(RANKS[rank]), child \"$(child)\" is listed with multiple distinct parents."
                )
                child = parent
                parent_id = grandparent_id
            end
        end
    end
    # Remove empty dicts in result
    vector = [[(k,v) for (k, v) in d] for d in [result[i] for i in 0:7]]
    while isempty(last(vector))
        pop!(vector)
        isempty(vector) && return vector
    end
    # Add top clade if  there isn't one
    last_parents = unique(map(last, last(vector)))
    if length(last_parents) > 1
        push!(vector, [(p, "top") for p in last_parents])
    end
    vector
end

function empty_taxmaps(x::GENOMES_JSON_T)::TAXMAPS_JSON_T
    [[(genome_name, "top") for (genome_name, _, _) in x]]
end



# qseqid     sseqid     qlen     sstart    send
function parse_sequences_blast(io::IO)::Dict{String, Tuple{Int, Vector{Tuple{String, Int, Int}}}}
    by_query = Dict{String, Tuple{Int, Vector{Tuple{String, Int, Int}}}}()
    for (lineno, line) in enumerate(eachline(io))
        stripped = rstrip(line)
        isempty(stripped) && continue
        m = match(r"^(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\d+)\t?", stripped)
        m === nothing && error(
            "On line $(lineno), expected line of the form QSEQID \t SSEQID \t QLEN \t SSTART \t SEND [\t FIELDS...] " *
            lazy"but got \"$(stripped)\""
        )
        (query, subject, qlen, sstart, send) = map(something, m.captures)
        len = parse(Int, qlen)
        from = parse(Int, sstart)
        to = parse(Int, send)
        (from, to) = minmax(from, to)
        (_, targets) = get!(by_query, query) do
            (len, Tuple{String, Int, Int}[])
        end
        push!(targets, (String(subject), from, to))
    end
    by_query
end

function parse_sequences_fasta(io::IO)
    map(FASTA.Reader(io; copy=false)) do record
        String(FASTA.identifier(record)), FASTA.seqsize(record)
    end
end

function parse_ref(
    taxmaps_path::Union{
        TAXMAPS_JSON_T,
        Union{Nothing, @NamedTuple{ncbi::String, taxmaps::String}}
    },
    blast_path::Union{SEQUENCES_JSON_T, String},
    genomes::Union{GENOMES_JSON_T, Vector{}},
)::Reference
    genomes_json = open_perhaps_gzipped(parse_genomes_fasta, genomes)
    taxmaps_json = if isnothing(taxmaps)
        empty_taxmaps(genomes_json)
    else
        ncbi_dict = isnothing(ncbi) ? nothing : open_perhaps_gzipped(parse_ncbi_taxes, ncbi)
        open_perhaps_gzipped(io -> parse_taxmaps(io, ncbi_dict), taxmaps)
    end
    seq_json = open_perhaps_gzipped(parse_sequences_blast, sequences)
    json = ReferenceJSON(1, genomes_json, seq_json, taxmaps_json)
    Reference(json)
end

# Sequences: BLAST + (name, len) of unblasted seqs??
# Taxmaps: A manually specified file.
#     Format: ??
# Genomes:
#     [(flagset, directory) ... ] OR
#     



# (NCBI + taxmaps_file) or genomes => TAXMAPS_JSON_T
# (
#    genomes_dir, virus_dir or
#    raw genomes_json_t?
# ) => GENOMES_JSON_T
# blast => SEQUENCES_JSON_T

sequences(args::SeqArgs)::SEQUENCES_JSON_T = sequences(args.paths)
sequences(x::@NamedTuple{json::String}) = open(io -> JSON3.read(io, SEQUENCES_JSON_T), x.json)

# TODO: Parse the two files concurrently
function sequences(x::@NamedTuple{blast::String, fasta::String})
    d = open(parse_sequences_blast, x.blast)
    for (identifier, len) in open_perhaps_gzipped(parse_sequences_fasta, x.fasta)
        get!(d, identifier, (len, Tuple{String, Int, Int}[]))
    end
    [(q, L, t) for (q, (L, t)) in d]
end

genomes(x::GenomeArgs)::GENOMES_JSON_T = genomes(x.x)

# TODO: This method may need some conversion to support the integer FlagSet, which users
# may not know how to set.
genomes(x::@NamedTuple{json::String}) = open(io -> JSON3.read(io, GENOMES_JSON_T), x.json)

function genomes(directories::Vector{Pair{FlagSet, String}})
    result = GENOMES_JSON_T()
    for (flagset, directory) in directories
        files = filter!(readdir(directory; join=true)) do path
            any(endswith(path, i) for i in [".fasta", ".fna", ".fa"])
        end
        subdir_result = GENOMES_JSON_T(undef, length(files))
        Threads.@threads for (subdir_index, path) in collect(enumerate(files))
            (genome_name, _) = splitext(basename(path))
            io = open(path; lock=false)
            io = endswith(path, ".gz") ? GzipDecompressorStream(io) : io
            sources = FASTA.Reader(io; copy=false) do reader
                map(reader) do record
                    (FASTA.identifier(record), FASTA.seqsize(record))
                end
            end
            subdir_result[subdir_index] = (genome_name, flagset.x, sources)
        end
        append!(result, subdir_result)
    end
    result
end