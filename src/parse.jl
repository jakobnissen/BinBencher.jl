# Iterate nonempty, nonheader lines of IO, checking the header is the argument passed 
function tsv_iterator(io::IO, name::AbstractString, header::String)
    lines = enumerate(eachline(io))
    lines = Iterators.filter(i -> !isempty(rstrip(last(i))), lines)
    peeled = Iterators.peel(lines)
    isnothing(peeled) && exitwith("Found no TSV header when parsing \"$(name)\"")
    ((_, obs_header), rest) = peeled
    if rstrip(header) != rstrip(obs_header)
        exitwith("In $(name), did not find expected header \"$(header)\"")
    end
    rest
end

function taxonomy(x::TaxArgs, genomes::GENOMES_JSON_T)::TAXMAPS_JSON_T
    paths = x.paths
    if isnothing(paths)
        @info "No taxonomy passed - skipping loading it"
        # BinBencher requires at least one Clade which is the ancestor of all genomes,
        # so we make an artificial one if none is passed.
        [[(genome_name, "top") for (genome_name, _, _) in genomes]]
    else
        @debug "Loading taxonomy"
        taxonomy(paths)
    end
end

function taxonomy(x::@NamedTuple{json::String})
    @info "Loading taxonomy from JSON file at \"$(x.json)\""
    open(io -> JSON3.read(io, TAXMAPS_JSON_T), x.json)
end

function taxonomy(x::@NamedTuple{tax::String, ncbi::Union{Nothing, String}})
    (; tax, ncbi) = x
    ncbi_dict = if !isnothing(ncbi)
        @info "Loading NCBI taxfile at \"$(ncbi)\""
        open_perhaps_gzipped(parse_ncbi_taxes, ncbi)
    else
        @info "No NCBI taxfile passed"
        nothing
    end
    @info "Loading taxonomy file from TSV at \"$(tax)\""
    open_perhaps_gzipped(io -> parse_taxmaps(io, ncbi_dict), tax)
end

# child_id => (rank, parent_id, name)
const NCBI_DICT_T = Dict{Int32, Tuple{Int32, Int32, String}}

const NCBI_HEADER = "child_id\trank\tparent_id\tname"

function parse_ncbi_taxes(io::IO)::NCBI_DICT_T
    @debug "Parsing NCBI taxonomy"
    result = NCBI_DICT_T()
    for (lineno, line) in tsv_iterator(io, "NCBI file", NCBI_HEADER)
        stripped = rstrip(line)
        m = match(r"^(\d+)\t(\d+)\t(\d+)\t([^\t\r\n]+)$", stripped)
        m === nothing && error(
            "On line $(lineno), expected line of the form " *
            "$(raw"^(\d+)\t(\d+)\t(\d+)\t([^\t\r\n]+)$"), " *
            lazy"but got \"$(stripped)\"",
        )
        (child_str, rank_str, parent_str, name) = eachsplit(stripped, '\t')
        result[parse(Int32, child_str)] =
            (parse(Int32, rank_str), parse(Int32, parent_str), name)
    end
    result
end

const TAX_TSV_HEADER = "rank\tchild\tparent"

function parse_taxmaps(io::IO, ncbi::Union{Nothing, NCBI_DICT_T})::TAXMAPS_JSON_T
    # Rank => (child => parent)
    result = Dict(i => Dict{String, String}() for i in 0:7)
    # If an NCBI id line (e.g. id=1234) is seen, we autofill all higher ranks.
    # but we don't want to fill in ranks higher than the highest manually specified
    # rank, because all genomes must be filled to the same rank.
    largest_parent_rank_seen = nothing

    for (line_number, line) in tsv_iterator(io, "Taxonomy TSV file", TAX_TSV_HEADER)
        fields = split(rstrip(line), '\t')
        length(fields) == 3 || error(
            lazy"In taxmaps file, on line $(line_number), " *
            "expected 3 tab-separated fields, got $(length(fields))",
        )
        (rank_str, child, parent) = fields
        # For now, we demand that the user uses the classical 7 taxonomic ranks
        # plus 'strain' (the latter for genomes).
        rank = get(RANK_BY_NAME, lowercase(rank_str), 8)
        rank == 8 && error(
            lazy"In taxmaps file, on line $(line_number), " *
            lazy"Invalid rank name: \"$(rank_str)\". Valid names are:\n" *
            lazy"$(join(RANKS[1:end-1], \", \"))",
        )
        if !startswith(parent, "id=")
            existing = get!(result[rank], child, parent)
            existing == parent || error(
                lazy"In taxmaps file, on line $(line_number), " *
                lazy"At rank \"$(rank_str)\", child \"$(child)\" has two distinct parents",
            )
            largest_parent_rank_seen = if isnothing(largest_parent_rank_seen)
                # The rank listed is the child rank, so parent rank is +1
                rank + 1
            else
                max(largest_parent_rank_seen, rank + 1)
            end
        else
            if isnothing(ncbi)
                error(
                    lazy"In tax file, on line $(line_number), the parent is given " *
                    "as an NCBI tax id. However, the NCBI tax file was not provided.",
                )
            end
            parent_id = parse(Int, parent[4:end])
            # We fill in all ranks, then we trim away the unused ones at the end
            # once we've traversed the whole file and know the highest rank manually
            # given
            for rank in rank:6
                lookup = get(ncbi, parent_id, nothing)
                isnothing(lookup) && error(
                    lazy"In taxmaps file, on line $(line_number), " *
                    lazy"NCBI taxid $(parent_id) is not a known NCBI taxid. " *
                    "Please note that we only accept taxids with a rank " *
                    lazy"of $(join(RANKS[1:end-1], \", \")). Verify the NCBI taxid, " *
                    "or manually include the taxonomic levels for child $(child)",
                ) # we have only the normal ranks (the seven from species to domain)
                (parent_rank, grandparent_id, parent) = lookup
                parent_rank == rank + 1 || error(
                    lazy"At rank $(rank_str), NCBI taxid $(parent_id) is listed as a parent. " *
                    lazy"However, the parent of a $(rank_str) must be a $(RANKS[rank+1]), " *
                    "whereas NCBI tax $(parent_id) is a $(RANKS[parent_rank])",
                )
                existing = get!(result[rank], child, parent)
                existing == parent || error(
                    lazy"At rank $(RANKS[rank]), child \"$(child)\" is listed " *
                    "with multiple distinct parents.",
                )
                child = parent
                parent_id = grandparent_id
            end
        end
    end
    # If NO ranks were specified manually, we keep all ranks.
    top_rank_to_keep = if isnothing(largest_parent_rank_seen)
        7
    else
        @assert largest_parent_rank_seen < 9
        largest_parent_rank_seen - 1
    end

    # Remove empty dicts in result
    vector = [[(k, v) for (k, v) in d] for d in [result[i] for i in 0:top_rank_to_keep]]
    while isempty(last(vector))
        pop!(vector)
        isempty(vector) && return vector
    end
    # Add top clade if there isn't one
    last_parents = unique(map(last, last(vector)))
    if length(last_parents) > 1
        push!(vector, [(p, "top") for p in last_parents])
    end
    vector
end

const BLAST_ROW = @NamedTuple{query::String, subject::String, sstart::Int, send::Int}

const EXPECTED_SEQ_TSV_HEADER = "sequence\tsource\tstart\tend"

function parse_sequences_tsv(io::IO)::Vector{BLAST_ROW}
    result = BLAST_ROW[]
    for (lineno, line) in tsv_iterator(io, "Sequence TSV file", EXPECTED_SEQ_TSV_HEADER)
        fields = split(line, '\t')
        if length(fields) != 4
            exitwith(
                lazy"On line $(lineno) in sequence TSV file, did not get 4 tab-sep columns",
            )
        end
        (sequence, target, start, stop) = fields
        sstart = tryparse(Int, start)
        send = tryparse(Int, stop)
        if isnothing(sstart) || isnothing(send)
            exitwith(
                lazy"On line $(lineno) in sequence TSV files, columns 3 and 4 " *
                "cannot be parsed to `Int`.",
            )
        end
        push!(result, (; query=String(sequence), subject=String(target), sstart, send))
    end
    result
end

function sequences(
    blast::Vector{BLAST_ROW},
    fasta::Vector{Tuple{String, Int}},
    subject_lengths::Dict{String, Int},
)::SEQUENCES_JSON_T
    @debug "Constructing sequence object from BLAST results and FASTA"
    has_warned = false
    dict = Dict{String, Tuple{Int, Vector{Tuple{String, Int, Int}}}}()
    for (identifier, len) in fasta
        dict[identifier] = (len, Tuple{String, Int, Int}[])
    end
    for (; query, subject, sstart, send) in blast
        slen = get(subject_lengths, subject) do
            exitwith(
                "BLAST file contains subject \"$(subject)\", but this was not found " *
                "when parsing the genomes. Make sure the BLAST file was created by BLASTing " *
                "against the genome sequences.",
            )
        end
        (qlen, targets) = get(dict, query) do
            exitwith(
                "BLAST file contains query \"$(query)\", but this was not found when " *
                "parsing the sequence FASTA file. Make sure the BLAST file was created " *
                "using the same sequence FASTA file as query",
            )
        end
        # Some users might copy the data in this file from a BLAST output file.
        # BLAST - and other similar tools - will often set sstart > send if the
        # sequences are reverse-complemented.
        (start, stop) = minmax(sstart, send)
        # We allow users to BLAST against a set of subjects where each BLAST subject
        # is the subjects concatenated with themselves. This enables alignments to
        # "wrap" around the original subjects, simulating an alignment around a circular
        # subject sequence.
        # We can detect this happens if stop > slen
        if stop ≤ slen
            push!(targets, (subject, start, stop))
        else
            stop > 2 * slen && exitwith(
                "In BLAST file, a BLAST hit spans $(start)-$(stop) even though " *
                "subject \"$(subject)\" only has length $(slen).",
            )
            if start > slen
                push!(targets, (subject, start - slen, stop - slen))
                continue
            end
            span_length = length(start:stop)
            qlen > 500 &&
                span_length > 1.25 * qlen &&
                exitwith(
                    "In BLAST file, query \"$(query)\" has a hit spanning $(start)-$(stop). " *
                    "This is more than 125% of the query length $(qlen), which indicates " *
                    "a bad alignment",
                )
            if !has_warned
                @info (
                    "BLAST hit spans $(start)-$(stop) even though subject \"$(subject)\"" *
                    "only has length $(slen). This will be interpreted as a wrapping alignment. " *
                    "Please see the documentation for more details."
                )
                has_warned = true
            end
            push!(targets, (subject, start, stop - slen))
        end
    end
    [(query, len, targets) for (query, (len, targets)) in dict]
end

function parse_sequences_fasta(io::IO)::Vector{Tuple{String, Int}}
    seen_identifiers = Set{String}()
    FASTA.Reader(io; copy=false) do reader
        map(reader) do record
            identifier = String(FASTA.identifier(record))
            # This is kind of an ugly hack in order to prevent a dual lookup in
            # the hash set. Ideally, Julia ought to support a way to do this check
            # using only one lookup.
            L = length(seen_identifiers)
            push!(seen_identifiers, identifier)
            length(seen_identifiers) == L &&
                exitwith("Duplicate query name in FASTA file: \"$(identifier)\"")
            (identifier, FASTA.seqsize(record))
        end
    end
end

function sequences(args::SeqArgs, genomes::GENOMES_JSON_T)::SEQUENCES_JSON_T
    @debug "Loading sequences"
    sequences(args.paths, genomes::GENOMES_JSON_T)
end

function sequences(x::@NamedTuple{json::String}, ::GENOMES_JSON_T)
    @info "Loading sequences from JSON file $(x.json)"
    open(io -> JSON3.read(io, SEQUENCES_JSON_T), x.json)
end

function sequences(
    x::@NamedTuple{seq_mapping::String, fasta::String},
    genomes::GENOMES_JSON_T,
)
    @info lazy"Loading sequence TSV file from seq mapping $(x.seq_mapping)"
    a = @spawn open(parse_sequences_tsv, x.seq_mapping)
    @info lazy"Loading FASTA file at $(x.fasta)"
    fasta = open_perhaps_gzipped(parse_sequences_fasta, x.fasta)
    blast = fetch(a)::Vector{BLAST_ROW}
    subject_lengths = Dict{String, Int}()
    for (_, _, subjects) in genomes, (subject, len) in subjects
        haskey(subject_lengths, subject) &&
            exitwith("Duplicate subject name in genomes: \"$(subject)\"")
        subject_lengths[subject] = len
    end
    sequences(blast, fasta, subject_lengths)
end

genomes(x::GenomeArgs)::GENOMES_JSON_T = genomes(x.x)

# TODO: This method may need some conversion to support the integer FlagSet, which users
# may not know how to set.
function genomes(x::@NamedTuple{json::String})
    @info lazy"Loading genomes from JSON file $(x.json)"
    open(io -> JSON3.read(io, GENOMES_JSON_T), x.json)
end

function genomes(directories::Vector{Pair{FlagSet, String}})
    @info "Loading genomes from directories"
    result = GENOMES_JSON_T()
    for (flagset, directory) in directories
        @debug "Loading genomes from directory \"$(directory)\""
        # TODO: Add warning if any file is not hidden and not of this format
        files = filter!(readdir(directory; join=true)) do path
            any(endswith(path, i) for i in [".fasta", ".fna", ".fa"])
        end
        @debug "Found $(length(files)) FASTA files in directory \"$(directory)\""
        subdir_result = GENOMES_JSON_T(undef, length(files))
        @threads for (subdir_index, path) in collect(enumerate(files))
            (genome_name, _) = splitext(basename(path))
            io = open(path; lock=false)
            io = endswith(path, ".gz") ? GzipDecompressorStream(io) : io
            sources = FASTA.Reader(io; copy=false) do reader
                map(reader) do record
                    (String(FASTA.identifier(record)), FASTA.seqsize(record))
                end
            end
            subdir_result[subdir_index] = (genome_name, flagset.x, sources)
        end
        append!(result, subdir_result)
    end
    result
end
