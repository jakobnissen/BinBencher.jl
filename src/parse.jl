function taxonomy(x::TaxArgs, genomes::GENOMES_JSON_T)::TAXMAPS_JSON_T
    paths = x.paths
    if isnothing(paths)
        [[(genome_name, "top") for (genome_name, _, _) in genomes]]
    else
        taxonomy(paths)
    end
end

taxonomy(x::@NamedTuple{json::String}) = open(io -> JSON3.read(io, TAXMAPS_JSON_T), x.json)

function taxonomy(x::@NamedTuple{tax::String, ncbi::String})
    (; tax, ncbi) = x
    ncbi_dict = open_perhaps_gzipped(parse_ncbi_taxes, ncbi)
    open_perhaps_gzipped(io -> parse_taxmaps(io, ncbi_dict), tax)
end

# child_id => (rank, parent_id, name)
const NCBI_DICT_T = Dict{Int32, Tuple{Int32, Int32, String}}

# child_id  rank    parent_id   name
function parse_ncbi_taxes(io::IO)::NCBI_DICT_T
    result = NCBI_DICT_T()
    for (lineno, line) in enumerate(eachline(io))
        stripped = rstrip(line)
        isempty(stripped) && continue
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

# OR: Maybe easier to set:
# genome species genus [...] id=3134
# Minimally two columns. Maximally X columns. Id must be last, and map to the clade above

# Max level: Find maximal level specified on line without id=. (non-ncbi lines)
# Check all non-ncbi lines have same level
# Check NCBI lines have at most same level
# Fill out NCBI lines to match that level

# strain   child   parent
# class    child   id=3143

# TODO: How should the user specify: I don't care about higher levels than this?
# TODO: If the user uses id=XXX, how do we prevent this function from taking too many levels
# when the user otherwise only uses X levels?


function parse_taxmaps(io::IO, ncbi::NCBI_DICT_T)::TAXMAPS_JSON_T
    result = Dict(i => Dict{String, String}() for i in 0:7)
    for (line_number, line) in enumerate(eachline(io))
        fields = split(rstrip(line), '\t')
        length(fields) == 3 || error(
            lazy"In taxmaps file, on line $(line_number), " *
            "expected 3 tab-separated fields, got $(length(fields))",
        )
        (rank_str, child, parent) = fields
        rank = get(RANK_BY_NAME, lowercase(rank_str), 8)
        rank == 8 && error(
            lazy"Invalid rank name: \"$(rank_str)\". Valid names are:\n" *
            lazy"$(join(RANKS[1:end-1], \", \"))",
        )
        if !startswith(parent, "id=")
            existing = get!(result[rank], child, parent)
            existing == parent || error(
                lazy"At rank \"$(rank_str)\", child \"$(child)\" is has two distinct parents",
            )
        else
            parent_id = parse(Int, parent[4:end])
            for rank in rank:6
                lookup = get(ncbi, parent_id, nothing)
                isnothing(lookup) && error(
                    lazy"NCBI taxid $(parent_id) is not a known NCBI taxid. " *
                    "Please note that we only accept taxids with a rank " *
                    lazy"of $(join(RANKS[1:end-1], \", \")). Verify the NCBI taxid, " *
                    "or manually include the taxonomic levels for child $(child)",
                ) # we have only the normal ranks
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
    # Remove empty dicts in result
    vector = [[(k, v) for (k, v) in d] for d in [result[i] for i in 0:7]]
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

const BLAST_ROW = @NamedTuple{query::String, subject::String, sstart::Int, send::Int}

const DEFAULT_ROW = @NamedTuple{
    query_index::Int,
    subject_index::Int,
    ident::Float64,
    qstart::Int,
    qend::Int,
    sstart::Int,
    send::Int,
}

"""
    parse_default_blast(io::IO)::IO

Parse the default format of BLAST -outfmt 6.
This also handles the case where the queries are rotated relative to the sources,
and also automatically filters the hits.
Output is an IO object you can call parse_sequences_blast on
"""
function parse_default_blast(
    io::IO,
    _query_lengths::Dict{<:AbstractString, <:Integer},
    # nothing, or ??? TODO
)::Vector{BLAST_ROW}
    # Map from queries and subjects to ints
    query_to_int = Dict(k => i for (i, k) in enumerate(keys(_query_lengths)))
    query_lengths = collect(values(_query_lengths))
    subject_to_int = Dict{String, Int}() # we build as we go

    # Load all hits in and sort them by (query, subject).
    # This allows us to group them by q/s pair by looping over them,
    # without having to store them in some nested dict or something.
    fields = SubString{String}[]
    hits = DEFAULT_ROW[]
    for line in eachline(io)
        row = parse_default_blast_row(line, fields, query_to_int, subject_to_int)
        # Some basic filtering criteria: Identity must be high
        row.ident ≥ 0.95 && push!(hits, row)
    end

    # No hits: Return early
    result = IOBuffer()
    isempty(hits) && return BLAST_ROW[]

    # Use these mappings to get back to strings from the indices
    int_to_query = Dict(v => k for (k, v) in query_to_int)
    int_to_subject = Dict(v => k for (k, v) in subject_to_int)
    empty!(subject_to_int) # save memory
    empty!(query_to_int)

    # It should be mostly sorted already, so this should be decently fast
    sort!(hits; by=i -> (i[1], -i[3]))

    # Group by query and process in that group in order to
    # get only the best hits per query.
    # (if a contig blasts to two subjects, one much better than the other,
    # then it belongs to the first)
    query_buffer = DEFAULT_ROW[]
    last_query = first(hits).query_index
    for row in hits
        if row.query_index != last_query
            process_query!(query_buffer, query_lengths, result, int_to_query, int_to_subject)
            empty!(query_buffer)
            last_query = row.query_index
        end
        push!(query_buffer, row)
    end
    process_query!(query_buffer, query_lengths, result, int_to_query, int_to_subject)
    return parse_sequences_blast(seekstart(result))
end

function process_query!(
    # NB: All hits from one query, in descending order by identity
    query_buffer::Vector{DEFAULT_ROW},
    query_lengths::Vector{Int},
    result::IOBuffer,
    int_to_query::Dict{Int, String},
    int_to_subject::Dict{Int, String},
)::Nothing
    # Filter for identity: At least 97%, no more than 1% ANI or 1.5x
    # more differences than best hit
    best_identity = first(query_buffer).ident
    min_identity = max(0.97, best_identity - 0.01, 1 - 1.5 * (1 - best_identity))
    filter!(i -> i.ident ≥ min_identity, query_buffer)
    isempty(query_buffer) && return nothing

    should_circularize = false
    if should_circularize # TODO: Make this into a non-dead branch!
        circularize!(query_buffer)
    end

    # Filter for query coverage
    filter!(query_buffer) do i
        (qs, qe) = minmax(i.qstart, i.qend)
        qcov = ((qe - qs) + 1) / query_lengths[i.query_index]
        qcov ≥ 0.9
    end

    # circularize will call this function
    if !should_circularize
        make_forward!(query_buffer)
    end

    # Flip the reverse oriented ones to be fw
    for i in eachindex(query_buffer)
        d = query_buffer[i]
        (ss, se) = minmax(d.sstart, d.send)
        query_buffer[i] = (;d..., sstart=ss, send=se)
    end

    # Write result out
    for i in query_buffer
        print(
            result,
            int_to_query[i.query_index], '\t',
            int_to_subject[i.subject_index], '\t',
            i.sstart, '\t',
            i.send, '\n'
        )
    end
    nothing
end

function circularize!(v::Vector{DEFAULT_ROW}) # TODO
    1 + "abc" # show up in static analysis
    # Return early if < 2
    # Split to fw and rv hits
    # Convert to forward
    # Then actually circularize
end

function make_forward!(v::Vector{DEFAULT_ROW})
    for i in eachindex(v)
        d = v[i]
        (ss, se) = minmax(d.sstart, d.send)
        v[i] = (;d..., sstart=ss, send=se)
    end
    v
end

function parse_default_blast_row(
    line::Union{String, SubString{String}},
    fields::Vector{SubString{String}},
    query_to_int::Dict{String, Int},
    subject_to_int::Dict{String, Int},
)::DEFAULT_ROW
    empty!(fields)
    for field in eachsplit(line, '\t')
        push!(fields, field)
    end
    subject = fields[2]
    subject_index = get!(subject_to_int, subject, length(subject_to_int) + 1)
    # Default fields in BLAST:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    if length(fields) != 12
        error("By default, tabular BLAST should have 12 fields")
    end
    (qs, qe) = (parse(Int, fields[7]), parse(Int, fields[8]))
    (qstart, qend) = minmax(qs, qe)
    (;
        query_index = query_to_int[fields[1]],
        subject_index,
        ident = parse(Float64, fields[3]) * 0.01, # convert from percent
        qstart,
        qend,
        sstart = parse(Int, fields[9]),
        send = parse(Int, fields[10]),
    )
end

# qseqid     sseqid    sstart    send
function parse_sequences_blast(io::IO)::Vector{BLAST_ROW}
    result = BLAST_ROW[]
    for (lineno, line) in enumerate(eachline(io))
        stripped = rstrip(line)
        isempty(stripped) && continue
        m = match(r"^(\S+)\t(\S+)\t(\d+)\t(\d+)\t?", stripped)
        m === nothing && error(
            "On line $(lineno), expected line of the form " *
            "QSEQID\\tSSEQID\\tSSTART\\tSEND[\\t FIELDS...] " *
            lazy"but got \"$(repr(stripped))\"",
        )
        (query, subject, sstart_str, send_str) = map(something, m.captures)
        sstart = parse(Int, sstart_str)
        send = parse(Int, send_str)
        push!(result, (; query, subject, sstart, send))
    end
    result
end

function sequences(
    blast::Vector{BLAST_ROW},
    fasta::Vector{Tuple{String, Int}},
    subject_lengths::Dict{String, Int},
)::SEQUENCES_JSON_T
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
                "In BLAST file, a BLAST hit spans spans $(start)-$(stop) even though " *
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
                    "This is more than 125% of the query length $(qlen).",
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
    sequences(args.paths, genomes::GENOMES_JSON_T)
end

function sequences(x::@NamedTuple{json::String}, ::GENOMES_JSON_T)
    open(io -> JSON3.read(io, SEQUENCES_JSON_T), x.json)
end

# TODO: Check if this is type stable
function sequences(
    x::@NamedTuple{default_blast::String, fasta::String},
    genomes::GENOMES_JSON_T
)
    fasta = open_perhaps_gzipped(parse_sequences_fasta, x.fasta)
    blast = open(io -> parse_default_blast(io, Dict(fasta)), x.default_blast)
    subject_lengths = Dict{String, Int}()
    for (_, _, subjects) in genomes, (subject, len) in subjects
        haskey(subject_lengths, subject) &&
            exitwith("Duplicate subject name in genomes: \"$(subject)\"")
        subject_lengths[subject] = len
    end
    sequences(blast, fasta, subject_lengths)
end

# TODO: Check if this is type stable
function sequences(
    x::@NamedTuple{custom_blast::String, fasta::String},
    genomes::GENOMES_JSON_T
)
    a = Threads.@spawn open(parse_sequences_blast, x.custom_blast)
    b = Threads.@spawn open_perhaps_gzipped(parse_sequences_fasta, x.fasta)
    blast = fetch(a)
    fasta = fetch(b)
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
genomes(x::@NamedTuple{json::String}) = open(io -> JSON3.read(io, GENOMES_JSON_T), x.json)

function genomes(directories::Vector{Pair{FlagSet, String}})
    result = GENOMES_JSON_T()
    for (flagset, directory) in directories
        # TODO: Add warning if any file is not hidden and not of this format
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
