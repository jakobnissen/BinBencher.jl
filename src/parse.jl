function taxonomy(x::TaxArgs, genomes::GENOMES_JSON_T)::TAXMAPS_JSON_T
    paths = x.paths
    if isnothing(paths)
        @info "No taxonomy passed - skipping loading it"
        [[(genome_name, "top") for (genome_name, _, _) in genomes]]
    else
        @info "Loading taxonomy"
        taxonomy(paths)
    end
end

function taxonomy(x::@NamedTuple{json::String})
    @info "Loading taxonomy from JSON file at \"$(x.json)\""
    open(io -> JSON3.read(io, TAXMAPS_JSON_T), x.json)
end

function taxonomy(x::@NamedTuple{tax::String, ncbi::String})
    (; tax, ncbi) = x
    @info "Loading taxonomy file from TSV and NCBI at \"$(tax)\" and \"$(ncbi)\""
    ncbi_dict = open_perhaps_gzipped(parse_ncbi_taxes, ncbi)
    open_perhaps_gzipped(io -> parse_taxmaps(io, ncbi_dict), tax)
end

# child_id => (rank, parent_id, name)
const NCBI_DICT_T = Dict{Int32, Tuple{Int32, Int32, String}}

# child_id  rank    parent_id   name
function parse_ncbi_taxes(io::IO)::NCBI_DICT_T
    @debug "Parsing NCBI taxonomy"
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
    @debug "Parsing taxonomy TSV file"
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

const EXPECTED_SEQ_TSV_HEADER = "sequence\tsource\tstart\tend"

# qseqid     sseqid    sstart    send
function parse_sequences_tsv(io::IO)::Vector{BLAST_ROW}
    result = BLAST_ROW[]
    lines = Iterators.filter(i -> !isempty(last(i)), enumerate(eachline(io)))
    lines = let
        peeled = Iterators.peel(lines)
        isnothing(peeled) && exitwith(lazy"Empty sequence TSV file")
        ((_, header), rest) = peeled
        header = rstrip(header)
        if header != EXPECTED_SEQ_TSV_HEADER
            exitwith(
                lazy"In sequence TSV file, expected header \"$(EXPECTED_SEQ_TSV_HEADER)\"" *
                lazy", got \"$(header)\"",
            )
        end
        rest
    end
    for (lineno, line) in lines
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
    @info "Loading sequences from JSON file $(x.json)"
    open(io -> JSON3.read(io, SEQUENCES_JSON_T), x.json)
end

function sequences(
    x::@NamedTuple{seq_mapping::String, fasta::String},
    genomes::GENOMES_JSON_T,
)
    @info "Loading sequences from seq mapping $(x.seq_mapping)"
    a = @spawn open(parse_sequences_tsv, x.seq_mapping)
    b = @spawn open_perhaps_gzipped(parse_sequences_fasta, x.fasta)
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
