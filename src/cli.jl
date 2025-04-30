# TODO
# Less whitespace in help string - https://github.com/comonicon/Comonicon.jl/issues/256
# In help string, improve the default value printing
# Help for the main function

const DEFAULT_RECALLS = [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
const DEFAULT_PRECISIONS = [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

# TODO: This Base piracy is shady as fuck
Base.tryparse(::Type{Union{Nothing, String}}, x::String) = x

function Base.tryparse(::Type{FlagSet}, x::String)
    Iterators.map(eachsplit(x, ',')) do str
        @something tryparse(Flag, str) parse_flag_error(str)
    end |> FlagSet
end

function parse_flag_error(s::AbstractString)
    return exitwith(
        "Cannot parse as flag: \"$(s)\". " *
            "Available flags are: $(join(instances(Flag), ", "))",
    )
end

function exitwith(s::AbstractString)
    @error String(s)
    return exit(1)
end

function check_file(s::AbstractString, name::AbstractString)
    @debug lazy"Checking existence of file \"$(s)\""
    return isfile(s) ? s : exitwith("$(name) does not point to an existing file: \"$(s)\"")
end

function check_dir(s::AbstractString, name::AbstractString)
    @debug lazy"Checking existence of dir \"$(s)\""
    return isdir(s) ? s : exitwith("$(name) does not point to an existing directory: \"$(s)\"")
end

parentdir(dir::AbstractString) = dirname(rstrip(dir, ['\\', '/']))

function mkdir_checked(path::String, exists_ok::Bool)
    @debug lazy"Attempting to make directory at \"$(path)\""
    if exists_ok
        isfile(path) && exitwith("Output directory \"$(path)\" is an existing file")
    else
        ispath(path) && exitwith("Output directory \"$(path)\" already exists")
    end
    parent = parentdir(path)
    !isempty(parent) &&
        !isdir(parent) &&
        exitwith("Outdir's parent \"$(parent)\" is not an existing directory")
    return mkdir(path)
end

struct Thresholds
    x::Union{Vector{Float64}, Nothing}
end

const default_thresholds = Thresholds(nothing)

function Base.tryparse(::Type{Thresholds}, s::String)
    v = map(eachsplit(s, ',')) do i
        @something tryparse(Float64, i) exitwith("Cannot parse as float: \"$(i)\"")
    end
    return Thresholds(v)
end


"""
# Intro
Benchmark a set of bins agains a reference

# Args
- `outdir`: Path of output directory to create (must not exist)
- `ref`: Path to reference JSON file (see the `makeref subcommand`)
- `bins`: Path to one or more TSV files of bins.

# Options
- `-s, --sep`: Binsplit separator; delimiter to split bins
- `--minsize`: Remove bins with shorter breadth
- `--minseqs`: Remove bins fewer sequences
- `--keep-flags`: Only keep genomes with all these flags set.
  Takes a comma-separated list of flags, e.g. "virus,plasmid"
- `--remove-flags`: Remove genomes with any of these flags set. See above.
- `--recalls`: Comma-separated list of recall thresholds to use
- `--precisions`: Comma-separated list of precision thresholds to use

# Flags
- `--intersect`: Allow the same sequence to be in multiple bins
- `--quiet`: Disable non-error logging to stderr
"""
Comonicon.@cast function bench(
        outdir::String,
        ref::String,
        bins...; # name=path strings, or a single path
        sep::Union{String, Nothing} = nothing,
        minsize::Int = 1,
        minseqs::Int = 1,
        keep_flags::FlagSet = FlagSet(),
        remove_flags::FlagSet = FlagSet(),
        recalls::Thresholds = default_thresholds,
        precisions::Thresholds = default_thresholds,
        intersect::Bool = false,
        quiet::Bool = false,
    )
    set_global_logger!(; quiet)
    @debug "Checking validity of input paths"
    pairs::Union{String, Vector{Pair{String, String}}} = if length(bins) == 1
        check_file(only(bins), "Binning file")
    else
        v = map(collect(bins)) do i
            p = split(i, '='; limit = 2)
            length(p) == 1 && exitwith(
                "If multiple bins files is passed, pass them as e.g. bin1=path_to_bin1.tsv bin2=path_to_bin2.tsv",
            )
            String(p[1]) => String(check_file(p[2], p[1]))
        end
        if !allunique(map(first, v))
            seen = Set{String}()
            for (name, _) in v
                if name in seen
                    exitwith(
                        lazy"Multiple binnings have been assigned the same name \"$(name)\"",
                    )
                end
                push!(seen, name)
            end
        end
        v
    end
    check_file(ref, "Reference")
    @debug "All files checked"
    outdir = abspath(expanduser(outdir))
    mkdir_checked(outdir, false)
    @debug "Creating log file"
    set_global_logger!(joinpath(outdir, "log.txt"); quiet)
    start_info()
    @info "Running subcommand `bench`"
    settings = OutputSettings(recalls.x, precisions.x)
    @info "Loading reference"
    refdata = open(read, ref)
    reference = Reference(IOBuffer(refdata))
    @debug "Reference loaded"
    @debug "Computing reference CRC32c"
    checksum = crc32c(refdata)
    @info "Reference file CRC32c: $(string(checksum, base = 16))"
    resize!(refdata, 0) # save memory
    bin_paths::Vector{String} = pairs isa String ? [pairs] : last.(pairs)
    binnings = Vector{Binning}(undef, length(bin_paths))
    @info "Loading binning(s) from TSV file(s)"
    @debug "Loading with $(nthreads()) threads"
    @threads :greedy for (i, path) in enumerate(bin_paths)
        @debug "Loading binning at path \"$(path)\""
        binnings[i] = Binning(
            path,
            reference;
            binsplit_separator = sep,
            min_size = minsize,
            min_seqs = minseqs,
            disjoint = !intersect,
            recalls = settings.recalls,
            precisions = settings.precisions,
            filter_genomes = g ->
            isdisjoint(flags(g), remove_flags) && issubset(keep_flags, flags(g)),
        )
    end
    @debug "Done loading binnings"
    if pairs isa String
        populate_output(outdir, only(binnings), settings)
    else
        names_binnings = collect(zip(Iterators.map(first, pairs), binnings))
        @threads for (name, binning) in names_binnings
            subdir = joinpath(outdir, name)
            mkdir_checked(subdir, false)
            populate_output(subdir, binning, settings)
        end
    end
    @info "Done"
end

struct SeqArgs
    paths::Union{@NamedTuple{json::String}, @NamedTuple{seq_mapping::String, fasta::String}}

    function SeqArgs(;
            json::Union{String, Nothing},
            seq_mapping::Union{String, Nothing},
            fasta::Union{String, Nothing},
        )
        return if !isnothing(json)
            (isnothing(seq_mapping) && isnothing(fasta)) ||
                exitwith("If seq-json is set, neither seq-blast nor seq-fasta can be")
            json = expanduser(json)
            check_file(json, "seq-json")
            new((; json))
        elseif isnothing(fasta) || isnothing(seq_mapping)
            exitwith("If seq-json is not set seq-fasta and seq-mapping must be")
        else
            fasta = expanduser(fasta)
            check_file(fasta, "seq-fasta")
            seq_mapping = expanduser(seq_mapping)
            check_file(seq_mapping, "seq-mapping")
            new((; seq_mapping, fasta))
        end
    end
end

# virus+plasmid=/path/to/these,organism=/path/to/those
function parse_genomes_dir(s::String)::Vector{Pair{FlagSet, String}}
    return map(split(strip(s), ',')) do segment
        p = findfirst(==(UInt8('=')), codeunits(segment))
        isnothing(p) && exitwith(
            "No \"=\" symbol found in --genome-directories cli argument \"$(segment)\"",
        )
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
        return if !isnothing(json)
            @debug "Validating genome args as JSON"
            isnothing(directories) ||
                exitwith("If genome-json is set, genome-directories must not be set")
            json = expanduser(json)
            check_file(json, "genome-json")
            new((; json))
        else
            @debug "Validating genome args from directories"
            isnothing(directories) &&
                exitwith("If genome-json is not set, genome-directories must be set")
            v = parse_genomes_dir(directories)
            seen_flagsets = Set{FlagSet}()
            for (flagset, path) in v
                flagset ∈ seen_flagsets &&
                    exitwith("Flagset passed twice: $(join(flagset, '+'))")
                push!(seen_flagsets, flagset)
                check_dir(path, "Genome subdirectory")
            end
            new(v)
        end
    end
end

struct TaxArgs
    paths::Union{
        Nothing,
        @NamedTuple{json::String},
        @NamedTuple{tax::String, ncbi::Union{Nothing, String}}
    }

    function TaxArgs(;
            json::Union{String, Nothing},
            tax::Union{String, Nothing},
            ncbi::Union{String, Nothing},
        )
        return if isnothing(json) && isnothing(tax) && isnothing(ncbi)
            # If no paths are passed, we create a dummy taxonomy
            new(nothing)
        elseif !(isnothing(json) ⊻ isnothing(tax))
            exitwith("Exactly one of tax-json or tax must be set")
        elseif !isnothing(json)
            isnothing(ncbi) || exitwith("If tax-json is set, tax-ncbi must not be set")
            json = expanduser(json)
            check_file(json, "tax-json")
            new((; json))
        else
            @assert !isnothing(tax)
            tax = expanduser(tax)
            check_file(tax, "tax")
            if !isnothing(ncbi)
                ncbi = expanduser(ncbi)
                check_file(ncbi, "tax-ncbi")
            end
            new(@NamedTuple{tax::String, ncbi::Union{Nothing, String}}((; tax, ncbi)))
        end
    end
end

"""
# Intro
Create a new reference JSON file. See more details in the documentation.

# Args
- `outdir`: Path to output directory

# Options
- `--seq-mapping`: Mapping of sequences to genomes in 4-col TSV file
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
- `--quiet`: Disable non-error logging to stderr
"""
Comonicon.@cast function makeref(
        outdir::String;
        # Mapping of sequences to genomes (four-column format, see docs)
        seq_mapping::Union{String, Nothing} = nothing,
        # This is for unmapped sequences, so they still appear in reference
        seq_fasta::Union{String, Nothing} = nothing,
        # ... or they can pass in the JSON directly
        seq_json::Union{String, Nothing} = nothing,

        # The special tax file
        tax::Union{String, Nothing} = nothing,
        # Path to dump of NCBI - used to parse the `id=412` elements in `tax`
        tax_ncbi::Union{String, Nothing} = nothing,
        # ... or they can pass in the JSON directly
        tax_json::Union{String, Nothing} = nothing,

        # A string of the format
        # virus+plasmid=/path/to/phasmids,organisms=/path/to/organisms
        # I.e. (flagset, path) pairs comma-sep, with = to delimit the two, and +
        # to delimit flags in the flagset.
        genome_directories::Union{String, Nothing} = nothing,
        # ... or they can pass in the JSON directly
        genome_json::Union{String, Nothing} = nothing,

        # Do not error if output directory already exists
        overwrite::Bool = false,

        # Disable logging
        quiet::Bool = false,
    )
    set_global_logger!(; quiet)
    @debug "Validating sequence args"
    seq_args = SeqArgs(; json = seq_json, seq_mapping = seq_mapping, fasta = seq_fasta)
    @debug "Validating genome args"
    genome_args = GenomeArgs(; json = genome_json, directories = genome_directories)
    @debug "Validating taxonomy args"
    tax_args = TaxArgs(; json = tax_json, ncbi = tax_ncbi, tax)
    @debug "All arguments validated"
    outdir = abspath(expanduser(outdir))
    mkdir_checked(outdir, overwrite)
    @debug "Creating log file"
    set_global_logger!(joinpath(outdir, "log.txt"); quiet)
    start_info()
    @info "Running subcommand `makeref`"

    gns = genomes(genome_args)
    @debug "Writing genomes.json file"
    isnothing(genome_json) &&
        open(io -> JSON3.write(io, gns), joinpath(outdir, "genomes.json"), "w")

    taxes = taxonomy(tax_args, gns)
    @debug "Writing taxonomy.json file"
    isnothing(tax) ||
        open(io -> JSON3.write(io, taxes), joinpath(outdir, "taxonomy.json"), "w")

    seqs = sequences(seq_args, gns)
    @debug "Writing seqs.json file"
    isnothing(seq_json) &&
        open(io -> JSON3.write(io, seqs), joinpath(outdir, "seqs.json"), "w")

    @info "Constructing Reference object"
    @debug "Constructing ReferenceJSON"
    json_struct = BBB.ReferenceJSON(BBB.JSON_VERSION, gns, seqs, taxes)
    @debug "Constructing Reference from ReferenceJSON"
    ref = Reference(json_struct)
    @debug "Writing reference object"
    open(joinpath(outdir, "reference.json"), "w") do io
        BBB.save(io, ref)
    end
    @info "Done"
end
