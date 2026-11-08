const DEFAULT_RECALLS = [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
const DEFAULT_PRECISIONS = [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

# Parse a comma-separated list of flags, e.g. "virus,plasmid"
function parse_flagset(s::AbstractString)::FlagSet
    return Iterators.map(eachsplit(s, ',')) do str
        @something tryparse(Flag, str) parse_flag_error(str)
    end |> FlagSet
end

function parse_flag_error(s::AbstractString)
    return cli_error(
        "Cannot parse as flag: \"$(s)\". " *
            "Available flags are: $(join(instances(Flag), ", "))",
    )
end

"""
    CLIError(msg::String) <: Exception

Exception thrown when the input to the command-line program is invalid, e.g. a
missing file or a malformed input file.
`main` catches it, logs `msg` as an error, and returns exit code 1.
"""
struct CLIError <: Exception
    msg::String
end

Base.showerror(io::IO, e::CLIError) = print(io, "CLIError: ", e.msg)

cli_error(s::AbstractString) = throw(CLIError(String(s)))

# The backend throws a `BinBencherException` when its input is invalid.
# Rethrow the expected exception type `T` as a `CLIError`, so it's reported as an error
# message instead of crashing with a stacktrace. Other exceptions are rethrown, so that
# new exception types in the backend are noticed and handled deliberately.
function rethrow_as_cli_error(f, ::Type{T}, description::AbstractString) where {T <: BinBencherException}
    return try
        f()
    catch e
        e isa T || rethrow()
        cli_error("Invalid $(description): $(sprint(showerror, e))")
    end
end

function check_file(s::AbstractString, name::AbstractString)
    @debug lazy"Checking existence of file \"$(s)\""
    return isfile(s) ? s : cli_error("$(name) does not point to an existing file: \"$(s)\"")
end

function check_dir(s::AbstractString, name::AbstractString)
    @debug lazy"Checking existence of dir \"$(s)\""
    return isdir(s) ? s : cli_error("$(name) does not point to an existing directory: \"$(s)\"")
end

parentdir(dir::AbstractString) = dirname(rstrip(dir, ['\\', '/']))

function mkdir_checked(path::String, exists_ok::Bool)
    @debug lazy"Attempting to make directory at \"$(path)\""
    if exists_ok
        isfile(path) && cli_error("Output directory \"$(path)\" is an existing file")
    else
        ispath(path) && cli_error("Output directory \"$(path)\" already exists")
    end
    parent = parentdir(path)
    !isempty(parent) &&
        !isdir(parent) &&
        cli_error("Outdir's parent \"$(parent)\" is not an existing directory")
    return mkdir(path)
end

# Parse a comma-separated list of recall or precision thresholds, e.g. "0.5,0.9"
function parse_thresholds(s::AbstractString)::Vector{Float64}
    return map(eachsplit(s, ',')) do i
        @something tryparse(Float64, i) cli_error("Cannot parse as float: \"$(i)\"")
    end
end


"""
    bench(outdir, ref, bins::Vector{String}; kwargs...)

Implementation of the `bench` subcommand. See `BenchArgs` for the meaning of
the arguments.
"""
function bench(
        outdir::String,
        ref::String,
        bins::Vector{String}; # name=path strings, or a single path
        sep::Union{String, Nothing} = nothing,
        minsize::Int = 1,
        minseqs::Int = 1,
        keep_flags::FlagSet = FlagSet(),
        remove_flags::FlagSet = FlagSet(),
        recalls::Union{Vector{Float64}, Nothing} = nothing,
        precisions::Union{Vector{Float64}, Nothing} = nothing,
        intersect::Bool = false,
        quiet::Bool = false,
    )
    length(bins) == 0 && cli_error("Must pass at least one path to bins, got zero")
    set_global_logger!(; quiet)
    @debug "Checking validity of input paths"
    pairs::Union{String, Vector{Pair{String, String}}} = if length(bins) == 1
        check_file(only(bins), "Binning file")
    else
        v = map(bins) do i
            p = split(i, '='; limit = 2)
            length(p) == 1 && cli_error(
                "If multiple bins files is passed, pass them as e.g. bin1=path_to_bin1.tsv bin2=path_to_bin2.tsv",
            )
            String(p[1]) => String(check_file(p[2], p[1]))
        end
        if !allunique(map(first, v))
            seen = Set{String}()
            for (name, _) in v
                if name in seen
                    cli_error(
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
    settings = OutputSettings(recalls, precisions)
    @info "Loading reference"
    refdata = open(read, ref)
    reference = rethrow_as_cli_error(ReferenceError, "reference file \"$(ref)\"") do
        Reference(IOBuffer(refdata))
    end
    @debug "Reference loaded"
    @debug "Computing reference CRC32c"
    checksum = crc32c(refdata)
    @info "Reference file CRC32c: $(string(checksum, base = 16))"
    refdata = nothing # free memory
    bin_paths::Vector{String} = pairs isa String ? [pairs] : last.(pairs)
    binnings = Vector{Binning}(undef, length(bin_paths))
    # Exceptions thrown in the threaded loop would be wrapped in a `TaskFailedException`,
    # so we store errors and check them after
    errors = Vector{Union{Nothing, BinningError}}(nothing, length(bin_paths))
    @info "Loading binning(s) from TSV file(s)"
    @debug "Loading with $(nthreads()) threads"
    @threads :greedy for (i, path) in collect(enumerate(bin_paths))
        @debug "Loading binning at path \"$(path)\""
        try
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
        catch e
            e isa BinningError || rethrow()
            errors[i] = e
        end
    end
    for (path, e) in zip(bin_paths, errors)
        isnothing(e) || cli_error("Invalid binning file \"$(path)\": $(sprint(showerror, e))")
    end
    @debug "Done loading binnings"
    if pairs isa String
        populate_output(outdir, only(binnings), settings)
    else
        # Create the subdirectories before the threaded loop, so a `CLIError` is not
        # thrown from within a task.
        subdirs = map(pairs) do (name, _)
            subdir = joinpath(outdir, name)
            mkdir_checked(subdir, false)
            subdir
        end
        @threads for (subdir, binning) in collect(zip(subdirs, binnings))
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
                cli_error("If seq-json is set, neither seq-mapping nor seq-fasta can be")
            json = expanduser(json)
            check_file(json, "seq-json")
            new((; json))
        elseif isnothing(fasta) || isnothing(seq_mapping)
            cli_error("If seq-json is not set seq-fasta and seq-mapping must be")
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
        isnothing(p) && cli_error(
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
                cli_error("If genome-json is set, genome-directories must not be set")
            json = expanduser(json)
            check_file(json, "genome-json")
            new((; json))
        else
            @debug "Validating genome args from directories"
            isnothing(directories) &&
                cli_error("If genome-json is not set, genome-directories must be set")
            v = parse_genomes_dir(directories)
            seen_flagsets = Set{FlagSet}()
            for (flagset, path) in v
                flagset ∈ seen_flagsets &&
                    cli_error("Flagset passed twice: $(join(flagset, '+'))")
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
            cli_error("Exactly one of tax-json or tax must be set")
        elseif !isnothing(json)
            isnothing(ncbi) || cli_error("If tax-json is set, tax-ncbi must not be set")
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
    makeref(outdir; kwargs...)

Implementation of the `makeref` subcommand. See `MakerefArgs` for the meaning of
the arguments.
"""
function makeref(
        outdir::String;
        seq_mapping::Union{String, Nothing} = nothing,
        seq_fasta::Union{String, Nothing} = nothing,
        seq_json::Union{String, Nothing} = nothing,
        tax::Union{String, Nothing} = nothing,
        tax_ncbi::Union{String, Nothing} = nothing,
        tax_json::Union{String, Nothing} = nothing,
        genome_directories::Union{String, Nothing} = nothing,
        genome_json::Union{String, Nothing} = nothing,
        overwrite::Bool = false,
        quiet::Bool = false,
    )
    set_global_logger!(; quiet)
    @debug "Validating sequence args"
    seq_args = SeqArgs(; json = seq_json, seq_mapping = seq_mapping, fasta = seq_fasta)
    @debug "Validating genome args"
    genome_args = GenomeArgs(; json = genome_json, directories = genome_directories)
    @debug "Validating taxonomy args"
    tax_args = TaxArgs(; json = tax_json, ncbi = tax_ncbi, tax = tax)
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
        open(io -> JSON.json(io, gns), joinpath(outdir, "genomes.json"), "w")

    taxes = taxonomy(tax_args, gns)
    @debug "Writing taxonomy.json file"
    isnothing(tax) ||
        open(io -> JSON.json(io, taxes), joinpath(outdir, "taxonomy.json"), "w")

    seqs = sequences(seq_args, gns)
    @debug "Writing seqs.json file"
    isnothing(seq_json) &&
        open(io -> JSON.json(io, seqs), joinpath(outdir, "seqs.json"), "w")

    @info "Constructing Reference object"
    @debug "Constructing ReferenceJSON"
    json_struct = BBB.ReferenceJSON(BBB.JSON_VERSION, gns, seqs, taxes)
    @debug "Constructing Reference from ReferenceJSON"
    ref = rethrow_as_cli_error(ReferenceError, "reference") do
        Reference(json_struct)
    end
    @debug "Writing reference object"
    open(joinpath(outdir, "reference.json"), "w") do io
        BBB.save(io, ref)
    end
    @info "Done"
end
