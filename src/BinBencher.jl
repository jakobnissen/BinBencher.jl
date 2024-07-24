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
using CodecZlib: GzipDecompressorStream, GzipCompressorStream
using LoggingExtras: LoggingExtras
using Logging: Logging, global_logger
using StyledStrings: @styled_str
using Dates: Dates
using TOML: TOML
using .Threads: @threads, @spawn, nthreads

const VERSION = let
    toml_path = joinpath(dirname(dirname(@__FILE__)), "Project.toml")
    VersionNumber(TOML.parse(read(toml_path, String))["version"])
end

const BBB = BinBencherBackend

# Hour, minutes, seconds, miliseconds, in 24H format.
const TIME_FORMAT = Dates.DateFormat("HH:mm:SS:sss")

function start_info()
    @info "BinBencher version v$(VERSION)"
    @info "Date is $(Dates.today())"
    @debug "Number of threads: $(nthreads())"
end

function is_binbencher_submodule(mod::Module)
    while mod !== BinBencher
        parent = parentmodule(mod)
        parent === mod && return false
    end
    return true
end

# TODO: Basic information should be printed to the files once created.
function set_global_logger!(paths::Vararg{String})
    stdout_sink = (
        LoggingExtras.FormatLogger() do io, args
            println(io, args.message)
        end
    )
    file_sinks = map(paths) do path
        LoggingExtras.FormatLogger(path) do io, args
            println(io, args.message)
        end
    end
    sinks = (stdout_sink, file_sinks...)
    logger = LoggingExtras.TeeLogger(sinks...)
    logger = LoggingExtras.TransformerLogger(logger) do log
        prefix = if log.level == Logging.Debug
            styled"{blue,bold:DEBUG}"
        elseif log.level == Logging.Info
            styled" {green,bold:INFO}"
        elseif log.level == Logging.Warn
            styled" {yellow,bold:WARN}"
        elseif log.level == Logging.Error
            styled"{red,bold:ERROR}"
        else
            "     "
        end
        merge(
            log,
            (;
                message=prefix *
                        " $(Dates.format(Dates.now(), TIME_FORMAT)) | $(log.message)"
            ),
        )
    end
    logger = LoggingExtras.EarlyFilteredLogger(logger) do log
        is_binbencher_submodule(log._module)
    end
    logger = LoggingExtras.MinLevelLogger(logger, Logging.Info)
    global_logger(logger)
    nothing
end

include("output.jl")
include("cli.jl")
include("parse.jl")

Comonicon.@main

end # module BinBencher
