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
using Dates: Dates
using .Threads: @threads, @spawn, nthreads
using CRC32c: crc32c

const VERSION = let
    toml_path = joinpath(dirname(dirname(@__FILE__)), "Project.toml")
    vline = open(toml_path) do io
        first(Iterators.filter(startswith("version = "), eachline(io)))
    end
    VersionNumber(strip(split(vline, " = ")[end], ['"']))
end

const BBB = BinBencherBackend

# Hour, minutes, seconds, miliseconds, in 24H format.
const TIME_FORMAT = Dates.DateFormat("HH:mm:SS:sss")

function start_info()
    @info "Running BinBencher version v$(VERSION)"
    @info "Date is $(Dates.today())"
    return @debug "Number of threads: $(nthreads())"
end

function is_binbencher_submodule(mod::Module)
    while mod !== BinBencher
        parent = parentmodule(mod)
        parent === mod && return false
    end
    return true
end
function metafmt(level::Logging.LogLevel, _module, group, id, file, line)
    if level == Logging.Debug
        return (:blue, "DEBUG", "")
    elseif level == Logging.Info
        return (:green, "INFO ", "")
    elseif level == Logging.Warn
        return (:yellow, "WARN ", "")
    elseif level == Logging.Error
        return (:red, "ERROR", "")
    else
        return (:default, "     ", "")
    end
end

function set_global_logger!(paths::Vararg{String}; quiet::Bool = false)
    stderr_loglevel = quiet ? Logging.Error : Logging.Info
    stderr_sink = Logging.ConsoleLogger(stderr, stderr_loglevel; meta_formatter=metafmt)
    stderr_sink = LoggingExtras.MinLevelLogger(stderr_sink, stderr_loglevel)
    file_sinks = map(paths) do path
        LoggingExtras.FormatLogger(path) do io, args
            println(io, args.message)
        end
    end
    sinks = (stderr_sink, file_sinks...)
    logger = LoggingExtras.TeeLogger(sinks...)
    logger = LoggingExtras.TransformerLogger(logger) do log
        merge(
            log,
            (;message = "$(Dates.format(Dates.now(), TIME_FORMAT)) | $(log.message)"),
        )
    end
    logger = LoggingExtras.EarlyFilteredLogger(logger) do log
        is_binbencher_submodule(log._module)
    end
    logger = LoggingExtras.MinLevelLogger(logger, Logging.Info)
    global_logger(logger)
    return nothing
end

include("output.jl")
include("cli.jl")
include("parse.jl")

Comonicon.@main

end # module BinBencher
