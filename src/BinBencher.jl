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
using Dates: now

const BBB = BinBencherBackend

function is_binbencher_submodule(mod::Module)
    while mod !== BinBencher
        parent = parentmodule(mod)
        parent === mod && return false
    end
    return true
end

function create_logger(logfile::String)
    stdout_logger = LoggingExtras.FormatLogger() do io, args
        println(io, args.message)
    end
    file_logger = LoggingExtras.FormatLogger(logfile) do io, args
        println(io, args.message)
    end
    logger = LoggingExtras.TeeLogger(stdout_logger, file_logger)
    logger = LoggingExtras.TransformerLogger(logger) do log
        prefix = if log.level == Logging.Debug
            styled"{blue,bold:DEBUG} "
        elseif log.level == Logging.Info
            styled" {green,bold:INFO} "
        elseif log.level == Logging.Warn
            styled" {yellow,bold:WARN} "
        elseif log.level == Logging.Error
            styled"{red,bold:ERROR} "
        else
            "      "
        end
        merge(log, (; message=prefix * " $(now()) | $(log.message)"))
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
