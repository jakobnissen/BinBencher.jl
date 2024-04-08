struct OutputSettings
    recalls::Vector{Float64}
    precisions::Vector{Float64}
end

function OutputSettings(
    recalls::Union{Vector{Float64}, Nothing},
    precisions::Union{Vector{Float64}, Nothing},
)
    OutputSettings(
        isnothing(recalls) ? DEFAULT_RECALLS : recalls,
        isnothing(precisions) ? DEFAULT_PRECISIONS : precisions,
    )
end

function populate_output(
    dir::AbstractString,
    binning::BBB.Binning,
    settings::OutputSettings,
)
    open(joinpath(dir, "stats.json"), "w") do json
        output_stats(json, binning)
    end
    open(joinpath(dir, "recovery.json"), "w") do json
        output_recovery(json, binning, settings)
    end
    open(GzipCompressorStream, joinpath(dir, "bins.json.gz"), "w") do json
        d = Dict{String, Any}()
        for bin in binning.bins
            d[bin.name] = json_bin(bin)
        end
        JSON3.write(json, d)
    end
end

function output_stats(io::IO, binning::BBB.Binning)
    b = binning
    d = Dict{String, Any}("n_bins" => BBB.n_bins(b))
    d["n_sequences"] = sum(BBB.n_seqs, binning.bins; init=0)
    add_mean_bin_stats!(d, b)
    JSON3.pretty(io, d)
    write(io, '\n')
end

function add_mean_bin_stats!(d::Dict{String, Any}, b::BBB.Binning, precision=5) # TODO: Collapse genomic and asm precision - it's only 1 thing
    for (stats, name) in [(b.bin_genome_stats, "genomic"), (b.bin_asm_stats, "asm")]
        for (metric, field) in [
            ("recall", :mean_bin_recall),
            ("precision", :mean_bin_precision),
            ("F1", :mean_bin_f1),
        ]
            d["mean_bin_$(name)_$(metric)"] =
                round(getproperty(stats, field); digits=precision)
        end
    end
end

function output_recovery(io::IO, binning::Binning, settings::OutputSettings)
    json = Dict{String, Any}(
        "recalls" => settings.recalls,
        "precisions" => settings.precisions,
    )
    for (name, Ms) in [
        ("bins_asm_recall", binning.bin_asms),
        ("bins_genomic_recall", binning.bin_genomes),
        ("genomes_asm_recall", binning.recovered_asms),
        ("genomes_genomic_recall", binning.recovered_genomes),
    ]
        json[name] = []
        for M in Ms
            push!(json[name], [])
            for row in eachrow(M)
                push!(last(json[name]), row)
            end
        end
    end
    JSON3.write(io, json)
end

function add_passing_bins!(
    d::Dict{String, Any},
    b::BBB.Binning,
    recall::AbstractFloat,
    precision::AbstractFloat,
)
    for rank in 0:(BBB.nranks(b.ref) - 1)
        rank_suffix = iszero(rank) ? "" : "_rank_" * string(rank)

        # TODO: Proper error if this is not an integer - see the functions
        d["recovered_genomes" * rank_suffix] =
            BBB.n_recovered(b, recall, precision; level=rank, assembly=false)::Integer
        d["recovered_asms" * rank_suffix] =
            BBB.n_recovered(b, recall, precision; level=rank, assembly=true)::Integer
        d["passing_bins_genomes" * rank_suffix] =
            BBB.n_passing_bins(b, recall, precision; level=rank, assembly=false)::Integer
        d["passing_bins_asms" * rank_suffix] =
            BBB.n_passing_bins(b, recall, precision; level=rank, assembly=true)::Integer
    end
end

function json_bin(bin::BBB.Bin)::Dict{String, Any}
    result = Dict{String, Any}()
    d = Dict{String, Any}()
    result["genomes"] = d
    for (; genome, recall, precision) in BBB.recalls_precisions(bin; assembly=false) # TODO: Can we do this in one loop?
        d[genome.name] = Dict(
            "genomic_recall" => round(recall; digits=5),
            "precision" => round(precision; digits=5),
        )
    end
    for (; genome, recall) in BBB.recalls_precisions(bin; assembly=true)
        d[genome.name]["asm_recall"] = round(recall; digits=5)
    end
    d = Dict{String, Any}()
    result["clades"] = d
    for (; clade, recall, precision) in
        BBB.recalls_precisions(BBB.Clade{BBB.Genome}, bin; assembly=false) # TODO: Also merge two loops here. Perhaps change API
        d[clade.name] = Dict(
            "genomic_recall" => round(recall; digits=5),
            "precision" => round(precision; digits=5),
        )
    end
    for (; clade, recall) in
        BBB.recalls_precisions(BBB.Clade{BBB.Genome}, bin; assembly=true)
        d[clade.name]["asm_recall"] = round(recall; digits=5)
    end
    result
end
