export pcacomp
export biplot
export screeplot
export npca

"""
    pcacomp(m; <keyword arguments>)

Calculate `n` first Primary Components (PCs).

# Arguments

- `m::Matrix{Float64}`: observations × variables
- `n=size(m, 2)::Int64`: number of PCs
- `zstd::Bool=true`: perform z score standardization before performing PCA

# Returns

Named tuple containing:
- `pc::DataFrame`: PC(1)..PC(n)
- `pcv::Vector{Float64}`: PC variances (fraction of total variance explained)
- `pcm::Vector{Float64}`: PC means
- `pcp::Matrix{Float64}`: PC projections
- `pc_model::MultivariateStats.PCA{Float64}`: PC model
"""
function pcacomp(m::Matrix{Float64}; n::Int64=size(m, 2), zstd::Bool=true)::@NamedTuple{pc::DataFrame, pcv::Vector{Float64}, pcm::Vector{Float64}, pcp::Matrix{Float64}, pc_model::MultivariateStats.PCA{Float64}}

    @assert n >= 1 "n must be ≥ 1."
    @assert n <= size(m, 2) "n must be ≤ $(size(m, 2))."

    if zstd
        for idx in axes(m, 2)
            if length(unique(m[:, idx])) > 1
                m[:, idx] = StatsBase.standardize(ZScoreTransform, m[:, idx])
            else
                m[:, idx] = zeros(size(m, 1))
            end
        end
    end

    # check maximum n
    pc_tmp = []
    n_tmp = n
    if zstd
        pc_tmp = @views MultivariateStats.fit(PCA, Matrix(m'), maxoutdim=n, pratio=1, mean=0)
    else
        pc_tmp = @views MultivariateStats.fit(PCA, Matrix(m'), maxoutdim=n, pratio=1)
    end
    size(pc_tmp.proj, 2) < n_tmp && (n_tmp = size(pc_tmp.proj, 2))
    n_tmp < n && _warn("Only $n_tmp PCs were generated.")
    n = n_tmp

    if zstd
        pc_model = MultivariateStats.fit(PCA, Matrix(m'), maxoutdim=n, pratio=1, mean=0)
    else
        pc_model = MultivariateStats.fit(PCA, Matrix(m'), maxoutdim=n, pratio=1)
    end
    pcv = MultivariateStats.principalvars(pc_model) ./ MultivariateStats.var(pc_model) * 100

    pc = DataFrame(Matrix(MultivariateStats.predict(pc_model, Matrix(m'))'), :auto)
    DataFrames.rename!(pc, ["PC" * string(x) for x in 1:n])

    return (pc=pc, pcv=pcv, pcm=pc_model.mean, pcp=pc_model.proj, pc_model=pc_model)

end

"""
    pcacomp(df, vars; <keyword arguments>)

Calculate `n` first Primary Components (PCs).

# Arguments

- `df::DataFrame`
- `vars::Union{Vector{String}, Vector{Symbol}}`: variable names
- `n::Int64`: number of PCs
- `zstd::Bool=true`: perform z score standardization before performing PCA

# Returns

Named tuple containing:
- `pc::DataFrame`: PC(1)..PC(n)
- `pcv::Vector{Float64}`: PC variances (fraction of total variance explained)
- `pcm::Vector{Float64}`: PC means
- `pcp::Matrix{Float64}`: PC projections
- `pc_model::MultivariateStats.PCA{Float64}`: PC model
"""
function pcacomp(df::DataFrame, vars::Union{Vector{String}, Vector{Symbol}}; n::Int64=length(vars), zstd::Bool=true)::@NamedTuple{pc::DataFrame, pcv::Vector{Float64}, pcm::Vector{Float64}, pcp::Matrix{Float64}, pc_model::MultivariateStats.PCA{Float64}}

    @assert length(vars) > 1 "vars must contain at least 2 variable names."
    @assert n >= 1 "n must be ≥ 1."
    @assert n <= length(vars) "n must be ≤ $(length(vars))."

    for idx in vars
        @assert idx in names(df) "$idx not found in df column names."
    end

    pca = pcacomp(Float64.(Matrix(df[!, vars])), n=n, zstd=zstd)

    return pca

end

"""
    biplot(df, vars; <keyword arguments>)

Plot PCA biplot.

# Arguments

- `df::DataFrame`
- `vars::Union{Vector{String}, Vector{Symbol}}`: variable names
- `n::Int64`: number of PCs
- `zstd::Bool=true`: perform z score standardization before performing PCA

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function biplot(df::DataFrame, vars::Union{Vector{String}, Vector{Symbol}}; n::Int64=length(vars), zstd::Bool=true, kwargs...)::Plots.Plot{Plots.GRBackend}

    pca = pcacomp(df, vars, n=n, zstd=zstd)
    n = length(pca.pc_model.prinvars)
    if n >= 2
        p = Plots.scatter(pca.pc[:, "PC1"],
                          pca.pc[:, "PC2"],
                          title="Biplot",
                          xlabel="PC1 ($(round(pca.pcv[1], digits=1))%)",
                          ylabel="PC2 ($(round(pca.pcv[2], digits=1))%)",
                          xlims=(-4, 4),
                          ylims=(-4, 4),
                          aspect_ratio=1,
                          palette=:darktest,
                          markersize=2,
                          markerstrokewidth=0,
                          markerstrokealpha=0,
                          framestyle=:box,
                          label=false;
                          kwargs...)
        for i in 1:n
            p = Plots.plot!([0, pca.pcp[i, 1] * 2],
                            [0, pca.pcp[i, 2] * 2],
                            arrow=true,
                            label=vars[i],
                            legend = :outerright)
        end

        Plots.plot(p; kwargs...)

        return p
    else
        _warn("Could not calculate ≥ 2 PCs.")
    end

end

"""
    screeplot(df, vars; <keyword arguments>)

Plot PCA scree plot.

# Arguments

- `df::DataFrame`
- `vars::Union{Vector{String}, Vector{Symbol}}`: variable names
- `n::Int64`: number of PCs
- `zstd::Bool=true`: perform z score standardization before performing PCA

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function screeplot(df::DataFrame, vars::Union{Vector{String}, Vector{Symbol}}; n::Int64=length(vars), zstd::Bool=true, kwargs...)::Plots.Plot{Plots.GRBackend}

    pca = pcacomp(df, vars, n=n, zstd=zstd)
    n = length(pca.pc_model.prinvars)

    xl = ["PC" * string(x) for x in 1:n]
    p1 = Plots.bar(1:n,
                   pca.pcv,
                   title="",
                   xticks=(1:n, xl),
                   ylabel="% of variance\nexplained",
                   ylims=(0, 100),
                   palette=:darktest,
                   framestyle=:box)
    p2 = Plots.plot(1:n,
                    pca.pc_model.prinvars,
                    ylims=(0, ceil(maximum(pca.pc_model.prinvars), digits=0)),
                    xticks=(1:n, xl),
                    ylabel="Eigenvalues",
                    palette=:darktest,
                    framestyle=:box)
    p2 = Plots.scatter!(1:n,
                        pca.pc_model.prinvars,
                        mc=1,
                        markerstrokewidth=0,
                        markerstrokealpha=0)
    p2 = hline!([1],
                ls=:dash,
                lc=:black)
    p = Plots.plot(p1, p2,
                   legend=false,
                   layout=(2, 1),
                   plot_title="Scree plot")
    Plots.plot(p; kwargs...)

    return p

end

"""
    npca(m; <keyword arguments>)

Calculate recommended number of Primary Components (PCs).

# Arguments

- `m::Matrix{Float64}`: observations × variables
- `zstd::Bool=true`: perform z score standardization before performing PCA
- `type::Symbol`
    - `:var`: use total variation
    - `:eig`: use eigenvalue
- `value::Real`: minimum % of total variation or threshold for eigenvalues (keep the components with eigenvalues greater than the threshold)

# Returns

- `n::Int64`: recommended number of PCs
"""
function npca(m::Matrix{Float64}; zstd::Bool=true, type::Symbol, value::Real)::Int64

    _check_var(type, [:var, :eig], "type")
    if type === :var
        _in(value, (0, 1), "value")
    else
        @assert value > 0 "For :eig type, value must be > 0."
    end

    if zstd
        for idx in axes(m, 2)
            if length(unique(m[:, idx])) > 1
                m[:, idx] = StatsBase.standardize(ZScoreTransform, m[:, idx])
            else
                m[:, idx] = zeros(size(m, 1))
            end
        end
    end

    n = size(m, 2)
    if zstd
        pc_model = MultivariateStats.fit(PCA, Matrix(m'), maxoutdim=n, pratio=1, mean=0)
    else
        pc_model = MultivariateStats.fit(PCA, Matrix(m'), maxoutdim=n, pratio=1)
    end
    if type === :var
        pcv = MultivariateStats.principalvars(pc_model) ./ MultivariateStats.var(pc_model)
        pcv = cumsum(pcv)
        if value > pcv[end]
            return length(pcv)
        else
            for idx in eachindex(pcv)
                value <= pcv[idx] && return idx
            end
        end
    else
        eig = pc_model.prinvars
        return count(x->x > value, eig)
    end

end
