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
function pcacomp(
    m::Matrix{Float64}; n::Int64 = size(m, 2), zstd::Bool = true
)::@NamedTuple{
    pc::DataFrame,
    pcv::Vector{Float64},
    pcm::Vector{Float64},
    pcp::Matrix{Float64},
    pc_model::MultivariateStats.PCA{Float64},
}

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
        pc_tmp = @views MultivariateStats.fit(PCA, Matrix(m'), maxoutdim = n, pratio = 1, mean = 0)
    else
        pc_tmp = @views MultivariateStats.fit(PCA, Matrix(m'), maxoutdim = n, pratio = 1)
    end
    size(pc_tmp.proj, 2) < n_tmp && (n_tmp = size(pc_tmp.proj, 2))
    n_tmp < n && _warn("Only $n_tmp PCs were generated.")
    n = n_tmp

    if zstd
        pc_model = MultivariateStats.fit(PCA, Matrix(m'); maxoutdim = n, pratio = 1, mean = 0)
    else
        pc_model = MultivariateStats.fit(PCA, Matrix(m'); maxoutdim = n, pratio = 1)
    end
    pcv = MultivariateStats.principalvars(pc_model) ./ MultivariateStats.var(pc_model) * 100

    pc = DataFrame(Matrix(MultivariateStats.predict(pc_model, Matrix(m'))'), :auto)
    DataFrames.rename!(pc, ["PC" * string(x) for x in 1:n])

    return (pc = pc, pcv = pcv, pcm = pc_model.mean, pcp = pc_model.proj, pc_model = pc_model)

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
function pcacomp(
    df::DataFrame, vars::Union{Vector{String}, Vector{Symbol}}; n::Int64 = length(vars), zstd::Bool = true
)::@NamedTuple{
    pc::DataFrame,
    pcv::Vector{Float64},
    pcm::Vector{Float64},
    pcp::Matrix{Float64},
    pc_model::MultivariateStats.PCA{Float64},
}

    @assert length(vars) > 1 "vars must contain at least 2 variable names."
    @assert n >= 1 "n must be ≥ 1."
    @assert n <= length(vars) "n must be ≤ $(length(vars))."

    for idx in vars
        @assert idx in names(df) "$idx not found in df column names."
    end

    pca = pcacomp(Float64.(Matrix(df[!, vars])); n = n, zstd = zstd)

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

  - `Union{Nothing, GLMakie.Figure}`
"""
function biplot(
    df::DataFrame, vars::Union{Vector{String}, Vector{Symbol}}; n::Int64 = length(vars), zstd::Bool = true
)::Union{Nothing, GLMakie.Figure}

    pca = pcacomp(df, vars; n = n, zstd = zstd)
    n = length(pca.pc_model.prinvars)
    if n >= 2
        p = GLMakie.Figure()
        ax = GLMakie.Axis(
            p[1, 1];
            aspect = 1,
            title = "Biplot",
            xlabel = "PC1 ($(round(pca.pcv[1], digits=1))%)",
            ylabel = "PC2 ($(round(pca.pcv[2], digits=1))%)",
        )
        GLMakie.xlims!(ax, (-4, 4))
        GLMakie.ylims!(ax, (-4, 4))
        GLMakie.scatter!(ax, pca.pc[:, "PC1"], pca.pc[:, "PC2"]; markersize = 10, colormap = :darktest)
        cmap = GLMakie.resample_cmap(:darktest, n)
        for idx in 1:n
            GLMakie.arrows2d!(
                (0, 0),
                (pca.pcp[idx, 1] * 2, pca.pcp[idx, 2] * 2);
                color = cmap[idx],
                colormap = :darktest,
                colorrange = 1:n,
                label = vars[idx],
            )
        end
        axislegend(; position = :rt, colormap = :darktest)
        return p
    else
        _warn("Could not calculate ≥ 2 PCs.")
        return nothing
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

  - `p::GLMakie.Figure`
"""
function screeplot(
    df::DataFrame, vars::Union{Vector{String}, Vector{Symbol}}; n::Int64 = length(vars), zstd::Bool = true
)::GLMakie.Figure

    pca = pcacomp(df, vars; n = n, zstd = zstd)
    n = length(pca.pc_model.prinvars)

    xl = ["PC" * string(x) for x in 1:n]
    p = GLMakie.Figure()
    ax1 = GLMakie.Axis(
        p[1, 1];
        title = "Scree plot",
        xticks = (1:n, xl),
        xlabel = "",
        ylabel = "% of variance explained",
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
    )
    GLMakie.xlims!(ax1, (0.5, n + 0.5))
    GLMakie.ylims!(ax1, (0, 100))
    cmap = GLMakie.resample_cmap(:darktest, n)
    for idx in 1:n
        GLMakie.barplot!(ax1, idx, pca.pcv[idx]; color = cmap[idx], colorrange = 1:n, colormap = :darktest)
    end
    ax2 = GLMakie.Axis(p[2, 1]; title = "", xticks = (1:n, xl), xlabel = "", ylabel = "Eigenvalues")
    GLMakie.xlims!(ax2, (0.5, n + 0.5))
    GLMakie.ylims!(ax2, (0, ceil(maximum(pca.pc_model.prinvars), digits = 0)))
    GLMakie.lines!(ax2, 1:n, pca.pc_model.prinvars; color = :black)
    GLMakie.scatter!(ax2, 1:n, pca.pc_model.prinvars; markersize = 10, color = :black)
    GLMakie.hlines!(ax2, 1; linestyle = :dash, color = :black)

    return p

end

"""
    npca(m; <keyword arguments>)

Calculate recommended number of Primary Components (PCs).

# Arguments

  - `m::Matrix{Float64}`: observations × variables
  - `zstd::Bool=true`: perform z score standardization before performing PCA
  - `type::Symbol`
      + `:var`: use total variation
      + `:eig`: use eigenvalue
  - `value::Real`: minimum % of total variation or threshold for eigenvalues (keep the components with eigenvalues greater than the threshold)

# Returns

  - `n::Int64`: recommended number of PCs
"""
function npca(m::Matrix{Float64}; zstd::Bool = true, type::Symbol, value::Real)::Int64

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
        pc_model = MultivariateStats.fit(PCA, Matrix(m'); maxoutdim = n, pratio = 1, mean = 0)
    else
        pc_model = MultivariateStats.fit(PCA, Matrix(m'); maxoutdim = n, pratio = 1)
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
