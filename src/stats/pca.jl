export pcacomp
export biplot
export screeplot
export npca

# ---------------------------------------------------------------------------
# internal helper: Z-score standardize each column of a matrix in-place
# constant columns (zero variance) are replaced with zeros
# ---------------------------------------------------------------------------
function _zstd_columns!(m::Matrix{Float64})
    for col in axes(m, 2)
        if length(unique(m[:, col])) > 1
            m[:, col] = StatsBase.standardize(ZScoreTransform, m[:, col])
        else
            m[:, col] .= 0.0
        end
    end
end

# ---------------------------------------------------------------------------
# internal helper: fit a PCA model, dispatching on the zstd flag
# ---------------------------------------------------------------------------
function _fit_pca(m::Matrix{Float64}, n::Int64, zstd::Bool)
    kwargs = zstd ? (maxoutdim=n, pratio=1, mean=0) : (maxoutdim=n, pratio=1)
    return MultivariateStats.fit(PCA, Matrix(m'); kwargs...)
end

"""
    pcacomp(m; <keyword arguments>)

Calculate the first `n` principal components (PCs) of a data matrix.

# Arguments

- `m::Matrix{Float64}`: data matrix of shape `(observations × variables)`; must have ≥ 2 rows and ≥ 1 column
- `n::Int64=size(m, 2)`: number of PCs to compute; must satisfy `1 ≤ n ≤ size(m, 2)`
- `zstd::Bool=true`: if `true`, Z-score standardize each variable before PCA

# Returns

Named tuple:

- `pc::DataFrame`: PC scores; columns named `PC1 … PCn`
- `pcv::Vector{Float64}`: percentage of total variance explained by each PC
- `pcm::Vector{Float64}`: column means of the (possibly standardised) data
- `pcp::Matrix{Float64}`: PC projection (loading) matrix
- `pc_model::MultivariateStats.PCA{Float64}`: fitted PCA model object

# Throws

- `ArgumentError`: if `n < 1`, `n > size(m, 2)`, or `size(m, 1) < 2`

# See also

[`pcacomp(::DataFrame, ::Union{Vector{String}, Vector{Symbol}})`](@ref), [`npca`](@ref), [`biplot`](@ref), [`screeplot`](@ref)
"""
function pcacomp(
    m::Matrix{Float64};
    n::Int64 = size(m, 2),
    zstd::Bool = true
)::@NamedTuple{
    pc::DataFrame,
    pcv::Vector{Float64},
    pcm::Vector{Float64},
    pcp::Matrix{Float64},
    pc_model::MultivariateStats.PCA{Float64},
}

    @assert size(m, 1) >= 2 "m must have at least 2 observations (rows)."
    @assert n >= 1 "n must be ≥ 1."
    @assert n <= size(m, 2) "n must be ≤ $(size(m, 2))."

    # work on a copy so the caller's matrix is not modified in-place
    m = copy(m)
    zstd && _zstd_columns!(m)

    if zstd
        for idx in axes(m, 2)
            if length(unique(m[:, idx])) > 1
                m[:, idx] = StatsBase.standardize(ZScoreTransform, m[:, idx])
            else
                m[:, idx] = zeros(size(m, 1))
            end
        end
    end

    # dry run to discover the actual number of PCs the solver can produce
    pc_tmp  = _fit_pca(m, n, zstd)
    n_actual = size(pc_tmp.proj, 2)
    if n_actual < n
        _warn("Only $n_actual PCs could be computed (requested $n).")
        n = n_actual
    end

    # final fit with the confirmed n
    pc_model = _fit_pca(m, n, zstd)
    pcv = MultivariateStats.principalvars(pc_model) ./
          MultivariateStats.var(pc_model) .* 100

    scores = Matrix(MultivariateStats.predict(pc_model, Matrix(m'))')
    pc = DataFrame(scores, ["PC$i" for i in 1:n])

    return (pc=pc, pcv=pcv, pcm=pc_model.mean, pcp=pc_model.proj, pc_model=pc_model)

end

"""
    pcacomp(df, vars; <keyword arguments>)

Calculate the first `n` principal components from selected columns of a DataFrame.

# Arguments

- `df::DataFrame`: input data
- `vars::Union{Vector{String}, Vector{Symbol}}`: column names to use; must all exist in `df` and there must be at least 2
- `n::Int64=length(vars)`: number of PCs; must satisfy `1 ≤ n ≤ length(vars)`
- `zstd::Bool=true`: if `true`, Z-score standardize each variable before PCA

# Returns

- `pc::DataFrame`: PC scores; columns named `PC1 … PCn`
- `pcv::Vector{Float64}`: percentage of total variance explained by each PC
- `pcm::Vector{Float64}`: column means of the (possibly standardised) data
- `pcp::Matrix{Float64}`: PC projection (loading) matrix
- `pc_model::MultivariateStats.PCA{Float64}`: fitted PCA model object

# Throws

- `ArgumentError`: if `vars` has fewer than 2 entries, any name is missing from `df`, or `n` is out of range

# See also

[`pcacomp(::Matrix{Float64})`](@ref), [`biplot`](@ref), [`screeplot`](@ref)
"""
function pcacomp(
    df::DataFrame,
    vars::Union{Vector{String}, Vector{Symbol}};
    n::Int64 = length(vars),
    zstd::Bool = true
)::@NamedTuple{
    pc::DataFrame,
    pcv::Vector{Float64},
    pcm::Vector{Float64},
    pcp::Matrix{Float64},
    pc_model::MultivariateStats.PCA{Float64},
}

    @assert length(vars) >= 2 "vars must contain at least 2 variable names."
    @assert n >= 1 "n must be ≥ 1."
    @assert n <= length(vars) "n must be ≤ $(length(vars))."

    for v in vars
        @assert string(v) in names(df) "Variable '$v' not found in df."
    end

    return pcacomp(Float64.(Matrix(df[!, vars])); n=n, zstd=zstd)

end

"""
    biplot(df, vars; <keyword arguments>)

Plot a PCA biplot (PC1 vs PC2 scores with loading vectors).

Requires at least 2 PCs; returns `nothing` with a warning otherwise.

# Arguments

- `df::DataFrame`: input data
- `vars::Union{Vector{String}, Vector{Symbol}}`: variable names for PCA
- `n::Int64=length(vars)`: number of PCs to compute (≥ 2 needed for the plot)
- `zstd::Bool=true`: if `true`, Z-score standardise before PCA

# Returns

- `GLMakie.Figure`: biplot figure, or `nothing` if fewer than 2 PCs result

# See also

[`screeplot`](@ref), [`pcacomp`](@ref)
"""
function biplot(
    df::DataFrame,
    vars::Union{Vector{String}, Vector{Symbol}};
    n::Int64 = length(vars),
    zstd::Bool = true
)::Union{Nothing, GLMakie.Figure}

    pca  = pcacomp(df, vars; n=n, zstd=zstd)
    n_pc = length(pca.pc_model.prinvars)

    if n_pc < 2
        _warn("Could not compute ≥ 2 PCs; biplot not generated.")
        return nothing
    end

    fig  = GLMakie.Figure()
    ax = GLMakie.Axis(
        fig[1, 1];
        aspect = 1,
        title = "Biplot",
        xlabel = "PC1 ($(round(pca.pcv[1], digits=1))%)",
        ylabel = "PC2 ($(round(pca.pcv[2], digits=1))%)",
    )
    GLMakie.xlims!(ax, (-4, 4))
    GLMakie.ylims!(ax, (-4, 4))
    GLMakie.scatter!(ax, pca.pc[:, "PC1"], pca.pc[:, "PC2"]; markersize=10)

    cmap = GLMakie.resample_cmap(:darktest, n_pc)
    for idx in 1:n_pc
        # scale loading vector by 2 for visibility; pcp is (variables × PCs)
        GLMakie.arrows2d!(
            (0, 0),
            (pca.pcp[1, idx] * 2, pca.pcp[2, idx] * 2);
            color = cmap[idx],
            label = string(vars[idx]),
        )
    end
    axislegend(; position=:rt)

    return fig

end

"""
    screeplot(df, vars; <keyword arguments>)

Plot a PCA scree plot showing variance explained and eigenvalues per PC.

# Arguments

- `df::DataFrame`: input data
- `vars::Union{Vector{String}, Vector{Symbol}}`: variable names for PCA
- `n::Int64=length(vars)`: number of PCs to compute
- `zstd::Bool=true`: if `true`, Z-score standardize before PCA

# Returns

- `GLMakie.Figure`: two-panel figure (% variance explained + eigenvalues)

# See also

[`biplot`](@ref), [`pcacomp`](@ref), [`npca`](@ref)
"""
function screeplot(
    df::DataFrame,
    vars::Union{Vector{String}, Vector{Symbol}};
    n::Int64 = length(vars),
    zstd::Bool = true
)::GLMakie.Figure

    pca = pcacomp(df, vars; n=n, zstd=zstd)
    n_pc = length(pca.pc_model.prinvars)
    xl = ["PC$i" for i in 1:n_pc]

    fig   = GLMakie.Figure()
    ax1 = GLMakie.Axis(
        fig[1, 1];
        title = "Scree plot",
        xticks = (1:n_pc, xl),
        ylabel = "% variance explained",
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
    )
    GLMakie.xlims!(ax1, (0.5, n_pc + 0.5))
    GLMakie.ylims!(ax1, (0, 100))

    cmap = GLMakie.resample_cmap(:darktest, n_pc)
    for idx in 1:n_pc
        GLMakie.barplot!(ax1, idx, pca.pcv[idx]; color=cmap[idx])
    end

    ax2 = GLMakie.Axis(
        fig[2, 1];
        xticks = (1:n_pc, xl),
        ylabel = "Eigenvalues",
    )
    GLMakie.xlims!(ax2, (0.5, n_pc + 0.5))
    GLMakie.ylims!(ax2, (0, ceil(maximum(pca.pc_model.prinvars), digits=0)))
    GLMakie.lines!(ax2, 1:n_pc, pca.pc_model.prinvars; color=:black)
    GLMakie.scatter!(ax2, 1:n_pc, pca.pc_model.prinvars; markersize=10, color=:black)
    GLMakie.hlines!(ax2, 1; linestyle=:dash, color=:black)

    return fig

end

"""
    npca(m; <keyword arguments>)

Calculate the recommended number of principal components (PCs).

# Arguments

- `m::Matrix{Float64}`: data matrix `(observations × variables)`; must have ≥ 2 rows and ≥ 1 column
- `zstd::Bool=true`: if `true`, Z-score standardise each variable before PCA. 
- `type::Symbol`: selection criterion:
    - `:var`: keep enough PCs to explain at least `value` fraction of total variance
    - `:eig`: keep PCs whose eigenvalue exceeds `value`
- `value::Real`: threshold; must be in `(0, 1)` for `:var`, or `> 0` for `:eig`

# Returns

- `Int64`: recommended number of PCs (≥ 1)

# Throws

- `ArgumentError`: if `type` is invalid or `value` is out of range

# See also

[`pcacomp`](@ref), [`screeplot`](@ref)
"""
function npca(m::Matrix{Float64}; zstd::Bool = true, type::Symbol, value::Real)::Int64

    _check_var(type, [:var, :eig], "type")
    if type === :var
        _in(value, (0, 1), "value")
    else
        @assert value > 0 "For :eig, value must be > 0."
    end
    @assert size(m, 1) >= 2 "m must have at least 2 observations (rows)."

    # work on a copy to avoid mutating the caller's matrix
    m = copy(m)
    zstd && _zstd_columns!(m)

    pc_model = _fit_pca(m, size(m, 2), zstd)

    if type === :var

        pcv = cumsum(
            MultivariateStats.principalvars(pc_model) ./
            MultivariateStats.var(pc_model)
        )
        # return the index of the first cumulative variance that meets the threshold
        idx = findfirst(>=(value), pcv)
        # if no single PC meets the threshold, return all PCs
        return isnothing(idx) ? length(pcv) : idx

    elseif type === :eig

        return max(1, count(>(value), pc_model.prinvars))

    end

end
