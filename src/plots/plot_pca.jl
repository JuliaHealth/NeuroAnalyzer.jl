export plot_pca
export plot_scree

"""
    plot_pca(df, vars; <keyword arguments>)

Plot a PCA biplot (PC1 vs PC2 scores with loading vectors).

Requires at least 2 PCs; returns `nothing` with a warning otherwise.

# Arguments

- `df::DataFrame`: input data
- `vars::Union{Vector{String}, Vector{Symbol}}`: variable names for PCA
- `n::Int64=length(vars)`: number of PCs to compute (≥ 2 needed for the plot)
- `zstd::Bool=true`: if `true`, Z-score standardise before PCA

# Returns

- `GLMakie.Figure`: biplot figure, or `nothing` if fewer than 2 PCs result
"""
function plot_pca(
    df::DataFrame,
    vars::Union{Vector{String}, Vector{Symbol}};
    n::Int64 = length(vars),
    zstd::Bool = true,
)::Union{Nothing, GLMakie.Figure}
    pca = pcacomp(df, vars; n = n, zstd = zstd)
    n_pc = length(pca.pc_model.prinvars)

    if n_pc < 2
        _warn("Could not compute ≥ 2 PCs; plot not generated.")
        return nothing
    end

    GLMakie.activate!(; title = "plot_pca()")
    fig = GLMakie.Figure()
    ax = GLMakie.Axis(
        fig[1, 1];
        aspect = 1,
        title = "Biplot",
        xlabel = "PC1 ($(round(pca.pcv[1], digits = 1))%)",
        ylabel = "PC2 ($(round(pca.pcv[2], digits = 1))%)",
    )
    GLMakie.scatter!(ax, pca.pc[:, "PC1"], pca.pc[:, "PC2"]; markersize = 10)

    # one loading vector per variable: (loading on PC1, loading on PC2);
    # pcp is (variables × PCs); scaled by 2 for visibility
    n_var = length(vars)
    cmap = GLMakie.resample_cmap(:darktest, max(2, n_var))
    for idx = 1:n_var
        GLMakie.arrows2d!(
            (0, 0),
            (pca.pcp[idx, 1] * 2, pca.pcp[idx, 2] * 2);
            color = cmap[idx],
            label = string(vars[idx]),
        )
    end
    axislegend(; position = :rt)

    return fig
end

"""
    plot_scree(df, vars; <keyword arguments>)

Plot a PCA scree plot showing variance explained and eigenvalues per PC.

# Arguments

- `df::DataFrame`: input data
- `vars::Union{Vector{String}, Vector{Symbol}}`: variable names for PCA
- `n::Int64=length(vars)`: number of PCs to compute
- `zstd::Bool=true`: if `true`, Z-score standardize before PCA

# Returns

- `GLMakie.Figure`: two-panel figure (% variance explained + eigenvalues)
"""
function plot_scree(
    df::DataFrame,
    vars::Union{Vector{String}, Vector{Symbol}};
    n::Int64 = length(vars),
    zstd::Bool = true,
)::GLMakie.Figure
    pca = pcacomp(df, vars; n = n, zstd = zstd)
    n_pc = length(pca.pc_model.prinvars)
    xl = ["PC$i" for i = 1:n_pc]

    GLMakie.activate!(; title = "plot_scree()")
    fig = GLMakie.Figure()
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
    for idx = 1:n_pc
        GLMakie.barplot!(ax1, idx, pca.pcv[idx]; color = cmap[idx])
    end

    ax2 = GLMakie.Axis(
        fig[2, 1];
        xticks = (1:n_pc, xl),
        ylabel = "Eigenvalues",
    )
    GLMakie.xlims!(ax2, (0.5, n_pc + 0.5))
    GLMakie.ylims!(ax2, (0, ceil(maximum(pca.pc_model.prinvars); digits = 0)))
    GLMakie.lines!(ax2, 1:n_pc, pca.pc_model.prinvars; color = :black)
    GLMakie.scatter!(ax2, 1:n_pc, pca.pc_model.prinvars; markersize = 10, color = :black)
    GLMakie.hlines!(ax2, 1; linestyle = :dash, color = :black)

    return fig
end
