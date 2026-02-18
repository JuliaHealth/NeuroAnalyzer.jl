export plot_coherence

"""
    plot_coherence(coh, f; <keyword arguments>)

Plot coherence.

# Arguments

  - `coh::Vector{Float64}`: coherence
  - `f::Vector{Float64}`: frequencies
  - `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit
  - `xlabel::String="Frequency [Hz]"`: x-axis label
  - `ylabel::String="Coherence"`: y-axis label
  - `title::String=""`: plot title
  - `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_coherence(
    coh::Vector{Float64},
    f::Vector{Float64};
    flim::Tuple{Real, Real} = (f[1], f[end]),
    xlabel::String = "Frequency [Hz]",
    ylabel::String = "Coherence",
    title::String = "",
    frq::Symbol = :lin,
    mono::Bool = false,
)::GLMakie.Figure

    @assert length(coh) == length(f) "Length of coherence vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, "flim")

    pal = mono ? :grays : :darktest

    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        flim = (sf[2], flim[2])
    end

    # prepare plot
    GLMakie.activate!(title = "plot_coherence()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        yticks = [0, 0.25, 0.5, 0.75, 1.0],
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xscale = frq === :lin ? identity : log,
        xautolimitmargin = (0, 0),
    )
    GLMakie.xlims!(ax, flim)
    GLMakie.ylims!(ax, -0.1, 1.1)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # draw coherence
    Makie.lines!(f, coh; linewidth = 2, color = :black)

    max_coh = maxat(coh, f)
    min_coh = minat(coh, f)
    _info("Minimum coherence $(round(coh[min_coh[2]], digits=3)) at $(round(min_coh[1], digits=2)) Hz")
    _info("Maximum coherence $(round(coh[max_coh[2]], digits=3)) at $(round(max_coh[1], digits=2)) Hz")

    return p

end

"""
    plot_coherence(coh, f; <keyword arguments>)

Plot multi-channel coherence.

# Arguments

  - `coh::Matrix{Float64}`: coherence
  - `f::Vector{Float64}`: frequencies
  - `clabels::Vector{String}=string.(1:size(coh, 1))`: channel pairs labels vector
  - `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit
  - `xlabel::String="Frequency [Hz]"`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  - `avg::Bool=false`: if true, plot averaged PSD
  - `ci95::Bool=false`: if true, plot mean and ±95% CI of averaged PSDs
  - `leg::Bool=true`: if true, add legend with channel labels
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_coherence(
    coh::Matrix{Float64},
    f::Vector{Float64};
    clabels::Vector{String} = string.(1:size(coh, 1)),
    flim::Tuple{Real, Real} = (f[1], f[end]),
    xlabel::String = "Frequency [Hz]",
    ylabel::String = "",
    title::String = "",
    frq::Symbol = :lin,
    avg::Bool = false,
    ci95::Bool = false,
    leg::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure

    ch_n = size(coh, 1)

    @assert size(coh, 2) == length(f) "Length of coherence vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, "flim")

    pal = mono ? :grays : :darktest

    # get mean and 95%CI
    if ci95
        coh_m, _, coh_u, coh_l = NeuroAnalyzer.msci95(coh)
    end

    # prepare plot
    GLMakie.activate!(title = "plot_coherence()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xscale = frq === :lin ? identity : log,
        xautolimitmargin = (0, 0),
    )
    GLMakie.xlims!(ax, flim)
    GLMakie.ylims!(ax, -0.1, 1.1)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    if ci95
        # draw upper 95% CI
        Makie.band!(f, coh_u, coh_l; alpha = 0.25, color = :grey, strokewidth = 0.5)

        # draw mean
        Makie.lines!(f, coh_m; color = :black, linewidth = 2)
    else
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            Makie.lines!(
                f,
                coh[idx, :];
                color = cmap[idx],
                colormap = pal,
                colorrange = 1:ch_n,
                linewidth = 2,
                label = clabels[idx],
            )
        end

        # draw averaged channels
        if avg
            coh_avg = mean(coh; dims = 1)[:]
            Makie.lines!(f, coh_avg; colormap = pal, linewidth = 4, color = :black)
        end

        (leg && ch_n < 30) && axislegend(; position = :rt, colormap = pal)

    end

    return p

end
