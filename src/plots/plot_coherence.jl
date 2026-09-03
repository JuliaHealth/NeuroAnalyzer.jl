export plot_coherence

"""
    plot_coherence(coh, f; <keyword arguments>)

Plot single-channel coherence as a function of frequency.

# Arguments

- `coh::Vector{Float64}`: coherence values (must match length of `f`)
- `f::Vector{Float64}`: frequency values (Hz)
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plot
- `xlabel::String="Frequency [Hz]"`: x-axis label
- `ylabel::String="Coherence"`: y-axis label
- `title::String=""`: plot title
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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
    # validate
    length(coh) == length(f) ||
        throw(
            ArgumentError(
                "Length of coherence vector must equal length of frequencies vector.",
            ),
        )
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(f), "flim")

    # set color palette
    pal = mono ? :grays : :darktest

    # adjust for log scale if flim[1] == 0
    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        flim = (sf[2], flim[2])
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_coherence()")
    plot_size = (900, 450)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        yticks = [0, 0.25, 0.5, 0.75, 1.0],
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xscale = frq === :lin ? identity : log,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )

    # set axis limits
    GLMakie.xlims!(ax, flim)
    GLMakie.ylims!(ax, -0.1, 1.1)
    _style_axis!(ax)

    # plot coherence
    GLMakie.lines!(
        f,
        coh;
        linewidth = 2,
        color = :black,
    )

    # find and log min/max coherence
    max_coh = maxat(coh, f)
    min_coh = minat(coh, f)
    _info(
        "Minimum coherence $(round(coh[min_coh[2]], digits = 3)) at $(round(min_coh[1], digits = 2)) Hz",
    )
    _info(
        "Maximum coherence $(round(coh[max_coh[2]], digits = 3)) at $(round(max_coh[1], digits = 2)) Hz",
    )

    return fig
end

"""
    plot_coherence(coh, f; <keyword arguments>)

Plot multi-channel coherence as a function of frequency.

# Arguments

- `coh::Matrix{Float64}`: coherence matrix, shape (channels, frequencies)
- `f::Vector{Float64}`: frequency values (Hz)
- `clabels::Vector{String}=string.(1:size(coh, 1))`: channel pair labels
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plot
- `xlabel::String="Frequency [Hz]"`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `avg::Bool=false`: if `true`, plot averaged coherence
- `ci95::Bool=false`: if `true`, plot mean and ±95% CI of averaged coherence
- `leg::Bool=true`: if `true`, add legend with channel labels
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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
    # validate
    size(coh, 2) == length(f) ||
        throw(
            ArgumentError(
                "Length of coherence vector must equal length of frequencies vector.",
            ),
        )
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(f), "flim")

    # number of channels
    ch_n = size(coh, 1)

    # set color palette
    pal = mono ? :grays : :darktest

    # prepare plot
    GLMakie.activate!(; title = "plot_coherence()")
    plot_size = (900, 450)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xscale = frq === :lin ? identity : log,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )

    # set axis limits
    GLMakie.xlims!(ax, flim)
    GLMakie.ylims!(ax, -0.1, 1.1)
    _style_axis!(ax)

    if ci95

        # calculate mean and 95% CI
        msci95_data = NeuroAnalyzer.msci95(coh)
        coh_m = msci95_data.sm
        coh_u = msci95_data.ul
        coh_l = msci95_data.ll

        # plot confidence interval band
        GLMakie.band!(
            f,
            coh_u,
            coh_l;
            alpha = 0.25,
            color = :grey,
            strokewidth = 0.5,
        )

        # plot mean coherence
        GLMakie.lines!(
            f,
            coh_m;
            color = :black,
            linewidth = 2,
        )

    else

        # set color palette
        cmap = GLMakie.resample_cmap(pal, ch_n)

        # plot each channel's coherence
        for idx = 1:ch_n
            GLMakie.lines!(
                f,
                coh[idx, :];
                color = cmap[idx],
                colormap = pal,
                colorrange = 1:ch_n,
                linewidth = 2,
                label = clabels[idx],
            )
        end

        # plot averaged coherence if requested
        if avg
            coh_avg = mean(coh; dims = 1)[:]
            Makie.lines!(f, coh_avg; colormap = pal, linewidth = 4, color = :black)
        end

        # add legend if requested and not too many channels
        if leg && ch_n < 30
            axislegend(ax; position = :rt, colormap = mono ? :grays : :darktest)
        end
    end

    return fig
end
