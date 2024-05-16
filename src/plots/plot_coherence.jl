export plot_coherence
export plot_coherence_avg
export plot_coherence_butterfly

"""
    plot_coherence(coh, f; <keyword arguments>)

Plot coherence.

# Arguments

- `coh::Vector{Float64}`: coherence
- `f::Vector{Float64}`: frequencies
- `frq_lim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit for the X-axis
- `xlabel::String="Frequency [Hz]"`: x-axis label
- `ylabel::String="Coherence"`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_coherence(coh::Vector{Float64}, f::Vector{Float64}; frq_lim::Tuple{Real, Real}=(f[1], f[end]), xlabel::String="Frequency [Hz]", ylabel::String="Coherence", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    @assert length(coh) == length(f) "Length of coherence vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin], "ax")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    if ax === :linlin
        xt = collect(_ticks(frq_lim))
        xsc = :identity
        ysc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.001, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.001 Hz")
            f[1] == 0 && (f[1] = 0.001)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        end
        xsc = :log10
        ysc = :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   ylim=(-0.1, 1.1),
                   yticks=[0, 0.25, 0.5, 0.75, 1.0],
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6)

    # plot coherence
    p = Plots.plot!(f,
                    coh,
                    linewidth=1,
                    alpha=0.5,
                    xticks=xt,
                    xscale=xsc,
                    yscale=ysc;
                    kwargs...)

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
- `clabels::Vector{String}=[""]`: channel pairs labels vector
- `frq_lim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit for the X-axis
- `xlabel::String="Frequency [Hz]"`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_coherence(coh::Matrix{Float64}, f::Vector{Float64}; clabels::Vector{String}=[""], frq_lim::Tuple{Real, Real}=(f[1], f[end]), xlabel::String="Frequency [Hz]", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    ch_n = size(coh, 1)
    @assert size(coh, 2) == length(f) "Length of coherence vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin], "ax")
    _check_tuple(frq_lim, "frq_lim")

    # reverse so 1st channel is on top
    coh_tmp = deepcopy(coh)
    coh = @views reverse(coh[:, eachindex(f)], dims = 1)
    # also, reverse colors if palette is not mono
    if mono
        pal = :grays
        channel_color = repeat([:black], ch_n)
    else
        pal = :darktest
        channel_color = ch_n:-1:1
    end

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(coh, 1)))

    # normalize and shift so all channels are visible
    # each channel is between 0 and +1.0
    for idx in 1:ch_n
        # scale by 0.5 so maxima do not overlap
        coh[idx, :] = @views normalize_n(coh[idx, :]) .+ (idx - 1)
    end

    if ax === :linlin
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.001, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.001 Hz")
            f[1] == 0 && (f[1] = 0.001)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        end
        xsc = :log10
        ysc = :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 800),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6)

    # plot zero line
    p = Plots.hline!(collect((ch_n - 1):-1:0),
                     color=:grey,
                     lw=0.5,
                     label="")
    # plot 0.25 line
    p = Plots.hline!(collect((ch_n - 0.25):-1:0),
                     color=:grey,
                     ls=:dot,
                     lw=0.5,
                     label="")
    # plot 0.5 line
    p = Plots.hline!(collect((ch_n - 0.5):-1:0),
                     color=:grey,
                     ls=:dot,
                     lw=0.5,
                     label="")
    # plot 0.75 line
    p = Plots.hline!(collect((ch_n - 0.75):-1:0),
                     color=:grey,
                     ls=:dot,
                     lw=0.5,
                     label="")
    # plot 1.0 line
    p = Plots.hline!([ch_n],
                     color=:grey,
                     ls=:dot,
                     lw=0.5,
                     label="")

    # plot channels
    for idx in 1:ch_n
        p = @views Plots.plot!(f,
                               coh[idx, :],
                               linewidth=1,
                               alpha=0.5,
                               label="",
                               xticks=xt,
                               xscale=xsc,
                               color=channel_color[idx])
    end

    # plot labels
    p = Plots.plot!(yticks=((ch_n - 1):-1:0, clabels))
    for idx in 1:size(coh, 1)
        max_coh = maxat(coh_tmp[idx, :], f)
        min_coh = minat(coh_tmp[idx, :], f)
        if clabels == repeat([""], size(coh, 1))
            _info("Channel pair $idx minimum coherence $(round(coh_tmp[idx, min_coh[2]], digits=3)) at $(round(min_coh[1], digits=2)) Hz")
            _info("Channel pair $idx maximum coherence $(round(coh_tmp[idx, max_coh[2]], digits=3)) at $(round(max_coh[1], digits=2)) Hz")
        else
            _info("Channel pair $(clabels[idx]) minimum coherence $(round(coh_tmp[idx, min_coh[2]], digits=3)) at $(round(min_coh[1], digits=2)) Hz")
            _info("Channel pair $(clabels[idx]) maximum coherence $(round(coh_tmp[idx, max_coh[2]], digits=3)) at $(round(max_coh[1], digits=2)) Hz")
        end
    end

    return p

end

"""
    plot_coherence_avg(coh, f; <keyword arguments>)

Plot coherence mean and ±95% CI of averaged channels.

# Arguments

- `coh::Matrix{Float64}`: coherence
- `f::Vector{Float64}`: frequencies
- `clabels::Vector{String}=[""]`: channel pairs labels vector
- `frq_lim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit for the X-axis
- `xlabel::String="Frequency [Hz]"`: x-axis label
- `ylabel::String="Coherence"`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_coherence_avg(coh::Matrix{Float64}, f::Vector{Float64}; clabels::Vector{String}=[""], frq_lim::Tuple{Real, Real}=(f[1], f[end]), xlabel::String="Frequency [Hz]", ylabel::String="Coherence", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    @assert size(coh, 2) == length(f) "Length of coherence vector must equal length of frequencies vector."
    _check_var(ax,[:linlin, :loglin], "ax")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = msci95(coh)

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(coh, 1)))

    if ax === :linlin
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.001, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.001 Hz")
            f[1] == 0 && (f[1] = 0.001)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        end
        xsc = :log10
        ysc = :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   ylim=(-0.1, 1.1),
                   xticks=xt,
                   yticks=[0, 0.25, 0.5, 0.75, 1.0],
                   xscale=xsc,
                   yscale=ysc,
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot upper 95% CI
    p = Plots.plot!(f,
                    s_u,
                    fillrange=s_l,
                    fillalpha=0.35,
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot lower 95% CI
    p = Plots.plot!(f,
                    s_l,
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot mean
    p = Plots.plot!(f,
                    s_m,
                    label=false,
                    t=:line,
                    c=:black,
                    lw=0.5)

    for idx in 1:size(coh, 1)
        max_coh = maxat(coh[idx, :], f)
        min_coh = minat(coh[idx, :], f)
        if clabels == repeat([""], size(coh, 1))
            _info("Channel pair $idx minimum coherence $(round(coh[idx, min_coh[2]], digits=3)) at $(round(min_coh[1], digits=2)) Hz")
            _info("Channel pair $idx maximum coherence $(round(coh[idx, max_coh[2]], digits=3)) at $(round(max_coh[1], digits=2)) Hz")
        else
            _info("Channel pair $(clabels[idx]) minimum coherence $(round(coh[idx, min_coh[2]], digits=3)) at $(round(min_coh[1], digits=2)) Hz")
            _info("Channel pair $(clabels[idx]) maximum coherence $(round(coh[idx, max_coh[2]], digits=3)) at $(round(max_coh[1], digits=2)) Hz")
        end
    end

    return p

end

"""
    plot_coherence_butterfly(coh, f; <keyword arguments>)

Butterfly PSD plot.

# Arguments

- `coh::Array{Float64, 3}`: coherence
- `f::Vector{Float64}`: frequencies
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `frq_lim::Tuple{Real, Real}=(f[1], f[end]): frequency limit for the x-axis
- `xlabel::String="Frequency [Hz]"`: x-axis label
- `ylabel::String="Coherence"`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_coherence_butterfly(coh::Matrix{Float64}, f::Vector{Float64}; clabels::Vector{String}=[""], frq_lim::Tuple{Real, Real}=(f[1], f[end]), xlabel::String="Frequency [Hz]", ylabel::String="Coherence", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    @assert size(coh, 2) == length(f) "Length of coherence vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin], "ax")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(coh, 1)))

    if ax === :linlin
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.001, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.001 Hz")
            f[1] == 0 && (f[1] = 0.001)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        end
        xsc = :log10
        ysc = :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   ylim=(-0.1, 1.1),
                   xticks=xt,
                   yticks=[0.0, 0.25, 0.5, 0.75, 1.0],
                   xscale=xsc,
                   yscale=ysc;
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6)

    # plot coherence
    for idx in 1:size(coh, 1)
        p = Plots.plot!(f,
                        coh[idx, :],
                        t=:line,
                        linecolor=idx,
                        linewidth=0.5,
                        label=clabels[idx],
                        legend=true;
                        kwargs...)
    end

    for idx in 1:size(coh, 1)
        max_coh = maxat(coh[idx, :], f)
        min_coh = minat(coh[idx, :], f)
        if clabels == repeat([""], size(coh, 1))
            _info("Channel pair $idx minimum coherence $(round(coh[idx, min_coh[2]], digits=3)) at $(round(min_coh[1], digits=2)) Hz")
            _info("Channel pair $idx maximum coherence $(round(coh[idx, max_coh[2]], digits=3)) at $(round(max_coh[1], digits=2)) Hz")
        else
            _info("Channel pair $(clabels[idx]) minimum coherence $(round(coh[idx, min_coh[2]], digits=3)) at $(round(min_coh[1], digits=2)) Hz")
            _info("Channel pair $(clabels[idx]) maximum coherence $(round(coh[idx, max_coh[2]], digits=3)) at $(round(max_coh[1], digits=2)) Hz")
        end
    end

    return p

end
