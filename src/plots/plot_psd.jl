export plot_psd
export plot_psd_avg
export plot_psd_butterfly
export plot_psd_3d
export plot_psd_topo

"""
    plot_psd(sf, sp; <keyword arguments>)

Plot PSD (power spectrum density).

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Vector{Float64}`: powers
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd(sf::Vector{Float64}, sp::Vector{Float64}; norm::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    @assert length(sp) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    frq_lim = tuple_order(frq_lim)

    pal = mono == true ? :grays : :darktest

    if ax === :linlin
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = norm == false ? :log10 : :identity
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
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6)

    # plot powers
    p = Plots.plot!(sf,
                    sp,
                    xticks=xticks,
                    xscale=xscale,
                    yscale=yscale;
                    kwargs...)

    return p

end

"""
    plot_psd(sf, sp; <keyword arguments>)

Plot multi-channel PSD (power spectrum density).

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Matrix{Float64}`: powers
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd(sf::Vector{Float64}, sp::Matrix{Float64}; clabels::Vector{String}=[""], norm::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    ch_n = size(sp, 1)
    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    frq_lim = tuple_order(frq_lim)

    # reverse so 1st channel is on top
    sp = @views reverse(sp[:, eachindex(sf)], dims = 1)
    # also, reverse colors if palette is not mono
    if mono == true
        pal = :grays
        channel_color = repeat([:black], ch_n)
    else
        pal = :darktest
        channel_color = ch_n:-1:1
    end

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    # get range of the original sp for the scale
    range = _get_range(sp)

    # normalize and shift so all channels are visible
    # each channel is between -1.0 and +1.0
    for idx in 1:ch_n
        # scale by 0.5 so maxima do not overlap
        sp[idx, :] = @views normalize(sp[idx, :], method=:minmax) .* 0.5 .+ (idx - 1)
    end

    if ax === :linlin
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        _warn("For multi-channel PSD plots, y-axis log-scale is ignored.")
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = :identity
    elseif ax === :loglog
        _warn("For multi-channel PSD plots, y-axis log-scale is ignored.")
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
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

    # plot channels
    for idx in 1:ch_n
        p = @views Plots.plot!(sf,
                               sp[idx, :],
                               linewidth=1,
                               label="",
                               xticks=xticks,
                               xscale=xscale,
                               color=channel_color[idx])
    end

    # plot labels
    p = Plots.plot!(yticks=((ch_n - 1):-1:0, clabels))

    return p

end

"""
    plot_psd_avg(sf, sp; <keyword arguments>)

Plot PSD mean and ±95% CI of averaged channels.

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd_avg(sf::Vector{Float64}, sp::Array{Float64, 2}; norm::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax,[:linlin, :loglin, :linlog, :loglog], "ax")

    frq_lim = tuple_order(frq_lim)

    pal = mono == true ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = msci95(sp)

    if ax === :linlin
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = norm == false ? :log10 : :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   xticks=xticks,
                   xscale=xscale,
                   yscale=yscale,
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
    p = Plots.plot!(sf,
                    s_u,
                    fillrange=s_l,
                    fillalpha=0.35, 
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot lower 95% CI
    p = Plots.plot!(sf,
                    s_l,
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot mean
    p = Plots.plot!(sf,
                    s_m,
                    label=false,
                    t=:line,
                    c=:black,
                    lw=0.5)

    return p

end

"""
    plot_psd_butterfly(sf, sp; <keyword arguments>)

Butterfly PSD plot.

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd_butterfly(sf::Vector{Float64}, sp::Array{Float64, 2}; clabels::Vector{String}=[""], norm::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    frq_lim = tuple_order(frq_lim)

    pal = mono == true ? :grays : :darktest
    
    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    if ax === :linlin
        xticks=_ticks(frq_lim)
        xscale=:identity
        yscale=:identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = norm == false ? :log10 : :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   xticks=xticks,
                   xscale=xscale,
                   yscale=yscale;
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

    # plot powers
    for idx in 1:size(sp, 1)
        p = Plots.plot!(sf,
                        sp[idx, :],
                        t=:line,
                        linecolor=idx,
                        linewidth=0.5,
                        label=clabels[idx],
                        legend=true;
                        kwargs...)
    end

    return p

end

"""
    plot_psd_w3d(sf, sp; <keyword arguments>)

Plot 3-d waterfall PSD plot.

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `zlabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `variant::Symbol`: waterfall (`:w`) or surface (`:s`)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}`
"""
function plot_psd_3d(sf::Vector{Float64}, sp::Array{Float64, 2}; clabels::Vector{String}=[""], norm::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", zlabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, variant::Symbol, kwargs...)

    _check_var(variant, [:w, :s], "variant")
    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    frq_lim = tuple_order(frq_lim)

    ch_n = size(sp, 1)

    pal = mono == true ? :grays : :darktest
    
    # channel labels
    clabels == [""] && (clabels = repeat([""], ch_n))

    if ax === :linlin
        xticks=_ticks(frq_lim)
        xscale=:identity
        zscale=:identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        zscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        zscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        zscale = norm == false ? :log10 : :identity
    end

    # prepare plot
    if variant === :w
        p = Plots.plot(sf,
                       ones(length(sf)),
                       sp[1, :],
                       xlabel=xlabel,
                       ylabel=ylabel,
                       zlabel=zlabel,
                       legend=false,
                       xlims=frq_lim,
                       xticks=xticks,
                       xscale=xscale,
                       zscale=zscale;
                       title=title,
                       palette=pal,
                       st=:line,
                       lc=:black,
                       size=(1200, 800),
                       margins=20Plots.px,
                       titlefontsize=8,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6)

        # plot powers
        for idx in 2:ch_n
            p = Plots.plot!(sf,
                            ones(length(sf)) .* idx,
                            sp[idx, :],
                            st=:line,
                            linecolor=idx,
                            linewidth=0.5,
                            kwargs...)
        end
    else
        f1 = vsearch(frq_lim[1], sf)
        f2 = vsearch(frq_lim[2], sf)
        p = Plots.plot(sf[f1:f2],
                       eachindex(clabels),
                       sp[:, f1:f2],
                       xlabel=xlabel,
                       ylabel=ylabel,
                       zlabel=zlabel,
                       legend=false,
                       xlims=frq_lim,
                       xticks=xticks,
                       xscale=xscale,
                       zscale=zscale;
                       title=title,
                       palette=pal,
                       st=:surface,
                       lc=:black,
                       size=(1200, 800),
                       margins=20Plots.px,
                       titlefontsize=8,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6)
    end

    p = Plots.plot!(yticks=(1:ch_n, clabels))

    return p

end

"""
    plot_psd_topo(locs, sf, sp; <keyword arguments>)

Plot topographical map PSDs.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `ch::Union{Vector{Int64}, AbstractRange}`: which channels to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `cart::Bool=false`: if true, use Cartesian x and y coordinates, otherwise use polar radius and theta coordinates
- `kwargs`: optional arguments for plot() function

# Returns

- `fig::GLMakie.Figure`
"""
function plot_psd_topo(locs::DataFrame, sf::Vector{Float64}, sp::Array{Float64, 2}; ch=Union{Vector{Int64}, AbstractRange}, clabels::Vector{String}=[""], norm::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, cart::Bool=false, kwargs...)

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    @assert length(ch) <= nrow(locs) "Some channels do not have locations."

    frq_lim = tuple_order(frq_lim)

    pal = mono == true ? :grays : :darktest
    
    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    if ax === :linlin
        xticks=_ticks(frq_lim)
        xscale=:identity
        yscale=:identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
        end
        sf[1] == 0 && (sf[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = norm == false ? :log10 : :identity
    end

    # plot parameters
    plot_size = 1200
    marker_size = (150, 75)
    
    # get locations
    if cart == false
        loc_x = zeros(size(locs, 1))
        loc_y = zeros(size(locs, 1))
        for idx in 1:size(locs, 1)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
        loc_x = loc_x[ch]
        loc_y = loc_y[ch]
    else
        loc_x = locs[ch, :loc_x]
        loc_y = locs[ch, :loc_y]
    end
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)
    # get marker centers
    loc_x .*= ((plot_size / 2) - marker_size[1] / 2)
    loc_y .*= ((plot_size / 2) - marker_size[2] / 2)

    fig = Figure(; resolution=(plot_size, plot_size))
    fig_axis = Axis(fig[1, 1])
    fig_axis.aspect = AxisAspect(1)
    fig_axis.title = title
    GLMakie.xlims!(fig_axis, [-plot_size / 1.75, plot_size / 1.75])
    GLMakie.ylims!(fig_axis, [-plot_size / 1.75, plot_size / 1.75])
    hidedecorations!(fig_axis, grid=true, ticks=true)

    for idx in 1:size(sp, 1)
        p = Plots.plot(sf,
                       sp[idx, :],
                       t=:line,
                       c=:black,
                       linewidth=0.5,
                       xlabel=xlabel,
                       ylabel=ylabel,
                       legend=false,
                       xlims=frq_lim,
                       xticks=false,
                       yticks=false,
                       xscale=xscale,
                       yscale=yscale,
                       title=clabels[idx],
                       palette=pal,
                       size=marker_size,
                       #left_margin=20Plots.px,
                       titlefontsize=8,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6;
                       kwargs...)
        marker_img = tempname() * ".png"
        savefig(p, marker_img)
        marker = GLMakie.load(marker_img)
        GLMakie.scatter!(fig_axis, (loc_x[idx], loc_y[idx]), marker=marker, markersize=marker_size)
        rm(marker_img)
    end

    return fig

end

"""
    plot_psd(obj; <keyword arguments>)

Plot power spectrum density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:welch`: method of calculating PSD:
    - `:welch`: Welch's periodogram
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: x-axis limit
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher` 
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is Power [dB] or Power [units^2/Hz]
- `zlabel::String="default"`: z-axis label for 3-d plots, default is Power [dB] or Power [units^2/Hz]
- `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or grey palette
- `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}`
"""
function plot_psd(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, ch::Union{Int64, Vector{Int64}, <:AbstractRange}, norm::Bool=true, method::Symbol=:welch, nt::Int64=8, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, kwargs...)

    _check_var(type, [:normal, :butterfly, :mean, :w3d, :s3d, :topo], "type")
    _check_var(method, [:welch, :mt, :mw], "method")
    _check_var(ref, [:abs, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "ref")
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")
    ref !== :abs && @assert method !== :mw "For relative PSD, method must be :welch or :mt."

    _check_channels(obj, ch)

    @assert seg[1] != seg[2] "Signal is too short for analysis."

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    if ep != 0
        _check_epochs(obj, ep)
        if epoch_n(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    # get time vector
    if seg[2] <= epoch_len(obj)
        signal = obj.data[ch, seg[1]:seg[2], 1]
    else
        signal = epoch(obj, ep_n=1).data[ch, seg[1]:seg[2], 1]
    end
    # t = _get_t(seg[1], seg[2], sr(obj))
    t = obj.time_pts[seg[1]:seg[2]]
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    # set units
    units = _set_units(obj, ch[1])

    clabels = labels(obj)[ch]
    length(ch) == 1 && (clabels = [clabels])

    ref !== :abs && (f = band_frq(obj, band=ref))

    # get frequency range
    fs = sr(obj)
    frq_lim = tuple_order(frq_lim)
    @assert !(frq_lim[1] < 0 || frq_lim[2] > fs / 2) "frq_lim must be ≥ 0 and ≤ $(fs / 2)."

    # # get time vector
    # _, t_s1, _, t_s2 = _convert_t(obj.epoch_time[1], obj.epoch_time[end])

    if ref === :abs
        if method === :welch
            sp, sf = psd(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $ch, epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd(signal, fs=fs, norm=norm, mt=true, nt=nt)
            title == "default" && (title = "Absolute PSD (multi-tapered) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $ch, epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_mw(signal, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc)
            title == "default" && (title = "Absolute PSD (Morlet wavelet convolution) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $ch, epoch: $ep, time window: $t_s1:$t_s2]")
        end
    elseif ref === :total
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "PSD (Welch's periodogram) relative to total power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $ch, epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, norm=norm, mt=true)
            title == "default" && (title = "PSD (multi-tapered) relative to total power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $ch, epoch: $ep, time window: $t_s1:$t_s2]")
        end
    else
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, norm=norm, mt=false, f=f)
            title == "default" && (title = "PSD (Welch's periodogram) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $ch, epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, norm=norm, mt=true, f=f)
            title == "default" && (title = "PSD (multi-tapered) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $ch, epoch: $ep, time window: $t_s1:$t_s2]")
        end
    end

    # set labels
    if type !== :w3d && type !== :s3d && type !== :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if ref !== :abs
            ylabel == "default" && (ylabel = "Power ratio")
        end
        if norm == true
            ylabel == "default" && (ylabel = "Power [dB]")
        else
            ylabel == "default" && (ylabel = "Power [$units^2/Hz]")
        end
    end

    if type === :normal
        ch_t = obj.header.recording[:channel_type]
        if length(ch) > 1
            ch_t_uni = unique(ch_t[ch])
            @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels should be of the same type."
        end
        if ndims(sp) == 1
            p = plot_psd(sf,
                         sp,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         norm=norm,
                         frq_lim=frq_lim,
                         ax=ax,
                         mono=mono;
                         kwargs...)
        else
            p = plot_psd(sf,
                         sp,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         clabels=clabels,
                         title=title,
                         norm=norm,
                         frq_lim=frq_lim,
                         ax=ax,
                         mono=mono;
                         kwargs...)
        end
    elseif type === :butterfly
        ch_t = obj.header.recording[:channel_type]
        ch_t_uni = unique(ch_t[ch])
        @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels should be of the same type."
        @assert ndims(sp) >= 2 "For type=:butterfly plot the signal must contain ≥ 2 channels."
        title = replace(title, "channel" => "channels")
        p = plot_psd_butterfly(sf,
                               sp,
                               clabels=clabels,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               norm=norm,
                               frq_lim=frq_lim,
                               ax=ax,
                               mono=mono;
                               kwargs...)
    elseif type === :mean
        ch_t = obj.header.recording[:channel_type]
        ch_t_uni = unique(ch_t[ch])
        @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels should be of the same type."
        @assert ndims(sp) >= 2 "For type=:mean plot the signal must contain ≥ 2 channels."
        title = replace(title, "PSD" => "PSD [mean ± 95%CI]")
        title = replace(title, "channel" => "averaged channels")
        p = plot_psd_avg(sf,
                         sp,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         norm=norm,
                         frq_lim=frq_lim,
                         ax=ax,
                         mono=mono;
                         kwargs...)
    elseif type === :w3d
        ch_t = obj.header.recording[:channel_type]
        ch_t_uni = unique(ch_t[ch])
        @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels should be of the same type."
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        zlabel == "default" && (zlabel = norm == true ? "Power [dB]" : "Power [$units^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(sf,
                        sp,
                        clabels=clabels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        norm=norm,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:w;
                        kwargs...)
    elseif type === :s3d
        ch_t = obj.header.recording[:channel_type]
        ch_t_uni = unique(ch_t[ch])
        @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels should be of the same type."
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        zlabel == "default" && (zlabel = norm == true ? "Power [dB]" : "Power [$units^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(sf,
                        sp,
                        clabels=clabels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        norm=norm,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:s;
                        kwargs...)
    elseif type === :topo
        ch_t = obj.header.recording[:channel_type]
        ch_t_uni = unique(ch_t[ch])
        @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels should be of the same type."
        @assert _has_locs(obj) "Electrode locations not available."
        ndims(sp) == 1 && (sp = reshape(sp, 1, length(sp)))
        xlabel == "default" && (xlabel = "")
        ylabel == "default" && (ylabel = "")
        title = replace(title, "channel" => "channels")
        p = plot_psd_topo(obj.locs,
                          sf,
                          sp,
                          ch=ch,
                          clabels=clabels,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          norm=norm,
                          frq_lim=frq_lim,
                          ax=ax,
                          mono=mono;
                          kwargs...)
    end

    if p isa Plots.Plot{Plots.GRBackend}
        Plots.plot(p)
    else
        p
    end

    return p

end

"""
    plot_psd(obj; <keyword arguments>)

Plot power spectrum density of embedded or external component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:welch`: method of calculating PSD:
    - `:welch`: Welch's periodogram
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: x-axis limit
- `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher` 
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is Power [dB] or Power [units^2/Hz]
- `zlabel::String="default"`: z-axis label for 3-d plots, default is Power [dB] or Power [units^2/Hz]
- `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or grey palette
- `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
- `units::String=""`
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, norm::Bool=true, method::Symbol=:welch, nt::Int64=8, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, units::String="", kwargs...)

    _check_var(type, [:normal, :butterfly, :mean, :w3d, :s3d, :topo], "type")
    _check_var(method, [:welch, :mt, :mw], "method")
    _check_var(ref, [:abs, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "ref")
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")
    ref !== :abs && @assert method !== :mw "For relative PSD, method must be :welch or :mt."

    @assert seg[1] != seg[2] "Signal is too short for analysis."
    
    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    if ep != 0
        _check_epochs(obj, ep)
        if epoch_n(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    # select component c_idxs, default is all c_idxs
    c isa Symbol && (c = _get_component(obj, c))
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    clabels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (clabels = [clabels])

    ref !== :abs && (f = band_frq(obj, band=ref))

    # get frequency range
    fs = sr(obj)
    frq_lim = tuple_order(frq_lim)
    @assert !(frq_lim[1] < 0 || frq_lim[2] > fs / 2) "frq_lim must be ≥ 0 and ≤ $(fs / 2)."

    # get time vector
    if seg[2] <= epoch_len(obj)
        signal = c[c_idx, seg[1]:seg[2], 1]
    else
        signal = reshape(c, size(c, 1), :, 1)[c_idx, seg[1]:seg[2], 1]
    end
    # t = _get_t(seg[1], seg[2], sr(obj))
    t = obj.time_pts[seg[1]:seg[2]]
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    if ref === :abs
        if method === :welch
            sp, sf = psd(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd(signal, fs=fs, norm=norm, mt=true, nt=nt)
            title == "default" && (title = "Absolute PSD (multi-tapered) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_psd(signal, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc)
            title == "default" && (title = "Absolute PSD (Morlet wavelet convolution) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end
    elseif ref === :total
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "PSD (Welch's periodogram) relative to total power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, norm=norm, mt=true)
            title == "default" && (title = "PSD (multi-tapered) relative to total power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end
    else
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, norm=norm, mt=false, f=f)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, norm=norm, mt=true, f=f)
            title == "default" && (title = "Absolute PSD (multi-tapered) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end
    end

    # set labels
    if type !== :w3d && type !== :s3d
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if norm == true
            ylabel == "default" && (ylabel = "Power [dB]")
        else
            ylabel == "default" && (ylabel = "Power [$units^2/Hz]")
        end
    end

    if type === :normal
        @assert ndims(sp) == 1 "For type=:normal the signal must contain 1 c_idx."
        p = plot_psd(sf,
                     sp,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     title=title,
                     norm=norm,
                     frq_lim=frq_lim,
                     ax=ax,
                     mono=mono;
                     kwargs...)
    elseif type === :butterfly
        @assert ndims(sp) >= 2 "For type=:butterfly plot the signal must contain ≥ 2 c_idxs."
        title = replace(title, "component" => "components")
        p = plot_psd_butterfly(sf,
                               sp,
                               clabels=clabels,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               norm=norm,
                               frq_lim=frq_lim,
                               ax=ax,
                               mono=mono;
                               kwargs...)
    elseif type === :mean
        @assert ndims(sp) >= 2 "For type=:mean plot the signal must contain ≥ 2 c_idxs."
        title = replace(title, "PSD" => "PSD [mean ± 95%CI]")
        title = replace(title, "component" => "averaged components")
        p = plot_psd_avg(sf,
                         sp,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         norm=norm,
                         frq_lim=frq_lim,
                         ax=ax,
                         mono=mono;
                         kwargs...)
    elseif type === :w3d
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        zlabel == "default" && (zlabel = norm == true ? "Power [dB]" : "Power [$units^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(sf,
                        sp,
                        clabels=clabels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        norm=norm,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:w;
                        kwargs...)
    elseif type === :s3d
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        zlabel == "default" && (zlabel = norm == true ? "Power [dB]" : "Power [$units^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(sf,
                        sp,
                        clabels=clabels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        norm=norm,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:s;
                        kwargs...)
    end

    Plots.plot(p)

    return p
    
end
