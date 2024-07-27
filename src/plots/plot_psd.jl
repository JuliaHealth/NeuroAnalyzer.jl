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
- `db::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
    - `:linlog`: linear-log10
    - `:loglog`: log10-log10
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd(sf::Vector{Float64}, sp::Vector{Float64}; db::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    @assert length(sp) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    if ax === :linlin
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        ysc = :identity
    elseif ax === :linlog
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = !db ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        ysc = !db ? :log10 : :identity
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
                   titlefontsize=10,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6)

    # plot powers
    p = Plots.plot!(sf,
                    sp,
                    xticks=xt,
                    xscale=xsc,
                    yscale=ysc;
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
- `db::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
    - `:linlog`: linear-log10
    - `:loglog`: log10-log10
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd(sf::Vector{Float64}, sp::Matrix{Float64}; clabels::Vector{String}=[""], db::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    ch_n = size(sp, 1)
    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")
    _check_tuple(frq_lim, "frq_lim")

    # reverse so 1st channel is on top
    sp = @views reverse(sp[:, eachindex(sf)], dims = 1)
    # also, reverse colors if palette is not mono
    if mono
        pal = :grays
        channel_color = repeat([:black], ch_n)
    else
        pal = :darktest
        channel_color = ch_n:-1:1
    end

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    # normalize and shift so all channels are visible
    # each channel is between -1.0 and +1.0
    for idx in 1:ch_n
        # scale by 0.5 so maxima do not overlap
        sp[idx, :] = @views normalize(sp[idx, :], method=:minmax) .* 0.5 .+ (idx - 1)
    end

    if ax === :linlin
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        ysc = :identity
    elseif ax === :linlog
        _warn("For multi-channel PSD plots, y-axis log-scale is ignored.")
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = :identity
    elseif ax === :loglog
        _warn("For multi-channel PSD plots, y-axis log-scale is ignored.")
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
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
                   size=size(sp, 1) <= 64 ? (1200, 800) : (1200, 1200),
                   margins=20Plots.px,
                   titlefontsize=10,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=size(sp, 1) <= 64 ? 6 : 5)

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
                               xticks=xt,
                               xscale=xsc,
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
- `sp::Matrix{Float64}`: powers
- `db::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
    - `:linlog`: linear-log10
    - `:loglog`: log10-log10
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd_avg(sf::Vector{Float64}, sp::Matrix{Float64}; db::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax,[:linlin, :loglin, :linlog, :loglog], "ax")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = msci95(sp)

    if ax === :linlin
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        ysc = :identity
    elseif ax === :linlog
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = !db ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        ysc = !db ? :log10 : :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   xticks=xt,
                   xscale=xsc,
                   yscale=ysc,
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=10,
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
- `db::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
    - `:linlog`: linear-log10
    - `:loglog`: log10-log10
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd_butterfly(sf::Vector{Float64}, sp::Matrix{Float64}; clabels::Vector{String}=[""], db::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    if ax === :linlin
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        ysc = :identity
    elseif ax === :linlog
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = !db ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        ysc = !db ? :log10 : :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   xticks=xt,
                   xscale=xsc,
                   yscale=ysc;
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=10,
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
                        legend=false;
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
- `db::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `zlabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
    - `:linlog`: linear-log10
    - `:loglog`: log10-log10
- `variant::Symbol`: waterfall (`:w`) or surface (`:s`)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd_3d(sf::Vector{Float64}, sp::Matrix{Float64}; clabels::Vector{String}=[""], db::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", zlabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, variant::Symbol, kwargs...)

    _check_var(variant, [:w, :s], "variant")
    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")
    _check_tuple(frq_lim, "frq_lim")

    ch_n = size(sp, 1)

    pal = mono ? :grays : :darktest

    # channel labels
    clabels == [""] && (clabels = repeat([""], ch_n))

    if ax === :linlin
        xt = _ticks(frq_lim)
        xsc = :identity
        zsc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        zsc = :identity
    elseif ax === :linlog
        xt = _ticks(frq_lim)
        xsc = :identity
        zsc = !db ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        zsc = !db ? :log10 : :identity
    end

    # prepare plot
    if variant === :w
        p = Plots.plot3d(sf,
                         ones(length(sf)),
                         sp[1, :],
                         xlabel=xlabel,
                         ylabel="",
                         zlabel=zlabel,
                         legend=false,
                         xlims=frq_lim,
                         xticks=xt,
                         xscale=xsc,
                         zscale=zsc;
                         title=title,
                         palette=pal,
                         st=:line,
                         lc=:black,
                         size=size(sp, 1) <= 64 ? (900, 900) : (1400, 1400),
                         margins=-50Plots.px,
                         xrotation=-15,
                         yrotation=-10,
                         titlefontsize=10,
                         xlabelfontsize=8,
                         ylabelfontsize=8,
                         zlabelfontsize=8,
                         xtickfontsize=6,
                         ytickfontsize=size(sp, 1) <= 64 ? 6 : 4,
                         ztickfontsize=6)

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
        p = Plots.plot3d(sf[f1:f2],
                         eachindex(clabels),
                         sp[:, f1:f2],
                         xlabel=xlabel,
                         ylabel="",
                         zlabel=zlabel,
                         legend=false,
                         xlims=frq_lim,
                         xticks=xt,
                         xscale=xsc,
                         zscale=zsc;
                         title=title,
                         palette=pal,
                         st=:surface,
                         lc=:black,
                         size=size(sp, 1) <= 64 ? (900, 800) : (1400, 1400),
                         margins=-10Plots.px,
                         xrotation=-15,
                         yrotation=-10,
                         titlefontsize=10,
                         xlabelfontsize=8,
                         ylabelfontsize=8,
                         zlabelfontsize=8,
                         xtickfontsize=6,
                         ytickfontsize=size(sp, 1) <= 64 ? 6 : 4,
                         ztickfontsize=6)
    end

    p = Plots.plot!(yticks=(1:ch_n, clabels))

    return p

end

"""
    plot_psd_topo(locs, sf, sp; <keyword arguments>)

Plot topographical map PSDs.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `ch::Union{Vector{Int64}, AbstractRange}`: which channels to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `db::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
    - `:linlog`: linear-log10
    - `:loglog`: log10-log10
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd_topo(locs::DataFrame, sf::Vector{Float64}, sp::Matrix{Float64}; ch=Union{Vector{Int64}, AbstractRange}, clabels::Vector{String}=[""], db::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, cart::Bool=false, kwargs...)

    @assert length(ch) == nrow(locs) "Some channels do not have locations."
    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    if ax === :linlin
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        ysc = :identity
    elseif ax === :linlog
        xt = _ticks(frq_lim)
        xsc = :identity
        ysc = !db ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.1 Hz")
            sf[1] == 0 && (sf[1] = 0.1)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)))
        end
        xsc = :log10
        ysc = !db ? :log10 : :identity
    end

    # plot parameters
    if size(sp, 1) <= 64
        plot_size = 800
        marker_size = (120, 80)
        rx = 1
        ry = 1
        offset = 0
    elseif size(sp, 1) <= 100
        size(sp, 1) > 64
        plot_size = 1200
        marker_size = (100, 60)
        rx = 0.9
        ry = 0.7
        offset = 150
    else
        plot_size = 2500
        marker_size = (90, 50)
        rx = 0.6
        ry = 0.4
        offset = 250
    end

    # get locations
    if !cart
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
    loc_x .*= ((plot_size / 2) - marker_size[1] / 2) .* rx
    loc_y .*= ((plot_size / 2) - marker_size[2] / 2) .* ry
    # origin is in the left top corner, convert positions
    loc_x = round.(Int64, loc_x .+ (plot_size / 2) .- marker_size[1] / 2)
    loc_y = (plot_size - marker_size[2]) .- round.(Int64, loc_y .+ (plot_size / 2) .- marker_size[2] / 2)

    c = CairoRGBSurface(plot_size, plot_size - 3 * offset)
    cr = CairoContext(c)
    Cairo.set_source_rgb(cr, 256, 256, 256)
    Cairo.rectangle(cr, 0.0, 0.0, plot_size, plot_size - 3 * offset)
    Cairo.fill(cr)
    for idx in 1:size(sp, 1)
        p = Plots.plot(sf,
                       sp[idx, :],
                       t=:line,
                       c=:black,
                       showaxis=false,
                       linewidth=0.5,
                       xlabel=xlabel,
                       ylabel=ylabel,
                       legend=false,
                       xlims=frq_lim,
                       xticks=false,
                       yticks=false,
                       xscale=xsc,
                       yscale=ysc,
                       title=clabels[idx],
                       palette=pal,
                       size=marker_size,
                       titlefontsize=10,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6;
                       kwargs...)
        p = Plots.hline!([minimum(sp[idx, :])], c=:black, linewidth=0.5, xlims=(0, sf[end]))
        p = Plots.vline!([0], c=:black, linewidth=0.5, ylims=(minimum(sp[idx, :]), maximum(sp[idx, :])))
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        Cairo.set_source_surface(cr, img, loc_x[idx], loc_y[idx] - 1.5 * offset)
        Cairo.paint(cr)
    end
    img_png = tempname() * ".png"
    Cairo.write_to_png(c, img_png)
    img = FileIO.load(img_png)
    p = Plots.plot(img,
                   size=(plot_size, plot_size - 3 * offset),
                   title=title,
                   titlefontsize=12,
                   border=:none)
    rm(img_png)

    return p

end

"""
    plot_psd(obj; <keyword arguments>)

Plot power spectrum density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `db::Bool=true`: normalize powers to dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher`
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
    - `:linlog`: linear-log10
    - `:loglog`: log10-log10
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is `Power [dB units^2/Hz] or Power [units^2/Hz]`
- `zlabel::String="default"`: z-axis label for 3-d plots, default is `Power [dB units^2/Hz] or Power [units^2/Hz]`
- `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or gray palette
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:butterfly`
    - `:mean`
    - `:w3d`: 3-d waterfall
    - `:s3d`: 3-d surface
    - `:topo`: topographical
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, ch::Union{String, Vector{String}}, db::Bool=true, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128), ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, kwargs...) where {T <: CWT}

    _check_var(type, [:normal, :butterfly, :mean, :w3d, :s3d, :topo], "type")
    _check_var(method, [:welch, :fft, :stft, :mt, :mw, :gh, :cwt], "method")
    _check_var(ref, [:abs, :total, :delta, :theta, :alpha, :alpha_lower, :alpha_higher, :beta, :beta_lower, :beta_higher, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "ref")
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    ch = get_channel(obj, ch=ch)
    _check_tuple(frq_lim, "frq_lim", (0, sr(obj) / 2))

    @assert seg[1] != seg[2] "Signal is too short for analysis."

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    if ep != 0
        _check_epochs(obj, ep)
        if nepochs(obj) == 1
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

    clabels = labels(obj)

    # set units
    units = _ch_units(obj, clabels[ch][1])

    ref !== :abs && (frq_lim = band_frq(obj, band=ref))

    # get frequency range
    fs = sr(obj)
    _check_tuple(frq_lim, "frq_lim", (0, sr(obj) / 2))

    # # get time vector
    # _, t_s1, _, t_s2 = _convert_t(obj.epoch_time[1], obj.epoch_time[end])

    if ref === :abs
        if method === :welch
            sp, sf = psd(signal, fs=fs, db=db, method=:welch, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (Welch's periodogram)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd(signal, fs=fs, db=db, method=:fft, w=w)
            title == "default" && (title = "Absolute PSD (fast Fourier transform)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (short-time Fourier transform)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (multi-tapered)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd(signal, fs=fs, db=db, method=:mw, ncyc=ncyc, w=w)
            title == "default" && (title = "Absolute PSD (Morlet wavelet convolution)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "Absolute PSD (Gaussian and Hilbert transform)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf = psd(signal, fs=fs, db=db, method=:cwt, wt=wt)
            title == "default" && (title = "Absolute PSD (CWT)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        end
    elseif ref === :total
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, db=db, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (Welch's periodogram) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:fft, w=w)
            title == "default" && (title = "PSD (fast Fourier transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (short-time Fourier transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (multi-tapered) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mw, ncyc=ncyc, w=w)
            title == "default" && (title = "PSD (Morlet wavelet convolution) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "PSD (Gaussian and Hilbert transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:cwt, wt=wt)
            title == "default" && (title = "PSD (CWT) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        end
    else
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, db=db, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (Welch's periodogram) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:fft, frq_lim=frq_lim, w=w)
            title == "default" && (title = "PSD (fast Fourier transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:stft, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (short-time Fourier transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mt, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (multi-tapered) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mw, frq_lim=frq_lim, ncyc=ncyc, w=w)
            title == "default" && (title = "PSD (Morlet wavelet convolution) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "PSD (Gaussian and Hilbert transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:cwt, wt=wt)
            title == "default" && (title = "PSD (CWT) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        end
    end

    # set labels
    if type !== :w3d && type !== :s3d && type !== :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if ref !== :abs
            ylabel == "default" && (ylabel = "Power ratio")
        end
        if method === :cwt
            ylabel == "default" && (ylabel = "Magnitude")
        else
            ylabel == "default" && (ylabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        end
    end

    if type === :normal
        ch_t = obj.header.recording[:channel_type]
        if length(ch) > 1
            ch_t_uni = unique(ch_t[ch])
            @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels should be of the same type."
        end
        if size(sp, 1) == 1
            p = plot_psd(sf,
                         sp[1, :],
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         db=db,
                         frq_lim=frq_lim,
                         ax=ax,
                         mono=mono;
                         kwargs...)
        else
            p = plot_psd(sf,
                         sp,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         clabels=clabels[ch],
                         title=title,
                         db=db,
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
                               clabels=clabels[ch],
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               db=db,
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
                         db=db,
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
        if method === :cwt
            zlabel == "default" && (zlabel = "Magnitude")
        else
            zlabel == "default" && (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        end
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(sf,
                        sp,
                        clabels=clabels[ch],
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        db=db,
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
        if method === :cwt
            zlabel == "default" && (zlabel = "Magnitude")
        else
            zlabel == "default" && (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        end
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(sf,
                        sp,
                        clabels=clabels[ch],
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        db=db,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:s;
                        kwargs...)
    elseif type === :topo
        @assert length(unique(obj.header.recording[:channel_type][ch])) == 1 "For multi-channel PSD plots all channels should be of the same type."
        _has_locs(obj)
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        @assert length(ch) == nrow(locs) "Some channels do not have locations."
        ndims(sp) == 1 && (sp = reshape(sp, 1, length(sp)))
        xlabel == "default" && (xlabel = "")
        ylabel == "default" && (ylabel = "")
        title = replace(title, "channel" => "channels")
        p = plot_psd_topo(locs,
                          sf,
                          sp,
                          ch=collect(1:nrow(locs)),
                          clabels=obj.locs[!, :label],
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          db=db,
                          frq_lim=frq_lim,
                          ax=ax,
                          mono=mono;
                          kwargs...)
    end

    Plots.plot(p)

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
- `db::Bool=true`: normalize powers to dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher`
- `ax::Symbol=:linlin`: type of axes scaling:
    - `:linlin`: linear-linear
    - `:loglin`: log10-linear
    - `:linlog`: linear-log10
    - `:loglog`: log10-log10
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is `Power [dB units^2/Hz] or Power [units^2/Hz]`
- `zlabel::String="default"`: z-axis label for 3-d plots, default is `Power [dB units^2/Hz] or Power [units^2/Hz]`
- `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or gray palette
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:butterfly`
    - `:mean`
    - `:w3d`: 3-d waterfall
    - `:s3d`: 3-d surface
    - `:topo`: topographical
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, db::Bool=true, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128), ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, kwargs...) where {T <: CWT}

    _check_var(type, [:normal, :butterfly, :mean, :w3d, :s3d, :topo], "type")
    _check_var(method, [:welch, :fft, :stft, :mt, :mw, :gh, :cwt], "method")
    _check_var(ref, [:abs, :total, :delta, :theta, :alpha, :alpha_lower, :alpha_higher, :beta, :beta_lower, :beta_higher, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "ref")
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    @assert seg[1] != seg[2] "Signal is too short for analysis."

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    if ep != 0
        _check_epochs(obj, ep)
        if nepochs(obj) == 1
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

    # select component channel, default is all channels
    c isa Symbol && (c = _get_component(obj, c))
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    clabels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (clabels = [clabels])

    ref !== :abs && (frq_lim = band_frq(obj, band=ref))

    # get frequency range
    fs = sr(obj)
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

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
            sp, sf = psd(signal, fs=fs, db=db, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (Welch's periodogram)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd(signal, fs=fs, db=db, method=:fft, w=w)
            title == "default" && (title = "Absolute PSD (fast Fourier transform)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (short-time Fourier transform)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (multi-tapered)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd(signal, fs=fs, db=db, method=:mw, ncyc=ncyc, w=w)
            title == "default" && (title = "Absolute PSD (Morlet wavelet convolution)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "Absolute PSD (Gaussian and Hilbert transform)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf = psd(signal, fs=fs, db=db, method=:cwt, wt=wt)
            title == "default" && (title = "Absolute PSD (CWT)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end
    elseif ref === :total
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, db=db, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (Welch's periodogram) relative to total power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:fft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (fast Fourier transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (short-time Fourier transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (multi-tapered) relative to total power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mw, ncyc=ncyc, w=w)
            title == "default" && (title = "PSD (Morlet wavelet convolution) relative to total power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "PSD (Gaussian and Hilbert transform) relative to total power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:cwt, wt=wt)
            title == "default" && (title = "PSD (CWT) relative to total power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end
    else
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, db=db, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) relative to $ref power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:fft, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (fast Fourier transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:stft, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (short-time Fourier transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mt, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (multi-tapered) relative to $(replace(string(ref), "_"=>" ")) power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mw, frq_lim=frq_lim, ncyc=ncyc, w=w)
            title == "default" && (title = "PSD (Morlet wavelet convolution) relative to $(replace(string(ref), "_"=>" ")) power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "PSD (Gaussian and Hilbert transform) relative to $(replace(string(ref), "_"=>" ")) power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:cwt, wt=wt)
            title == "default" && (title = "PSD (CWT) relative to $(replace(string(ref), "_"=>" ")) power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end
    end

    # set labels
    if type !== :w3d && type !== :s3d && type !== :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if ref !== :abs
            ylabel == "default" && (ylabel = "Power ratio")
        end
        if method === :cwt
            ylabel == "default" && (ylabel = "Magnitude")
        else
            ylabel == "default" && (ylabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        end
    end

    if type === :normal
        @assert ndims(sp) == 1 "For type=:normal the signal must contain 1 c_idx."
        p = plot_psd(sf,
                     sp,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     title=title,
                     db=db,
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
                               db=db,
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
                         db=db,
                         frq_lim=frq_lim,
                         ax=ax,
                         mono=mono;
                         kwargs...)
    elseif type === :w3d
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        if method === :cwt
            zlabel == "default" && (zlabel = "Magnitude")
        else
            zlabel == "default" && (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        end
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(sf,
                        sp,
                        clabels=clabels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        db=db,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:w;
                        kwargs...)
    elseif type === :s3d
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        if method === :cwt
            zlabel == "default" && (zlabel = "Magnitude")
        else
            zlabel == "default" && (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        end
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(sf,
                        sp,
                        clabels=clabels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        db=db,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:s;
                        kwargs...)
    end

    Plots.plot(p)

    return p

end
