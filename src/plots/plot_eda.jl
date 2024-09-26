export plot_eda
export plot_eda_butterfly
export plot_eda_avg

"""
    plot_eda(t, s, bad; <keyword arguments>)

Plot EDA.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractVector`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_eda(t::Union{AbstractVector, AbstractRange}, s::AbstractVector; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :darktest

    # get limits
    ylim = (0, ceil(maximum(s) * 1.1, digits=0))
    yticks = [0, ylim[2]]

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_ticks(t), string.(_ticks(t))),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    p = Plots.plot!(t,
                    s,
                    linewidth=1,
                    label="",
                    color=:black)

    return p

end

"""
    plot_eda(t, s, bad; <keyword arguments>)

Plot EDA.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractMatrix`: data to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_eda(t::Union{AbstractVector, AbstractRange}, s::AbstractMatrix; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    ch_n = size(s, 1)
    @assert size(s, 2) == length(t) "Length of EDA values must equal length of time points vector."

    pal = mono ? :grays : :darktest

    # reverse so 1st channel is on top
    s = @views reverse(s[:, eachindex(t)], dims = 1)
    # also, reverse colors if palette is not mono
    if mono
        pal = :grays
        channel_color = repeat([:black], ch_n)
    else
        pal = :darktest
        channel_color = ch_n:-1:1
    end

    # normalize and shift so all channels are visible
    # each channel is between 0 and +1.0
    for idx in 1:ch_n
        # scale by 0.5 so maxima do not overlap
        s[idx, :] = @views normalize_n(s[idx, :]) .+ (idx - 1)
    end

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(s, 1)))

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   xlims=_xlims(t),
                   xticks=(_ticks(t), string.(_ticks(t))),
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot channels
    for idx in 1:ch_n
        p = @views Plots.plot!(t,
                               s[idx, :],
                               label="",
                               color=channel_color[idx])
    end

    # plot labels
    p = Plots.plot!(yticks=((ch_n - 1):-1:0, clabels))

    return p

end

"""
    plot_eda_butterfly(t, s; <keyword arguments>)

Butterfly plot of EDA.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `avg::Bool=false`: plot average EDA
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_eda_butterfly(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, avg::Bool=true, kwargs...)::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :darktest

    ch_n = size(s, 1)

    # get limits
    ylim = (0, ceil(maximum(s) * 1.1, digits=0))
    yticks = [0, ylim[2]]

    # plot channels
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_ticks(t), string.(_ticks(t))),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot signals
    for idx in 1:ch_n
        if clabels == [""]
            p = Plots.plot!(t,
                            s[idx, :],
                            t=:line,
                            linecolor=idx,
                            linewidth=0.2,
                            alpha=0.2,
                            legend=false)
        else
            if clabels == repeat([""], ch_n)
                p = Plots.plot!(t,
                                s[idx, :],
                                t=:line,
                                legend=false,
                                linecolor=idx,
                                linewidth=0.5,
                                alpha=0.5)
            else
                p = Plots.plot!(t,
                                s[idx, :],
                                t=:line,
                                label=clabels[idx],
                                linecolor=idx,
                                linewidth=0.5,
                                alpha=0.5)
            end
        end
    end

    # plot averaged EDA
    if avg
        if ch_n == 1
            s = mean(s, dims=2)[:]
        else
            s = mean(s, dims=1)[:]
        end
        p = Plots.plot!(t,
                        s,
                        linewidth=2,
                        linecolor=:black,
                        label=false)
    end

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    return p

end

"""
    plot_eda_avg(t, s; <keyword arguments>)

Plot EDA amplitude mean and ±95% CI.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_eda_avg(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = msci95(s)

    # get limits
    ylim = (0, round(minimum(s) * 1.5, digits=-2))
    yticks = [0, ylim[2]]

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_ticks(t), string.(_ticks(t))),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot upper 95% CI
    p = Plots.plot!(t,
                    s_u,
                    fillrange=s_l,
                    fillalpha=0.35,
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot lower 95% CI
    p = Plots.plot!(t,
                    s_l,
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot mean
    p = Plots.plot!(t,
                    s_m,
                    label=false,
                    t=:line,
                    c=:black,
                    lw=0.5)

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    return p

end
