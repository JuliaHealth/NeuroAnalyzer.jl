export plot_mep
export plot_mep_stack

"""
    plot_mep(t, s, bad; <keyword arguments>)

Plot MEP (single channel).

# Arguments

  - `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  - `s::AbstractVector`: data to plot
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `zl::Bool=true`: draw line at t = 0
  - `yrev::Bool=false`: reverse y-axis
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_mep(
        t::Union{AbstractVector, AbstractRange},
        s::AbstractVector;
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        zl::Bool = true,
        yrev::Bool = false,
        mono::Bool = false,
    )::GLMakie.Figure

    # prepare plot
    GLMakie.activate!(title = "plot_mep()")
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
        yticks = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        yreversed = yrev,
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yrev ? reverse(_ylims(s) .* 1.5) : (_ylims(s) .* 1.5))
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot 0 v-line
    if zl
        GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)
    end

    # plot MEP
    GLMakie.lines!(ax, t, s; color = :black, linewidth = 1)

    return p

end

"""
    plot_mep(t, s; <keyword arguments>)

Plot MEP (multi-channel).

# Arguments

  - `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  - `s::AbstractMatrix`: data to plot
  - `clabels::Vector{String}=string.(1:size(s, 1))`: signal channel labels vector
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `yrev::Bool=false`: reverse y-axis
  - `avg::Bool=true`: if true, plot averaged MEP
  - `ci95::Bool=false`: if true, plot mean and ±95% CI
  - `leg::Bool=true`: if true, add legend with channel labels
  - `zl::Bool=true`: draw line at t = 0
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_mep(
        t::Union{AbstractVector, AbstractRange},
        s::AbstractMatrix;
        clabels::Vector{String} = string.(1:size(s, 1)),
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        yrev::Bool = false,
        avg::Bool = true,
        ci95::Bool = false,
        leg::Bool = true,
        zl::Bool = true,
        mono::Bool = false,
    )::GLMakie.Figure

    pal = mono ? :grays : :darktest

    ch_n = size(s, 1)

    # prepare plot
    GLMakie.activate!(title = "plot_mep()")
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
        yticks = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(10),
        yreversed = yrev,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yrev ? reverse(_ylims(s) .* 1.5) : (_ylims(s) .* 1.5))
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot 0 v-line
    if zl
        GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)
    end

    # plot MEPs
    if ci95
        avg = false
        leg = false
        s_m, _, s_u, s_l = NeuroAnalyzer.msci95(s)
        # draw 95% CI
        Makie.band!(ax, t, s_u, s_l; alpha = 0.25, color = :grey, strokewidth = 0.5)

        # draw mean
        Makie.lines!(ax, t, s_m; color = :black, linewidth = 2)
    else
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            GLMakie.lines!(
                ax,
                t,
                s[idx, :];
                color = cmap[idx],
                colormap = pal,
                colorrange = 1:ch_n,
                linewidth = 1,
                alpha = avg ? 0.25 : 1.0,
                label = clabels[idx],
            )
        end
    end

    # plot averaged MEP
    if avg
        if ch_n == 1
            s = mean(s; dims = 2)[:]
        else
            s = mean(s; dims = 1)[:]
        end
        GLMakie.lines!(ax, t, s; color = :black, linewidth = 2)
    end

    (leg && ch_n < 30) && axislegend(; position = :rt, colormap = pal)

    return p

end

"""
    plot_mep_stack(s; <keyword arguments>)

Plot MEPs stacked by channels or by epochs.

# Arguments

  - `t::AbstractVector`: x-axis values
  - `s::AbstractMatrix`
  - `clabels::Vector{String}=string.(1:size(s, 1))`: signal channel labels vector
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `cb::Bool=true`: plot color bar
  - `cb_title::String=""`: color bar title
  - `smooth::Bool=false`: smooth the image using Gaussian blur
  - `ks::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
  - `zl::Bool=true`: draw line at t = 0
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_mep_stack(
        t::AbstractVector,
        s::AbstractArray;
        clabels::Vector{String} = string.(1:size(s, 1)),
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        cb::Bool = true,
        cb_title::String = "",
        smooth::Bool = false,
        ks::Int64 = 3,
        zl::Bool = true,
        mono::Bool = false,
    )::GLMakie.Figure

    @assert length(t) == size(s, 2) "Number of s columns ($(size(s, 2))) must equal length of t ($(length(t)))."

    pal = mono ? :grays : :darktest

    if smooth
        s = imfilter(s, Kernel.gaussian(ks))
    end

    # prepare plot
    GLMakie.activate!(title = "plot_mep()")
    plot_size = size(s, 1) <= 64 ? (1200, 800) : (1200, 1200)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        yticks = (axes(s, 1), clabels),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        yticklabelsize = size(s, 1) <= 64 ? 8 : 5,
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    hm = GLMakie.heatmap!(ax, t, axes(s, 1), rotr90(s); colormap = pal)

    # plot 0 v-line
    if zl
        GLMakie.vlines!(ax, 0; color = :white, linestyle = :dash, linewidth = 2)
    end

    # draw colorbar
    if cb
        Colorbar(p[1, 2], hm; label = cb_title, labelsize = 16)
    end

    return p

end

"""
    plot_mep(obj; <keyword arguments>)

Plot MEP.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  - `ch::Union{String, Vector{String}, Regex}`: channel name(s)
  - `xlabel::String="default"`: x-axis label
  - `ylabel::String="default"`: y-axis label
  - `title::String="default"`: plot title
  - `cb::Bool=true`: plot color bar
  - `cb_title::String="default"`: color bar title
  - `peaks::Bool=true`: draw peaks
  - `leg::Bool=true`: if true, add legend with channel labels
  - `type::Symbol=:normal`: multi-channel plot type:
      + `:normal`: butterfly or mean and ±95% CI
      + `:stack`: stacked channels
  - `yrev::Bool=false`: reverse y-axis
  - `avg::Bool=true`: if true, plot averaged MEP
  - `ci95::Bool=false`: if true, plot mean and ±95% CI
  - `smooth::Bool=false`: smooth the image using Gaussian blur
  - `ks::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
  - `zl::Bool=true`: draw line at t = 0
  - `mono::Bool=false`: use color or gray palette
  - `gui::Bool=false`: ignored

# Returns

  - `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_mep(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        xlabel::String = "default",
        ylabel::String = "default",
        title::String = "default",
        cb::Bool = true,
        cb_title::String = "default",
        peaks::Bool = true,
        leg::Bool = true,
        type::Symbol = :normal,
        yrev::Bool = false,
        avg::Bool = true,
        ci95::Bool = false,
        smooth::Bool = false,
        ks::Int64 = 3,
        zl::Bool = true,
        mono::Bool = false,
        gui::Bool = false,
    )::GLMakie.Figure

    _check_datatype(obj, "mep")
    _check_var(type, [:normal, :stack], "type")

    # check channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    @assert !(length(ch) > 1 && length(unique(obj.header.recording[:channel_type][ch])) > 1) "All channels must be of the same type."

    # set units
    units = _ch_units(obj, labels(obj)[ch[1]])

    # get data
    ep_n = nepochs(obj)
    if length(ch) == 1
        s = obj.data[ch, :, 1][:]
    else
        s = obj.data[ch, :, 1]
    end

    # get labels
    clabels = labels(obj)[ch]

    # get time vector
    t = obj.epoch_time

    if length(ch) == 1
        xl, yl, tt = NeuroAnalyzer._set_defaults(
            xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "MEP amplitude, channel: $(clabels[1])"
        )
        p = plot_mep(t, s; xlabel = xl, ylabel = yl, title = tt, mono = mono, yrev = yrev, zl = zl)
    elseif type === :normal
        xl, yl, tt = NeuroAnalyzer._set_defaults(
            xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "MEP amplitude, $(length(ch)) channels"
        )
        p = plot_mep(
            t,
            s;
            xlabel = xl,
            ylabel = yl,
            title = tt,
            clabels = clabels,
            mono = mono,
            yrev = yrev,
            avg = avg,
            ci95 = ci95,
            leg = leg,
            zl = zl,
        )
    elseif type === :stack
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "MEP amplitude, $(length(ch)) channels")
        cb_title == "default" && (cb_title = "Amplitude [$units]")
        p = plot_mep_stack(
            t,
            s;
            xlabel = xl,
            ylabel = yl,
            title = tt,
            clabels = clabels,
            cb = cb,
            cb_title = cb_title,
            mono = mono,
            ks = ks,
            smooth = smooth,
            zl = zl,
        )
    end

    # draw peaks
    if peaks
        if length(ch) == 1
            pp = erp_peaks(obj)
            GLMakie.scatter!(
                p[1, 1],
                t[pp[ch, 1]][1],
                obj.data[ch, pp[ch, 1], 1][1];
                marker = :xcross,
                color = mono ? :black : :red,
                markersize = 15,
            )
            GLMakie.scatter!(
                p[1, 1],
                t[pp[ch, 2]][1],
                obj.data[ch, pp[ch, 2], 1][1];
                marker = :xcross,
                color = mono ? :black : :blue,
                markersize = 15,
            )
            _info("Positive peak time: $(round(t[pp[ch, 1]][1] * 1000, digits = 0)) ms")
            _info("Positive peak amplitude: $(round(obj.data[ch, pp[ch, 1], 1][1], digits = 2)) $units")
            _info("Negative peak time: $(round(t[pp[ch, 2]][1] * 1000, digits = 0)) ms")
            _info("Negative peak amplitude: $(round(obj.data[ch, pp[ch, 2], 1][1], digits = 2)) $units")
        elseif length(ch) > 1 && type === :normal
            mep_tmp = mean(obj.data[ch, :, 1]; dims = 1)[:, :, :]
            obj_tmp = keep_channel(obj, ch = labels(obj)[1])
            obj_tmp.data = mep_tmp
            pp = erp_peaks(obj_tmp)
            GLMakie.scatter!(
                p[1, 1], t[pp[1, 1]], mep_tmp[pp[1, 1]]; marker = :xcross, color = mono ? :black : :red, markersize = 15
            )
            GLMakie.scatter!(
                p[1, 1],
                t[pp[1, 2]],
                mep_tmp[pp[1, 2]];
                marker = :xcross,
                color = mono ? :black : :blue,
                markersize = 15,
            )
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits = 0)) ms")
            _info("Positive peak amplitude: $(round(mep_tmp[pp[1, 1]], digits = 2)) $units")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits = 0)) ms")
            _info("Negative peak amplitude: $(round(mep_tmp[pp[1, 2]], digits = 2)) $units")
        end
    end

    return p

end
