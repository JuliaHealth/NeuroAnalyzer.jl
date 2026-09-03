export plot_mep
export plot_mep_stack

"""
    plot_mep(t, s, bad; <keyword arguments>)

Plot a single-channel Motor Evoked Potential (MEP) waveform with customizable visualization.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: vector of time points in seconds
- `s::AbstractVector`: signal amplitude values
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `yrev::Bool=false`: if `true`, reverse the y-axis
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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
    GLMakie.activate!(; title = "plot_mep()")
    fig = GLMakie.Figure(; size = (900, 450))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xlabel,
        ylabel             = ylabel,
        title              = title,
        xticks             = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        yticks             = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks        = IntervalsBetween(10),
        yreversed          = yrev,
        xautolimitmargin   = (0, 0),
        yautolimitmargin   = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yrev ? reverse(_ylims(s) .* 1.5) : (_ylims(s) .* 1.5))
    _style_axis!(ax)

    # draw zero line if requested
    zl && GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)

    GLMakie.lines!(ax, t, s; color = :black, linewidth = 1)

    return fig
end

"""
    plot_mep(t, s; <keyword arguments>)

Plot multi-channel Motor Evoked Potentials (MEPs) with customizable visualization and averaging options.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: vector of time points in seconds
- `s::AbstractMatrix`: signal amplitude values, shape (channel, samples)
- `clabels::Vector{String}=string.(1:size(s, 1))`: channel labels (default: auto-generated)
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `yrev::Bool=false`: if `true`, reverse the y-axis
- `avg::Bool=true`: if `true`, plot averaged MEP across channels
- `ci95::Bool=false`: if `true`, plot mean and ±95% confidence interval of averaged MEPs
- `leg::Bool=true`: if `true`, show legend with channel labels
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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
    # validate
    avg && ci95 && throw(ArgumentError("avg and ci95 cannot both be true."))

    # set color palette
    pal = mono ? :grays : :darktest

    # number of channels
    ch_n = size(s, 1)

    # prepare plot
    GLMakie.activate!(; title = "plot_mep()")
    fig = GLMakie.Figure(; size = (900, 450))
    ax  = GLMakie.Axis(
    fig[1, 1];
    xlabel             = xlabel,
    ylabel             = ylabel,
    title              = title,
    xticks             = LinearTicks(10),
    xminorticksvisible = true,
    xminorticks        = IntervalsBetween(10),
    yticks             = LinearTicks(10),
    yminorticksvisible = true,
    yminorticks        = IntervalsBetween(10),
    yreversed          = yrev,
    xautolimitmargin   = (0, 0),
    yautolimitmargin   = (0, 0),
    _AXIS_LOCK_KWARGS...
)
    GLMakie.ylims!(ax, yrev ? reverse(_ylims(s) .* 1.5) : (_ylims(s) .* 1.5))
    _style_axis!(ax)

    # draw zero line if requested
    zl && GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)

    if ci95
        # get mean and 95%CI
        msci95_data = NeuroAnalyzer.msci95(s)
        s_m = msci95_data.sm
        s_u = msci95_data.ul
        s_l = msci95_data.ll
        # draw 95% CI
        GLMakie.band!(ax, t, s_u, s_l; alpha = 0.25, color = :grey, strokewidth = 0.5)
        # draw mean
        GLMakie.lines!(ax, t, s_m; color = :black, linewidth = 2)
    else
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx = 1:ch_n
            GLMakie.lines!(
                ax, t, s[idx, :];
                color      = cmap[idx],
                colormap   = pal,
                colorrange = 1:ch_n,
                linewidth  = 1,
                alpha      = avg ? 0.25 : 1.0,
                label      = clabels[idx],
            )
        end

        # draw averaged channels
        if avg
            s_avg = mean(s; dims = 1)[:]
            GLMakie.lines!(ax, t, s_avg; color = :black, linewidth = 2)
        end

        # add legend if requested
        (leg && ch_n < 30) && axislegend(; position = :rt, colormap = pal)
    end

    return fig
end

"""
    plot_mep_stack(s; <keyword arguments>)

Plot Motor Evoked Potentials (MEPs) stacked by channels or epochs with customizable visualization.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: vector of time points in seconds
- `s::AbstractMatrix`: signal amplitude values, shape (channel, samples)
- `clabels::Vector{String}=string.(1:size(s, 1))`: channel labels (default: auto-generated)
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `cb::Bool=true`: if `true`, show color bar
- `cb_title::String=""`: color bar title
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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
    # validate
    length(t) == size(s, 2) || throw(
        ArgumentError(
            "Number of s columns ($(size(s, 2))) must equal length of t ($(length(t))).",
        ),
    )

    # set color palette
    pal = mono ? :grays : :darktest

    # apply Gaussian filter if requested
    smooth && (s = imfilter(s, Kernel.gaussian(ks)))

    ytick_size = size(s, 1) <= 64 ? 8 : 5

    # prepare plot
    GLMakie.activate!(; title = "plot_mep_stack()")
    plot_size = size(s, 1) <= 64 ? (1200, 800) : (1200, 1200)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xlabel,
        ylabel             = ylabel,
        title              = title,
        xticks             = LinearTicks(10),
        yticks             = (axes(s, 1), clabels),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        xautolimitmargin   = (0, 0),
        yautolimitmargin   = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)
    ax.yticklabelsize = ytick_size  # conditional, overwrite _style_axis!

    hm = GLMakie.heatmap!(ax, t, axes(s, 1), rotr90(s); colormap = pal)

    # draw zero line if requested
    zl && GLMakie.vlines!(ax, 0; color = :white, linestyle = :dash, linewidth = 2)

    # draw zero line if requested
    cb && GLMakie.Colorbar(fig[1, 2], hm; label = cb_title, labelsize = 16)

    return fig
end

"""
    plot_mep(obj; <keyword arguments>)

Plot Motor Evoked Potentials (MEPs) from a NEURO object with customizable visualization options.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `cb::Bool=true`: if `true`, show color bar
- `cb_title::String="default"`: color bar title
- `peaks::Symbol=:detect`: method for drawing peaks:
    - `:detect`: detect and draw peaks automatically
    - `:embed`: embed peaks in the plot
    - `:off`: do not draw peaks
- `leg::Bool=true`: if `true`, show legend with channel labels
- `type::Symbol=:normal`: multi-channel plot type:
    - `:normal`: butterfly or mean and ±95% CI
    - `:stack`: stacked channels
- `yrev::Bool=false`: if `true`, reverse the y-axis
- `avg::Bool=true`: if `true`, plot averaged MEP
- `ci95::Bool=false`: if `true`, plot mean and ±95% CI
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette
- `gui::Bool=false`: ignored parameter (kept for compatibility)

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_mep(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    cb::Bool = true,
    cb_title::String = "default",
    peaks::Symbol = :detect,
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
    # validate
    _check_datatype(obj, "mep")
    _check_var(type, [:normal, :stack], "type")
    _check_var(peaks, [:detect, :embed, :off], "peaks")

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    (length(ch) > 1 && length(unique(obj.header.recording[:channel_type][ch])) > 1) &&
        throw(ArgumentError("All channels must be of the same type."))

    # set units
    units = _ch_units(obj, labels(obj)[ch[1]])

    # get data
    s = length(ch) == 1 ? obj.data[ch, :, 1][:] : obj.data[ch, :, 1]

    # get labels
    clabels = labels(obj)[ch]

    # get time vector
    t = obj.epoch_time

    if length(ch) == 1
        xl, yl, tt = NeuroAnalyzer._set_defaults(
            xlabel, ylabel, title,
            "Time [ms]", "Amplitude [$units]",
            "MEP amplitude, channel: $(clabels[1])",
        )
        fig = plot_mep(t, s; xlabel = xl, ylabel = yl, title = tt,
            mono = mono, yrev = yrev, zl = zl)

    elseif type === :normal
        xl, yl, tt = NeuroAnalyzer._set_defaults(
            xlabel, ylabel, title,
            "Time [ms]", "Amplitude [$units]",
            "MEP amplitude, $(length(ch)) channels",
        )
        fig = plot_mep(t, s; xlabel = xl, ylabel = yl, title = tt,
            clabels = clabels, mono = mono, yrev = yrev,
            avg = avg, ci95 = ci95, leg = leg, zl = zl)

    elseif type === :stack
        xl, yl, tt = NeuroAnalyzer._set_defaults(
            xlabel, ylabel, title,
            "Time [ms]", "",
            "MEP amplitude, $(length(ch)) channels",
        )
        cb_title == "default" && (cb_title = "Amplitude [$units]")
        fig = plot_mep_stack(t, s; xlabel = xl, ylabel = yl, title = tt,
            clabels = clabels, cb = cb, cb_title = cb_title,
            mono = mono, ks = ks, smooth = smooth, zl = zl)
    end

    # draw peaks - single-channel only
    if peaks !== :off
        if length(ch) == 1
            pp = if peaks === :detect
                mep_peaks(obj)
            elseif peaks === :embed
                hcat(
                    obj.header.recording[:markers_pos],
                    obj.header.recording[:markers_neg],
                )
            end
            GLMakie.scatter!(
                fig[1, 1], t[pp[ch, 1]][1], obj.data[ch, pp[ch, 1], 1][1];
                marker     = :xcross,
                color      = mono ? :black : :red,
                markersize = 15,
            )
            GLMakie.scatter!(
                fig[1, 1], t[pp[ch, 2]][1], obj.data[ch, pp[ch, 2], 1][1];
                marker     = :xcross,
                color      = mono ? :black : :blue,
                markersize = 15,
            )
            _info("Positive peak time: $(round(t[pp[ch, 1]][1] * 1000; digits = 0)) ms")
            _info(
                "Positive peak amplitude: $(round(obj.data[ch, pp[ch, 1], 1][1]; digits = 2)) $units",
            )
            _info("Negative peak time: $(round(t[pp[ch, 2]][1] * 1000; digits = 0)) ms")
            _info(
                "Negative peak amplitude: $(round(obj.data[ch, pp[ch, 2], 1][1]; digits = 2)) $units",
            )

        elseif length(ch) > 1 && type === :normal
            mep_tmp = mean(obj.data[ch, :, 1]; dims = 1)
            obj_tmp = keep_channel(obj; ch = labels(obj)[1])
            obj_tmp.data = reshape(mep_tmp, 1, :, 1)
            pp = mep_peaks(obj_tmp)
            GLMakie.scatter!(
                fig[1, 1], t[pp[1, 1]], mep_tmp[pp[1, 1]];
                marker     = :xcross,
                color      = mono ? :black : :red,
                markersize = 15,
            )
            GLMakie.scatter!(
                fig[1, 1], t[pp[1, 2]], mep_tmp[pp[1, 2]];
                marker     = :xcross,
                color      = mono ? :black : :blue,
                markersize = 15,
            )
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000; digits = 0)) ms")
            _info("Positive peak amplitude: $(round(mep_tmp[pp[1, 1]]; digits = 2)) $units")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000; digits = 0)) ms")
            _info("Negative peak amplitude: $(round(mep_tmp[pp[1, 2]]; digits = 2)) $units")
        end
    end

    return fig
end
