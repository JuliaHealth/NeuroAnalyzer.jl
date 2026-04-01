# TO DO:

# types by colors
# add marker start : end
# time format (SS:MS HH:MM:SS)
# change scaling

export plot_cont

"""
    plot_cont(obj; <keyword arguments>)

Plot continuous signal with interactive editing.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name(s)
- `seg::Tuple{Real, Real}=(0, 10)`: time segment to display in seconds (from, to), default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `markers::Bool`: if `true`, draw markers if available
- `scale::Bool=true`: if `true`, draw scale reference
- `group_ch::Bool=true`: if `true`, group channels by type (e.g. EEG, EOG, ECG)
- `type::Symbol=:normal`: plot type:
    - `:normal`: standard multi-channel plot
    - `:butterfly`: butterfly plot showing all channels overlaid
- `avg::Bool=false`: if `true`, plot averaged channel in butterfly plot
- `ci95::Bool=false`: if `true`, plot mean and ±95% confidence interval of averaged channels in butterfly plot
- `n_channels::Int64=20`: maximum number of visible channels
- `res::Int64=1`: resampling factor (draw every `res`-nth sample)
- `snap::Bool=true`: if `true`, snap to grid when placing markers
- `gui::Bool=true`: if `true`, keep window open and interactive

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_cont(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = "all",
    seg::Tuple{Real, Real} = (0, 10),
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    mono::Bool = false,
    markers::Bool = true,
    scale::Bool = true,
    group_ch::Bool = true,
    type::Symbol = :normal,
    avg::Bool = true,
    ci95::Bool = false,
    n_channels::Int64 = 20,
    res::Int64 = 1,
    snap::Bool = true,
    gui::Bool = true,
)::GLMakie.Figure
    # validate
    res >= 1 || throw(ArgumentError("res must be ≥ 1."))
    res > 10 && _warn("At res > 10 plot will be inaccurate.")
    n_channels >= 1 || throw(ArgumentError("n_channels must be ≥ 1."))
    n_channels <= nchannels(obj) ||
        throw(ArgumentError("n_channels must be ≤ $(nchannels(obj))."))
    _check_var(type, [:normal, :butterfly], "type")
    avg && ci95 && throw(ArgumentError("avg and ci95 cannot both be true."))
    !_has_markers(obj) && (markers = false)

    # set color palette
    pal = mono ? :grays : :darktest

    # set maximum length to display
    if signal_len(obj) <= seg[2] * sr(obj)
        seg = (obj.time_pts[1], obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end

    # check channels and meta data
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    ch_n = length(ch)
    if group_ch
        ch_order = _sort_channels(obj.header.recording[:channel_type][ch])
    else
        ch_order = collect(1:ch_n)
    end
    ctypes  = obj.header.recording[:channel_type][ch][ch_order]
    cunits  = obj.header.recording[:unit][ch][ch_order]

    # order by ctypes
    # and markers for ax3
    ctypes_uni = unique(ctypes)
    ctypes_pos = zeros(Int64, length(ctypes_uni))
    for idx in eachindex(ctypes_uni)
        ctypes_pos[idx] = findfirst(isequal(ctypes_uni[idx]), ctypes)
    end
    # ctypes_uni_pos contains a list of ticks where scale has to be drawn
    ctypes_uni_pos = zeros(Int64, ch_n)
    ctypes_uni_pos[ctypes_pos] .= 1

    # get time points vector
    t = obj.time_pts[1:res:end]
    # get signal matrix
    s = obj.data[ch, :, 1][ch_order, 1:res:end]

    # set defaults
    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "")

    # list of bad channels
    bad_ch = Observable(obj.header.recording[:bad_channel])

    # displayed segment
    seg_pos = Observable(Float64(seg[1]))
    seg_len = Float64(seg[2]) - Float64(seg[1])

    if type === :normal
        nch      = Observable(n_channels)
        ch1      = Observable(1)
        ch2_init = gui && ch_n > nch[] ? ch1[] + nch[] - 1 : ch_n
        clabels = labels(obj)[ch][ch_order]
    else
        ch_n     = length(ctypes_uni)
        ch1      = Observable(1)
        ch2_init = length(ctypes_uni)
        nch      = Observable(ch_n)
        clabels  = uppercase.(ctypes_uni)
    end

    # get ranges of the original signal for the scales
    # normalize in groups by channel type
    # between -0.5 and +0.5 and shift so all channels are visible
    r = Observable(Float64[])
    for idx in eachindex(ctypes_uni)
        group = s[ctypes .== ctypes_uni[idx], :]
        push!(r[], round(_get_range(group)))
        # remove per-channel DC offset
        group = group .- mean(group; dims = 2)
        # map to [-0.5, 0.5]
        s[ctypes .== ctypes_uni[idx], :] = normalize_minmax(group, 0.5; bych = true)
    end
    if type === :normal
        s .+= collect(1:ch_n)
    elseif type === :butterfly
        for idx in eachindex(ctypes_uni)
            s[ctypes .== ctypes_uni[idx], :] .+= idx
        end
    end

    # y-axis labels colors
    if type === :normal
        ytc = repeat([:black], nchannels(obj))
        ytc[bad_ch[]] .= :lightgray
    else
        ytc = repeat([:black], ch_n)
    end

    # prepare markers
    if markers
        markers_pos  = obj.markers[!, :start]
        markers_id   = obj.markers[!, :id]
        markers_desc = obj.markers[!, :value]
    end

    # prepare plot
    plot_size = if gui
        type === :normal ? (1250, 700) : (1200, 700)
    else
        (1200, 650)
    end
    GLMakie.activate!(; title = "plot_cont()")
    fig = GLMakie.Figure(; size = plot_size, figure_padding = (10, 20, 10, 10))

    # create axis with customizable properties
    ax1 = GLMakie.Axis(
        fig[1, 1];
        xlabel             = "",
        ylabel             = yl,
        title              = tt,
        xticks             = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        yticks             = (1:ch_n, clabels),
        # TO DO: yticklabelcolor = ytc
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
        yticklabelspace = let ml = maximum(length, clabels)
            ml <= 5 ? 60.0 : ml >= 10 ? 100.0 : 80.0
        end,
    )
    GLMakie.xlims!(ax1, seg)
    if gui && ch_n > nch[]
        GLMakie.ylims!(ax1, ch2_init + 0.5, ch2_init - nch[] + 0.5)
    else
        GLMakie.ylims!(ax1, ch_n + 0.5, 0.5)
    end
    _style_axis!(ax1)

    # draw channels
    if type === :normal
        for idx in 1:ch_n
            line_color = @lift($bad_ch[idx] ? :lightgray : :black)
            GLMakie.lines!(ax1, t, s[idx, :]; linewidth = 1.5, color = line_color)
        end

    elseif type === :butterfly
        if ci95
            for idx in eachindex(ctypes_uni)
                msci95_data = NeuroAnalyzer.msci95(s[ctypes .== ctypes_uni[idx], :])
                s_m = msci95_data.sm
                s_u = msci95_data.ul
                s_l = msci95_data.ll
                GLMakie.band!(
                    ax1, t, s_u, s_l;
                    alpha = 0.25,
                    color = :grey,
                    strokewidth = 0.5,
                )
                GLMakie.lines!(ax1, t, s_m; color = :black, linewidth = 2)
            end

        else
            !mono && (cmap = GLMakie.resample_cmap(pal, size(s, 1)))
            for idx in axes(s, 1)
                GLMakie.lines!(
                    ax1, t, s[idx, :];
                    color      = mono ? :black : cmap[idx],
                    colormap   = pal,
                    colorrange = 1:size(s, 1),
                    linewidth  = 0.5,
                    alpha = 1.0,
                )
            end
            if avg
                for idx in eachindex(ctypes_uni)
                    s_avg = mean(s[ctypes .== ctypes_uni[idx], :]; dims = 1)[:]
                    GLMakie.lines!(ax1, t, s_avg; linewidth = 2, color = :black)
                end
            end
        end
    end

    # draw scale bars
    # TO DO: place scale values on the left side, below channel label
    if scale
        if type === :normal
            idx2 = 1
            for idx1 in 1:ch_n
                if ctypes_uni_pos[idx1] == 1
                    s_rectangle = lift(seg_pos) do sp
                        return Rect(sp, (idx1 - 0.49), 0.01, 0.98)
                    end
                    l_pos = lift(seg_pos) do sp
                        return (sp + 0.01, idx1 + 0.49)
                    end
                    GLMakie.poly!(
                        ax1,
                        s_rectangle;
                        color = :red,
                        strokecolor = :red,
                        strokewidth = 2,
                    )
                    GLMakie.text!(
                        ax1, l_pos;
                        markerspace = :pixel,
                        text        = string(r[][idx2]) * " " * cunits[idx1],
                        fontsize    = 10,
                        color       = :red,
                        align       = (:left, :bottom),
                        offset      = (5, 0),
                    )
                    idx2 += 1
                end
            end
        elseif type === :butterfly
            for idx in 1:ch_n
                s_rectangle = lift(seg_pos) do sp
                    return Rect(sp, (idx - 0.475), 0.01, 0.975)
                end
                l_pos = lift(seg_pos) do sp
                    return (sp, idx + 0.5)
                end
                GLMakie.poly!(
                    ax1,
                    s_rectangle;
                    color = :red,
                    strokecolor = :red,
                    strokewidth = 2,
                )
                GLMakie.text!(
                    ax1, l_pos;
                    text        = string(r[][idx]) * " " * cunits[ctypes .== ctypes_uni[idx]][1],
                    markerspace = :pixel,
                    fontsize    = 10,
                    color       = :red,
                    align       = (:left, :bottom),
                    offset      = (5, 0),
                )
            end
        end
    end

    # draw event markers
    if markers
        GLMakie.vlines!(ax1, markers_pos; linestyle = :dash, linewidth = 1, color = :black)
        for idx in eachindex(markers_pos)
            markers_ypos = lift(ch1, nch) do v1, v2
                return (markers_pos[idx], v1 + (v2 - 1) + 0.5)
            end
            GLMakie.textlabel!(
                ax1, markers_ypos;
                text = "$(markers_id[idx]) / $(markers_desc[idx])",
                text_align = (:left, :center),
                fontsize = 8,
                cornerradius = 0,
                cornervertices = 2,
                padding = 2,
                strokewidth = 1,
                offset = (0, 5),
                text_rotation = pi / 2,
            )
        end
    end

    vmarker1     = Observable(NaN)
    vmarker2     = Observable(NaN)
    marker_range = Observable([NaN, NaN])

    if gui
        # time bar
        ax2 = GLMakie.Axis(
            fig[2, 1];
            xlabel           = xl,
            ylabel           = "",
            title            = "",
            xticks           = LinearTicks(25),
            yticksvisible    = false,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            backgroundcolor  = :white,
            xzoomlock        = true,
            yzoomlock        = true,
            xpanlock         = true,
            ypanlock         = true,
            xrectzoom        = false,
            yrectzoom        = false,
        )
        GLMakie.xlims!(ax2, t[1], t[end])
        GLMakie.ylims!(ax2, 0, 1)
        hideydecorations!(ax2)
        hidexdecorations!(ax2; label = false, ticks = false, ticklabels = false)
        ax2.xticklabelsize = 12

        markers && GLMakie.vlines!(
            ax2,
            markers_pos;
            linestyle = :dash,
            linewidth = 1,
            color = :black,
        )

        t_rectangle = lift(seg_pos) do v
            return Rect(v, 0, seg_len, 1)
        end
        GLMakie.poly!(
            ax2,
            t_rectangle;
            color = :darkgrey,
            strokecolor = :black,
            strokewidth = 2,
            alpha = 0.5,
        )

        # channel bar (normal type only)
        if type === :normal
            ax3 = GLMakie.Axis(
                fig[1, 2];
                xlabel           = "",
                ylabel           = "",
                title            = "",
                yticks           = 1:ch_n,
                xticksvisible    = false,
                yticksvisible    = false,
                yreversed        = true,
                xautolimitmargin = (0, 0),
                yautolimitmargin = (0, 0),
                backgroundcolor  = :white,
                xzoomlock        = true,
                yzoomlock        = true,
                xpanlock         = true,
                ypanlock         = true,
                xrectzoom        = false,
                yrectzoom        = false,
            )
            ch_n > 1 && GLMakie.ylims!(ax3, ch_n, 1)
            hidedecorations!(ax3)

            if group_ch
                for idx in eachindex(ctypes_pos)
                    GLMakie.hlines!(ax3, ctypes_pos[idx]; linewidth = 5, color = :black)
                end
            end

            ch_rectangle = @lift(Rect(0, $ch1, 1, $nch - 1))
            GLMakie.poly!(
                ax3, ch_rectangle;
                color       = :darkgrey,
                strokecolor = :black,
                strokewidth = 2,
                alpha       = 0.25,
            )
        end

        GLMakie.vlines!(ax1, vmarker1; color = (:blue, 0.8), linewidth = 1)
        GLMakie.vlines!(ax1, vmarker2; color = (:blue, 0.8), linewidth = 1)
        GLMakie.band!(ax1, marker_range, 0.5, ch_n + 0.5; color = (:blue, 0.1))

        # mouse events
        on(events(fig).mousebutton) do event
            ax1_x = mouseposition(ax1)[1]
            ax1_y = mouseposition(ax1)[2]
            ax2_x = mouseposition(ax2)[1]
            ax2_y = mouseposition(ax2)[2]
            if type === :normal
                ax3_x = mouseposition(ax3)[1]
                ax3_y = mouseposition(ax3)[2]
            end

            if event.action == Mouse.press
                if event.button == Mouse.right
                    if type === :normal
                        # mark/unmark channel as bad
                        if ax1_x < ax1.limits[][1][1]
                            bad_ch[][round(Int64, ax1_y)] = !bad_ch[][round(Int64, ax1_y)]
                            obj.header.recording[:bad_channel][
                                get_channel(obj; ch = clabels[round(Int64, ax1_y)])[1],
                            ] =
                                !obj.header.recording[:bad_channel][
                                    get_channel(obj; ch = clabels[round(Int64, ax1_y)])[1],
                                ]
                            notify(bad_ch)
                        end
                        # clear markers
                        if ax1_x >= ax1.limits[][1][1] &&
                           ax1_x <= ax1.limits[][1][2] &&
                           ax1_y >= ax1.limits[][2][1] &&
                           ax1_y <= ax1.limits[][2][2]
                            vmarker1[] = NaN
                            vmarker2[] = NaN
                            marker_range[] = [NaN, NaN]
                            notify(vmarker1);
                            notify(vmarker2);
                            notify(marker_range)
                        end
                    end

                elseif event.button == Mouse.left
                    if type === :normal
                        # get channel info
                        ax1_x < ax1.limits[][1][1] &&
                            channel_info(obj; ch = clabels[round(Int64, ax1_y)])

                        # place marker
                        if ax1_x >= ax1.limits[][1][1] &&
                           ax1_x <= ax1.limits[][1][2] &&
                           ax1_y >= ax1.limits[][2][1] &&
                           ax1_y <= ax1.limits[][2][2]
                            vmarker_pos = snap ? round(ax1_x; digits = 1) : ax1_x
                            if isnan(vmarker1[])
                                vmarker1[] = vmarker_pos
                            else
                                vmarker2[] = vmarker_pos
                            end
                            vmarker1[] > vmarker2[] &&
                                ((vmarker1[], vmarker2[]) = (vmarker2[], vmarker1[]))
                            vmarker1[] > t[end] && (vmarker1[] = t[end])
                            vmarker2[] > t[end] && (vmarker2[] = t[end])
                            marker_range[] = [vmarker1[], vmarker2[]]
                            notify(vmarker1);
                            notify(vmarker2);
                            notify(marker_range)
                        end
                    end

                    # change time window
                    if ax2_x >= 0 && ax2_y >= 0 && ax2_y <= 1
                        if ax2_x <= ax2.limits[][1][2] - seg_len
                            seg = (round(Int64, ax2_x), round(Int64, ax2_x) + seg_len)
                            ax1.limits[] = (seg, ax1.limits[][2])
                            seg_pos[] = round(Int64, ax2_x)
                        else
                            seg = (ceil(t[end]) - seg_len, ceil(t[end]))
                            ax1.limits[] = (seg, ax1.limits[][2])
                            seg_pos[] = seg[1]
                        end
                    end

                    # change channels window
                    if type === :normal
                        if ch_n > n_channels
                            if ax3_x >= 0 && ax3_x <= 1 && ax3_y >= 0 &&
                               ax3_y <= ax3.limits[][2][2]
                                ch1[] = floor(Int64, ax3_y)
                                ch1[] > ch_n - nch[] + 1 && (ch1[] = ch_n - nch[] + 1)
                                ax1.limits[] =
                                    (ax1.limits[][1], (ch1[] - 0.5, ch1[] + nch[] - 0.5))
                            end
                        end
                    end
                end
            end
        end

        # keyboard events
        on(events(fig).keyboardbutton) do event
            update_ax2 = false
            update_ax3 = false
            if event.action == Keyboard.press || event.action == Keyboard.repeat
                if type === :normal
                    if event.key == Keyboard.d
                        if !isnan(vmarker1[]) && !isnan(vmarker2[])
                            trim!(obj; seg = (marker_range[][1], marker_range[][2]))
                            screen = display(fig)
                            close(screen)
                            NeuroAnalyzer.plot(
                                obj;
                                ch         = ch,
                                seg        = (ax1.limits[][1][1], ax1.limits[][1][1] + seg_len),
                                xlabel     = xlabel,
                                ylabel     = ylabel,
                                title      = title,
                                markers    = markers,
                                scale      = scale,
                                group_ch   = group_ch,
                                n_channels = n_channels,
                                mono       = mono,
                                res        = res,
                                gui        = gui,
                            )
                        end
                    end

                    event.key == Keyboard.s && (snap = !snap)

                    if event.key == Keyboard.down
                        if ch1[] < ch_n - nch[] + 1
                            ch1[] += 1
                            update_ax3 = true
                        end
                    end
                    if event.key == Keyboard.up
                        if ch1[] > 1
                            ch1[] -= 1
                            update_ax3 = true
                        end
                    end
                    if ispressed(fig, Keyboard.page_down)
                        if ch_n > 1 && nch[] > 1
                            nch[] -= 1
                            update_ax3 = true
                        end
                    end
                    if ispressed(fig, Keyboard.page_up)
                        if ch_n > 1 && nch[] < ch_n && ch1[] + (nch[] - 1) < ch_n
                            nch[] += 1
                            update_ax3 = true
                        end
                    end
                end

                if event.key == Keyboard.home
                    seg_pos[] = 0
                    update_ax2 = true
                end
                if event.key == Keyboard._end
                    seg_pos[] = ceil(Int64, t[end] - seg_len)
                    update_ax2 = true
                end
                if event.key == Keyboard.left && seg_pos[] > 0
                    seg_pos[] -= 1
                    update_ax2 = true
                end
                if ispressed(fig, Keyboard.left_shift & Keyboard.left) && seg_pos[] >= 9
                    seg_pos[] -= 9
                    update_ax2 = true
                end
                if event.key == Keyboard.right && seg_pos[] < t[end] - seg_len
                    seg_pos[] += 1
                    update_ax2 = true
                end
                if ispressed(fig, Keyboard.left_shift & Keyboard.right) &&
                   seg_pos[] <= t[end] - seg_len - (seg_len - 1)
                    seg_pos[] += (seg_len - 1)
                    update_ax2 = true
                end

                update_ax2 &&
                    (ax1.limits[] = ((seg_pos[], seg_pos[] + seg_len), ax1.limits[][2]))
                update_ax3 &&
                    (ax1.limits[] = (ax1.limits[][1], (ch1[] - 0.5, ch1[] + nch[] - 0.5)))
            end
        end

        type === :normal && colsize!(fig.layout, 2, GLMakie.Fixed(20))
        rowsize!(fig.layout, 2, GLMakie.Fixed(20))

        wait(display(fig))
    end

    return fig
end
