export plot_cont

"""
    plot_cont(obj1, obj2; <keyword arguments>)

Plot two continuous signals overlayed for comparison.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object, will be drawn in blue
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object, will be drawn in black
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name(s)
- `seg::Tuple{Real, Real}=(0, 10)`: time segment to display in seconds (from, to), default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `scale::Bool=true`: if `true`, draw scale reference
- `group_ch::Bool=true`: if `true`, group channels by type (e.g. EEG, EOG, ECG)
- `n_channels::Int64=20`: maximum number of visible channels
- `res::Int64=1`: resampling factor (draw every `res`-nth sample)
- `gui::Bool=true`: if `true`, keep window open and interactive

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_cont(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = "all",
    seg::Tuple{Real, Real} = (0, 10),
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    scale::Bool = true,
    group_ch::Bool = true,
    n_channels::Int64 = 20,
    res::Int64 = 1,
    gui::Bool = true,
)::GLMakie.Figure
    # validate
    res >= 1 || throw(ArgumentError("res must be ≥ 1."))
    res > 10 && _warn("At res > 10 plot will be inaccurate.")
    n_channels >= 1 || throw(ArgumentError("n_channels must be ≥ 1."))
    n_channels <= nchannels(obj1) ||
        throw(ArgumentError("n_channels must be ≤ $(nchannels(obj1))."))
    size(obj1) == size(obj2) ||
        throw(ArgumentError("Size of OBJ1 and OBJ2 must equal."))
    sr(obj1) == sr(obj2) ||
        throw(ArgumentError("Sampling rate of OBJ1 and OBJ2 must equal."))
    labels(obj1) == labels(obj2) ||
        throw(ArgumentError("Labels of OBJ1 and OBJ2 must equal."))
    obj1.header.recording[:channel_type] == obj2.header.recording[:channel_type] ||
        throw(ArgumentError("Channel types of OBJ1 and OBJ2 must equal."))
    obj1.header.recording[:unit] == obj2.header.recording[:unit] ||
        throw(ArgumentError("Channel units of OBJ1 and OBJ2 must equal."))

    # set maximum length to display
    if signal_len(obj1) <= seg[2] * sr(obj1)
        seg = (obj1.time_pts[1], obj1.time_pts[end])
    else
        _check_segment(obj1, seg)
    end

    # check channels and meta data
    ch = get_channel(obj1; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    ch_n = length(ch)
    if group_ch
        ch_order = _sort_channels(obj1.header.recording[:channel_type][ch])
    else
        ch_order = collect(1:ch_n)
    end
    clabels = labels(obj1)[ch][ch_order]
    ctypes  = obj1.header.recording[:channel_type][ch][ch_order]
    cunits  = obj1.header.recording[:unit][ch][ch_order]

    # order by ctypes
    # and markers for ax3
    ctypes_uni = unique(ctypes)
    ctypes_pos = zeros(Int64, length(ctypes_uni))
    for idx in eachindex(ctypes_uni)
        ctypes_pos[idx] = findfirst(isequal(ctypes_uni[idx]), ctypes)
    end
    ctypes_uni_pos = zeros(Int64, ch_n)
    ctypes_uni_pos[ctypes_pos] .= 1

    # get time points vector
    t = obj1.time_pts
    # get signal matrices
    s1 = obj1.data[ch, :, 1][ch_order, 1:res:end]
    s2 = obj2.data[ch, :, 1][ch_order, 1:res:end]

    # set defaults
    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "")

    # displayed segment
    seg_pos = Observable(seg[1])
    seg_len = (seg[2] - seg[1])

    nch      = Observable(n_channels)
    ch1      = Observable(1)
    ch2_init = gui && ch_n > nch[] ? ch1[] + nch[] - 1 : ch_n

    # get ranges of the original signal for the scales
    # normalize in groups by channel type
    # between -0.5 and +0.5 and shift so all channels are visible
    r = Observable(Float64[])
    for idx in eachindex(ctypes_uni)
        group1 = s1[ctypes .== ctypes_uni[idx], :]
        group2 = s2[ctypes .== ctypes_uni[idx], :]
        push!(r[], round(_get_range(group1)))
        # remove per-channel DC offset
        group1 = group1 .- mean(group1; dims = 2)
        group2 = group2 .- mean(group2; dims = 2)
        # map to [-0.5, 0.5]
        s1[ctypes .== ctypes_uni[idx], :] = normalize_minmax(group1, 0.5; bych = true)
        s2[ctypes .== ctypes_uni[idx], :] = normalize_minmax(group2, 0.5; bych = true)
    end
    s1 .+= collect(1:ch_n)
    s2 .+= collect(1:ch_n)

    # prepare plot
    if gui
        plot_size = (1250, 700)
    else
        plot_size = (1200, 650)
    end
    GLMakie.activate!(; title = "plot()")

    # create axis with customizable properties
    fig = GLMakie.Figure(;
        size = plot_size,
        figure_padding = (10, 20, 10, 10), # L R B T
    )
    ax1 = GLMakie.Axis(
        fig[1, 1];
        xlabel = "",
        ylabel = yl,
        title = tt,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        yticks = (1:ch_n, clabels),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
        yticklabelspace = let ml = maximum(length, clabels)
            ml <= 5 ? 60.0 : ml >= 10 ? 100.0 : 80.0
        end,
    )
    GLMakie.xlims!(ax1, seg)
    if gui
        if ch_n > nch[]
            GLMakie.ylims!(ax1, ch2_init + 0.5, ch2_init - nch[] + 0.5)
        else
            GLMakie.ylims!(ax1, ch_n + 0.5, 0.5)
        end
    else
        GLMakie.ylims!(ax1, ch_n + 0.5, 0.5)
    end
    _style_axis!(ax1)

    # draw channels
    for idx in 1:ch_n
        GLMakie.lines!(
            ax1,
            t[1:res:end],
            s1[idx, 1:res:end];
            linewidth = 1.5,
            alpha = 0.4,
            color = :blue,
        )
        GLMakie.lines!(
            ax1,
            t[1:res:end],
            s2[idx, 1:res:end];
            linewidth = 1.5,
            color = :black,
        )
    end

    # draw scale bars
    # TO DO: place scale values on the left side, below channel label
    if scale
        idx2 = 1
        for idx1 in 1:ch_n
            if ctypes_uni_pos[idx1] == 1
                s_rectangle = lift(seg_pos) do seg_pos
                    return Rect(seg_pos, (idx1 - 0.49), 0.01, 0.98)
                end
                l_pos = lift(seg_pos) do seg_pos
                    return (seg_pos + 0.01, idx1 + 0.49)
                end
                GLMakie.poly!(
                    ax1, s_rectangle;
                    color = :red,
                    strokecolor = :red,
                    strokewidth = 2,
                )
                GLMakie.text!(
                    ax1, l_pos;
                    markerspace = :pixel,
                    text = string(r[][idx2]) * " " * cunits[idx1],
                    fontsize = 10,
                    color = :red,
                    align = (:left, :bottom),
                    offset = (5, 0),
                )
                idx2 += 1
            end
        end
    end

    if gui

        # time bar
        ax2 = GLMakie.Axis(
            fig[2, 1];
            xlabel = xl,
            ylabel = "",
            title = "",
            xticks = LinearTicks(25),
            yticksvisible = false,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            backgroundcolor = :white,
            _AXIS_LOCK_KWARGS...,
        )
        GLMakie.xlims!(ax2, t[1], t[end])
        GLMakie.ylims!(ax2, 0, 1)
        hideydecorations!(ax2)
        hidexdecorations!(ax2; label = false, ticks = false, ticklabels = false)
        ax2.xticklabelsize = 12

        # time line marker
        # define a square: Rect(x, y, width, height)
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

        # channel bar
        ax3 = GLMakie.Axis(
            fig[1, 2];
            xlabel = "",
            ylabel = "",
            title = "",
            yticks = 1:ch_n,
            xticksvisible = false,
            yticksvisible = false,
            yreversed = true,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            _AXIS_LOCK_KWARGS...,
        )
        ch_n > 1 && (GLMakie.ylims!(ax3, ch_n, 1))
        hidedecorations!(ax3)

        # mark channel types
        if group_ch
            for idx in eachindex(ctypes_pos)
                GLMakie.hlines!(ax3, ctypes_pos[idx]; linewidth = 5, color = :black)
            end
        end

        # channel marker
        # define a square: Rect(x, y, width, height)
        ch_rectangle = @lift(Rect(0, $ch1, 1, $nch - 1))
        GLMakie.poly!(
            ax3,
            ch_rectangle;
            color = :darkgrey,
            strokecolor = :black,
            strokewidth = 2,
            alpha = 0.25,
        )

        # mouse events
        on(events(fig).mousebutton) do event
            ax1_x = mouseposition(ax1)[1]
            ax1_y = mouseposition(ax1)[2]
            ax2_x = mouseposition(ax2)[1]
            ax2_y = mouseposition(ax2)[2]
            ax3_x = mouseposition(ax3)[1]
            ax3_y = mouseposition(ax3)[2]

            if event.action == Mouse.press
                if event.button == Mouse.left

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

        # keyboard events
        on(events(fig).keyboardbutton) do event
            update_ax2 = false
            update_ax3 = false
            if event.action == Keyboard.press || event.action == Keyboard.repeat
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

        colsize!(fig.layout, 2, GLMakie.Fixed(20))
        rowsize!(fig.layout, 2, GLMakie.Fixed(20))

        wait(display(fig))
    end

    return fig
end
