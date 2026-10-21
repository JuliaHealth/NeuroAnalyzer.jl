# TO DO:

# types by colors
# time format (SS:MS HH:MM:SS)
# change scaling

export plot_ep

"""
    plot_ep(t, s; <keyword arguments>)

Plot epoched signal with interactive editing.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name or list of channel names
- `ep::Int64=1`: first epoch to plot
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `markers::Bool`: draw markers if available
- `scale::Bool=true`: draw scale
- `group_ch::Bool=true`: group channels by type
- `type::Symbol=:normal`: plot type:
    - `:normal`: standard multi-channel plot
    - `:butterfly`: butterfly plot showing all channels overlaid
- `avg::Bool=false`: plot averaged channel in butterfly plot
- `ci95::Bool=false`: plot averaged channels and 95% CI in butterfly plot
- `n_channels::Int64=20`: number of visible channels
- `n_epochs::Int64=5`: number of visible epochs
- `res::Int64=1`: resampling factor (draw every `res`-nth sample)
- `gui::Bool=true`: if `true`, keep window open and interactive

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_ep(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = "all",
    ep::Int64 = 1,
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
    n_epochs::Int64 = 5,
    res::Int64 = 1,
    gui::Bool = true,
)::GLMakie.Figure
    # validate
    res >= 1 || throw(ArgumentError("res must be ≥ 1."))
    res > 10 && _warn("At res > 10 plot will be inaccurate.")
    n_channels >= 1 || throw(ArgumentError("n_channels must be ≥ 1."))
    # at most all channels are shown at once
    n_channels = min(n_channels, nchannels(obj))
    _check_var(type, [:normal, :butterfly], "type")
    avg && ci95 && throw(ArgumentError("avg and ci95 cannot both be true."))
    !_has_markers(obj) && (markers = false)

    # set color palette
    pal = mono ? :grays : :darktest

    _check_epochs(obj, ep)
    ep_len = epoch_duration(obj)
    ep_n = Observable(nepochs(obj))
    ep_n[] > 1 || throw(ArgumentError("Use plot_cont() for continuous object."))
    # displayed segment starts at epoch `ep`; the last valid position (in epochs) leaves
    # no empty space after the last epoch
    seg_pos_max() = Float64(max(ep_n[] - n_epochs, 0))
    seg_pos = Observable(clamp(Float64(ep - 1), 0.0, seg_pos_max()))
    seg = (seg_pos[] * ep_len, (seg_pos[] + n_epochs) * ep_len)
    epmarkers = @lift([(idx - 1) * ep_len for idx = 1:$ep_n])
    ep_selected = Observable(zeros(Bool, ep_n[]))

    # check channels and meta data
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    ch_n = length(ch)
    if group_ch
        ch_order = _sort_channels(obj.header.recording[:channel_type][ch])
    else
        ch_order = collect(1:ch_n)
    end
    clabels = labels(obj)[ch][ch_order]
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

    # get time points vector (resampled per epoch, matching the signal matrix)
    t = reshape(obj.time_pts, :, ep_n[])[1:res:end, :][:]
    # get signal matrix
    s = obj.data[ch, :, :][ch_order, 1:res:end, :]
    s = reshape(s, :, size(s, 2) * size(s, 3), 1)
    s = dropdims(s; dims = 3)

    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Epochs", "", "")

    # list of bad channels
    bad_ch = Observable(obj.header.recording[:bad_channel])

    # displayed segment length
    seg_len = (seg[2] - seg[1])

    if type === :normal
        nch = Observable(n_channels)
        ch1 = Observable(1)
        ch2_init = ch_n > nch[] ? ch1[] + nch[] - 1 : ch_n
        clabels = labels(obj)[ch][ch_order]
    else
        ch_n = length(ctypes_uni)
        ch1 = Observable(1)
        ch2_init = length(ctypes_uni)
        nch = Observable(ch_n)
        clabels = uppercase.(ctypes_uni)
    end

    # common amplitude scale of each channel type (one row = r signal units), robust to artifacts;
    # each channel is centered at its median
    r = Observable(Float64[])
    for idx in eachindex(ctypes_uni)
        g_idx = ctypes .== ctypes_uni[idx]
        group = s[g_idx, :]
        # butterfly: all channels of the type share one row, ±6 SD fit in it
        sc = _plot_scale(group) * (type === :butterfly ? 2 : 1)
        push!(r[], round(sc; sigdigits = 3))
        s[g_idx, :] = (group .- median(group; dims = 2)) ./ sc
    end
    if type === :normal
        s .+= collect(1:ch_n)
    elseif type === :butterfly
        for idx in eachindex(ctypes_uni)
            s[ctypes .== ctypes_uni[idx], :] .+= idx
        end
    end

    # signal is observable so that deleted epochs can be removed from the plot
    ep_cols = size(s, 2) ÷ ep_n[]
    s = Observable(s)
    # time points matching the (possibly shortened) signal
    _pts(v) = Point2.(view(t, 1:length(v)), v)

    # y-axis labels colors
    if type === :normal
        ytc = repeat([:black], nchannels(obj))
        ytc[bad_ch[]] .= :lightgray
    else
        ytc = repeat([:black], ch_n)
    end

    # prepare markers
    if markers
        markers_pos = Observable(obj.markers[!, :start])
    end

    # prepare plot
    if gui
        if type === :normal
            plot_size = (1250, 700)
        else
            plot_size = (1200, 700)
        end
    else
        plot_size = (1200, 650)
    end
    GLMakie.activate!(; title = "plot_ep()")
    fig = GLMakie.Figure(;
        size = plot_size,
        figure_padding = (10, 20, 10, 10), # L R B T
    )

    # create axis with customizable properties
    ax1 = GLMakie.Axis(
        fig[1, 1];
        xlabel = "",
        ylabel = yl,
        title = tt,
        _TIME_AXIS_KWARGS...,
        yticks = (1:ch_n, clabels),
        # TO DO: yticklabelcolor=ytc[1:end],
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
        yticklabelspace = let ml = maximum(length, clabels)
            ml <= 5 ? 60.0 : ml >= 10 ? 100.0 : 80.0
        end,
    )
    GLMakie.xlims!(ax1, seg)
    if ch_n > nch[]
        GLMakie.ylims!(ax1, ch2_init + 0.5, ch2_init - nch[] + 0.5)
    else
        GLMakie.ylims!(ax1, ch_n + 0.5, 0.5)
    end
    _style_axis!(ax1)

    # highlight selected epochs (drawn first so that signal stays on top)
    ep_sel_color = mono ? (:black, 0.15) : (:orange, 0.25)
    ep_sel_rect1 = lift(ep_selected) do v
        return Rect2d[Rect2d((idx - 1) * ep_len, 0, ep_len, ch_n + 1) for idx in findall(v)]
    end
    GLMakie.poly!(ax1, ep_sel_rect1; color = ep_sel_color)

    # draw channels
    if type === :normal
        for idx = 1:ch_n
            line_color = @lift($bad_ch[idx] ? :lightgray : :black)
            GLMakie.lines!(ax1, @lift(_pts($s[idx, :])); linewidth = 1.5, color = line_color)
        end

    elseif type === :butterfly
        if ci95
            for idx in eachindex(ctypes_uni)
                msci95_data = @lift(NeuroAnalyzer._msci95($s[ctypes .== ctypes_uni[idx], :]))
                GLMakie.band!(
                    ax1, @lift(_pts($msci95_data.ul)), @lift(_pts($msci95_data.ll));
                    alpha = 0.25,
                    color = :grey,
                    strokewidth = 0.5,
                )
                GLMakie.lines!(
                    ax1, @lift(_pts($msci95_data.sm)); color = :black, linewidth = 2,
                )
            end

        else
            !mono && (cmap = GLMakie.resample_cmap(pal, max(2, size(s[], 1))))
            for idx in axes(s[], 1)
                GLMakie.lines!(
                    ax1, @lift(_pts($s[idx, :]));
                    color = mono ? :black : cmap[idx],
                    colormap = pal,
                    linewidth = 0.5,
                    alpha = 1.0,
                )
            end
            if avg
                for idx in eachindex(ctypes_uni)
                    s_avg = @lift(_pts(mean($s[ctypes .== ctypes_uni[idx], :]; dims = 1)[:]))
                    GLMakie.lines!(ax1, s_avg; linewidth = 2, color = :black)
                end
            end
        end
    end

    # draw epochs markers
    # TO DO: draw epoch numbers
    GLMakie.vlines!(
        ax1,
        epmarkers;
        linestyle = :dot,
        linewidth = 0.5,
        color = mono ? :black : :blue,
    )

    # draw scale bars
    # TO DO: place scale values on the left side, below channel label
    if scale
        # scale bar width as a fixed fraction of the displayed segment, so its
        # on-screen width does not depend on the segment length
        # (equals the former 0.01 s for a 10 s segment)
        sb_w = seg_len / 1000
        # seg_pos is in epochs, ax1 x-axis is in seconds
        if type === :normal
            idx2 = 1
            for idx1 = 1:ch_n
                if ctypes_uni_pos[idx1] == 1
                    s_rectangle = lift(seg_pos) do sp
                        return Rect(sp * ep_len, (idx1 - 0.49), sb_w, 0.98)
                    end
                    l_pos = lift(seg_pos) do sp
                        return (sp * ep_len + sb_w, idx1 + 0.49)
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
                        text = string(r[][idx2]) * " " * cunits[idx1],
                        fontsize = 10,
                        color = :red,
                        align = (:left, :bottom),
                        offset = (5, 0),
                    )
                    idx2 += 1
                end
            end
        elseif type === :butterfly
            for idx in eachindex(ctypes_uni)
                s_rectangle = lift(seg_pos) do sp
                    return Rect(sp * ep_len, (idx - 0.475), sb_w, 0.975)
                end
                l_pos = lift(seg_pos) do sp
                    return (sp * ep_len, idx + 0.5)
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
                    text = string(r[][idx]) * " " * cunits[ctypes .== ctypes_uni[idx]][1],
                    markerspace = :pixel,
                    fontsize = 10,
                    color = :red,
                    align = (:left, :bottom),
                    offset = (5, 0),
                )
            end
        end
    end

    # draw event markers
    if markers
        GLMakie.vlines!(ax1, markers_pos; linestyle = :dash, linewidth = 1, color = :black)
        markers_ypos = lift(markers_pos, ch1, nch) do mp, v1, v2
            return [Point2(p, v1 + (v2 - 1) + 0.5) for p in mp]
        end
        # labels are redrawn (not updated) after epoch deletion:
        # textlabel!() leaves stale backgrounds when its positions become empty
        draw_marker_labels() = isempty(obj.markers) ? nothing : GLMakie.textlabel!(
            ax1, markers_ypos;
            text = string.(obj.markers[!, :id], " / ", obj.markers[!, :value]),
            text_align = (:left, :center),
            fontsize = 8,
            cornerradius = 0,
            cornervertices = 2,
            padding = 2,
            strokewidth = 1,
            offset = (0, 5),
            text_rotation = pi / 2,
        )
        markers_labels = draw_marker_labels()
    end

    if gui
        println()

        # time/epoch bar
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
        GLMakie.xlims!(ax2, 0, ep_n[])
        GLMakie.ylims!(ax2, 0, 1)
        hideydecorations!(ax2)
        hidexdecorations!(ax2; label = false, ticks = false, ticklabels = false)
        ax2.xticklabelsize = 12

        # selected epochs
        ep_sel_rect2 = lift(ep_selected) do v
            return Rect2d[Rect2d(idx - 1, 0, 1, 1) for idx in findall(v)]
        end
        GLMakie.poly!(ax2, ep_sel_rect2; color = ep_sel_color)

        # epoch markers
        GLMakie.vlines!(
            ax2, @lift(collect(1:$ep_n)); linestyle = :dash, linewidth = 1, color = :black,
        )

        # time line marker
        # define a square: Rect(x, y, width, height)
        t_rectangle = lift(seg_pos) do v
            return Rect(v, 0, n_epochs, 1)
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
        if type === :normal
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
                backgroundcolor = :white,
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
        end

        # mouse events
        on(events(fig).mousebutton) do event
            if event.action == Mouse.press
                ax1_x = mouseposition(ax1)[1]
                ax1_y = mouseposition(ax1)[2]
                ax2_x = mouseposition(ax2)[1]
                ax2_y = mouseposition(ax2)[2]
                if type === :normal
                    ax3_x = mouseposition(ax3)[1]
                    ax3_y = mouseposition(ax3)[2]
                end

                if event.button == Mouse.right

                    # mark/unmark channel as bad
                    if type === :normal
                        if ax1_x < 0
                            bad_ch[][round(Int64, ax1_y)] = !bad_ch[][round(Int64, ax1_y)]
                            obj.header.recording[:bad_channel][
                                get_channel(obj; ch = clabels[round(Int64, ax1_y)])[1],
                            ] =
                                !obj.header.recording[:bad_channel][
                                    get_channel(obj; ch = clabels[round(Int64, ax1_y)])[1],
                                ]
                            notify(bad_ch)
                        end
                    end

                elseif event.button == Mouse.left

                    # get channel info
                    if ax1_x < 0
                        channel_info(obj; ch = clabels[round(Int64, ax1_y)])
                    end

                    # select / deselect epochs
                    if is_mouseinside(ax1)
                        nep = ceil(Int64, ax1_x / ep_len)
                        if 1 <= nep <= ep_n[]
                            ep_selected[][nep] = !ep_selected[][nep]
                            notify(ep_selected)
                        end
                    end

                    # change displayed epoch window
                    if ax2_x >= 0 && ax2_x <= ax2.limits[][1][2] && ax2_y >= 0 && ax2_y <= 1
                        seg_pos[] = clamp(Float64(round(Int64, ax2_x) - 1), 0.0, seg_pos_max())
                        seg = (seg_pos[] * ep_len, (seg_pos[] + n_epochs) * ep_len)
                        ax1.limits[] = (seg, ax1.limits[][2])
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
                # delete selected epochs
                if event.key == Keyboard.d && event.action == Keyboard.press &&
                   any(ep_selected[])
                    ep_del = findall(ep_selected[])
                    if length(ep_del) == ep_n[]
                        _warn("Cannot delete all epochs.")
                    else
                        delete_epoch!(obj; ep = ep_del)
                        ep_keep = setdiff(1:ep_n[], ep_del)
                        cols = reduce(vcat, [((e - 1) * ep_cols + 1):(e * ep_cols) for e in ep_keep])
                        ep_selected[] = zeros(Bool, length(ep_keep))
                        ep_n[] = length(ep_keep)
                        s[] = s[][:, cols]
                        if markers
                            # delete labels before updating positions, so they never receive an empty vector
                            !isnothing(markers_labels) && delete!(ax1, markers_labels)
                            markers_pos[] = obj.markers[!, :start]
                            markers_labels = draw_marker_labels()
                        end
                        GLMakie.xlims!(ax2, 0, ep_n[])
                        seg_pos[] = clamp(seg_pos[], 0.0, seg_pos_max())
                        update_ax2 = true
                    end
                end
                if event.key == Keyboard.home
                    seg_pos[] = 0.0
                    update_ax2 = true
                end
                if event.key == Keyboard._end
                    seg_pos[] = seg_pos_max()
                    update_ax2 = true
                end
                if event.key == Keyboard.left && seg_pos[] > 0
                    seg_pos[] -= 1.0
                    update_ax2 = true
                end
                if ispressed(fig, Keyboard.left_shift & Keyboard.left)
                    seg_pos[] =
                        clamp(seg_pos[] - (n_epochs - 1), 0.0, seg_pos_max())
                    update_ax2 = true
                end
                if event.key == Keyboard.right && seg_pos[] < seg_pos_max()
                    seg_pos[] += 1.0
                    update_ax2 = true
                end
                if ispressed(fig, Keyboard.left_shift & Keyboard.right)
                    seg_pos[] =
                        clamp(seg_pos[] + (n_epochs - 1), 0.0, seg_pos_max())
                    update_ax2 = true
                end
                if update_ax2
                    seg = (seg_pos[] * ep_len, (seg_pos[] + n_epochs) * ep_len)
                    ax1.limits[] = (seg, ax1.limits[][2])
                end
                if update_ax3
                    ax1.limits[] = (ax1.limits[][1], (ch1[] - 0.5, ch1[] + nch[] - 0.5))
                end
            end
        end

        type === :normal && colsize!(fig.layout, 2, GLMakie.Fixed(20))
        rowsize!(fig.layout, 2, GLMakie.Fixed(20))

        screen = display(fig)
        GLMakie.GLFW.MaximizeWindow(screen.glscreen)
        wait(screen)
    end

    return fig
end
