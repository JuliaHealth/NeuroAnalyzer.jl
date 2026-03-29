# TO DO:

# types by colors
# add marker start : end
# select region
# select epoch
# time format (SS:MS HH:MM:SS)
# delete epoch
# change scaling
# plot(obj1, obj2)

export plot_ep

"""
    plot_ep(t, s; <keyword arguments>)

Plot epoched signal.

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
    - `:normal`
    - `:butterfly`: butterfly plot
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
    n_channels <= nchannels(obj) ||
        throw(ArgumentError("n_channels must be ≤ $(nchannels(obj))."))
    _check_var(type, [:normal, :butterfly], "type")
    avg && ci95 && throw(ArgumentError("avg and ci95 cannot both be true."))
    !_has_markers(obj) && (markers = false)

    # set color palette
    pal = mono ? :grays : :darktest

    _check_epochs(obj, ep)
    ep_len = epoch_duration(obj)
    ep_n = Observable(nepochs(obj))
    ep_n[] > 1 || throw(ArgumentError("Use plot_cont() for continuous object."))
    seg = (0, n_epochs * ep_len)
    epmarkers = [(idx - 1) * (epoch_len(obj) / sr(obj)) for idx in 1:ep_n[]]
    ep_selected = zeros(Bool, ep_n[])

    # check channels and meta data
    _ = get_channel(obj; ch = ch)
    obj_tmp = deepcopy(obj)
    keep_channel!(obj_tmp; ch = ch)
    obj_tmp.data = reshape(
        obj_tmp.data,
        size(obj_tmp.data, 1),
        size(obj_tmp.data, 2) * size(obj_tmp.data, 3),
        1,
    )

    ch_n = nchannels(obj_tmp)
    if group_ch
        ch_order = _sort_channels(obj_tmp.header.recording[:channel_type])
    else
        ch_order = collect(1:ch_n)
    end
    clabels = labels(obj_tmp)[ch_order]
    ctypes = obj_tmp.header.recording[:channel_type][ch_order]
    cunits = obj_tmp.header.recording[:unit][ch_order]

    # order by ctypes
    # and markers for ax3
    ctypes_uni = unique(ctypes)
    ctypes_pos = zeros(Int64, length(ctypes_uni))
    for idx in eachindex(ctypes_uni)
        ctypes_pos[idx] = findfirst(isequal(ctypes_uni[idx]), ctypes)
    end
    ctypes_uni_pos = zeros(Int64, ch_n)
    ctypes_uni_pos[ctypes_pos] .= 1

    t = obj_tmp.time_pts
    s = obj_tmp.data[ch_order, :, 1]

    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Epochs", "", "")

    # list of bad channels
    bad_ch = Observable(obj_tmp.header.recording[:bad_channel])

    # displayed segment
    seg_pos = Observable(seg[1])
    seg_len = (seg[2] - seg[1])

    if type === :normal
        nch      = Observable(n_channels)
        ch1      = Observable(1)
        ch2_init = gui && ch_n > nch[] ? ch1[] + nch[] - 1 : ch_n
    else
        ch1      = Observable(1)
        ch2_init = length(ctypes_uni)
        nch      = Observable(ch_n)
    end

    # get ranges of the original signal for the scales
    # normalize in groups by channel type
    # between -0.5 and +0.5 and shift so all channels are visible
    r = Observable(Float64[])
    for idx in eachindex(ctypes_uni)
        group = s[ctypes .== ctypes_uni[idx], :]
        push!(r[], round(_get_range(group)))
        # remove per-channel DC offset
        group = group .- mean(group; dims=2)
        # map to [-0.5, 0.5]
        s[ctypes .== ctypes_uni[idx], :] = normalize_minmax(group, 0.5; bych=true)
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
        ytc = repeat([:black], nchannels(obj_tmp))
        ytc[bad_ch[]] .= :lightgray
    else
        ytc = repeat([:black], ch_n)
    end

    # prepare markers
    if markers
        markers_pos = obj.markers[!, :start]
        markers_id = obj.markers[!, :id]
        markers_desc = obj.markers[!, :value]
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
        xlabel             = "",
        ylabel             = yl,
        title              = tt,
        xticks             = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        yticks             = (1:ch_n, clabels),
        # TO DO: yticklabelcolor=ytc[1:end],
        xautolimitmargin   = (0, 0),
        yautolimitmargin   = (0, 0),
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
                GLMakie.band!(ax1, t, s_u, s_l; alpha = 0.25, color = :grey, strokewidth = 0.5)
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
                    GLMakie.poly!(ax1, s_rectangle; color = :red, strokecolor = :red, strokewidth = 2)
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
                GLMakie.poly!(ax1, s_rectangle; color = :red, strokecolor = :red, strokewidth = 2)
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
                text           = "$(markers_id[idx]) / $(markers_desc[idx])",
                text_align     = (:left, :center),
                fontsize       = 8,
                cornerradius   = 0,
                cornervertices = 2,
                padding        = 2,
                strokewidth    = 1,
                offset         = (0, 5),
                text_rotation  = pi / 2,
            )
        end
    end

    if gui
        println()

        # time/epoch bar
        ax2 = GLMakie.Axis(
            fig[2, 1];
            xlabel             = xl,
            ylabel             = "",
            title              = "",
            xticks             = LinearTicks(25),
            yticksvisible      = false,
            xautolimitmargin   = (0, 0),
            yautolimitmargin   = (0, 0),
            backgroundcolor    = :white,
            _AXIS_LOCK_KWARGS...,
        )
        GLMakie.xlims!(ax2, 0, ep_n[])
        GLMakie.ylims!(ax2, 0, 1)
        hideydecorations!(ax2)
        hidexdecorations!(ax2; label = false, ticks = false, ticklabels = false)
        ax2.xticklabelsize = 12

        # epoch markers
        GLMakie.vlines!(ax2, 1:ep_n[]; linestyle = :dash, linewidth = 1, color = :black)

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
                ax3_x = mouseposition(ax3)[1]
                ax3_y = mouseposition(ax3)[2]
                if event.button == Mouse.right

                    # mark/unmark channel as bad
                    if type === :normal
                        if ax1_x < 0
                            bad_ch[][round(Int64, ax1_y)] = !bad_ch[][round(Int64, ax1_y)]
                            obj.header.recording[:bad_channel][
                                get_channel(obj; ch = clabels[round(Int64, ax1_y)])[1],
                            ] = !obj.header.recording[:bad_channel][
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
                    if ax1_y >= ax1.limits[][2][1] && ax1_y <= ax1.limits[][2][2]
                        nep = ceil(Int64, ax1_x / ep_len)
                        if 1 <= nep <= ep_n[]
                            ep_selected[nep] = !ep_selected[nep]
                        end
                    end

                    # change displayed epoch window
                    if ax2_x >= 0 && ax2_x <= ax2.limits[][1][2] && ax2_y >= 0 && ax2_y <= 1
                        nep = clamp(round(Int64, ax2_x), 1, ep_n[])
                        seg = ((nep - 1) * ep_len, (nep + n_epochs - 1) * ep_len)
                        ax1.limits[] = (seg, ax1.limits[][2])
                        seg_pos[] = Float64(nep - 1)
                    end

                    # change channels window
                    if type === :normal
                        ax3_x = mouseposition(ax3)[1]
                        ax3_y = mouseposition(ax3)[2]
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
                if event.key == Keyboard.home
                    seg_pos[] = 0.0
                    update_ax2 = true
                end
                if event.key == Keyboard._end
                    seg_pos[] = Float64(ep_n[] - n_epochs)
                    update_ax2 = true
                end
                if event.key == Keyboard.left && seg_pos[] > 0
                    seg_pos[] -= 1.0
                    update_ax2 = true
                end
                if ispressed(fig, Keyboard.left_shift & Keyboard.left)
                    seg_pos[] = clamp(seg_pos[] - (n_epochs - 1), 0.0, Float64(ep_n[] - n_epochs))
                    update_ax2 = true
                end
                if event.key == Keyboard.right && seg_pos[] < ep_n[] - n_epochs
                    seg_pos[] += 1.0
                    update_ax2 = true
                end
                if ispressed(fig, Keyboard.left_shift & Keyboard.right)
                    seg_pos[] = clamp(seg_pos[] + (n_epochs - 1), 0.0, Float64(ep_n[] - n_epochs))
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

        wait(display(fig))
    end

    return fig
end
