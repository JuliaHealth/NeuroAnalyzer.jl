# TO DO:

# types by colors
# add marker start : end
# select region
# select epoch
# time format (SS:MS HH:MM:SS)
# delete region
# delete epoch
# click channel to mark as bad
# click label to get info / position
# change scaling
# plot(obj1, obj2)

export plot_cont

"""
    plot_cont(obj; <keyword arguments>)

Plot continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name or list of channel names
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title
- `mono::Bool=false`: use color or gray palette
- `markers::Bool`: draw markers if available
- `scale::Bool=true`: draw scale
- `group_ch::Bool=true`: group channels by type
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:butterfly`: butterfly plot
- `avg::Bool=false`: plot averaged channel in butterfly plot
- `ci95::Bool=false`: plot averaged channels and 95% CI in butterfly plot
- `n_channels::Int64=20`: number of visible channels
- `res::Int64=1`: resampling factor (draw every res-nth sample)
- `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

- `p::GLMakie.Figure`
"""
function plot_cont(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}="all", seg::Tuple{Real, Real}=(0, 10), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, scale::Bool=true, group_ch::Bool=true, type::Symbol=:normal, avg::Bool=true, ci95::Bool=false, n_channels::Int64=20, res::Int64=1, gui::Bool=true)::GLMakie.Figure

    @assert res >= 1 "res must be ≥ 1."
    res > 10 && _warn("At res > 10 plot will be inaccurate.")
    @assert n_channels >= 1 "n_channels must be ≥ 1."
    @assert n_channels <= nchannels(obj) "n_channels must be ≤ $(nchannels(obj))."
    _check_var(type, [:normal, :butterfly], "type")
    !_has_markers(obj) && (markers = false)

    pal = mono ? :grays : :darktest

    if signal_len(obj) <= seg[2] * sr(obj)
        seg = (obj.time_pts[1], obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end

    # check channels and meta data
    _ = get_channel(obj, ch=ch)
    obj_tmp = deepcopy(obj)
    keep_channel!(obj_tmp, ch=ch)
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

    t = Observable(obj_tmp.time_pts)
    s = Observable(obj_tmp.data[ch_order, :, 1])

    xl, yl, tt = _set_defaults(xlabel,
                               ylabel,
                               title,
                               "Time [s]",
                               "",
                               "")

    # list of bad channels
    bad_ch = Observable(obj_tmp.header.recording[:bad_channel])

    # displayed segment
    seg_pos = Observable(seg[1])
    seg_len = (seg[2] - seg[1])

    if type === :normal
        nch = Observable(n_channels)
        if gui
            if ch_n > nch[]
                ch1 = Observable(1)
                ch2 = ch1[] + nch[] - 1
            else
                ch1 = Observable(1)
                ch2 = ch_n
            end
        else
            ch1 = Observable(1)
            ch2 = ch_n
        end
    else
        ch1 = Observable(1)
        ch2 = length(ctypes_uni)
        clabels = uppercase.(ctypes_uni)
        ch_n = length(ctypes_uni)
        nch = Observable(ch_n)
    end
    # get ranges of the original signal for the scales
    # normalize in groups by channel type
    # between -1.0 and +1.0 and shift so all channels are visible
    r = Observable(Float64[])
    for idx in eachindex(ctypes_uni)
        push!(r[], round(_get_range(s[][ctypes .== ctypes_uni[idx], :])))
        s[][ctypes .== ctypes_uni[idx], :] = normalize_minmax(s[][ctypes .== ctypes_uni[idx], :])
    end
    if type === :normal
        s[] .+= collect(1:ch_n)
    else
        for idx in eachindex(ctypes_uni)
            s[][ctypes .== ctypes_uni[idx], :] .+= idx
        end
    end

    # y-axis labels colors
    if type === :normal
        ytc = repeat([:black], nchannels(obj_tmp))
        ytc[bad_ch[]] .= :lightgray
    else
        ytc = repeat([:black], ch_n)
    end
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
    p = GLMakie.Figure(size=plot_size,
                       figure_padding=(10, 20, 10, 10)) # L R B T)
    ax1 = GLMakie.Axis(p[1, 1],
                       xlabel="",
                       ylabel=yl,
                       title=tt,
                       xticks=LinearTicks(10),
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(10),
                       yticks=(1:ch_n, clabels),
                       # TO DO: yticklabelcolor=ytc[1:end],
                       xautolimitmargin=(0, 0),
                       yautolimitmargin=(0, 0),
                       xzoomlock=true,
                       yzoomlock=true,
                       xpanlock=true,
                       ypanlock=true,
                       xrectzoom=false,
                       yrectzoom=false,
                       yticklabelspace=60.0)
    GLMakie.xlims!(ax1, seg)
    if gui
        if ch_n > nch[]
            GLMakie.ylims!(ax1, ch2 + 0.5, ch2 - nch[] + 0.5)
        else
            GLMakie.ylims!(ax1, ch_n + 0.5, 0.5)
        end
    else
        GLMakie.ylims!(ax1, ch_n + 0.5, 0.5)
    end
    ax1.titlesize = 20
    ax1.xlabelsize = 12
    ax1.ylabelsize = 12
    ax1.xticklabelsize = 12
    ax1.yticklabelsize = 12

    # draw channels
    if type === :normal
        for idx in 1:ch_n
            GLMakie.lines!(ax1,
                           t[][1:res:end],
                           @lift($s[idx, 1:res:end]),
                           linewidth=1.5,
                           color=bad_ch[][idx] ? :lightgray : :black)
        end
    else
        if ci95
            for idx in eachindex(ctypes_uni)
                s_m, _, s_u, s_l = NeuroAnalyzer.msci95(s[][ctypes .== ctypes_uni[idx], :])
                # draw 95% CI
                Makie.band!(ax1,
                            t[][1:res:end],
                            s_u[1:res:end],
                            s_l[1:res:end],
                            alpha=0.25,
                            color=:grey,
                            strokewidth=0.5)
                # draw mean
                Makie.lines!(ax1,
                             t[][1:res:end],
                             s_m[1:res:end],
                             color=:black,
                             linewidth=2)
            end
        else
            !mono && (cmap = GLMakie.resample_cmap(pal, size(s[], 1)))
            for idx in axes(s[], 1)
                GLMakie.lines!(ax1,
                               t[][1:res:end],
                               @lift($s[idx, 1:res:end]),
                               color=mono ? :black : cmap[idx],
                               colormap=pal,
                               colorrange=1:size(s[], 1),
                               linewidth=0.5)
            end

            # plot averaged channels
            if avg
                for idx in eachindex(ctypes_uni)
                    s_avg = mean(s[][ctypes .== ctypes_uni[idx], :], dims=1)[:]
                    GLMakie.lines!(ax1,
                                   t[][1:res:end],
                                   s_avg[1:res:end],
                                   linewidth=2,
                                   color=:black)
                end
            end
        end
    end
    # draw epochs markers
    # TODO: draw epoch numbers
#    if !isnothing(emarkers)
#        GLMakie.vlines!(ax1,
#                        emarkers,
#                        linestyle=:dot,
#                        linewidth=0.5,
#                        color=:blue)
#    end

    # draw scale bars
    # TO DO: place scale values on the left side, below channel label
    if scale
        if type === :normal
            idx2 = 1
            for idx1 in 1:ch_n
                if ctypes_uni_pos[idx1] == 1
                    s_rectangle = lift(seg_pos) do seg_pos
                        Rect(seg_pos, (idx1 - 0.49), 0.01, 0.98)
                    end
                    l_pos = lift(seg_pos) do seg_pos
                        (seg_pos + 0.01, idx1 + 0.49)
                    end
                    GLMakie.poly!(ax1,
                                  s_rectangle,
                                  color=:red,
                                  strokecolor=:red,
                                  strokewidth=2)
                    GLMakie.text!(ax1,
                                  l_pos,
                                  markerspace=:pixel,
                                  text=string(r[][idx2]) * " " * cunits[idx1],
                                  fontsize=10,
                                  color=:red,
                                  align=(:left, :bottom),
                                  #rotation=pi/2,
                                  offset=(5, 0))
                    idx2 += 1
                end
            end
        else
            for idx in 1:ch_n
                s_rectangle = lift(seg_pos) do seg_pos
                    Rect(seg_pos, (idx - 0.475), 0.01, 0.975)
                end
                l_pos = lift(seg_pos) do seg_pos
                    (seg_pos, idx + 0.5)
                end
                GLMakie.poly!(ax1,
                              s_rectangle,
                              color=:red,
                              strokecolor=:red,
                              strokewidth=2)
                GLMakie.text!(ax1,
                              l_pos,
                              text=string(r[][idx]) * " " * cunits[ctypes .== ctypes_uni[idx]][1],
                              markerspace=:pixel,
                              fontsize=10,
                              color=:red,
                              align=(:left, :bottom),
                              #rotation=pi/2,
                              offset=(5, 0))
            end
        end
    end

    # plot markers if available
    if markers
        GLMakie.vlines!(ax1,
                        markers_pos,
                        linestyle=:dash,
                        linewidth=1,
                        color=:black)
        for idx in eachindex(markers_pos)
            markers_ypos = lift(ch1, nch) do v1, v2
                (markers_pos[idx], v1 + (v2 - 1) + 0.5)
            end
            GLMakie.textlabel!(ax1,
                               markers_ypos,
                               text="$(markers_id[idx]) / $(markers_desc[idx])",
                               text_align=(:left, :center),
                               fontsize=8,
                               cornerradius=0,
                               cornervertices=2,
                               padding=2,
                               strokewidth=1,
                               offset=(0, 5),
                               text_rotation=pi/2)
        end
    end

    if gui

        # selection
        # define a square: Rect(x, y, width, height)
#        square = Rect(2, ax1.limits[][2][1],
#                      6, (ax1.limits[][2][2] - ax1.limits[][2][1]))
#        GLMakie.poly!(ax1,
#                      square,
#                      alpha=0.25,
#                      color=:darkgrey,
#                      strokecolor=:black,
#                      strokewidth=2)

        # time bar
        ax2 = GLMakie.Axis(p[2, 1],
                           xlabel=xl,
                           ylabel="",
                           title="",
                           xticks=LinearTicks(25),
                           yticksvisible=false,
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0),
                           backgroundcolor = :white,
                           xzoomlock=true,
                           yzoomlock=true,
                           xpanlock=true,
                           ypanlock=true,
                           xrectzoom=false,
                           yrectzoom=false)
        GLMakie.xlims!(ax2, t[][1], t[][end])
        GLMakie.ylims!(ax2, 0, 1)
        hideydecorations!(ax2)
        hidexdecorations!(ax2, label=false, ticks=false, ticklabels=false)
        ax2.xticklabelsize = 12

        if markers
            GLMakie.vlines!(ax2,
                            markers_pos,
                            linestyle=:dash,
                            linewidth=1,
                            color=:black)
        end

        # time line marker
        # define a square: Rect(x, y, width, height)
        t_rectangle = lift(seg_pos) do v
            Rect(v, 0, seg_len, 1)
        end
        poly!(ax2,
              t_rectangle,
              color=:darkgrey,
              strokecolor=:black,
              strokewidth=2,
              alpha=0.5)

        # channel bar
        if type === :normal
            ax3 = GLMakie.Axis(p[1, 2],
                               xlabel="",
                               ylabel="",
                               title="",
                               yticks=1:ch_n,
                               xticksvisible=false,
                               yticksvisible=false,
                               yreversed=true,
                               xautolimitmargin=(0, 0),
                               yautolimitmargin=(0, 0),
                               backgroundcolor = :white,
                               xzoomlock=true,
                               yzoomlock=true,
                               xpanlock=true,
                               ypanlock=true,
                               xrectzoom=false,
                               yrectzoom=false)
            ch_n > 1 && (GLMakie.ylims!(ax3, ch_n, 1))
            hidedecorations!(ax3)

            # mark channel types
            if group_ch
                for idx in eachindex(ctypes_pos)
                    GLMakie.hlines!(ax3,
                                    ctypes_pos[idx],
                                    linewidth=5,
                                    color=:black)
                end
            end

            # channel marker
            # define a square: Rect(x, y, width, height)
            ch_rectangle = @lift(Rect(0, $ch1, 1, $nch - 1))
            GLMakie.poly!(ax3,
                          ch_rectangle,
                          color=:darkgrey,
                          strokecolor=:black,
                          strokewidth=2,
                          alpha=0.25)

            # mark bad channels
            for idx in 1:ch_n
                if bad_ch[][idx]
                    GLMakie.hlines!(ax3,
                                    idx,
                                    linestyle=:dots,
                                    linewidth=5,
                                    color=:black)
                end
            end

        end

        on(events(p).mousebutton) do event
            if event.button == Mouse.left
                if event.action == Mouse.press # || event.action == Mouse.release

                    # change time
                    ax2_x = mouseposition(ax2)[1]
                    ax2_y = mouseposition(ax2)[2]
                    seg = (round(Int64, ax2_x), round(Int64, ax2_x) + seg_len)
                    if ax2_x >= 0 && ax2_x <= ax2.limits[][1][2] && ax2_y >= 0 && ax2_y <= 1
                        ax1.limits[] = (seg, ax1.limits[][2])
                        seg_pos[] = round(Int64, ax2_x)
                    end

                    # change channels
                    if type === :normal
                        ax3_x = mouseposition(ax3)[1]
                        ax3_y = mouseposition(ax3)[2]
                        if ax3_x >= 0 && ax3_x <= 1 && ax3_y >= 0 && ax3_y <= ax3.limits[][2][2]
                            ch1[] = floor(Int64, ax3_y)
                            ch1[] > ch_n - nch[] + 1 && (ch1[] = ch_n - nch[] + 1)
                            ax1.limits[] = (ax1.limits[][1], (ch1[] - 0.5, ch1[] + nch[] - 0.5))
                        end
                    end

                end
            end
        end

        on(events(p).keyboardbutton) do event
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
                    if ispressed(p, Keyboard.page_down)
                        if ch_n > 1 && nch[] > 1
                            nch[] -= 1
                            update_ax3 = true
                        end
                    end
                    if ispressed(p, Keyboard.page_up)
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
                    seg_pos[] = ceil(Int64, t[][end] - seg_len)
                    update_ax2 = true
                end
                if event.key == Keyboard.left
                    if seg_pos[] > 0
                        seg_pos[] -= 1
                        update_ax2 = true
                    end
                end
                if ispressed(p, Keyboard.left_shift & Keyboard.left)
                    if seg_pos[] >= 9
                        seg_pos[] -= 9
                        update_ax2 = true
                    end
                end
                if event.key == Keyboard.right
                    if seg_pos[] <= t[][end] - seg_len
                        seg_pos[] += 1
                        update_ax2 = true
                    end
                end
                if ispressed(p, Keyboard.left_shift & Keyboard.right)
                    if seg_pos[] <= t[][end] - seg_len - 9
                        seg_pos[] += 9
                        update_ax2 = true
                    end
                end
                if update_ax2
                    seg = (seg_pos[], seg_pos[] + seg_len)
                    ax1.limits[] = (seg, ax1.limits[][2])
                end
                if update_ax3
                    ax1.limits[] = (ax1.limits[][1], (ch1[] - 0.5, ch1[] + nch[] - 0.5))
                end
            end
        end

        type === :normal && colsize!(p.layout, 2, GLMakie.Fixed(20))
        rowsize!(p.layout, 2, GLMakie.Fixed(20))

    end

    return p

end

"""
    plot_signal(t, s; <keyword arguments>)

Plot epoched signal.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (time points or samples)
- `s::AbstractArray`: data to plot
- `seg::Tuple{Real, Real}`: segment (from, to) in seconds to display
- `clabels::Vector{String}=string.(1:size(s, 1))`: channel labels
- `ctypes:::Vector{String}=repeat([""], size(s, 1))`: channel types
- `cunits::Vector{String}=repeat([""], size(s, 1))`: channel units
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `scale::Bool=true`: draw scale
- `bad::Vector{Bool}=zeros(Bool, size(s, 1))`: list of bad channels

# Returns

- `p::GLMakie.Figure`
"""
function plot_ep(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; seg::Tuple{Real, Real}, clabels::Vector{String}=string.(1:size(s, 1)), ctypes::Vector{String}=repeat([""], size(s, 1)), cunits::Vector{String}=repeat([""], size(s, 1)), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, bad::Vector{Bool}=zeros(Bool, size(s, 1)))::GLMakie.Figure

    ch_n = size(s, 1)

    if ch_n == 1

        fs = round(Int64, 1/(t[2] - t[1]))
        ep_len = size(s, 2) / fs
        ep_n = size(s, 3)
        seg = (0, 5 * ep_len)
        seg_pos = Observable(seg[1])
        seg_len = 5 * ep_len

        # prepare plot
        plot_size = (900, 450)
        p = GLMakie.Figure(size=plot_size)
        ax1 = GLMakie.Axis(p[1, 1],
                           xlabel="",
                           ylabel=ylabel,
                           title=title,
                           xticks=LinearTicks(10),
                           xminorticksvisible=true,
                           xminorticks=IntervalsBetween(10),
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0),
                           xzoomlock=true,
                           yzoomlock=true,
                           xpanlock=true,
                           ypanlock=true,
                           xrectzoom=false,
                           yrectzoom=false)
        GLMakie.xlims!(ax1, seg)
        if minimum(s) == 0
            GLMakie.ylims!(ax1, 0, maximum(s)[2])
        else
            GLMakie.ylims!(ax1, extrema(s))
        end
        ax1.titlesize = 20
        ax1.xlabelsize = 12
        ax1.ylabelsize = 12
        ax1.xticklabelsize = 12
        ax1.yticklabelsize = 12

        # plot signal
        GLMakie.lines!(ax1,
                       t,
                       s[1, :, :][:],
                       linewidth=1.5,
                       color=bad ? :lightgray : :black)

    else

        ctypes_uni = unique(ctypes)
        t_pos = zeros(Int64, length(ctypes_uni))
        for idx in eachindex(ctypes_uni)
             t_pos[idx] = findfirst(isequal(ctypes_uni[idx]), ctypes)
        end
        ctypes_uni_pos = zeros(Int64, length(ctypes))
        ctypes_uni_pos[t_pos] .= 1

        # get ranges of the original signal for the scales
        # normalize and shift so all channels are visible
        r = Float64[]
        @inbounds for ch_idx in eachindex(ctypes_uni)
            push!(r, round(_get_range(s[ctypes .== ctypes_uni[ch_idx], :])))
            s[ctypes .== ctypes_uni[ch_idx], :] = @views normalize_minmax(s[ctypes .== ctypes_uni[ch_idx], :])
        end
        # reverse so 1st channel is on top
        s = @views reverse(s[:, eachindex(t)], dims = 1)
        bad = reverse(bad)
        # each channel is between -1.0 and +1.0
        # scale by 0.5 so maxima do not overlap
        @inbounds for idx in 1:ch_n
            s[idx, :] = @views (s[idx, :]) .+ (idx - 1)
        end

        # prepare plot
        plot_size = (1200, 100 + 40 * ch_n)
        p = GLMakie.Figure(size=plot_size)
        ax = GLMakie.Axis(p[1, 1],
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          xticks=LinearTicks(10),
                          xminorticksvisible=true,
                          xminorticks=IntervalsBetween(10),
                          yticks=((ch_n - 1):-1:0, clabels),
                          yticksvisible=false,
                          xautolimitmargin=(0, 0),
                          yautolimitmargin=(0, 0))
        GLMakie.xlims!(ax, seg)
        GLMakie.ylims!(ax, -1, ch_n)
        ax.titlesize = 20
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = ch_n <= 64 ? 12 : 10;

        # plot channels
        if length(ctypes_uni) > 1
            cmap = GLMakie.resample_cmap(pal, length(ctypes_uni))
        else
            cmap = reverse(GLMakie.resample_cmap(pal, ch_n))
        end

        for idx in 1:ch_n
            if !bad[idx]
                if mono
                    GLMakie.lines!(ax,
                                   t,
                                   s[idx, :],
                                   linewidth=0.75,
                                   color=:black)
                else
                    if length(ctypes_uni) > 1
                        GLMakie.lines!(ax,
                                       t,
                                       s[idx, :],
                                       linewidth=0.75,
                                       color=cmap[channel_color[idx]],
                                       colormap=pal,
                                       colorrange=1:length(ctypes_uni))
                    else
                        GLMakie.lines!(ax,
                                       t,
                                       s[idx, :],
                                       linewidth=0.75,
                                       color=cmap[idx],
                                       colormap=pal,
                                       colorrange=1:ch_n)
                    end
                end
            else
                GLMakie.lines!(ax,
                               t,
                               s[idx, :],
                               linewidth=0.75,
                               alpha=0.2,
                               color=:black)
            end
        end

        # draw scales
        if scale
            idx2 = 1
            for idx1 in 1:ch_n
                if ctypes_uni_pos[idx1] == 1
                    s_pos = ch_n - idx1 + 1
                    GLMakie.lines!(ax,
                                   #[_xlims(t)[1], _xlims(t)[1]],
                                   [ax.limits[][1][1], ax.limits[][1][1]],
                                   [(s_pos - 1.5), (s_pos - 0.5)],
                                   color=:red,
                                   linewidth=5)
                    GLMakie.text!(ax,
                                  (_xlims(t)[1], s_pos - 1),
                                  text="$(r[idx2]) $(cunits[idx1])  ",
                                  fontsize=8,
                                  align=(:center, :top),
                                  rotation=pi/2,
                                  offset=(5, 0))
                    idx2 += 1
                end
            end
        end
    end

    if gui

        ax2 = GLMakie.Axis(p[2, 1],
                           xlabel="Epochs",
                           ylabel="",
                           title="",
                           xticks=LinearTicks(25),
                           yticksvisible=false,
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0),
                           backgroundcolor = :white,
                           xzoomlock=true,
                           yzoomlock=true,
                           xpanlock=true,
                           ypanlock=true,
                           xrectzoom=false,
                           yrectzoom=false)
        GLMakie.xlims!(ax2, 0, ep_n)
        GLMakie.ylims!(ax2, 0, 1)
        hideydecorations!(ax2)
        ax2.xticklabelsize = 12

        # epoch lengths
        for idx in 1:ep_n
            GLMakie.vlines!(ax1,
                            idx * ep_len,
                            linewidth=0.5,
                            linestyle=:dash,
                            color=:black)
        end

        # epoch marker
        # define a square: Rect(x, y, width, height)
        rectangle = lift(seg_pos) do seg_pos
            Rect(seg_pos, 0, 5, 1)
        end
        poly!(ax2,
              rectangle,
              color=:darkgrey,
              strokecolor=:black,
              strokewidth=2)

        ep_selected = zeros(Bool, ep_n)

        on(events(p).mousebutton) do event
            if event.button == Mouse.left
                if event.action == Mouse.press
                    ax1_x = mouseposition(ax1)[1]
                    ax1_y = mouseposition(ax1)[2]
                    nep = ceil(Int64, ax1_x / ep_len)
                    if ax1_y >= ax1.limits[][2][1] && ax1_y <= ax1.limits[][2][2]
                        ep_selected[nep] = !ep_selected[nep]
                    end

                    ax2_x = mouseposition(ax2)[1]
                    ax2_y = mouseposition(ax2)[2]
                    nep = round(Int64, ax2_x)
                    nep < 1 && (nep = 1)
                    seg = ((nep - 1) * ep_len, (nep + 4) * ep_len)
                    if ax2_x >= 0 && ax2_x <= ax2.limits[][1][2] && ax2_y >= 0 && ax2_y <= 1
                        ax1.limits[] = (seg, ax1.limits[][2])
                        seg_pos[] = nep
                    end
                end
            end
        end

        on(events(p).keyboardbutton) do event
            if event.action == Keyboard.press || event.action == Keyboard.repeat
                if event.key == Keyboard.home
                    seg_pos[] = 0
                end
                if event.key == Keyboard._end
                    seg_pos[] = ep_n - 5
                end
                if event.key == Keyboard.left
                    if seg_pos[] >= 1
                        seg_pos[] -= 1
                    end
                end
                if event.key == Keyboard.right
                    if seg_pos[] <= ep_n - 5
                        seg_pos[] += 1
                    end
                end
                seg = (seg_pos[] * ep_len, (seg_pos[] + 5) * ep_len)
                ax1.limits[] = (seg, ax1.limits[][2])
            end
        end

        rowsize!(p.layout, 2, GLMakie.Fixed(20))

    end

    return p

end