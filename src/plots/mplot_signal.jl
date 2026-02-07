export mplot_signal
export mplot_signal_avg
export mplot_signal_butterfly
export mplot

"""
    mplot_signal(t, s; <keyword arguments>)

Plot amplitude of single-channel continuous signal.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (time points or samples)
- `s::AbstractVector`: data to plot
- `seg::Tuple{Real, Real}=(t[1], t[end])`: segment (from, to) in seconds to display
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `bad::Bool=false`: is this a bad channel
- `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

- `p::GLMakie.Figure`
"""
function mplot_signal(t::Union{AbstractVector, AbstractRange}, s::AbstractVector; seg::Tuple{Real, Real}=(t[1], t[end]), xlabel::String="", ylabel::String="", title::String="", bad::Bool=false, gui::Bool=false)::Union{GLMakie.Figure, Nothing}

    seg_len = (seg[2] - seg[1])
    seg_pos = Observable(seg[1])

    # prepare plot
    plot_size = (1600, 450)
    p = GLMakie.Figure(size=plot_size)
    ax1 = GLMakie.Axis(p[1, 1],
                       xlabel=gui ? "" : xlabel,
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
    ax1.xlabelsize = 18
    ax1.ylabelsize = 18
    ax1.xticklabelsize = 12
    ax1.yticklabelsize = 12

    # plot signal
    GLMakie.lines!(ax1,
                   t,
                   s,
                   linewidth=1.5,
                   color=bad ? :lightgray : :black)

    if gui

        # time line
        ax2 = GLMakie.Axis(p[2, 1],
                           xlabel=xlabel,
                           ylabel="",
                           title="",
                           xticks=LinearTicks(25),
                           yticksvisible=false,
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0),
                           backgroundcolor = :white,
                           xzoomlock=false,
                           yzoomlock=true,
                           xpanlock=true,
                           ypanlock=true,
                           xrectzoom=false,
                           yrectzoom=false)
        GLMakie.xlims!(ax2, t[1], t[end])
        GLMakie.ylims!(ax2, 0, 1)
        hideydecorations!(ax2)
        ax2.xticklabelsize = 12

        # time line marker
        # define a square: Rect(x, y, width, height)
        rectangle = lift(seg_pos) do seg_pos
            Rect(seg_pos, 0, seg_len, 1)
        end
        poly!(ax2,
              rectangle,
              color=:darkgrey,
              strokecolor=:black,
              strokewidth=2)

        on(events(p).mousebutton) do event
            if event.button == Mouse.left
                if event.action == Mouse.press
                    ax2_x = mouseposition(ax2)[1]
                    ax2_y = mouseposition(ax2)[2]
                    seg = (round(Int64, ax2_x), round(Int64, ax2_x) + seg_len)
                    if ax2_x >= 0 && ax2_x <= ax2.limits[][1][2] && ax2_y >= 0 && ax2_y <= 1
                        ax1.limits[] = (seg, ax1.limits[][2])
                        seg_pos[] = round(Int64, ax2_x)
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
                    seg_pos[] = ceil(Int64, t[end] - seg_len)
                end
                if event.key == Keyboard.left
                    if seg_pos[] > 0
                        seg_pos[] -= 1
                    end
                end
                if ispressed(p, Keyboard.left_shift & Keyboard.left)
                    if seg_pos[] >= 9
                        seg_pos[] -= 9
                    end
                end
                if event.key == Keyboard.right
                    if seg_pos[] <= t[end] - seg_len
                        seg_pos[] += 1
                    end
                end
                if ispressed(p, Keyboard.left_shift & Keyboard.right)
                    if seg_pos[] <= t[end] - seg_len - 9
                        seg_pos[] += 9
                    end
                end
                seg = (seg_pos[], seg_pos[] + seg_len)
                ax1.limits[] = (seg, ax1.limits[][2])
            end
        end

        rowsize!(p.layout, 2, GLMakie.Fixed(20))

#        wait(display(p))

#        return nothing

    end

    return p

end

"""
    mplot_signal(t, s; <keyword arguments>)

Plot amplitude of single-channel epoched signal.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (time points or samples)
- `s::AbstractArray`: data to plot
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `bad::Bool=false`: is this a bad channel
- `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

- `p::GLMakie.Figure`
"""
function mplot_signal(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; ylabel::String="", title::String="", bad::Bool=false, gui::Bool=true)::GLMakie.Figure

    @assert size(s, 1) == 1 "Signal must contain only one channel."

    fs = round(Int64, 1/(t[2] - t[1]))
    ep_len = size(s, 2) / fs
    ep_n = size(s, 3)
    seg = (0, 5 * ep_len)
    seg_pos = Observable(seg[1])
    seg_len = 5 * ep_len

    # prepare plot
    plot_size = (1600, 450)
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
    ax1.xlabelsize = 18
    ax1.ylabelsize = 18
    ax1.xticklabelsize = 12
    ax1.yticklabelsize = 12

    # plot signal
    GLMakie.lines!(ax1,
                   t,
                   s[1, :, :][:],
                   linewidth=1.5,
                   color=bad ? :lightgray : :black)

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
            Rect(seg_pos, 0, seg_len, 1)
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
                        @show ep_selected
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
                if event.key == Keyboard.left
                    if seg_pos[] >= 1
                        @show seg_pos[]
                        seg_pos[] -= 1
                        seg = (seg_pos[] * ep_len, (seg_pos[] + 5) * ep_len)
                        ax1.limits[] = (seg, ax1.limits[][2])
                    end
                end
                if event.key == Keyboard.right
                    if seg_pos[] <= ep_n - 5
                        @show seg_pos[]
                        seg_pos[] += 1
                        seg = (seg_pos[] * ep_len, (seg_pos[] + 5) * ep_len)
                        ax1.limits[] = (seg, ax1.limits[][2])
                    end
                end
            end
        end

        rowsize!(p.layout, 2, GLMakie.Fixed(20))

        wait(display(p))

    end

    return p

end

"""
    mplot_signal(t, s; <keyword arguments>)

Plot amplitude of multi-channel continuous signal.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (time points or samples)
- `s::AbstractMatrix`: data to plot
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display
- `clabels::Vector{String}=string.(1:size(s, 1))`: channel labels
- `ctypes:::Vector{String}=repeat([""], size(s, 1))`: channel types
- `cunits::Vector{String}=repeat([""], size(s, 1))`: channel units
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `scale::Bool=true`: draw scale
- `bad::Vector{Bool}=zeros(Bool, size(s, 1))`: list of bad channels
- `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

- `p::GLMakie.Figure`
"""
function mplot_signal(t::Union{AbstractVector, AbstractRange}, s::AbstractMatrix; seg::Tuple{Real, Real}=(0, 10), clabels::Vector{String}=string.(1:size(s, 1)), ctypes::Vector{String}=repeat([""], size(s, 1)), cunits::Vector{String}=repeat([""], size(s, 1)), xlabel::String="", ylabel::String="", title::String="", scale::Bool=true, bad::Vector{Bool}=zeros(Bool, size(s, 1)), gui::Bool=true)::GLMakie.Figure

    if clabels != string.(1:size(s, 1))
        l = length.(clabels)
        ml = maximum(l)
        for idx in eachindex(clabels)
            length(clabels[idx]) < ml && (clabels[idx] = lpad(clabels[idx], ml, ' '))
        end
    end

    # zooming factor
    zoom = Observable(1.0)

    seg_len = (seg[2] - seg[1])
    seg_pos = Observable(seg[1])

    ch_n = size(s, 1)

    if gui
        if ch_n > 20
            ch1 = Observable(1)
            ch2 = ch1[] + 19
        else
            ch1 = Observable(1)
            ch2 = ch_n
        end
    else
        ch1 = Observable(1)
        ch2 = ch_n
    end

    # order by ctypes
    # and markers for ax3

    ctypes_uni = unique(ctypes)

    t_pos = zeros(Int64, length(ctypes_uni))
    for idx in eachindex(ctypes_uni)
         t_pos[idx] = findfirst(isequal(ctypes_uni[idx]), ctypes)
    end
    ctypes_uni_pos = zeros(Int64, length(ctypes))
    ctypes_uni_pos[t_pos] .= 1

    # get ranges of the original signal for the scales
    # normalize in groups by channel type
    # between -1.0 and +1.0 and shift so all channels are visible
    r = Float64[]
    @inbounds for idx in eachindex(ctypes_uni)
        push!(r, round(NeuroAnalyzer._get_range(s[ctypes .== ctypes_uni[idx], :])))
        s[ctypes .== ctypes_uni[idx], :] = @views normalize_minmax(s[ctypes .== ctypes_uni[idx], :])
#        s = normalize_minmax(s, bych=true)
    end

    # reverse so 1st channel is on top
    # reverse!(clabels)
    # reverse!(ctypes)
    # reverse!(s, dims=1)
    # bad = reverse(bad)

    # each channel is 
    # scale by 0.5 so maxima do not overlap
    #    @inbounds for idx in 1:ch_n
    #        s[idx, :] = @views (s[idx, :] .* 0.5) .+ idx
    #    end

    # y-axis labels colors
    ytc = repeat([:black], ch_n)
    ytc[bad] .= :lightgray

    # prepare plot
    if gui
        plot_size = (1650, 950)
    else
        plot_size = (1650, ch_n * 50)
    end
    p = GLMakie.Figure(size=plot_size)
    ax1 = GLMakie.Axis(p[1, 1],
                       xlabel="",
                       ylabel=ylabel,
                       title=title,
                       xticks=LinearTicks(10),
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(10),
                       yticks=clabels==string.(1:ch_n) ? (1:ch_n, string.(1:ch_n)) : (1:ch_n, clabels),
                       # TO DO: yticklabelcolor=ytc[1:end],
                       yreversed=true,
                       xautolimitmargin=(0, 0),
                       yautolimitmargin=(0, 0),
                       xzoomlock=true,
                       yzoomlock=true,
                       xpanlock=true,
                       ypanlock=true,
                       xrectzoom=false,
                       yrectzoom=false)
    GLMakie.xlims!(ax1, seg)
    if gui
        if ch_n > 20
            GLMakie.ylims!(ax1, ch2 + 0.5, ch2 - 20 + 0.5)
        else
            GLMakie.ylims!(ax1, ch_n + 0.5, 0.5)
        end
    else
        GLMakie.ylims!(ax1, ch_n + 0.5, 0.5)
    end
    ax1.titlesize = 20
    ax1.xlabelsize = 18
    ax1.ylabelsize = 18
    ax1.xticklabelsize = 12
    ax1.yticklabelsize = 12

    # draw channels
    for idx in ch_n:-1:1
        GLMakie.lines!(ax1,
                       t,
                       @lift((s[idx, :] .* $zoom .* 0.5) .+ idx),
                       linewidth=1.5,
                       color=bad[idx] ? :lightgray : :black)
    end

    # draw scale bars
    if scale
        idx2 = 1
        for idx1 in 1:ch_n
            if ctypes_uni_pos[idx1] == 1
                x = ax1.limits[][1][1]
                GLMakie.lines!(ax1,
                               [x, x],
                               [(idx1 - 0.5), (idx1 + 0.5)],
                               color=:red,
                               linewidth=5)
                GLMakie.text!(ax1,
                              x,
                              idx1,
                              text=cunits == repeat([""], ch_n) ? "$(r[idx2])" : "$(r[idx2]) $(cunits[idx1])",
                              fontsize=8,
                              color=:red,
                              align=(:center, :top),
                              rotation=pi/2,
                              offset=(5, 0))
                idx2 += 1
            end
        end
    end

    if gui

        # selection
        # define a square: Rect(x, y, width, height)
        square = Rect(2, ax1.limits[][2][1],
                      6, (ax1.limits[][2][2] - ax1.limits[][2][1]))
        poly!(ax1,
              square,
              alpha=0.25,
              color=:darkgrey,
              strokecolor=:black,
              strokewidth=2)

        # time bar
        ax2 = GLMakie.Axis(p[2, 1],
                           xlabel=xlabel,
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
        GLMakie.xlims!(ax2, t[1], t[end])
        GLMakie.ylims!(ax2, 0, 1)
        hideydecorations!(ax2)
        ax2.xticklabelsize = 12

        # time line marker
        # define a square: Rect(x, y, width, height)
        t_rectangle = lift(seg_pos) do seg_pos
            Rect(seg_pos, 0, seg_len, 1)
        end
        poly!(ax2,
              t_rectangle,
              color=:darkgrey,
              strokecolor=:black,
              strokewidth=2)

        # channel bar
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
        GLMakie.ylims!(ax3, ch_n, 1)
        hidedecorations!(ax3)

        # channel marker
        # define a square: Rect(x, y, width, height)
        ch_rectangle = lift(ch1) do ch1
            Rect(0, ch1, 1, 19)
        end
        poly!(ax3,
              ch_rectangle,
              color=:darkgrey,
              strokecolor=:black,
              strokewidth=1)

        on(events(p).mousebutton) do event
            if event.button == Mouse.left
                if event.action == Mouse.press # || event.action == Mouse.release
                    ax2_x = mouseposition(ax2)[1]
                    ax2_y = mouseposition(ax2)[2]
                    seg = (round(Int64, ax2_x), round(Int64, ax2_x) + seg_len)
                    if ax2_x >= 0 && ax2_x <= ax2.limits[][1][2] && ax2_y >= 0 && ax2_y <= 1
                        ax1.limits[] = (seg, ax1.limits[][2])
                        seg_pos[] = round(Int64, ax2_x)
                    end
                    ax3_x = mouseposition(ax3)[1]
                    ax3_y = mouseposition(ax3)[2]
                    if ax3_x >= 0 && ax3_x <= 1 && ax3_y >= 0 && ax3_y <= ax3.limits[][2][2]
                        ch1[] = floor(Int64, ax3_y)
                        ch1[] > ch_n - 19 && (ch1[] = ch_n - 19)
                        ax1.limits[] = (ax1.limits[][1], (ch1[] - 0.5, ch1[] + 20 - 0.5))
                    end
                end
            end
        end

        on(events(p).keyboardbutton) do event
            if event.action == Keyboard.press || event.action == Keyboard.repeat
                if event.key == Keyboard.kp_add
                    if zoom[] < 20
                        zoom[] += 0.25
                    end
                end
                if event.key == Keyboard.kp_subtract
                    if zoom[] > 1
                        zoom[] -= 0.25
                    end
                end
                if event.key == Keyboard.home
                    seg_pos[] = 0
                end
                if event.key == Keyboard._end
                    seg_pos[] = ceil(Int64, t[end] - seg_len)
                end
                if event.key == Keyboard.left
                    if seg_pos[] > 0
                        seg_pos[] -= 1
                    end
                end
                if ispressed(p, Keyboard.left_shift & Keyboard.left)
                    if seg_pos[] >= 9
                        seg_pos[] -= 9
                    end
                end
                if event.key == Keyboard.right
                    if seg_pos[] <= t[end] - seg_len
                        seg_pos[] += 1
                    end
                end
                if ispressed(p, Keyboard.left_shift & Keyboard.right)
                    if seg_pos[] <= t[end] - seg_len - 9
                        seg_pos[] += 9
                    end
                end
                seg = (seg_pos[], seg_pos[] + seg_len)
                ax1.limits[] = (seg, ax1.limits[][2])
            end
        end

        colsize!(p.layout, 2, GLMakie.Fixed(20))
        rowsize!(p.layout, 2, GLMakie.Fixed(20))

#        wait(display(p))

    end

    return p

end

"""
    mplot_signal(t, s; <keyword arguments>)

Plot amplitude of multi-channel epoched signal.

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
function mplot_signal(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; seg::Tuple{Real, Real}, clabels::Vector{String}=string.(1:size(s, 1)), ctypes::Vector{String}=repeat([""], size(s, 1)), cunits::Vector{String}=repeat([""], size(s, 1)), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, bad::Vector{Bool}=zeros(Bool, size(s, 1)))::GLMakie.Figure

    ch_n = size(s, 1)

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
        s[idx, :] = @views (s[idx, :] .* 0.5) .+ (idx - 1)
    end

    # prepare plot
    plot_size = (1600, 100 + 40 * ch_n)
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
@show ax.limits[][1]
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

    return p

end

"""
    mplot_signal(t, s1, s2; <keyword arguments>)

Plot amplitude of single-channel signal.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s1::AbstractVector`: data to plot
- `s2::AbstractVector`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title

# Returns

- `p::GLMakie.Figure`
"""
function mplot_signal(t::Union{AbstractVector, AbstractRange}, s1::AbstractVector, s2::AbstractVector; xlabel::String="", ylabel::String="", title::String="")::GLMakie.Figure

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # prepare plot
    plot_size = (1600, 400)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(10),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0))
    GLMakie.xlims!(ax, _xlims(t))
    if minimum(s1) == 0 || minimum(s2) == 0
        if maximum(s1) > maximum(s2)
            GLMakie.ylims!(ax, 0, _ylims(s1)[2])
        else
            GLMakie.ylims!(ax, 0, _ylims(s2)[2])
        end
    else
        if maximum(s1) < maximum(s2)
            if minimum(s1) < minimum(s2)
                GLMakie.ylims!(ax, _ylims(s1)[1], _ylims(s2)[2])
            else
                GLMakie.ylims!(ax, _ylims(s2)[1], _ylims(s2)[2])
            end
        else
            if minimum(s1) < minimum(s2)
                GLMakie.ylims!(ax, _ylims(s1)[1], _ylims(s1)[2])
            else
                GLMakie.ylims!(ax, _ylims(s2)[1], _ylims(s1)[2])
            end
        end
    end
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot signals
    GLMakie.lines!(t,
                   s1,
                   linewidth=1,
                   alpha=0.5,
                   color=:black)
    GLMakie.lines!(t,
                   s2,
                   linewidth=1,
                   alpha=0.5,
                   color=:blue)

    return p

end

"""
    mplot_signal_avg(t, signal; <keyword arguments>)

Plot amplitude mean and ±95% CI of averaged `signal` channels.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette

# Returns

- `p::GLMakie.Figure`
"""
function mplot_signal_avg(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false)::GLMakie.Figure

    pal = mono ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = NeuroAnalyzer.msci95(s)

    # get limits
    ylim = (round(minimum(s_l) * 1.5, digits=0), round(maximum(s_u) * 1.5, digits=0))
    ylim = _tuple_max(ylim)

    # prepare plot
    plot_size = (1600, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(10),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0))
    GLMakie.xlims!(ax, _xlims(t))
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot upper 95% CI
    GLMakie.band!(t,
                  s_u,
                  s_l,
                  alpha=0.25,
                  color=:grey,
                  strokewidth=0.5)
    # plot mean
    GLMakie.lines!(t,
                   s_m,
                   color=:black,
                   linewidth=1)

    return p

end

"""
    mplot_signal_butterfly(t, s; <keyword arguments>)

Butterfly plot of `s` channels.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `avg::Bool=false`: plot average channels
- `mono::Bool=false`: use color or gray palette

# Returns

- `p::GLMakie.Figure`
"""
function mplot_signal_butterfly(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", avg::Bool=true, mono::Bool=false)::GLMakie.Figure

    pal = mono ? :grays : :darktest

    ch_n = size(s, 1)

    # get limits
    ylim = (round(minimum(s) * 1.5, digits=0), round(maximum(s) * 1.5, digits=0))
    ylim = _tuple_max(ylim)

    # channel labels
    clabels == [""] && (clabels = repeat([""], ch_n))

    # plot channels
    plot_size = (1600, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(10),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0))
    GLMakie.xlims!(ax, _xlims(t))
    cmap = GLMakie.resample_cmap(pal, ch_n)
    for idx in 1:ch_n
        GLMakie.lines!(t,
                       s[idx, :],
                       color=cmap[idx],
                       colormap=pal,
                       colorrange=1:ch_n,
                       linewidth=0.5,
                       label=clabels[idx])
    end
    ch_n < 40 && axislegend(position = :rb)

    # plot averaged channels
    if avg
        s = mean(s, dims=1)[:]
        GLMakie.lines!(t,
                       s,
                       linewidth=2,
                       color=:black)
    end

    return p

end

"""
    mplot_signal(t, s1, s2; <keyword arguments>)

Plot amplitude of multi-channel signals.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s1::AbstractArray`: data to plot
- `s2::AbstractArray`: data to plot
- `clabels::Vector{String}=string.(1:size(s, 1))`: channel labels
- `ctypes:::Vector{String}=repeat([""], size(s1, 1))`: channel types
- `cunits::Vector{String}=repeat([""], size(s1, 1))`: channel units
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `scale::Bool=true`: draw scale

# Returns

- `p::GLMakie.Figure`
"""
function mplot_signal(t::Union{AbstractVector, AbstractRange}, s1::AbstractArray, s2::AbstractArray; clabels::Vector{String}=string.(1:size(s, 1)), ctypes::Vector{String}=repeat([""], size(s1, 1)), cunits::Vector{String}=repeat([""], size(s1, 1)), xlabel::String="", ylabel::String="", title::String="", scale::Bool=true)::GLMakie.Figure

    ch_n = size(s1, 1)

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
        push!(r, round(_get_range(s1[ctypes .== ctypes_uni[ch_idx], :])))
        s1[ctypes .== ctypes_uni[ch_idx], :] = @views normalize_minmax(s1[ctypes .== ctypes_uni[ch_idx], :])
        s2[ctypes .== ctypes_uni[ch_idx], :] = @views normalize_minmax(s2[ctypes .== ctypes_uni[ch_idx], :])
    end
    # reverse so 1st channel is on top
    s1 = @views reverse(s1[:, eachindex(t)], dims = 1)
    s2 = @views reverse(s2[:, eachindex(t)], dims = 1)
    # each channel is between -1.0 and +1.0
    # scale by 0.5 so maxima do not overlap
    @inbounds for idx in 1:ch_n
        s1[idx, :] = @views (s1[idx, :] .* 0.5) .+ (idx - 1)
        s2[idx, :] = @views (s2[idx, :] .* 0.5) .+ (idx - 1)
    end

    # prepare plot
    plot_size = (1600, 100 + 40 * ch_n)
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
    GLMakie.xlims!(ax, _xlims(t))
    GLMakie.ylims!(ax, -1, ch_n)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot channels
    for idx in 1:ch_n
        GLMakie.lines!(t,
                       s1[idx, :],
                       linewidth=1,
                       color=:black,
                       alpha=0.5)
        GLMakie.lines!(t,
                       s2[idx, :],
                       linewidth=1,
                       color=:blue,
                       alpha=0.5)
    end

    # draw scale
    if scale
        idx2 = 1
        for idx1 in 1:ch_n
            if ctypes_uni_pos[idx1] == 1
                s_pos = ch_n - idx1 + 1
                GLMakie.lines!(ax,
                               [_xlims(t)[1], _xlims(t)[1]],
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

    return p

end

"""
    mplot(obj; <keyword arguments>)

Plot signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ep::Int64=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name or list of channel names
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title
- `mono::Bool=false`: use color or gray palette
- `emarkers::Bool`: draw epoch markers if available
- `markers::Bool`: draw markers if available
- `scale::Bool=true`: draw scale
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:mean`: mean ± 95%CI
    - `:butterfly`: butterfly plot
- `avg::Bool=false`: plot averaged channel in butterfly plot
- `bad::Bool=false`: plot bad channels
- `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

- `p::GLMakie.Figure`
"""
function mplot(obj::NeuroAnalyzer.NEURO; ep::Int64=0, ch::Union{String, Vector{String}, Regex}="all", seg::Tuple{Real, Real}=(0, 10), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, markers::Bool=true, scale::Bool=true, type::Symbol=:normal, avg::Bool=true, bad::Bool=true, gui::Bool=true)::GLMakie.Figure

    datatype(obj) == "erp" && _warn("For ERP objects, use plot_erp()")
    datatype(obj) == "erf" && _warn("For ERF objects, use plot_erp()")
    datatype(obj) == "mep" && _warn("For MEP objects, use plot_mep()")

    if signal_len(obj) <= seg[2] * sr(obj)
        seg = (obj.time_pts[1], obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end

    _check_var(type, [:normal, :butterfly, :mean], "type")

    ep1 = nothing
    ep2 = nothing
    if ep != 0
        _check_epochs(obj, ep)
        @assert nepochs(obj) > 1 "To use ep the signal must be epoched."
        ep1 = (ep - 1) * epoch_len(obj) + 1
        ep2 = ep * epoch_len(obj)
        seg = e2t(obj, ep)
    end

    # do not show epoch markers if there are no epochs
    nepochs(obj) == 1 && (emarkers = false)
    if emarkers
        epoch_markers = _get_epoch_markers(obj)
    end

    # check channels
    ch = get_channel(obj, ch=ch)
    ch_n = length(ch)

    length(ch) == 1 && (ch = ch[1])

    xl, yl, tt = "", "", ""

    bm = obj.header.recording[:bad_channel]
    ctypes = obj.header.recording[:channel_type]
    clabels = labels(obj)

    # sort channels by their type
    if !isa(ch, Int64)
        ctypes = ctypes[ch]
        clabels = clabels[ch]
        cunits = obj.header.recording[:unit][ch]
    end
    if type === :normal
        if isa(ch, Int64)
            ch_name = _ch_rename(ctypes[ch])
            xl, yl, tt = _set_defaults(xlabel,
                                       ylabel,
                                       title,
                                       "Time [s]",
                                       "",
                                       "Channel $ch: $(clabels[ch]) ($ch_name)")
            if datatype(obj) == "eda"
                ylabel == "default" && (yl = "Impedance [μS]")
            else
                ylabel == "default" && (yl = "Amplitude [$(_ch_units(obj, clabels[ch]))]")
            end
            if ep == 0

                if nepochs(obj) == 1
                    p = mplot_signal(obj.time_pts,
                                    obj.data[ch, :, :][:],
                                    seg=seg,
                                    xlabel=xl,
                                    ylabel=yl,
                                    title=tt,
                                    bad=bm[ch],
                                    gui=gui)
                else
                    p = mplot_signal(obj.time_pts,
                                    reshape(obj.data[ch, :, :], 1, :, nepochs(obj)),
                                    ylabel=yl,
                                    title=tt,
                                    bad=bm[ch],
                                    gui=gui)
                end
            else
                p = mplot_signal(obj.time_pts[ep1:ep2],
                                obj.data[ch, :, ep],
                                seg=round.(Int64, seg),
                                xlabel=xl,
                                ylabel=yl,
                                title=tt,
                                bad=bm[ch],
                                gui=false)
            end
        else
            xl, yl, tt = _set_defaults(xlabel,
                                       ylabel,
                                       title,
                                       "Time [s]",
                                       "",
                                       "")
            p = mplot_signal(obj.time_pts,
                            obj.data[ch, :, :][:, :],
                            seg=seg,
                            ctypes=ctypes,
                            clabels=clabels,
                            cunits=cunits,
                            xlabel=xl,
                            ylabel=yl,
                            title=tt,
                            bad=bm[ch],
                            scale=scale,
                            gui=gui)
        end
    end

    if type === :butterfly
        @assert length(unique(ctypes)) == 1 "For plot type=:butterfly all channels must be of the same type."
        @assert ch_n > 1 "For plot type=:butterfly the signal must contain ≥ 2 channels."
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [s]",
                                   datatype(obj) == "eda" ? "Impedance [μS]" : "Amplitude [$(obj.header.recording[:unit][ch[1]])]",
                                   "$ch_n $(uppercase(unique(ctypes)[1])) channels")
        p = mplot_signal_butterfly(obj.time_pts,
                                  obj.data[ch, :, :][:, :],
                                  clabels=clabels,
                                  xlabel=xl,
                                  ylabel=yl,
                                  title=tt,
                                  avg=avg,
                                  mono=mono)
    end

    if type === :mean
        @assert length(unique(ctypes)) == 1 "For plot type=:mean all channels must be of the same type."
        @assert ch_n > 1 "For plot type=:mean the signal must contain ≥ 2 channels."
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [s]",
                                   datatype(obj) == "eda" ? "Impedance [μS]" : "Amplitude [$(obj.header.recording[:unit][ch[1]])]",
                                   "$ch_n $(uppercase(unique(ctypes)[1])) channels")
        p = mplot_signal_avg(obj.time_pts,
                            obj.data[ch, :, :][:, :],
                            xlabel=xl,
                            ylabel=yl,
                            title=tt,
                            mono=mono)
    end

    # add epochs markers
    # TODO: draw epoch numbers
    if emarkers
        GLMakie.vlines!(p[1, 1],
                        epoch_markers,
                        linestyle=:dot,
                        linewidth=0.5,
                        color=:blue)
    end

    # plot markers if available
    if markers && _has_markers(obj)
        markers_pos = obj.markers[!, :start]
        markers_id = obj.markers[!, :id]
        markers_desc = obj.markers[!, :value]
        if gui
            GLMakie.vlines!(p[2, 1],
                            markers_pos,
                            linestyle=:dash,
                            linewidth=1,
                            color=:black)
        end
        for idx in eachindex(markers_pos)
            if _in(markers_pos[idx], (obj.time_pts[1], obj.time_pts[end]))
                GLMakie.vlines!(p[1, 1],
                                markers_pos[idx],
                                linestyle=:dash,
                                linewidth=1,
                                color=:black)
                if length(ch) > 1
                    GLMakie.textlabel!(p[1, 1],
                                       (markers_pos[idx] + 0.07, ch_n > 20 ? 20.4 : ch_n + 0.4),
                                       text="$(markers_id[idx]) / $(markers_desc[idx])",
                                       text_align=(:left, :center),
                                       fontsize=8,
                                       text_rotation=pi/2)
                else
                    GLMakie.textlabel!(p[1, 1],
                                       (markers_pos[idx] + 0.07, 0.97 * minimum(obj.data[ch, :, :])),
                                       text="$(markers_id[idx]) / $(markers_desc[idx])",
                                       text_align=(:left, :center),
                                       fontsize=8,
                                       text_rotation=pi/2)
                end
            end
        end
    end

    if gui
#        time_s = Slider(p[2, 1],
#                        range = LinRange(t[1], t[end], length(t)),
#                        horizontal = true)
#
#        wait(display(p))
    end

    return p

end

"""
    plot(obj1, obj2; <keyword arguments>)

Plot signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name or list of channel names
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title
- `scale::Bool=true`: draw scale

# Returns

- `p::GLMakie.Figure`
"""
function mplot(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ep::Union{Int64, AbstractRange}=0, ch::Union{String, Vector{String}, Regex}="all", seg::Tuple{Real, Real}=(0, 10), xlabel::String="default", ylabel::String="default", title::String="default", scale::Bool=true)::GLMakie.Figure

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert size(obj1.data) == size(obj2.data) "Signals of OBJ1 and OBJ2 must have the same size."
    @assert datatype(obj1) == obj2.header.recording[:data_type] "OBJ1 and OBJ2 must have the same data type."

    if signal_len(obj1) <= 10 * sr(obj1) && seg == (0, 10)
        seg = (obj1.time_pts[1], obj1.time_pts[end])
    else
        _check_segment(obj1, seg)
    end
    seg = (vsearch(seg[1], obj1.time_pts), vsearch(seg[2], obj1.time_pts))

    if ep != 0
        _check_epochs(obj1, ep)
        if nepochs(obj1) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj1) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj1) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj1)))
            end
            ep = 0
        end
    end

    # check channels
    @assert labels(obj1)[get_channel(obj1, ch=ch)] == labels(obj2)[get_channel(obj1, ch=ch)] "OBJ1 and OBJ2 channel labels must be the same."
    _ = get_channel(obj2, ch=ch)
    ch = get_channel(obj1, ch=ch)
    clabels = labels(obj1)

    # get time vector
    if seg[2] <= epoch_len(obj1)
        s1 = obj1.data[:, seg[1]:seg[2], 1]
        s2 = obj2.data[:, seg[1]:seg[2], 1]
    else
        s1 = epoch(obj1, ep_n=1).data[:, seg[1]:seg[2], 1]
        s2 = epoch(obj2, ep_n=1).data[:, seg[1]:seg[2], 1]
    end
    t = obj1.time_pts[seg[1]:seg[2]]

    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj1, seg[1], seg[2])

    (ch isa(Vector{Int64}) && length(ch) == 1) && (ch = ch[1])

    xl, yl, tt = "", "", ""

    # sort channels by their type
    ctypes = obj1.header.recording[:channel_type]
    if !isa(ch, Int64)
        s1 = @views s1[ch, :]
        s2 = @views s2[ch, :]
        ctypes = ctypes[ch]
        clabels = clabels[ch]
        cunits = obj1.header.recording[:unit][ch]
    end

    if isa(ch, Int64)
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [s]",
                                   "",
                                   "")
        ylabel == "default" && (yl = "Amplitude [$(_ch_units(obj1, labels(obj1)[ch]))]")
        p = mplot_signal(t,
                        vec(s1[ch, :]),
                        vec(s2[ch, :]),
                        xlabel=xl,
                        ylabel=yl,
                        title=tt)
    else
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [s]",
                                   "",
                                   "")
        p = mplot_signal(t,
                        s1,
                        s2,
                        ctypes=ctypes,
                        clabels=clabels,
                        cunits=cunits,
                        xlabel=xl,
                        ylabel=yl,
                        title=tt,
                        scale=scale)
    end

    return p

end
