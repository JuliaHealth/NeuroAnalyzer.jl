export plot_ica

# scalp map of every independent component: the component weights at the channels (columns
# of the mixing matrix), interpolated on a grid; returns the maps (NaN outside the head),
# grid coordinates, electrode positions and axis limits
function _ica_maps(
    obj::NeuroAnalyzer.NEURO,
    chn::Vector{Int64},
    ic_mw::Matrix{Float64},
    imethod::Symbol;
    iter::Int64 = 64,
)
    # locations are stored in their original units: scale the montage to the unit sphere
    locs = _locs_unit(obj.locs)
    rows = Int64[]
    for l in labels(obj)[chn]
        r = "label" in names(locs) && DataFrames.nrow(locs) > 0 ? _find_bylabel(locs, l) : Int64[]
        r isa Int64 || throw(ArgumentError("Channel $l has no location; add locations first."))
        push!(rows, r)
    end
    loc_x = zeros(length(rows))
    loc_y = zeros(length(rows))
    for (idx, r) in enumerate(rows)
        loc_x[idx], loc_y[idx] = pol2cart(locs[r, :loc_radius], locs[r, :loc_theta])
    end

    head12 = maximum(abs, loc_x) < 1.2 && maximum(abs, loc_y) < 1.2
    lim = head12 ? 1.2 : 1.4
    # the maps fill the head circle (radius 1), or reach the outermost electrodes if they lie
    # outside it
    rmask = min(lim, max(1.0, maximum(hypot.(loc_x, loc_y))))

    maps = Matrix{Float64}[]
    gx = Float64[]
    gy = Float64[]
    for idx in axes(ic_mw, 2)
        m, gx, gy = _interpolate2d(ic_mw[:, idx], loc_x, loc_y, iter, imethod, :none)
        m = Matrix(m'[:, end:-1:1])
        for i in eachindex(gx), j in eachindex(gy)
            hypot(gx[i], gy[j]) > rmask && (m[i, j] = NaN)
        end
        push!(maps, m)
    end

    return (; maps, gx, gy, loc_x, loc_y, lim)
end

# symmetric color range: the sign of a component is arbitrary, zero is always in the middle
function _ica_crange(m::AbstractMatrix)::Tuple{Float64, Float64}
    v = filter(!isnan, m)
    a = isempty(v) ? 1.0 : maximum(abs, v)
    a == 0 && (a = 1.0)
    return (-a, a)
end

# build the plot_ica() window (without displaying it); the widgets are returned for testing
function _plot_ica_build(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_var::Union{Nothing, Vector{Float64}} = nothing,
    imethod::Symbol = :tp,
    flim::Tuple{Real, Real} = (0, min(sr(obj) / 2, 80)),
    wlen::Real = 10,
)
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("plot_ica() works for continuous (single-epoch) objects only."))
    _check_var(imethod, _IMETHODS, "imethod")
    chn = get_channel(obj; ch = ch)
    isempty(chn) && throw(ArgumentError("No channels selected."))
    size(ic_mw, 1) == length(chn) ||
        throw(
            ArgumentError("ic_mw has $(size(ic_mw, 1)) rows, $(length(chn)) channels are selected."),
        )
    size(ic, 1) == size(ic_mw, 2) ||
        throw(ArgumentError("Dimension mismatch between ic $(size(ic)) and ic_mw $(size(ic_mw))."))
    size(ic, 2) == epoch_len(obj) ||
        throw(ArgumentError("ic has $(size(ic, 2)) samples, OBJ has $(epoch_len(obj))."))
    isnothing(ic_var) ||
        length(ic_var) == size(ic, 1) ||
        throw(ArgumentError("ic_var must contain $(size(ic, 1)) values."))
    _check_tuple(flim, (0, sr(obj) / 2), "flim")
    wlen > 0 || throw(ArgumentError("wlen must be > 0."))

    ic_n = size(ic, 1)
    fs = sr(obj)
    t = obj.time_pts
    t_first = t[1]
    t_last = t[end]
    clabels = labels(obj)[chn]
    units = obj.header.recording[:unit][chn]

    sm = _ica_maps(obj, chn, ic_mw, imethod)
    # large map of the selected component (finer grid)
    lm = _ica_maps(obj, chn, ic_mw[:, 1:1], imethod; iter = 128)
    big_map(i) = i == 1 ? lm.maps[1] : _ica_maps(obj, chn, ic_mw[:, i:i], imethod; iter = 128).maps[1]

    # state
    cur_ic = Observable(1)
    cur_ch = Observable(1)
    marked = Observable(falses(ic_n))
    applied = Ref(false)
    psd_cache = Dict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
    ctl = Dict{Symbol, Any}()

    function ic_psd(i::Int64)
        get!(psd_cache, i) do
            r = _psd(ic[i, :]; fs = fs, db = true)
            keep = flim[1] .<= r.f .<= flim[2]
            (r.f[keep], r.p[keep])
        end
    end

    ic_title(i) =
        isnothing(ic_var) ? "IC $i" : "IC $i ($(round(ic_var[i]; digits = 1))% variance)"

    # ---------------------------------------------------------------- #
    # layout                                                            #
    # ---------------------------------------------------------------- #
    fig = GLMakie.Figure(; size = (1600, 1000))
    left = fig[1, 1] = GLMakie.GridLayout()
    right = fig[1, 2] = GLMakie.GridLayout()
    GLMakie.colsize!(fig.layout, 1, GLMakie.Relative(0.38))

    GLMakie.Label(
        left[1, 1],
        "Components: click to select, right-click to mark for removal";
        font = :bold,
        halign = :left,
        tellwidth = false,
    )
    maps_grid = left[2, 1] = GLMakie.GridLayout(; valign = :top, tellheight = false)
    ncol = ceil(Int64, sqrt(ic_n))
    nrow = ceil(Int64, ic_n / ncol)
    map_axes = GLMakie.Axis[]
    for i = 1:ic_n
        ax = GLMakie.Axis(
            maps_grid[(i - 1) ÷ ncol + 1, (i - 1) % ncol + 1];
            aspect = 1,
            title = "IC $i",
            titlesize = 12,
            backgroundcolor = :transparent,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            _AXIS_LOCK_KWARGS...,
        )
        GLMakie.hidedecorations!(ax)
        GLMakie.hidespines!(ax)
        GLMakie.heatmap!(
            ax, sm.gx, sm.gy, sm.maps[i];
            colormap = :bluesreds,
            colorrange = _ica_crange(sm.maps[i]),
            nan_color = :transparent,
        )
        _draw_head_outline!(ax; lw = 1)
        GLMakie.xlims!(ax, -sm.lim, sm.lim)
        GLMakie.ylims!(ax, -sm.lim, sm.lim)
        push!(map_axes, ax)
    end
    # square cells (with room for the title)
    for r = 1:nrow
        GLMakie.rowsize!(maps_grid, r, GLMakie.Aspect(1, 1.15))
    end
    GLMakie.rowgap!(maps_grid, 4)
    GLMakie.colgap!(maps_grid, 4)

    ax_ic = GLMakie.Axis(right[1, 1]; xlabel = "Time [s]", ylabel = "Amplitude [A.U.]")
    ax_map = GLMakie.Axis(
        right[1:2, 2];
        aspect = 1,
        width = 320,
        backgroundcolor = :transparent,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.hidedecorations!(ax_map)
    GLMakie.hidespines!(ax_map)
    ax_psd = GLMakie.Axis(right[2, 1]; xlabel = "Frequency [Hz]", ylabel = "Power [dB]")
    ch_grid = right[3, 1:2] = GLMakie.GridLayout()
    ax_ch = GLMakie.Axis(ch_grid[2, 1]; xlabel = "Time [s]")

    # plots updated through observables
    ic_pts = Observable(GLMakie.Point2f[])
    psd_pts = Observable(GLMakie.Point2f[])
    orig_pts = Observable(GLMakie.Point2f[])
    clean_pts = Observable(GLMakie.Point2f[])
    big = Observable(lm.maps[1])
    big_range = Observable(_ica_crange(lm.maps[1]))
    GLMakie.lines!(ax_ic, ic_pts; color = :black, linewidth = 1)
    GLMakie.lines!(ax_psd, psd_pts; color = :black, linewidth = 1.5)
    GLMakie.lines!(ax_ch, orig_pts; color = :gray60, linewidth = 1, label = "original")
    GLMakie.lines!(
        ax_ch,
        clean_pts;
        color = :dodgerblue3,
        linewidth = 1,
        label = "marked components removed",
    )
    GLMakie.Legend(
        ch_grid[1, 1],
        ax_ch;
        orientation = :horizontal,
        framevisible = false,
        tellwidth = false,
        halign = :right,
    )
    GLMakie.rowgap!(ch_grid, 0)
    GLMakie.heatmap!(
        ax_map, lm.gx, lm.gy, big;
        colormap = :bluesreds,
        colorrange = big_range,
        nan_color = :transparent,
    )
    _draw_head_outline!(ax_map; lw = 2)
    GLMakie.scatter!(ax_map, lm.loc_x, lm.loc_y; color = :black, markersize = 4)
    GLMakie.xlims!(ax_map, -lm.lim, lm.lim)
    GLMakie.ylims!(ax_map, -lm.lim, lm.lim)

    # ---------------------------------------------------------------- #
    # controls                                                          #
    # ---------------------------------------------------------------- #
    ctrl = right[4, 1:2] = GLMakie.GridLayout(; halign = :left)
    GLMakie.Label(ctrl[1, 1], "Component"; halign = :right)
    bt_ic_prev = GLMakie.Button(ctrl[1, 2]; label = "◀", width = 36)
    lab_ic = GLMakie.Label(ctrl[1, 3], "1"; width = 40)
    bt_ic_next = GLMakie.Button(ctrl[1, 4]; label = "▶", width = 36)
    bt_mark = GLMakie.Button(ctrl[1, 5]; label = "Mark for removal", width = 150)
    bt_clear = GLMakie.Button(ctrl[1, 6]; label = "Clear marks", width = 110)
    GLMakie.Label(ctrl[1, 7], "Channel"; halign = :right)
    bt_ch_prev = GLMakie.Button(ctrl[1, 8]; label = "◀", width = 36)
    menu_ch = GLMakie.Menu(ctrl[1, 9]; options = clabels, width = 110)
    bt_ch_next = GLMakie.Button(ctrl[1, 10]; label = "▶", width = 36)

    nav = right[5, 1:2] = GLMakie.GridLayout(; halign = :left)
    GLMakie.Label(nav[1, 1], "Time"; halign = :right)
    bt_start = GLMakie.Button(nav[1, 2]; label = "⇤", width = 36)
    bt_back = GLMakie.Button(nav[1, 3]; label = "←", width = 36)
    max_start() = max(t_first, t_last - wlen_s[])
    wlen_s = Observable(Float64(min(wlen, t_last - t_first)))
    # window start positions in 0.1 s steps (at least two positions)
    time_range() = t_first:0.1:max(max_start(), t_first + 0.1)
    sl_time = GLMakie.Slider(
        nav[1, 4];
        range = time_range(),
        startvalue = t_first,
        width = 500,
    )
    bt_fwd = GLMakie.Button(nav[1, 5]; label = "→", width = 36)
    bt_end = GLMakie.Button(nav[1, 6]; label = "⇥", width = 36)
    lab_time = GLMakie.Label(nav[1, 7], ""; width = 150)
    GLMakie.Label(nav[1, 8], "Window [s]"; halign = :right)
    wlen_lab(w) = isinteger(w) ? string(Int64(w)) : string(w)
    wlen_opts = sort(unique(Float64.(vcat([1, 2, 5, 10, 20, 30, 60], wlen))))
    menu_wlen = GLMakie.Menu(
        nav[1, 9];
        options = wlen_lab.(wlen_opts),
        default = wlen_lab(Float64(wlen)),
        width = 70,
    )

    status = Observable("Keys: ← → time, ↑ ↓ component, m mark/unmark.")
    bottom = right[6, 1:2] = GLMakie.GridLayout(; halign = :left)
    bt_apply = GLMakie.Button(bottom[1, 1]; label = "Apply", width = 90)
    bt_close = GLMakie.Button(bottom[1, 2]; label = "Close", width = 90)
    GLMakie.Label(bottom[1, 3], status; halign = :left, color = :gray25)

    GLMakie.rowsize!(right, 1, GLMakie.Relative(0.27))
    GLMakie.rowsize!(right, 2, GLMakie.Relative(0.27))

    # ---------------------------------------------------------------- #
    # drawing                                                           #
    # ---------------------------------------------------------------- #
    msg!(s) = (status[] = s; nothing)
    marked_idx() = findall(marked[])

    # time window in samples
    function window()
        t0 = clamp(sl_time.value[], t_first, max_start())
        i1 = clamp(searchsortedfirst(t, t0), 1, length(t))
        i2 = clamp(i1 + round(Int64, wlen_s[] * fs) - 1, i1, length(t))
        return i1, i2
    end

    function draw_signals!()
        i1, i2 = window()
        tt = t[i1:i2]
        i = cur_ic[]
        c = cur_ch[]
        ic_pts[] = GLMakie.Point2f.(tt, ic[i, i1:i2])
        x = obj.data[chn[c], i1:i2, 1]
        m = marked_idx()
        y = isempty(m) ? x : x .- vec(ic_mw[c:c, m] * ic[m, i1:i2])
        orig_pts[] = GLMakie.Point2f.(tt, x)
        clean_pts[] = GLMakie.Point2f.(tt, y)
        for ax in (ax_ic, ax_ch)
            GLMakie.autolimits!(ax)
            GLMakie.xlims!(ax, tt[1], tt[end])
        end
        ax_ch.title[] = "Channel: $(clabels[c])" * (isempty(m) ? "" : " (removed: IC $(_v2s(m)))")
        ax_ch.ylabel[] = "Amplitude [$(units[c])]"
        lab_time.text[] = "$(round(tt[1]; digits = 2))–$(round(tt[end]; digits = 2)) s"
        return nothing
    end

    function draw_ic!()
        i = cur_ic[]
        f, p = ic_psd(i)
        psd_pts[] = GLMakie.Point2f.(f, p)
        GLMakie.autolimits!(ax_psd)
        m = big_map(i)
        big_range[] = _ica_crange(m)
        big[] = m
        ax_ic.title[] = ic_title(i) * (marked[][i] ? " — marked for removal" : "")
        ax_psd.title[] = "IC $i: power spectrum"
        lab_ic.text[] = "$i / $ic_n"
        bt_mark.label[] = marked[][i] ? "Unmark" : "Mark for removal"
        for (k, ax) in enumerate(map_axes)
            ax.backgroundcolor[] = k == i ? GLMakie.RGBAf(0.85, 0.85, 0.85, 1) : GLMakie.RGBAf(1, 1, 1, 0)
            ax.titlecolor[] = GLMakie.to_color(marked[][k] ? :red : :black)
            ax.title[] = marked[][k] ? "IC $k ✕" : "IC $k"
        end
        return nothing
    end

    function toggle!(i::Int64)
        v = copy(marked[])
        v[i] = !v[i]
        marked[] = v
        return msg!(
            isempty(marked_idx()) ? "No components marked." :
            "Marked for removal: IC $(_v2s(marked_idx())).",
        )
    end

    # ---------------------------------------------------------------- #
    # events                                                            #
    # ---------------------------------------------------------------- #
    on(_ -> draw_ic!(), cur_ic)
    on(_ -> draw_signals!(), cur_ic)
    on(_ -> (draw_ic!(); draw_signals!()), marked)
    on(_ -> draw_signals!(), cur_ch)
    on(_ -> draw_signals!(), sl_time.value)

    step_ic!(d) = (cur_ic[] = mod1(cur_ic[] + d, ic_n))
    step_ch!(d) = (cur_ch[] = mod1(cur_ch[] + d, length(chn)))
    move_time!(t0) = GLMakie.set_close_to!(sl_time, clamp(t0, t_first, max_start()))

    on(_ -> step_ic!(-1), bt_ic_prev.clicks)
    on(_ -> step_ic!(1), bt_ic_next.clicks)
    on(_ -> toggle!(cur_ic[]), bt_mark.clicks)
    on(bt_clear.clicks) do _
        marked[] = falses(ic_n)
        return msg!("Marks cleared.")
    end
    on(_ -> step_ch!(-1), bt_ch_prev.clicks)
    on(_ -> step_ch!(1), bt_ch_next.clicks)
    on(cur_ch) do c
        menu_ch.i_selected[] != c && (menu_ch.i_selected[] = c)
    end
    on(menu_ch.selection) do s
        isnothing(s) && return nothing
        c = findfirst(==(s), clabels)
        !isnothing(c) && c != cur_ch[] && (cur_ch[] = c)
        return nothing
    end

    on(_ -> move_time!(t_first), bt_start.clicks)
    on(_ -> move_time!(sl_time.value[] - wlen_s[]), bt_back.clicks)
    on(_ -> move_time!(sl_time.value[] + wlen_s[]), bt_fwd.clicks)
    on(_ -> move_time!(t_last), bt_end.clicks)
    on(menu_wlen.selection) do s
        isnothing(s) && return nothing
        wlen_s[] = min(parse(Float64, s), t_last - t_first)
        t0 = sl_time.value[]
        sl_time.range[] = time_range()
        move_time!(t0)
        draw_signals!()
        return nothing
    end

    # clicking a component map: left button selects it, right button marks it
    on(GLMakie.events(fig.scene).mousebutton; priority = 2) do ev
        ev.action == GLMakie.Mouse.press || return GLMakie.Consume(false)
        for (k, ax) in enumerate(map_axes)
            if GLMakie.is_mouseinside(ax.scene)
                if ev.button == GLMakie.Mouse.left
                    cur_ic[] = k
                elseif ev.button == GLMakie.Mouse.right
                    toggle!(k)
                end
                return GLMakie.Consume(true)
            end
        end
        return GLMakie.Consume(false)
    end

    on(GLMakie.events(fig.scene).keyboardbutton) do ev
        ev.action in (GLMakie.Keyboard.press, GLMakie.Keyboard.repeat) || return nothing
        if ev.key == GLMakie.Keyboard.left
            move_time!(sl_time.value[] - wlen_s[])
        elseif ev.key == GLMakie.Keyboard.right
            move_time!(sl_time.value[] + wlen_s[])
        elseif ev.key == GLMakie.Keyboard.up
            step_ic!(-1)
        elseif ev.key == GLMakie.Keyboard.down
            step_ic!(1)
        elseif ev.key == GLMakie.Keyboard.m && ev.action == GLMakie.Keyboard.press
            toggle!(cur_ic[])
        end
        return nothing
    end

    on(_ -> (applied[] = true), bt_apply.clicks)

    merge!(
        ctl,
        Dict(
            :ic_prev => bt_ic_prev,
            :ic_next => bt_ic_next,
            :mark => bt_mark,
            :clear => bt_clear,
            :ch_prev => bt_ch_prev,
            :ch_next => bt_ch_next,
            :menu_ch => menu_ch,
            :start => bt_start,
            :back => bt_back,
            :fwd => bt_fwd,
            :end => bt_end,
            :time => sl_time,
            :menu_wlen => menu_wlen,
            :apply => bt_apply,
            :close => bt_close,
            :status => status,
            :clean => clean_pts,
            :orig => orig_pts,
        ),
    )

    draw_ic!()
    draw_signals!()

    return (; fig, ctl, cur_ic, cur_ch, marked, applied)
end

"""
    plot_ica(obj; <keyword arguments>)

Interactive browser of independent components (ICA). For every component its scalp map (the component weights at the channels), time course and power spectrum are shown, together with a channel signal before and after removing the components marked for removal.

Components are selected by clicking their maps (or with ↑ ↓) and marked for removal by right-clicking (or with the Mark button or the M key); ← → move through the signal. Apply closes the window and returns the marked components, which can then be removed with `ica_remove()`:

```julia
ic, ic_mw, ic_var = ica_decompose(eeg; ch = "eeg")
bad = plot_ica(eeg; ch = "eeg", ic = ic, ic_mw = ic_mw, ic_var = ic_var)
ica_remove!(eeg; ch = "eeg", ic = ic, ic_mw = ic_mw, ic_idx = bad)
```

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (continuous, with channel locations)
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) used for the decomposition
- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: mixing matrix (component weights at the channels), shape (channels, n)
- `ic_var::Union{Nothing, Vector{Float64}}=nothing`: variance explained by each component (%), as returned by `ica_decompose()`
- `imethod::Symbol=:tp`: interpolation method of the scalp maps, see `plot_topo()`
- `flim::Tuple{Real, Real}=(0, min(sr(obj) / 2, 80))`: frequency limits of the power spectrum
- `wlen::Real=10`: length of the displayed signal window in seconds
- `gui::Bool=true`: if `true`, open the interactive window and wait until it is closed

# Returns

- `Vector{Int64}`: components marked for removal when the window is closed with Apply (empty if it is closed with Close), if `gui = true`
- `GLMakie.Figure`: the figure (not displayed), if `gui = false`
"""
function plot_ica(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_var::Union{Nothing, Vector{Float64}} = nothing,
    imethod::Symbol = :tp,
    flim::Tuple{Real, Real} = (0, min(sr(obj) / 2, 80)),
    wlen::Real = 10,
    gui::Bool = true,
)::Union{Vector{Int64}, GLMakie.Figure}
    pv = _plot_ica_build(
        obj;
        ch = ch,
        ic = ic,
        ic_mw = ic_mw,
        ic_var = ic_var,
        imethod = imethod,
        flim = flim,
        wlen = wlen,
    )
    gui || return pv.fig

    screen = GLMakie.Screen(; title = "NeuroAnalyzer: plot_ica()")
    on(_ -> close(screen), pv.ctl[:apply].clicks)
    on(_ -> close(screen), pv.ctl[:close].clicks)
    display(screen, pv.fig)
    GLMakie.GLFW.MaximizeWindow(screen.glscreen)
    wait(screen)

    return pv.applied[] ? findall(pv.marked[]) : Int64[]
end
