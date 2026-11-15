export rt_plotter

# window of rt_plotter(): one axis showing a segment of `l` seconds
function _rt_plotter_build(fs::Int64, l::Int64, title::String)
    n = l * fs # samples per displayed segment
    fig = GLMakie.Figure(; size = (900, 450))
    ax = GLMakie.Axis(
        fig[1, 1];
        title = title,
        xlabel = "Time [s]",
        ylabel = "Amplitude",
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, -1, 1)
    GLMakie.xlims!(ax, 0, l)
    pts = Observable(GLMakie.Point2f[])
    GLMakie.lines!(ax, pts; color = :black, linewidth = 1)

    # show the first k samples of segment seg_idx (0-based)
    function show!(y::Vector{Float64}, k::Int64, seg_idx::Int64)
        t0 = seg_idx * l
        pts[] = GLMakie.Point2f.(t0 .+ (0:(k - 1)) ./ fs, y[1:k])
        GLMakie.xlims!(ax, t0, t0 + l)
        return nothing
    end

    return (; fig, ax, pts, n, show!)
end

"""
    rt_plotter(; <keyword arguments>)

Plot a recorded signal in real time.

The displayed window of `l` seconds is filled from left to right and then starts again. Plotting stops after `duration` seconds or when the window is closed. Work in progress: the signal is a random placeholder (no acquisition source yet).

# Arguments

- `fs::Int64`: sampling rate in Hz; must be > 0
- `l::Int64=5`: displayed segment length in seconds; must be > 0
- `duration::Int64=20`: duration of recording in seconds; must be > 0
- `title::String="Real-time signal"`: plot title

# Returns

- `GLMakie.Figure`: the plot (showing the last segment)
"""
function rt_plotter(;
    fs::Int64,
    l::Int64 = 5,
    duration::Int64 = 20,
    title::String = "Real-time signal",
)::GLMakie.Figure
    _wip()

    # validate
    fs > 0 || throw(ArgumentError("fs must be > 0."))
    l > 0 || throw(ArgumentError("l must be > 0."))
    duration > 0 || throw(ArgumentError("duration must be > 0."))

    w = _rt_plotter_build(fs, l, title)
    y = fill(NaN, w.n)
    seg_idx = 0
    # the plot is refreshed at most 10 times per second
    t_refresh_every = 0.1

    screen = GLMakie.Screen(; title = "NeuroRecorder: rt_plotter()")
    display(screen, w.fig)

    t_s = time()
    t_refresh = t_s
    n_sample = 0
    while time() < t_s + duration && isopen(screen)
        # new samples due since the start (placeholder: random values)
        while n_sample < floor(Int64, (time() - t_s) * fs)
            idx = n_sample % w.n + 1
            if idx == 1 && n_sample > 0
                seg_idx += 1
                fill!(y, NaN)
            end
            y[idx] = rand(-1:0.1:1)
            n_sample += 1
        end
        if time() - t_refresh >= t_refresh_every
            w.show!(y, (n_sample - 1) % w.n + 1, seg_idx)
            t_refresh = time()
        end
        sleep(1 / (2 * fs))
    end

    return w.fig
end
