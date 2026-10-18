export rt_plotter

"""
    rt_plotter(; <keyword arguments>)

Plot a recorded signal in real time.

The displayed window of `l` seconds is filled from left to right and then starts again. Work in progress: the signal is a random placeholder (no acquisition source yet).

# Arguments

- `fs::Int64`: sampling rate in Hz; must be > 0
- `l::Int64=5`: displayed segment length in seconds; must be > 0
- `duration::Int64=20`: duration of recording in seconds; must be > 0
- `title::String="Real-time signal"`: plot title

# Returns

- `Plots.Plot{Plots.GRBackend}`: last displayed plot
"""
function rt_plotter(;
    fs::Int64,
    l::Int64 = 5,
    duration::Int64 = 20,
    title::String = "Real-time signal",
)::Plots.Plot{Plots.GRBackend}
    _wip()

    # validate
    fs > 0 || throw(ArgumentError("fs must be > 0."))
    l > 0 || throw(ArgumentError("l must be > 0."))
    duration > 0 || throw(ArgumentError("duration must be > 0."))

    n = l * fs # samples per displayed segment
    y = fill(NaN, n)
    seg_idx = 0
    # the plot is refreshed at most 10 times per second
    t_refresh_every = 0.1

    _rt_plot() = Plots.plot(
        (seg_idx * l) .+ (0:(n - 1)) ./ fs,
        y;
        ylims = (-1, 1),
        xlims = (seg_idx * l, (seg_idx + 1) * l),
        legend = false,
        title = title,
        palette = :darktest,
        size = (800, 400),
        margins = 20Plots.PlotMeasures.px,
        titlefontsize = 10,
        xlabelfontsize = 8,
        ylabelfontsize = 8,
        xtickfontsize = 8,
        ytickfontsize = 8,
    )

    p = _rt_plot()
    display(p)

    t_s = time()
    t_refresh = t_s
    n_sample = 0
    while time() < t_s + duration
        # new samples due since the start (placeholder: random values)
        while n_sample < floor(Int64, (time() - t_s) * fs)
            idx = n_sample % n + 1
            if idx == 1 && n_sample > 0
                seg_idx += 1
                fill!(y, NaN)
            end
            y[idx] = rand(-1:0.1:1)
            n_sample += 1
        end
        if time() - t_refresh >= t_refresh_every
            p = _rt_plot()
            display(p)
            t_refresh = time()
        end
        sleep(1 / (2 * fs))
    end

    return p
end
