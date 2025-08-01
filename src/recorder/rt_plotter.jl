export rt_plotter

"""
    rt_plotter(; <keyword arguments>)

Plot recorded signal in real time.

# Arguments

- `fs::Int64`
- `l::Int64=5`: displayed segment length (in seconds)
- `duration::Int64=20`: duration of recording (in seconds)

# Returns

- `Plots.Plot{Plots.GRBackend}`
"""
function rt_plotter(; fs::Int64, l::Int64=5, duration::Int64=20)::Plots.Plot{Plots.GRBackend}

    _wip()

    t = 0:1/fs:l
    y = zeros(length(t))
    t_idx = 1
    seg_idx = 0
    sleep(0.1)
    t_s = time()

    while true
        if time() >= t_s + duration
            break
        else
            t_refresh = time()
        end
        while true
            if time() >= (t_refresh + 1/fs)
                y[t_idx] = rand(-1:0.1:1, 1)[1]
                t_idx += 1
                if t_idx > length(y)
                    t_idx = 1
                    seg_idx += 1
                end
            end
            if time() >= (t_refresh + 1/fs)
                t = (l * seg_idx):1/fs:(l * (seg_idx + 1))
                p = Plots.plot(t,
                               y,
                               ylims=(-1, 1),
                               xticks=(t[1]:t[end]),
                               legend=false,
                               title="dsfdsf",
                               palette=:darktest,
                               size=(800, 400),
                               margins=20Plots.px,
                               titlefontsize=10,
                               xlabelfontsize=8,
                               ylabelfontsize=8,
                               xtickfontsize=8,
                               ytickfontsize=8)
                display(p)
                break
            end
        end
    end

end