# export plot_locs_ecog

"""
    plot_locs_ecog(locs; <keyword arguments>)

Draw ECOG grid electrode and channel numbers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `grid::Bool=false`: draw grid, useful for locating positions
- `mono::Bool=false`: use color or grey palette
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs_nirs(obj::NeuroAnalyzer.NEURO; grid::Bool=false, mono::Bool=false, plot_size::Int64=800)

    _check_datatype(obj, :ecog)

    if plot_size > 400
        marker_size = plot_size ÷ 100
        font_size = plot_size ÷ 50
    else
        marker_size = plot_size ÷ 200
        font_size = plot_size ÷ 100
    end

    if grid == true
        p = Plots.plot(grid=true,
                       xlim=(-1.22, 1.23),
                       ylim=(-1.1, 1.2),
                       ratio=1,
                       legend=false,
                       xticks=-1:0.1:1,
                       yticks=-1:0.1:1,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       right_margin=-20*Plots.px,
                       bottom_margin=-10*Plots.px,
                       top_margin=-20*Plots.px,
                       left_margin=-10*Plots.px,
                       size=(plot_size, plot_size))
    else
        p = Plots.plot(border=:none,
                       grid=false,
                       # xlim=(-1.22, 1.23),
                       # ylim=(-1.1, 1.2),
                       ratio=1,
                       legend=false,
                       right_margin=-20*Plots.px,
                       bottom_margin=-10*Plots.px,
                       top_margin=-20*Plots.px,
                       left_margin=-10*Plots.px,
                       size=(plot_size, plot_size))
    end

#=
    locs = Matrix{Tuple{Real, Real}}(undef, 4, 4)
    x = linspace(-1, 1, size(locs, 1))
    y = linspace(-1, 1, size(locs, 2))
    for idx_x in 1:size(locs, 1), idx_y in 1:size(locs, 2)
        locs[idx_x, idx_y] = (x[idx_x], y[idx_y])
        Plots.scatter!(p, locs[idx_x, idx_y])
    end
=#

    p = Plots.plot!()

    return p

end