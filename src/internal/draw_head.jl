function _draw_head(p::Plots.Plot{Plots.GRBackend}; head_labels::Bool=true, head_details::Bool=true, topo::Bool=false, kwargs...)
    _deprecated("_draw_head")
    # Draw head over a topographical plot `p`.
    # - `p::Plots.Plot{Plots.GRBackend}`: electrodes plot
    # - `loc_x::Vector{<:Real}`: vector of x electrode position
    # - `loc_y::Vector{<:Real}`: vector of y electrode position
    # - `head_labels::Bool=true`: add text labels to the plot
    # - `topo::Bool=false`: if true, perform peripheral erasing for topo plots
    # - `kwargs`: optional arguments for plot() function
    # loc_x, loc_y = loc_y, loc_x
    pts = Plots.partialcircle(0, 2π, 100, 1.1)
    x, y = Plots.unzip(pts)
    maxx = maximum(x)
    maxy = maximum(y)
    minx = minimum(x)
    miny = minimum(y)
    head = Plots.Shape(x, y)
    if head_details
        nose = Plots.Shape([(-0.2, maxy - 0.015),
                            (0, maxy + 0.08),
                            (0.2, maxy - 0.015),
                            (-0.01, maxy)
                           ])
        ear_r = Plots.Shape([(maxx, -0.05),
                             (maxx - 0.005, 0.09),
                             (maxx + 0.02, 0.125),
                             (maxx + 0.04, 0.13),
                             (maxx + 0.06, 0.115),
                             (maxx + 0.075, 0.085),
                             (maxx + 0.07, 0),
                             (maxx + 0.08, -0.155),
                             (maxx + 0.05, -0.215),
                             (maxx + 0.015, -0.225),
                             (maxx - 0.016, -0.19)
                            ])
        ear_l = Plots.Shape([(minx, -0.05),
                             (minx + 0.005, 0.09),
                             (minx - 0.02, 0.125),
                             (minx - 0.04, 0.13),
                             (minx - 0.06, 0.115),
                             (minx - 0.075, 0.085),
                             (minx - 0.07, 0),
                             (minx - 0.08, -0.155),
                             (minx - 0.05, -0.215),
                             (minx - 0.015, -0.225),
                             (minx + 0.016, -0.19)
                            ])
        p = Plots.plot!(nose, fill=nothing, label="")
        p = Plots.plot!(ear_l, fill=nothing, label="")
        p = Plots.plot!(ear_r, fill=nothing, label="")
    end

    p = Plots.plot!(p, head, fill=nothing, label="")

    if head_labels
        p = Plots.plot!(annotations=(0, 1.05, Plots.text("IN", pointsize=4, halign=:center, valign=:center)))
        p = Plots.plot!(annotations=(0, -1.05, Plots.text("NAS", pointsize=4, halign=:center, valign=:center)))
        p = Plots.plot!(annotations=(-1.05, 0, Plots.text("LPA", pointsize=4, halign=:center, valign=:center, rotation=90)))
        p = Plots.plot!(annotations=(1.05, 0, Plots.text("RPA", pointsize=4, halign=:center, valign=:center, rotation=-90)))
    end

    if topo
        pts = Plots.partialcircle(0, 2π, 100, 1.4)
        x, y = Plots.unzip(pts)
        for idx in 1:0.001:1.7
            peripheral = Shape(x .* idx, y .* idx)
            p = Plots.plot!(p, peripheral, label="", fill=nothing, lc=:white)
        end
        p = Plots.plot!(xlims=(-1.4, 1.4), ylims=(-1.4, 1.4); kwargs...)
    end

    return p
end

