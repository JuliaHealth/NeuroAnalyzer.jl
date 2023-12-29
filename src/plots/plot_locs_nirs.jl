export plot_locs_nirs

"""
    plot_locs_nirs(locs; <keyword arguments>)

Preview of NIRS optodes and channel locations. It uses Cartesian `:loc_x` and `:loc_y` locations.

# Arguments

- `locs::DataFrame`: columns: labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `opt_pairs::Matrix{Int64}`: pairs of source and detector
- `src_n::Int64`: number of sources
- `det_n::Int64`: number of detectors
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
- `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
- `src_labels::Bool=false`: plot source labels
- `det_labels::Bool=false`: plot detector labels
- `opt_labels::Bool=false`: plot optode type (S for source, D for detector) and number
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: Use color or gray palette
- `grid::Bool=false`: draw grid, useful for locating positions
- `plot_size::Int64=400`: plot dimensions in pixels (size × size)

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs_nirs(locs::DataFrame, opt_pairs::Matrix{Int64}, src_n::Int64, det_n::Int64; src_labels::Bool=false, det_labels::Bool=false, opt_labels::Bool=false, ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1, selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0, head::Bool=true, head_labels::Bool=true, mono::Bool=false, grid::Bool=false, plot_size::Int64=400, large::Bool=true)

    # TO DO: plot channel numbers

    if plot_size > 400
        marker_size = plot_size ÷ 100
        font_size = plot_size ÷ 50
    else
        marker_size = plot_size ÷ 200
        font_size = plot_size ÷ 100
    end

    if large
        marker_size = 4
        font_size = 6
        # loc_x = @. round(origin[1] + (loc_x * 250), digits=2)
        # loc_y = @. round(origin[2] - (loc_y * 250), digits=2)
    else
        marker_size = 4
        font_size = 4
        src_labels = false
        det_labels = false
        grid = false
        # loc_x = @. round(origin[1] + (loc_x * 100), digits=2)
        # loc_y = @. round(origin[2] - (loc_y * 100), digits=2)
    end

    if grid
        p = Plots.plot(grid=true,
                       xlim=(-1.2, 1.2),
                       ylim=(-1.2, 1.2),
                       ratio=1,
                       legend=false,
                       xticks=-1.2:0.1:1.2,
                       yticks=-1.2:0.1:1.2,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       right_margin=-20*Plots.px,
                       bottom_margin=-10*Plots.px,
                       top_margin=-20*Plots.px,
                       left_margin=-5*Plots.px,
                       size=(plot_size, plot_size))
    else
        p = Plots.plot(border=:none,
                       grid=false,
                       xlim=(-1.2, 1.2),
                       ylim=(-1.2, 1.2),
                       ratio=1,
                       legend=false,
                       xticks=-1.2:0.1:1.2,
                       yticks=-1.2:0.1:1.2,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       right_margin=-20*Plots.px,
                       bottom_margin=-10*Plots.px,
                       top_margin=-20*Plots.px,
                       left_margin=-5*Plots.px,
                       size=(plot_size, plot_size))
    end

    x = locs[!, :loc_x]
    y = locs[!, :loc_y]

    for idx in 1:size(opt_pairs, 1)
        xs = x[opt_pairs[idx, 1]]
        xd = x[src_n + opt_pairs[idx, 2]]
        ys = y[opt_pairs[idx, 1]]
        yd = y[src_n + opt_pairs[idx, 2]]
        if mono == true
            p = Plots.plot!([xs, xd], [ys, yd], lc=:gray, lw=1, alpha=0.5)
        else
            p = Plots.plot!([xs, xd], [ys, yd], lc=:blue, lw=1, alpha=0.5)
        end
    end
    
    if src_labels
        for idx in 1:src_n
            p = Plots.plot!(annotations=(x[idx], y[idx], Plots.text(locs[!, :labels][idx], pointsize=font_size)))
        end
    elseif opt_labels == false
        if mono == true
            p = Plots.scatter!(x[1:src_n], y[1:src_n], c=:black, msc=:black, ms=marker_size, msa=1)
        else
            p = Plots.scatter!(x[1:src_n], y[1:src_n], c=:red, msc=:red, ms=marker_size, msa=1)
        end
    end
    
    if det_labels
        for idx in (src_n + 1):(src_n + det_n)
            p = Plots.plot!(annotations=(x[idx], y[idx], Plots.text(locs[!, :labels][idx], pointsize=font_size)))
        end
    elseif opt_labels == false
        if mono == true
            p = Plots.scatter!(x[(src_n + 1):end], y[(src_n + 1):end], c=:white, msc=:black, ms=marker_size, msa=1)
        else
            p = Plots.scatter!(x[(src_n + 1):end], y[(src_n + 1):end], c=:green, msc=:green, ms=marker_size, msa=1)
        end
    end
    
    if opt_labels
        for idx in 1:src_n
            p = Plots.plot!(annotations=(x[idx], y[idx], Plots.text("S" * string(idx), pointsize=font_size)))
        end
        # for idx in (src_n + 1):(src_n + det_n)
        for idx in 1:det_n
            p = Plots.plot!(annotations=(x[src_n + idx], y[src_n + idx], Plots.text("D" * string(idx), pointsize=font_size)))
        end
    end

    if head == true
        _warn("TO DO: add head outline.")
    end

    p = Plots.plot!()

    return p

end