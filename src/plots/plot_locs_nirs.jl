export plot_locs_nirs

"""
    plot_locs_nirs(locs; <keyword arguments>)

Preview of NIRS optodes and channel locations. It uses Cartesian `:loc_x` and `:loc_y` locations.

# Arguments

  - `locs::DataFrame`: columns: labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
  - `opt_pairs::Matrix{Int64}`: pairs of source and detector
  - `src_n::Int64`: number of sources
  - `det_n::Int64`: number of detectors
  - `src_labels::Bool=false`: plot source labels
  - `det_labels::Bool=false`: plot detector labels
  - `opt_labels::Bool=false`: plot optode type (S for source, D for detector) and number
  - `head::Bool=true`: draw head
  - `head_labels::Bool=false`: plot head labels
  - `mono::Bool=false`: use color or gray palette
  - `grid::Bool=false`: draw grid, useful for locating positions
  - `ps::Symbol=:l`: plot size (`:l`: large (800×800 px), `:m`: medium (300×300 px), `:s`: small (100×100 px))
  - `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  - `plane::Symbol=:xy`: which plane to plot:
      + `:xy`: horizontal (top)
      + `:xz`: coronary (front)
      + `:yz`: sagittal (side)
  - `ch_info::Vector{String}=string.(1:DataFrames.nrow(locs))`: channels info

# Returns

  - `p::GLMakie.Figure`
"""
function plot_locs_nirs(
        locs::DataFrame,
        opt_pairs::Matrix{Int64},
        src_n::Int64,
        det_n::Int64;
        src_labels::Bool = false,
        det_labels::Bool = false,
        opt_labels::Bool = false,
        head::Bool = true,
        head_labels::Bool = true,
        mono::Bool = false,
        grid::Bool = false,
        ps::Symbol = :l,
        cart::Bool = false,
        plane::Symbol = :xy,
        ch_info::Vector{String} = string.(1:DataFrames.nrow(locs)),
    )::GLMakie.Figure

    # TO DO: plot channel numbers

    ch = 1:DataFrames.nrow(locs)

    _check_var(ps, [:l, :m, :s], "ps")
    _check_var(plane, [:xy, :yz, :xz], "plane")
    pal = mono ? :grays : :darktest

    if plane === :xy
        if !cart
            loc_x = zeros(length(ch))
            loc_y = zeros(length(ch))
            for idx in 1:length(ch)
                loc_x[idx], loc_y[idx] = pol2cart(locs[ch, :loc_radius][idx], locs[ch, :loc_theta][idx])
            end
        else
            loc_x = locs[ch, :loc_x]
            loc_y = locs[ch, :loc_y]
        end
    elseif plane === :xz
        if !cart
            loc_x = zeros(length(ch))
            loc_y = zeros(length(ch))
            for idx in 1:length(ch)
                loc_x[idx], _, loc_y[idx] = sph2cart(
                    locs[ch, :loc_radius_sph][idx], locs[ch, :loc_theta_sph][idx], locs[ch, :loc_phi_sph][idx]
                )
            end
        else
            loc_x = locs[ch, :loc_x]
            loc_y = locs[ch, :loc_z]
        end
    elseif plane === :yz
        if !cart
            loc_x = zeros(length(ch))
            loc_y = zeros(length(ch))
            for idx in 1:length(ch)
                _, loc_x[idx], loc_y[idx] = sph2cart(
                    locs[ch, :loc_radius_sph][idx], locs[ch, :loc_theta_sph][idx], locs[ch, :loc_phi_sph][idx]
                )
            end
        else
            loc_x = locs[ch, :loc_y]
            loc_y = locs[ch, :loc_z]
        end
    end

    xl = (-1.2, 1.2)
    yl = (-1.2, 1.2)

    if ps === :l
        plot_size = (800, 800)
        marker_size = length(ch) > 64 ? 10 : 20
        font_size = 14
        length(ch) > 64 && (ch_labels = false)
    elseif ps === :m
        plot_size = (300, 300)
        marker_size = length(ch) > 64 ? 5 : 10
        font_size = 8
        src_labels = false
        det_labels = false
        grid = false
    elseif ps === :s
        plot_size = (100, 100)
        marker_size = length(ch) > 64 ? 4 : 8
        font_size = 8
        head_labels = false
        src_labels = false
        det_labels = false
        grid = false
    end

    # prepare plot
    GLMakie.activate!(title = "plot_locs_nirs()")
    p = GLMakie.Figure(
        size = plot_size,
        figure_padding = 0,
    )
    if grid
        ax = GLMakie.Axis(
            p[1, 1];
            aspect = 1,
            xlabel = "",
            ylabel = "",
            title = "",
            xticks = xt,
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(2),
            yticks = yt,
            yminorticksvisible = true,
            yminorticks = IntervalsBetween(2),
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            backgroundcolor = :transparent,
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
    else
        ax = GLMakie.Axis(
            p[1, 1];
            aspect = 1,
            xlabel = "",
            ylabel = "",
            title = "",
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            backgroundcolor = :transparent,
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        hidedecorations!(ax; grid = true)
        hidespines!(ax)
    end
    GLMakie.xlims!(ax, xl)
    GLMakie.ylims!(ax, yl)

    # draw head
    if head
        ps === :l && (lw = 3)
        ps === :m && (lw = 2)
        ps === :s && (lw = 1)
        if plane === :xy
            # nose
            GLMakie.lines!(ax, [-0.2, 0], [0.98, 1.08]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [0.2, 0], [0.98, 1.08]; linewidth = lw, color = :black)

            # ears
            # left
            GLMakie.lines!(ax, [-0.995, -1.03], [0.1, 0.15]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.03, -1.06], [0.15, 0.16]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.06, -1.1], [0.16, 0.14]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.1, -1.12], [0.14, 0.05]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.12, -1.1], [0.05, -0.1]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.1, -1.13], [-0.1, -0.3]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.13, -1.09], [-0.3, -0.37]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.09, -1.02], [-0.37, -0.39]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.02, -0.98], [-0.39, -0.33]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-0.98, -0.975], [-0.33, -0.22]; linewidth = lw, color = :black)
            # right
            GLMakie.lines!(ax, [0.995, 1.03], [0.1, 0.15]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.03, 1.06], [0.15, 0.16]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.06, 1.1], [0.16, 0.14]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.1, 1.12], [0.14, 0.05]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.12, 1.1], [0.05, -0.1]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.1, 1.13], [-0.1, -0.3]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.13, 1.09], [-0.3, -0.37]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.09, 1.02], [-0.37, -0.39]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.02, 0.98], [-0.39, -0.33]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [0.98, 0.975], [-0.33, -0.22]; linewidth = lw, color = :black)

            # head
            GLMakie.arc!(ax, (0, 0), 1, 0, 2pi; linewidth = lw, color = :black)
        elseif plane === :yz
            # head
            GLMakie.arc!(ax, (0, 0), 1, 0, pi; linewidth = lw, color = :black)
        elseif plane === :xz
            # head
            GLMakie.arc!(ax, (0, 0), 1, 0, pi; linewidth = lw, color = :black)
        end
    end

    ch_n = length(ch)
    cmap = GLMakie.resample_cmap(pal, ch_n)

    ps === :l && (sw = 2)
    ps === :m && (sw = 1)
    ps === :s && (sw = 0.5)

    for idx in axes(opt_pairs, 1)
        xs = loc_x[opt_pairs[idx, 1]]
        xd = loc_x[src_n + opt_pairs[idx, 2]]
        ys = loc_y[opt_pairs[idx, 1]]
        yd = loc_y[src_n + opt_pairs[idx, 2]]
        GLMakie.lines!([xs, xd], [ys, yd]; color = mono ? :gray : :blue, alpha = 0.5)
    end

    label_offset_x = 0.0
    label_offset_y = -0.08

    if src_labels
        for idx in 1:src_n
            GLMakie.text!(
                loc_x[idx] + label_offset_x,
                loc_y[idx] + label_offset_y;
                text = locs[!, :label][idx],
                align = (:center, :bottom),
                fontsize = font_size,
            )
        end
    elseif !opt_labels
        GLMakie.scatter!(
            loc_x[1:src_n],
            loc_y[1:src_n];
            markersize = marker_size,
            color = mono ? :black : :red,
            strokewidth = sw,
            strokecolor = :black,
        )
    end

    if det_labels
        for idx in (src_n + 1):(src_n + det_n)
            GLMakie.text!(
                loc_x[idx] + label_offset_x,
                loc_y[idx] + label_offset_y;
                text = locs[!, :label][idx],
                align = (:center, :bottom),
                fontsize = font_size,
            )
        end
    elseif !opt_labels
        GLMakie.scatter!(
            loc_x[(src_n + 1):end],
            loc_y[(src_n + 1):end];
            markersize = marker_size,
            color = mono ? :white : :green,
            strokewidth = sw,
            strokecolor = :black,
        )
    end

    if opt_labels
        for idx in 1:src_n
            GLMakie.text!(
                loc_x[idx] + label_offset_x,
                loc_y[idx] + label_offset_y;
                text = "S" * string(idx),
                align = (:center, :bottom),
                fontsize = font_size,
            )
        end
        for idx in 1:det_n
            GLMakie.text!(
                loc_x[idx] + label_offset_x,
                loc_y[idx] + label_offset_y;
                text = "D" * string(idx),
                align = (:center, :bottom),
                fontsize = font_size,
            )
        end
    end

    if head_labels
        fid_names = ["NAS", "IN", "LPA", "RPA"]
        for idx in 1:length(NeuroAnalyzer.fiducial_points)
            if plane === :xy
                fid_loc_x = NeuroAnalyzer.fiducial_points[idx][1]
                fid_loc_y = NeuroAnalyzer.fiducial_points[idx][2]
            elseif plane === :xz
                fid_loc_x = NeuroAnalyzer.fiducial_points[idx][1]
                fid_loc_y = NeuroAnalyzer.fiducial_points[idx][3]
            elseif plane === :yz
                fid_loc_x = NeuroAnalyzer.fiducial_points[idx][2]
                fid_loc_y = NeuroAnalyzer.fiducial_points[idx][3]
            end
            GLMakie.text!(fid_loc_x, fid_loc_y; text = fid_names[idx], fontsize = font_size, align = (:center, :center))
        end
    end

    return p

end
