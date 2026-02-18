export plot_locs

"""
    plot_locs(locs; <keyword arguments>)

Preview channel locations.

# Arguments

  - `locs::DataFrame`: columns: `channel`, `labels`, `loc_radius`, `loc_theta`, `loc_x`, `loc_y`, `loc_z`, `loc_radius_sph`, `loc_theta_sph`, `loc_phi_sph`
  - `ch::Union{Int64, Vector{Int64}, AbstractRange}=1:DataFrames.nrow(locs)`: list of locations to plot, default is all locations
  - `sch::Union{Int64, Vector{Int64}, AbstractRange}=0`: which channels are selected
  - `ch_labels::Bool=true`: plot locations labels
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
  - `connections::Matrix{<:Real}=[0 0; 0 0]`: matrix of connections weights (channels by channels)
  - `threshold::Real=0`: threshold for plotting, see below
  - `threshold_type::Symbol=:neq`: rule for thresholding:
      + `:eq`: draw region is values are equal to threshold
      + `:neq`: draw region is values are not equal to threshold
      + `:geq`: draw region is values are ≥ to threshold
      + `:leq`: draw region is values are ≤ to threshold
      + `:g`: draw region is values are > to threshold
      + `:l`: draw region is values are < to threshold
  - `weights::Union{Bool, Vector{<:Real}}=true`: weight line widths and alpha based on connection value, if false connections values will be drawn or vector of weights
  - `ch_info::Vector{String}=string.(1:DataFrames.nrow(locs))`: channels info

# Returns

  - `p::GLMakie.Figure`
"""
function plot_locs(
    locs::DataFrame;
    ch::Union{Int64, Vector{Int64}, AbstractRange} = 1:DataFrames.nrow(locs),
    sch::Union{Int64, Vector{Int64}, AbstractRange} = 0,
    ch_labels::Bool = true,
    head::Bool = true,
    head_labels::Bool = false,
    mono::Bool = false,
    grid::Bool = false,
    ps::Symbol = :l,
    cart::Bool = false,
    plane::Symbol = :xy,
    connections::Matrix{<:Real} = [0 0; 0 0],
    threshold::Real = 0,
    threshold_type::Symbol = :neq,
    weights::Union{Bool, Vector{<:Real}} = true,
    ch_info::Vector{String} = string.(1:DataFrames.nrow(locs)),
    gui::Bool = true,
)::GLMakie.Figure

    _check_var(ps, [:l, :m, :s], "ps")
    _check_var(plane, [:xy, :yz, :xz], "plane")
    pal = mono ? :grays : :darktest
    sch_labels = ch_labels

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

    loc_x = _n2v(loc_x)
    loc_y = _n2v(loc_y)

    head12 = false
    maximum(abs.(locs[:, :loc_x])) <= 1.2 &&
        maximum(abs.(locs[:, :loc_y])) <= 1.2 &&
        maximum(abs.(locs[:, :loc_z])) <= 1.5 &&
        (head12 = true)

    if head12
        xl = (-1.2, 1.2)
        yl = (-1.2, 1.2)
    else
        xl = (-1.6, 1.6)
        yl = (-1.6, 1.6)
    end

    if ps === :l
        plot_size = (800, 800)
        marker_size = length(ch) > 64 ? 10 : 20
        font_size = 14
        length(ch) > 64 && (ch_labels = false)
    elseif ps === :m
        plot_size = (300, 300)
        marker_size = length(ch) > 64 ? 5 : 10
        font_size = 8
        ch_labels = false
        sch_labels = false
        grid = false
    elseif ps === :s
        plot_size = (100, 100)
        marker_size = length(ch) > 64 ? 4 : 8
        font_size = 8
        head_labels = false
        ch_labels = false
        sch_labels = false
        grid = false
    end

    # prepare plot
    GLMakie.activate!(title = "plot_locs()")
    p = GLMakie.Figure(
                    size = plot_size,
                    figure_padding = (0, 0, 0, 0),
                ) # L R B T
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
            GLMakie.lines!(ax, [-0.2, 0], [0.980, 1.08]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [0.2, 0], [0.980, 1.08]; linewidth = lw, color = :black)

            # ears
            # left
            GLMakie.lines!(ax, [-0.995, -1.03], [0.1, 0.15]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.03, -1.06], [0.15, 0.16]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.06, -1.1], [0.16, 0.14]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.1, -1.12], [0.14, 0.05]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.12, -1.10], [0.05, -0.1]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.10, -1.13], [-0.1, -0.3]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.13, -1.09], [-0.3, -0.37]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.09, -1.02], [-0.37, -0.39]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-1.02, -0.98], [-0.39, -0.33]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [-0.98, -0.975], [-0.33, -0.22]; linewidth = lw, color = :black)
            # right
            GLMakie.lines!(ax, [0.995, 1.03], [0.1, 0.15]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.03, 1.06], [0.15, 0.16]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.06, 1.1], [0.16, 0.14]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.1, 1.12], [0.14, 0.05]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.12, 1.10], [0.05, -0.1]; linewidth = lw, color = :black)
            GLMakie.lines!(ax, [1.10, 1.13], [-0.1, -0.3]; linewidth = lw, color = :black)
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

    # draw connections
    if connections != [0 0; 0 0]
        sch = ""
        @assert size(connections, 1) == length(ch) "Length of channel and number of connections rows must be equal."
        _check_var(threshold_type, [:eq, :neq, :geq, :leq, :g, :l, :in, :bin], "threshold_type")
        if threshold_type in [:eq, :neq, :geq, :leq, :g, :l]
            @assert length(threshold) == 1 "threshold must contain a single value."
        else
            @assert length(threshold) == 2 "threshold must contain two values."
            _check_tuple(threshold, "threshold")
        end
        m_tmp = normalize_n(abs.(connections))
        for idx1 in axes(connections, 1)
            for idx2 in 2:size(connections, 1)
                if idx1 != idx2
                    if threshold_type === :g
                        if connections[idx1, idx2] > threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :red,
                                        )
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                            linestyle = :dot,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :blue,
                                        )
                                    end
                                end
                            else
                                GLMakie.lines!(
                                    [loc_x[idx1], loc_x[idx2]],
                                    [loc_y[idx1], loc_y[idx2]];
                                    linewidth = 0.2,
                                    color = :black,
                                )
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(
                                        l_pos[1], l_pos[2]; text = string(connections[idx1, idx2]), fontsize = font_size
                                    )
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :red,
                                        )
                                    else
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :blue,
                                        )
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :l
                        if connections[idx1, idx2] < threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :red,
                                        )
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                            linestyle = :dot,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :blue,
                                        )
                                    end
                                end
                            else
                                GLMakie.lines!(
                                    [loc_x[idx1], loc_x[idx2]],
                                    [loc_y[idx1], loc_y[idx2]];
                                    linewidth = 0.2,
                                    color = :black,
                                )
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(
                                        l_pos[1], l_pos[2]; text = string(connections[idx1, idx2]), fontsize = font_size
                                    )
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :red,
                                        )
                                    else
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :blue,
                                        )
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :eq
                        if connections[idx1, idx2] == threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :red,
                                        )
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                            linestyle = :dot,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :blue,
                                        )
                                    end
                                end
                            else
                                GLMakie.lines!(
                                    [loc_x[idx1], loc_x[idx2]],
                                    [loc_y[idx1], loc_y[idx2]];
                                    linewidth = 0.2,
                                    color = :black,
                                )
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(
                                        l_pos[1], l_pos[2]; text = string(connections[idx1, idx2]), fontsize = font_size
                                    )
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :red,
                                        )
                                    else
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :blue,
                                        )
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :neq
                        if connections[idx1, idx2] != threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :red,
                                        )
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                            linestyle = :dot,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :blue,
                                        )
                                    end
                                end
                            else
                                GLMakie.lines!(
                                    [loc_x[idx1], loc_x[idx2]],
                                    [loc_y[idx1], loc_y[idx2]];
                                    linewidth = 0.2,
                                    color = :black,
                                )
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(
                                        l_pos[1], l_pos[2]; text = string(connections[idx1, idx2]), fontsize = font_size
                                    )
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :red,
                                        )
                                    else
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :blue,
                                        )
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :leq
                        if connections[idx1, idx2] <= threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :red,
                                        )
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                            linestyle = :dot,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :blue,
                                        )
                                    end
                                end
                            else
                                GLMakie.lines!(
                                    [loc_x[idx1], loc_x[idx2]],
                                    [loc_y[idx1], loc_y[idx2]];
                                    linewidth = 0.2,
                                    color = :black,
                                )
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(
                                        l_pos[1], l_pos[2]; text = string(connections[idx1, idx2]), fontsize = font_size
                                    )
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :red,
                                        )
                                    else
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :blue,
                                        )
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :geq
                        if connections[idx1, idx2] >= threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :red,
                                        )
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                            linestyle = :dot,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :blue,
                                        )
                                    end
                                end
                            else
                                GLMakie.lines!(
                                    [loc_x[idx1], loc_x[idx2]],
                                    [loc_y[idx1], loc_y[idx2]];
                                    linewidth = 0.2,
                                    color = :black,
                                )
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(
                                        l_pos[1], l_pos[2]; text = string(connections[idx1, idx2]), fontsize = font_size
                                    )
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :red,
                                        )
                                    else
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :blue,
                                        )
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :in
                        if connections[idx1, idx2] >= threshold[1] && connections[idx1, idx2] <= threshold[2]
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :red,
                                        )
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                            linestyle = :dot,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :blue,
                                        )
                                    end
                                end
                            else
                                GLMakie.lines!(
                                    [loc_x[idx1], loc_x[idx2]],
                                    [loc_y[idx1], loc_y[idx2]];
                                    linewidth = 0.2,
                                    color = :black,
                                )
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(
                                        l_pos[1], l_pos[2]; text = string(connections[idx1, idx2]), fontsize = font_size
                                    )
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :red,
                                        )
                                    else
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :blue,
                                        )
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :bin
                        if connections[idx1, idx2] > threshold[1] && connections[idx1, idx2] < threshold[2]
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :red,
                                        )
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :black,
                                            linestyle = :dot,
                                        )
                                    else
                                        GLMakie.lines!(
                                            [loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]];
                                            linewidth = 6 * m_tmp[idx1, idx2],
                                            alpha = 0.25 * m_tmp[idx1, idx2],
                                            color = :blue,
                                        )
                                    end
                                end
                            else
                                GLMakie.lines!(
                                    [loc_x[idx1], loc_x[idx2]],
                                    [loc_y[idx1], loc_y[idx2]];
                                    linewidth = 0.2,
                                    color = :black,
                                )
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(
                                        l_pos[1], l_pos[2]; text = string(connections[idx1, idx2]), fontsize = font_size
                                    )
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :red,
                                        )
                                    else
                                        GLMakie.text!(
                                            l_pos[1],
                                            l_pos[2];
                                            text = string(connections[idx1, idx2]),
                                            fontsize = font_size,
                                            color = :blue,
                                        )
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    ch_n = length(ch)
    cmap = GLMakie.resample_cmap(pal, ch_n)
    ch = setdiff(ch, sch)

    ps === :l && (sw = 2)
    ps === :m && (sw = 1)
    ps === :s && (sw = 0.0)

    for idx in 1:ch_n
        if idx in sch
            if mono
                GLMakie.scatter!(
                    loc_x[idx],
                    loc_y[idx];
                    markersize = marker_size,
                    color = :gray,
                    strokewidth = sw,
                    strokecolor = :black,
                )

            else
                GLMakie.scatter!(
                    loc_x[idx],
                    loc_y[idx];
                    markersize = marker_size,
                    color = cmap[idx],
                    colormap = pal,
                    colorrange = 1:ch_n,
                    strokewidth = sw,
                    strokecolor = :black,
                )
            end
        else
            GLMakie.scatter!(
                loc_x[idx], loc_y[idx]; markersize = marker_size, color = :gray, strokewidth = sw, strokecolor = :black
            )
        end
    end

    label_offset_x = 0.0
    label_offset_y = -0.08

    if ch_labels
        for idx in eachindex(locs[!, :label])
            if idx in ch
                GLMakie.text!(
                    loc_x[idx] + label_offset_x,
                    loc_y[idx] + label_offset_y;
                    text = locs[!, :label][idx],
                    align = (:center, :bottom),
                    fontsize = font_size,
                )
            end
        end
    end
    if sch_labels
        for idx in eachindex(locs[!, :label])
            if idx in sch
                GLMakie.text!(
                    loc_x[idx] + label_offset_x,
                    loc_y[idx] + label_offset_y;
                    text = locs[!, :label][idx],
                    align = (:center, :bottom),
                    fontsize = font_size,
                )
            end
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

    # draw weights
    if typeof(weights) <: Vector
        label_offset_x = 0.0
        label_offset_y = 0.07
        @assert length(weights) <= length(ch) "Number of weights must be ≤ number of channels to plot ($(length(ch)))."
        @assert length(weights) >= 1 "weights must contain at least one value."
        for idx in eachindex(locs[ch, :label])
            if idx in ch
                if mono
                    GLMakie.text!(
                        loc_x[idx] + label_offset_x,
                        loc_y[idx] + label_offset_y;
                        text = string(weights[idx]),
                        fontsize = font_size,
                        align = (:center, :top),
                    )
                else
                    if weights[idx] >= 0
                        GLMakie.text!(
                            loc_x[idx] + label_offset_x,
                            loc_y[idx] + label_offset_y;
                            text = string(weights[idx]),
                            fontsize = font_size,
                            color = :red,
                            align = (:center, :top),
                        )
                    else
                        GLMakie.text!(
                            loc_x[idx] + label_offset_x,
                            loc_y[idx] + label_offset_y;
                            text = string(weights[idx]),
                            fontsize = font_size,
                            color = :blue,
                            align = (:center, :top),
                        )
                    end
                end
            end
        end
    end

    loc_x_range = Tuple{Float64, Float64}[]
    loc_y_range = Tuple{Float64, Float64}[]
    for idx in eachindex(loc_x)
        push!(loc_x_range, (loc_x[idx] - 0.02, loc_x[idx] + 0.02))
        push!(loc_y_range, (loc_y[idx] - 0.02, loc_y[idx] + 0.02))
    end

    if gui
        println()

        on(events(p).mousebutton) do event
            if event.button == Mouse.left
                if event.action == Mouse.press
                    ax_x = mouseposition(ax)[1]
                    ax_y = mouseposition(ax)[2]
                    for idx in eachindex(loc_x)
                        if ax_x >= loc_x_range[idx][1] &&
                            ax_x <= loc_x_range[idx][2] &&
                            ax_y >= loc_y_range[idx][1] &&
                            ax_y <= loc_y_range[idx][2]
                            println(ch_info[idx])
                            break
                        end
                    end

                end
            end
        end

        wait(display(p))

    end

    return p

end

"""
    plot_locs(obj; <keyword arguments>)

Preview of channel locations.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `sch::Union{String, Vector{String}, Regex}`: which channels are selected
  - `ch_labels::Bool=true`: plot channel labels
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
  - `connections::Matrix{<:Real}=[0 0; 0 0]`: matrix of connections weights (channels by channels)
  - `threshold::Real=0`: threshold for plotting, see below
  - `threshold_type::Symbol=:neq`: rule for thresholding:
      + `:eq`: draw region is values are equal to threshold
      + `:neq`: draw region is values are not equal to threshold
      + `:geq`: draw region is values are ≥ to threshold
      + `:leq`: draw region is values are ≤ to threshold
      + `:g`: draw region is values are > to threshold
      + `:l`: draw region is values are < to threshold
  - `weights::Union{Bool, Vector{<:Real}}=true`: weight line widths and alpha based on connection value, if false connections values will be drawn or vector of weights
  - `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

  - `Union{GLMakie.Figure, Nothing}`
"""
function plot_locs(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    sch::Union{String, Vector{String}, Regex} = "",
    ch_labels::Bool = true,
    src_labels::Bool = false,
    det_labels::Bool = false,
    opt_labels::Bool = false,
    head::Bool = true,
    head_labels::Bool = false,
    mono::Bool = false,
    grid::Bool = false,
    ps::Symbol = :l,
    cart::Bool = false,
    plane::Symbol = :xy,
    connections::Matrix{<:Real} = [0 0; 0 0],
    threshold::Real = 0,
    threshold_type::Symbol = :neq,
    weights::Union{Bool, Vector{<:Real}} = true,
    gui::Bool = true,
)::Union{GLMakie.Figure, Nothing}

    @assert datatype(obj) != "ecog" "Use plot_locs_ecog() for ECoG data."

    ch_info = String[]
    ch = exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") : get_channel(obj; ch = ch, exclude = "")
    [push!(ch_info, channel_info(obj; ch = labels(obj)[ch[idx]], pr = false)) for idx in eachindex(ch)]
    chs = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    ch = collect(1:DataFrames.nrow(locs))

    if sch == ""
        sch = 0
    else
        sch = exclude_bads ? get_channel(obj; ch = sch, exclude = "bad") : get_channel(obj; ch = sch, exclude = "")
        sch = intersect(locs[!, :label], labels(obj)[sch])
        sch = _find_bylabel(locs, sch)
    end

    if datatype(obj) in ["eeg", "meg", "csd", "erp", "erf"]
        p = plot_locs(
            locs;
            ch = ch,
            sch = sch,
            ch_labels = ch_labels,
            head = head,
            head_labels = head_labels,
            grid = grid,
            ps = ps,
            mono = mono,
            cart = cart,
            plane = plane,
            connections = connections,
            threshold = threshold,
            threshold_type = threshold_type,
            weights = weights,
            ch_info = ch_info,
            gui = gui,
        )
    elseif datatype(obj) == "nirs"
        opt_pairs = obj.header.recording[:optode_pairs]
        src_n = length(source_labels(obj))
        det_n = length(detector_labels(obj))
        p = plot_locs_nirs(
            obj.locs,
            opt_pairs,
            src_n,
            det_n;
            src_labels = src_labels,
            det_labels = det_labels,
            opt_labels = opt_labels,
            ps = ps,
            head = head,
            head_labels = head_labels,
            cart = cart,
            grid = grid,
            mono = mono,
            plane = plane,
            ch_info = ch_info,
        )
    elseif datatype(obj) == "ecog"
        _warn("ECOG locs are not supported yet.")
        return nothing
    elseif datatype(obj) == "seeg"
        _warn("SEEG locs are not supported yet.")
        return nothing
    elseif datatype(obj) == "ieeg"
        _warn("iEEG locs are not supported yet.")
        return nothing
    elseif datatype(obj) in ["sensors", "eda", "mep", "tpt"]
        _warn("For $(datatype(obj)) object type locs are not available.")
    end

    return p

end
