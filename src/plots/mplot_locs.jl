export mplot_locs
export mplot_gridlocs

"""
    mplot_locs(locs; <keyword arguments>)

Preview channel locations.

# Arguments

- `locs::DataFrame`: columns: `channel`, `labels`, `loc_radius`, `loc_theta`, `loc_x`, `loc_y`, `loc_z`, `loc_radius_sph`, `loc_theta_sph`, `loc_phi_sph`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=1:DataFrames.nrow(locs)`: list of locations to plot, default is all locations
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: which channels should be highlighted
- `ch_labels::Bool=true`: plot locations labels
- `head::Bool=true`: draw head
- `head_labels::Bool=false`: plot head labels
- `mono::Bool=false`: use color or gray palette
- `grid::Bool=false`: draw grid, useful for locating positions
- `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `plane::Symbol=:xy`: which plane to plot:
    - `:xy`: horizontal (top)
    - `:xz`: coronary (front)
    - `:yz`: sagittal (side)
- `transparent::Bool=false`: if true, do not paint the background
- `connections::Matrix{<:Real}=[0 0; 0 0]`: matrix of connections weights (channels by channels)
- `threshold::Real=0`: threshold for plotting, see below
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
- `weights::Union{Bool, Vector{<:Real}}=true`: weight line widths and alpha based on connection value, if false connections values will be drawn or vector of weights

# Returns

- `p::GLMakie.Figure`
"""
function mplot_locs(locs::DataFrame; ch::Union{Int64, Vector{Int64}, AbstractRange}=1:DataFrames.nrow(locs), selected::Union{Int64, Vector{Int64}, AbstractRange}=0, ch_labels::Bool=true, head::Bool=true, head_labels::Bool=false, mono::Bool=false, grid::Bool=false, large::Bool=true, cart::Bool=false, plane::Symbol=:xy, transparent::Bool=false, connections::Matrix{<:Real}=[0 0; 0 0], threshold::Real=0, threshold_type::Symbol=:neq, weights::Union{Bool, Vector{<:Real}}=true)::GLMakie.Figure

    _check_var(plane, [:xy, :yz, :xz], "plane")
    pal = mono ? :grays : :darktest
    sch_labels = ch_labels

    if plane === :xy
        head_shape = large ? rotr90(FileIO.load(joinpath(res_path, "head_t_large.png"))) : rotr90(FileIO.load(joinpath(res_path, "head_t_small.png")))
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
        head_shape = large ? rotr90(FileIO.load(joinpath(res_path, "head_f_large.png"))) : rotr90(FileIO.load(joinpath(res_path, "head_f_small.png")))
        if !cart
            loc_x = zeros(length(ch))
            loc_y = zeros(length(ch))
            for idx in 1:length(ch)
                loc_x[idx], _, loc_y[idx] = sph2cart(locs[ch, :loc_radius_sph][idx], locs[ch, :loc_theta_sph][idx], locs[ch, :loc_phi_sph][idx])
            end
        else
            loc_x = locs[ch, :loc_x]
            loc_y = locs[ch, :loc_z]
        end
    elseif plane === :yz
        head_shape = large ? rotr90(FileIO.load(joinpath(res_path, "head_s_large.png"))) : rotr90(FileIO.load(joinpath(res_path, "head_s_small.png")))
        if !cart
            loc_x = zeros(length(ch))
            loc_y = zeros(length(ch))
            for idx in 1:length(ch)
                _, loc_x[idx], loc_y[idx] = sph2cart(locs[ch, :loc_radius_sph][idx], locs[ch, :loc_theta_sph][idx], locs[ch, :loc_phi_sph][idx])
            end
        else
            loc_x = locs[ch, :loc_y]
            loc_y = locs[ch, :loc_z]
        end
    end

    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)
    loc_y .= -loc_y

    head12 = false
    maximum(abs.(locs[:, :loc_x])) <= 1.2 && maximum(abs.(locs[:, :loc_y])) <= 1.2 && maximum(abs.(locs[:, :loc_z])) <= 1.5 && (head12 = true)

    if head12
        xt = round.(linspace(0, size(head_shape, 1), 25))
        yt = round.(linspace(0, size(head_shape, 2), 25))
        xl = (0, size(head_shape, 1))
        yl = (0, size(head_shape, 2))
        origin = size(head_shape) .÷ 2
        if large
            marker_size = length(ch) > 64 ? 10 : 20
            font_size = 12
            loc_x = @. round(origin[1] + (loc_x * 250), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 250), digits=2)
            length(ch) > 64 && (ch_labels = false)
        else
            marker_size = length(ch) > 64 ? 5 : 10
            font_size = 6
            ch_labels = false
            sch_labels = false
            grid = false
            loc_x = @. round(origin[1] + (loc_x * 100), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 100), digits=2)
        end
    else
        m = zeros(RGBA{FixedPointNumbers.N0f8}, size(head_shape) .+ 200)
        m[101:100+size(head_shape, 1), 101:100+size(head_shape, 2)] .= head_shape
        head_shape = m
        xt = (round.(linspace(0, size(head_shape, 1), 25)), string.(round.(linspace(-1.6, 1.6, 25), digits=1)))
        yt = (round.(linspace(0, size(head_shape, 2), 25)), string.(round.(linspace(1.6, -1.6, 25), digits=1)))
        xl = (0, size(head_shape, 1))
        yl = (0, size(head_shape, 2))
        origin = size(head_shape) ./ 2 .+ 1
        if large
            marker_size = length(ch) > 64 ? 20 : 10
            font_size = 12
            loc_x = @. round(origin[1] + (loc_x * 250), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 250), digits=2)
            length(ch) > 64 && (ch_labels = false)
        else
            marker_size = length(ch) > 64 ? 10 : 5
            font_size = 6
            ch_labels = false
            sch_labels = false
            grid = false
            loc_x = @. round(origin[1] + (loc_x * 100), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 100), digits=2)
        end
    end

    ma = ch_labels ? 0.5 : 1.0

    # prepare plot
    plot_size = size(head_shape)
    p = GLMakie.Figure(size=plot_size)
    if grid
        ax = GLMakie.Axis(p[1, 1],
                          aspect=1,
                          xlabel="",
                          ylabel="",
                          title="",
                          xticks=xt,
                          xminorticksvisible=true,
                          xminorticks=IntervalsBetween(2),
                          yticks=yt,
                          yminorticksvisible=true,
                          yminorticks=IntervalsBetween(2),
                          xautolimitmargin=(0, 0),
                          yautolimitmargin=(0, 0),
                          backgroundcolor=transparent ? :transparent : :white)
    else
        # prepare plot
        ax = GLMakie.Axis(p[1, 1],
                          aspect=1,
                          xlabel="",
                          ylabel="",
                          title="",
                          xautolimitmargin=(0, 0),
                          yautolimitmargin=(0, 0),
                          backgroundcolor=transparent ? :transparent : :white)
        hidedecorations!(ax, grid=true)
        hidespines!(ax)
    end
    GLMakie.xlims!(ax, xl)
    GLMakie.ylims!(ax, yl)

    head && GLMakie.image!(head_shape)

    # draw connections
    if connections != [0 0; 0 0]
        selected = ""
        @assert size(connections, 1) == length(ch) "Length of channel and number of connections rows must be equal."
        _check_var(threshold_type, [:eq, :neq, :geq, :leq, :g, :l], "threshold_type")
        m_tmp = normalize_n(abs.(connections))

        for idx1 in axes(connections, 1)
            for idx2 in 2:size(connections, 1)
                if idx1 != idx2
                    if threshold_type === :g
                        if connections[idx1, idx2] > threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                       [loc_y[idx1], loc_y[idx2]],
                                                        linewidth=6 * m_tmp[idx1, idx2],
                                                        alpha=0.25 * m_tmp[idx1, idx2],
                                                        color=:black)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                       [loc_y[idx1], loc_y[idx2]],
                                                       linewidth=6 * m_tmp[idx1, idx2],
                                                       alpha=0.25 * m_tmp[idx1, idx2],
                                                       color=:red)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                       [loc_y[idx1], loc_y[idx2]],
                                                       linewidth=6 * m_tmp[idx1, idx2],
                                                       alpha=0.25 * m_tmp[idx1, idx2],
                                                       color=:black,
                                                       linestyle=:dot)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                       [loc_y[idx1], loc_y[idx2]],
                                                       linewidth=6 * m_tmp[idx1, idx2],
                                                       alpha=0.25 * m_tmp[idx1, idx2],
                                                       color=:blue)
                                    end
                                end
                            else
                                GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                               [loc_y[idx1], loc_y[idx2]],
                                               linewidth=0.2,
                                               color=:black)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(l_pos[1],
                                                  l_pos[2],
                                                  text=string(connections[idx1, idx2]),
                                                  fontsize=font_size+2)
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(l_pos[1],
                                                      l_pos[2],
                                                      text=string(connections[idx1, idx2]),
                                                      fontsize=font_size+2,
                                                      color=:red)
                                    else
                                        GLMakie.text!(l_pos[1],
                                                      l_pos[2],
                                                      text=string(connections[idx1, idx2]),
                                                      fontsize=font_size+2,
                                                      color=:blue)
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :l
                        if connections[idx1, idx2] < threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:red)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black,
                                                    linestyle=:dot)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:blue)
                                    end
                                end
                            else
                                GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]],
                                            linewidth=0.2,
                                            color=:black)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(l_pos[1], l_pos[2],
                                                             text=string(connections[idx1, idx2]),
                                                             fontsize=font_size+2)
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(l_pos[1], l_pos[2],
                                                                 text=string(connections[idx1, idx2]),
                                                                 fontsize=font_size+2,
                                                                 color=:red)
                                    else
                                        GLMakie.text!(l_pos[1], l_pos[2],
                                                                 text=string(connections[idx1, idx2]),
                                                                 fontsize=font_size+2,
                                                                 color=:blue)
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :eq
                        if connections[idx1, idx2] == threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:red)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black,
                                                    linestyle=:dot)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:blue)
                                    end
                                end
                            else
                                GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]],
                                            linewidth=0.2,
                                            color=:black)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(l_pos[1],
                                                             l_pos[2],
                                                             text=string(connections[idx1, idx2]),
                                                             fontsize=font_size+2)
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(l_pos[1],
                                                                 l_pos[2],
                                                                 text=string(connections[idx1, idx2]),
                                                                 fontsize=font_size+2,
                                                                 color=:red)
                                    else
                                        GLMakie.text!(l_pos[1],
                                                                 l_pos[2],
                                                                 text=string(connections[idx1, idx2]),
                                                                 fontsize=font_size+2,
                                                                 color=:blue)
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :neq
                        if connections[idx1, idx2] != threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:red)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black,
                                                    linestyle=:dot)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:blue)
                                    end
                                end
                            else
                                GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]],
                                            linewidth=0.2,
                                            color=:black)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(l_pos[1],
                                                             l_pos[2],
                                                             text=string(connections[idx1, idx2]),
                                                             fontsize=font_size+2)
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(l_pos[1],
                                                                 l_pos[2],
                                                                 text=string(connections[idx1, idx2]),
                                                                 fontsize=font_size+2,
                                                                 color=:red)
                                    else
                                        GLMakie.text!(l_pos[1],
                                                                 l_pos[2],
                                                                 text=string(connections[idx1, idx2]),
                                                                 fontsize=font_size+2,
                                                                 color=:blue)
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :leq
                        if connections[idx1, idx2] <= threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:red)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black,
                                                    linestyle=:dot)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:blue)
                                    end
                                end
                            else
                                GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]],
                                            linewidth=0.2, color=:black)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(l_pos[1],
                                                             l_pos[2],
                                                             text=string(connections[idx1, idx2]),
                                                             fontsize=font_size+2)
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(l_pos[1],
                                                                 l_pos[2],
                                                                 text=string(connections[idx1, idx2]),
                                                                 fontsize=font_size+2,
                                                                 color=:red)
                                    else
                                        GLMakie.text!(l_pos[1],
                                                                 l_pos[2],
                                                                 text=string(connections[idx1, idx2]),
                                                                 fontsize=font_size+2,
                                                                 color=:blue)
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :geq
                        if connections[idx1, idx2] >= threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:red)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:black,
                                                    linestyle=:dot)
                                    else
                                        GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                                    [loc_y[idx1], loc_y[idx2]],
                                                    linewidth=6 * m_tmp[idx1, idx2],
                                                    alpha=0.25 * m_tmp[idx1, idx2],
                                                    color=:blue)
                                    end
                                end
                            else
                                GLMakie.lines!([loc_x[idx1], loc_x[idx2]],
                                            [loc_y[idx1], loc_y[idx2]],
                                            linewidth=0.2,
                                            color=:black)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    GLMakie.text!(l_pos[1],
                                                  l_pos[2],
                                                  text=string(connections[idx1, idx2]),
                                                  fontsize=font_size+2)
                                else
                                    if connections[idx1, idx2] >= 0
                                        GLMakie.text!(l_pos[1],
                                                      l_pos[2],
                                                      text=string(connections[idx1, idx2]),
                                                      fontsize=font_size+2,
                                                      color=:red)
                                    else
                                        GLMakie.text!(l_pos[1],
                                                      l_pos[2],
                                                      text=string(connections[idx1, idx2]),
                                                      fontsize=font_size+2,
                                                      color=:blue)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    cmap = GLMakie.resample_cmap(pal, length(ch))
    ch = setdiff(ch, selected)

    for idx in eachindex(locs[!, :label])
        if idx in ch
            GLMakie.scatter!((loc_x[idx], loc_y[idx]),
                             color=:darkgrey,
                             markersize=marker_size,
                             alpha=ma)
        end
    end

    for idx in eachindex(locs[!, :label])
        if idx in selected
            if mono != true
                GLMakie.scatter!((loc_x[idx], loc_y[idx]),
                                 color=cmap[idx],
                                 colormap=pal,
                                 colorrange=length(selected),
                                 markersize=marker_size,
                                 alpha=ma)
            else
                GLMakie.scatter!((loc_x[idx], loc_y[idx]),
                                 color=:lightgrey,
                                 colormap=pal,
                                 colorrange=length(selected),
                                 markersize=marker_size)
            end
        end
    end

    label_offset_x = 10
    label_offset_y = -10

    if ch_labels
        for idx in eachindex(locs[!, :label])
            if idx in ch
                GLMakie.text!(loc_x[idx] + label_offset_x,
                              loc_y[idx] + label_offset_y,
                              text=locs[!, :label][idx],
                              fontsize=font_size)
            end
        end
    end
    if sch_labels
        for idx in eachindex(locs[!, :label])
            if idx in selected
                GLMakie.text!(loc_x[idx] + label_offset_x,
                              loc_y[idx] + label_offset_y,
                              text=locs[!, :label][idx],
                              fontsize=font_size)
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
            if large
                fid_loc_x = @. origin[1] + (fid_loc_x * 260)
                fid_loc_y = @. origin[2] + (fid_loc_y * 260)
            else
                fid_loc_x = @. origin[1] + (fid_loc_x * 100)
                fid_loc_y = @. origin[2] + (fid_loc_y * 100)
            end
            GLMakie.text!(fid_loc_x,
                          fid_loc_y,
                          text=fid_names[idx],
                          fontsize=font_size + 2,
                          align = (:center, :center))
        end
    end

    # draw weights
    if typeof(weights) <: Vector
        if ch_labels
            label_offset_x = -10
            label_offset_y = 10
        else
            label_offset_x = 10
            label_offset_y = -10
        end
        @assert length(weights) <= length(ch) "Number of weights must be ≤ number of channels to plot ($(length(ch)))."
        @assert length(weights) >= 1 "weights must contain at least one value."
        for idx in eachindex(locs[ch, :label])
            if idx in ch
                if mono
                    GLMakie.text!(loc_x[idx] + label_offset_x,
                                  loc_y[idx] + label_offset_y,
                                  text=string(weights[idx]),
                                  fontsize=font_size)
                else
                    if weights[idx] >= 0
                        GLMakie.text!(loc_x[idx] + label_offset_x,
                                      loc_y[idx] + label_offset_y,
                                      text=string(weights[idx]),
                                      fontsize=font_size,
                                      color=:red)
                    else
                        GLMakie.text!(loc_x[idx] + label_offset_x,
                                      loc_y[idx] + label_offset_y,
                                      text=string(weights[idx]),
                                      fontsize=font_size,
                                      color=:blue)
                    end
                end
            end
        end
    end

    return p

end

"""
    mplot_locs(obj; <keyword arguments>)

Preview of channel locations.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `selected::Union{String, Vector{String}, Regex}`: which channels should be highlighted
- `ch_labels::Bool=true`: plot channel labels
- `src_labels::Bool=false`: plot source labels
- `det_labels::Bool=false`: plot detector labels
- `opt_labels::Bool=false`: plot optode type (S for source, D for detector) and number
- `head::Bool=true`: draw head
- `head_labels::Bool=false`: plot head labels
- `mono::Bool=false`: use color or gray palette
- `grid::Bool=false`: draw grid, useful for locating positions
- `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `plane::Symbol=:xy`: which plane to plot:
    - `:xy`: horizontal (top)
    - `:xz`: coronary (front)
    - `:yz`: sagittal (side)
- `transparent::Bool=false`: if true, do not paint the background
- `connections::Matrix{<:Real}=[0 0; 0 0]`: matrix of connections weights (channels by channels)
- `threshold::Real=0`: threshold for plotting, see below
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
- `weights::Union{Bool, Vector{<:Real}}=true`: weight line widths and alpha based on connection value, if false connections values will be drawn or vector of weights

# Returns

- `Union{GLMakie.Figure, Nothing}`
"""
function mplot_locs(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, selected::Union{String, Vector{String}, Regex}="", ch_labels::Bool=true, src_labels::Bool=false, det_labels::Bool=false, opt_labels::Bool=false, head::Bool=true, head_labels::Bool=false, mono::Bool=false, grid::Bool=false, large::Bool=true, cart::Bool=false, plane::Symbol=:xy, transparent::Bool=false, connections::Matrix{<:Real}=[0 0; 0 0], threshold::Real=0, threshold_type::Symbol=:neq, weights::Union{Bool, Vector{<:Real}}=true)::Union{GLMakie.Figure, Nothing}

    @assert datatype(obj) != "ecog" "Use mplot_locs_ecog() for ECoG data."

    ch = get_channel(obj, ch=ch)
    chs = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    ch = collect(1:DataFrames.nrow(locs))

    if selected == ""
        selected = 0
    else
        selected = get_channel(obj, ch=selected)
        selected = intersect(locs[!, :label], labels(obj)[selected])
        selected = _find_bylabel(locs, selected)
    end

    if datatype(obj) == "nirs"
        opt_pairs = obj.header.recording[:optode_pairs]
        src_n = length(source_labels(obj))
        det_n = length(detector_labels(obj))
        p = mplot_locs_nirs(obj.locs,
                           opt_pairs,
                           src_n,
                           det_n,
                           src_labels=src_labels,
                           det_labels=det_labels,
                           opt_labels=opt_labels,
                           head=head,
                           head_labels=head_labels,
                           grid=grid,
                           mono=mono)
        return p
    else
        p = mplot_locs(locs,
                      ch=ch,
                      selected=selected,
                      ch_labels=ch_labels,
                      head=head,
                      head_labels=head_labels,
                      grid=grid,
                      large=large,
                      mono=mono,
                      cart=cart,
                      plane=plane,
                      transparent=transparent,
                      connections=connections,
                      threshold=threshold,
                      threshold_type=threshold_type,
                      weights=weights)
    end

    return p

end

"""
    mplot_gridlocs()

Plot a simplified plot of 10-20 EEG channels on a grid.

# Arguments

- `mono::Bool=false`: use color or gray palette

# Returns

- `p::GLMakie.Figure`
"""
function mplot_gridlocs(; mono::Bool=false)::GLMakie.Figure

    pal = mono ? :grays : :darktest

    plot_size=(800, 800)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      aspect=1,
                      xlabel="",
                      ylabel="",
                      title="",
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0))
    hidedecorations!(ax, grid=true)
    hidespines!(ax)
    GLMakie.xlims!(ax, (-1.2, 1.2))
    GLMakie.ylims!(ax, (-1.2, 1.2))

    GLMakie.lines!([-1, 1], [-1, -1], color=:black, linewidth=0.2)
    GLMakie.lines!([-1, 1], [1, 1], color=:black, linewidth=0.2)

    GLMakie.lines!([-1, -1], [-1, 1], color=:black, linewidth=0.2)
    GLMakie.lines!([1, 1], [-1, 1], color=:black, linewidth=0.2)

    GLMakie.lines!([-1, -0.5], [0.5, 1], color=:black, linewidth=0.5)
    GLMakie.lines!([0.5, 1], [1, 0.5], color=:black, linewidth=0.5)
    GLMakie.lines!([-1, -0.5], [-0.5, -1], color=:black, linewidth=0.5)
    GLMakie.lines!([0.5, 1], [-1, -0.5], color=:black, linewidth=0.5)

    GLMakie.lines!([-0.5, 0.5], [-1, -1], color=:black, linewidth=0.5)
    GLMakie.lines!([-1, 1], [-0.5, -0.5], color=:black, linewidth=0.5)
    GLMakie.lines!([-1, 1], [0, 0], color=:black, linewidth=0.5)
    GLMakie.lines!([-1, 1], [0.5, 0.5], color=:black, linewidth=0.5)
    GLMakie.lines!([-0.5, 0.5], [1, 1], color=:black, linewidth=0.5)

    GLMakie.lines!([-1, -1], [-0.5, 0.5], color=:black, linewidth=0.5)
    GLMakie.lines!([-0.5, -0.5], [-1, 1], color=:black, linewidth=0.5)
    GLMakie.lines!([0, 0], [-1, 1], color=:black, linewidth=0.5)
    GLMakie.lines!([0.5, 0.5], [-1, 1], color=:black, linewidth=0.5)
    GLMakie.lines!([1, 1], [-0.5, 0.5], color=:black, linewidth=0.5)

    loc_x = [-0.5, 0, 0.5, -1, -0.5, 0, 0.5, 1, -1, -0.5, 0, 0.5, 1, -1, -0.5, 0, 0.5, 1, -0.5, 0, 0.5]
    loc_y = [1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, -0.5, -0.5, -0.5, -0.5, -0.5, -1, -1, -1]
    loc_lab = ["Fp1", "Fpz", "Fp2", "F7", "F3", "Fz", "F4", "F8", "T3", "C3", "Cz", "C4", "T4", "T5", "P3", "Pz", "P4", "T6", "O1", "Oz", "O2"]
    font_size = 16
    label_offset_x = 0.01
    label_offset_y = 0.01
    for idx in eachindex(loc_x)
        GLMakie.scatter!(loc_x[idx],
                         loc_y[idx],
                         markersize=3.0,
                         color=:black)
        GLMakie.text!(loc_x[idx] + label_offset_x,
                      loc_y[idx] + label_offset_y,
                      text=loc_lab[idx],
                      fontsize=font_size)
    end

    return p

end
