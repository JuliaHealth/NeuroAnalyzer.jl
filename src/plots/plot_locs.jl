export plot_locs
export plot_locs3d

"""
    plot_locs(locs; <keyword arguments>)

Preview channel locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=1:nrow(locs)`: list of locations to plot, default is all locations
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

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs(locs::DataFrame; ch::Union{Int64, Vector{Int64}, AbstractRange}=1:nrow(locs), selected::Union{Int64, Vector{Int64}, AbstractRange}=0, ch_labels::Bool=true, head::Bool=true, head_labels::Bool=false, mono::Bool=false, grid::Bool=false, large::Bool=true, cart::Bool=false, plane::Symbol=:xy, transparent::Bool=false, connections::Matrix{<:Real}=[0 0; 0 0], threshold::Real=0, threshold_type::Symbol=:neq, weights::Union{Bool, Vector{<:Real}}=true)::Plots.Plot{Plots.GRBackend}

    _check_var(plane, [:xy, :yz, :xz], "plane")
    pal = mono ? :grays : :darktest
    sch_labels = ch_labels

    if plane === :xy
        head_shape = large ? FileIO.load(joinpath(res_path, "head_t_large.png")) : FileIO.load(joinpath(res_path, "head_t_small.png"))
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
        head_shape = large ? FileIO.load(joinpath(res_path, "head_f_large.png")) : FileIO.load(joinpath(res_path, "head_f_small.png"))
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
        head_shape = large ? FileIO.load(joinpath(res_path, "head_s_large.png")) : FileIO.load(joinpath(res_path, "head_s_small.png"))
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

    head12 = false
    maximum(abs.(locs[:, :loc_x])) <= 1.2 && maximum(abs.(locs[:, :loc_y])) <= 1.2 && maximum(abs.(locs[:, :loc_z])) <= 1.5 && (head12 = true)

    if head12
        xt = (round.(linspace(0, size(head_shape, 1), 25)), string.(round.(linspace(-1.2, 1.2, 25), digits=1)))
        yt = (round.(linspace(0, size(head_shape, 2), 25)), string.(round.(linspace(1.2, -1.2, 25), digits=1)))
        xl = (0, size(head_shape, 1))
        yl = (0, size(head_shape, 2))
        !head && (loc_y .= -loc_y)
        origin = size(head_shape) ./ 2 .+ 1
        if large
            marker_size = length(ch) > 64 ? 3 : 6
            font_size = 6
            loc_x = @. round(origin[1] + (loc_x * 250), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 250), digits=2)
            length(ch) > 64 && (ch_labels = false)
        else
            marker_size = length(ch) > 64 ? 2 : 4
            font_size = 4
            ch_labels = false
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
        !head && (loc_y .= -loc_y)
        origin = size(head_shape) ./ 2 .+ 1
        if large
            marker_size = length(ch) > 64 ? 3 : 6
            font_size = 6
            loc_x = @. round(origin[1] + (loc_x * 250), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 250), digits=2)
            length(ch) > 64 && (ch_labels = false)
        else
            marker_size = length(ch) > 64 ? 2 : 4
            font_size = 4
            ch_labels = false
            sch_labels = false
            grid = false
            loc_x = @. round(origin[1] + (loc_x * 100), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 100), digits=2)
        end
    end

    ma = ch_labels ? 0.5 : 1.0

    if grid
        p = Plots.plot(grid=true,
                       framestyle=:grid,
                       palette=pal,
                       aspect_ratio=1,
                       size=size(head_shape) .+ 50,
                       right_margin=10*Plots.px,
                       bottom_margin=0*Plots.px,
                       top_margin=10*Plots.px,
                       left_margin=0*Plots.px,
                       xtickfontsize=font_size,
                       ytickfontsize=font_size,
                       xticks=xt,
                       yticks=yt,
                       xlims=xl,
                       ylims=yl,
                       background_color=transparent ? :transparent : :white,
                       foreground_color=:black)
    else
        if large
            p = Plots.plot(grid=false,
                           framestyle=:none,
                           border=:none,
                           palette=pal,
                           aspect_ratio=1,
                           size=head12 ? size(head_shape) : size(head_shape) .+ 50,
                           right_margin=head12 ? -100*Plots.px : 10*Plots.px,
                           bottom_margin=head12 ? -100*Plots.px : 0*Plots.px,
                           top_margin=head12 ? -100*Plots.px : 10*Plots.px,
                           left_margin=head12 ? -100*Plots.px : 0*Plots.px,
                           ticks_fontsize=font_size,
                           xticks=xt,
                           yticks=yt,
                           xlims=xl,
                           ylims=yl,
                           background_color=transparent ? :transparent : :white,
                           foreground_color=:black)
        else
            p = Plots.plot(grid=false,
                           framestyle=:none,
                           border=:none,
                           palette=pal,
                           aspect_ratio=1,
                           size=head12 ? size(head_shape) : size(head_shape) .+ 7,
                           right_margin=head12 ? -10*Plots.px : -10*Plots.px,
                           bottom_margin=head12 ? -30*Plots.px : -30*Plots.px,
                           top_margin=head12 ? -20*Plots.px : -20*Plots.px,
                           left_margin=head12 ? -40*Plots.px : -30*Plots.px,
                           xticks=xt,
                           yticks=yt,
                           xlims=xl,
                           ylims=yl,
                           background_color=transparent ? :transparent : :white,
                           foreground_color=:black)
        end
    end

    head && (p = Plots.plot!(head_shape))

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
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:red, legend=false)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, ls=:dot, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:blue, legend=false)
                                    end
                                end
                            else
                                p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2)))
                                else
                                    if connections[idx1, idx2] >= 0
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :red)))
                                    else
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :blue)))
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :l
                        if connections[idx1, idx2] < threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:red, legend=false)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, ls=:dot, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:blue, legend=false)
                                    end
                                end
                            else
                                p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2)))
                                else
                                    if connections[idx1, idx2] >= 0
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :red)))
                                    else
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :blue)))
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :eq
                        if connections[idx1, idx2] == threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:red, legend=false)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, ls=:dot, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:blue, legend=false)
                                    end
                                end
                            else
                                p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2)))
                                else
                                    if connections[idx1, idx2] >= 0
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :red)))
                                    else
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :blue)))
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :neq
                        if connections[idx1, idx2] != threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:red, legend=false)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, ls=:dot, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:blue, legend=false)
                                    end
                                end
                            else
                                p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2)))
                                else
                                    if connections[idx1, idx2] >= 0
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :red)))
                                    else
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :blue)))
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :leq
                        if connections[idx1, idx2] <= threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:red, legend=false)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, ls=:dot, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:blue, legend=false)
                                    end
                                end
                            else
                                p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2)))
                                else
                                    if connections[idx1, idx2] >= 0
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :red)))
                                    else
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :blue)))
                                    end
                                end
                            end
                        end
                    elseif threshold_type === :geq
                        if connections[idx1, idx2] >= threshold
                            if weights
                                if connections[idx1, idx2] > 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:red, legend=false)
                                    end
                                elseif connections[idx1, idx2] < 0
                                    if mono
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, ls=:dot, legend=false)
                                    else
                                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:blue, legend=false)
                                    end
                                end
                            else
                                p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                                l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
                                if mono
                                    p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2)))
                                else
                                    if connections[idx1, idx2] >= 0
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :red)))
                                    else
                                        p = Plots.plot!(annotations=(l_pos[1], l_pos[2], Plots.text(connections[idx1, idx2], pointsize=font_size+2, :blue)))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    ch = setdiff(ch, selected)

    for idx in eachindex(locs[!, :label])
        if idx in ch
            p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor=Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
        end
    end

    for idx in eachindex(locs[!, :label])
        if idx in selected
            if mono != true
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=idx,
                                markerstrokecolor=Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markeralpha=ma,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            else
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor=Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            end
        end
    end

    label_offset_y = -10

    if ch_labels
        for idx in eachindex(locs[!, :label])
            if idx in ch
                Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + label_offset_y, Plots.text(locs[!, :label][idx], pointsize=font_size)))
            end
#=            if idx in selected
                Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + 1, Plots.text(locs[!, :label][idx], pointsize=font_size)))
            end
=#
        end
    end
    if sch_labels
        for idx in eachindex(locs[!, :label])
            if idx in selected
                Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + label_offset_y, Plots.text(locs[!, :label][idx], pointsize=font_size)))
            end
        end
    end

    if head_labels
        fid_names = ["NAS", "IN", "LPA", "RPA"]
        for idx in eachindex(NeuroAnalyzer.fiducial_points)
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
                fid_loc_x = @. origin[1] + (fid_loc_x * 250)
                fid_loc_y = @. origin[2] - (fid_loc_y * 250)
            else
                fid_loc_x = @. origin[1] - (fid_loc_x * 100)
                fid_loc_y = @. origin[2] - (fid_loc_y * 100)
            end
            p = Plots.plot!(annotations=(fid_loc_x, fid_loc_y, Plots.text(fid_names[idx], pointsize=font_size + 2)))
        end
    end

    # draw weights
    if typeof(weights) <: Vector
        @assert length(weights) <= length(ch) "Number of weights must be ≤ number of channels to plot ($(length(ch)))."
        @assert length(weights) >= 1 "weights must contain at least one value."
        for idx in eachindex(locs[ch, :label])
            if idx in ch
                if mono
                    Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + round(p.attr[:size][2] * 0.03), Plots.text(string(weights[idx]), pointsize=font_size)))
                else
                    if weights[idx] >= 0
                        Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + round(p.attr[:size][2] * 0.03), Plots.text(string(weights[idx]), :red, pointsize=font_size)))
                    else
                        Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + round(p.attr[:size][2] * 0.03), Plots.text(string(weights[idx]), :blue, pointsize=font_size)))
                    end
                end
            end
        end
    end

    Plots.plot!(p)

    return p

end

"""
    plot_locs3d(locs; <keyword arguments>)

3D preview of channel locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Union{Int64, Vector{Int64}}=1:nrow(locs)`: list of channels, default is all channels
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: which channel should be highlighted
- `ch_labels::Bool=true`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or gray palette
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use spherical coordinates
- `camera::Tuple{Real, Real}=(20, 45)`: camera position -- (XY plane angle, XZ plane angle)

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs3d(locs::DataFrame; ch::Union{Int64, Vector{Int64}, AbstractRange}=1:nrow(locs), selected::Union{Int64, Vector{Int64}, AbstractRange}=0, ch_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, cart::Bool=false, camera::Tuple{Real, Real}=(20, 45))::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :darktest

    if !cart
        loc_x = zeros(nrow(locs))
        loc_y = zeros(nrow(locs))
        loc_z = zeros(nrow(locs))
        for idx in 1:nrow(locs)
            loc_x[idx], loc_y[idx], loc_z[idx] = sph2cart(locs[idx, :loc_radius_sph], locs[idx, :loc_theta_sph], locs[idx, :loc_phi_sph])
        end
    else
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
        loc_z = locs[!, :loc_z]
    end

    if maximum(locs[:, :loc_x]) <= 1.2 && maximum(locs[:, :loc_y]) <= 1.2 && maximum(locs[:, :loc_z]) <= 1.5
        x_lim = (-1.5, 1.5)
        y_lim = (-1.5, 1.5)
        z_lim = (-1.5, 1.5)
        plot_size = 640
    else
        x_lim = (-2.0, 2.0)
        y_lim = (-2.0, 2.0)
        z_lim = (-2.0, 2.0)
        plot_size = 850
    end

    marker_size = 6
    font_size = 6

    ch = setdiff(ch, selected)

    p = Plots.scatter3d(grid=true,
                        palette=pal,
                        size=(plot_size, plot_size),
                        aspect_ratios=:equal,
                        right_margin=-20*Plots.px,
                        bottom_margin=-20*Plots.px,
                        top_margin=-20*Plots.px,
                        left_margin=-20*Plots.px,
                        legend=false,
                        camera=camera,
                        xticks=([-1, 0, 1]),
                        yticks=([-1, 0, 1]),
                        zticks=([-1, 0, 1]),
                        xlabel="X",
                        ylabel="Y",
                        zlabel="Z",
                        xlim=x_lim,
                        ylim=y_lim,
                        zlim=z_lim)

    p = Plots.scatter3d!((loc_x, loc_y, loc_z),
                         markercolor=:gray,
                         markerstrokecolor=Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                         markershape=:circle,
                         markersize=marker_size,
                         markerstrokewidth=0,
                         markerstrokealpha=0)

    if selected != 0
        if mono
            p = Plots.scatter3d!((loc_x[selected], loc_y[selected], loc_z[selected]),
                                 markercolor=:gray,
                                 markerstrokecolor=Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                 markershape=:circle,
                                 markersize=marker_size,
                                 markerstrokewidth=0,
                                 markerstrokealpha=0)
        else
            for idx in selected
                p = Plots.scatter3d!((loc_x[idx], loc_y[idx], loc_z[idx]),
                                     markercolor=idx,
                                     markerstrokecolor=Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                     markershape=:circle,
                                     markersize=marker_size,
                                     markerstrokewidth=0,
                                     markerstrokealpha=0)
            end
        end
    end

    if ch_labels
        for idx in eachindex(locs[!, :label])
            idx in ch && Plots.annotate!(loc_x[idx] * 1.1, loc_y[idx] * 1.1, loc_z[idx] * 1.1, Plots.text(locs[!, :label][idx], font_size))
            idx in selected && Plots.annotate!(loc_x[idx] * 1.1, loc_y[idx] * 1.1, loc_z[idx] * 1.1, Plots.text(locs[!, :label][idx], font_size))
        end
    end

    if head_labels
        fid_names = ["NAS", "IN", "LPA", "RPA"]
        for idx in 1:length(NeuroAnalyzer.fiducial_points)
            Plots.annotate!(NeuroAnalyzer.fiducial_points[idx][1],
                            NeuroAnalyzer.fiducial_points[idx][2],
                            NeuroAnalyzer.fiducial_points[idx][3],
                            Plots.text(fid_names[idx], font_size))
        end
    end

    p = Plots.plot!(p)

    return p

end

"""
    plot_locs(obj; <keyword arguments>)

Preview of channel locations.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `selected::Union{String, Vector{String}}`: which channels should be highlighted
- `ch_labels::Bool=true`: plot channel labels
- `src_labels::Bool=false`: plot source labels
- `det_labels::Bool=false`: plot detector labels
- `opt_labels::Bool=false`: plot optode type (S for source, D for detector) and number
- `head::Bool=true`: draw head
- `head_labels::Bool=false`: plot head labels
- `d::Int64=2`: 2- or 3-dimensional plot
- `mono::Bool=false`: use color or gray palette
- `grid::Bool=false`: draw grid, useful for locating positions
- `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `plane::Symbol=:xy`: which plane to plot:
    - `:xy`: horizontal (top)
    - `:xz`: coronary (front)
    - `:yz`: sagittal (side)
- `interactive::Bool=true`: if true, use interactive 3-dimensional plot
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

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, selected::Union{String, Vector{String}}="", ch_labels::Bool=true, src_labels::Bool=false, det_labels::Bool=false, opt_labels::Bool=false, head::Bool=true, head_labels::Bool=false, d::Int64=2, mono::Bool=false, grid::Bool=false, large::Bool=true, cart::Bool=false, plane::Symbol=:xy, interactive::Bool=true, transparent::Bool=false, connections::Matrix{<:Real}=[0 0; 0 0], threshold::Real=0, threshold_type::Symbol=:neq, weights::Union{Bool, Vector{<:Real}}=true)::Plots.Plot{Plots.GRBackend}

    @assert datatype(obj) != "ecog" "Use plot_locs_ecog() for ECoG data."
    @assert (d == 2 || d == 3) "d must be 2 or 3."

    ch = get_channel(obj, ch=ch)
    chs = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    ch = collect(1:nrow(locs))

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
        p = plot_locs_nirs(obj.locs,
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
    end

    if d == 3
        if interactive
            iplot_locs3d(locs,
                         ch=ch,
                         selected=selected,
                         ch_labels=ch_labels,
                         head_labels=head_labels,
                         mono=mono)
            return nothing
        else
            p = plot_locs3d(locs,
                            ch=ch,
                            selected=selected,
                            ch_labels=ch_labels,
                            head_labels=head_labels,
                            mono=mono)
        end
    else
        p = plot_locs(locs,
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
