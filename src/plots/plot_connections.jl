export plot_connections

"""
    plot_connections(obj; <keyword arguments>)

Plot connections between channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `connections::Matrix{<:Real}`: matrix of connections weights (channels by channels)
- `threshold::Real=0`: threshold for plotting, see below
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: plot if connection weight is equal to threshold
    - `:neq`: plot if connection weight is not equal to threshold
    - `:geq`: plot if connection weight is ≥ to threshold
    - `:leq`: plot if connection weight is ≤ to threshold
    - `:g`: plot if connection weight is > to threshold
    - `:l`: plot if connection weight is < to threshold
- `weights::Bool=true`: weight line widths and alpha based on connection value, if false connections values will be drawn
- `ch::Union{String, Vector{String}}=1:nrow(locs)`: list of channels, default is all channels
- `ch_labels::Bool=false`: plot channel labels
- `head::Bool=true`: draw head
- `head_labels::Bool=false`: plot head labels
- `mono::Bool=false`: use color or gray palette
- `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `plane::Symbol=:xy`: which plane to plot:
    - `:xy`: horizontal (top)
    - `:xz`: coronary (front)
    - `:yz`: sagittal (side)
- `title::String=""`: plot title

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_connections(locs::DataFrame; connections::Matrix{<:Real}, threshold::Real=0, threshold_type::Symbol=:neq, weights::Bool=true, ch::Union{String, Vector{String}}=1:nrow(locs), ch_labels::Bool=false, head::Bool=true, head_labels::Bool=false, mono::Bool=false, large::Bool=true, cart::Bool=false, plane::Symbol=:xy, title::String="")

    @assert size(connections, 1) == length(ch) "Length of channel and number of connections rows must be equal."
    _check_var(threshold_type, [:eq, :neq, :geq, :leq, :g, :l], "threshold_type")

    _check_var(plane, [:xy, :yz, :xz], "plane")

    pal = mono ? :grays : :darktest

    locs = locs[ch, :]

    if plane === :xy
        if large
            head_shape = FileIO.load(joinpath(res_path, "head_t_large.png"))
        else
            head_shape = FileIO.load(joinpath(res_path, "head_t_small.png"))
        end
        if !cart
            loc_x = zeros(nrow(locs))
            loc_y = zeros(nrow(locs))
            for idx in 1:nrow(locs)
                loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
            end
        else
            loc_x = locs[!, :loc_x]
            loc_y = locs[!, :loc_y]
        end
    elseif plane === :xz
        if large
            head_shape = FileIO.load(joinpath(res_path, "head_f_large.png"))
        else
            head_shape = FileIO.load(joinpath(res_path, "head_f_small.png"))
        end
        if !cart
            loc_x = zeros(nrow(locs))
            loc_y = zeros(nrow(locs))
            for idx in 1:nrow(locs)
                loc_x[idx], _, loc_y[idx] = sph2cart(locs[!, :loc_radius_sph][idx], locs[!, :loc_theta_sph][idx], locs[!, :loc_phi_sph][idx])
            end
        else
            loc_x = locs[!, :loc_x]
            loc_y = locs[!, :loc_z]
        end
    elseif plane === :yz
        if large
            head_shape = FileIO.load(joinpath(res_path, "head_s_large.png"))
        else
            head_shape = FileIO.load(joinpath(res_path, "head_s_small.png"))
        end
        if !cart
            loc_x = zeros(nrow(locs))
            loc_y = zeros(nrow(locs))
            for idx in 1:nrow(locs)
                _, loc_x[idx], loc_y[idx] = sph2cart(locs[!, :loc_radius_sph][idx], locs[!, :loc_theta_sph][idx], locs[!, :loc_phi_sph][idx])
            end
        else
            loc_x = locs[!, :loc_y]
            loc_y = locs[!, :loc_z]
        end
    end

    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    if head
        xt = (linspace(0, size(head_shape, 1), 25), string.(-1.2:0.1:1.2))
        yt = (linspace(0, size(head_shape, 2), 25), string.(1.2:-0.1:-1.2))
        xl = (0, size(head_shape, 1))
        yl = (0, size(head_shape, 2))
    else
        xt = (-1.2:0.1:1.2)
        yt = (1.2:-0.1:-1.2)
        xl = (-1.2, 1.2)
        yl = (-1.2, 1.2)
    end

    origin = size(head_shape) ./ 2
    if large
        marker_size = 10
        font_size = 6
        loc_x = @. round(origin[1] + (loc_x * 250), digits=2)
        loc_y = @. round(origin[2] - (loc_y * 250), digits=2)
    else
        marker_size = 4
        font_size = 4
        ch_labels = false
        grid = false
        loc_x = @. round(origin[1] + (loc_x * 100), digits=2)
        loc_y = @. round(origin[2] - (loc_y * 100), digits=2)
    end

    ma = 1.0
    ch_labels && (ma = 0.75)

    if large
        p = Plots.plot(grid=false,
                       framestyle=:none,
                       border=:none,
                       palette=pal,
                       aspect_ratio=1,
                       size=title == "" ? size(head_shape) .+ 95 : size(head_shape) .+ 75,
                       right_margin=-50*Plots.px,
                       bottom_margin=5*Plots.px,
                       top_margin=title == "" ? 50*Plots.px : 10*Plots.px,
                       left_margin=-50*Plots.px,
                       titlefontsize=font_size,
                       xlims=xl,
                       ylims=yl,
                       title=title)
    else
        p = Plots.plot(grid=false,
                       framestyle=:none,
                       border=:none,
                       palette=pal,
                       aspect_ratio=1,
                       size=size(head_shape) .+ 2,
                       right_margin=-50*Plots.px,
                       bottom_margin=-50*Plots.px,
                       top_margin=-50*Plots.px,
                       left_margin=-100*Plots.px,
                       titlefontsize=font_size,
                       xlims=xl,
                       ylims=yl,
                       title=title)
    end

    head && (p = Plots.plot!(head_shape))

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
    if ch_labels
        for idx in eachindex(locs[!, :label])
            if idx in ch
                Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + 1, Plots.text(locs[!, :label][idx], pointsize=font_size)))
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
                fid_loc_x = @. origin[1] + (fid_loc_x * 250)
                fid_loc_y = @. origin[2] - (fid_loc_y * 250)
            else
                fid_loc_x = @. origin[1] - (fid_loc_x * 100)
                fid_loc_y = @. origin[2] - (fid_loc_y * 100)
            end
            p = Plots.plot!(annotations=(fid_loc_x, fid_loc_y, Plots.text(fid_names[idx], pointsize=font_size+2 + 2)))
        end
    end

    m_tmp = normalize_n(abs.(connections))

    for idx1 in 1:size(connections, 1)
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

    !large && (title = "")
    Plots.plot!(p,
                title=title,
                titlefontsize=10)

    Plots.plot!(p)

    return p

end

"""
    plot_connections(obj; <keyword arguments>)

Plot weights at electrode positions.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `connections::Matrix{<:Real}`: matrix of connections weights
- `threshold::Real=0`: threshold value
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: plot if connection weight is equal to threshold
    - `:neq`: plot if connection weight is not equal to threshold
    - `:geq`: plot if connection weight is ≥ to threshold
    - `:leq`: plot if connection weight is ≤ to threshold
    - `:g`: plot if connection weight is > to threshold
    - `:l`: plot if connection weight is < to threshold
- `weights::Bool=true`: weight line widths and alpha based on connection value
- `ch::Union{String, Vector{String}}=get_channel(obj, type=datatype(obj))`: list of channels, default is all channels
- `ch_labels::Bool=false`: plot ch_labels
- `head::Bool=true`: draw head
- `head_labels::Bool=false`: plot head labels
- `mono::Bool=false`: use color or gray palette
- `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `plane::Symbol=:xy`: which plane to plot:
    - `:xy`: horizontal (top)
    - `:xz`: coronary (front)
    - `:yz`: sagittal (side)
- `title::String=""`: plot title

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_connections(obj::NeuroAnalyzer.NEURO; connections::Matrix{<:Real}, threshold::Real=0, threshold_type::Symbol=:neq, weights::Bool=true, ch::Union{String, Vector{String}}=get_channel(obj, type=datatype(obj)), ch_labels::Bool=false, head::Bool=true, head_labels::Bool=false, mono::Bool=false, large::Bool=true, cart::Bool=false, plane::Symbol=:xy, title::String="")

    @assert _has_locs(obj) "Electrode locations not available, use load_locs() or add_locs() first."
    _check_var(threshold_type, [:eq, :neq, :geq, :leq, :g, :l], "threshold_type")

    _check_channels(obj, ch, datatype(obj))
    ch = _ch_idx(obj, ch)

    p = plot_connections(obj.locs, connections=connections, threshold=threshold, threshold_type=threshold_type, weights=weights, ch=ch, ch_labels=ch_labels, head=head, head_labels=head_labels, large=large, mono=mono, cart=cart, plane=plane, title=title)

    return p

end