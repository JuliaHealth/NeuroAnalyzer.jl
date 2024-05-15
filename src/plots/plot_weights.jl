export plot_weights

"""
    plot_weights(locs; <keyword arguments>)

Plot weights at electrode positions.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `weights::Vector{<:Real}=[]`: weights vector
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
- `ch_labels::Bool=true`: plot channel labels
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
function plot_weights(locs::DataFrame; weights::Vector{<:Real}=[], ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs), ch_labels::Bool=true, head::Bool=true, head_labels::Bool=false, mono::Bool=false, large::Bool=true, cart::Bool=false, plane::Symbol=:xy, title::String="")

    @assert length(weights) <= length(ch) "Number of weights must be ≤ number of channels to plot ($(length(ch)))."
    @assert length(weights) >= 1 "weights must contain at least one value."

    _check_var(plane, [:xy, :yz, :xz], "plane")

    pal = mono ? :grays : :darktest

    locs = locs[ch, :]

    if plane === :xy
        if large
            head_shape = FileIO.load(joinpath(res_path, "head_t_large.png"))
        else
            head_shape = FileIO.load(joinpath(res_path, "head_t_small.png"))
        end
        if cart == false
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
        if cart == false
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
        if cart == false
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

    for idx in eachindex(locs[!, :labels])
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
        for idx in eachindex(locs[!, :labels])
            if idx in ch
                Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + 1, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
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
            p = Plots.plot!(annotations=(fid_loc_x, fid_loc_y, Plots.text(fid_names[idx], pointsize=font_size + 2)))
        end
    end
    for idx in eachindex(locs[ch, :labels])
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

    large == false && (title = "")
    Plots.plot!(p,
                title=title,
                titlefontsize=10)


    Plots.plot!(p)

    return p

end

"""
    plot_weights(obj; <keyword arguments>)

Plot weights at electrode positions.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `weights::Matrix{<:Real}`: matrix of weights
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
function plot_weights(obj::NeuroAnalyzer.NEURO; weights::Vector{<:Real}, ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type=datatype(obj)), ch_labels::Bool=true, head::Bool=true, head_labels::Bool=false, mono::Bool=false, large::Bool=true, cart::Bool=false, plane::Symbol=:xy, title::String="")

    @assert _has_locs(obj) "Electrode locations not available, use load_locs() or add_locs() first."

    _check_channels(obj, ch)

    p = plot_weights(obj.locs, weights=weights, ch=ch, ch_labels=ch_labels, head=head, head_labels=head_labels, large=large, mono=mono, cart=cart, plane=plane, title=title)

    return p
    
end
