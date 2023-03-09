export plot_locs
export plot_locs3d

"""
    plot_locs(locs; <keyword arguments>)

Preview of channel locations. It uses polar `:loc_radius` and `:loc_theta` locations, which are translated into Cartesian x and y positions.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
- `channel_labels::Bool=true`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `grid::Bool=false`: draw grid, useful for locating positions
- `plot_size::Int64=400`: plot dimensions in pixels (size × size)

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs(locs::DataFrame; channel::Union{Int64, Vector{Int64}, AbstractRange}, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, channel_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, head_details::Bool=true, grid::Bool=false, plot_size::Int64=400)

    pal = mono == true ? :grays : :darktest

    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    if plot_size > 300
        marker_size = plot_size ÷ 75
        font_size = plot_size ÷ 75
    else
        marker_size = plot_size ÷ 50
        font_size = plot_size ÷ 50
        channel_labels = false
    end

    length(channel) > 64 && (font_size = plot_size ÷ 100)

    if grid == false
        p = Plots.plot(grid=false,
                       framestyle=:none,
                       palette=pal,
                       size=(plot_size, plot_size),
                       border=:none,
                       aspect_ratio=1,
                       right_margin=-30 * Plots.px,
                       bottom_margin=-20 * Plots.px,
                       top_margin=-30 * Plots.px,
                       left_margin=-50 * Plots.px,
                       xlim=(-1.22, 1.23),
                       ylim=(-1.1, 1.2))
    else
        p = Plots.plot(grid=true,
                       palette=pal,
                       size=(plot_size, plot_size),
                       aspect_ratio=1,
                       right_margin=-30 * Plots.px,
                       bottom_margin=-50 * Plots.px,
                       top_margin=-50 * Plots.px,
                       left_margin=-5 * Plots.px,
                       xticks=-1:0.1:1,
                       yticks=-1:0.1:1,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       xlim=(-1.22, 1.23),
                       ylim=(-1.1, 1.2))
    end
    hd = _draw_head(p, head_labels=head_labels, head_details=head_details)
    p = Plots.plot!(hd)

    for idx in eachindex(locs[!, :labels])
        if idx in channel
            if selected != 0
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                #seriestype=:scatter,
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            else
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                #seriestype=:scatter,
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            end
        end
        if idx in selected
            if mono != true
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=idx,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                #seriestype=:scatter,
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            else
                #p = Plots.plot!((loc_x[idx], loc_y[idx]),
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                #seriestype=:scatter,
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            end
        end
    end

    if channel_labels == true
        for idx in eachindex(locs[!, :labels])
            if idx in channel
                Plots.plot!(annotation=(loc_x[idx], loc_y[idx] + 0.075, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
            if idx in selected
                Plots.plot!(annotation=(loc_x[idx], loc_y[idx] + 0.075, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
        end
    end

    Plots.plot(p)

    return p
end

"""
    plot_locs3d(locs; <keyword arguments>)

3D interactive preview of channel locations. It uses spherical :loc_radius_sph, :loc_theta_sph and :loc_phi_sph locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
- `channel_labels::Bool=true`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `plot_size::Int64=800`: plot dimensions in pixels (plot_size×plot_size)
c

# Returns

- `fig::GLMakie.Figure`
"""
function plot_locs3d(locs::DataFrame; channel::Union{Int64, Vector{Int64}, AbstractRange}, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, channel_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, plot_size::Int64=800)

    # selected != 0 && length(intersect(channel, selected)) < length(selected) && throw(ArgumentError("channel must include selected."))
    # channel = setdiff(channel, selected)

    pal = mono == true ? :grays : :darktest

    loc_x = zeros(nrow(locs))
    loc_y = zeros(nrow(locs))
    loc_z = zeros(nrow(locs))

    for idx in 1:nrow(locs)
        loc_x[idx], loc_y[idx], loc_z[idx] = sph2cart(locs[idx, :loc_radius_sph], locs[idx, :loc_theta_sph], locs[idx, :loc_phi_sph])
    end

    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    x_lim = (-1.1, 1.1)
    y_lim = (-1.1, 1.1)
    z_lim = extrema(loc_z)

    marker_size = plot_size ÷ 40
    font_size = plot_size ÷ 40

    fig = Figure(; resolution=(plot_size, plot_size))
    ax = Axis3(fig[1, 1]; aspect=(1, 1, 0.5), perspectiveness=0.5, limits = (x_lim, y_lim, z_lim))
    # hidedecorations!(ax, grid=true, ticks=true)

    GLMakie.scatter!(ax, loc_x[channel], loc_y[channel], loc_z[channel], markersize=marker_size, color=:gray)
    if selected != 0
        if mono == true
            GLMakie.scatter!(ax, loc_x[selected], loc_y[selected], loc_z[selected], markersize=marker_size, color=:gray)
        else
            GLMakie.scatter!(ax, loc_x[selected], loc_y[selected], loc_z[selected], markersize=marker_size, color=:red)
        end
    end

    if channel_labels == true
        for idx in eachindex(locs[!, :labels])
            if idx in channel
                GLMakie.text!(ax, locs[!, :labels][idx], position=(loc_x[idx], loc_y[idx], loc_z[idx]), fontsize=font_size)
            end
            if idx in selected
                GLMakie.text!(ax, locs[!, :labels][idx], position=(loc_x[idx], loc_y[idx], loc_z[idx]), fontsize=font_size)
            end
        end
    end

    if head_labels == true
        GLMakie.text!(ax, "Nz", position=(0, 1.025, 0), fontsize = font_size)
        GLMakie.text!(ax, "Iz", position=(0, -1.025, 0), fontsize = font_size)
        GLMakie.text!(ax, "LPA", position=(-1.025, 0, 0), fontsize = font_size)
        GLMakie.text!(ax, "RPA", position=(1.025, 0, 0), fontsize = font_size)
        GLMakie.text!(ax, "top", position=(0, 0, 1.025), fontsize = font_size)
    end
    fig

    return fig
end

"""
    plot_locs(obj; <keyword arguments>)

Preview of channel locations.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: which channel should be highlighted
- `channel_labels::Bool=true`: plot channel labels
- `head::Bool`=true: plot head
- `head_labels::Bool=false`: plot head labels
- `plot_size::Int64=400`: plot dimensions in pixels (plot_size×plot_size)
- `head_details::Bool=true`: draw nose and ears
- `mono::Bool=false`: use color or grey palette
- `threed::Bool=false`: 3-dimensional plot
- `grid::Bool=false`: draw grid, useful for locating positions
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), selected::Union{Int64, Vector{Int64}, AbstractRange}=0, channel_labels::Bool=true, head::Bool=true, head_labels::Bool=false, plot_size::Int64=400, head_details::Bool=true, mono::Bool=false, threed::Bool=false, grid::Bool=false, kwargs...)

    obj.header.has_locs == false && throw(ArgumentError("Channel locations not available, use load_locs() or add_locs() first."))

    # select channels, default is all channels
    _check_channels(obj, channel, Symbol(obj.header.recording[:data_type]))
    selected != 0 && _check_channels(obj, selected)

    if threed == false
        p = plot_locs(obj.locs, channel=channel, selected=selected, channel_labels=channel_labels, head_labels=head_labels, mono=mono, head_details=head_details, plot_size=plot_size, grid=grid)
    else
        p = plot_locs3d(obj.locs, channel=channel, selected=selected, channel_labels=channel_labels, head_labels=head_labels, mono=mono, plot_size=plot_size)
    end

    return p
end
