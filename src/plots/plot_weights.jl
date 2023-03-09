export plot_weights
export plot_weights

"""
    plot_weights(locs; <keyword arguments>)

Plot weights at electrode positions. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
- `weights::Vector{<:Real}=[]`: weights vector
- `channel_labels::Bool=true`: plot channel_labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `plot_size::Int64=400`: plot dimensions in pixels (size × size)

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_weights(locs::DataFrame; channel::Union{Int64, Vector{Int64}, AbstractRange}, weights::Vector{<:Real}=[], channel_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, head_details::Bool=true, plot_size::Int64=400)

    length(weights) > length(channel) && throw(ArgumentError("Number of weights must be ≤ number of channels to plot ($(length(channel)))."))
    length(weights) < 1 && throw(ArgumentError("weights must contain at least one value."))

    # selected != 0 && length(intersect(channel, selected)) < length(selected) && throw(ArgumentError("channel must include selected."))
    # channel = setdiff(channel, selected)

    pal = mono == true ? :grays : :darktest

    marker_size = plot_size ÷ 100
    font_size = plot_size ÷ 100

    p = Plots.plot(grid=true,
                   framestyle=:none,
                   palette=pal,
                   size=(plot_size, plot_size),
                   markerstrokewidth=0,
                   border=:none,
                   aspect_ratio=1,
                   margins=-plot_size * Plots.px,
                   titlefontsize=plot_size ÷ 50)

    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    for idx in eachindex(locs[!, :labels])
        if idx in channel
            p = Plots.plot!((loc_x[idx], loc_y[idx]),
                            color=:black,
                            seriestype=:scatter,
                            grid=true,
                            label="",
                            markersize=marker_size,
                            markerstrokewidth=0,
                            markerstrokealpha=0)
        end
    end

    if channel_labels == true
        for idx in eachindex(locs[!, :labels])
            if idx in channel
                Plots.plot!(annotation=(loc_x[idx], loc_y[idx] + 0.05, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
        end
    end

    for idx in eachindex(locs[!, :labels])
        if idx in channel
            Plots.plot!(annotation=(loc_x[idx], loc_y[idx], Plots.text(string(weights[idx]), pointsize=font_size)))
        end
    end

    hd = _draw_head(p, head_labels=head_labels, head_details=head_details)
    p = Plots.plot!(hd)

    Plots.plot(p)

    return p
end

"""
    plot_weights(obj; <keyword arguments>)

Plot weights at electrode positions. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `weights::Matrix{<:Real}`: matrix of weights
- `channel_labels::Bool=false`: plot channel_labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `title::String=""`: plot title
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_weights(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), weights::Vector{<:Real}, channel_labels::Bool=true, head_labels::Bool=false, mono::Bool=false, head_details::Bool=true, plot_size::Int64=800, title::String="", kwargs...)

    obj.header.has_locs == false && throw(ArgumentError("Electrode locations not available, use load_electrodes() or add_electrodes() first."))

    # remove non-signal channels
    obj_tmp = deepcopy(obj)
    keep_channel_type!(obj_tmp, type=Symbol(obj_tmp.header.recording[:data_type]))

    _check_channels(obj, channel, Symbol(obj.header.recording[:data_type]))
    typeof(channel) == Int64 && throw(ArgumentError("≥ 2 channels are required."))

    p = plot_weights(obj_tmp.locs, weights=weights, channel=channel, channel_labels=channel_labels, head_labels=head_labels, mono=mono, plot_size=plot_size, head_details=head_details)

    Plots.plot!(p, title=title; kwargs)

    return p
end
