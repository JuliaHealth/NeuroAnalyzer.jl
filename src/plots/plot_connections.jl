export plot_connections

"""
    plot_connections(obj; <keyword arguments>)

Plot weights at electrode positions.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `connections::Matrix{<:Real}`: matrix of connections weights
- `threshold::Real`: plot all connection above threshold
- `threshold_type::Symbol=:g`: rule for thresholding: = (`:eq`), ≥ (`:geq`), ≤ (`:leq`), > (`:g`), < (`:l`)
- `weights::Bool=true`: weight line widths and alpha based on connection value
- `channel_labels::Bool=false`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `title::String=""`: plot title
- `cart::Bool=false`: if true, use Cartesian x and y coordinates, otherwise use polar radius and theta coordinates
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_connections(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), connections::Matrix{<:Real}, threshold::Real, threshold_type::Symbol=:g, weights::Bool=true, channel_labels::Bool=true, head_labels::Bool=false, mono::Bool=false, head_details::Bool=true, plot_size::Int64=800, title::String="", cart::Bool=false, kwargs...)

    _has_locs(obj) == false && throw(ArgumentError("Electrode locations not available, use load_locs() or add_locs() first."))
    _check_var(threshold_type, [:eq, :geq, :leq, :g, :l], "threshold_type")

    # remove non-signal channels
    obj_tmp = deepcopy(obj)
    keep_channel_type!(obj_tmp, type=Symbol(obj_tmp.header.recording[:data_type]))

    _check_channels(obj, ch, Symbol(obj.header.recording[:data_type]))
    typeof(ch) == Int64 && throw(ArgumentError("≥ 2 channels are required."))

    p = plot_connections(obj_tmp.locs, connections=connections, ch=ch, threshold=threshold, threshold_type=threshold_type, weights=weights, channel_labels=channel_labels, head_labels=head_labels, mono=mono, plot_size=plot_size, head_details=head_details, cart=cart)

    Plots.plot!(p, title=title; kwargs)

    return p

end

"""
    plot_connections(obj; <keyword arguments>)

Plot connections between channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Vector{Int64}, AbstractRange}`: channel(s) to plot
- `connections::Matrix{<:Real}`: matrix of connections weights
- `threshold::Real`: plot all connection above threshold
- `threshold_type::Symbol=:g`: rule for thresholding: = (`:eq`), ≥ (`:geq`), ≤ (`:leq`), > (`:g`), < (`:l`)
- `weights::Bool=true`: weight line widths and alpha based on connection value
- `channel_labels::Bool=false`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `cart::Bool=false`: if true, use Cartesian x and y coordinates, otherwise use polar radius and theta coordinates
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_connections(locs::DataFrame; ch::Union{Vector{Int64}, AbstractRange}, connections::Matrix{<:Real}, threshold::Real, threshold_type::Symbol=:g, weights::Bool=true, channel_labels::Bool=true, head_labels::Bool=false, mono::Bool=false, head_details::Bool=true, plot_size::Int64=800, cart::Bool=false, kwargs...)

    size(connections, 1) == length(ch) || throw(ArgumentError("Length of channel and number of connections rows must be equal."))
    _check_var(threshold_type, [:eq, :geq, :leq, :g, :l], "threshold_type")

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

    if cart == false
        loc_x = zeros(size(locs, 1))
        loc_y = zeros(size(locs, 1))
        for idx in 1:size(locs, 1)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
    else
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
    end
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    for idx in eachindex(locs[!, :labels])
        if idx in ch
            p = Plots.plot!((loc_x[idx], loc_y[idx]),
                            color=:gray,
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
            if idx in ch
                Plots.plot!(annotation=(loc_x[idx], loc_y[idx] + 0.05, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
        end
    end

    hd = _draw_head(p, head_labels=head_labels, head_details=head_details)
    p = Plots.plot!(hd)

    m_tmp = normalize_n(connections)

    for idx1 in 1:size(connections, 1)
        for idx2 in 1:size(connections, 1)
            if threshold_type === :g
                if connections[idx1, idx2] > threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            elseif threshold_type === :l
                if connections[idx1, idx2] < threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            elseif threshold_type === :eq
                if connections[idx1, idx2] == threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            elseif threshold_type === :leq
                if connections[idx1, idx2] <= threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            elseif threshold_type === :geq
                if connections[idx1, idx2] >= threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            end
        end
    end

    return p
    
end
