export ampdiff

"""
    ampdiff(s; <keyword arguments>)

Calculate amplitude difference to reference mean: amplitude difference between each channel and mean amplitude of reference channels.

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)
- `ch::Union{Int64, Vector{Int64}}=size(s, 1)`: indices of reference channels; default is all channels, for each analyzed channel, that channel itself is excluded from the reference mean

# Returns

- `ad::Array{Float64, 3}`: amplitude difference, shape `(channels, samples, epochs)`
"""
function ampdiff(
    s::AbstractArray;
    ch::Union{Int64, Vector{Int64}} = _c(size(s, 1))
)::Array{Float64, 3}

    # validate shape and that all requested channel indices are in bounds
    _chk3d(s)
    _check_channels(s, ch)

    # pre-allocate output
    ad = similar(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ref_ch = setdiff(ch, ch_idx)
        amp_ref = dropdims(mean(@view(s[ref_ch, :, ep_idx]), dims = 1), dims = 1)
        ad[ch_idx, :, ep_idx] .= @view(s[ch_idx, :, ep_idx]) .- amp_ref
    end

    return ad

end

"""
    ampdiff(obj; <keyword arguments>)

Calculate amplitude difference to reference mean: amplitude difference between each channel and mean amplitude of reference channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: reference channel name(s)

# Returns

- `ad::Array{Float64, 3}`: amplitude difference, shape `(channels, samples, epochs)`
"""
function ampdiff(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    ad = ampdiff(@view(obj.data[ch, :, :]), ch = _c(length(ch)))

    return ad

end
