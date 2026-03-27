export ampdiff

"""
    ampdiff(s)

Calculate amplitude difference to reference mean: amplitude difference between each channel and mean amplitude of reference channels.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

- `Array{Float64, 3}`: amplitude difference, shape (channels, samples, epochs)
"""
function ampdiff(
    s::AbstractArray,
)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    amp_diff = similar(s, Float64)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ref_ch = setdiff(ch, ch_idx)
        amp_ref = dropdims(mean(@view(s[ref_ch, :, ep_idx]), dims = 1), dims = 1)
        amp_diff[ch_idx, :, ep_idx] .= @view(s[ch_idx, :, ep_idx]) .- amp_ref
    end

    return amp_diff
end

"""
    ampdiff(obj; <keyword arguments>)

Calculate amplitude difference to reference mean: amplitude difference between each channel and mean amplitude of reference channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: reference channel name(s)

# Returns

- `Array{Float64, 3}`: amplitude difference, shape (channels, samples, epochs)
"""
function ampdiff(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")

    return ampdiff(@view(obj.data[ch, :, :]))
end
