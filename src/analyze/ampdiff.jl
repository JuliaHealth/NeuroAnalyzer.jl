export ampdiff

"""
    ampdiff(s; <keyword arguments>)

Calculate amplitude difference between each channel and mean amplitude of reference channels.

# Arguments

- `s::AbstractArray`
- `ch::Union{Int64, Vector{Int64}}=size(s, 1)`: index of reference channels, default is all channels except the analyzed one

# Returns

- `ad::Array{Float64, 3}`
"""
function ampdiff(s::AbstractArray; ch::Union{Int64, Vector{Int64}}=_c(size(s, 1)))::Array{Float64, 3}

    _chk3d(s)
    _check_channels(s, ch)

    ad = similar(s)

    @inbounds for ep_idx in axes(s, 3)
        for ch_idx in axes(s, 1)
            ref_ch = setdiff(ch, ch_idx)
            amp_ref = @views vec(mean(s[ref_ch, :, ep_idx], dims=1))
            ad[ch_idx, :, ep_idx] = @views s[ch[ch_idx], :, ep_idx] - amp_ref
        end
    end

    return ad

end

"""
    ampdiff(obj; <keyword arguments>)

Calculate amplitude difference between each channel and mean amplitude of reference channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: index of reference channels

# Returns

- `ad::Array{Float64, 3}`
"""
function ampdiff(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")

    ad = @views ampdiff(obj.data[ch, :, :], ch=ch)

    return ad

end
