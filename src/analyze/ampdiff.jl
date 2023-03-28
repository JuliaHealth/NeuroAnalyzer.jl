export ampdiff

"""
    ampdiff(s; ch)

Calculate amplitude difference between each channel and mean amplitude of reference channels `ch`.

# Arguments

- `s::AbstractArray`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=size(s, 1)`: index of reference channels, default is all channels except the analyzed one

# Returns
 
- `ad::Array{Float64, 3}`
"""
function ampdiff(s::AbstractArray; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=size(s, 1))

    _check_channels(s, ch)

    ch_n = length(ch)
    ep_n = size(s, 3)

    ad = similar(s)

    @inbounds @simd for ep_idx in 1:ep_n
        for ch_idx in 1:ch_n
            ref_ch = setdiff(ch, ch_idx)
            amp_ref = @views vec(mean(s[ref_ch, :, ep_idx], dims=1))
            ad[ch_idx, :, ep_idx] = @views s[ch[ch_idx], :, ep_idx] - amp_ref
        end
    end

    return ad
end


"""
    ampdiff(obj; ch)

Calculate amplitude difference between each channel and mean amplitude of reference channels `ch`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of reference channels, default is all signal channels except the analyzed one

# Returns
 
- `ad::Array{Float64, 3}`
"""
function ampdiff(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)

    ad = @views ampdiff(obj.data[ch, :, :], ch=ch)

    return ad

end
