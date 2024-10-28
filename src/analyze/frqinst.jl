export frqinst

"""
    frqinst(s; <keyword arguments>)

Calculate instantaneous frequency.

# Arguments

- `s::Vector{<:Real}`

# Returns

- `f::Vector{Float64}`
"""
function frqinst(s::Vector{<:Real})::Vector{Float64}

    _, _, _, pha = hspectrum(s)
    f = 1 / (2 * pi) * derivative(DSP.unwrap(pha))

    return f

end

"""
    frqinst(s; <keyword arguments>)

Calculate instantaneous frequency.

# Arguments

- `s::Array{<:Real, 3}`

# Returns

- `f::Array{Float64, 3}`
"""
function frqinst(s::Array{<:Real, 3})::Array{Float64, 3}

    _warn("frqinst() uses Hilbert transform, the signal should be narrowband for best results.")

    f = similar(s)
    @inbounds for ep_idx in axes(s, 3)
        Threads.@threads for ch_idx in axes(s, 1)
            f[ch_idx, :, ep_idx] = frqinst(s[ch_idx, :, ep_idx])
        end
    end

    return f

end

"""
    frqinst(obj; <keyword arguments>)

Calculate instantaneous frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

- `f::Array{Float64, 3}`
"""
function frqinst(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})::Array{Float64, 3}

    ch = get_channel(obj, ch=ch)
    f = frqinst(obj.data[ch, :, :])

    return f

end
