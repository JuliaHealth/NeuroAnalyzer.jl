export frqinst

"""
    frqinst(s; <keyword arguments>)

Calculate instantaneous frequency.

# Arguments

- `s::AbstractVector`

# Returns

- `f::Vector{Float64}`
"""
function frqinst(s::AbstractVector)::Vector{Float64}

    _, _, _, pha = hspectrum(s)
    f = 1 / (2 * pi) * derivative(DSP.unwrap(pha))

    return f

end

"""
    frqinst(s; <keyword arguments>)

Calculate instantaneous frequency.

# Arguments

- `s::AbstractVector`

# Returns

- `f::Array{Float64, 3}`
"""
function frqinst(s::AbstractArray)::Array{Float64, 3}

    _warn("frqinst() uses Hilbert transform, the signal should be narrowband for best results.")

    _chk3d(s)

    f = similar(s)
    @inbounds for ep_idx in axes(s, 3)
        Threads.@threads :greedy for ch_idx in axes(s, 1)
            f[ch_idx, :, ep_idx] = @views frqinst(s[ch_idx, :, ep_idx])
        end
    end

    return f

end

"""
    frqinst(obj; <keyword arguments>)

Calculate instantaneous frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `f::Array{Float64, 3}`
"""
function frqinst(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    f = @views frqinst(obj.data[ch, :, :])

    return f

end
