export frqinst

"""
    frqinst(s; <keyword arguments>)

Calculate instantaneous frequencies.

# Arguments

- `s::AbstractVector`

# Returns

- `f::Vector{Float64}`: instantaneous frequencies (in rad/s)
"""
function frqinst(s::AbstractVector)::Vector{Float64}

    _, _, _, ph = htransform(s)

    f = NeuroAnalyzer.derivative(DSP.unwrap(ph)) / (2 * pi)

    return f

end

"""
    frqinst(s; <keyword arguments>)

Calculate instantaneous frequencies.

# Arguments

- `s::AbstractVector`

# Returns

- `f::Array{Float64, 3}`: instantaneous frequencies (in rad/s)
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

Calculate instantaneous frequencies.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `f::Array{Float64, 3}`: instantaneous frequencies (in Hz)
"""
function frqinst(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    f = @views frqinst(obj.data[ch, :, :]) .* sr(obj)

    return f

end
