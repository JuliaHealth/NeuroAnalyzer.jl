export frqinst

"""
    frqinst(s; fs)

Calculate instantaneous frequency.

# Arguments

- `s::AbstractVector`
- `fs::Int64`

# Returns

- `f::Vector{Float64}`
"""
function frqinst(s::AbstractVector; fs::Int64)

    @assert fs >= 1 "fs must be ≥ 1."

    _, _, _, pha = hspectrum(s)
    f = fs / (2 * pi) * derivative(DSP.unwrap(pha))

    return f

end

"""
    frqinst(s; fs)

Calculate instantaneous frequency.

# Arguments

- `s::AbstractVector`
- `fs::Int64`

# Returns

- `f::Array{Float64, 2}`
"""
function frqinst(s::AbstractArray; fs::Int64)

    _warn("frqinst() uses Hilbert transform, the signal should be narrowband for best results.")

    ch_n = size(s, 1)
    ep_len = size(s, 2)
    ep_n = size(s, 3)
    
    f = zeros(ch_n, ep_len, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            f[ch_idx, :, ep_idx] = @views frqinst(s[ch_idx, :, ep_idx], fs=fs)
        end
    end

    return f

end

"""
    frqinst(obj; channel)

Calculate instantaneous frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `f::Array{Float64, 3}`
"""
function frqinst(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)

    f = @views frqinst(obj.data[channel, :, :], fs=sr(obj))

    return f

end
