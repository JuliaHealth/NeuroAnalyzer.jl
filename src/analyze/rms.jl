export rms
export msa
export amp

"""
    rms(s)

Calculate Root Mean Square.

# Arguments

- `s::AbstractVector`

# Returns

- `rms::Float64`
"""
function rms(s::AbstractVector)

    # return sqrt(1/length(s) * sum(s.^2))
    return norm(s) / sqrt(length(s))

end

"""
    rms(s)

Calculate Root Mean Square.

# Arguments

- `s::AbstractArray`

# Returns

- `r::Array{Float64, 2}`
"""
function rms(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    
    r = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            r[ch_idx, ep_idx] = @views rms(s[ch_idx, :, ep_idx])
        end
    end

    return r

end

"""
    rms(obj; channel)

Calculate Root Mean Square.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `r::Array{Float64, 2}`
"""
function rms(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)
    length(ch) == 1 && (ch = [ch])

    r = @views rms(obj.data[ch, :, :])

    return r

end

"""
    msa(s)

Calculate Mean Square Amplitude.

# Arguments

- `s::AbstractVector`

# Returns

- `msa::Float64`
"""
function msa(s::AbstractVector)

    return 1/length(s) * sum(s.^2)

end

"""
    msa(s)

Calculate Mean Square Amplitude.

# Arguments

- `s::AbstractArray`

# Returns

- `r::Array{Float64, 2}`
"""
function msa(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    
    r = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            r[ch_idx, ep_idx] = @views msa(s[ch_idx, :, ep_idx])
        end
    end

    return r

end

"""
    msa(obj; channel)

Calculate Mean Square Amplitude.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `r::Array{Float64, 2}`
"""
function msa(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)
    length(ch) == 1 && (ch = [ch])

    r = @views msa(obj.data[ch, :, :])

    return r

end

export amp

"""
    amp(s)

Calculate amplitudes.

# Arguments

- `s::AbstractVector`

# Returns
 
Named tuple containing:
- `p::Float64`: peak amplitude
- `r::Float64`: RMS amplitude
- `p2p::Float64`: peak-to-peak amplitude
- `semi_p2p::Float64`: half of the peak-to-peak amplitude
- `rmsa::Float64`: root mean square amplitude
"""
function amp(s::AbstractVector)

    p = maximum(abs.(s))
    r = p / sqrt(2)
    p2p = abs(maximum(s)) + abs(minimum(s))
    semi_p2p = p2p / 2
    rmsa = p2p / sqrt(2) 

    return (p=p, r=r, p2p=p2p, semi_p2p=semi_p2p, rmsa=rmsa)

end

"""
    amp(s)

Calculate amplitudes.

# Arguments

- `s::AbstractArray`

# Returns
 
Named tuple containing:
- `p::Matrix{Float64}`: peak amplitude
- `r::Matrix{Float64}`: RMS amplitude
- `p2p::Matrix{Float64}`: peak-to-peak amplitude
- `semi_p2p::Matrix{Float64}`: half of the peak-to-peak amplitude
- `rmsa::Matrix{Float64}`: root mean square amplitude
"""
function amp(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    p = zeros(ch_n, ep_n)
    r = zeros(ch_n, ep_n)
    p2p = zeros(ch_n, ep_n)
    semi_p2p = zeros(ch_n, ep_n)
    rmsa = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        for ch_idx in 1:ch_n
            p[ch_idx, ep_idx], r[ch_idx, ep_idx], p2p[ch_idx, ep_idx], semi_p2p[ch_idx, ep_idx], rmsa[ch_idx, ep_idx] = @views amp(s[ch_idx, :, ep_idx])
        end
    end

    return (p=p, r=r, p2p=p2p, semi_p2p=semi_p2p, rmsa=rmsa)

end

"""
    amp(obj; ch)

Calculate amplitudes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of reference channels, default is all signal channels except the analyzed one

# Returns
 
Named tuple containing:
- `p::Matrix{Float64}`: peak amplitude
- `r::Matrix{Float64}`: RMS amplitude
- `p2p::Matrix{Float64}`: peak-to-peak amplitude
- `semi_p2p::Matrix{Float64}`: half of the peak-to-peak amplitude
- `rmsa::Matrix{Float64}`: root mean square amplitude
"""
function amp(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)
    length(ch) == 1 && (ch = [ch])

    p, r, p2p, semi_p2p, rmsa = @views amp(obj.data[ch, :, :])

    return (p=p, r=r, p2p=p2p, semi_p2p=semi_p2p, rmsa=rmsa)

end
