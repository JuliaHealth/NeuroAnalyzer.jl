export psd

"""
    psd(s; fs, norm, mt, nt)

Calculate power spectrum density.

# Arguments
- `s::Vector{Float64}`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `pw::Vector{Float64}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd(s::Vector{Float64}; fs::Int64, norm::Bool=false, mt::Bool=false, nt::Int64=8)

    nt < 1 && throw(ArgumentError("nt must be ≥ 1."))
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    # for short signals use multi-tapered periodogram
    length(s) < 4 * fs && (mt = true)

    if mt == true
        p = mt_pgram(s, fs=fs, nw=(nt÷2+1), ntapers=nt)
    else
        p = welch_pgram(s, 4*fs, fs=fs)
    end

    pw = power(p)
    pf = Vector(freq(p))
    pw = pw[1:length(pf)]

    # replace powers at extreme frequencies
    pw[1] = pw[2]
    pw[end] = pw[end - 1]
    
    norm == true && (pw = pow2db.(pw))

    return (pw=pw, pf=pf)

end

"""
    psd(s; fs, norm, mt, nt)

Calculate power spectrum density.

# Arguments

- `s::AbstractMatrix`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `pw::Array{Float64, 2}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd(s::AbstractMatrix; fs::Int64, norm::Bool=false, mt::Bool=false, nt::Int64=8)

    # for short signals use multi-tapered periodogram
    if size(s, 2) < 4 * fs
        mt = true
        _info("Using multi-tapered periodogram.")
    end

    ch_n = size(s, 1)
    _, pf = psd(s[1, :], fs=fs, norm=norm, mt=mt, nt=nt)

    pw = zeros(ch_n, length(pf))

    @inbounds @simd for ch_idx in 1:ch_n
        pw[ch_idx, :], _ = psd(s[ch_idx, :], fs=fs, norm=norm, mt=mt, nt=nt)
    end
    
    return (pw=pw, pf=pf)

end

"""
    psd(s; fs, norm, mt, nt)

Calculate power spectrum density.

# Arguments
- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `pw::Array{Float64, 3}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd(s::AbstractArray; fs::Int64, norm::Bool=false, mt::Bool=false, nt::Int64=8)

    # for short signals use multi-tapered periodogram
    if size(s, 2) < 4 * fs
        mt = true
        _info("Using multi-tapered periodogram.")
    end

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, pf = psd(s[1, :, 1], fs=fs, norm=norm, mt=mt, nt=nt)

    pw = zeros(ch_n, length(pf), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pw[ch_idx, :, ep_idx], _ = psd(s[ch_idx, :, ep_idx], fs=fs, norm=norm, mt=mt, nt=nt)
        end
    end
    
    return (pw=pw, pf=pf)

end

"""
    psd(obj; ch, norm, mt)

Calculate power spectrum density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `pw::Array{Float64, 3}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), norm::Bool=false, mt::Bool=false, nt::Int64=8)

    _check_channels(obj, ch)

    if length(ch) == 1
        pw, pf = psd(reshape(obj.data[ch, :, :], length(ch), :, size(obj.data[ch, :, :], 2)), fs=sr(obj), norm=norm, mt=mt, nt=nt)
    else
        pw, pf = psd(obj.data[ch, :, :], fs=sr(obj), norm=norm, mt=mt, nt=nt)
    end

    return (pw=pw, pf=pf)

end
