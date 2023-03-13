export psd_rel

"""
    psd_rel(s; fs, norm, mt, nt, f)

Calculate relative power spectrum density.

# Arguments
- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers
- `f::Union(Tuple{Real, Real}, Nothing)=nothing`: calculate power relative to frequency range or total power

# Returns

Named tuple containing:
- `pw::Vector{Float64}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd_rel(s::AbstractVector; fs::Int64, norm::Bool=false, mt::Bool=false, nt::Int64=8, f::Union{Tuple{Real, Real}, Nothing}=nothing)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    if f !== nothing
        f = tuple_order(f)
        f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0.")) 
        f[2] > fs / 2 && throw(ArgumentError("Lower frequency bound must be ≤ $(fs / 2).")) 
    end

    # for short signals use multi-tapered periodogram
    length(s) < 4 * fs && (mt = true)

    ref_pw = f === nothing ? total_power(s, fs=fs, mt=mt) : band_power(s, fs=fs, mt=mt, nt=nt, f=f)

    if mt == true
        p = mt_pgram(s, fs=fs, nw=(nt÷2+1), ntapers=nt)
    else
        p = welch_pgram(s, 4*fs, fs=fs)
    end

    pw = power(p)
    pf = Vector(freq(p))
    pw = pw[1:length(pf)]

    # replace power at 0 Hz
    pw[1] = pw[2]

    pw = pw / ref_pw

    norm == true && (pw = pow2db.(pw))

    return (pw=pw, pf=pf)

end

"""
    psd_rel(s; fs, norm, mt, nt, f)

Calculate relative power spectrum density.

# Arguments
- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers
- `f::Union(Tuple{Real, Real}, Nothing)=nothing`: calculate power relative to frequency range or total power

# Returns

Named tuple containing:
- `pw::Array{Float64, 3}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd_rel(s::AbstractArray; fs::Int64, norm::Bool=false, mt::Bool=false, nt::Int64=8, f::Union{Tuple{Real, Real}, Nothing}=nothing)

    # for short signals use multi-tapered periodogram
    if size(s, 2) < 4 * fs
        mt = true
        _info("Using multi-tapered periodogram.")
    end

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, pf = psd_rel(s[1, :, 1], fs=fs, norm=norm, mt=mt, nt=nt, f=f)

    pw = zeros(ch_n, length(pf), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pw[ch_idx, :, ep_idx], _ = psd_rel(s[ch_idx, :, ep_idx], fs=fs, norm=norm, mt=mt, nt=nt, f=f)
        end
    end
    
    return (pw=pw, pf=pf)

end

"""
    psd_rel(obj; ch, norm, mt, nt, f)

Calculate relative power spectrum density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers
- `f::Union(Tuple{Real, Real}, Nothing)=nothing`: calculate power relative to frequency range or total power

# Returns

Named tuple containing:
- `pw::Array{Float64, 3}`: powers
- `pf::Array{Float64, 3}`: frequencies
"""
function psd_rel(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), norm::Bool=false, mt::Bool=false, nt::Int64=8, f::Union{Tuple{Real, Real}, Nothing}=nothing)

    _check_channels(obj, ch)

    pw, pf = @views psd_rel(obj.data[ch, :, :], fs=sr(obj), norm=norm, mt=mt, nt=nt, f=f)

    return (pw=pw, pf=pf)

end
