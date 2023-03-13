export total_power

"""
    total_power(signal; fs, mt)

Calculate total power.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`: sampling rate
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

- `tp::Float64`: signal total power
"""
function total_power(signal::AbstractVector; fs::Int64, mt::Bool=false, nt::Int64=8)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

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

    # dx: frequency resolution
    dx = pf[2] - pf[1]
    tp = simpson(psd_pow, dx=dx)

    return tp

end

"""
    total_power(signal; fs, mt)

Calculate total power.

`# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

- `tp::Matrix{Float64}`: total power

"""
function total_power(signal::AbstractArray; fs::Int64, mt::Bool=false, nt::Int64=8)

    ch_n = size(signal, 1)
    ep_n = size(signal, 3)

    tp = zeros(ch_n, ep_n)
    
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            @views tp[ch_idx, ep_idx] = total_power(signal[ch_idx, :, ep_idx], fs=fs, mt=mt, nt=nt)
        end
    end

    return tp
end

"""
    total_power(obj, channel, mt)

Calculate total power.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(record)`: index of channels, default is all signal channels
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns
 
- `tp::Matrix{Float64}`: total power
"""
function total_power(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), mt::Bool=false, nt::Int64=8)

    _check_channels(obj, channel)

    tp = @views total_power(obj.data[channel, :, :], fs=sr(obj), mt=mt, nt=nt)

    return tp

end
