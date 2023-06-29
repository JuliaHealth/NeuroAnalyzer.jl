export total_power

"""
    total_power(s; fs, mt)

Calculate total power.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

- `tp::Float64`: total power
"""
function total_power(s::AbstractVector; fs::Int64, mt::Bool=false, nt::Int64=8)

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

    # dx: frequency resolution
    dx = pf[2] - pf[1]
    tp = simpson(pw, dx=dx)

    return tp

end

"""
    total_power(s; fs, mt)

Calculate total power.

`# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

- `tp::Matrix{Float64}`: total power

"""
function total_power(s::AbstractArray; fs::Int64, mt::Bool=false, nt::Int64=8)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    tp = zeros(ch_n, ep_n)
    
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            @views tp[ch_idx, ep_idx] = total_power(s[ch_idx, :, ep_idx], fs=fs, mt=mt, nt=nt)
        end
    end

    return tp
end

"""
    total_power(obj, ch, mt)

Calculate total power.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(record)`: index of channels, default is all signal channels
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns
 
- `tp::Matrix{Float64}`: total power
"""
function total_power(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), mt::Bool=false, nt::Int64=8)

    _check_channels(obj, ch)

    tp = @views total_power(obj.data[ch, :, :], fs=sr(obj), mt=mt, nt=nt)

    return tp

end
