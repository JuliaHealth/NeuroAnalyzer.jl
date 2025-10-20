export linspace
export log10space
export log2space
export cmax
export cmin
export tuple_order
export cums
export f_nearest
export ntapers
export trtm
export fir_order_bw
export fir_order_f
export iir_order

"""
    linspace(start, stop, length)

Generates a sequence of evenly spaced numbers between `start` and `stop`.

# Arguments

- `start::Number`
- `stop::Number`
- `n::Int64`: sequence length

# Returns

- `range::Vector{Float64}`
"""
function linspace(start::Number, stop::Number, n::Int64)::Vector{Float64}

    @assert n >= 2 "n must be ≥ 2."

    return collect(range(start, stop, n))

end

"""
    log10space(start, stop, n)

Generates a sequence of log10-spaced numbers between `start` and `stop`.

# Arguments

- `start::Number`
- `stop::Number`
- `n::Int64`: sequence length

# Returns

- `range::Vector{Float64}`
"""
function log10space(start::Number, stop::Number, n::Int64)::Vector{Float64}

    @assert n >= 2 "n must be ≥ 2."

    return collect(exp10.(range(start, stop, n)))

end

"""
    log2space(start, stop, n)

Generates a sequence of log2-spaced numbers between `start` and `stop`.

# Arguments

- `start::Number`
- `stop::Number`
- `n::Int64`: sequence length

# Returns

- `range::Vector{Float64}`
"""
function log2space(start::Number, stop::Number, n::Int64)::Vector{Float64}

    @assert n >= 2 "n must be ≥ 2."

    return collect(exp2.(range(start, stop, n)))

end

"""
    cmax(x)

Return maximum value of the complex vector.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmax::ComplexF64`
"""
function cmax(x::Vector{<:Complex})::ComplexF64

    return argmax(abs2, x)

end

"""
    cmin(x)

Return minimum value of the complex vector.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmin::ComplexF64`
"""
function cmin(x::Vector{<:Complex})::ComplexF64

    return argmin(abs2, x)

end

"""
    cextrema(x)

Return extreme values of the complex vector.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

Tuple containing:
- `cmax::ComplexF64`
- `cmin::ComplexF64`
"""
function cextrema(x::Vector{<:Complex})::Tuple{ComplexF64, ComplexF64}

    return (cmax(x), cmin(x))

end

"""
    tuple_order(t, rev)

Order tuple elements in ascending or descending (`rev=true`) order.

# Arguments

- `t::Tuple{Real, Real}`
- `rev::Bool=false`

# Returns

- `t::Tuple{Real, Real}`
"""
function tuple_order(t::Tuple{Real, Real}, rev::Bool=false)::Tuple{Real, Real}

    (!rev && t[1] > t[2]) && (t = (t[2], t[1]))
    (rev && t[1] < t[2]) && (t = (t[2], t[1]))

    return t

end

"""
    cums(signal)

Calculate cumulative sum of a 3-dimensional array.

# Arguments

- `signal::Array{<:Real, 3}`

# Returns

- `signal_cs::Array{Float64, 3}`
"""
function cums(signal::Array{<:Real, 3})::Array{Float64, 3}

    ch_n, _, ep_n = size(signal)
    signal_cs = similar(signal)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            signal_cs[ch_idx, :, ep_idx] = @views cumsum(signal[ch_idx, :, ep_idx])
        end
    end

    return signal_cs

end

"""
    f_nearest(m, pos)

Find nearest position tuple in a matrix of positions.

# Arguments

- `m::Matrix{Tuple{Float64, Float64}}`: matrix of positions
- `p::Tuple{Float64, Float64}`: position tuple

# Returns

- `pos::Tuple{Int64, Int64}`: row and column in m
"""
function f_nearest(m::Matrix{Tuple{Float64, Float64}}, p::Tuple{Float64, Float64})::Tuple{Int64, Int64}

    d = zeros(size(m))

    @inbounds for idx1 in axes(m, 1)
        for idx2 in axes(m, 2)
            d[idx1, idx2] = euclidean(m[idx1, idx2], p)
        end
    end

    return (findmin(d)[2][1], findmin(d)[2][2])

end

"""
    ntapers(obj; df)

Return recommended number of Slepian tapers for multi-taper power spectrum analysis.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `df::Real`: frequency resolution (bandwidth); smallest distance between frequency peaks that we want to observe (e.g. 1 Hz)

# Returns

- `nt::Int64`
"""
function ntapers(obj::NeuroAnalyzer.NEURO; df::Real)::Int64

    _bin(df, (0, sr(obj) / 2))

    n = epoch_len(obj) / sr(obj)
    nt = round(Int64, df * n) - 1

    return nt

end

"""
    trtm(obj; <keyword arguments>)

Return signal channel in the form trials × time.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `ep::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj))`: epoch numbers; default use all epochs

# Returns

- `s::Matrix{Float64}`
"""
function trtm(obj::NeuroAnalyzer.NEURO; ch::String, ep::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj)))::Matrix{Float64}

    _check_epochs(obj, ep)
    ch = get_channel(obj, ch=ch)

    s = zeros(length(ep), epoch_len(obj))
    @inbounds for ep_idx in eachindex(ep)
        s[ep_idx, :] = @views obj.data[ch, :, ep_idx]
    end

    return s

end

"""
    fir_order_bw(; <keyword arguments>)

Calculate order of FIR filter using Harris formula.

# Arguments

- `bw::Real`: transition band width
- `a::Real`: attenuation in dB
- `fs::Int64`: sampling rate

# Returns

- `n::Int64`
"""
function fir_order_bw(; bw::Real, a::Real=60, fs::Int64)::Int64

    n = round(Int64, (a * fs)/(22 * bw))

    return n

end

"""
    fir_order_bw(obj; <keyword arguments>)

Calculate order of FIR filter using Harris formula.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `bw::Real`: transition band width
- `a::Real`: attenuation in dB

# Returns

- `n::Int64`
"""
function fir_order_bw(obj::NeuroAnalyzer.NEURO; bw::Real, a::Real=60)::Int64

    n = fir_order_bw(bw=bw, a=a, fs=sr(obj))

    return n

end

"""
    fir_order_f(obj; <keyword arguments>)

Calculate order of FIR filter using lower frequency bound.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `f::Real`: lower frequency bound for analyzed range

# Returns

- `n::Tuple{Int64, Int64}`: recommended order range
"""
function fir_order_f(obj::NeuroAnalyzer.NEURO; f::Real)::Tuple{Int64, Int64}

    n = fir_order_f(f=f, fs=sr(obj))

    return n

end

"""
    fir_order_f(; <keyword arguments>)

Calculate order of FIR filter using lower frequency bound.

# Arguments

- `f::Real`: lower frequency bound for analyzed range
- `fs::Int64`: sampling rate

# Returns

- `n::Tuple{Int64, Int64}`: recommended order range
"""
function fir_order_f(; fs::Int64, f::Real)::Tuple{Int64, Int64}

    # single cycle length in s
    sc = (1 / f)
    
    # convert to samples
    n = (4 * t2s(sc, fs), 5 * t2s(sc, fs))

    return n

end

"""
    iir_order(; <keyword arguments>)

Calculate order of IIR filter.

# Arguments

- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (a pair of frequencies for band pass and band stop filters)
- `bw::Real`: transition band width
- `rp::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `fs::Int64`: sampling rate

# Returns

- `n::Int64`
"""
function iir_order(; fprototype::Symbol, ftype::Symbol, cutoff::Union{Real, Tuple{Real, Real}}, bw::Real, rp::Union{Nothing, Real}=nothing, rs::Union{Nothing, Real}=nothing, fs::Int64)::Int64

    _check_var(fprototype, [:butterworth, :chebyshev1, :chebyshev2, :elliptic], "fprototype")
    _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")

    if rp isa Nothing
        if fprototype === :elliptic
            rp = 0.0025
        else
            rp = 2
        end
    end

    if rs isa Nothing
        if fprototype === :elliptic
            rs = 40
        else
            rs = 20
        end
    end

    if ftype === :lp
        @assert length(cutoff) == 1 "cutoff must specify only one frequency."
        wp = (cutoff[1] - (bw /2)) / (fs / 2)
        ws = (cutoff[1] + (bw /2)) / (fs / 2)
    elseif ftype === :hp
        @assert length(cutoff) == 1 "cutoff must specify only one frequency."
        ws = (cutoff[1] - (bw /2)) / (fs / 2)
        wp = (cutoff[1] + (bw /2)) / (fs / 2)
    elseif ftype === :bp
        @assert length(cutoff) == 2 "cutoff must specify two frequencies."
        wp = ((cutoff[1] + (bw /2)) / (fs / 2), (cutoff[2] - (bw /2)) / (fs / 2))
        ws = ((cutoff[1] - (bw /2)) / (fs / 2), (cutoff[2] + (bw /2)) / (fs / 2))
    elseif ftype === :bs
        @assert length(cutoff) == 2 "cutoff must specify two frequencies."
        ws = ((cutoff[1] + (bw /2)) / (fs / 2), (cutoff[2] - (bw /2)) / (fs / 2))
        wp = ((cutoff[1] - (bw /2)) / (fs / 2), (cutoff[2] + (bw /2)) / (fs / 2))
    end

    if fprototype === :butterworth
        n, _ = buttord(wp, ws, rp, rs)
    elseif fprototype === :chebyshev1
        n, _ = cheb1ord(wp, ws, rp, rs)
    elseif fprototype === :chebyshev2
        n, _ = cheb2ord(wp, ws, rp, rs)
    elseif fprototype === :elliptic
        n, _ = ellipord(wp, ws, rp, rs)
    end

    return n

end

"""
    iir_order(obj; <keyword arguments>)

Calculate order of FIR filter using harris' formula.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (a pair of frequencies for band pass and band stop filters)
- `bw::Real`: transition band width
- `rp::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others

# Returns

- `n::Int64`
"""
function iir_order(obj::NeuroAnalyzer.NEURO; fprototype::Symbol, ftype::Symbol, cutoff::Union{Real, Tuple{Real, Real}}, bw::Real, rp::Union{Nothing, Real}=nothing, rs::Union{Nothing, Real}=nothing)::Int64

    n = iir_order(fprototype=fprototype, ftype=ftype, cutoff=cutoff, bw=bw, rp=rp, rs=rs, fs=sr(obj))

    return n

end