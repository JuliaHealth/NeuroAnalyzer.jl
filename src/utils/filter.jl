export fir_order_bw
export fir_order_f
export iir_order

"""
    fir_order_bw(; <keyword arguments>)

Calculate order of FIR filter using Harris formula.

# Arguments

  - `bw::Real`: transition band width in Hz
  - `a::Real`: attenuation in dB
  - `fs::Int64`: sampling rate

# Returns

  - `n::Int64`
"""
function fir_order_bw(; bw::Real, a::Real = 60, fs::Int64)::Int64

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
function fir_order_bw(obj::NeuroAnalyzer.NEURO; bw::Real, a::Real = 60)::Int64

    n = fir_order_bw(; bw = bw, a = a, fs = sr(obj))

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
    fir_order_f(obj; <keyword arguments>)

Calculate order of FIR filter using lower frequency bound.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `f::Real`: lower frequency bound for analyzed range

# Returns

  - `n::Tuple{Int64, Int64}`: recommended order range
"""
function fir_order_f(obj::NeuroAnalyzer.NEURO; f::Real)::Tuple{Int64, Int64}

    n = fir_order_f(; f = f, fs = sr(obj))

    return n

end

"""
    iir_order(; <keyword arguments>)

Calculate order of IIR filter.

# Arguments

  - `fprototype::Symbol`: filter prototype:
      + `:butterworth`
      + `:chebyshev1`
      + `:chebyshev2`
      + `:elliptic`
  - `ftype::Symbol`: filter type:
      + `:lp`: low pass
      + `:hp`: high pass
      + `:bp`: band pass
      + `:bs`: band stop
  - `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (a pair of frequencies for band pass and band stop filters)
  - `bw::Real`: transition band width
  - `rp::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  - `rs::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  - `fs::Int64`: sampling rate

# Returns

  - `n::Int64`
"""
function iir_order(;
    fprototype::Symbol,
    ftype::Symbol,
    cutoff::Union{Real, Tuple{Real, Real}},
    bw::Real,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    fs::Int64,
)::Int64

    _check_var(fprototype, [:butterworth, :chebyshev1, :chebyshev2, :elliptic], "fprototype")
    _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")

    nqf = fs / 2

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
        wp = (cutoff[1] - (bw / 2)) / nqf
        ws = (cutoff[1] + (bw / 2)) / nqf
    elseif ftype === :hp
        @assert length(cutoff) == 1 "cutoff must specify only one frequency."
        ws = (cutoff[1] - (bw / 2)) / nqf
        wp = (cutoff[1] + (bw / 2)) / nqf
    elseif ftype === :bp
        @assert length(cutoff) == 2 "cutoff must specify two frequencies."
        wp = ((cutoff[1] + (bw / 2)) / nqf, (cutoff[2] - (bw / 2)) / nqf)
        ws = ((cutoff[1] - (bw / 2)) / nqf, (cutoff[2] + (bw / 2)) / nqf)
    elseif ftype === :bs
        @assert length(cutoff) == 2 "cutoff must specify two frequencies."
        ws = ((cutoff[1] + (bw / 2)) / nqf, (cutoff[2] - (bw / 2)) / nqf)
        wp = ((cutoff[1] - (bw / 2)) / nqf, (cutoff[2] + (bw / 2)) / nqf)
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

Calculate order of IIR filter.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `fprototype::Symbol`: filter prototype:
      + `:butterworth`
      + `:chebyshev1`
      + `:chebyshev2`
      + `:elliptic`
  - `ftype::Symbol`: filter type:
      + `:lp`: low pass
      + `:hp`: high pass
      + `:bp`: band pass
      + `:bs`: band stop
  - `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (a pair of frequencies for band pass and band stop filters)
  - `bw::Real`: transition band width
  - `rp::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  - `rs::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others

# Returns

  - `n::Int64`
"""
function iir_order(
    obj::NeuroAnalyzer.NEURO;
    fprototype::Symbol,
    ftype::Symbol,
    cutoff::Union{Real, Tuple{Real, Real}},
    bw::Real,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
)::Int64

    n = iir_order(;
                fprototype = fprototype,
                ftype = ftype,
                cutoff = cutoff,
                bw = bw,
                rp = rp,
                rs = rs,
                fs = sr(obj),
            )

    return n

end
