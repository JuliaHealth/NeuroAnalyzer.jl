export fir_order_bw
export fir_order_f
export iir_order

"""
    fir_order_bw(; <keyword arguments>)

Calculate the order of a FIR filter using the Harris formula.

The formula is: `n = round((a × fs) / (22 × bw))`, where 22 is the Harris empirical constant relating attenuation, bandwidth, and sample rate.

# Arguments

- `bw::Real`: transition band width in Hz; must be > 0.
- `a::Real=60`: attenuation in dB; must be > 0
- `fs::Int64`: sampling rate in Hz; must be > 0

# Returns

- `Int64`: estimated FIR filter order

# Throws

- `ArgumentError`: if `bw ≤ 0`, `a ≤ 0`, or `fs ≤ 0`.

# See also

[`fir_order_f`](@ref), [`iir_order`](@ref)
"""
function fir_order_bw(; bw::Real, a::Real = 60, fs::Int64)::Int64

    !(bw > 0) && throw(ArgumentError("bw must be > 0."))
    !(a > 0) && throw(ArgumentError("a must be > 0."))
    !(fs > 0) && throw(ArgumentError("fs must be > 0."))

    return round(Int64, (a * fs) / (22 * bw))

end

"""
    fir_order_bw(obj; <keyword arguments>)

Calculate the order of a FIR filter using the Harris formula.

Convenience wrapper that reads the sampling rate from `obj`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `bw::Real`: transition band width in Hz; must be > 0
- `a::Real=60`: attenuation in dB; must be > 0

# Returns

- `Int64`: estimated FIR filter order

# Throws

- `ArgumentError`: if `bw ≤ 0` or `a ≤ 0`

# See also

[`fir_order_f`](@ref), [`iir_order`](@ref)
"""
function fir_order_bw(obj::NeuroAnalyzer.NEURO; bw::Real, a::Real = 60)::Int64

    return fir_order_bw(bw = bw, a = a, fs = sr(obj))

end

"""
    fir_order_f(; <keyword arguments>)

Calculate a recommended FIR filter order range from the lower frequency bound.

The rule of thumb is that the filter should span 4–5 full cycles of the lowest frequency of interest. The range is returned as `(4 × n_samples, 5 × n_samples)`, where `n_samples` is the number of samples in one cycle of `f`.

# Arguments

- `fs::Int64`: sampling rate in Hz; must be > 0
- `f::Real`: lower frequency bound of the analyzed range in Hz; must be > 0

# Returns

- `Tuple{Int64, Int64}`: `(lower_order, upper_order)` recommended filter order range

# Throws

- `ArgumentError`: if `fs ≤ 0` or `f ≤ 0`.

# See also

[`fir_order_bw`](@ref)
"""
function fir_order_f(; fs::Int64, f::Real)::Tuple{Int64, Int64}

    !(fs > 0) && throw(ArgumentError("fs must be > 0."))
    !(f  > 0) && throw(ArgumentError("f must be > 0."))

    # samples per one full cycle of the lowest frequency of interest
    cycle_samples = t2s(1 / f, fs)

    # convert to samples
    return (4 * cycle_samples, 5 * cycle_samples)

end

"""
    fir_order_f(obj; <keyword arguments>)

Calculate a recommended FIR filter order range from the lower frequency bound.

Convenience wrapper that reads the sampling rate from `obj`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `f::Real`: lower frequency bound of the analysed range in Hz; must be > 0

# Returns

- `Tuple{Int64, Int64}`: `(lower_order, upper_order)` recommended filter order range
"""
function fir_order_f(obj::NeuroAnalyzer.NEURO; f::Real)::Tuple{Int64, Int64}

    return fir_order_f(; f = f, fs = sr(obj))

end

"""
    iir_order(; <keyword arguments>)

Calculate the minimum order of an IIR filter meeting the given specifications.

Pass-band and stop-band edges are derived from `cutoff` and `bw`:

- For `:lp` / `:hp`: pass edge = `cutoff ± bw/2`, stop edge = `cutoff ∓ bw/2`.
- For `:bp` / `:bs`: inner edges shifted inward by `bw/2`, outer edges shifted outward by `bw/2`.

The order is estimated via the appropriate `DSP.jl` design function (`buttord`, `cheb1ord`, `cheb2ord`, or `ellipord`).

# Arguments

- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
- `ftype::Symbol`: filter type:
    - `:lp`: low-pass
    - `:hp`: high-pass
    - `:bp`: band-pass
    - `:bs`: band-stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: cutoff frequency in Hz; scalar for `:lp` / `:hp`; two-element tuple for `:bp` / `:bs`
- `bw::Real`: transition band width in Hz; must be > 0
- `rp::Union{Nothing, Real}=nothing`: pass-band ripple in dB; defaults to 0.5 dB
- `rs::Union{Nothing, Real}=nothing`: stop-band ripple in dB; defaults to 20 dB
- `fs::Int64`: sampling rate in Hz; must be > 0

# Returns

- `order::Int64`: minimum filter order satisfying the specifications

# Throws

- `ArgumentError`: if `fprototype` or `ftype` is invalid, `bw ≤ 0`, `fs ≤ 0`, or `cutoff` has wrong length for the chosen `ftype`.

# See also

[`fir_order_bw`](@ref), [`fir_order_f`](@ref)
"""
function iir_order(;
    fprototype::Symbol,
    ftype::Symbol,
    cutoff::Union{Real, Tuple{Real, Real}},
    bw::Real,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    fs::Int64
)::Int64

    _check_var(
        fprototype, [:butterworth, :chebyshev1, :chebyshev2, :elliptic], "fprototype"
    )
    _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")
    !(bw > 0) && throw(ArgumentError("bw must be > 0."))
    !(fs > 0) && throw(ArgumentError("fs must be > 0."))

    # nyquist frequency; used to normalise cutoff edges to [0, 1]
    nqf = fs / 2

    # apply prototype-specific ripple defaults if not supplied by the caller
    isnothing(rp) && (rp = 0.5)
    isnothing(rs) && (rs = 20.0)

    # compute normalized pass-band (wp) and stop-band (ws) edges
    # the transition band is centered on `cutoff`; half-width = bw/2.
    if ftype === :lp
        !(length(cutoff) == 1) && throw(ArgumentError("cutoff must specify exactly one frequency for :lp."))
        wp = (cutoff[1] - bw / 2) / nqf  # pass edge (below cutoff)
        ws = (cutoff[1] + bw / 2) / nqf  # stop edge (above cutoff)
    elseif ftype === :hp
        !(length(cutoff) == 1) && throw(ArgumentError("cutoff must specify exactly one frequency for :hp."))
        ws = (cutoff[1] - bw / 2) / nqf  # stop edge (below cutoff)
        wp = (cutoff[1] + bw / 2) / nqf  # pass edge (above cutoff)
    elseif ftype === :bp
        !(length(cutoff) == 2) && throw(ArgumentError("cutoff must specify exactly two frequencies for :bp."))
        wp = ((cutoff[1] + bw / 2) / nqf, (cutoff[2] - bw / 2) / nqf)  # inner pass edges
        ws = ((cutoff[1] - bw / 2) / nqf, (cutoff[2] + bw / 2) / nqf)  # outer stop edges
    elseif ftype === :bs
        !(length(cutoff) == 2) && throw(ArgumentError("cutoff must specify exactly two frequencies for :bs."))
        ws = ((cutoff[1] + bw / 2) / nqf, (cutoff[2] - bw / 2) / nqf)  # inner stop edges
        wp = ((cutoff[1] - bw / 2) / nqf, (cutoff[2] + bw / 2) / nqf)  # outer pass edges
    end

    # delegate order estimation to the appropriate DSP.jl function
    if fprototype === :butterworth
        order = buttord(wp, ws, rp, rs)[1]
    elseif fprototype === :chebyshev1
        order = cheb1ord(wp, ws, rp, rs)[1]
    elseif fprototype === :chebyshev2
        order = cheb2ord(wp, ws, rp, rs)[1]
    elseif fprototype === :elliptic
        order = ellipord(wp, ws, rp, rs)[1]
    end

    return order

end

"""
    iir_order(obj; <keyword arguments>)

Calculate the minimum order of an IIR filter meeting the given specifications.

Convenience wrapper that reads the sampling rate from `obj`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
- `ftype::Symbol`: filter type:
    - `:lp`: low-pass
    - `:hp`: high-pass
    - `:bp`: band-pass
    - `:bs`: band-stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: cutoff frequency in Hz; scalar for `:lp` / `:hp`; two-element tuple for `:bp` / `:bs`
- `bw::Real`: transition band width in Hz; must be > 0
- `rp::Union{Nothing, Real}=nothing`: pass-band ripple in dB; defaults to 0.5 dB
- `rs::Union{Nothing, Real}=nothing`: stop-band ripple in dB; defaults to 20 dB

# Returns

- `Int64`: minimum filter order satisfying the specifications

# See also

[`fir_order_bw`](@ref), [`fir_order_f`](@ref)
"""
function iir_order(
    obj::NeuroAnalyzer.NEURO;
    fprototype::Symbol,
    ftype::Symbol,
    cutoff::Union{Real, Tuple{Real, Real}},
    bw::Real,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing
)::Int64

    return iir_order(
        fprototype = fprototype,
        ftype = ftype,
        cutoff = cutoff,
        bw = bw,
        rp = rp,
        rs = rs,
        fs = sr(obj),
    )

end
