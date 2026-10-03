export filter_create
export filter_apply
export filter_apply!
export filter
export filter!
export filter_order
export filter_report

# window-method FIR: transition width (pass edge -> stop edge) ≈ k · fs / N,
# and minimum stop-band attenuation (dB) of the window
const _FIR_WINDOWS = Dict(
    :rect => (k = 0.9, att = 21.0),
    :hann => (k = 3.1, att = 44.0),
    :hamming => (k = 3.3, att = 53.0),
    :blackman => (k = 5.5, att = 74.0),
)

const _FLT_TYPES = Union{
    Vector{Float64},
    ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64},
    Biquad{:z, Float64},
}

function _fir_window(window::Symbol, n::Int64)::Vector{Float64}
    window === :hamming && return DSP.hamming(n)
    window === :hann && return DSP.hanning(n)
    window === :blackman && return DSP.blackman(n)
    window === :rect && return DSP.rect(n)
    throw(ArgumentError("Unknown window: $window"))
end

# transition band(s) centred on the cutoff(s)
function _transition_edges(ftype::Symbol, cutoff, bw::Real)
    if ftype in (:lp, :hp)
        return [(cutoff[1] - bw / 2, cutoff[1] + bw / 2)]
    else
        return [(cutoff[1] - bw / 2, cutoff[1] + bw / 2), (cutoff[2] - bw / 2, cutoff[2] + bw / 2)]
    end
end

# transition bands must lie within (0, Nyquist) and must not overlap
function _check_edges(ftype::Symbol, cutoff, bw::Real, nqf::Real; strict::Bool)::Nothing
    tr = _transition_edges(ftype, cutoff, bw)
    msg = String[]
    tr[1][1] <= 0 && push!(msg, "transition band ($(tr[1][1])–$(tr[1][2]) Hz) reaches 0 Hz")
    tr[end][2] >= nqf && push!(msg, "transition band ($(tr[end][1])–$(tr[end][2]) Hz) reaches Nyquist ($nqf Hz)")
    length(tr) == 2 && tr[1][2] >= tr[2][1] && push!(msg, "transition bands overlap")
    isempty(msg) && return nothing
    m = Base.join(msg, "; ") * " (cutoff=$cutoff Hz, bw=$bw Hz)"
    if strict
        throw(ArgumentError(m * ". Reduce bw."))
    else
        _warn(m * ": the nominal stop-band attenuation will not be reached; reduce bw (increase order). Check filter_report().")
    end
    return nothing
end

"""
    filter_order(; <keyword arguments>)

Estimate the FIR filter order (number of taps) for a given transition band width.

# Arguments

- `fprototype::Symbol`: `:fir`, `:firls` or `:remez`
- `fs::Int64`: sampling rate in Hz
- `bw::Real`: transition band width in Hz (pass-band edge to stop-band edge)
- `window::Symbol=:hamming`: window for `:fir` (`:rect`, `:hann`, `:hamming`, `:blackman`)
- `rs::Union{Nothing, Real}=nothing`: target stop-band attenuation in dB for `:firls`/`:remez` (default 53 dB, as Hamming)

# Returns

- `Int64`: number of taps, rounded up to an odd number (type-I linear-phase FIR, valid for all filter types)

# Notes

- `:fir`: N = ⌈k · fs / bw⌉, k = 0.9 (rect), 3.1 (hann), 3.3 (hamming), 5.5 (blackman)
- `:firls`, `:remez`: harris' rule N = ⌈(A / 22) · fs / bw⌉, A = stop-band attenuation in dB
"""
function filter_order(;
    fprototype::Symbol,
    fs::Int64,
    bw::Real,
    window::Symbol = :hamming,
    rs::Union{Nothing, Real} = nothing,
)::Int64
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    bw > 0 || throw(ArgumentError("bw must be > 0."))
    n = if fprototype === :fir
        _check_var(window, collect(keys(_FIR_WINDOWS)), "window")
        ceil(Int64, _FIR_WINDOWS[window].k * fs / bw)
    elseif fprototype in (:firls, :remez)
        a = isnothing(rs) ? 53.0 : rs
        ceil(Int64, a / 22 * fs / bw)
    else
        throw(ArgumentError("Automatic order is available only for :fir, :firls and :remez."))
    end
    return isodd(n) ? n : n + 1
end

"""
    filter_create(; <keyword arguments>)

Create a FIR or IIR filter object.

# Arguments

- `fprototype::Symbol`: filter prototype:
    - `:fir`: FIR filter (window method)
    - `:firls`: weighted least-squares FIR filter
    - `:remez`: Remez (Parks-McClellan) FIR filter
    - `:butterworth`: Butterworth IIR filter
    - `:chebyshev1`: Chebyshev type-I IIR filter
    - `:chebyshev2`: Chebyshev type-II IIR filter
    - `:elliptic`: elliptic IIR filter
    - `:iirnotch`: second-order IIR notch filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz
    - for `:lp`/`:hp`: single frequency
    - for `:bp`/`:bs`: frequency range (f1, f2)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `order::Union{Nothing, Int64}=nothing`: filter order (number of taps for FIR, filter order for IIR); for FIR prototypes, if `nothing` it is calculated from `bw` (see `filter_order`); required for IIR prototypes
- `rp::Union{Nothing, Real}=nothing`: pass-band ripple in dB (default 0.5 dB)
- `rs::Union{Nothing, Real}=nothing`: stop-band attenuation in dB (default 20 dB for IIR; target for automatic `:firls`/`:remez` order, default 53 dB)
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz, centred on the cutoff (required for `:firls`, `:remez`, `:iirnotch`; for `:fir` required unless `order` or `w` is given)
- `w::Union{Nothing, AbstractVector}=nothing`: window vector for `:fir` (overrides `window` and `order`) or weight vector for `:firls`
- `window::Symbol=:hamming`: window for `:fir` when `w` is not given (`:rect`, `:hann`, `:hamming`, `:blackman`)

# Returns

- `Vector{Float64}`: FIR filter coefficients (for `:fir`, `:firls`, `:remez`)
- `ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}`: IIR filter in zero-pole-gain form (for `:butterworth`, `:chebyshev1`, `:chebyshev2`, `:elliptic`)
- `Biquad{:z, Float64}`: second-order biquad filter (for `:iirnotch`)
"""
function filter_create(;
    fprototype::Symbol,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Real, Tuple{Real, Real}},
    fs::Int64,
    order::Union{Nothing, Int64} = nothing,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    bw::Union{Nothing, Real} = nothing,
    w::Union{Nothing, AbstractVector} = nothing,
    window::Symbol = :hamming,
)::_FLT_TYPES
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    nqf = fs / 2

    _check_var(
        fprototype,
        [:fir, :firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch],
        "fprototype",
    )
    !isnothing(ftype) && _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")
    _check_var(window, collect(keys(_FIR_WINDOWS)), "window")
    !isnothing(order) && order < 1 && throw(ArgumentError("order must be ≥ 1."))
    !isnothing(bw) && bw <= 0 && throw(ArgumentError("bw must be > 0."))

    # --- ftype required (all except :iirnotch) ---
    if fprototype !== :iirnotch
        isnothing(ftype) && throw(ArgumentError("ftype must be specified for $fprototype."))
    end

    # --- cutoff arity, values and normalization ---
    if fprototype === :iirnotch || ftype in (:lp, :hp)
        length(cutoff) == 1 || throw(ArgumentError("For :$(something(ftype, fprototype)), cutoff must be a scalar."))
        cutoff = cutoff[1]
        cutoff > 0 || throw(ArgumentError("cutoff must be > 0 Hz."))
        cutoff < nqf || throw(ArgumentError("cutoff must be < $nqf Hz (Nyquist)."))
    else
        length(cutoff) == 2 || throw(ArgumentError("For :$ftype, cutoff must specify two frequencies."))
        if cutoff[1] > cutoff[2]
            cutoff = (cutoff[2], cutoff[1])
            _warn("cutoff frequencies swapped to $cutoff Hz.")
        end
        cutoff[1] == cutoff[2] && throw(ArgumentError("cutoff frequencies must differ."))
        cutoff[1] > 0 || throw(ArgumentError("cutoff must be > 0 Hz."))
        cutoff[2] < nqf || throw(ArgumentError("cutoff must be < $nqf Hz (Nyquist)."))
    end

    # --- :fir: order from bw (automatic) or manual; window ---
    if fprototype === :fir
        custom_w = !isnothing(w)
        if custom_w
            length(w) >= 1 || throw(ArgumentError("Length of w must be ≥ 1."))
            !isnothing(order) && order != length(w) &&
                throw(ArgumentError("Length of w ($(length(w))) must equal order ($order)."))
            order = length(w)
            _info("Custom window: order = length(w) = $order taps; transition width not estimated, check filter_report()")
        else
            if isnothing(order)
                isnothing(bw) && throw(ArgumentError("bw or order must be specified for :fir."))
                order = filter_order(; fprototype = :fir, fs = fs, bw = bw, window = window)
                _info("order calculated from bw=$bw Hz ($window window): $order taps")
            else
                bw_eff = _FIR_WINDOWS[window].k * fs / order
                if !isnothing(bw) && abs(bw_eff - bw) / bw > 0.1
                    _warn(
                        "order=$order gives a transition width of ≈$(round(bw_eff; digits = 2)) Hz " *
                        "($window window), not bw=$bw Hz; bw=$bw Hz needs ≈$(filter_order(; fprototype = :fir, fs = fs, bw = bw, window = window)) taps.",
                    )
                end
            end
            w = _fir_window(window, order)
        end
        ftype in (:hp, :bp, :bs) && iseven(order) &&
            throw(ArgumentError("order must be odd for :hp/:bp/:bs filters."))
        # edge check with the transition width the design actually has (standard windows only)
        custom_w || _check_edges(ftype, cutoff, _FIR_WINDOWS[window].k * fs / order, nqf; strict = false)
    end

    # --- :firls / :remez: bw required, order from bw (automatic) or manual ---
    if fprototype in (:firls, :remez)
        isnothing(bw) && throw(ArgumentError("bw must be specified for $fprototype."))
        _check_edges(ftype, cutoff, bw, nqf; strict = true)
        if isnothing(order)
            order = filter_order(; fprototype = fprototype, fs = fs, bw = bw, rs = rs)
            _info("order calculated from bw=$bw Hz (target attenuation $(isnothing(rs) ? 53 : rs) dB): $order taps")
        end
        ftype in (:hp, :bs) && iseven(order) &&
            throw(ArgumentError("order must be odd for :hp/:bs filters."))
    end

    # --- :iirnotch ---
    if fprototype === :iirnotch
        isnothing(bw) && throw(ArgumentError("bw must be specified for :iirnotch."))
        !isnothing(ftype) && _info("For :iirnotch filter ftype is ignored")
        !isnothing(order) && _info("For :iirnotch filter order is ignored")
        (cutoff - bw / 2 > 0 && cutoff + bw / 2 < nqf) ||
            throw(ArgumentError("Notch band ($(cutoff - bw / 2)–$(cutoff + bw / 2) Hz) must lie within (0, $nqf) Hz."))
    end

    # --- IIR: order required; bw not used ---
    if fprototype in (:butterworth, :chebyshev1, :chebyshev2, :elliptic)
        isnothing(order) && throw(ArgumentError("order must be specified for $fprototype."))
        !isnothing(bw) && _info("bw is not used for $fprototype (order must be set manually)")
    end

    # --- :firls weight vector defaults ---
    if fprototype === :firls
        nw = ftype in (:bp, :bs) ? 6 : 4
        if !isnothing(w)
            length(w) == nw || throw(ArgumentError("Length of w must be $nw for :$ftype filter."))
        else
            w = ones(nw)
        end
    end

    # --- ripple defaults for equiripple IIR prototypes ---
    if fprototype in (:chebyshev1, :chebyshev2, :elliptic)
        if isnothing(rp)
            rp = 0.5
            _info("rp set at $rp dB")
        end
        if isnothing(rs)
            rs = 20
            _info("rs set at $rs dB")
        end
    end

    # -----------------------------------------------------------------------
    # FIR filters
    # -----------------------------------------------------------------------

    if fprototype === :fir
        responsetype = if ftype === :lp
            Lowpass(cutoff)
        elseif ftype === :hp
            Highpass(cutoff)
        elseif ftype === :bp
            Bandpass(cutoff[1], cutoff[2])
        elseif ftype === :bs
            Bandstop(cutoff[1], cutoff[2])
        end
        _info("Creating $(uppercase(string(ftype))) FIR filter ($(order) taps)")
        return digitalfilter(responsetype, FIRWindow(w); fs = fs)
    end

    if fprototype === :firls
        if ftype === :bp
            f1_stop, f1_pass = cutoff[1] - bw / 2, cutoff[1] + bw / 2
            f2_pass, f2_stop = cutoff[2] - bw / 2, cutoff[2] + bw / 2
            flt_shape = [0, 0, 1, 1, 0, 0]
            flt_frq = [0, f1_stop, f1_pass, f2_pass, f2_stop, nqf]
            _info("Creating BP FIRLS filter ($order taps, bw=$bw Hz)")
            _info(" Bands: stop=[0,$f1_stop], pass=[$f1_pass,$f2_pass], stop=[$f2_stop,$nqf]")
        elseif ftype === :bs
            f1_pass, f1_stop = cutoff[1] - bw / 2, cutoff[1] + bw / 2
            f2_stop, f2_pass = cutoff[2] - bw / 2, cutoff[2] + bw / 2
            flt_shape = [1, 1, 0, 0, 1, 1]
            flt_frq = [0, f1_pass, f1_stop, f2_stop, f2_pass, nqf]
            _info("Creating BS FIRLS filter ($order taps, bw=$bw Hz)")
            _info(" Bands: pass=[0,$f1_pass], stop=[$f1_stop,$f2_stop], pass=[$f2_pass,$nqf]")
        elseif ftype === :lp
            f_pass, f_stop = cutoff - bw / 2, cutoff + bw / 2
            flt_shape = [1, 1, 0, 0]
            flt_frq = [0, f_pass, f_stop, nqf]
            _info("Creating LP FIRLS filter ($order taps, bw=$bw Hz, pass=$f_pass, stop=$f_stop)")
        elseif ftype === :hp
            f_stop, f_pass = cutoff - bw / 2, cutoff + bw / 2
            flt_shape = [0, 0, 1, 1]
            flt_frq = [0, f_stop, f_pass, nqf]
            _info("Creating HP FIRLS filter ($order taps, bw=$bw Hz, stop=$f_stop, pass=$f_pass)")
        end
        return FIRLSFilterDesign.firls_design(order - 1, flt_frq, flt_shape, w, true; fs = fs)
    end

    if fprototype === :remez
        if ftype === :bp
            f1_stop, f1_pass = cutoff[1] - bw / 2, cutoff[1] + bw / 2
            f2_pass, f2_stop = cutoff[2] - bw / 2, cutoff[2] + bw / 2
            w = [(0, f1_stop) => 0, (f1_pass, f2_pass) => 1, (f2_stop, nqf) => 0]
            _info("Creating BP Remez filter ($order taps, bw=$bw Hz)")
            _info(" Bands: stop=[0,$f1_stop], pass=[$f1_pass,$f2_pass], stop=[$f2_stop,$nqf]")
        elseif ftype === :bs
            f1_pass, f1_stop = cutoff[1] - bw / 2, cutoff[1] + bw / 2
            f2_stop, f2_pass = cutoff[2] - bw / 2, cutoff[2] + bw / 2
            w = [(0, f1_pass) => 1, (f1_stop, f2_stop) => 0, (f2_pass, nqf) => 1]
            _info("Creating BS Remez filter ($order taps, bw=$bw Hz)")
            _info(" Bands: pass=[0,$f1_pass], stop=[$f1_stop,$f2_stop], pass=[$f2_pass,$nqf]")
        elseif ftype === :lp
            f_pass, f_stop = cutoff - bw / 2, cutoff + bw / 2
            w = [(0, f_pass) => 1, (f_stop, nqf) => 0]
            _info("Creating LP Remez filter ($order taps, bw=$bw Hz, pass=$f_pass, stop=$f_stop)")
        elseif ftype === :hp
            f_stop, f_pass = cutoff - bw / 2, cutoff + bw / 2
            w = [(0, f_stop) => 0, (f_pass, nqf) => 1]
            _info("Creating HP Remez filter ($order taps, bw=$bw Hz, stop=$f_stop, pass=$f_pass)")
        end
        return remez(order, w; Hz = fs, maxiter = 100)
    end

    # -----------------------------------------------------------------------
    # IIR filters
    # -----------------------------------------------------------------------

    if fprototype in (:butterworth, :chebyshev1, :chebyshev2, :elliptic)
        responsetype = if ftype === :lp
            Lowpass(cutoff)
        elseif ftype === :hp
            Highpass(cutoff)
        elseif ftype === :bp
            Bandpass(cutoff[1], cutoff[2])
        elseif ftype === :bs
            Bandstop(cutoff[1], cutoff[2])
        end
        prototype = if fprototype === :butterworth
            Butterworth(order)
        elseif fprototype === :chebyshev1
            Chebyshev1(order, rp)
        elseif fprototype === :chebyshev2
            Chebyshev2(order, rs)
        elseif fprototype === :elliptic
            Elliptic(order, rp, rs)
        end
        _info("Creating $(uppercase(string(ftype))) $(fprototype) filter (order=$order)")
        return digitalfilter(responsetype, prototype; fs = fs)
    end

    if fprototype === :iirnotch
        _info("Creating IIR notch filter (cutoff=$cutoff Hz, bw=$bw Hz)")
        return iirnotch(cutoff, bw; fs = fs)
    end
end

"""
    filter_report(flt; <keyword arguments>)

Measure the frequency response of a filter created by `filter_create` and report its parameters.

Gains are those of the filter as applied: single pass, or squared magnitude (dB doubled) for `dir=:twopass`.

# Arguments

- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter object
- `fs::Int64`: sampling rate in Hz
- `dir::Symbol=:twopass`: filtering direction (`:twopass`, `:onepass`, `:reverse`)
- `ftype::Union{Nothing, Symbol}=nothing`: filter type; with `cutoff` and `bw`, pass-band deviation and stop-band attenuation are measured
- `cutoff::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: filter cutoff in Hz
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz (for `:fir` use the effective width, k · fs / taps)
- `fcheck::AbstractVector{<:Real}=Float64[]`: frequencies (Hz) at which to report the gain
- `n::Int64=2^16`: number of frequency points between 0 and Nyquist
- `verbose::Bool=true`: print the report

# Returns

Named tuple:

- `fs`, `dir`, `kind` ("FIR", "IIR", "IIR biquad"), `taps` (FIR only), `order` (FIR: taps; IIR: number of poles)
- `kernel_s` (FIR length in s), `delay_s` (group delay in s; 0 for two-pass)
- `crossings::Dict{Float64, Vector{Float64}}`: frequencies (Hz) where the gain crosses -0.1, -1, -3, -6, -20 and -40 dB
- `gain_dc`, `gain_nyquist`: gain (dB) at 0 Hz and at Nyquist
- `gain_at::Dict{Float64, Float64}`: gain (dB) at `fcheck` frequencies
- `pass_dev`: maximum absolute pass-band deviation (dB)
- `stop_att`: minimum stop-band attenuation (dB)
- `summary::String`: one-line summary (stored in the object history by `filter_apply` and `filter`)
- `f`, `gain_db`: frequency grid and gain
"""
function filter_report(
    flt::_FLT_TYPES;
    fs::Int64,
    dir::Symbol = :twopass,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    bw::Union{Nothing, Real} = nothing,
    fcheck::AbstractVector{<:Real} = Float64[],
    n::Int64 = 2^16,
    verbose::Bool = true,
)
    _check_var(dir, [:twopass, :onepass, :reverse], "dir")
    nqf = fs / 2
    f = collect(range(0, nqf; length = n))
    ω = 2π .* f ./ fs
    h = flt isa Vector{Float64} ? freqresp(PolynomialRatio(flt, [1.0]), ω) : freqresp(flt, ω)
    g = 20 .* log10.(max.(abs.(h), floatmin(Float64)))
    dir === :twopass && (g .*= 2)

    # threshold crossings (linear interpolation)
    thr = [-0.1, -1.0, -3.0, -6.0, -20.0, -40.0]
    crossings = Dict{Float64, Vector{Float64}}()
    for t in thr
        x = Float64[]
        for i in 1:(n - 1)
            if (g[i] - t) * (g[i + 1] - t) < 0
                push!(x, f[i] + (t - g[i]) * (f[i + 1] - f[i]) / (g[i + 1] - g[i]))
            end
        end
        crossings[t] = round.(x; digits = 3)
    end

    gain_at = Dict{Float64, Float64}()
    for fc in fcheck
        (0 <= fc <= nqf) || throw(ArgumentError("fcheck frequencies must be within [0, $nqf] Hz."))
        gain_at[fc] = round(g[clamp(round(Int64, fc / nqf * (n - 1)) + 1, 1, n)]; digits = 2)
    end

    # pass-band deviation and stop-band attenuation
    pass_dev = nothing
    stop_att = nothing
    if !isnothing(ftype) && !isnothing(cutoff) && !isnothing(bw)
        tr = _transition_edges(ftype, cutoff, bw)
        pass, stop = if ftype === :lp
            f .<= tr[1][1], f .>= tr[1][2]
        elseif ftype === :hp
            f .>= tr[1][2], f .<= tr[1][1]
        elseif ftype === :bp
            (f .>= tr[1][2]) .& (f .<= tr[2][1]), (f .<= tr[1][1]) .| (f .>= tr[2][2])
        else
            (f .<= tr[1][1]) .| (f .>= tr[2][2]), (f .>= tr[1][2]) .& (f .<= tr[2][1])
        end
        any(pass) && (pass_dev = round(maximum(abs.(g[pass])); digits = 3))
        any(stop) && (stop_att = round(-maximum(g[stop]); digits = 1))
    end

    taps = flt isa Vector{Float64} ? length(flt) : nothing
    # IIR order = number of poles (Biquad: 2)
    iir_order = flt isa Vector{Float64} ? nothing : (flt isa Biquad ? 2 : length(flt.p))
    kind = flt isa Vector{Float64} ? "FIR" : (flt isa Biquad ? "IIR biquad" : "IIR")
    kernel_s = isnothing(taps) ? nothing : round((taps - 1) / fs; digits = 3)
    delay_s = dir === :twopass ? 0.0 : (isnothing(taps) ? nothing : round((taps - 1) / 2 / fs; digits = 3))

    if verbose
        _info("$kind filter response ($(dir === :twopass ? "two-pass, |H|²" : String(dir))), fs=$fs Hz")
        !isnothing(taps) && _info("  taps: $taps; kernel length: $kernel_s s; group delay: $delay_s s")
        !isnothing(iir_order) && _info("  order: $iir_order$(dir === :twopass ? " (effective $(2 * iir_order) two-pass)" : "")")
        for t in thr
            _info("  $t dB at: $(isempty(crossings[t]) ? "not reached" : Base.join(crossings[t], ", ") * " Hz")")
        end
        _info("  gain at 0 Hz: $(round(g[1]; digits = 2)) dB; at Nyquist: $(round(g[end]; digits = 2)) dB")
        for fc in fcheck
            _info("  gain at $fc Hz: $(gain_at[fc]) dB")
        end
        !isnothing(pass_dev) && _info("  max pass-band deviation: $pass_dev dB")
        !isnothing(stop_att) && _info("  min stop-band attenuation: $stop_att dB")
    end

    summary =
        "$kind" * (isnothing(taps) ? ", order=$iir_order" : ", taps=$taps") * ", dir=:$dir" *
        "; -3 dB at $(isempty(crossings[-3.0]) ? "n/a" : Base.join(crossings[-3.0], ", ")) Hz" *
        "; -6 dB at $(isempty(crossings[-6.0]) ? "n/a" : Base.join(crossings[-6.0], ", ")) Hz" *
        "; 0 Hz: $(round(g[1]; digits = 2)) dB" *
        (isnothing(pass_dev) ? "" : "; pass-band dev ≤ $pass_dev dB") *
        (isnothing(stop_att) ? "" : "; stop-band att ≥ $stop_att dB") *
        Base.join(["; $fc Hz: $(gain_at[fc]) dB" for fc in fcheck])

    return (
        fs = fs,
        dir = dir,
        kind = kind,
        taps = taps,
        order = isnothing(taps) ? iir_order : taps,
        kernel_s = kernel_s,
        delay_s = delay_s,
        crossings = crossings,
        gain_dc = g[1],
        gain_nyquist = g[end],
        gain_at = gain_at,
        pass_dev = pass_dev,
        stop_att = stop_att,
        summary = summary,
        f = f,
        gain_db = g,
    )
end

"""
    filter_apply(s; <keyword arguments>)

Apply a pre-designed IIR or FIR filter to a signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter object returned by `filter_create`
- `dir:Symbol=:twopass`: filtering direction:
    - `:twopass`: forward pass followed by reverse pass (zero phase distortion; magnitude response is squared)
    - `:onepass`: single forward pass (introduces phase delay)
    - `:reverse`: single reverse pass
- `fs::Union{Nothing, Int64}=nothing`: sampling rate in Hz; required if `report=true`
- `report::Bool=false`: print the filter properties and measured frequency response (see `filter_report`)
- `ftype::Union{Nothing, Symbol}=nothing`, `cutoff::Union{Nothing, Real, Tuple{Real, Real}}=nothing`, `bw::Union{Nothing, Real}=nothing`: optional design parameters; if all three are given, pass-band deviation and stop-band attenuation are also reported
- `fcheck::AbstractVector{<:Real}=Float64[]`: frequencies (Hz) at which to report the gain

# Returns

- `Vector{Float64}`: filtered signal of the same length as `s`
"""
function filter_apply(
    s::AbstractVector;
    flt::_FLT_TYPES,
    dir::Symbol = :twopass,
    fs::Union{Nothing, Int64} = nothing,
    report::Bool = false,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    bw::Union{Nothing, Real} = nothing,
    fcheck::AbstractVector{<:Real} = Float64[],
)::Vector{Float64}
    # validate
    _check_var(dir, [:twopass, :onepass, :reverse], "dir")
    if report
        isnothing(fs) && throw(ArgumentError("fs must be specified when report=true."))
        filter_report(flt; fs = fs, dir = dir, ftype = ftype, cutoff = cutoff, bw = bw, fcheck = fcheck, verbose = true)
    end

    if flt isa Vector{Float64} && length(s) <= 3 * (length(flt) - 1)
        _warn(
            "Signal ($(length(s)) samples) is short relative to the filter ($(length(flt)) taps): " *
            "edge transients will affect a large part of it. Filter the continuous signal before epoching.",
        )
    end

    if dir === :onepass
        return filt(flt, s)
    elseif dir === :twopass
        return filtfilt(flt, s)
    elseif dir === :reverse
        return reverse(filt(flt, reverse(s)))
    end
end

"""
    filter_apply(obj; <keyword arguments>)

Apply a pre-designed filter to selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter object returned by `filter_create`
- `dir:Symbol=:twopass`: filtering direction:
    - `:twopass`: forward pass followed by reverse pass (zero phase distortion; magnitude response is squared)
    - `:onepass`: single forward pass (introduces phase delay)
    - `:reverse`: single reverse pass
- `report::Bool=true`: print the filter properties and measured frequency response (see `filter_report`)
- `ftype::Union{Nothing, Symbol}=nothing`, `cutoff::Union{Nothing, Real, Tuple{Real, Real}}=nothing`, `bw::Union{Nothing, Real}=nothing`: optional design parameters; if all three are given, pass-band deviation and stop-band attenuation are also reported
- `fcheck::AbstractVector{<:Real}=Float64[]`: frequencies (Hz) at which to report the gain

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object with filtered data; the measured filter properties are stored in `obj.history`

# Notes

- For best results apply to a continuous (single-epoch) signal. A warning is issued when `nepochs(obj) > 1`.
- Taper the signal before filtering to reduce edge artifacts.
"""
function filter_apply(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    flt::_FLT_TYPES,
    dir::Symbol = :twopass,
    report::Bool = true,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    bw::Union{Nothing, Real} = nothing,
    fcheck::AbstractVector{<:Real} = Float64[],
)::NeuroAnalyzer.NEURO
    # validate
    _check_var(dir, [:twopass, :onepass, :reverse], "dir")

    # filter properties: always measured (for history), printed if report=true
    r = filter_report(
        flt;
        fs = sr(obj),
        dir = dir,
        ftype = ftype,
        cutoff = cutoff,
        bw = bw,
        fcheck = fcheck,
        n = 2^14,
        verbose = report,
    )

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    ch_n = length(ch)
    ep_n = nepochs(obj)

    ep_n > 1 && _warn("filter_apply() should preferably be used on a continuous signal.")
    _info("Taper the signal before filtering to reduce edge artifacts")
    dir === :twopass && _info("Two-pass filtering: magnitude response is squared (gains in dB doubled)")
    if flt isa Vector{Float64} && size(obj.data, 2) <= 3 * (length(flt) - 1)
        _warn(
            "Epoch length ($(size(obj.data, 2)) samples) is short relative to the filter ($(length(flt)) taps): " *
            "edge transients will affect a large part of each epoch.",
        )
    end

    # create new dataset
    obj_tmp = deepcopy(obj)

    # initialize progress bar
    progbar = Progress(ep_n * ch_n; dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        obj_tmp.data[ch[ch_idx], :, ep_idx] = _filter_one(obj, ch[ch_idx], ep_idx, flt, dir)
        # update progress bar
        progress_bar && next!(progbar)
    end

    push!(obj_tmp.history, "filter_apply(obj; ch=$ch, dir=:$dir) [$(r.summary)]")

    return obj_tmp
end

# per-channel/epoch worker (avoids repeating the length warning in every thread)
function _filter_one(obj::NeuroAnalyzer.NEURO, ch_idx::Int64, ep_idx::Int64, flt::_FLT_TYPES, dir::Symbol)
    s = @view(obj.data[ch_idx, :, ep_idx])
    dir === :onepass && return filt(flt, s)
    dir === :twopass && return filtfilt(flt, s)
    return reverse(filt(flt, reverse(s)))
end

"""
    filter_apply!(obj; <keyword arguments>)

Apply a pre-designed filter in-place to selected channels of a NEURO object.

Delegates to `filter_apply` and copies the result back.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter object returned by `filter_create`
- `dir:Symbol=:twopass`: filtering direction (`:twopass`, `:onepass`, `:reverse`)
- `report::Bool=true`: print the filter properties and measured frequency response (see `filter_report`)
- `ftype::Union{Nothing, Symbol}=nothing`, `cutoff::Union{Nothing, Real, Tuple{Real, Real}}=nothing`, `bw::Union{Nothing, Real}=nothing`: optional design parameters; if all three are given, pass-band deviation and stop-band attenuation are also reported
- `fcheck::AbstractVector{<:Real}=Float64[]`: frequencies (Hz) at which to report the gain

# Returns

- `Nothing`
"""
function filter_apply!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    flt::_FLT_TYPES,
    dir::Symbol = :twopass,
    report::Bool = true,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    bw::Union{Nothing, Real} = nothing,
    fcheck::AbstractVector{<:Real} = Float64[],
)::Nothing
    obj_tmp = filter_apply(
        obj;
        ch = ch,
        flt = flt,
        dir = dir,
        report = report,
        ftype = ftype,
        cutoff = cutoff,
        bw = bw,
        fcheck = fcheck,
    )
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end

"""
    filter(obj; <keyword arguments>)

Design and apply a digital filter to selected channels of a NEURO object in a single call.

Combines `filter_create`, `filter_report` and `filter_apply`. When `preview=true`, the filter frequency response is plotted without modifying the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `fprototype::Symbol`: filter prototype (`:fir`, `:firls`, `:remez`, `:butterworth`, `:chebyshev1`, `:chebyshev2`, `:elliptic`, `:iirnotch`)
- `ftype::Union{Nothing, Symbol}=nothing`: filter type (`:lp`, `:hp`, `:bp`, `:bs`)
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz
    - for `:lp`/`:hp`: single frequency
    - for `:bp`/`:bs`: frequency range (f1, f2)
- `order::Union{Nothing, Int64}=nothing`: filter order (number of taps for FIR, filter order for IIR); FIR: calculated from `bw` if `nothing`
- `rp::Union{Nothing, Real}=nothing`: pass-band ripple in dB (default 0.5 dB)
- `rs::Union{Nothing, Real}=nothing`: stop-band attenuation in dB (default 20 dB for IIR)
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz
- `w::Union{Nothing, AbstractVector}=nothing`: window vector for `:fir` or weight vector for `:firls`
- `window::Symbol=:hamming`: window for `:fir` (`:rect`, `:hann`, `:hamming`, `:blackman`)
- `dir:Symbol=:twopass`: filtering direction (`:twopass`, `:onepass`, `:reverse`)
- `report::Bool=true`: print the filter parameters and measured frequency response
- `fcheck::AbstractVector{<:Real}=Float64[]`: frequencies (Hz) at which to report the gain
- `preview::Bool=false`: if `true`, plot the filter frequency response and return the figure without filtering the signal

# Returns

- `NeuroAnalyzer.NEURO`: filtered object (when `preview=false`)
- `GLMakie.Figure`: filter frequency-response plot (when `preview=true`)
"""
function filter(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    fprototype::Symbol,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Real, Tuple{Real, Real}},
    order::Union{Nothing, Int64} = nothing,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    bw::Union{Nothing, Real} = nothing,
    w::Union{Nothing, AbstractVector} = nothing,
    window::Symbol = :hamming,
    dir::Symbol = :twopass,
    report::Bool = true,
    fcheck::AbstractVector{<:Real} = Float64[],
    preview::Bool = false,
)::Union{NeuroAnalyzer.NEURO, GLMakie.Figure}
    fs = sr(obj)

    if preview
        _info("Previewing filter response, signal will not be filtered")
        return plot_filter(;
            fs = fs,
            fprototype = fprototype,
            ftype = ftype,
            cutoff = cutoff,
            order = order,
            rp = rp,
            rs = rs,
            bw = bw,
            w = w,
            window = window,
            dir = dir,
        )
    end

    # resolve FIR order here so that its source can be reported
    order_src = isnothing(order) ? "auto" : "manual"
    if fprototype in (:fir, :firls, :remez) && isnothing(order) && isnothing(w) && !isnothing(bw)
        order = filter_order(; fprototype = fprototype, fs = fs, bw = bw, window = window, rs = rs)
        order_src = "calculated from bw"
    elseif fprototype === :fir && !isnothing(w)
        order = length(w)
        order_src = "length(w)"
    end

    flt = filter_create(;
        fprototype = fprototype,
        ftype = ftype,
        cutoff = cutoff,
        fs = fs,
        order = order,
        rp = rp,
        rs = rs,
        bw = bw,
        w = w,
        window = window,
    )

    # effective transition width
    bw_eff = if fprototype === :fir && isnothing(w)
        _FIR_WINDOWS[window].k * fs / length(flt)
    elseif fprototype in (:firls, :remez, :iirnotch)
        bw
    else
        nothing
    end

    desc =
        "filter(obj; ch=$ch, fprototype=:$fprototype" *
        (isnothing(ftype) ? "" : ", ftype=:$ftype") *
        ", cutoff=$cutoff, fs=$fs" *
        (isnothing(order) ? "" : ", order=$order ($order_src)") *
        (isnothing(bw) ? "" : ", bw=$bw") *
        (isnothing(bw_eff) ? "" : ", bw_effective=$(round(bw_eff; digits = 3))") *
        (fprototype === :fir ? ", window=$(isnothing(w) ? ":$window" : "custom")" : "") *
        (isnothing(rp) ? "" : ", rp=$rp") *
        (isnothing(rs) ? "" : ", rs=$rs") *
        ", dir=:$dir)"

    # design parameters, then measured response (printed by filter_apply)
    # report && _info(desc)
    obj_tmp = filter_apply(
        obj;
        ch = ch,
        flt = flt,
        dir = dir,
        report = report,
        ftype = fprototype === :iirnotch ? nothing : ftype,
        cutoff = cutoff,
        bw = bw_eff,
        fcheck = fcheck,
    )
    # history: design line followed by the filter_apply line with the measured summary
    insert!(obj_tmp.history, length(obj_tmp.history), desc)

    return obj_tmp
end

"""
    filter!(obj; <keyword arguments>)

Design and apply a digital filter in-place to selected channels of a NEURO object.

Arguments as in `filter`. When `preview=true`, the filter frequency response is plotted and returned without modifying the signal.

# Returns

- `Nothing` when `preview=false`
- `GLMakie.Figure`: filter frequency-response plot (when `preview=true`)
"""
function filter!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    fprototype::Symbol,
    ftype::Union{Symbol, Nothing} = nothing,
    cutoff::Union{Real, Tuple{Real, Real}},
    order::Union{Nothing, Int64} = nothing,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    bw::Union{Nothing, Real} = nothing,
    w::Union{Nothing, AbstractVector} = nothing,
    window::Symbol = :hamming,
    dir::Symbol = :twopass,
    report::Bool = true,
    fcheck::AbstractVector{<:Real} = Float64[],
    preview::Bool = false,
)::Union{Nothing, GLMakie.Figure}
    out = NeuroAnalyzer.filter(
        obj;
        ch = ch,
        fprototype = fprototype,
        ftype = ftype,
        cutoff = cutoff,
        order = order,
        rp = rp,
        rs = rs,
        bw = bw,
        w = w,
        window = window,
        dir = dir,
        report = report,
        fcheck = fcheck,
        preview = preview,
    )
    preview && return out

    obj.data = out.data
    obj.history = out.history

    return nothing
end