export fooof

# to be appended to test/analyze.jl
# reference values from fooof 1.1.0 (Python) on the same noise-free spectra:
#   fixed: ap = [1.00846, 1.50140], peaks = [10.0016 0.78872 2.91857; 22.0024 0.38969 3.80960], R² = 0.99986
#   knee:  ap = [1.0, 25.0, 2.0], no peaks
 
@info "Test: fooof()"
f = collect(1:0.5:40)
p = 10 .^ (1 .- 1.5 .* log10.(f) .+ 0.8 .* exp.(-(f .- 10) .^ 2 ./ (2 * 1.5^2)) .+ 0.4 .* exp.(-(f .- 22) .^ 2 ./ (2 * 2.0^2)))
r = fooof(f, p)
@test isapprox(r.ap, [1.00846, 0.0, 1.50140], atol = 1e-3)
@test size(r.pe) == (2, 3)
@test isapprox(r.pe, [10.0016 0.78872 2.91857; 22.0024 0.38969 3.80960], atol = 1e-3)
@test isapprox(r.r2, 0.99986, atol = 1e-4)
@test length(r.model) == length(r.f) == length(f)
pk = 10 .^ (1 .- log10.(25 .+ f .^ 2))
r = fooof(f, pk, aperiodic_mode = :knee)
@test isapprox(r.ap, [1.0, 25.0, 2.0], rtol = 1e-3)
@test size(r.pe) == (0, 3)
r = fooof(f, p, flim = (2, 30), max_n_peaks = 1)
@test r.f[1] == 2.0
@test r.f[end] == 30.0
@test size(r.pe, 1) == 1
@test_throws AssertionError fooof(f, pow2db.(p))
@test_throws AssertionError fooof(vcat(0.0, f), vcat(1.0, p), flim = (0, 40))
r = fooof(rand(10, 1000, 10), fs = 100)
@test size(r.ap) == (10, 3, 10)
@test size(r.pe) == (10, 10)
@test size(r.r2) == (10, 10)
@test r.f[1] == 1.0
@test r.f[end] == 50.0
r = fooof(rand(10, 1000, 10), fs = 100, method = :stft)
@test size(r.ap) == (10, 3, 10)
@test r.f[1] == 1.0
r = fooof(e10, ch = "Fp1")
@test size(r.ap) == (1, 3, 10)
@test size(r.model) == (1, length(r.f), 10)

# FOOOF / specparam: parameterizing neural power spectra into aperiodic and periodic components
#
# Reference:
#   Donoghue T, Haller M, Peterson EJ, et al. Parameterizing neural power spectra into periodic
#   and aperiodic components. Nat Neurosci. 2020;23(12):1655-1665. doi:10.1038/s41593-020-00744-x
#
# This is a re-implementation of the fitting algorithm of fooof 1.1.0 / specparam 2.0.0rc
# (the core algorithm and its default settings are identical in both), with these differences:
#   - curve fitting uses an internal box-constrained Levenberg-Marquardt solver instead of
#     scipy.optimize.curve_fit (MINPACK LM for the unbounded aperiodic fit, TRF for the bounded
#     peak fit); results agree with the Python reference to within fitting tolerance, but are
#     not bit-identical
#   - for `aperiodic_mode=:fixed` the knee is reported as 0 (the fixed model is the knee model
#     with knee = 0), so that the output always has 3 aperiodic parameters

# private settings (as in fooof/specparam)
const _FOOOF_AP_PERCENTILE_THRESH = 0.025   # percentile (not fraction!), passed as in np.percentile()
const _FOOOF_BW_STD_EDGE = 1.0              # drop peaks closer to the edge than this × gaussian SD
const _FOOOF_GAUSS_OVERLAP_THRESH = 0.75    # drop the lower of two peaks overlapping by this × SD
const _FOOOF_CF_BOUND = 1.5                 # CF bounds: ± 2 × cf_bound × gaussian SD
const _FOOOF_MAXFEV = 5000                  # maximum number of model evaluations per fit
const _FOOOF_TOL = 1e-5                     # convergence tolerance

# model functions (power in log10)

# log10 returning NaN instead of throwing DomainError (as numpy does); needed because the solver may
# try knee < -f^exponent, which must be rejected as a bad step, not abort the fit
_log10nan(x::Real)::Float64 = x > 0 ? log10(x) : NaN

function _fooof_ap(f::AbstractVector{<:Real}, p::AbstractVector{<:Real})::Vector{Float64}
    # fixed: [offset, exponent]; knee: [offset, knee, exponent]
    if length(p) == 2
        return p[1] .- p[2] .* log10.(f)
    else
        return p[1] .- _log10nan.(p[2] .+ f .^ p[3])
    end
end

function _fooof_ap_jac(f::AbstractVector{<:Real}, p::AbstractVector{<:Real})::Matrix{Float64}
    if length(p) == 2
        J = zeros(length(f), 2)
        J[:, 1] .= 1.0
        J[:, 2] .= -log10.(f)
    else
        J = zeros(length(f), 3)
        fx = f .^ p[3]
        d = (p[2] .+ fx) .* log(10)
        J[:, 1] .= 1.0
        J[:, 2] .= -1.0 ./ d
        J[:, 3] .= -fx .* log.(f) ./ d
    end
    return J
end

function _fooof_gauss(f::AbstractVector{<:Real}, p::AbstractVector{<:Real})::Vector{Float64}
    # p = [cf1, h1, sd1, cf2, h2, sd2, ...]
    y = zeros(length(f))
    @inbounds for idx in 1:3:length(p)
        y .+= p[idx + 1] .* exp.(-(f .- p[idx]) .^ 2 ./ (2 * p[idx + 2]^2))
    end
    return y
end

function _fooof_gauss_jac(f::AbstractVector{<:Real}, p::AbstractVector{<:Real})::Matrix{Float64}
    J = zeros(length(f), length(p))
    @inbounds for idx in 1:3:length(p)
        c, h, w = p[idx], p[idx + 1], p[idx + 2]
        dx = f .- c
        e = exp.(-dx .^ 2 ./ (2 * w^2))
        J[:, idx] .= h .* e .* dx ./ w^2
        J[:, idx + 1] .= e
        J[:, idx + 2] .= h .* e .* dx .^ 2 ./ w^3
    end
    return J
end

# box-constrained Levenberg-Marquardt (bounds enforced by projection)

function _fooof_lm(
        fun::Function,
        jac::Function,
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        p0::AbstractVector{<:Real},
        lb::AbstractVector{<:Real},
        ub::AbstractVector{<:Real};
        maxfev::Int64 = _FOOOF_MAXFEV,
        tol::Real = _FOOOF_TOL,
    )::Vector{Float64}

    p = clamp.(Float64.(p0), lb, ub)
    r = y .- fun(x, p)
    sse = sum(abs2, r)
    @assert isfinite(sse) "Model is not finite at the initial guess."
    nfev = 1
    λ = 1e-3

    while nfev < maxfev
        J = jac(x, p)
        g = J' * r
        A = J' * J
        # gradient convergence
        maximum(abs, g) <= tol && return p
        d = max.(diag(A), eps())
        converged = false
        while true
            nfev >= maxfev && error("Maximum number of function evaluations ($maxfev) exceeded.")
            δ = try
                Symmetric(A + λ * Diagonal(d)) \ g
            catch
                nothing
            end
            if isnothing(δ) || !all(isfinite, δ)
                λ *= 10
                λ > 1e16 && return p
                continue
            end
            p_new = clamp.(p .+ δ, lb, ub)
            r_new = y .- fun(x, p_new)
            nfev += 1
            sse_new = sum(abs2, r_new)
            if isfinite(sse_new) && sse_new < sse
                converged = (sse - sse_new) <= tol * sse ||
                    norm(p_new .- p) <= tol * (norm(p) + tol)
                p, r, sse = p_new, r_new, sse_new
                λ = max(λ / 10, 1e-12)
                break
            else
                λ *= 10
                # no descent direction left: p is a (constrained) minimum
                λ > 1e16 && return p
            end
        end
        converged && return p
    end

    error("Maximum number of function evaluations ($maxfev) exceeded.")

end

# aperiodic fitting

function _fooof_ap_guess(f::AbstractVector{<:Real}, p_fit::AbstractVector{<:Real}, p_full::AbstractVector{<:Real}, knee::Bool)::Vector{Float64}
    # note: as in fooof, the exponent guess is always taken from the full (original) spectrum
    off = p_fit[1]
    ex = abs((p_full[end] - p_full[1]) / (log10(f[end]) - log10(f[1])))
    return knee ? [off, 0.0, ex] : [off, ex]
end

function _fooof_ap_fit(f::AbstractVector{<:Real}, p::AbstractVector{<:Real}, p0::AbstractVector{<:Real})::Vector{Float64}
    n = length(p0)
    return _fooof_lm(_fooof_ap, _fooof_ap_jac, f, p, p0, fill(-Inf, n), fill(Inf, n))
end

function _fooof_robust_ap_fit(f::AbstractVector{<:Real}, p::AbstractVector{<:Real}, knee::Bool)::Vector{Float64}
    # initial fit
    popt = _fooof_ap_fit(f, p, _fooof_ap_guess(f, p, p, knee))
    # flatten, clip negative values, keep points at or below the percentile threshold
    flat = p .- _fooof_ap(f, popt)
    flat[flat .< 0] .= 0
    # Julia quantile() default (type 7) == numpy percentile() default (linear)
    thr = quantile(flat, _FOOOF_AP_PERCENTILE_THRESH / 100)
    mask = flat .<= thr
    @assert count(mask) >= length(popt) "Too few points for the robust aperiodic fit."
    return _fooof_ap_fit(f[mask], p[mask], popt)
end

# periodic fitting

function _fooof_fit_peaks(
        f::AbstractVector{<:Real},
        flat::AbstractVector{<:Real},
        peak_width_limits::Tuple{Real, Real},
        max_n_peaks::Int64,
        min_peak_height::Real,
        peak_threshold::Real,
    )::Matrix{Float64}

    sd_lim = peak_width_limits ./ 2
    fres = f[2] - f[1]
    frange = (f[1], f[end])
    flat_iter = copy(flat)
    guess = zeros(0, 3)

    # find candidate peaks iteratively
    while size(guess, 1) < max_n_peaks
        max_idx = argmax(flat_iter)
        max_h = flat_iter[max_idx]
        # relative (to SD of the flattened spectrum) and absolute height thresholds
        max_h <= peak_threshold * std(flat_iter, corrected = false) && break
        max_h > min_peak_height || break
        # data-driven SD guess from the shorter half-width at half-maximum
        # note: as in fooof, the left-side search excludes the first frequency bin
        hh = 0.5 * max_h
        le_idx = nothing
        for i in (max_idx - 1):-1:2
            if flat_iter[i] <= hh
                le_idx = i
                break
            end
        end
        ri_idx = nothing
        for i in (max_idx + 1):length(flat_iter)
            if flat_iter[i] <= hh
                ri_idx = i
                break
            end
        end
        sides = [abs(i - max_idx) for i in (le_idx, ri_idx) if !isnothing(i)]
        if isempty(sides)
            g_sd = mean(peak_width_limits)
        else
            fwhm = minimum(sides) * 2 * fres
            g_sd = fwhm / (2 * sqrt(2 * log(2)))
        end
        g_sd = clamp(g_sd, sd_lim[1], sd_lim[2])
        guess = vcat(guess, [f[max_idx] max_h g_sd])
        flat_iter = flat_iter .- _fooof_gauss(f, [f[max_idx], max_h, g_sd])
    end

    # drop peaks too close to the edges
    if size(guess, 1) > 0
        keep = (abs.(guess[:, 1] .- frange[1]) .> guess[:, 3] .* _FOOOF_BW_STD_EDGE) .&
            (abs.(guess[:, 1] .- frange[2]) .> guess[:, 3] .* _FOOOF_BW_STD_EDGE)
        guess = guess[keep, :]
    end

    # drop the lower of two overlapping peaks (pairs evaluated on the sorted list, then dropped together)
    if size(guess, 1) > 1
        guess = guess[sortperm(guess[:, 1]), :]
        lo = guess[:, 1] .- guess[:, 3] .* _FOOOF_GAUSS_OVERLAP_THRESH
        hi = guess[:, 1] .+ guess[:, 3] .* _FOOOF_GAUSS_OVERLAP_THRESH
        drop = Int64[]
        for idx in 1:(size(guess, 1) - 1)
            hi[idx] > lo[idx + 1] && push!(drop, guess[idx, 2] <= guess[idx + 1, 2] ? idx : idx + 1)
        end
        guess = guess[setdiff(1:size(guess, 1), drop), :]
    end

    size(guess, 1) == 0 && return zeros(0, 3)

    # fit all gaussians together
    lb = Float64[]
    ub = Float64[]
    for idx in axes(guess, 1)
        cf, _, sd = guess[idx, :]
        append!(lb, [max(cf - 2 * _FOOOF_CF_BOUND * sd, frange[1]), 0.0, sd_lim[1]])
        append!(ub, [min(cf + 2 * _FOOOF_CF_BOUND * sd, frange[2]), Inf, sd_lim[2]])
    end
    gp = _fooof_lm(_fooof_gauss, _fooof_gauss_jac, f, flat, vec(permutedims(guess)), lb, ub)
    gp = permutedims(reshape(gp, 3, :))

    return gp[sortperm(gp[:, 1]), :]

end

"""
    fooof(f, p; <keyword arguments>)

Parameterize a power spectrum into aperiodic and periodic (peak) components using the FOOOF / specparam algorithm (Donoghue et al., 2020). Power is modeled in log10 space as `L(f) = offset - log10(knee + f^exponent) + Σ gaussians`.

# Arguments

  - `f::AbstractVector`: frequencies (in Hz, evenly spaced, > 0)
  - `p::AbstractVector`: powers (linear scale, not dB)
  - `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency range to fit
  - `peak_width_limits::Tuple{Real, Real}=(0.5, 12.0)`: limits on peak bandwidth (in Hz)
  - `max_n_peaks::Int64=typemax(Int64)`: maximum number of peaks
  - `min_peak_height::Real=0.0`: absolute threshold for peak height (in log10 power above the aperiodic component)
  - `peak_threshold::Real=2.0`: relative threshold for peak detection (in SD of the flattened spectrum)
  - `aperiodic_mode::Symbol=:fixed`: aperiodic model:
      + `:fixed`: no knee, `offset - log10(f^exponent)`
      + `:knee`: with knee, `offset - log10(knee + f^exponent)`

# Returns

Named tuple containing:

  - `ap::Vector{Float64}`: aperiodic parameters: offset, knee (0 for `:fixed`), exponent
  - `pe::Matrix{Float64}`: peak parameters (one row per peak): center frequency (Hz), power above the aperiodic component (log10), bandwidth (Hz, 2 × gaussian SD)
  - `gp::Matrix{Float64}`: gaussian parameters (one row per peak): center, height, SD
  - `r2::Float64`: R² of the model fit
  - `err::Float64`: mean absolute error of the model fit
  - `f::Vector{Float64}`: frequencies of the fit
  - `p::Vector{Float64}`: fitted spectrum (log10)
  - `model::Vector{Float64}`: full model (log10)
  - `ap_fit::Vector{Float64}`: aperiodic component (log10)

If fitting fails, all parameters are `NaN` and `pe` / `gp` are empty.
"""
function fooof(
        f::AbstractVector{<:Real},
        p::AbstractVector{<:Real};
        flim::Tuple{Real, Real} = (f[1], f[end]),
        peak_width_limits::Tuple{Real, Real} = (0.5, 12.0),
        max_n_peaks::Int64 = typemax(Int64),
        min_peak_height::Real = 0.0,
        peak_threshold::Real = 2.0,
        aperiodic_mode::Symbol = :fixed,
    )::@NamedTuple{
        ap::Vector{Float64},
        pe::Matrix{Float64},
        gp::Matrix{Float64},
        r2::Float64,
        err::Float64,
        f::Vector{Float64},
        p::Vector{Float64},
        model::Vector{Float64},
        ap_fit::Vector{Float64},
    }

    _check_var(aperiodic_mode, [:fixed, :knee], "aperiodic_mode")
    @assert length(f) == length(p) "f and p must have the same length."
    @assert flim[1] > 0 "flim[1] must be > 0 (log10 of 0 Hz is undefined)."
    _check_tuple(flim, (f[1], f[end]), "flim")
    _check_tuple(peak_width_limits, (0, Inf), "peak_width_limits")
    @assert peak_width_limits[1] > 0 "peak_width_limits[1] must be > 0."
    @assert max_n_peaks >= 0 "max_n_peaks must be ≥ 0."
    @assert peak_threshold > 0 "peak_threshold must be > 0."

    # inclusive frequency range, as in fooof
    f_idx = (f .>= flim[1]) .& (f .<= flim[2])
    f = Float64.(f[f_idx])
    p = Float64.(p[f_idx])
    _fooof_check(f, p)
    _fooof_width_warn(f, peak_width_limits)

    return _fooof(
        f,
        p;
        peak_width_limits = peak_width_limits,
        max_n_peaks = max_n_peaks,
        min_peak_height = min_peak_height,
        peak_threshold = peak_threshold,
        aperiodic_mode = aperiodic_mode,
        verbose = true,
    )

end

function _fooof_nan(f::Vector{Float64}, p::Vector{Float64})
    nans = fill(NaN, length(f))
    pl = all(p .> 0) ? log10.(p) : copy(nans)
    return (
        ap = fill(NaN, 3),
        pe = zeros(0, 3),
        gp = zeros(0, 3),
        r2 = NaN,
        err = NaN,
        f = f,
        p = pl,
        model = nans,
        ap_fit = copy(nans),
    )
end

function _fooof_check(f::Vector{Float64}, p::Vector{Float64})::Nothing
    @assert length(f) >= 4 "At least 4 frequency bins are required within flim."
    @assert all(x -> isapprox(x, f[2] - f[1]), diff(f)) "Frequencies must be evenly spaced."
    @assert all(isfinite, p) "p must not contain NaN or Inf values."
    @assert all(p .> 0) "p must be > 0 (power in linear scale, not dB)."
    return nothing
end

function _fooof_width_warn(f::Vector{Float64}, peak_width_limits::Tuple{Real, Real})::Nothing
    fres = f[2] - f[1]
    peak_width_limits[1] < 2 * fres && _warn(
        "Lower peak width limit ($(peak_width_limits[1]) Hz) is < 2 × frequency resolution ($(round(2 * fres, digits = 3)) Hz); this may lead to overfitting of narrow peaks."
    )
    return nothing
end

function _fooof(
        f::Vector{Float64},
        p::Vector{Float64};
        peak_width_limits::Tuple{Real, Real},
        max_n_peaks::Int64,
        min_peak_height::Real,
        peak_threshold::Real,
        aperiodic_mode::Symbol,
        verbose::Bool = true,
    )::@NamedTuple{
        ap::Vector{Float64},
        pe::Matrix{Float64},
        gp::Matrix{Float64},
        r2::Float64,
        err::Float64,
        f::Vector{Float64},
        p::Vector{Float64},
        model::Vector{Float64},
        ap_fit::Vector{Float64},
    }

    pl = log10.(p)
    knee = aperiodic_mode === :knee

    try
        # robust initial aperiodic fit
        ap = _fooof_robust_ap_fit(f, pl, knee)
        # flatten and fit peaks
        flat = pl .- _fooof_ap(f, ap)
        gp = _fooof_fit_peaks(f, flat, peak_width_limits, max_n_peaks, min_peak_height, peak_threshold)
        pe_fit = _fooof_gauss(f, vec(permutedims(gp)))
        # final aperiodic fit on the peak-removed spectrum
        p_rm = pl .- pe_fit
        ap = _fooof_ap_fit(f, p_rm, _fooof_ap_guess(f, p_rm, pl, knee))
        ap_fit = _fooof_ap(f, ap)
        model = ap_fit .+ pe_fit
        # peak parameters: CF, power above the aperiodic fit at the nearest bin, bandwidth = 2 × SD
        pe = zeros(size(gp, 1), 3)
        for idx in axes(gp, 1)
            f_idx = argmin(abs.(f .- gp[idx, 1]))
            pe[idx, :] = [gp[idx, 1], model[f_idx] - ap_fit[f_idx], 2 * gp[idx, 3]]
        end
        r2 = cor(pl, model)^2
        err = mean(abs.(pl .- model))
        ap = knee ? ap : [ap[1], 0.0, ap[2]]
        return (ap = ap, pe = pe, gp = gp, r2 = r2, err = err, f = f, p = pl, model = model, ap_fit = ap_fit)
    catch e
        verbose && _warn("FOOOF model fitting failed: $(sprint(showerror, e))")
        return _fooof_nan(f, p)
    end

end

"""
    fooof(s; <keyword arguments>)

Parameterize power spectra into aperiodic and periodic (peak) components using the FOOOF / specparam algorithm (Donoghue et al., 2020). PSD is calculated first (default: Welch's periodogram).

# Arguments

  - `s::AbstractArray`: channels × samples × epochs
  - `fs::Int64`: sampling rate
  - `flim::Tuple{Real, Real}=(1, fs / 2)`: frequency range to fit
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
      + `:gh`: Gaussian and Hilbert transform
  - `nt::Int64=7`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet
  - `gw::Real=5`: Gaussian width in Hz
  - `peak_width_limits::Tuple{Real, Real}=(0.5, 12.0)`: limits on peak bandwidth (in Hz)
  - `max_n_peaks::Int64=typemax(Int64)`: maximum number of peaks
  - `min_peak_height::Real=0.0`: absolute threshold for peak height (log10 power)
  - `peak_threshold::Real=2.0`: relative threshold for peak detection (in SD)
  - `aperiodic_mode::Symbol=:fixed`: aperiodic model (`:fixed` or `:knee`)

# Returns

Named tuple containing:

  - `ap::Array{Float64, 3}`: aperiodic parameters (channels × [offset, knee, exponent] × epochs)
  - `pe::Matrix{Matrix{Float64}}`: peak parameters (channels × epochs), each a peaks × [CF, PW, BW] matrix
  - `r2::Matrix{Float64}`: R² (channels × epochs)
  - `err::Matrix{Float64}`: mean absolute error (channels × epochs)
  - `f::Vector{Float64}`: frequencies of the fit
  - `model::Array{Float64, 3}`: full models (log10, channels × frequencies × epochs)
  - `ap_fit::Array{Float64, 3}`: aperiodic components (log10, channels × frequencies × epochs)
"""
function fooof(
        s::AbstractArray;
        fs::Int64,
        flim::Tuple{Real, Real} = (1, fs / 2),
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        gw::Real = 5,
        peak_width_limits::Tuple{Real, Real} = (0.5, 12.0),
        max_n_peaks::Int64 = typemax(Int64),
        min_peak_height::Real = 0.0,
        peak_threshold::Real = 2.0,
        aperiodic_mode::Symbol = :fixed,
    )::@NamedTuple{
        ap::Array{Float64, 3},
        pe::Matrix{Matrix{Float64}},
        r2::Matrix{Float64},
        err::Matrix{Float64},
        f::Vector{Float64},
        model::Array{Float64, 3},
        ap_fit::Array{Float64, 3},
    }

    _check_var(method, [:welch, :fft, :mt, :stft, :mw, :gh], "method")
    _chk3d(s)
    _check_tuple(flim, (0, fs / 2), "flim")
    @assert flim[1] > 0 "flim[1] must be > 0 (log10 of 0 Hz is undefined)."
    _check_var(aperiodic_mode, [:fixed, :knee], "aperiodic_mode")
    _check_tuple(peak_width_limits, (0, Inf), "peak_width_limits")
    @assert peak_width_limits[1] > 0 "peak_width_limits[1] must be > 0."
    @assert max_n_peaks >= 0 "max_n_peaks must be ≥ 0."
    @assert peak_threshold > 0 "peak_threshold must be > 0."
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    pw, pf = psd(
        s; fs = fs, db = false, method = method, nt = nt, wlen = wlen, woverlap = woverlap, w = w, ncyc = ncyc, gw = gw
    )

    f_idx = (pf .>= flim[1]) .& (pf .<= flim[2])
    f = Float64.(pf[f_idx])
    _fooof_check(f, ones(length(f)))
    _fooof_width_warn(f, peak_width_limits)

    kw = (
        peak_width_limits = peak_width_limits,
        max_n_peaks = max_n_peaks,
        min_peak_height = min_peak_height,
        peak_threshold = peak_threshold,
        aperiodic_mode = aperiodic_mode,
        verbose = false,
    )

    ap = zeros(ch_n, 3, ep_n)
    pe = Matrix{Matrix{Float64}}(undef, ch_n, ep_n)
    r2 = zeros(ch_n, ep_n)
    err = zeros(ch_n, ep_n)
    model = zeros(ch_n, length(f), ep_n)
    ap_fit = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            p = Float64.(pw[ch_idx, f_idx, ep_idx])
            # invalid spectra (e.g. flat / zero-power channels) are not fitted, results are NaN
            r = all(isfinite, p) && all(p .> 0) ? _fooof(f, p; kw...) : _fooof_nan(f, p)
            ap[ch_idx, :, ep_idx] = r.ap
            pe[ch_idx, ep_idx] = r.pe
            r2[ch_idx, ep_idx] = r.r2
            err[ch_idx, ep_idx] = r.err
            model[ch_idx, :, ep_idx] = r.model
            ap_fit[ch_idx, :, ep_idx] = r.ap_fit
        end
    end

    n_fail = count(isnan, r2)
    n_fail > 0 && _warn("FOOOF model fitting failed for $n_fail of $(length(r2)) spectra (parameters set to NaN).")

    return (ap = ap, pe = pe, r2 = r2, err = err, f = f, model = model, ap_fit = ap_fit)

end

"""
    fooof(obj; <keyword arguments>)

Parameterize power spectra into aperiodic and periodic (peak) components using the FOOOF / specparam algorithm (Donoghue et al., 2020). PSD is calculated first (default: Welch's periodogram).

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `flim::Tuple{Real, Real}=(1, sr(obj) / 2)`: frequency range to fit
  - `method::Symbol=:welch`: method used to calculate PSD (`:welch`, `:fft`, `:mt`, `:stft`, `:mw`, `:gh`)
  - `nt::Int64=7`: number of Slepian tapers
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet
  - `gw::Real=5`: Gaussian width in Hz
  - `peak_width_limits::Tuple{Real, Real}=(0.5, 12.0)`: limits on peak bandwidth (in Hz)
  - `max_n_peaks::Int64=typemax(Int64)`: maximum number of peaks
  - `min_peak_height::Real=0.0`: absolute threshold for peak height (log10 power)
  - `peak_threshold::Real=2.0`: relative threshold for peak detection (in SD)
  - `aperiodic_mode::Symbol=:fixed`: aperiodic model (`:fixed` or `:knee`)

# Returns

Named tuple containing:

  - `ap::Array{Float64, 3}`: aperiodic parameters (channels × [offset, knee, exponent] × epochs)
  - `pe::Matrix{Matrix{Float64}}`: peak parameters (channels × epochs), each a peaks × [CF, PW, BW] matrix
  - `r2::Matrix{Float64}`: R² (channels × epochs)
  - `err::Matrix{Float64}`: mean absolute error (channels × epochs)
  - `f::Vector{Float64}`: frequencies of the fit
  - `model::Array{Float64, 3}`: full models (log10, channels × frequencies × epochs)
  - `ap_fit::Array{Float64, 3}`: aperiodic components (log10, channels × frequencies × epochs)
"""
function fooof(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        flim::Tuple{Real, Real} = (1, sr(obj) / 2),
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = sr(obj),
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        gw::Real = 5,
        peak_width_limits::Tuple{Real, Real} = (0.5, 12.0),
        max_n_peaks::Int64 = typemax(Int64),
        min_peak_height::Real = 0.0,
        peak_threshold::Real = 2.0,
        aperiodic_mode::Symbol = :fixed,
    )::@NamedTuple{
        ap::Array{Float64, 3},
        pe::Matrix{Matrix{Float64}},
        r2::Matrix{Float64},
        err::Matrix{Float64},
        f::Vector{Float64},
        model::Array{Float64, 3},
        ap_fit::Array{Float64, 3},
    }

    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    _log_off()
    r = fooof(
        obj.data[ch, :, :];
        fs = sr(obj),
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        peak_width_limits = peak_width_limits,
        max_n_peaks = max_n_peaks,
        min_peak_height = min_peak_height,
        peak_threshold = peak_threshold,
        aperiodic_mode = aperiodic_mode,
    )
    _log_on()

    return r

end
