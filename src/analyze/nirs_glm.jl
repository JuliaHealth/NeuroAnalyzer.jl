export nirs_hrf
export nirs_glm

"""
    nirs_hrf(; <keyword arguments>)

Canonical hemodynamic response function (HRF): a difference of two gamma functions (the SPM canonical HRF) with a peak at about 5 s and an undershoot at about 15 s:

`h(t) = t^(a1 − 1) e^(−t) / Γ(a1) − c · t^(a2 − 1) e^(−t) / Γ(a2)`, with `a1 = 6`, `a2 = 16`, `c = 1/6`

The HRF is normalized to unit sum, so that a stimulus block longer than the response gives a plateau of 1 in the convolved regressor.

# Arguments

- `fs::Int64`: sampling rate in Hz
- `duration::Real=32`: length of the HRF in seconds

# Returns

Named tuple:

- `t::Vector{Float64}`: time points in seconds
- `h::Vector{Float64}`: HRF
"""
function nirs_hrf(;
    fs::Int64,
    duration::Real = 32,
)::@NamedTuple{
    t::Vector{Float64},
    h::Vector{Float64},
}
    # validate
    fs > 0 || throw(ArgumentError("fs must be > 0."))
    duration > 0 || throw(ArgumentError("duration must be > 0."))

    t = collect(0:(1 / fs):duration)
    # gamma density with integer shape a: Γ(a) = (a - 1)!
    g(a) = t .^ (a - 1) .* exp.(-t) ./ factorial(a - 1)
    h = g(6) .- g(16) ./ 6
    h ./= sum(h)

    return (; t, h)
end

# discrete cosine drift basis (SPM): cosines with periods longer than 1 / cutoff
function _dct_drift(n::Int64, fs::Int64, cutoff::Real)::Matrix{Float64}
    k = floor(Int64, 2 * (n / fs) * cutoff + 1)
    k <= 1 && return zeros(n, 0)
    t = 0:(n - 1)
    return reduce(hcat, [cos.(pi .* (2 .* t .+ 1) .* j ./ (2 * n)) for j = 1:(k - 1)])
end

# AR(p) coefficients of x (Yule-Walker, Levinson-Durbin) and the innovation variance
function _ar_yw(x::AbstractVector, p::Int64)::Tuple{Vector{Float64}, Float64}
    n = length(x)
    r = [sum(x[1:(n - k)] .* x[(1 + k):n]) / n for k = 0:p]
    a = zeros(p)
    e = r[1]
    for k = 1:p
        acc = r[k + 1] - sum(a[j] * r[k - j + 1] for j = 1:(k - 1); init = 0.0)
        kk = acc / e
        a_new = copy(a)
        a_new[k] = kk
        for j = 1:(k - 1)
            a_new[j] = a[j] - kk * a[k - j]
        end
        a = a_new
        e *= (1 - kk^2)
    end

    return a, e
end

"""
    nirs_glm(obj; <keyword arguments>)

Estimate the hemodynamic activation of NIRS channels using the general linear model (GLM).

For each condition (marker value), a boxcar function of the stimulation blocks (marker onsets, block durations from the marker lengths or `duration`) is convolved with the canonical HRF (`nirs_hrf()`). The design matrix contains these regressors, a constant and a discrete cosine basis modeling slow drifts (periods longer than 1 / `drift` seconds). The model is fitted to each channel by ordinary least squares; the serially correlated noise of NIRS signals is accounted for by prewhitening with an autoregressive model of the residuals (AR(p), order selected by BIC up to `ar`, as in the AR-IRLS method of Barker et al. 2013, without the robust weighting).

The object must be continuous (one epoch) and contain markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj; type = ["nirs_hbo", "nirs_hbr"])`: channel name(s), default is all HbO and HbR channels
- `cond::Union{Nothing, String, Vector{String}}=nothing`: conditions (marker values); `nothing`: all marker values
- `duration::Union{Nothing, Real}=nothing`: block duration in seconds; `nothing`: the marker lengths
- `drift::Real=0.01`: cut-off frequency of the drift model in Hz (`0`: no drift regressors)
- `ar::Int64=2 * sr(obj)`: maximum order of the AR prewhitening model (`0`: no prewhitening)

# Returns

- `DataFrame`: columns `channel`, `condition`, `beta` (in the units of the channel, for a block plateau), `se`, `t`, `df`, `p`

# References

1. Barker JW, Aarabi A, Huppert TJ. Autoregressive model based algorithm for correcting motion and serially correlated errors in fNIRS. Biomedical Optics Express. 2013;4(8):1366–1379.
2. Friston KJ, et al. Statistical parametric maps in functional imaging: a general linear approach. Human Brain Mapping. 1994;2(4):189–210.
"""
function nirs_glm(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = get_channel(obj; type = ["nirs_hbo", "nirs_hbr"]),
    cond::Union{Nothing, String, Vector{String}} = nothing,
    duration::Union{Nothing, Real} = nothing,
    drift::Real = 0.01,
    ar::Int64 = 2 * sr(obj),
)::DataFrame
    # validate
    _check_datatype(obj, "nirs")
    nepochs(obj) == 1 || throw(ArgumentError("OBJ must be continuous (one epoch)."))
    _has_markers(obj) || throw(ArgumentError("OBJ does not contain markers."))
    isnothing(duration) || duration > 0 || throw(ArgumentError("duration must be > 0."))
    drift >= 0 || throw(ArgumentError("drift must be ≥ 0."))
    ar >= 0 || throw(ArgumentError("ar must be ≥ 0."))
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    fs = sr(obj)
    n = epoch_len(obj)
    mrk = obj.markers
    conds = isnothing(cond) ? unique(string.(mrk[!, :value])) :
            (cond isa String ? [cond] : cond)
    for c in conds
        c in string.(mrk[!, :value]) || throw(ArgumentError("No markers with value $c."))
    end

    # task regressors: boxcars convolved with the canonical HRF
    h = nirs_hrf(; fs = fs).h
    X_task = zeros(n, length(conds))
    for (k, c) in enumerate(conds)
        u = zeros(n)
        for r in eachrow(mrk)
            string(r[:value]) == c || continue
            dur = isnothing(duration) ? r[:length] : duration
            dur > 0 || throw(
                ArgumentError("Marker $c has length 0: block duration must be given."),
            )
            i1 = clamp(round(Int64, r[:start] * fs) + 1, 1, n)
            i2 = clamp(round(Int64, (r[:start] + dur) * fs), 1, n)
            u[i1:i2] .= 1
        end
        X_task[:, k] = DSP.conv(u, h)[1:n]
    end
    X = hcat(X_task, ones(n), _dct_drift(n, fs, drift))
    k = size(X, 2)

    df_out = DataFrame(;
        channel = String[], condition = String[], beta = Float64[], se = Float64[],
        t = Float64[], df = Int64[], p = Float64[],
    )
    for c in ch
        y = obj.data[c, :, 1]
        b = X \ y
        res = y .- X * b
        Xw, yw = X, y
        p_ar = 0
        if ar > 0
            # AR order by BIC
            best = (Inf, 0, Float64[])
            for p = 1:min(ar, n ÷ 4)
                a, e = _ar_yw(res .- mean(res), p)
                e > 0 || break
                bic = n * log(e) + p * log(n)
                bic < best[1] && (best = (bic, p, a))
            end
            p_ar, a = best[2], best[3]
            if p_ar > 0
                # whitening filter: x[t] - Σ a[j] x[t - j]
                wf(z) = [z[i] - sum(a[j] * z[i - j] for j = 1:p_ar) for i = (p_ar + 1):n]
                yw = wf(y)
                Xw = reduce(hcat, [wf(X[:, j]) for j = 1:k])
                b = Xw \ yw
            end
        end
        res_w = yw .- Xw * b
        dfr = length(yw) - k
        sigma2 = sum(abs2, res_w) / dfr
        cov_b = sigma2 .* inv(Xw' * Xw)
        for (j, cn) in enumerate(conds)
            se = sqrt(cov_b[j, j])
            tv = b[j] / se
            pv = 2 * ccdf(TDist(dfr), abs(tv))
            push!(df_out, (labels(obj)[c], cn, b[j], se, tv, dfr, pv))
        end
    end

    return df_out
end
