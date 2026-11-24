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

# AR(p) models of x for all orders p = 1:pmax (Yule-Walker, Levinson-Durbin): the coefficients
# and the innovation variance of each order; the recursion of order p contains all lower orders,
# so the autocorrelation is calculated once
function _ar_yw(
    x::AbstractVector,
    pmax::Int64,
)::Tuple{Vector{Vector{Float64}}, Vector{Float64}}
    n = length(x)
    r = [dot(view(x, 1:(n - k)), view(x, (1 + k):n)) / n for k = 0:pmax]
    as = Vector{Vector{Float64}}(undef, pmax)
    es = zeros(pmax)
    a = Float64[]
    e = r[1]
    for k = 1:pmax
        acc = r[k + 1] - sum(a[j] * r[k - j + 1] for j = 1:(k - 1); init = 0.0)
        kk = acc / e
        a_new = zeros(k)
        a_new[k] = kk
        for j = 1:(k - 1)
            a_new[j] = a[j] - kk * a[k - j]
        end
        a = a_new
        e *= (1 - kk^2)
        as[k] = a
        es[k] = e
    end

    return as, es
end

"""
    nirs_glm(obj; <keyword arguments>)

Estimate the hemodynamic activation of NIRS channels using the general linear model (GLM).

For each condition (by default a marker value; see `condition`), a boxcar function of the stimulation blocks (marker onsets, block durations from the marker lengths or `duration`) is convolved with the canonical HRF (`nirs_hrf()`). The design matrix contains these regressors, a constant and a discrete cosine basis modeling slow drifts (periods longer than 1 / `drift` seconds). The model is fitted to each channel by ordinary least squares; the serially correlated noise of NIRS signals is accounted for by prewhitening with an autoregressive model of the residuals (AR(p), order selected by BIC up to `ar`, as in the AR-IRLS method of Barker et al. 2013, without the robust weighting).

The object must be continuous (one epoch) and contain markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj; type = ["nirs_hbo", "nirs_hbr"])`: channel name(s), default is all HbO and HbR channels
- `cond::Union{Nothing, String, Vector{String}}=nothing`: conditions (condition labels, by default marker values); `nothing`: all conditions
- `marker_type::String=""`: if not empty, only markers of this type are used (e.g. `"stim"`)
- `condition::Union{String, AbstractDict}="value"`: condition label of each marker: the name of a marker column (`"value"`, `"type"` or an extra column, e.g. `"trial_type"` of BIDS events) or a dictionary of marker values and condition labels (e.g. `Dict("1" => "task", "2" => "task")`: one regressor for both markers; values not in the dictionary are their own labels)
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
    marker_type::String = "",
    condition::Union{String, AbstractDict} = "value",
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
    # condition label of each marker (markers of other types and without a label are not used)
    m_cond = _marker_conditions(mrk, condition)
    m_use = _marker_select(mrk, "", marker_type) .& (m_cond .!= "")
    avail = unique(m_cond[m_use])
    conds = isnothing(cond) ? avail : (cond isa String ? [cond] : cond)
    isempty(conds) && throw(ArgumentError("No markers to model."))
    for c in conds
        c in avail || throw(ArgumentError("No markers of condition $c."))
    end

    # task regressors: boxcars convolved with the canonical HRF
    h = nirs_hrf(; fs = fs).h
    X_task = zeros(n, length(conds))
    for (k, c) in enumerate(conds)
        u = zeros(n)
        for (i, r) in enumerate(eachrow(mrk))
            (m_use[i] && m_cond[i] == c) || continue
            dur = isnothing(duration) ? r[:length] : duration
            dur > 0 || throw(
                ArgumentError("Markers of condition $c have length 0: block duration must be given."),
            )
            i1 = clamp(round(Int64, r[:start] * fs) + 1, 1, n)
            i2 = clamp(round(Int64, (r[:start] + dur) * fs), 1, n)
            u[i1:i2] .= 1
        end
        X_task[:, k] = DSP.conv(u, h)[1:n]
    end
    X = hcat(X_task, ones(n), _dct_drift(n, fs, drift))
    k = size(X, 2)

    # ordinary least squares of all channels at once (the design is the same for all channels)
    B_ols = X \ Matrix(obj.data[ch, :, 1]')

    # one row per channel and condition: (beta, se, t, df, p)
    rows = Vector{Vector{Tuple{Float64, Float64, Float64, Int64, Float64}}}(undef, length(ch))
    Threads.@threads :static for ch_idx in eachindex(ch)
        y = obj.data[ch[ch_idx], :, 1]
        b = B_ols[:, ch_idx]
        res = y .- X * b
        Xw, yw = X, y
        p_ar = 0
        if ar > 0
            # AR order by BIC (models of all orders from one Levinson-Durbin recursion)
            pmax = min(ar, n ÷ 4)
            as, es = _ar_yw(res .- mean(res), pmax)
            best = (Inf, 0)
            for p = 1:pmax
                es[p] > 0 || break
                bic = n * log(es[p]) + p * log(n)
                bic < best[1] && (best = (bic, p))
            end
            p_ar = best[2]
            if p_ar > 0
                # whitening filter: x[t] - Σ a[j] x[t - j], for t > p_ar
                wf = vcat(1.0, -as[p_ar])
                yw = DSP.filt(wf, y)[(p_ar + 1):n]
                Xw = DSP.filt(wf, X)[(p_ar + 1):n, :]
            end
        end
        # normal equations: XᵀX is needed for the covariance of the estimates anyway
        XtX = Symmetric(Xw' * Xw)
        p_ar > 0 && (b = XtX \ (Xw' * yw))
        res_w = yw .- Xw * b
        dfr = length(yw) - k
        sigma2 = sum(abs2, res_w) / dfr
        cov_b = sigma2 .* inv(XtX)
        rows[ch_idx] = map(eachindex(conds)) do j
            se = sqrt(cov_b[j, j])
            tv = b[j] / se
            (b[j], se, tv, dfr, 2 * ccdf(TDist(dfr), abs(tv)))
        end
    end

    df_out = DataFrame(;
        channel = String[], condition = String[], beta = Float64[], se = Float64[],
        t = Float64[], df = Int64[], p = Float64[],
    )
    for (ch_idx, c) in enumerate(ch), (j, cn) in enumerate(conds)
        push!(df_out, (labels(obj)[c], cn, rows[ch_idx][j]...))
    end

    return df_out
end
