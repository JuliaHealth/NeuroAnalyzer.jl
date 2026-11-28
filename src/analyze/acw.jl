export acw

"""
    _acw(r, lags, level)

Return the first lag at which an auto-correlation function drops to or below `level`, linearly interpolated between the two neighbouring lags.

# Arguments

- `r::AbstractVector`: auto-correlation function for lags ≥ 0 (`r[1]` is lag 0)
- `lags::AbstractVector`: lags (same length as `r`)
- `level::Real`: level, e.g. `0.5` (ACW-50), `0` (ACW-0) or `1/e` (ACW-e)

# Returns

- `Float64`: interpolated lag (in the units of `lags`); `NaN` if the ACF does not reach `level`
"""
function _acw(r::AbstractVector, lags::AbstractVector, level::Real)::Float64
    length(r) == length(lags) || throw(ArgumentError("r and lags must have the same length."))
    i = findfirst(<=(level), r)
    isnothing(i) && return NaN
    i == 1 && return Float64(lags[1])
    # interpolate between lags i - 1 and i
    dl = lags[i] - lags[i - 1]
    return lags[i - 1] + (r[i - 1] - level) / (r[i - 1] - r[i]) * dl
end

"""
    _acw_auc(r, lags)

Return the area under an auto-correlation function from lag 0 to its first zero crossing (ACW-0), using the trapezoidal rule; the last segment ends at the interpolated zero crossing.

# Arguments

- `r::AbstractVector`: auto-correlation function for lags ≥ 0 (`r[1]` is lag 0)
- `lags::AbstractVector`: lags (same length as `r`)

# Returns

- `Float64`: area (in the units of `lags`); `NaN` if the ACF does not reach 0
"""
function _acw_auc(r::AbstractVector, lags::AbstractVector)::Float64
    i = findfirst(<=(0), r)
    (isnothing(i) || i == 1) && return NaN
    auc = 0.0
    for k in 2:(i - 1)
        auc += (r[k - 1] + r[k]) / 2 * (lags[k] - lags[k - 1])
    end
    # last segment: from lag i - 1 to the zero crossing (a triangle)
    auc += r[i - 1] / 2 * (_acw(r, lags, 0) - lags[i - 1])
    return auc
end

"""
    acw(obj; <keyword arguments>)

Calculate auto-correlation windows (ACW), measures of intrinsic neural timescales: the lags at which the auto-correlation function (ACF, see `acor()`) first drops to a given level, linearly interpolated between neighbouring lags.

- ACW-50: the lag at which the ACF first drops to 0.5
- ACW-0: the lag at which the ACF first reaches 0
- ACW-e: the lag at which the ACF first drops to 1/e
- AUC: the area under the ACF from lag 0 to ACW-0

Measures that are not reached within the maximum lag `l` are `NaN`. For ERP the first epoch is the trial-averaged ACF (see `acor()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `l::Real=1.0`: maximum lag in seconds (limited to the epoch length)
- `demean::Bool=true`: subtract the mean before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance
- `method::Symbol=:sum`: method of calculating auto-correlation (see `acor()`)

# Returns

Named tuple (all values in seconds, shape (channels, epochs)):

- `acw50::Matrix{Float64}`: ACW-50
- `acw0::Matrix{Float64}`: ACW-0
- `acwe::Matrix{Float64}`: ACW-e
- `auc::Matrix{Float64}`: area under the ACF up to ACW-0

# References

1. Honey CJ, et al. Slow cortical dynamics and the accumulation of information over long timescales. Neuron. 2012;76(2):423-434.
2. Golesorkhi M, et al. The brain and its time: intrinsic neural timescales are key for input processing. Commun Biol. 2021;4:970.
"""
function acw(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    l::Real = 1.0,
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)::@NamedTuple{
    acw50::Matrix{Float64},
    acw0::Matrix{Float64},
    acwe::Matrix{Float64},
    auc::Matrix{Float64},
}
    # validate
    l > 0 || throw(ArgumentError("l must be > 0."))

    # maximum lag in samples, limited to the epoch length
    l_smp = min(round(Int64, l * sr(obj)), epoch_len(obj) - 1)

    ac = acor(
        obj;
        ch = ch,
        l = l_smp,
        demean = demean,
        biased = biased,
        method = method,
        exclude_bad_channels = exclude_bad_channels,
    )

    # lags ≥ 0
    z = l_smp + 1
    lags = ac.lags[z:end]
    r = @view ac.autocor[:, z:end, :]

    ch_n = size(r, 1)
    ep_n = size(r, 3)
    acw50 = zeros(ch_n, ep_n)
    acw0 = zeros(ch_n, ep_n)
    acwe = zeros(ch_n, ep_n)
    auc = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n, ch_idx in 1:ch_n
        rr = @view r[ch_idx, :, ep_idx]
        acw50[ch_idx, ep_idx] = _acw(rr, lags, 0.5)
        acw0[ch_idx, ep_idx] = _acw(rr, lags, 0)
        acwe[ch_idx, ep_idx] = _acw(rr, lags, exp(-1))
        auc[ch_idx, ep_idx] = _acw_auc(rr, lags)
    end

    return (; acw50, acw0, acwe, auc)
end
