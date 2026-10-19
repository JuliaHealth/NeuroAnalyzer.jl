export filter_mmed
export filter_mmed!

"""
    _filter_mmed(s; <keyword arguments>)

Filter a signal using a weighted moving median with an optional threshold.

Samples within the threshold band `[mean(s) − t×std(s), mean(s) + t×std(s)]` are left unchanged; only outliers outside that band are replaced by the local weighted median. When `t = 0` (default) every sample is filtered.

# Arguments

- `s::AbstractVector`: signal vector
- `k::Int64=8`: half-window length; full window is `N = 2k + 1` samples; for a desired normalized cutoff `F = f/fs` of the (unweighted) moving average, `N ≈ sqrt(0.196202 + F^2) / F`, i.e. `k = round(Int, (N - 1) / 2)`; must satisfy `1 ≤ k < length(s)`
- `t::Real=0`: threshold multiplier (≥ 0). `t = 0` filters all samples
- `ww::AbstractVector=ones(2k+1)`: weighting window of length `2k + 1`; weights must be ≥ 0 (weighted median)

# Returns

- `Vector{Float64}`: filtered signal of the same length as `s`

# References

1. https://dsp.stackexchange.com/questions/9966/what-is-the-cutoff-frequency-of-a-moving-average-filter
"""
function _filter_mmed(
    s::AbstractVector;
    k::Int64 = 8,
    t::Real = 0,
    ww::AbstractVector = ones(2 * k + 1),
)::Vector{Float64}
    # check k
    _in(k, (1, length(s) - 1), "k")
    # check weighting window
    length(ww) == 2 * k + 1 ||
        throw(ArgumentError("length(ww) must be 2k + 1 ($(2k + 1))."))

    all(>=(0), ww) && any(>(0), ww) ||
        throw(ArgumentError("ww must be non-negative with at least one positive weight."))

    # cache threshold bounds once
    s_mean = mean(s)
    s_std = std(s)
    lo = s_mean - t * s_std
    hi = s_mean + t * s_std

    # unfiltered samples stay at their original value
    s_filtered = Vector{Float64}(s)

    # centered window [idx-k … idx+k], truncated at the signal edges (with the matching
    # part of the weighting window)
    n = length(s)
    @inbounds for idx = 1:n
        (t <= 0 || s[idx] < lo || s[idx] > hi) || continue
        i1 = max(1, idx - k)
        i2 = min(n, idx + k)
        s_tmp = @view s[i1:i2]
        w_tmp = @view ww[(i1 - idx + k + 1):(i2 - idx + k + 1)]
        s_filtered[idx] = _wmedian_cols(reshape(Vector{Float64}(s_tmp), :, 1), Vector{Float64}(w_tmp))[1]
    end

    return s_filtered
end

"""
    _filter_mmed(s; <keyword arguments>)

Apply a weighted moving median filter to every channel × epoch slice of a 3-D signal array.

# Arguments

- `s::AbstractArray`: 3-D signal array, shape (channels, samples, epochs)
- `k::Int64=8`: half-window length; full window is `N = 2k + 1` samples; for a desired normalized cutoff `F = f/fs` of the (unweighted) moving average, `N ≈ sqrt(0.196202 + F^2) / F`, i.e. `k = round(Int, (N - 1) / 2)`; must satisfy `1 ≤ k < length(s)`
- `t::Real=0`: threshold multiplier (≥ 0). `t = 0` filters all samples
- `ww::AbstractVector=ones(2k+1)`: weighting window of length `2k + 1`; weights must be ≥ 0 (weighted median)

# Returns

- `Array{Float64, 3}`: filtered array of the same shape as `s`
"""
function _filter_mmed(
    s::AbstractArray;
    k::Int64 = 8,
    t::Real = 0,
    ww::AbstractVector = ones(2 * k + 1),
)::Array{Float64, 3}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_filtered = similar(s, Float64)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_filtered[ch_idx, :, ep_idx] =
            _filter_mmed(@view(s[ch_idx, :, ep_idx]), k = k, t = t, ww = ww)
    end

    return s_filtered
end

"""
    filter_mmed(obj; <keyword arguments>)

Apply a weighted moving median filter to selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `k::Int64=8`: half-window length; full window is `N = 2k + 1` samples; for a desired normalized cutoff `F = f/fs` of the (unweighted) moving average, `N ≈ sqrt(0.196202 + F^2) / F`, i.e. `k = round(Int, (N - 1) / 2)`; must satisfy `1 ≤ k < length(s)`
- `t::Real=0`: threshold multiplier (≥ 0). `t = 0` filters all samples
- `ww::AbstractVector=ones(2k+1)`: weighting window of length `2k + 1`; weights must be ≥ 0 (weighted median)

# Returns

- `NeuroAnalyzer.NEURO`: new object with filtered channels
"""
function filter_mmed(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    k::Int64 = 8,
    t::Real = 0,
    ww::AbstractVector = ones(2 * k + 1),
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # window length
    wlen = 2 * k + 1
    _info("Window length: $wlen samples")

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.data[ch, :, :] = _filter_mmed(
        @view(obj.data[ch, :, :]); k = k, t = t, ww = ww,
    )
    push!(obj_tmp.history, "filter_mmed(obj; ch=$ch, k=$k, t=$t, ww=$ww)")

    return obj_tmp
end

"""
    filter_mmed!(obj; <keyword arguments>)

Apply a weighted moving median filter in-place to selected channels of a NEURO object

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `k::Int64=8`: half-window length; full window is `N = 2k + 1` samples; for a desired normalized cutoff `F = f/fs` of the (unweighted) moving average, `N ≈ sqrt(0.196202 + F^2) / F`, i.e. `k = round(Int, (N - 1) / 2)`; must satisfy `1 ≤ k < length(s)`
- `t::Real=0`: threshold multiplier (≥ 0). `t = 0` filters all samples
- `ww::AbstractVector=ones(2k+1)`: weighting window of length `2k + 1`; weights must be ≥ 0 (weighted median)

# Returns

- `Nothing`
"""
function filter_mmed!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    k::Int64 = 8,
    t::Real = 0,
    ww::AbstractVector = ones(2 * k + 1),
)::Nothing
    obj_tmp = filter_mmed(obj; ch = ch, k = k, t = t, ww = ww)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
