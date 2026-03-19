export filter_mmed
export filter_mmed!

"""
    filter_mmed(s; <keyword arguments>)

Filter a signal using a weighted moving median with an optional threshold.

Samples within the threshold band `[mean(s) − t×std(s), mean(s) + t×std(s)]` are left unchanged; only outliers outside that band are replaced by the local weighted median. When `t = 0` (default) every sample is filtered.

# Arguments

- `s::AbstractVector`: signal vector
- `k::Int64=8`: half-window length; full window is `2k + 1` samples; for a desired normalized cutoff `F = f/fs`, choose
  `k = round(Int, sqrt(0.196202 + F^2) / F)`; must satisfy `1 ≤ k < length(s)`
- `t::Real=0`: threshold multiplier (≥ 0). `t = 0` filters all samples
- `ww::Union{Nothing, AbstractVector}=ones(2k+1)`: weighting window of length `2k + 1`

# Returns

- `Vector{Float64}`: filtered signal of the same length as `s`

# Throws

- `ArgumentError`: if `k` is out of range or `length(ww) ≠ 2k + 1`

# References

1. https://dsp.stackexchange.com/questions/9966/what-is-the-cutoff-frequency-of-a-moving-average-filter

# See also

[`filter_mmed(::AbstractArray)`](@ref), [`filter_mmed(::NeuroAnalyzer.NEURO)`](@ref)
"""
function filter_mmed(
    s::AbstractVector;
    k::Int64 = 8,
    t::Real = 0,
    ww::AbstractVector = ones(2 * k + 1)
)::Vector{Float64}

    # check k
    _in(k, (1, length(s) - 1), "k")
    # check weighting window
    length(ww) == 2 * k + 1 || throw(ArgumentError("length(ww) must be 2k + 1 ($(2k + 1))."))

    # cache threshold bounds once
    s_mean = mean(s)
    s_std  = std(s)
    lo = s_mean - t * s_std
    hi = s_mean + t * s_std

    # helper: return true when sample x is outside the threshold band
    needs_filter(x) = t <= 0 || x < lo || x > hi

    # unfiltered samples stay at their original value
    s_filtered = copy(s)

    # left edge: truncated window [1 … idx]
    @inbounds for idx in 1:k
        needs_filter(s[idx]) || continue
        s_tmp = @view s[1:idx]
        w_tmp = @view ww[1:idx]
        s_filtered[idx] = median(s_tmp .* w_tmp)
    end

    # interior: full centered window [idx-k … idx+k]
    @inbounds for idx in (k + 1):(length(s) - k)
        needs_filter(s[idx]) || continue
        s_tmp = @view s[(idx - k):(idx + k)]
        s_filtered[idx] = median(s_tmp .* ww)
    end

    # right edge: truncated window [idx … end]
    @inbounds for idx in (length(s) - k + 1):length(s)
        needs_filter(s[idx]) || continue
        s_tmp = @view s[idx:end]
        w_tmp = @view ww[(end - length(s_tmp) + 1):end]
        s_filtered[idx] = median(s_tmp .* w_tmp)
    end

    return s_filtered

end

"""
    filter_mmed(s; <keyword arguments>)

Apply a weighted moving median filter to every channel × epoch slice of a 3-D signal array. Delegates to [`filter_mavg(::AbstractVector)`](@ref).

# Arguments

- `s::AbstractArray`: 3-D signal array, shape `(channels, samples, epochs)`
- `k::Int64=8`: half-window length; full window is `2k + 1` samples; for a desired normalized cutoff `F = f/fs`, choose
  `k = round(Int, sqrt(0.196202 + F^2) / F)`; must satisfy `1 ≤ k < length(s)`
- `t::Real=0`: threshold multiplier (≥ 0). `t = 0` filters all samples
- `ww::Union{Nothing, AbstractVector}=ones(2k+1)`: weighting window of length `2k + 1`

# Returns

- `Array{Float64, 3}`: filtered array of the same shape as `s`

# See also

[`filter_mmed(::AbstractVector)`](@ref), [`filter_mmed(::NeuroAnalyzer.NEURO)`](@ref)
"""
function filter_mmed(
    s::AbstractArray;
    k::Int64 = 8,
    t::Real = 0, ww::AbstractVector = ones(2 * k + 1)
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
        s_filtered[ch_idx, :, ep_idx] = filter_mmed(@view(s[ch_idx, :, ep_idx]), k=k, t=t, ww=ww)
    end

    return s_filtered

end

"""
    filter_mmed(obj; <keyword arguments>)

Apply a weighted moving median filter to selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `k::Int64=8`: half-window length; full window is `2k + 1` samples; for a desired normalized cutoff `F = f/fs`, choose
  `k = round(Int, sqrt(0.196202 + F^2) / F)`; must satisfy `1 ≤ k < length(s)`
- `t::Real=0`: threshold multiplier (≥ 0). `t = 0` filters all samples
- `ww::Union{Nothing, AbstractVector}=ones(2k+1)`: weighting window of length `2k + 1`

# Returns

- `NeuroAnalyzer.NEURO`: new object with filtered channels

# See also

[`filter_mmed!`](@ref), [`filter_mmed(::AbstractArray)`](@ref)
"""
function filter_mmed(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    k::Int64 = 8,
    t::Real = 0,
    ww::AbstractVector = ones(2 * k + 1),
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch=ch)

    wlen = 2 * k + 1
    _info("Window length: $wlen samples")

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = filter_mmed(
        @view(obj.data[ch, :, :]), k=k, t=t, ww=ww
    )
    push!(obj_new.history, "filter_mmed(OBJ, ch=$ch, k=$k, t=$t, ww=$ww)")
    return obj_new

end

"""
    filter_mmed!(obj; <keyword arguments>)

Apply a weighted moving median filter in-place to selected channels of a NEURO object. Delegates to [`filter_mavg`](@ref) and copies the result back.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `k::Int64=8`: half-window length; full window is `2k + 1` samples; for a desired normalized cutoff `F = f/fs`, choose
  `k = round(Int, sqrt(0.196202 + F^2) / F)`; must satisfy `1 ≤ k < length(s)`
- `t::Real=0`: threshold multiplier (≥ 0). `t = 0` filters all samples
- `ww::Union{Nothing, AbstractVector}=ones(2k+1)`: weighting window of length `2k + 1`

# Returns

- `Nothing`

# See also

[`filter_mmed`](@ref)
"""
function filter_mmed!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    k::Int64 = 8,
    t::Real = 0,
    ww::AbstractVector = ones(2 * k + 1),
)::Nothing

    obj_new = filter_mmed(obj, ch = ch, k = k, t = t, ww = ww)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
