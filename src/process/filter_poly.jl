export filter_poly
export filter_poly!

"""
    filter_poly(s; <keyword arguments>)

Filter a signal using a piecewise polynomial (Savitzky-Golay-style) filter.

The signal is split into non-overlapping windows. A polynomial of degree `order` is fitted to each window and used to replace the raw samples. Window boundaries are then smoothed with a Loess regression to suppress discontinuities between adjacent windows.

# Arguments

- `s::AbstractVector`: signal vector
- `order::Int64=8`: polynomial degree; must be ≥ 2 and < `window`
- `window::Int64=10`: window length in samples; must satisfy `1 ≤ window ≤ length(s)`

# Returns

- `Vector{Float64}`: filtered signal of the same length as `s`

# Throws

- `ArgumentError`: if `order < 2`, `order ≥ window`, or `window` is out of range

# See also

[`filter_poly(::AbstractArray)`](@ref), [`filter_poly(::NeuroAnalyzer.NEURO)`](@ref)
"""
function filter_poly(
    s::AbstractVector;
    order::Int64 = 8,
    window::Int64 = 10
)::Vector{Float64}

    order >= 2 || throw(ArgumentError("order must be ≥ 2."))
    (window >= 1 && window <= length(s)) || throw(ArgumentError("window must be in [1, $(length(s))]."))
    order < window || throw(ArgumentError("order must be < window ($window)."))

    s_filtered = float(copy(s))

    # TO DO: smooth spikes between windows
    window_n = length(s) ÷ window
    # remaining samples after full windows
    window_last = length(s) - window_n * window

    # --- fit and replace interior windows (all except the last) ---
    @inbounds for window_idx in 1:(window_n - 1)
        i1 = (window_idx - 1) * window + 1
        i2 = window_idx * window
        s_tmp = s[i1:i2]
        t = eachindex(s_tmp)
        p = Polynomials.fit(t, s_tmp, order)
        s_filtered[i1:i2] = [p(ti) for ti in t]
    end

    # --- fit and replace the last (possibly longer) window ---
    i1 = (window_n - 1) * window + 1
    i2 = window_n * window + window_last
    s_tmp = s[i1:i2]
    t = eachindex(s_tmp)
    p = Polynomials.fit(t, s_tmp, order)
    s_filtered[i1:i2] = [p(ti) for ti in t]

    # --- smooth window junctions with Loess to suppress discontinuities ---
    # half-width of the junction smoothing region
    half = window ÷ 4
    @inbounds for junction in window:window:(window * (window_n - 1))
        j1 = max(1, junction - half)
        j2 = min(length(s), junction + half - 1)   # was: no bounds check — could exceed signal length
        s_tmp = s_filtered[j1:j2]
        t = collect(1.0:length(s_tmp))
        model = Loess.loess(t, Vector{Float64}(s_tmp); span=1.0)
        smoothed = Loess.predict(model, t)
        # preserve edge values to avoid boundary artifacts
        smoothed[1] = s_tmp[1]
        smoothed[end] = s_tmp[end]
        s_filtered[j1:j2] = smoothed
    end

    return s_filtered

end

"""
    filter_poly(s; <keyword arguments>)

Apply a piecewise polynomial filter to every channel × epoch slice of a 3-D signal array. Delegates to [`filter_poly(::AbstractVector)`](@ref).

# Arguments

- `s::AbstractArray`: 3-D signal array, shape `(channels, samples, epochs)`
- `order::Int64=8`: polynomial degree; must be ≥ 2 and < `window`
- `window::Int64=10`: window length in samples; must satisfy `1 ≤ window ≤ size(s, 2)`

# Returns

# Returns

- `Array{Float64, 3}`: filtered array of the same shape as `s`

# See also

[`filter_poly(::AbstractVector)`](@ref), [`filter_poly(::NeuroAnalyzer.NEURO)`](@ref)
"""
function filter_poly(
    s::AbstractArray;
    order::Int64 = 8,
    window::Int64 = 10
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
        s_filtered[ch_idx, :, ep_idx] = filter_poly(
                                            @view(s[ch_idx, :, ep_idx]),
                                            order=order,
                                            window=window
                                        )
    end

    return s_filtered

end

"""
    filter_poly(obj; <keyword arguments>)

Apply a piecewise polynomial filter to selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `order::Int64=8`: polynomial degree; must be ≥ 2 and < `window`
- `window::Int64=10`: window length in samples; must satisfy `1 ≤ window ≤ epoch_len(obj)`

# Returns

- `NeuroAnalyzer.NEURO`: new object with filtered channels

# See also

[`filter_poly!`](@ref), [`filter_poly(::AbstractArray)`](@ref)
"""
function filter_poly(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    order::Int64 = 8,
    window::Int64 = 10
)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch = ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = filter_poly(
                                @view(obj.data[ch, :, :]),
                                order=order,
                                window=window
                             )
    push!(obj_new.history, "filter_poly(OBJ, ch=$ch, order=$order, window=$window)")

    return obj_new

end

"""
    filter_poly!(obj; <keyword arguments>)

Apply a piecewise polynomial filter in-place to selected channels of a NEURO object. Delegates to [`filter_poly`](@ref) and copies the result back.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `order::Int64=8`: polynomial degree; must be ≥ 2 and < `window`
- `window::Int64=10`: window length in samples; must satisfy `1 ≤ window ≤ epoch_len(obj)`

# Returns

- `Nothing`

# See also
[`filter_poly`](@ref)
"""
function filter_poly!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    order::Int64 = 8,
    window::Int64 = 10
)::Nothing

    obj_new = filter_poly(obj, ch = ch, order = order, window = window)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
