export filter_sg
export filter_sg!

"""
    filter_sg(s; <keyword arguments>)

Filter a signal using a Savitzky-Golay smoothing filter.

Fits a polynomial of degree `order` to successive overlapping windows of `window` samples using least squares, then uses the fitted polynomial to replace each sample. Preserves peak heights and widths better than a simple moving average.

# Arguments

- `s::AbstractVector`: signal vector; must have at least `window` elements
- `order::Int64=6`: polynomial degree; must satisfy `2 ≤ order < window`
- `window::Int64=11`: filter window length (number of coefficients); must be odd and satisfy `1 ≤ window ≤ length(s)`

# Returns

- `Vector{Float64}`: filtered signal of the same length as `s`

# Throws

- `ArgumentError`: if `window` is out of range, even, `order < 2`, or `order ≥ window`

# See also

[`filter_sg(::AbstractArray)`](@ref), [`filter_sg(::NeuroAnalyzer.NEURO)`](@ref)
"""
function filter_sg(
    s::AbstractVector;
    order::Int64 = 6,
    window::Int64 = 11
)::Vector{Float64}

    (window >= 1 && window <= length(s)) || throw(ArgumentError("window must be in [1, $(length(s))]."))
    isodd(window) || throw(ArgumentError("window must be odd."))
    order >= 2 || throw(ArgumentError("order must be ≥ 2."))
    order < window || throw(ArgumentError("order must be < window ($window)."))

    return savitzky_golay(s, window, order).y

end

"""
    filter_sg(s; <keyword arguments>)

Apply a Savitzky-Golay filter to every channel × epoch slice of a 3-D signal array. Delegates to [`filter_sg(::AbstractVector)`](@ref).

# Arguments

- `s::AbstractArray`: 3-D signal array, shape `(channels, samples, epochs)`
- `order::Int64=6`: polynomial degree; must satisfy `2 ≤ order < window`
- `window::Int64=11`: filter window length (number of coefficients); must be odd and satisfy `1 ≤ window ≤ size(s, 2)`

# Returns

- `Array{Float64, 3}`: filtered array of the same shape as `s`

# See also

[`filter_sg(::AbstractVector)`](@ref), [`filter_sg(::NeuroAnalyzer.NEURO)`](@ref)
"""
function filter_sg(
    s::AbstractArray;
    order::Int64 = 6,
    window::Int64 = 11
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
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_filtered[ch_idx, :, ep_idx] = filter_sg(@view(s[ch_idx, :, ep_idx]); order=order, window=window)
    end

    return s_filtered

end

"""
    filter_sg(obj; <keyword arguments>)

Apply a Savitzky-Golay filter to selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `order::Int64=6`: polynomial degree; must satisfy `2 ≤ order < window`
- `window::Int64=11`: filter window length (number of coefficients); must be odd and satisfy `1 ≤ window ≤ epoch_len(obj)`

# Returns

- `NeuroAnalyzer.NEURO`: new object with filtered channels

# See also

[`filter_sg!`](@ref), [`filter_sg(::AbstractArray)`](@ref)
"""
function filter_sg(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    order::Int64 = 6,
    window::Int64 = 11
)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch = ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = filter_sg(obj.data[ch, :, :]; order = order, window = window)
    push!(obj_new.history, "filter_sg(OBJ, ch=$ch, order=$order, window=$window)")

    return obj_new

end

"""
    filter_sg!(obj; <keyword arguments>)

Apply a Savitzky-Golay filter in-place to selected channels of a NEURO object. Delegates to [`filter_sg`](@ref) and copies the result back.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `order::Int64=6`: order of the polynomial used to fit the samples; must be less than `window`
- `window::Int64=11`: length of the filter window (i.e., the number of coefficients); must be an odd number

# Returns

- `Nothing`

# See also
[`filter_sg`](@ref)
"""
function filter_sg!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    order::Int64 = 6,
    window::Int64 = 11
)::Nothing

    obj_new = filter_sg(obj, ch = ch, order = order, window = window)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
