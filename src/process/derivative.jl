export derivative
export derivative!

"""
    derivative(s)

Return the derivative of a discrete signal using the symmetric difference quotient.

For interior samples: `ds[i] = (s[i+1] − s[i−1]) / 2`.

For the boundary samples the one-sided half-difference is used so that the output has the same length as the input.

# Arguments

- `s::AbstractVector`: signal vector; must contain at least 3 elements

# Returns

- `Vector{Float64}`: derivative signal of the same length as `s`

# Throws

- `ArgumentError`: if `length(s) ≤ 2`

# See also

[`derivative(::AbstractArray)`](@ref), [`derivative(::NeuroAnalyzer.NEURO)`](@ref)
"""
function derivative(s::AbstractVector)::AbstractVector

    length(s) > 2 || throw(ArgumentError("Signal length must be > 2."))

    # half-differences: length n-1
    dv = diff(s) / 2
    s_new = Vector{Float64}(undef, length(s))
    # forward half-difference at left boundary
    s_new[1] = dv[1]
    # symmetric: left + right half-differences
    s_new[2:end-1] = dv[1:end-1] .+ dv[2:end]
    # backward half-difference at right boundary
    s_new[end] = dv[end]

    return s_new

end

"""
    derivative(s)

Return the derivative of a discrete signal using the symmetric difference quotient.

For interior samples: `ds[i] = (s[i+1] − s[i−1]) / 2`.

For the boundary samples the one-sided half-difference is used so that the output has the same length as the input.

# Arguments

- `s::AbstractArray`: signal array, shape `(channels, samples, epochs)`

# Returns

- `Array{Float64, 3}`: derivative array of the same shape as `s`

# Throws

- `ArgumentError`: if `s` is not 3-dimensional

# See also

[`derivative(::AbstractVector)`](@ref), [`derivative(::NeuroAnalyzer.NEURO)`](@ref)
"""
function derivative(s::AbstractArray)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = derivative(@view(s[ch_idx, :, ep_idx]))
    end

    return s_new

end

"""
    derivative(obj; <keyword arguments>)

Return the derivative of a discrete signal using the symmetric difference quotient of selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `NeuroAnalyzer.NEURO`: new object with differentiated channels

# See also

[`derivative!`](@ref), [`derivative(::AbstractArray)`](@ref)
"""
function derivative(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex}
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = derivative(obj.data[ch, :, :])
    push!(obj_new.history, "derivative(OBJ, ch=$ch)")

    return obj_new

end

"""
    derivative!(obj; <keyword arguments>)

Return the derivative of a discrete signal using the symmetric difference quotient in-place of selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Nothing`

# See also

[`derivative`](@ref)
"""
function derivative!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Nothing

    obj_new = derivative(obj, ch = ch)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
