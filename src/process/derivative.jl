export derivative
export derivative!

"""
    derivative(s)

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `s_new::AbstractVector`
"""
function derivative(s::AbstractVector)::AbstractVector

    !(length(s) > 2) && throw(ArgumentError("Signal length must be > 2."))

    dv = diff(s) / 2       # half the derivative
    s_new = [dv[1]; dv]     # copies first element
    s_new .+= [dv; dv[end]] # copies last element, add both results to compute average

    return s_new

end

"""
    derivative(s)

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

# Arguments

- `s::AbstractArray`

# Returns

- `s_new::Array{Float64, 3}`
"""
function derivative(s::AbstractArray)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :dynamic for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views derivative(s[ch_idx, :, ep_idx])
        end
    end

    return s_new

end

"""
    derivative(obj; <keyword arguments>)

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: output NEURO object
"""
function derivative(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch = ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = derivative(obj.data[ch, :, :])
    push!(obj_new.history, "derivative(OBJ, ch=$ch)")

    return obj_new

end

"""
    derivative!(obj; <keyword arguments>)

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Nothing`
"""
function derivative!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Nothing

    obj_new = derivative(obj, ch = ch)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
