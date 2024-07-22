export derivative
export derivative!

"""
    derivative(s)

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

# Arguments

- `s::AbstractVector`

# Returns

- `s_new::AbstractVector`
"""
function derivative(s::AbstractVector)

    @assert length(s) > 2 "Signal length must be > 2."

    dv  = diff(s) / 2       # half the derivative
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
function derivative(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views derivative(s[ch_idx, :, ep_idx])
        end
    end

    return s_new

end

"""
    derivative(obj; <keyword arguments>)

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function derivative(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = derivative(obj.data[ch, :, :])
    reset_components!(obj_new)
    push!(obj_new.history, "derivative(OBJ, ch=$ch)")

    return obj_new

end

"""
    derivative!(obj; <keyword arguments>)

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels
"""
function derivative!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})

    obj_new = derivative(obj, ch=ch)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
