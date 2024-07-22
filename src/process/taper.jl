export taper
export taper!

"""
    taper(s; <keyword arguments>)

Taper the signal.

# Arguments

- `s::AbstractVector`
- `t::Union{AbstractVector, Vector{ComplexF64}}`

# Returns

- `s_new::Vector{Union{Float64, ComplexF64}}`
"""
function taper(s::AbstractVector; t::Union{AbstractVector, Vector{ComplexF64}})

    @assert length(t) == length(s) "Taper and signal lengths must be equal."

    return s .* t

end

"""
    taper(s; <keyword arguments>)

Taper the signal.

# Arguments

- `s::AbstractArray`
- `t::Union{Vector{Real, Vector{ComplexF64}}`

# Returns

- `s_new::Array{Float64, 3`
"""
function taper(s::AbstractArray; t::Union{Vector{<:Real}, Vector{ComplexF64}})

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views taper(s[ch_idx, :, ep_idx], t=t)
        end
    end

    return s_new

end

"""
    taper(obj; <keyword arguments>)

Taper the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels
- `t::Union{Vector{Real, Vector{ComplexF64}}`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function taper(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, t::Union{Vector{<:Real}, Vector{ComplexF64}})

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = taper(obj.data[ch, :, :], t=t)
    reset_components!(obj_new)
    push!(obj_new.history, "taper(OBJ, ch=$ch), t=$t")

    return obj_new

end

"""
    taper!(obj; <keyword arguments>)

Taper the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels
- `t::Union{Vector{<:Real}, Vector{ComplexF64}}`
"""
function taper!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    obj_new = taper(obj, ch=ch, t=t)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
