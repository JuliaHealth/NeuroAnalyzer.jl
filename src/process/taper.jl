export taper
export taper!

"""
    taper(signal; taper)

Taper the signal.

# Arguments

- `signal::AbstractVector`
- `t::Union{AbstractVector, Vector{ComplexF64}}`

# Returns

- `t::Vector{Union{Float64, ComplexF64}}`
"""
function taper(signal::AbstractVector; t::Union{AbstractVector, Vector{ComplexF64}})

    length(t) == length(signal) || throw(ArgumentError("Taper and signal lengths must be equal."))

    return signal .* t

end

"""
    taper(s; t)

Taper the signal.

# Arguments

- `s::AbstractArray`
- `t::Union{Vector{Real, Vector{ComplexF64}}`

# Returns

- `s::Array{Float64, 3`
"""
function taper(s::AbstractArray; t::Union{Vector{<:Real}, Vector{ComplexF64}})

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views taper(s[ch_idx, :, ep_idx], t=t)
        end
    end

    return s_new

end

"""
    taper(obj; channel, t)

Taper the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `t::Union{Vector{Real, Vector{ComplexF64}}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function taper(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), t::Union{Vector{<:Real}, Vector{ComplexF64}})

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = taper(obj.data[ch, :, :], t=t)
    reset_components!(obj_new)
    push!(obj_new.history, "taper(OBJ, ch=$ch), t=$t")

    return obj_new

end

"""
    taper!(obj; ch, t)

Taper the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `t::Union{Vector{<:Real}, Vector{ComplexF64}}`
"""
function taper!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    obj_new = taper(obj, ch=ch, t=t)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
