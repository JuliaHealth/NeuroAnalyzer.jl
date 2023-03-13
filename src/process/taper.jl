export taper
export taper!

"""
    taper(signal; taper)

Taper the signal.

# Arguments

- `signal::AbstractVector`
- `t::Union{AbstractVector, Vector{ComplexF64}}`

# Returns

- `taper::Vector{Union{Float64, ComplexF64}}`
"""
function taper(signal::AbstractVector; t::Union{AbstractVector, Vector{ComplexF64}})
    length(t) == length(signal) || throw(ArgumentError("Taper and signal lengths must be equal."))
    return signal .* t
end

"""
    taper(obj; channel, taper)

Taper channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `taper::Union{Vector{Real, Vector{ComplexF64}}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function taper(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), t::Union{Vector{<:Real}, Vector{ComplexF64}})

    _check_channels(obj, channel)

    ep_n = epoch_n(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
            @views obj_new.data[channel[ch_idx], :, ep_idx] = taper(obj_new.data[channel[ch_idx], :, ep_idx], t=t)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "taper(OBJ, channel=$channel), t=$t")

    return obj_new
end

"""
    taper!(obj; channel, taper)

Taper channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `t::Union{Vector{<:Real}, Vector{ComplexF64}}`
"""
function taper!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    obj_tmp = taper(obj, channel=channel, t=t)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
