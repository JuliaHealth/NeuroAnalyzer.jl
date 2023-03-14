export add_noise
export add_noise!

"""
    add_noise(s, n)

Adds noise.

# Arguments

- `s::AbstractVector`
- `n::AbstractVector`

# Returns

- `s_noisy::AbstractVector`
"""
function add_noise(s::AbstractVector, n::AbstractVector)

    length(s) == length(n) || throw(ArgumentError("Length of signal and noise must be equal."))

    return s .+ n

end

"""
    add_noise(obj; ch, n)

Add noise.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `n::AbstractVector`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_noise(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), n::AbstractVector)

    _check_channels(obj, ch)

    ch_n = length(ch)
    ep_n = epoch_n(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            @views obj_new.data[ch[ch_idx], :, ep_idx] = add_noise(obj_new.data[ch[ch_idx], :, ep_idx], n=n)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "add_noise(OBJ, ch=$ch)")

    return obj_new
end

"""
    add_noise!(obj; ch, n)

Add noise.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `n::AbstractVector`
"""
function add_noise!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), noise::AbstractVector)

    obj_tmp = add_noise(obj, ch=ch, n=n)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing
end
