export add_signal
export add_signal!

"""
    add_signal(s1, s2)

Adds noise.

# Arguments

- `s1::AbstractVector`: target signal
- `s2::AbstractVector`: signal to be added

# Returns

- `s_noisy::AbstractVector`
"""
function add_signal(s1::AbstractVector, s2::AbstractVector)

    length(s1) == length(s2) || throw(ArgumentError("s1 and s2 must have the same length."))

    return s1 .+ s2

end

"""
    add_signal(obj; ch, s)

Add signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `s::AbstractVector`: signal to be added to each channel

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_signal(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), s::AbstractVector)

    _check_channels(obj, ch)

    ch_n = length(ch)
    ep_n = epoch_n(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            @views obj_new.data[ch[ch_idx], :, ep_idx] = add_signal(obj_new.data[ch[ch_idx], :, ep_idx], s)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.history, "add_signal(OBJ, ch=$ch)")

    return obj_new

end

"""
    add_signal!(obj; ch, n)

Add signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `s::AbstractVector`: signal to be added to each channel
"""
function add_signal!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), s::AbstractVector)

    obj_new = add_signal(obj, ch=ch, s=s)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end