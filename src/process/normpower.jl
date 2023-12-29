export normpower
export normpower!

"""
    normpower(s)

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

# Arguments

- `s::AbstractVector`

# Returns

- `s_new::Vector{Float64}`
"""
function normpower(s::AbstractVector)

    return s .* rms(s)

end

"""
    normpower(s)

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

# Arguments

- `s::AbstractArray`

# Returns

- `s_new::Array{Float64, 3}`
"""
function normpower(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views normpower(s[ch_idx, :, ep_idx])
        end
    end

    return s_new

end

"""
    normpower(obj; channel, t)

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function normpower(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)))

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = normpower(obj.data[ch, :, :])
    reset_components!(obj_new)
    push!(obj_new.history, "normpower(OBJ, ch=$ch)")

    return obj_new

end

"""
    normpower!(obj; ch)

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
"""
function normpower!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)))

    obj_new = normpower(obj, ch=ch)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
