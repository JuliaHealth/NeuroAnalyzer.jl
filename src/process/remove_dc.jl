export remove_dc
export remove_dc!

"""
    remove_dc(s, n)

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

# Arguments

- `s::AbstractVector`
- `n::Int64=0`: baseline is the first `n` samples

# Returns

- `s_new::Vector{Float64}`
"""
function remove_dc(s::AbstractVector, n::Int64=0)
    
    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > length(s) && throw(ArgumentError("n must be ≤ $(length(s))."))

    if n == 0
        s_new = s .- mean(s)
    else
        s_new = s .- mean(s[1:n])
    end

    return s_new

end

"""
    remove_dc(s; n)

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

# Arguments

- `s::AbstractArray`
- `n::Int64=0`: baseline is the first `n` samples

# Returns

- `s::Array{Float64, 3}`
"""
function remove_dc(s::AbstractArray, n::Int64=0)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = remove_dc(s[ch_idx, :, ep_idx], n)
        end
    end

    return s_new

end

"""
    remove_dc(obj; ch, n)

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `n::Int64=0`: baseline is the first `n` samples

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function remove_dc(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), n::Int64=0)

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views remove_dc(obj.data[ch, :, :], n)
    reset_components!(obj_new)
    push!(obj_new.history, "remove_dc(OBJ, ch=$ch, n=$n)")

    return obj_new

end

"""
    remove_dc!(obj; ch, n)

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `n::Int64=0`: baseline is the first `n` samples
"""
function remove_dc!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), n::Int64=0)

    obj_new = remove_dc(obj, ch=ch, n=n)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
