export remove_dc
export remove_dc!

"""
    remove_dc(s, n)

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

# Arguments

- `s::AbstractVector`
- `n::Int64`: baseline is the first `n` samples

# Returns

- `s_new::Vector{Float64}`
"""
function remove_dc(s::AbstractVector, n::Int64)

    @assert n >=0 "n must be ≥ 0."
    @assert n <= length(s) "n must be ≤ $(length(s))."

    s_new = n == 0 ? s .- mean(s) : s .- mean(s[1:n])

    return s_new

end

"""
    remove_dc(s, n)

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

# Arguments

- `s::AbstractArray`
- `n::Int64`: baseline is the first `n` samples

# Returns

- `s::Array{Float64, 3}`
"""
function remove_dc(s::AbstractArray, n::Int64)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views remove_dc(s[ch_idx, :, ep_idx], n)
        end
    end

    return s_new

end

"""
    remove_dc(obj; ch, n)

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `n::Int64=0`: baseline is the first `n` samples

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function remove_dc(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), n::Int64=0)

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
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `n::Int64=0`: baseline is the first `n` samples
"""
function remove_dc!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), n::Int64=0)

    obj_new = remove_dc(obj, ch=ch, n=n)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    remove_dc(s, n)

Remove mean value (DC offset). If `n` is greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples.

# Arguments

- `s::AbstractVector`
- `n::Tuple{Int64, Int64}`: baseline is the `n[1]` to `n[2]` samples

# Returns

- `s_new::Vector{Float64}`
"""
function remove_dc(s::AbstractVector, n::Tuple{Int64, Int64})

    n != (0, 0) && _check_tuple(n, "n", (1, length(s)))

    s_new = n == (0, 0) ? s .- mean(s) : s .- mean(s[n[1]:n[2]])

    return s_new

end

"""
    remove_dc(s; n)

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

# Arguments

- `s::AbstractArray`
- `n::Tuple{Int64, Int64}`: baseline is the `n[1]` to `n[2]` samples

# Returns

- `s::Array{Float64, 3}`
"""
function remove_dc(s::AbstractArray, n::Tuple{Int64, Int64})

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    n != (0, 0) && _info("Baseline at $(n[1]):$(n[2]) samples")

    s_new = similar(s)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views remove_dc(s[ch_idx, :, ep_idx], n)
        end
    end

    return s_new

end

"""
    remove_dc(obj; ch, n)

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `n::Tuple{Int64, Int64}=(0, 0)`: baseline is the `n[1]` to `n[2]` samples

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function remove_dc(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), n::Tuple{Int64, Int64}=(0, 0))

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
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `n::Tuple{Int64, Int64}=(0, 0)`: baseline is the `n[1]` to `n[2]` samples
"""
function remove_dc!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), n::Tuple{Int64, Int64}=(0, 0))

    obj_new = remove_dc(obj, ch=ch, n=n)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
