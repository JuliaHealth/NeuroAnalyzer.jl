export remove_dc
export remove_dc!

"""
    remove_dc(s, n)

Remove mean value (DC offset).

# Arguments

- `s::AbstractVector`
- `n::Union{Int64, Tuple{Int64, Int64}}=0`: if `n` is greater than 0, mean value is calculated for the first `n` samples or if `n` is a tuple greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples

# Returns

- `s_new::Vector{Float64}`
"""
function remove_dc(s::AbstractVector, n::Union{Int64, Tuple{Int64, Int64}}=0)

    if isa(n, Int64)
        @assert n >=0 "n must be ≥ 0."
        @assert n <= length(s) "n must be ≤ $(length(s))."

        s_new = n == 0 ? s .- mean(s) : s .- mean(s[1:n])
    else
        n != (0, 0) && _check_tuple(n, "n", (1, length(s)))

        s_new = n == (0, 0) ? s .- mean(s) : s .- mean(s[n[1]:n[2]])
    end

    return s_new

end

"""
    remove_dc(s, n)

Remove mean value (DC offset).

# Arguments

- `s::AbstractArray`
- `n::Union{Int64, Tuple{Int64, Int64}}=0`: if `n` is greater than 0, mean value is calculated for the first `n` samples or if `n` is a tuple greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples

# Returns

- `s::Array{Float64, 3}`
"""
function remove_dc(s::AbstractArray, n::Union{Int64, Tuple{Int64, Int64}}=0)

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
    remove_dc(obj; <keyword arguments>)

Remove mean value (DC offset).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels
- `n::Union{Int64, Tuple{Int64, Int64}}=0`: if `n` is greater than 0, mean value is calculated for the first `n` samples or if `n` is a tuple greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function remove_dc(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, n::Union{Int64, Tuple{Int64, Int64}}=0)

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views remove_dc(obj.data[ch, :, :], n)
    reset_components!(obj_new)
    push!(obj_new.history, "remove_dc(OBJ, ch=$ch, n=$n)")

    return obj_new

end

"""
    remove_dc!(obj; <keyword arguments>)

Remove mean value (DC offset).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels
- `n::Union{Int64, Tuple{Int64, Int64}}=0`: if `n` is greater than 0, mean value is calculated for the first `n` samples or if `n` is a tuple greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples
"""
function remove_dc!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, n::Union{Int64, Tuple{Int64, Int64}}=0)

    obj_new = remove_dc(obj, ch=ch, n=n)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
