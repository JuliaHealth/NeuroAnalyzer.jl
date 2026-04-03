export remove_dc
export remove_dc!

"""
    remove_dc(s, n)

Remove mean value (DC offset).

# Arguments

- `s::AbstractVector`: signal vector
- `n::Union{Int64, Tuple{Int64, Int64}}=0`: if `n` is greater than 0, mean value is calculated for the first `n` samples or if `n` is a tuple greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples

# Returns

- `Vector{Float64}`
"""
function remove_dc(
    s::AbstractVector,
    n::Union{Int64, Tuple{Int64, Int64}} = 0,
)::Vector{Float64}
    # validate
    if isa(n, Int64)
        n >= 0 || throw(ArgumentError("n must be ≥ 0."))
        n <= length(s) || throw(ArgumentError("n must be ≤ $(length(s))."))
        return n == 0 ? s .- mean(s) : s .- mean(s[1:n])

    else
        n != (0, 0) && _check_tuple(n, (1, length(s)), "n")
        return n == (0, 0) ? s .- mean(s) : s .- mean(s[n[1]:n[2]])
    end
end

"""
    remove_dc(s, n)

Remove mean value (DC offset).

# Arguments

- `s::AbstractMatrix`: signal matrix, shape (channels, samples)
- `n::Union{Int64, Tuple{Int64, Int64}}=0`: if `n` is greater than 0, mean value is calculated for the first `n` samples or if `n` is a tuple greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples

# Returns

- `Matrix{Float64}`
"""
function remove_dc(
    s::AbstractMatrix,
    n::Union{Int64, Tuple{Int64, Int64}} = 0,
)::Matrix{Float64}
    # validate
    ch_n = size(s, 1)

    s_new = similar(s, Float64)
    Threads.@threads :static for ch_idx in 1:ch_n
        s_new[ch_idx, :] = remove_dc(@view(s[ch_idx, :]), n)
    end

    return s_new
end

"""
    remove_dc(s, n)

Remove mean value (DC offset) for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `n::Union{Int64, Tuple{Int64, Int64}}=0`: if `n` is greater than 0, mean value is calculated for the first `n` samples or if `n` is a tuple greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples

# Returns

- `Array{Float64, 3}`
"""
function remove_dc(
    s::AbstractArray,
    n::Union{Int64, Tuple{Int64, Int64}} = 0,
)::Array{Float64, 3}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    # calculate over channels and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = remove_dc(@view(s[ch_idx, :, ep_idx]), n)
    end

    return s_new
end

"""
    remove_dc(obj; <keyword arguments>)

Remove mean value (DC offset).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `n::Union{Int64, Tuple{Int64, Int64}}=0`: if `n` is greater than 0, mean value is calculated for the first `n` samples or if `n` is a tuple greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function remove_dc(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    n::Union{Int64, Tuple{Int64, Int64}} = 0,
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] = remove_dc(@view(obj.data[ch, :, :]), n)
    push!(obj_new.history, "remove_dc(obj; ch=$ch, n=$n)")

    return obj_new
end

"""
    remove_dc!(obj; <keyword arguments>)

Remove mean value (DC offset).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `n::Union{Int64, Tuple{Int64, Int64}}=0`: if `n` is greater than 0, mean value is calculated for the first `n` samples or if `n` is a tuple greater than (0, 0), mean value is calculated for `n[1]` to `n[2]` samples

# Returns

- `Nothing`
"""
function remove_dc!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    n::Union{Int64, Tuple{Int64, Int64}} = 0,
)::Nothing
    obj_new = remove_dc(obj; ch = ch, n = n)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
