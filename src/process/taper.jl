export taper
export taper!

"""
    taper(s; <keyword arguments>)

Taper the signal.

# Arguments

- `s::AbstractVector`: signal vector
- `t::Vector{<:Real}`

# Returns

- `Vector{Float64}`
"""
function taper(s::AbstractVector; t::Vector{<:Real})::Vector{Float64}

    # validate
    length(t) == length(s) ||
        throw(ArgumentError("Taper and signal lengths must be equal."))

    return s .* t
end

"""
    taper(s; <keyword arguments>)

Taper a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `t::Vector{<:Real}`

# Returns

- `Array{Float64, 3}`
"""
function taper(s::AbstractArray; t::Vector{<:Real})::Array{Float64, 3}
    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s, Float64)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = taper(@view(s[ch_idx, :, ep_idx]), t = t)
    end

    return s_new
end

"""
    taper(obj; <keyword arguments>)

Taper the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `t::Vector{<:Real}`

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function taper(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        t::Vector{<:Real},
    )::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] = taper(obj.data[ch, :, :]; t = t)
    push!(obj_new.history, "taper(obj; ch=$ch), t=$t")

    return obj_new
end

"""
    taper!(obj; <keyword arguments>)

Taper the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `t::Vector{<:Real}`

# Returns

- `Nothing`
"""
function taper!(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        t::Vector{<:Real},
    )::Nothing
    obj_new = taper(obj; ch = ch, t = t)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
