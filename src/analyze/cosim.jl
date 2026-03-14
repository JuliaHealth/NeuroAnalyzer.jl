export cosim

"""
    cosim(s1, s2)

Measures the cosine of the angle between two signal vectors:

CS = (s1 · s2) / (‖s1‖ · ‖s2‖) ∈ [-1, 1]

- CS = 1 → identical direction (perfectly similar)
- CS = 0 → orthogonal (no similarity)
- CS = -1 → opposite direction

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

- `cs::Float64`: cosine similarity value ∈ [-1, 1]
"""
function cosim(s1::AbstractVector, s2::AbstractVector)::Float64

    !(length(s1) == length(s2)) && throw(ArgumentError("Both signals must have the same length."))

    # CS = (s1 · s2) / (‖s1‖ · ‖s2‖)
    cs = dot(s1, s2) / (norm(s1) * norm(s2))

    return cs

end

"""
    cosim(obj1, obj2; <keyword arguments>)

Measures the cosine of the angle between paired channels across two objects.

CS = (s1 · s2) / (‖s1‖ · ‖s2‖) ∈ [-1, 1]

- CS = 1 → identical direction (perfectly similar)
- CS = 0 → orthogonal (no similarity)
- CS = -1 → opposite direction

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

- `cs::Matrix{Float64}`: cosine similarity values, shape `(channels, epochs)`
"""
function cosim(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    !(length(ch1) == length(ch2)) && throw(ArgumentError("Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."))

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    !(length(ep1) == length(ep2)) && throw(ArgumentError("Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."))
    !(epoch_len(obj1) == epoch_len(obj2)) && throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    # number of channels
    ch_n = length(ch1)
    # number of epochs
    ep_n = length(ep1)

    # pre-allocate output
    cs = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        cs[ch_idx, ep_idx] = cosim(
            @view(obj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(obj2.data[ch2[ch_idx], :, ep2[ep_idx]]),
        )
    end

    return cs

end

"""
    cosim(obj; <keyword arguments>)

Measures the cosine of the angle between all channel pairs within one object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `cs::Array{Float64, 3}`: cosine similarity values, shape `(channels, channels, epochs)`
"""
function cosim(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = length(ep)

    # pre-allocate output
    cs = zeros(ch_n, ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 in 1:ch_idx1
            cs[ch_idx1, ch_idx2, ep_idx] = cosim(
                @view(obj.data[ch[ch_idx1], :, ep_idx]),
                @view(obj.data[ch[ch_idx2], :, ep_idx]),
            )
        end
    end

    # mirror the lower triangle to the upper triangle
    # to produce the full symmetric similarity matrix
    cs = _copy_lt2ut(cs)

    return cs

end
