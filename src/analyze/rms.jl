export rms
export rmse

"""
    rms(s)

alculate Root Mean Square (RMS) of a 1-D signal vector.

RMS = √( mean(s²) ) = norm(s) / √length(s)

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Float64`: RMS value
"""
function rms(s::AbstractVector)::Float64

    # equivalent to sqrt(mean(s.^2)) but avoids an intermediate allocation
    return norm(s) / sqrt(length(s))
end

"""
    rms(s)

Calculate Root Mean Square (RMS) for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

- `Matrix{Float64}`: RMS values, shape (channels, epochs)
"""
function rms(s::AbstractArray)::Matrix{Float64}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    r = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        r[ch_idx, ep_idx] = rms(@view(s[ch_idx, :, ep_idx]))
    end

    return r
end

"""
    rms(obj; <keyword arguments>)

Calculate Root Mean Square (RMS) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}: list of channels
- `ep::Union{Int64, Vector{Int64}, UnitRange{Int64}}=_c(nepochs(obj))`: default use all epochs

# Returns

- `Matrix{Float64}`: RMS
"""
function rms(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ep::Union{Int64, Vector{Int64}, UnitRange{Int64}} = _c(nepochs(obj)),
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    _check_epochs(obj, ep)
    isa(ep, Int64) && (ep = [ep])

    return rms(@view(obj.data[ch, :, ep]))
end

"""
    rmse(s1, s2)

Calculate Root Mean Square Error (RMSE).

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

- `Float64`: RMSE
"""
function rmse(s1::AbstractVector, s2::AbstractVector)::Float64

    # validate
    length(s1) == length(s2) || throw(ArgumentError("s1 and s2 must have the same length."))

    return sqrt(mean((s2 .- s1) .^ 2))
end

"""
    rmse(s1, s2)

Calculate Root Mean Square Error (RMSE) for two 3-D signal arrays.

# Arguments

- `s1::AbstractArray`: signal array, shape (channels, samples, epochs)
- `s2::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

- `Matrix{Float64}`: RMSE
"""
function rmse(s1::AbstractArray, s2::AbstractArray)::Matrix{Float64}

    # validate
    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))
    _chk3d(s1)
    _chk3d(s2)

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    r = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        r[ch_idx, ep_idx] = rmse(
            @view(s1[ch_idx, :, ep_idx]),
            @view(s2[ch_idx, :, ep_idx])
        )
    end

    return r
end

"""
    rmse(obj1, obj2; <keyword arguments>)

Calculate Root Mean Square Error (RMSE) for two NEURO objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s) in `obj1`
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s) in `obj2`
- `ep1::Union{Int64, Vector{Int64}, UnitRange{Int64}}=_c(nepochs(obj1))`: epoch number(s) in `obj1`
- `ep2::Union{Int64, Vector{Int64}, UnitRange{Int64}}=_c(nepochs(obj2))`: epoch number(s) in `obj2`

# Returns

- `Matrix{Float64}`: RMSE
"""
function rmse(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, UnitRange{Int64}} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, UnitRange{Int64}} = _c(nepochs(obj2)),
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 =
        exclude_bads ? get_channel(obj1; ch = ch1, exclude = "bad") :
        get_channel(obj1; ch = ch1, exclude = "")
    ch2 =
        exclude_bads ? get_channel(obj2; ch = ch2, exclude = "bad") :
        get_channel(obj2; ch = ch2, exclude = "")
    isempty(ch1) && throw(ArgumentError("No channels selected."))
    isempty(ch2) && throw(ArgumentError("No channels selected."))

    # validate
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    length(ch1) == length(ch2) ||
        throw(
            ArgumentError(
                "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal.",
            ),
        )
    length(ep1) == length(ep2) ||
        throw(
            ArgumentError(
                "Lengths of ep1 ($(length(ep1)) and ep2 ($(length(ep2)) must be equal.",
            ),
        )
    epoch_len(obj1) == epoch_len(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    return rmse(@view(obj1.data[ch1, :, ep1]), @view(obj2.data[ch2, :, ep2]))
end
