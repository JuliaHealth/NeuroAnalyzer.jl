export mutual_information

"""
    mutual_information(s1, s2)

Calculate mutual information between two signals.

Wraps `InformationMeasures.get_mutual_information()` for mutual information estimation using the maximum likelihood estimator.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

- `mi::Float64`: mutual information
"""
function mutual_information(s1::AbstractVector, s2::AbstractVector)::Float64

    mi = InformationMeasures.get_mutual_information(s1, s2)

    return mi

end

"""
    mutual_information(s1, s2)

Calculate mutual information between matched channel pairs of two 3-D arrays.

Wraps `InformationMeasures.get_mutual_information()` for mutual information estimation using the maximum likelihood estimator.

# Arguments

- `s1::AbstractArray`: signal array (channels × samples × epochs)
- `s2::AbstractArray`: signal array (channels × samples × epochs)

# Returns

- `mi::Matrix{Float64}`: MI matrix, shape `(channels, epochs)`
"""
function mutual_information(s1::AbstractArray, s2::AbstractArray)::Matrix{Float64}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s1)
    _chk3d(s2)

    # number of channels
    ch_n = size(s1, 1)
    # number of epochs
    ep_n = size(s1, 3)

    # pre-allocate output
    mi = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        mi[ch_idx, ep_idx] = mutual_information(
            @view(s1[ch_idx, :, ep_idx]),
            @view(s2[ch_idx, :, ep_idx]),
        )
    end

    return mi

end

"""
    mutual_information(s)

Calculate mutual information between all channel pairs.

Wraps `InformationMeasures.get_mutual_information()` for mutual information estimation using the maximum likelihood estimator.

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)

# Returns

- `mi::Array{Float64, 3}`: symmetric MI matrix, shape `(channels, channels, epochs)`
"""
function mutual_information(s::AbstractArray)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # initialize progress bar
    progbar = Progress(ep_n * ch_n, dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    mi = zeros(ch_n, ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads for ep_idx in 1:ep_n
        for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1 - 1
                mi[ch_idx1, ch_idx2, ep_idx] = mutual_information(
                    @view(s[ch_idx1, :, ep_idx]),
                    @view(s[ch_idx2, :, ep_idx]),
                )
            end

            progress_bar && next!(progbar)
        end
    end

    # mirror lower triangle to upper triangle
    mi = _copy_lt2ut(mi)

    return mi

end

"""
    mutual_information(obj; <keyword arguments>)

Calculate mutual information between all channel pairs.

Wraps `InformationMeasures.get_mutual_information()` for mutual information estimation using the maximum likelihood estimator.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `mi::Array{Float64, 3}`: symmetric MI matrix, shape `(channels, channels, epochs)`
"""
function mutual_information(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    mi = mutual_information(@view(obj.data[ch, :, :]))

    return mi

end

"""
    mutual_information(obj1, obj2; <keyword arguments>)

Calculate mutual information between matched channel pairs.

Wraps `InformationMeasures.get_mutual_information()` for mutual information estimation using the maximum likelihood estimator.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

- `m::Matrix{Float64}`: MI matrix, shape `(channels, epochs)`
"""
function mutual_information(
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
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    mi = mutual_information(
        @view(obj1.data[ch1, :, ep1]),
        @view(obj2.data[ch2, :, ep2]),
    )

    return mi

end
