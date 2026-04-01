export sumsim

"""
    sumsim(s1, s2; theta)

Calculate summed similarity using an exponential decay model between two 1-D signal vectors.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `theta::Real`: decay parameter

# Returns

- `Float64`: summed similarity

# Notes

Values of `ss` are in the range [0, 1]; higher value indicates larger similarity.
"""
function sumsim(s1::AbstractVector, s2::AbstractVector; theta::Real)::Float64

    # validate
    length(s1) == length(s2) ||
        throw(
            ArgumentError(
                "Lengths of s1 ($(length(s1))) and s2 ($(length(s2))) must be equal.",
            ),
        )

    ss = exp(-theta * sqrt(sum((s1 .- s2) .^ 2)))

    return ss
end

"""
    sumsim(s1, s2; theta)

Calculate summed similarity using an exponential decay model between two 3-D signal arrays.

# Arguments

- `s1::AbstractArray`: signal array, shape (channels, samples, epochs)
- `s2::AbstractArray`: signal array, shape (channels, samples, epochs)
- `theta::Real`: decay parameter

# Returns

- `Matrix{Float64}`: summed similarity

# Notes

Values of `ss` are in the range [0, 1]; higher value indicates larger similarity.
"""
function sumsim(s1::AbstractArray, s2::AbstractArray; theta::Real)::Matrix{Float64}

    # validate
    size(s1) == size(s2) || throw(
        ArgumentError("Sizes of s1 ($(size(s1))) and s2 ($(size(s2))) must be equal."),
    )
    _chk3d(s1)
    _chk3d(s2)

    ch_n = size(s1, 1)
    ep_n = size(s2, 3)

    ss = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ss[ch_idx, ep_idx] = sumsim(
            @view(s1[ch_idx, :, ep_idx]),
            @view(s2[ch_idx, :, ep_idx]),
            theta = theta,
        )
    end

    return ss
end

"""
    sumsim(obj; <keyword arguments>)

Calculate summed similarity using an exponential decay model between two NEURO objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s) in `obj1`
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s) in `obj2`
- `ep1::Union{Int64, Vector{Int64}, UnitRange{Int64}}=_c(nepochs(obj1))`: epoch number(s) in `obj1`
- `ep2::Union{Int64, Vector{Int64}, UnitRange{Int64}}=_c(nepochs(obj2))`: epoch number(s) in `obj2`
- `theta::Real`: decay parameter

# Returns

- `Matrix{Float64}`: summed similarity

# Notes

Values of `ss` are in the range [0, 1]; higher value indicates larger similarity.
"""
function sumsim(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, UnitRange{Int64}} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, UnitRange{Int64}} = _c(nepochs(obj2)),
    theta::Real,
)::Matrix{Float64}

    # validate
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

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 =
        exclude_bads ? get_channel(obj1; ch = ch1, exclude = "bad") :
        get_channel(obj1; ch = ch1, exclude = "")
    ch2 =
        exclude_bads ? get_channel(obj2; ch = ch2, exclude = "bad") :
        get_channel(obj2; ch = ch2, exclude = "")
    isempty(ch1) && throw(ArgumentError("No channels selected."))
    isempty(ch2) && throw(ArgumentError("No channels selected."))
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    ep1 = _n2v(ep1)
    ep2 = _n2v(ep2)

    return sumsim(
        @view(obj1.data[ch1, :, ep1]),
        @view(obj2.data[ch2, :, ep2]);
        theta = theta,
    )
end
