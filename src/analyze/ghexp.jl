export ghexp

# convenience alias kept local to reduce repetition in signatures
const _QRange = Union{
    Nothing,
    StepRangeLen{
        Float64,
        Base.TwicePrecision{Float64},
        Base.TwicePrecision{Float64},
        Int64,
    },
}

"""
    _ghexp(s; <keyword arguments>)

Calculate the Generalised Hurst Exponent(s) of a 1-D signal vector by analysing how the q-th moment of the absolute increments |s(t+τ) − s(t)|^q scales with lag τ.

Two modes:
- q_range = nothing → standard Hurst exponent via hurst_exponent(); output shape: (1, 2) → (exponent, goodness-of-fit)
- q_range provided →  generalised exponents via generalised_hurst_range(); output shape: (length(q_range), 2)

# Arguments

- `s::AbstractVector`: signal vector
- `tau_range::UnitRange{Int64}`: lag range over which the q-th moment of absolute increments is estimated
- `q_range::Union{Nothing, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}=nothing`: moments at which GHEs are estimated; if `nothing`, the standard Hurst exponent is returned

# Returns

- `Matrix{Float64}`: shape `(1, 2)` when `q_range=nothing`, otherwise `(length(q_range), 2)` - columns are (exponent, goodness-of-fit)
"""
function _ghexp(
    s::AbstractVector;
    tau_range::UnitRange{Int64},
    q_range::_QRange = nothing,
)::Matrix{Float64}
    # validate
    tau_range[end] < length(s) ||
        throw(
            ArgumentError(
                "End of tau_range ($(tau_range[end])) must be < length of s ($(length(s))).",
            ),
        )

    if isnothing(q_range)
        ghe = hurst_exponent(Vector(s), tau_range)
    else
        ghe = generalised_hurst_range(Vector(s), tau_range, q_range)
    end

    return ghe
end

"""
    _ghexp(s)

Calculate the Generalised Hurst Exponents (GHEs) for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `tau_range::UnitRange{Int64}`: lag range over which the q-th moment of absolute increments is estimated
- `q_range::Union{Nothing, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}=nothing`: moments at which GHEs are estimated; if `nothing`, the standard Hurst exponent is returned

# Returns

- `Array{Float64, 4}`: shape (channels, q, 2, epochs) where `q` is 1 when `q_range=nothing`, otherwise `length(q_range)`
"""
function _ghexp(
    s::AbstractArray;
    tau_range::UnitRange{Int64},
    q_range::_QRange = nothing,
)::Array{Float64, 4}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # determine the q dimension from the range, or 1 for the standard exponent
    q_n = isnothing(q_range) ? 1 : length(q_range)
    ghe = zeros(ch_n, q_n, 2, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ghe[ch_idx, :, :, ep_idx] = _ghexp(
            @view(s[ch_idx, :, ep_idx]),
            tau_range = tau_range,
            q_range = q_range,
        )
    end

    return ghe
end

"""
    ghexp(obj; <keyword arguments>)

Calculate the Generalised Hurst Exponents (GHEs) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `exclude_bad_segments::Bool=NeuroAnalyzer.exclude_bad_segments`: if `true`, bad segments (markers of type `"BAD"` of all channels or of the selected channels) of a continuous signal are not used: their samples are removed and the remaining parts are joined (default: the `exclude_bad_segments` preference, see `set_exclude_bad_segments()`)
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `tau_range::UnitRange{Int64}`: lag range over which the q-th moment of absolute increments is estimated
- `q_range::Union{Nothing, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}=nothing`: moments at which GHEs are estimated; if `nothing`, the standard Hurst exponent is returned

# Returns

- `Array{Float64, 4}`: shape `(channels, q, 2, epochs)`
"""
function ghexp(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    tau_range::UnitRange{Int64},
    q_range::_QRange = nothing,
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
    exclude_bad_segments::Bool = NeuroAnalyzer.exclude_bad_segments,
)::Array{Float64, 4}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch = _get_ch(obj, ch, exclude_bad_channels)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _ghexp(_seg_data(obj, ch, exclude_bad_segments); tau_range = tau_range, q_range = q_range)
end
