export scovm

"""
    scovm(s1, s2; <keyword arguments>)

Calculate corrected or uncorrected sample covariance matrix of two signals (S = κ (Xc' * Xc), where Xc = vcat(s1', s2')).

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `corrected::Bool=false`: if `true`, calculate corrected covariance (κ = number of samples, otherwise κ = number of samples - 1)

# Returns

- `Matrix{Float64}`: 2×2 covariance matrix
"""
function scovm(s1::AbstractVector, s2::AbstractVector; corrected::Bool = false)::Matrix{Float64}
    # validate
    length(s1) == length(s2) || throw(ArgumentError("s1 and s2 must have the same length."))

    # compute the 2×2 channels-vs-channels covariance matrix
    cm = cov(SimpleCovariance(corrected = corrected), hcat(s1, s2))

    return cm
end

"""
    scovm(s; <keyword arguments>)

Calculate sample covariance matrix of a matrix (S = κ (Xc' * Xc)).

# Arguments

- `s::AbstractMatrix`: signal matrix (channels, samples)
- `corrected::Bool=false`: if `true`, calculate corrected covariance (κ = number of samples, otherwise κ = number of samples - 1)

# Returns

- `Matrix{Float64}`: covariance matrix, shape `(channels, channels)`
"""
function scovm(s::AbstractMatrix; corrected::Bool = false)::Matrix{Float64}
    # transpose so columns are channels
    cm = cov(SimpleCovariance(corrected = corrected), s')

    return cm
end

"""
    scovm(s; <keyword arguments>)

Calculate sample covariance matrix of a 3-D signal array (S = κ (Xc' * Xc)).

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `corrected::Bool=false`: if `true`, calculate corrected covariance (κ = number of samples, otherwise κ = number of samples - 1)

# Returns

- `Array{Float64, 3}`: covariance matrix, shape `(channels, channels, epochs)`
"""
function scovm(s::AbstractArray; corrected::Bool = false)::Array{Float64, 3}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    cm = zeros(ch_n, ch_n, ep_n)

    # calculate over epochs
    @inbounds Threads.@threads :static for ep_idx in 1:ep_n
        cm[:, :, ep_idx] = scovm(@view(s[:, :, ep_idx]), corrected = corrected)
    end

    return cm
end

"""
    scovm(obj; <keyword arguments>)

Calculate sample covariance matrix between all channel pairs within a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}: channel name(s)
- `corrected::Bool=false`: if `true`, calculate corrected covariance matrix

# Returns

- `Array{Float64, 3}`: covariance matrix, shape `(channels, channels, epochs)`
"""
function scovm(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    corrected::Bool = false,
)::Array{Float64, 3}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return scovm(@view(obj.data[ch, :, :]); corrected = corrected)
end
