export sumsim

"""
    sumsim(s1, s2; theta)

Calculate summed similarity using an exponential decay model between two signals.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `theta::Real`: decay parameter

# Returns

- `ss::Float64`: summed similarity

# Notes

Values of `ss` are in the range [0, 1]; higher value indicates larger similarity.
"""
function sumsim(s1::AbstractVector, s2::AbstractVector; theta::Real)::Float64

    @assert length(s1) == length(s2) "Lengths of s1 ($(length(s1))) and s2 ($(length(s2))) must be equal."
    ss = exp(-theta * sqrt(sum((s1 .- s2).^2)))

    return ss

end

"""
    sumsim(s1, s2; theta)

Calculate summed similarity using an exponential decay model between two signals.

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `theta::Real`: decay parameter

# Returns

- `ss::Matrix{Float64}`: signal symmetry

# Notes

Values of `ss` are in the range [0, 1]; higher value indicates larger similarity.
"""
function sumsim(s1::AbstractArray, s2::AbstractArray; theta::Real)::Matrix{Float64}

    @assert size(s1) == size(s2) "Sizes of s1 ($(size(s1))) and s2 ($(size(s2))) must be equal."
    _chk3d(s1)
    _chk3d(s2)

    ch_n = size(s1, 1)
    ep_n = size(s2, 3)

    ss = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ss[ch_idx, ep_idx] = @views sumsim(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], theta=theta)
        end
    end

    return ss

end

"""
    sumsim(obj; <keyword arguments>)

Calculate signal symmetry (ratio of positive to negative amplitudes). Perfectly symmetrical signal has symmetry of 1.0. Symmetry above 1.0 indicates there are more positive amplitudes.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}`: list of channels
- `ch2::Union{String, Vector{String}}`: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `theta::Real`: decay parameter

# Returns

- `ss::Matrix{Float64}`: signal symmetry

# Notes

Values of `ss` are in the range [0, 1]; higher value indicates larger similarity.
"""
function sumsim(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)), theta::Real)::Matrix{Float64}

    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal."
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1)) and ep2 ($(length(ep2)) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    ss = @views sumsim(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2], theta=theta)

    return ss

end
