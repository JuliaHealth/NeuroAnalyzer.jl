export rmse

"""
    rmse(s1, s2)

Calculate Root Mean Square Error (RMSE).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

- `rmse::Float64`: RMSE
"""
function rmse(s1::AbstractVector, s2::AbstractVector)

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # r = sum(s1 .* s2) ./ (sqrt(sum(s1.^2)) .* sqrt(sum(s2.^2)))
    return sqrt(mean(s2 - s1)^2)

end

"""
    rmse(s1, s2)

Calculate Root Mean Square Error (RMSE).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`

# Returns

- `r::Array{Float64, 2}`: RMSE
"""
function rmse(s1::AbstractArray, s2::AbstractArray)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)
    
    r = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            r[ch_idx, ep_idx] = @views rmse(s1[ch_idx, :, ep_idx], s1[ch_idx, :, ep_idx])
        end
    end

    return r

end

"""
    rmse(obj1, obj2; ch1, ch2, ep1, ep2)

Calculate Root Mean Square Error (RMSE).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `r::Array{Float64, 3}`: RMSE
- `cps_ph::Array{Float64, 3}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64, 3}`: cross power spectrum frequencies
"""
function rmse(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2)))

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."
    
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    length(ch1) == 1 && (ch1 = [ch1])
    length(ch2) == 1 && (ch2 = [ch2])

    r = @views rmse(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return r

end
