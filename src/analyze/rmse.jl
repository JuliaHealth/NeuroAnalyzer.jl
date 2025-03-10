export rmse

"""
    rmse(s1, s2)

Calculate Root Mean Square Error (RMSE).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

- `r::Float64`: RMSE
"""
function rmse(s1::AbstractVector, s2::AbstractVector)::Float64

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # r = sum(s1 .* s2) ./ (sqrt(sum(s1.^2)) .* sqrt(sum(s2.^2)))
    r = sqrt(mean(s2 - s1)^2)

    return r

end

"""
    rmse(s1, s2)

Calculate Root Mean Square Error (RMSE).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`

# Returns

- `r::Matrix{Float64}`: RMSE
"""
function rmse(s1::AbstractArray, s2::AbstractArray)::Matrix{Float64}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _chk3d(s1)
    _chk3d(s2)

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    r = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            r[ch_idx, ep_idx] = @views rmse(s1[ch_idx, :, ep_idx], s1[ch_idx, :, ep_idx])
        end
    end

    return r

end

"""
    rmse(obj1, obj2; <keyword arguments>)

Calculate Root Mean Square Error (RMSE).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

- `r::Matrix{Float64}`: RMSE
"""
function rmse(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))::Matrix{Float64}

    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    r = @views rmse(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return r

end
