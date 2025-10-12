export rms
export rmse

"""
    rms(s)

Calculate Root Mean Square (RMS).

# Arguments

- `s::AbstractVector`

# Returns

- `r::Float64`: RMS
"""
function rms(s::AbstractVector)::Float64

    # r = sqrt(1/length(s) * sum(s.^2))
    # r = sqrt(mean(s.^2))
    r = norm(s) / sqrt(length(s))

    return r

end

"""
    rms(s)

Calculate Root Mean Square (RMS).

# Arguments

- `s::AbstractArray`

# Returns

- `r::Matrix{Float64}`: RMS
"""
function rms(s::AbstractArray)::Matrix{Float64}

    _chk3d(s)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    r = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            r[ch_idx, ep_idx] = @views rms(s[ch_idx, :, ep_idx])
        end
    end

    return r

end

"""
    rms(obj; <keyword arguments>)

Calculate Root Mean Square (RMS).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}: list of channels
- `ep::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj))`: default use all epochs

# Returns

- `r::Matrix{Float64}`: RMS
"""
function rms(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, ep::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj)))::Matrix{Float64}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")

    _check_epochs(obj, ep)
    isa(ep, Int64) && (ep = [ep])

    r = @views rms(obj.data[ch, :, ep])

    return r

end

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
    r = sqrt(mean((s2 .- s1).^2))

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
        Threads.@threads for ch_idx in 1:ch_n
            r[ch_idx, ep_idx] = @views rmse(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx])
        end
    end

    return r

end

"""
    rmse(obj1, obj2; <keyword arguments>)

Calculate Root Mean Square Error (RMSE).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

- `r::Matrix{Float64}`: RMSE
"""
function rmse(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))::Matrix{Float64}

    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal."

    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1)) and ep2 ($(length(ep2)) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    r = @views rmse(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return r

end
