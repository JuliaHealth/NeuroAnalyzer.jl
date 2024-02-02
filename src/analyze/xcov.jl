export xcov

"""
   xcov(s1, s2; l, demean, n, biased)

Calculate cross-covariance.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing cross-covariance
- `biased::Bool=true`: calculate biased or unbiased cross-covariance

# Returns

- `xc::Matrix{Float64}`
"""
function xcov(s1::AbstractVector, s2::AbstractVector; l::Int64=round(Int64, min(length(s) - 1, 10 * log10(length(s)))), demean::Bool=true, biased::Bool=true)

    @assert length(s1) == length(s2) "Both signals must have the same length."

    xc = zeros(l + 1)
    xc_neg = zeros(l + 1)

    ms1 = mean(s1)
    ms2 = mean(s2)
    if demean == true
        s1_tmp = s1 .- ms1
        s2_tmp = s2 .- ms2
    else
        s1_tmp = s1
        s2_tmp = s2
    end

    for idx in 0:l
        xc[idx + 1] = @views sum(s1_tmp[(1 + idx):end] .* s2_tmp[1:(end - idx)])
        if biased == true 
            xc[idx + 1] /= length(s1)
        else
            xc[idx + 1] /= (length(s1) - idx)
        end
    end
    
    for idx in 0:l
        xc_neg[idx + 1] = @views sum(s1_tmp[1:(end - idx)] .* s2_tmp[(1 + idx):end])
        if biased == true 
            xc_neg[idx + 1] /= length(s1)
        else
            xc_neg[idx + 1] /= (length(s1) - idx)
        end
    end

    xc = vcat(reverse(xc_neg), xc[2:end])
    xc = round.(xc, digits=3)

    return reshape(xc, 1, :, 1)

end

"""
   xcov(s1, s2; l, demean)

Calculate cross-covariance.

# Arguments

- `s1::AbstractMatrix`
- `s2::AbstractMatrix`
- `l::Int64=round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing cross-covariance
- `biased::Bool=true`: calculate biased or unbiased cross-covariance

# Returns

- `xc::Array{Float64, 3}`
"""
function xcov(s1::AbstractMatrix, s2::AbstractMatrix; l::Int64=round(Int64, min(size(s1, 1), 10 * log10(size(s1, 1)))), demean::Bool=true, biased::Bool=true)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ep_n = size(s1, 2)

    xc = zeros(1, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        xc[1, :, ep_idx] = @views xcov(s1[1, :, ep_idx], s2[1, :, ep_idx], l=l, demean=demean, biased=biased)
    end

    return xc

end

"""
   xcov(s1, s2; l, demean)

Calculate cross-covariance.

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `l::Int64=round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing cross-covariance
- `biased::Bool=true`: calculate biased or unbiased cross-covariance

# Returns

- `xc::Array{Float64, 3}`
"""
function xcov(s1::AbstractArray, s2::AbstractArray; l::Int64=round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2)))), demean::Bool=true, biased::Bool=true)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    xc = zeros(ch_n, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            xc[ch_idx, :, ep_idx] = @views xcov(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], l=l, demean=demean, biased=biased)
        end
    end

    return xc

end

"""
    xcov(obj1, obj2; ch1, ch2, ep1, ep2, l, norm)

Calculate cross-covariance. For ERP return trial-averaged cross-covariance.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `l::Real=1`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing cross-covariance
- `biased::Bool=true`: calculate biased or unbiased cross-covariance

# Returns

Named tuple containing:
- `xc::Array{Float64, 3}`: cross-covariance
- `l::Vector{Float64}`: lags [s]
"""
function xcov(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)), l::Real=1, demean::Bool=true, biased::Bool=true)

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."

    # check channels
    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."
    
    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    @assert l <= size(obj1, 2) "l must be ≤ $(size(obj1, 2))."
    @assert l >= 0 "l must be ≥ 0."

    if obj1.header.recording[:data_type] == "erp" && obj2.header.recording[:data_type] == "erp"
        xc = @views xcov(reshape(obj1.data[ch1, :, 2:end], length(ch1), :, (nepochs(obj1) - 1)), reshape(obj2.data[ch2, :, 2:end], length(ch2), :, (nepochs(obj2) - 1)), l=l, demean=demean, biased=biased)
        xc = cat(mean(xc, dims=3), xc, dims=3)
    else
        xc = @views xcov(reshape(obj1.data[ch1, :, ep1], length(ch1), :, length(ep1)), reshape(obj2.data[ch2, :, ep2], length(ch2), :, length(ep2)), l=l, demean=demean, biased=biased)
    end

    return (xc=xc, l=round.(collect(-l:l) .* (1/sr(obj1)), digits=3))

end
