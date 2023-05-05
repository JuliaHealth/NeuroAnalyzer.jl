export xcor

"""
   xcor(s1, s2; l, demean)

Calculate cross-correlation (a measure of similarity of two signals as a function of the displacement of one relative to the other).

# Arguments

- `s1::AbstractMatrix`
- `s2::AbstractMatrix`
- `l::Int64=round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2))))`
- `demean::Bool=true`: demean signal before computing cross-correlation

# Returns

- `xc::Array{Float64, 3}`
"""
function xcor(s1::AbstractMatrix, s2::AbstractMatrix; l::Int64=round(Int64, min(size(s1, 1), 10 * log10(size(s1, 1)))), demean::Bool=true)

    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))

    ep_n = size(s1, 2)

    xc = zeros(1, length(-l:l), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            xc[1, :, ep_idx] = @views crosscor(s1[1, :, ep_idx], s2[1, :, ep_idx], -l:l, demean=demean)
        end
    end

    return xc

end

"""
   xcor(s1, s2; l, demean)

Calculate cross-correlation (a measure of similarity of two signals as a function of the displacement of one relative to the other).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `l::Int64=round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2))))`
- `demean::Bool=true`: demean signal before computing cross-correlation

# Returns

- `xc::Array{Float64, 3}`
"""
function xcor(s1::AbstractArray, s2::AbstractArray; l::Int64=round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2)))), demean::Bool=true)

    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    xc = zeros(ch_n, length(-l:l), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            xc[ch_idx, :, ep_idx] = @views crosscor(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], -l:l, demean=demean)
        end
    end

    return xc

end

"""
    xcor(obj1, obj2; ch1, ch2, ep1, ep2, lag, norm)

Calculate cross-correlation (a measure of similarity of two signals as a function of the displacement of one relative to the other).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `lag::Real=1`: lags range is `-lag:lag` [s]
- `demean::Bool=true`: demean signal before computing cross-correlation

# Returns

Named tuple containing:
- `xc::Array{Float64, 3}`: cross-correlation
- `lag::Vector{Float64}`: lags [s]
"""
function xcor(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)), lag::Real=1, demean::Bool=true)

    # check channels
    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    length(ch1) == length(ch2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    lag > obj1.epoch_time[end] && throw(ArgumentError("lag must be ≤ $(obj1.epoch_time[end])."))

    l = vsearch(lag, obj1.epoch_time)
    xc = @views xcor(reshape(obj1.data[ch1, :, ep1], length(ch1), :, length(ep1)), reshape(obj2.data[ch2, :, ep2], length(ch2), :, length(ep2)), l=l, demean=demean)

    return (xc=xc, lag=vcat(-reverse(obj1.epoch_time[1:l]), 0, obj1.epoch_time[1:l]))

end
