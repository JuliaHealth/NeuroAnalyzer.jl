export xcov

"""
   xcov(s1, s2; lag, norm)

Calculate cross-covariance.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `lag::Int64`: lags range is `-lag:lag`
- `norm::Bool`: normalize cross-covariance

# Returns

Named tuple containing:
- `xc::Vector{Float64}`: cross-covariance
- `l::Vector{Int64}`: lags
"""
function xcov(s1::AbstractVector, s2::AbstractVector; lag::Int64=1, norm::Bool=false)

    length(s1) == length(s2) || throw(ArgumentError("Both signals must be of the same as length."))
    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

    l = collect(-lag:lag)

    xc = zeros(length(l))
    
    @inbounds @fastmath @simd for idx in eachindex(l)
        # no lag
        l[idx] == 0 && (xc[idx] = sum(s1 .* s2))
        # positive lag
        l[idx] > 0 && (xc[idx] = @views sum(s1[(1 + l[idx]):end] .* s2[1:(end - l[idx])]))
        # negative lag
        l[idx] < 0 && (xc[idx] = @views sum(s1[1:(end - abs(l[idx]))] .* s2[(1 + abs(l[idx])):end]))
    end
    norm == true && (xc ./ length(s1))

    return (xc=xc, l=l)
    
end

"""
    xcov(s1, s2; lag, norm)

Calculate cross-covariance.

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `lag::Int64=1`: lags range is `-lag:lag`
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `xc::Array{Float64, 3}`: cross-covariance
- `l::Vector{Float64}`: lags
"""
function xcov(s1::AbstractArray, s2::AbstractArray; lag::Int64=1, norm::Bool=false)

    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    # create vector of lags
    l = collect(-lag:lag)

    xc = zeros(ch_n, (2 * lag + 1), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            xc[ch_idx, :, ep_idx], _ = @views xcov(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], lag=lag, norm=norm)
        end
    end

    return (xc=xc, l=l)

end


"""
    xcov(obj; ch, lag, norm)

Calculate cross-covariance.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `lag::Int64=1`: lags range is `-lag:lag`
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `xc::Matrix{Float64}`: cross-covariance
- `l::Vector{Float64}`: lags in ms
"""
function xcov(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), lag::Int64=1, norm::Bool=false)

    #_check_channels(obj, ch)
    ch_n = length(ch)
    ep_n = epoch_n(obj)

    # create vector of lags
    l = 1/sr(obj) .* collect(-lag:lag) .* 1000
    xc = zeros(ch_n^2, length(l), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n

        # create half of the cross-covariance matrix
        xcov_packed = Array{Vector{Float64}}(undef, ch_n, ch_n)
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                xcov_packed[ch_idx1, ch_idx2], _ = @views xcov(obj.data[ch[ch_idx1], :, ep_idx], obj.data[ch[ch_idx2], :, ep_idx], lag=lag, norm=norm)
            end
        end
        
        # copy to the other half
        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                xcov_packed[ch_idx1, ch_idx2] = @views xcov_packed[ch_idx2, ch_idx1]
            end
        end

        # unpack by channels
        Threads.@threads for ch_idx in 1:ch_n^2
            xc[ch_idx, :, ep_idx] = @views xcov_packed[ch_idx]
        end
    end

    return (xc=xc, l=l)

end

"""
    xcov(obj1, obj2; ch1, ch2, ep1, ep2, lag, norm)

Calculate cross-covariance between two NeuroAnalyzer NEURO objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `lag::Int64=1`: lags range is `-lag:lag`
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `xc::Array{Float64, 3}`: cross-covariance
- `l::Vector{Float64}`: lags in ms
"""
function xcov(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)), lag::Int64=1, norm::Bool=false)

    # check channels
    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    length(ch1) == length(ch2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    # create vector of lags
    l = 1/sr(obj1) .* collect(-lag:lag) .* 1000

    xc, _ = @views xcov(reshape(obj1.data[ch1, :, ep1], length(ch1), :, length(ep1)), reshape(obj2.data[ch2, :, ep2], length(ch2), :, length(ep2)), lag=lag, norm=norm)

    return (xc=xc, l=l)

end
