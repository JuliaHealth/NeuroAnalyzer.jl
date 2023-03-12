export vartest

"""
    vartest(obj; channel)

Calculate variance F-test.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function vartest(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    f = zeros(ch_n, ch_n, ep_n)
    p = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
       Threads.@threads for ch_idx1 in 1:ch_n
            # create half of the matrix
            for ch_idx2 in 1:ch_idx1
                ftest = @views VarianceFTest(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx])
                f[ch_idx1, ch_idx2, ep_idx] = ftest.F
                p[ch_idx1, ch_idx2, ep_idx] = pvalue(ftest)
            end
        end
        # copy to the other half
        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                f[ch_idx1, ch_idx2, ep_idx] = @views f[ch_idx2, ch_idx1, ep_idx]
                p[ch_idx1, ch_idx2, ep_idx] = @views p[ch_idx2, ch_idx1, ep_idx]
            end
        end
    end

    return (f=f, p=p)
end

"""
    vartest(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculate variance F-test.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function vartest(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), channel2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    f = zeros(ch_n, ch_n, ep_n)
    p = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
       Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_n
                ftest = @views VarianceFTest(obj1.data[channel1[ch_idx1], :, epoch1[ep_idx]], obj2.data[channel2[ch_idx2], :, epoch2[ep_idx]])
                f[ch_idx1, ch_idx2, ep_idx] = ftest.F
                p[ch_idx1, ch_idx2, ep_idx] = pvalue(ftest)
            end
        end
    end

    return (f=f, p=p)
end
