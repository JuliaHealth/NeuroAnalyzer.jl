export xcorr

"""
   xcorr(s1, s2)

Calculate cross-correlation.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

- `xc::Vector{Float64}`
"""
function xcorr(s1::AbstractVector, s2::AbstractVector)

    xc = @views DSP.xcorr(s1, s2)

    return xc

end


"""
   xcorr(s1, s2)

Calculate cross-correlation.

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`

# Returns

- `xc::Array{Float64, 3}`
"""
function xcorr(s1::AbstractArray, s2::AbstractArray)

    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    xc = zeros(ch_n, (2 * size(s1, 2) - 1), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            xc[ch_idx, :, ep_idx] = @views DSP.xcorr(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx])
        end
    end

    return xc

end

"""
    xcorr(obj1, obj2; ch1, ch2, ep1, ep2)

Calculate cross correlation.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

- `xc::Array{Float64, 3}`
"""
function xcorr(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj2)))

    # check channels
    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    length(ch1) == length(ch2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    xc = @views NeuroAnalyzer.xcorr(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return xc

end

"""
   xcorr(s)

Calculate cross-correlation.

# Arguments

- `s::AbstractArray`

# Returns

- `xc::Array{Float64, 3}`
"""
function xcorr(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    xc = zeros(ch_n, ch_n, (2 * size(s, 2) - 1), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_n
                xc[ch_idx1, ch_idx2, :, ep_idx] = @views DSP.xcorr(s[ch_idx1, :, ep_idx], s[ch_idx2, :, ep_idx])
            end
        end
    end

    return xc

end

"""
   xcorr(obj; ch, lag, norm)

Calculate cross-correlation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `xc::Array{Float64, 3}`
"""
function xcorr(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), lag::Int64=1, norm::Bool=false)

    _check_channels(obj, ch)

    xc = @views NeuroAnalyzer.xcorr(obj.data[ch, :, :])

    return xc

end
