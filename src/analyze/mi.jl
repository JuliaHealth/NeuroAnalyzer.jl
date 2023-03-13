export mutual_information

"""
    mutual_information(s1, s2)

Calculate mutual information.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

- `mutual_information::Float64`
"""
function mutual_information(s1::AbstractVector, s2::AbstractVector)

    return get_mutual_information(s1, s2)

end

"""
    mutual_information(s1, s2)

Calculate mutual information (channels of `s1` vs channels of `s2`).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`

# Returns

- `mutual_information::Array{Float64}`
"""
function mutual_information(s1::AbstractArray, s2::AbstractArray)

    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    m = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            m[ch_idx, ep_idx] = @views mutual_information(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx])
        end
    end

    return m

end

"""
    mutual_information(s)

Calculate mutual information (channels vs channels).

# Arguments

- `s::AbstractArray`

# Returns

"""
function mutual_information(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    
    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    m = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
           for ch_idx2 in 1:ch_idx1
                m[ch_idx1, ch_idx2, ep_idx] = @views mutual_information(s[ch_idx1, :, ep_idx], s[ch_idx2, :, ep_idx])
            end

        # update progress bar
        progress_bar == true && next!(p)
        end
    end

    # copy lower triangle to upper triangle
    Threads.@threads for ep_idx in 1:ep_n
        @inbounds m[:, :, ep_idx] = _copy_lt2ut(m[:, :, ep_idx])
    end

    return m

end

"""
    mutual_information(obj; channel)

Calculate mutual information between channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `mutual_information::Array{Float64, 3}`
"""
function mutual_information(obj::NeuroAnalyzer.NEURO; channel::Union{Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    m = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n

        # create half of the matrix
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                m[ch_idx1, ch_idx2, ep_idx] = @views mutual_information(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx])
            end
        end

        # copy to the other half
        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                m[ch_idx1, ch_idx2, ep_idx] = @views m[ch_idx2, ch_idx1, ep_idx]
            end
        end
    end

    return m

end

"""
    mutual_information(obj1, obj2; ch1, ch2, ep1, ep2)

Calculate mutual information between two channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

- `m::Array{Float64, 3}`
"""
function mutual_information(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj2)))

    # check channels
    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    length(ch1) == length(ch2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    m = @views mutual_information(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    ch_n = length(ch1)
    ep_n = length(ep1)

    m = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_n
                m[ch_idx1, ch_idx2, ep_idx] = @views mutual_information(obj1.data[ch1[ch_idx1], :, ep1[ep_idx]], obj2.data[ch2[ch_idx2], :, ep2[ep_idx]])
            end
        end
    end

    return m

end