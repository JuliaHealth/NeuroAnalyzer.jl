export aecor
export escor

"""
    aecor(s1, s2)

Calculate Amplitude Envelope Correlation (AEC).

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

- `aec::Float64`: AEC value

# References

1. Bruns, A., Eckhorn, R., Jokeit, H., & Ebner, A. (2000). Amplitude envelope correlation detects coupling among incoherent brain signals. Neuroreport, 11(7), 1509-1514.
"""
function aecor(s1::AbstractVector, s2::AbstractVector)::Float64

    @assert length(s1) == length(s2) "Both signals must have the same length."

    # instantaneous amplitude envelopes via Hilbert transform
    ht1 = htransform(s1)
    ht2 = htransform(s2)
    e1 = ht1.a
    e2 = ht2.a

    # AEC is the Pearson correlation of the two amplitude envelopes
    aec = cor(e1, e2)

    return aec

end

"""
    aecor(obj1, obj2; <keyword arguments>)

Calculate Amplitude Envelope Correlation (AEC).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

- `aec::Matrix{Float64}`: AEC value, shape `(channels, epochs)`
"""
function aecor(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    # number of channels
    ch_n = length(ch1)
    # number of epochs
    ep_n = length(ep1)

    # pre-allocate output
    aec = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        aec[ch_idx, ep_idx] = aecor(
            @view(obj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(obj2.data[ch2[ch_idx], :, ep2[ep_idx]])
        )
    end

    return aec

end

"""
    aecor(obj; <keyword arguments>)

Calculate Amplitude Envelope Correlation (AEC).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `aec::Array{Float64, 3}`: AEC value, shape `(channels, channels, epochs)`
"""
function aecor(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex}
)::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch1)
    # number of epochs
    ep_n = nepochs(obj)

    # pre-allocate output
    aec = zeros(ch_n, ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 in 1:ch_idx1
            aec[ch_idx1, ch_idx2, ep_idx] = aecor(
                @view(obj.data[ch[ch_idx1], :, ep_idx]),
                @view(obj.data[ch[ch_idx2], :, ep_idx])
            )
        end
    end

    # mirror the lower triangle to the upper triangle to produce the full
    # symmetric matrix
    aec = _copy_lt2ut(aec)

    return aec

end

"""
    escor(s1, s2)

Calculate Envelope-to-Signal Correlation (ESC).

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

- `esc::Float64`: ESC value

# References

Bruns, A., & Eckhorn, R. (2004). Task-related coupling from high-to low-frequency signals among visual cortical areas in human subdural recordings. International Journal of Psychophysiology, 51(2), 97-116.
"""
function escor(s1::AbstractVector, s2::AbstractVector)::Float64

    @assert length(s1) == length(s2) "Both signals must have the same length."

    # instantaneous amplitude envelope via Hilbert transform
    # only s2's envelope is needed; s1 enters the correlation as the raw signal
    ht2 = htransform(s2)
    e2 = ht2.a

    # ESC: correlation of the raw s1 signal against s2's amplitude envelope
    esc = cor(s1, e2)

    return esc

end

"""
    escor(obj1, obj2; <keyword arguments>)

Calculate Envelope-to-Signal Correlation (ESC).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

- `esc::Matrix{Float64}`: ESC value, shape `(channels, epochs)`
"""
function escor(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    # number of channels
    ch_n = length(ch1)
    # number of epochs
    ep_n = length(ep1)

    # pre-allocate output
    esc = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        esc[ch_idx, ep_idx] = @views escor(
            obj1.data[ch1[ch_idx], :, ep1[ep_idx]],
            obj2.data[ch2[ch_idx], :, ep2[ep_idx]]
        )
    end

    return esc

end

"""
    escor(obj; <keyword arguments>)

Calculate Envelope-to-Signal Correlation (ESC).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `esc::Array{Float64, 3}`: ESC value, shape `(channels, channels, epochs)`
"""
function escor(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    
    # number of channels
    ch_n = length(ch1)
    # number of epochs
    ep_n = nepochs(obj)

    # pre-allocate output
    esc = zeros(ch_n, ch_n, ep_n)

    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 in 1:ch_idx1
            esc[ch_idx1, ch_idx2, ep_idx] = escor(
                @view(obj.data[ch[ch_idx1], :, ep_idx]),
                @view(obj.data[ch[ch_idx2], :, ep_idx])
            )
        end
    end

    # mirror the lower triangle to the upper triangle to produce the full symmetric matrix
    esc = _copy_lt2ut(esc)

    return esc

end
