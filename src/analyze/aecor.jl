export aecor
export escor

"""
    aecor(s1, s2)

Calculate Amplitude Envelope Correlation (AEC).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

- `aec::Float64`: AEC value

# Source

1. Bruns, A., Eckhorn, R., Jokeit, H., & Ebner, A. (2000). Amplitude envelope correlation detects coupling among incoherent brain signals. Neuroreport, 11(7), 1509-1514.
"""
function aecor(s1::AbstractVector, s2::AbstractVector)::Float64

    @assert length(s1) == length(s2) "Both signals must have the same length."

    # envelopes made of instantaneous amplitudes
    _, e1, _, _ = htransform(s1)
    _, e2, _, _ = htransform(s2)

    aec = cor(e1, e2)

    return aec

end

"""
    aecor(obj1, obj2; <keyword arguments>)

Calculate Amplitude Envelope Correlation (AEC).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

- `aec::Matrix{Float64}`: AEC value
"""
function aecor(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))::Matrix{Float64}

    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1)) and ep2 ($(length(ep2)) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    ep_n = length(ep1)
    ch_n = length(ch1)

    aec = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            aec[ch_idx, ep_idx] = @views aecor(obj1.data[ch1[ch_idx], :, ep1[ep_idx]], obj2.data[ch2[ch_idx], :, ep2[ep_idx]])
        end
    end

    return aec

end

"""
    aecor(obj; <keyword arguments>)

Calculate Amplitude Envelope Correlation (AEC).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Named tuple containing:
- `aec::Array{Float64, 3}`: AEC value
"""
function aecor(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    isa(ch, Int64) && (ch = [ch])

    aec = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                aec[ch_idx1, ch_idx2, ep_idx] = @views aecor(obj.data[ch[ch_idx1], :, ep_idx], obj.data[ch[ch_idx2], :, ep_idx])
            end
        end
    end

    aec = _copy_lt2ut(aec)

    return aec

end

"""
    escor(s1, s2)

Calculate Envelope-to-Signal Correlation (ESC).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

- `esc::Float64`: ESC value

# Source

1. Bruns, A., & Eckhorn, R. (2004). Task-related coupling from high-to low-frequency signals among visual cortical areas in human subdural recordings. International Journal of Psychophysiology, 51(2), 97-116.
"""
function escor(s1::AbstractVector, s2::AbstractVector)::Float64

    @assert length(s1) == length(s2) "Both signals must have the same length."

    # envelopes made of instantaneous amplitudes
    _, e2, _, _ = htransform(s2)

    esc = cor(s1, e2)

    return esc

end

"""
    escor(obj1, obj2; <keyword arguments>)

Calculate Envelope-to-Signal Correlation (ESC).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

- `aec::Matrix{Float64}`: ESC value
"""
function escor(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))::Matrix{Float64}

    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1)) and ep2 ($(length(ep2)) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    ep_n = length(ep1)
    ch_n = length(ch1)

    esc = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            esc[ch_idx, ep_idx] = @views escor(obj1.data[ch1[ch_idx], :, ep1[ep_idx]], obj2.data[ch2[ch_idx], :, ep2[ep_idx]])
        end
    end

    return esc

end

"""
    escor(obj; <keyword arguments>)

Calculate Envelope-to-Signal Correlation (ESC).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Named tuple containing:
- `pv::Array{Float64, 3}`: ESC value
"""
function escor(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    isa(ch, Int64) && (ch = [ch])

    esc = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                esc[ch_idx1, ch_idx2, ep_idx] = @views escor(obj.data[ch[ch_idx1], :, ep_idx], obj.data[ch[ch_idx2], :, ep_idx])
            end
        end
    end

    esc = _copy_lt2ut(esc)

    return esc

end
