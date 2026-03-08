export frqinst

"""
    frqinst(s; <keyword arguments>)

Estimate the instantaneous frequency of a signal via the Hilbert transform:

1. Compute the analytic signal with htransform().
2. Extract the instantaneous phase (index 3 of htransform output).
3. Unwrap the phase to remove 2π discontinuities.
4. Differentiate to obtain dφ/dt in radians per sample.
5. Divide by 2π to convert to cycles per sample.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `f::Vector{Float64}`: instantaneous frequencies in cycles per sample (multiply by the sampling rate to obtain Hz)

# Notes

Uses the Hilbert transform; best results for narrowband signals. Broadband signals produce meaningless instantaneous frequencies.
"""
function frqinst(s::AbstractVector)::Vector{Float64}

    h = htransform(s)
    ph = h.ph
    f = NeuroAnalyzer.derivative(DSP.unwrap(ph)) / (2 * π)

    return f

end

"""
    frqinst(s; <keyword arguments>)

Estimate the instantaneous frequency of a signal via the Hilbert transform:

1. Compute the analytic signal with htransform().
2. Extract the instantaneous phase (index 3 of htransform output).
3. Unwrap the phase to remove 2π discontinuities.
4. Differentiate to obtain dφ/dt in radians per sample.
5. Divide by 2π to convert to cycles per sample.

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)

# Returns

- `f::Array{Float64, 3}`: instantaneous frequencies in cycles per sample (multiply by the sampling rate to obtain Hz), shape `(channels, samples, epochs)`

# Notes

Uses the Hilbert transform; best results for narrowband signals. Broadband signals produce meaningless instantaneous frequencies.
"""
function frqinst(s::AbstractArray)::Array{Float64, 3}

    _warn("frqinst() uses Hilbert transform, the signal should be narrowband for best results.")

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    f = similar(s)
    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        f[ch_idx, :, ep_idx] = frqinst(@view(s[ch_idx, :, ep_idx]))
    end

    return f

end

"""
    frqinst(obj; <keyword arguments>)

Estimate the instantaneous frequency of a signal via the Hilbert transform:

1. Compute the analytic signal with htransform().
2. Extract the instantaneous phase (index 3 of htransform output).
3. Unwrap the phase to remove 2π discontinuities.
4. Differentiate to obtain dφ/dt in radians per sample.
5. Divide by 2π to convert to cycles per sample.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `f::Array{Float64, 3}`: instantaneous frequencies in Hz, shape `(channels, samples, epochs)`

# Notes

Uses the Hilbert transform; best results for narrowband signals. Broadband signals produce meaningless instantaneous frequencies.
"""
function frqinst(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    f = frqinst(@view(obj.data[ch, :, :])) .* sr(obj)

    return f

end
