export resample
export resample!
export upsample
export upsample!
export downsample
export downsample!

"""
    resample(s; <keyword arguments>)

Resample to `new_sr` sampling frequency.

# Arguments

- `s::AbstractVector`: signal vector
- `old_sr::Int64`: old sampling rate
- `new_sr::Int64`: new sampling rate

# Returns

- `Vector{Float64}`
"""
function resample(s::AbstractVector; old_sr::Int64, new_sr::Int64)::Vector{Float64}
    # validate
    old_sr >= 1 || throw(ArgumentError("old_sr must be ≥ 1."))
    new_sr >= 1 || throw(ArgumentError("new_sr must be ≥ 1."))

    new_sr == old_sr && return (s)

    # resample
    sr_ratio = new_sr / old_sr
    s_new = DSP.resample(s, sr_ratio)

    return s_new
end

"""
    resample(s; <keyword arguments>)

Resamples all channels and time vector `t` to `new_sr` sampling frequency for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `old_sr::Int64`: old sampling rate
- `new_sr::Int64`: new sampling rate

# Returns

- `Array{Float64, 3}`
"""
function resample(s::AbstractArray; old_sr::Int64, new_sr::Int64)::Array{Float64, 3}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # validate
    new_sr >= 1 || throw(ArgumentError("new_sr must be ≥ 1."))

    # number of channels and epochs
    ch_n, _, ep_n = size(s)

    # dry run
    s_tmp = NeuroAnalyzer.resample(s[1, :, 1]; old_sr = old_sr, new_sr = new_sr)

    # pre-allocate output
    s_new = zeros(ch_n, length(s_tmp), ep_n)

    # calculate over channels and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = NeuroAnalyzer.resample(
            @view(s[ch_idx, :, ep_idx]),
            old_sr = old_sr,
            new_sr = new_sr,
        )
    end

    return s_new
end

"""
    resample(obj; <keyword arguments>)

Resample all channels to `new_sr` sampling frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `new_sr::Int64`: new sampling rate

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function resample(obj::NeuroAnalyzer.NEURO; new_sr::Int64)::NeuroAnalyzer.NEURO
    # validate
    new_sr >= 1 || throw(ArgumentError("new_sr must be ≥ 1."))

    if new_sr > sr(obj)
        return upsample(obj; new_sr = new_sr)
    elseif new_sr < sr(obj)
        return downsample(obj; new_sr = new_sr)
    else
        return obj
    end
end

"""
    resample!(obj; <keyword arguments>)

Resample all channels to `new_sr` sampling frequency in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `new_sr::Int64`: new sampling rate

# Returns

- `Nothing`
"""
function resample!(obj::NeuroAnalyzer.NEURO; new_sr::Int64)::Nothing
    obj_new = resample(obj; new_sr = new_sr)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing
end

"""
    upsample(obj; <keyword arguments>)

Upsample all channels to `new_sr` sampling frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `new_sr::Int64`: new sampling rate

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function upsample(obj::NeuroAnalyzer.NEURO; new_sr::Int64)::NeuroAnalyzer.NEURO
    # validate
    new_sr / sr(obj) != new_sr ÷ sr(obj) && _warn(
        "New sampling rate should be easily captured by integer fractions, e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.",
    )

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data = NeuroAnalyzer.resample(obj.data; old_sr = sr(obj), new_sr = new_sr)
    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    obj_new.header.recording[:sampling_rate] = new_sr
    push!(obj_new.history, "upsample(obj, new_sr=$new_sr)")

    return obj_new
end

"""
    upsample!(obj; <keyword arguments>)

Upsample all channels to `new_sr` sampling frequency in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `new_sr::Int64`: new sampling rate

# Returns

- `Nothing`
"""
function upsample!(obj::NeuroAnalyzer.NEURO; new_sr::Int64)::Nothing
    obj_new = upsample(obj; new_sr = new_sr)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing
end

"""
    downsample(obj; <keyword arguments>)

Downsample all channels to `new_sr` sampling frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `new_sr::Int64`: new sampling rate

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function downsample(obj::NeuroAnalyzer.NEURO; new_sr::Int64)::NeuroAnalyzer.NEURO
    # validate
    new_sr < sr(obj) && _warn(
        "To prevent aliasing due to down-sampling, a low-pass filter should be applied before removing data points. The filter cutoff should be the Nyquist frequency of the new down-sampled rate, ($(new_sr / 2) Hz), not the original Nyquist frequency ($(sr(obj) / 2) Hz).",
    )

    new_sr / sr(obj) != new_sr ÷ sr(obj) && _warn(
        "New sampling rate should be easily captured by integer fractions e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.",
    )

    # create new dataset
    obj_new = deepcopy(obj)
    s_new = NeuroAnalyzer.resample(obj.data; old_sr = sr(obj), new_sr = new_sr)

    obj_new.data = s_new

    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    obj_new.header.recording[:sampling_rate] = new_sr
    push!(obj_new.history, "downsample(obj, new_sr=$new_sr)")

    return obj_new
end

"""
    downsample!(obj; <keyword arguments>)

Downsample all channels to `new_sr` sampling frequency in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `new_sr::Int64`: new sampling rate

# Returns

- `Nothing`
"""
function downsample!(obj::NeuroAnalyzer.NEURO; new_sr::Int64)::Nothing
    obj_new = downsample(obj; new_sr = new_sr)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing
end
