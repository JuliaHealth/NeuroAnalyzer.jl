export resample
export resample!
export upsample
export upsample!
export downsample
export downsample!

"""
    _resample(s; <keyword arguments>)

Resample to `new_sr` sampling frequency.

# Arguments

- `s::AbstractVector`: signal vector
- `old_sr::Int64`: old sampling rate
- `new_sr::Int64`: new sampling rate

# Returns

- `Vector{Float64}`
"""
function _resample(s::AbstractVector; old_sr::Int64, new_sr::Int64)::Vector{Float64}
    # validate
    old_sr >= 1 || throw(ArgumentError("old_sr must be ≥ 1."))
    new_sr >= 1 || throw(ArgumentError("new_sr must be ≥ 1."))

    new_sr == old_sr && return Vector{Float64}(s)

    # exact rational ratio: polyphase resampling with an anti-aliasing low-pass filter
    # (a floating-point ratio would select DSP's approximate arbitrary-rate filter)
    return DSP.resample(Vector{Float64}(s), new_sr // old_sr)
end

"""
    _resample(s; <keyword arguments>)

Resamples all channels and time vector `t` to `new_sr` sampling frequency for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `old_sr::Int64`: old sampling rate
- `new_sr::Int64`: new sampling rate

# Returns

- `Array{Float64, 3}`
"""
function _resample(s::AbstractArray; old_sr::Int64, new_sr::Int64)::Array{Float64, 3}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # validate
    new_sr >= 1 || throw(ArgumentError("new_sr must be ≥ 1."))

    # number of channels and epochs
    ch_n, _, ep_n = size(s)

    # dry run
    s_tmp = NeuroAnalyzer._resample(s[1, :, 1]; old_sr = old_sr, new_sr = new_sr)

    # pre-allocate output
    s_new = zeros(ch_n, length(s_tmp), ep_n)

    # calculate over channels and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = NeuroAnalyzer._resample(
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
        # always a copy, never the input object
        return deepcopy(obj)
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
    obj_tmp = resample(obj; new_sr = new_sr)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj_tmp = nothing

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
    new_sr > sr(obj) ||
        throw(ArgumentError("new_sr must be > $(sr(obj)) Hz (the current sampling rate)."))

    return _resample_obj(obj, new_sr, "upsample")
end

# resampled copy of obj: data, sampling rate and time axes (markers are stored in seconds
# and do not change)
function _resample_obj(
    obj::NeuroAnalyzer.NEURO,
    new_sr::Int64,
    fname::String,
)::NeuroAnalyzer.NEURO
    r = new_sr // sr(obj)
    max(numerator(r), denominator(r)) > 1000 && _warn(
        "The resampling ratio $(numerator(r))/$(denominator(r)) is complex and slow; a sampling rate with a simpler ratio to $(sr(obj)) Hz is faster.",
    )

    # create new dataset
    obj_tmp = deepcopy(obj)
    obj_tmp.data = NeuroAnalyzer._resample(obj.data; old_sr = sr(obj), new_sr = new_sr)
    # the new sampling rate must be set before the time axes are rebuilt
    obj_tmp.header.recording[:sampling_rate] = new_sr
    obj_tmp.time_pts, obj_tmp.epoch_time = _get_t(obj_tmp)
    push!(obj_tmp.history, "$fname(obj, new_sr=$new_sr)")

    return obj_tmp
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
    obj_tmp = upsample(obj; new_sr = new_sr)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj_tmp = nothing

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
    1 <= new_sr < sr(obj) ||
        throw(ArgumentError("new_sr must be in [1, $(sr(obj) - 1)] Hz (below the current sampling rate)."))

    # (the resampling filter is an anti-aliasing low-pass filter, cutoff: the new Nyquist frequency)
    return _resample_obj(obj, new_sr, "downsample")
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
    obj_tmp = downsample(obj; new_sr = new_sr)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj_tmp = nothing

    return nothing
end
