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
    obj_new = deepcopy(obj)
    resample!(obj_new; new_sr = new_sr)

    return obj_new
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
    # validate
    new_sr >= 1 || throw(ArgumentError("new_sr must be ≥ 1."))

    if new_sr > sr(obj)
        upsample!(obj; new_sr = new_sr)
    elseif new_sr < sr(obj)
        downsample!(obj; new_sr = new_sr)
    end

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
    obj_new = deepcopy(obj)
    upsample!(obj_new; new_sr = new_sr)

    return obj_new
end

# resampled copy of obj: data, sampling rate and time axes (markers are stored in seconds
# and do not change)
function _resample_obj!(
    obj::NeuroAnalyzer.NEURO,
    new_sr::Int64,
    fname::String,
)::Nothing
    r = new_sr // sr(obj)
    max(numerator(r), denominator(r)) > 1000 && _warn(
        "The resampling ratio $(numerator(r))/$(denominator(r)) is complex and slow; a sampling rate with a simpler ratio to $(sr(obj)) Hz is faster.",
    )

    obj.data = NeuroAnalyzer._resample(obj.data; old_sr = sr(obj), new_sr = new_sr)
    # anti-aliasing: the signal contains no frequencies above the new Nyquist frequency
    rec = obj.header.recording
    _rec_filters!(rec)
    rec[:lowpass] = min.(rec[:lowpass], new_sr / 2)
    rec[:highpass] = min.(rec[:highpass], new_sr / 2)
    # notches at or above the new Nyquist frequency are outside the signal band
    rec[:notch] = Vector{Float64}[Base.filter(<(new_sr / 2), n) for n in rec[:notch]]
    # the new sampling rate must be set before the time axes are rebuilt
    rec[:sampling_rate] = new_sr
    obj.time_pts, obj.epoch_time = _get_t(obj)
    push!(obj.history, "$fname(obj, new_sr=$new_sr)")

    return nothing
end

# out-of-place version of _resample_obj!()
function _resample_obj(
    obj::NeuroAnalyzer.NEURO,
    new_sr::Int64,
    fname::String,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    _resample_obj!(obj_new, new_sr, fname)

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
    # validate
    new_sr > sr(obj) ||
        throw(ArgumentError("new_sr must be > $(sr(obj)) Hz (the current sampling rate)."))

    _resample_obj!(obj, new_sr, "upsample")

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
    obj_new = deepcopy(obj)
    downsample!(obj_new; new_sr = new_sr)

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
    # validate
    1 <= new_sr < sr(obj) ||
        throw(
            ArgumentError(
                "new_sr must be in [1, $(sr(obj) - 1)] Hz (below the current sampling rate).",
            ),
        )

    # (the resampling filter is an anti-aliasing low-pass filter, cutoff: the new Nyquist frequency)
    _resample_obj!(obj, new_sr, "downsample")

    return nothing
end
