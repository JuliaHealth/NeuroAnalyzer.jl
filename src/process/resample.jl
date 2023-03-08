export resample
export resample!
export upsample
export upsample!
export downsample
export downsample!

"""
    resample(obj; new_sr)

Resample (up- or down-sample).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `new_sr::Int64`: new sampling rate

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function resample(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    new_sr < 1 && throw(ArgumentError("new_sr must be ≥ 1."))
    new_sr > sr(obj) && return upsample(obj, new_sr=new_sr)
    new_sr < sr(obj) && return downsample(obj, new_sr=new_sr)
    new_sr == sr(obj) && return obj

end

"""
    resample!(obj; new_sr)

Resample (up- or down-sample).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `new_sr::Int64`: new sampling rate
"""
function resample!(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    obj_tmp = resample(obj, new_sr=new_sr)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end

"""
    upsample(obj; new_sr)

Upsample.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `new_sr::Int64`: new sampling rate

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function upsample(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    new_sr / sr(obj) != new_sr ÷ sr(obj) && _info("New sampling rate should be easily captured by integer fractions, e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.")
    
    obj_new = deepcopy(obj)

    t = obj.time_pts[1]:(1 / obj.header.recording[:sampling_rate]):obj.time_pts[end]
    s_upsampled, t_upsampled = s_resample(obj_new.data, t=t, new_sr=new_sr)

    t_upsampled = collect(t_upsampled)
    obj_new.data = s_upsampled
    obj_new.time_pts = t_upsampled
    obj_new.epoch_time = linspace(obj_new.epoch_time[1], obj_new.epoch_time[end], size(s_upsampled, 2))
    obj_new.header.recording[:duration_samples] = size(s_upsampled, 2) * size(s_upsampled, 3)
    obj_new.header.recording[:duration_seconds] = (size(s_upsampled, 2) * size(s_upsampled, 3)) / new_sr
    obj_new.header.recording[:epoch_duration_samples] = size(s_upsampled, 2)
    obj_new.header.recording[:epoch_duration_seconds] = size(s_upsampled, 2) / new_sr
    obj_new.header.recording[:sampling_rate] = new_sr
    reset_components!(obj_new)
    push!(obj_new.header[:history], "upsample(OBJ, new_sr=$new_sr)")

    return obj_new
end

"""
    upsample!(obj; new_sr)

Upsample.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `new_sr::Int64`: new sampling rate

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function upsample!(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    obj_tmp = upsample(obj, new_sr=new_sr)
    obj.data = obj_tmp.data
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end

"""
    downsample(obj; new_sr)

Downsample.
.
# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `new_sr::Int64`: new sampling rate

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function downsample(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    new_sr < sr(obj) && _info("To prevent aliasing due to down-sampling, a low-pass filter should be applied before removing data points. The filter cutoff should be the Nyquist frequency of the new down-sampled rate, ($(new_sr / 2) Hz), not the original Nyquist frequency ($(sr(obj) / 2) Hz).")

    new_sr / sr(obj) != new_sr ÷ sr(obj) && _info("New sampling rate should be easily captured by integer fractions e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.")

    obj_new = deepcopy(obj)

    t = obj.time_pts[1]:(1 / obj.header.recording[:sampling_rate]):obj.time_pts[end]
    s_downsampled, t_downsampled = s_resample(obj_new.data, t=t, new_sr=new_sr)

    t_downsampled = collect(t_downsampled)
    obj_new.time_pts = t_downsampled
    obj_new.data = s_downsampled
    obj_new.epoch_time = linspace(obj_new.epoch_time[1], obj_new.epoch_time[end], size(s_downsampled, 2))
    obj_new.header.recording[:duration_samples] = size(s_downsampled, 2) * size(s_downsampled, 3)
    obj_new.header.recording[:duration_seconds] = (size(s_downsampled, 2) * size(s_downsampled, 3)) / new_sr
    obj_new.header[:epoch_duration_samples] = size(s_downsampled, 2)
    obj_new.header[:epoch_duration_seconds] = size(s_downsampled, 2) / new_sr
    obj_new.header[:sampling_rate] = new_sr
    reset_components!(obj_new)
    push!(obj_new.header[:history], "downsample(OBJ, new_sr=$new_sr)")

    return obj_new
end

"""
    downsample!(obj; new_sr)

Downsample.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `new_sr::Int64`: new sampling rate
"""
function downsample!(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    obj_tmp = downsample(obj, new_sr=new_sr)
    obj.data = obj_tmp.data
    obj.time_pts = obj_tmp.time_pts
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
