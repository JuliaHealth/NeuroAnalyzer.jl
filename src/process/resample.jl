export resample
export resample!
export upsample
export upsample!
export downsample
export downsample!

"""
    resample(s; t, new_sr)

Resample to `new_sr` sampling frequency.

# Arguments

- `s::AbstractVector`
- `old_sr::Int64`: old sampling rate
- `new_sr::Int64`: new sampling rate

# Returns

- `s_new::Vector{Float64}`
"""
function resample(s::AbstractVector; old_sr::Int64, new_sr::Int64)

    @assert old_sr >= 1 "old_sr must be ≥ 1."
    @assert new_sr >= 1 "new_sr must be ≥ 1."

    new_sr == old_sr && return(s)

    # resample
    sr_ratio = new_sr / old_sr
    s_new = DSP.resample(s, sr_ratio)

    return s_new

end

"""
    resample(s; old_sr::Int64, new_sr::Int64)

Resamples all channels and time vector `t` to `new_sr` sampling frequency.

# Arguments

- `s::AbstractArray`
- `old_sr::Int64`: old sampling rate
- `new_sr::Int64`: new sampling rate

# Returns

- `s_new::Array{Float64, 3}`
"""
function resample(s::AbstractArray; old_sr::Int64, new_sr::Int64)

    @assert new_sr >= 1 "new_sr must be ≥ 1."

    ch_n, _, ep_n = size(s)

    s_new = NeuroAnalyzer.resample(s[1, :, 1], old_sr=old_sr, new_sr=new_sr)
    s_new = zeros(ch_n, length(s_new), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views NeuroAnalyzer.resample(s[ch_idx, :, ep_idx], old_sr=old_sr, new_sr=new_sr)
        end
    end

    return s_new

end

"""
    resample(obj; new_sr)

Resample (up- or down-sample).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `old_sr::Int64`: old sampling rate - `new_sr::Int64`: new sampling rate

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function resample(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    @assert new_sr >= 1 "new_sr must be ≥ 1."

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

    obj_new = resample(obj, new_sr=new_sr)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing

end

"""
    upsample(obj; new_sr)

Upsample.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `new_sr::Int64`: new sampling rate

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function upsample(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    new_sr / sr(obj) != new_sr ÷ sr(obj) && _warn("New sampling rate should be easily captured by integer fractions, e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.")

    obj_new = deepcopy(obj)

    s_new = NeuroAnalyzer.resample(obj.data, old_sr=sr(obj), new_sr=new_sr)

    obj_new.data = s_new

    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    obj_new.header.recording[:sampling_rate] = new_sr
    reset_components!(obj_new)
    push!(obj_new.history, "upsample(OBJ, new_sr=$new_sr)")

    return obj_new

end

"""
    upsample!(obj; new_sr)

Upsample.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `new_sr::Int64`: new sampling rate
"""
function upsample!(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    obj_new = upsample(obj, new_sr=new_sr)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing

end

"""
    downsample(obj; new_sr)

Downsample.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `new_sr::Int64`: new sampling rate

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function downsample(obj::NeuroAnalyzer.NEURO; new_sr::Int64)

    new_sr < sr(obj) && _warn("To prevent aliasing due to down-sampling, a low-pass filter should be applied before removing data points. The filter cutoff should be the Nyquist frequency of the new down-sampled rate, ($(new_sr / 2) Hz), not the original Nyquist frequency ($(sr(obj) / 2) Hz).")

    new_sr / sr(obj) != new_sr ÷ sr(obj) && _warn("New sampling rate should be easily captured by integer fractions e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.")

    obj_new = deepcopy(obj)

    s_new = NeuroAnalyzer.resample(obj.data, old_sr=sr(obj), new_sr=new_sr)

    obj_new.data = s_new

    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    obj_new.header.recording[:sampling_rate] = new_sr
    reset_components!(obj_new)
    push!(obj_new.history, "downsample(OBJ, new_sr=$new_sr)")

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

    obj_new = downsample(obj, new_sr=new_sr)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing

end
