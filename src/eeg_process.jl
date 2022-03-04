"""
    eeg_reference_channel(eeg; channel)

Reference the `eeg` to specific channel `channel`.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_reference_channel(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    typeof(channel) <: AbstractRange && (channel = collect(channel))

    s_referenced = signal_reference_channel(eeg.eeg_signals, channel=channel)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_referenced, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:reference] = "channel: $channel"
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_reference_channel(EEG, channel=$channel)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_reference_channel!(eeg; channel)

Reference the `eeg` to specific channel `channel`.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
"""
function eeg_reference_channel!(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    eeg.eeg_signals = signal_reference_channel(eeg.eeg_signals, channel=channel)
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_reference_channel!(EEG, channel=$channel)")

    return
end

"""
    eeg_reference_car(eeg)

Reference the `eeg` to common average reference.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_reference_car(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_referenced = signal_reference_car(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_referenced, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:reference] = "CAR"
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_reference_car(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_reference_car!(eeg)

Reference the `eeg` to common average reference.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_reference_car!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    eeg.eeg_signals = signal_reference_car(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_reference_car!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_derivative(eeg)

Return the derivative of the `eeg` with length same as the signal.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_derivative(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_der = signal_derivative(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_der, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_derivative(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_derivative!(eeg)

Return the derivative of the `eeg` with length same as the signal.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_derivative!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    eeg.eeg_signals = signal_derivative(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_derivative!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_detrend(eeg; type)

Remove linear trend from the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `type::Symbol=:linear`, optional
    - `:linear`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:constant`: the mean of `signal` is subtracted

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_detrend(eeg::NeuroJ.EEG; type::Symbol=:linear)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_det = signal_detrend(eeg.eeg_signals, type=type)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_det, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_detrend(EEG, type=$type)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_detrend!(eeg; type)

Remove linear trend from the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `type::Symbol=:linear`, optional
    - `:linear`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:constant`: the mean of `signal` is subtracted
"""
function eeg_detrend!(eeg::NeuroJ.EEG; type::Symbol=:linear)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    eeg.eeg_signals = signal_detrend(eeg.eeg_signals, type=type)
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_detrend!(EEG, type=$type)")

    return
end

"""
    eeg_taper(eeg; taper)

Taper `eeg` with `taper`.

# Arguments

- `eeg::NeuroJ.EEG`
- `taper::Vector`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_taper(eeg::NeuroJ.EEG; taper::Vector)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_tapered = signal_taper(eeg.eeg_signals, taper=taper)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_tapered, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_taper(EEG, taper=$taper)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_taper!(eeg; taper)

Taper `eeg` with `taper`.

# Arguments

- `eeg::NeuroJ.EEG`
- `taper::Vector`
"""
function eeg_taper!(eeg::NeuroJ.EEG; taper::Vector)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    eeg.eeg_signals = signal_taper(eeg.eeg_signals, taper=taper)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_taper!(EEG, taper=$taper)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_demean(eeg)

Remove mean value (DC offset).

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_demean(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_demeaned = signal_demean(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_demeaned, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_demean(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_demean!(eeg)

Remove mean value (DC offset).

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_demean!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    eeg.eeg_signals = signal_demean(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_demean!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_normalize_zscore(eeg)

Normalize by z-score.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_normalize_zscore(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_normalized = signal_normalize_zscore(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_normalized, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_normalize_zscore(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_normalize_zscore!(eeg)

Normalize by z-score.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_normalize_zscore!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    eeg.eeg_signals = signal_normalize_zscore(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_normalize_zscore!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_normalize_minmax(eeg)

Normalize to 0...1

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_normalize_minmax(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_normalized = signal_normalize_minmax(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_normalized, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_normalize_minmax(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_normalize_minmax!(eeg)

Normalize to 0...1

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_normalize_minmax!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    eeg.eeg_signals = signal_normalize_minmax(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_normalize_minmax!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_upsample(eeg; new_sr)

Upsample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroJ.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_upsample(eeg::NeuroJ.EEG; new_sr::Int64)

    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    s_upsampled, t_upsampled = signal_upsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    # create new dataset
    eeg_time = collect(t_upsampled)
    eeg_new = EEG(deepcopy(eeg.eeg_header), eeg_time, s_upsampled, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_upsampled, 2) * size(s_upsampled, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_upsampled, 2) * size(s_upsampled, 3)) / new_sr
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_upsampled, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_upsampled, 2) / new_sr
    eeg_new.eeg_header[:sampling_rate] = repeat([new_sr], eeg_channel_n(eeg_new))

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_upsample(EEG, new_sr=$new_sr)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_upsample!(eeg; new_sr)

Upsample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroJ.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_upsample!(eeg::NeuroJ.EEG; new_sr::Int64)

    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    eeg.eeg_signals, t_upsampled = signal_upsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    eeg.eeg_time = collect(t_upsampled)
    eeg.eeg_header[:eeg_duration_samples] = eeg_signal_len(eeg) * eeg_epoch_n(eeg)
    eeg.eeg_header[:eeg_duration_seconds] = (eeg_signal_len(eeg) * eeg_epoch_n(eeg)) / new_sr
    eeg.eeg_header[:epoch_duration_seconds] = eeg_signal_len(eeg)
    eeg.eeg_header[:epoch_duration_seconds] = eeg_signal_len(eeg) / new_sr
    eeg.eeg_header[:sampling_rate] = repeat([new_sr], eeg.eeg_header[:channel_n])

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_upsample!(EEG, new_sr=$new_sr)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_tconv(eeg; kernel)

Perform convolution in the time domain.

# Arguments

- `eeg::NeuroJ.EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel used for convolution

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_tconv(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_convoluted = signal_tconv(eeg.eeg_signals, kernel=kernel)

    ## EEG signal can only store Float64
    typeof(kernel) == Vector{ComplexF64} && (s_convoluted = abs.(s_convoluted))

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_convoluted, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_convoluted, 2) * size(s_convoluted, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_convoluted, 2) * size(s_convoluted, 3)) / eeg_sr(eeg_new)
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_convoluted, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_convoluted, 2) / eeg_sr(eeg_new)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_tconv(EEG, kernel=$kernel)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_tconv!(eeg; kernel)

Perform convolution in the time domain.

# Arguments

- `eeg::NeuroJ.EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel used for convolution
"""
function eeg_tconv!(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    # EEG signal can only store Float64
    if typeof(kernel) == Vector{ComplexF64}
        eeg.eeg_signals = abs.(signal_tconv(eeg.eeg_signals, kernel=kernel))
    else
        eeg.eeg_signals = signal_tconv(eeg.eeg_signals, kernel=kernel)
    end

    # create new dataset
    eeg.eeg_header[:eeg_duration_samples] = size(s_convoluted, 2) * size(s_convoluted, 3)
    eeg.eeg_header[:eeg_duration_seconds] = (size(s_convoluted, 2) * size(s_convoluted, 3)) / eeg_sr(eeg)
    eeg.eeg_header[:epoch_duration_samples] = size(s_convoluted, 2)
    eeg.eeg_header[:epoch_duration_seconds] = size(s_convoluted, 2) / eeg_sr(eeg)
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_tconv!(EEG, kernel=$kernel)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_filter(eeg; <keyword arguments>)

Filter `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:mavg`: moving average (with threshold and/or weight window)
    - `:mmed`: moving median (with threshold and/or weight window)
    - `:poly`: polynomial of `order` order
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Int64, Float64, Tuple}=0`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64`=0: filter order
- `rp::Union{Int64, Float64}=-1`: dB ripple in the passband
- `rs::Union{Int64, Float64}=-1`: dB attentuation in the stopband
- `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass)
- `d::Int64=1`: window length for mean average and median average filter
- `t::Union{Int64, Float64}=0`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{Float64}, Nothing}=nothing`: window, required for FIR filter

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_filter(eeg::NeuroJ.EEG; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

    s_filtered = signal_filter(eeg.eeg_signals,
                                    fprototype=fprototype,
                                    ftype=ftype,
                                    cutoff=cutoff,
                                    fs=eeg_sr(eeg),
                                    order=order,
                                    rp=rp,
                                    rs=rs,
                                    dir=dir,
                                    d=d,
                                    t=t,
                                    window=window)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_filtered, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_filter(EEG, fprototype=$fprototype, ftype=$ftype, cutoff=$cutoff, order=$order, rp=$rp, rs=$rs, dir=$dir, d=$d, window=$window)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_filter!(eeg; <keyword arguments>)

Filter `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:mavg`: moving average (with threshold and/or weight window)
    - `:mmed`: moving median (with threshold and/or weight window)
    - `:poly`: polynomial of `order` order
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Int64, Float64, Tuple}=0`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64`=0: filter order
- `rp::Union{Int64, Float64}=-1`: dB ripple in the passband
- `rs::Union{Int64, Float64}=-1`: dB attentuation in the stopband
- `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass)
- `d::Int64=1`: window length for mean average and median average filter
- `t::Union{Int64, Float64}=0`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{Float64}, Nothing}=nothing`: window, required for FIR filter
"""
function eeg_filter!(eeg::NeuroJ.EEG; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

    eeg.eeg_signals = signal_filter(eeg.eeg_signals,
                                    fprototype=fprototype,
                                    ftype=ftype,
                                    cutoff=cutoff,
                                    fs=eeg_sr(eeg),
                                    order=order,
                                    rp=rp,
                                    rs=rs,
                                    dir=dir,
                                    d=d,
                                    t=t,
                                    window=window)

    # create new dataset
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_filter!(EEG, fprototype=$fprototype, ftype=$ftype, cutoff=$cutoff, order=$order, rp=$rp, rs=$rs, dir=$dir, d=$d, window=$window)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_downsample(eeg; new_sr)

Downsamples all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroJ.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_downsample(eeg::NeuroJ.EEG; new_sr::Int64)

    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    s_downsampled, t_downsampled = signal_downsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    # create new dataset
    eeg_time = collect(t_downsampled)
    eeg_new = EEG(deepcopy(eeg.eeg_header), eeg_time, s_downsampled, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_downsampled, 2) * size(s_downsampled, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_downsampled, 2) * size(s_downsampled, 3)) / new_sr
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_downsampled, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_downsampled, 2) / new_sr
    eeg_new.eeg_header[:sampling_rate] = repeat([new_sr], eeg_channel_n(eeg_new))

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_downsample(EEG, new_sr=$new_sr)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_downsample!(eeg; new_sr)

Downsamples all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroJ.EEG`
- `new_sr::Int64`: new sampling rate
"""
function eeg_downsample!(eeg::NeuroJ.EEG; new_sr::Int64)

    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    eeg.eeg_signals, t_downsampled = signal_downsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    # create new dataset
    eeg.eeg_time = collect(t_downsampled)
    eeg.eeg_header[:eeg_duration_samples] = eeg_signal_len(eeg) * eeg_epoch_n(eeg)
    eeg.eeg_header[:eeg_duration_seconds] = (eeg_signal_len(eeg) * eeg_epoch_n(eeg)) / new_sr
    eeg.eeg_header[:epoch_duration_samples] = eeg_signal_len(eeg)
    eeg.eeg_header[:epoch_duration_seconds] = eeg_signal_len(eeg) / new_sr
    eeg.eeg_header[:sampling_rate] = repeat([new_sr], eeg.eeg_header[:channel_n])

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_downsample!(EEG, new_sr=$new_sr)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_pca(eeg; n)

Calculates `n` first PCs for `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `n::Int64`: number of PCs

# Returns

- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}`: PC_VAR(1)..PC_VAR(n) × epoch
"""
function eeg_pca(eeg::NeuroJ.EEG; n::Int64)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    pc, pc_var, m = signal_pca(eeg.eeg_signals, n=n)

    return pc, pc_var, m
end

"""
    eeg_pca!(eeg; n)

Calculates `n` first PCs for `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `n::Int64`: number of PCs
"""
function eeg_pca!(eeg::NeuroJ.EEG; n::Int64)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :pca in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:pca)
    :pca_var in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:pca_var)
    :pca_m in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:pca_m)
    pc, pc_var, m = eeg_pca(eeg, n=n)
    push!(eeg.eeg_components, pc)
    push!(eeg.eeg_components, pc_var)
    push!(eeg.eeg_components, m)
    push!(eeg.eeg_header[:components], :pca)
    push!(eeg.eeg_header[:components], :pca_var)
    push!(eeg.eeg_header[:components], :pca_m)
    push!(eeg.eeg_header[:history], "eeg_pca!(EEG, n=$n)")

    return
end

"""
    eeg_fconv(eeg, kernel)

Performs convolution in the time domain.

# Arguments

- `eeg::NeuroJ.EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel for convolution

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_fconv(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_convoluted = signal_fconv(eeg.eeg_signals, kernel=kernel)

    ## EEG signal can only store Float64
    s_convoluted = abs.(s_convoluted)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_convoluted, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_convoluted, 2) * size(s_convoluted, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_convoluted, 2) * size(s_convoluted, 3)) / eeg_sr(eeg_new)
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_convoluted, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_convoluted, 2) / eeg_sr(eeg_new)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_fconv(EEG, kernel=$kernel)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_fconv!(eeg, kernel)

Performs convolution in the time domain.

# Arguments

- `eeg::NeuroJ.EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel for convolution
"""
function eeg_fconv!(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    # EEG signal can only store Float64
    eeg.eeg_signals = abs.(signal_fconv(eeg.eeg_signals, kernel=kernel))

    eeg.eeg_header[:eeg_duration_samples] = size(s_convoluted, 2) * size(s_convoluted, 3)
    eeg.eeg_header[:eeg_duration_seconds] = (size(s_convoluted, 2) * size(s_convoluted, 3)) / eeg_sr(eeg)
    eeg.eeg_header[:epoch_duration_samples] = size(s_convoluted, 2)
    eeg.eeg_header[:epoch_duration_seconds] = size(s_convoluted, 2) / eeg_sr(eeg)
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_fconv!(EEG, kernel=$kernel)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_ica(eeg; <keyword arguments>)

Calculates `n` first ICs for `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `n::Int64`: number of ICs
- `tol::Float64=1.0e-6`: tolerance for ICA
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor: :tanh, :gaus

# Returns

- `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
"""
function eeg_ica(eeg::NeuroJ.EEG; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    ic, ic_mw = signal_ica(eeg.eeg_signals, n=n, tol=tol, iter=iter, f=f)

    return ic, ic_mw
end

"""
    eeg_ica!(eeg; <keyword arguments>)

Calculates `n` first ICs for `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `n::Int64`: number of ICs
- `tol::Float64=1.0e-6`: tolerance for ICA
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor: :tanh, :gaus
"""
function eeg_ica!(eeg::NeuroJ.EEG; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :ica in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:ica)
    :ica_mw in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:ica_mw)
    ic, ic_mw = signal_ica(eeg.eeg_signals, n=n, tol=tol, iter=iter, f=f)
    push!(eeg.eeg_components, ic)
    push!(eeg.eeg_components, ic_mw)
    push!(eeg.eeg_header[:components], :ica)
    push!(eeg.eeg_header[:components], :ica_mw)
    push!(eeg.eeg_header[:history], "eeg_ica!(EEG, n=$n, tol=$tol, iter=$iter, f=$f))")

    return
end

"""
    eeg_average(eeg)

Returns the average signal of all `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_average(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    # create new dataset
    eeg_new = deepcopy(eeg)
    eeg_keep_channel!(eeg_new, channel=1)
    eeg_new.eeg_signals = signal_average(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_average(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_average!(eeg)

Returns the average signal of all `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_average!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_avg = signal_average(eeg.eeg_signals)
    eeg_delete_channel!(eeg, channel=2:eeg_channel_n(eeg))
    eeg.eeg_signals = s_avg

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_average!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_ica_reconstruct(eeg; ica)

Reconstructs `eeg` signals using removal of `ica` ICA components.

# Arguments

- `eeg::NeuroJ.EEG`
- `ica::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_ica_reconstruct(eeg::NeuroJ.EEG; ica::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :ica in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain ICA. Perform eeg_ica!(EEG) first."))
    :ica_mw in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain ICA. Perform eeg_ica!(EEG) first."))

    eeg_new = deepcopy(eeg)
    ica_a_idx = findfirst(isequal(:ica), eeg.eeg_header[:components])
    ica_mw_idx = findfirst(isequal(:ica_mw), eeg.eeg_header[:components])
    eeg_new.eeg_signals = signal_ica_reconstruct(eeg_new.eeg_signals, ic_activations=eeg_new.eeg_components[ica_a_idx], ic_mw=eeg_new.eeg_components[ica_mw_idx], ic_v=ica)
    
    push!(eeg_new.eeg_header[:history], "eeg_ica_reconstruct(EEG, ica=$ica")

    return eeg_new
end

"""
    eeg_ica_reconstruct!(eeg; ica)

Reconstructs `eeg` signals using removal of `ica` ICA components.

# Arguments

- `eeg::NeuroJ.EEG`
- `ica::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove
"""
function eeg_ica_reconstruct!(eeg::NeuroJ.EEG; ica::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :ica in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain ICA. Perform eeg_ica!(EEG) first."))
    :ica_mw in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain ICA. Perform eeg_ica!(EEG) first."))

    ica_a_idx = findfirst(isequal(:ica), eeg.eeg_header[:components])
    ica_mw_idx = findfirst(isequal(:ica_mw), eeg.eeg_header[:components])
    eeg.eeg_signals = signal_ica_reconstruct(eeg.eeg_signals, ic_activations=eeg.eeg_components[ica_a_idx], ic_mw=eeg.eeg_components[ica_mw_idx], ic_v=ica)

    push!(eeg.eeg_header[:history], "eeg_ica_reconstruct!(EEG, ica=$ica")

    return
end

"""
    eeg_resample(eeg; new_sr)

Resample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroJ.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_resample(eeg::NeuroJ.EEG; new_sr::Int64)

    new_sr > eeg_sr(eeg) && (eeg_new = eeg_upsample(eeg, new_sr=new_sr))
    new_sr < eeg_sr(eeg) && (eeg_new = eeg_downsample(eeg, new_sr=new_sr))
    new_sr = eeg_sr(eeg) && (eeg_new = eeg)

    return eeg_new
end

"""
    eeg_resample!(eeg; new_sr)

Resample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroJ.EEG`
- `new_sr::Int64`: new sampling rate
"""
function eeg_resample!(eeg::NeuroJ.EEG; new_sr::Int64)

    new_sr > eeg_sr(eeg) && eeg_upsample!(eeg, new_sr=new_sr)
    new_sr < eeg_sr(eeg) && eeg_downsample!(eeg, new_sr=new_sr)

    return
end