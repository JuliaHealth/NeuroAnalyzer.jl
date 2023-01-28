"""
    eeg_reference_ch(eeg; channel, med)

Reference to selected channel(s). Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_ch(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}, med::Bool=false)

    # keep EEG channels
    _check_channels(eeg, channel)
    eeg_new = deepcopy(eeg)
    channels = eeg_signal_channels(eeg)
    signal = @view eeg_new.eeg_signals[channels, :, :]

    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if length(channel) == 1
                reference_channel = @views vec(signal[channel, :, epoch_idx])
                if channel_idx != channel
                    @views signal[channel_idx, :, epoch_idx] .-= reference_channel
                end
            else
                if med == false
                    reference_channel = @views vec(mean(signal[channel, :, epoch_idx], dims=1))
                else
                    reference_channel = @views vec(median(signal[channel, :, epoch_idx], dims=1))
                end
                @views signal[channel_idx, :, epoch_idx] .-= reference_channel
            end
        end
    end

    eeg_new.eeg_signals[channels, :, :] = signal
    eeg_new.eeg_header[:reference] = "channel: $channel"
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_reference_ch(EEG, channel=$channel, med=$med")

    return eeg_new
end

"""
    eeg_reference_ch!(eeg; channel, med)

Reference to selected channel(s). Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_ch!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}, med::Bool=false)

    eeg_tmp = eeg_reference_ch(eeg, channel=channel, med=med)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_reference_car(eeg; exclude_fpo, exclude_current, med)

Reference to common average reference. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `exclude_fpo::Bool=false`: exclude Fp1, Fp2, O1, O2 from CAR calculation
- `exclude_current::Bool=true`: exclude current channel from CAR calculation
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_car(eeg::NeuroAnalyzer.EEG; exclude_fpo::Bool=false, exclude_current::Bool=true, med::Bool=false)

    # keep EEG channels
    eeg_new = deepcopy(eeg)
    channels = eeg_signal_channels(eeg)
    signal = @view eeg_new.eeg_signals[channels, :, :]

    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            channels2exclude = Vector{Int64}()
            if exclude_fpo == true
                l = lowercase.(eeg_labels(eeg_new))
                "fp1" in l && push!(channels2exclude, findfirst(isequal("fp1"), l))
                "fp2" in l && push!(channels2exclude, findfirst(isequal("fp2"), l))
                "o1" in l && push!(channels2exclude, findfirst(isequal("o1"), l))
                "o2" in l && push!(channels2exclude, findfirst(isequal("o2"), l))
            end
            exclude_current == true && push!(channels2exclude, channel_idx)
            reference_channels = @view signal[setdiff(1:channel_n, unique(channels2exclude)), :, epoch_idx]
            if med == false
                reference_channel = vec(mean(reference_channels, dims=1))
            else
                reference_channel = vec(median(reference_channels, dims=1))
            end
            @views signal[channel_idx, :, epoch_idx] .-= reference_channel
        end
    end

    eeg_new.eeg_signals[channels, :, :] = signal
    eeg_new.eeg_header[:reference] = "CAR"
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_reference_car(EEG, exclude_fpo=$exclude_fpo, exclude_current=$exclude_current, med=$med))")

    return eeg_new
end

"""
    eeg_reference_car!(eeg; exclude_fpo, exclude_current, med)

Reference to common average reference. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `exclude_fpo::Bool=false`: exclude Fp1, Fp2, O1, O2 from CAR mean calculation
- `exclude_current::Bool=true`: exclude current channel from CAR mean calculation
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_car!(eeg::NeuroAnalyzer.EEG; exclude_fpo::Bool=false, exclude_current::Bool=true, med::Bool=false)

    eeg_tmp = eeg_reference_car(eeg, exclude_fpo=exclude_fpo, exclude_current=exclude_current, med=med)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_derivative(eeg; channel)

Return the derivative of EEG channel(s) with length same as the signal.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_derivative(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    epoch_n = eeg_epoch_n(eeg)
    
    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in eachindex(channel)
            @views eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = s_derivative(eeg.eeg_signals[channel[channel_idx], :, epoch_idx])
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_derivative(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_derivative!(eeg; channel)

Return the derivative of EEG channel(s) with length same as the signal.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_derivative!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_derivative(eeg, channel=channel)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_detrend(eeg; channel, type, offset, order, f)

Perform piecewise detrending of EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `type::Symbol=:linear`: detrending method
    - `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:linear`: linear trend is subtracted from `signal`
    - `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` is subtracted
    - `:loess`: fit and subtract LOESS approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for `:constant` detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_detrend(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)

    _check_var(type, [:ls, :linear, :constant, :poly, :loess, :hp], "type")

    epoch_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in eachindex(channel)
            @views eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = s_detrend(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], type=type, offset=offset, order=order, f=f, fs=fs)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_detrend(EEG, channel=$channel, type=$type, offset=$offset, order=$order, f=$f)")

    return eeg_new
end

"""
    eeg_detrend!(eeg; channel, type, offset, order, span)

Perform piecewise detrending of EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `type::Symbol=:linear`: detrending method
    - `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:linear`: linear trend is subtracted from `signal`
    - `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` order is subtracted
    - `:loess`: fit and subtract LOESS approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`
"""
function eeg_detrend!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)

    eeg_tmp = eeg_detrend(eeg, channel=channel, type=type, offset=offset, order=order, f=f)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_taper(eeg; channel, taper)

Taper EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `taper::Union{Vector{Real, Vector{ComplexF64}}`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_taper(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    _check_channels(eeg, channel)

    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in eachindex(channel)
            @views eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = s_taper(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], taper=taper)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_taper(EEG, taper=$taper, channel=$channel)")

    return eeg_new
end

"""
    eeg_taper!(eeg; channel, taper)

Taper EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `taper::Union{Vector{<:Real}, Vector{ComplexF64}}`
"""
function eeg_taper!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    eeg_tmp = eeg_taper(eeg, channel=channel, taper=taper)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_demean(eeg; channel)

Remove mean value (DC offset).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_demean(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)

    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in eachindex(channel)
            @views eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = s_demean(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx])
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_demean(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_demean!(eeg; channel)

Remove mean value (DC offset).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_demean!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_demean(eeg, channel=channel)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_normalize(eeg; channel, method)

Normalize EEG channel(s)

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `method::Symbol`: method for normalization, see `s_normalize()` for details

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_normalize(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), method::Symbol)

    _check_channels(eeg, channel)
    channel_n = length(channel)
    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            @views eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = s_normalize(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], method=method)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_normalize(EEG, channel=$channel, method=$method)")

    return eeg_new
end

"""
    eeg_normalize!(eeg; channel, method)

Normalize EEG channel(s)

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `method::Symbol`: method for normalization, see `s_normalize()` for details
"""
function eeg_normalize!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), method::Symbol)

    eeg_tmp = eeg_normalize(eeg, channel=channel, method=method)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_add_noise(eeg; channel)

Add random noise to EEG channel(s)

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_add_noise(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)

    channel_n = length(channel)
    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in eachindex(channel)
            @views eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = s_add_noise(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx])
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_add_noise(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_add_noise!(eeg; channel)

Add random noise to EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_add_noise!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_add_noise(eeg, channel=channel)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_filter(eeg; <keyword arguments>)

Apply filtering to EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`: second-order IIR notch filter
    - `:remez`: Remez FIR filter
    - `:mavg`: moving average (with threshold)
    - `:mmed`: moving median (with threshold)
    - `:poly`: polynomial of `order`
    - `:conv`: convolution
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
- `order::Int64=8`: filter order, number of taps for :remez filter, k-value for :mavg and :mmed (window length = 2 × k + 1)
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Real=-1`: bandwidth for `:iirnotch` and :remez filters
- `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, `:twopass`), for causal filter use `:onepass`
- `t::Real`: threshold for `:mavg` and `:mmed` filters; threshold = threshold * std(signal) + mean(signal) for `:mavg` or threshold = threshold * std(signal) + median(signal) for `:mmed` filter
- `window::Union{AbstractVector, Nothing} - window, required for FIR filter, weighting window for `:mavg` and `:mmed` 

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_filter(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple{Real, Real}}=0, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, dir::Symbol=:twopass, d::Int64=1, t::Real=0, window::Union{Vector{<:Real}, Nothing}=nothing, preview::Bool=false)

    _check_channels(eeg, channel)

    epoch_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)

    if preview == true
        verbose == true && @info "When `preview=true`, signal is not being filtered."
        fprototype === :iirnotch && (ftype = :bs)    
        p = plot_filter_response(fs=fs, fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, window=window)
        Plots.plot(p)
        return p
    end

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in eachindex(channel)
            eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = @views s_filter(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], fprototype=fprototype, ftype=ftype, cutoff=cutoff, fs=fs, order=order, rp=rp, rs=rs, bw=bw, dir=dir, t=t, window=window)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_filter(EEG, channel=$channel, fprototype=$fprototype, ftype=$ftype, cutoff=$cutoff, order=$order, rp=$rp, rs=$rs, dir=$dir, t=$t, window=$window)")

    return eeg_new
end

"""
    eeg_filter!(eeg; <keyword arguments>)

Apply filtering to EEG channel(s).

# Arguments

- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`: second-order IIR notch filter
    - `:remez`: Remez FIR filter
    - `:mavg`: moving average (with threshold)
    - `:mmed`: moving median (with threshold)
    - `:poly`: polynomial of `order`
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
- `order::Int64=8`: filter order, number of taps for `:remez` filter, k-value for `:mavg` and `:mmed` (window length = 2 × k + 1)
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Real=-1`: bandwidth for `:iirnotch` and `:remez filters
- `dir:Symbol=:twopass`: filter direction (`:onepass`, :onepass_reverse, `:twopass`), for causal filter use `:onepass`
- `t::Real`: threshold for `:mavg` and `:mmed` filters; threshold = threshold * std(signal) + mean(signal) for `:mavg` or threshold = threshold * std(signal) + median(signal) for `:mmed` filter
- `window::Union{Vector{<:Real}, Nothing} - window, required for `:fir` filter
- `preview::Bool=false`: plot filter response
"""
function eeg_filter!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple{Real, Real}}=0, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, dir::Symbol=:twopass, t::Real=0, window::Union{Vector{<:Real}, Nothing}=nothing, preview::Bool=false)

    if preview == true
        verbose == true && @info "When `preview=true`, signal is not being filtered."
        fprototype === :iirnotch && (ftype = :bs)
        p = plot_filter_response(fs=eeg_sr(eeg), fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, window=window)
        Plots.plot(p)
        return p
    end

    eeg_tmp = eeg_filter(eeg,
                         channel=channel,
                         fprototype=fprototype,
                         ftype=ftype,
                         cutoff=cutoff,
                         order=order,
                         rp=rp,
                         rs=rs,
                         bw=bw,
                         dir=dir,
                         t=t,
                         window=window)

    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_pca(eeg; channel, n)

Perform principal component analysis (PCA).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `n::Int64`: number of PCs to calculate

# Returns

Named tuple containing:
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
- `pc_m::Vector{Float64}`: PC mean
- `pca::MultivariateStats.PCA{Float64}`: PC model
"""
function eeg_pca(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), n::Int64)

    _check_channels(eeg, channel)

    pc, pc_var, pc_m, pca = @views s_pca(eeg.eeg_signals[channel, :, :], n=n)

    return (pc=pc, pc_var=pc_var, pc_m=pc_m, pca=pca)
end

"""
    eeg_pca_reconstruct(eeg; channel)

Reconstruct EEG signals using embedded PCA components (`:pc`) and model (`:pca`).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_pca_reconstruct(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    :pc in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain :pc component. Perform eeg_pca(EEG) first."))
    :pca in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain :pca component. Perform eeg_pca(EEG) first."))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    pc_idx = findfirst(isequal(:pc), eeg.eeg_header[:components])
    pc_m_idx = findfirst(isequal(:pca), eeg.eeg_header[:components])
    eeg_new.eeg_signals[channel, :, :] = @views s_pca_reconstruct(eeg.eeg_signals[channel, :, :], pc=eeg_new.eeg_components[pc_idx], pca=eeg_new.eeg_components[pc_m_idx])
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_pca_reconstruct(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_pca_reconstruct!(eeg; channel)

Reconstruct EEG signals using embedded PCA components (`:pc`) and model (`:pca`).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_pca_reconstruct!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_pca_reconstruct(eeg, channel=channel)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_pca_reconstruct(eeg, pc, pca; channel)

Reconstruct EEG signals using external PCA components (`pc` and `pca`).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pca::MultivariateStats.PCA{Float64}`: PC model
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_pca_reconstruct(eeg::NeuroAnalyzer.EEG, pc::Array{Float64, 3}, pca::MultivariateStats.PCA{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channel, :, :] = @views s_pca_reconstruct(eeg_new.eeg_signals[channel, :, :], pc=pc, pca=pca)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_pca_reconstruct(EEG, channel=$channel, pc=$pc, pca=$pca)")

    return eeg_new
end

"""
    eeg_pca_reconstruct!(eeg, pc, pca; channel)

Reconstruct EEG signals using external PCA components (`pc` and `pca`).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pca::MultivariateStats.PCA{Float64}`: PC model
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_pca_reconstruct!(eeg::NeuroAnalyzer.EEG, pc::Array{Float64, 3}, pca::MultivariateStats.PCA{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_pca_reconstruct(eeg, pc, pca, channel=channel)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_ica(eeg; <keyword arguments>)

Perform independent component analysis (ICA).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `n::Int64`: number of ICs
- `tol::Float64=1.0e-6`: tolerance for ICA
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor: :tanh, :gaus

# Returns

Named tuple containing:
- `ica::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ica_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
"""
function eeg_ica(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    _check_channels(eeg, channel)

    ica, ica_mw = @views s_ica(eeg.eeg_signals[channel, :, :], n=n, tol=tol, iter=iter, f=f)

    return (ica=ica, ica_mw=ica_mw)
end

"""
    eeg_ica_reconstruct(eeg; channel, ic)

Reconstruct EEG signals using embedded ICA components (`:ica` and `:ica_mw`).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `ic::Union{Int64, Vector{Int64}, AbstractRange}`: list of ICs to remove

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_ica_reconstruct(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), ic::Union{Int64, Vector{Int64}, AbstractRange})

    :ica in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain :ica component. Perform eeg_ica(EEG) first."))
    :ica_mw in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain :ica_mw component. Perform eeg_ica(EEG) first."))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    ica_a_idx = findfirst(isequal(:ica), eeg.eeg_header[:components])
    ica_mw_idx = findfirst(isequal(:ica_mw), eeg.eeg_header[:components])
    eeg_new.eeg_signals[channel, :, :] = @views s_ica_reconstruct(eeg_new.eeg_signals[channel, :, :], ica=eeg_new.eeg_components[ica_a_idx], ica_mw=eeg_new.eeg_components[ica_mw_idx], ic=ic)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_ica_reconstruct(EEG, channel=$channel, ic=$ic)")

    return eeg_new
end

"""
    eeg_ica_reconstruct!(eeg; channel, ic)

Reconstruct EEG signals using embedded ICA components (`:ica` and `:ica_mw`).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove
"""
function eeg_ica_reconstruct!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), ic::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_tmp = eeg_ica_reconstruct(eeg, channel=channel, ic=ic)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_ica_reconstruct(eeg, ica, ica_mw; channel, ic)

Reconstruct EEG signals using external ICA components (`ica` and `ica_mw`).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `ica::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ica_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_ica_reconstruct(eeg::NeuroAnalyzer.EEG, ica::Array{Float64, 3}, ica_mw::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), ic::Union{Int64, Vector{Int64}, AbstractRange})

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channel, :, :] = @views s_ica_reconstruct(eeg_new.eeg_signals[channel, :, :], ica=ica, ica_mw=ica_mw, ic=ic)
    
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_ica_reconstruct(EEG, ic=$ica, ica_mw=$ica_mw, channel=$channel, ic=$ic)")

    return eeg_new
end

"""
    eeg_ica_reconstruct!(eeg, ica, ica_mw; channel, ic)

Reconstruct EEG signals using external ICA components (`ica` and `ica_mw`).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `ica::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ica_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove
"""
function eeg_ica_reconstruct!(eeg::NeuroAnalyzer.EEG, ica::Array{Float64, 3}, ica_mw::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), ic::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_tmp = eeg_ica_reconstruct(eeg, ica, ica_mw, channel=channel, ic=ic)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_average(eeg; channel)

Return the average signal of EEG channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_average(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    eeg_keep_channel!(eeg_new, channel=1)
    eeg_new.eeg_signals = @views s_average(eeg.eeg_signals[channel, :, :])
    eeg_new.eeg_header[:labels]=["averaged channel"]
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_average(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_average!(eeg; channel)

Return the average signal of EEG channels.  

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_average!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_average(eeg, channel=channel)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_average(eeg1, eeg2)

Return the average signal of all `eeg1` and `eeg2` channels.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_average(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG)

    size(eeg1.eeg_signals) == size(eeg2.eeg_signals) || throw(ArgumentError("Both signals must have the same size."))
    channel_n = eeg_channel_n(eeg1)
    epoch_n = eeg_epoch_n(eeg1)

    eeg_new = deepcopy(eeg1)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            eeg_new.eeg_signals[channel_idx, :, epoch_idx] = @views s2_average(signal1[channel_idx, :, epoch_idx], signal2[channel_idx, :, epoch_idx])
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg.eeg_header[:history], "eeg_average(EEG1, EEG2)")

    return eeg_new
end

"""
    eeg_invert_polarity(eeg; channel)

Invert polarity.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_invert_polarity(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)
    
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channel, :, :] = .- eeg_new.eeg_signals[channel, :, :]
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_invert_polarity(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_invert_polarity!(eeg; channel)

Invert polarity.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: channel(s) to invert
"""
function eeg_invert_polarity!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_invert_polarity(eeg, channel=channel)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_resample(eeg; new_sr)

Resample (up- or down-sample).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_resample(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    new_sr < 1 && throw(ArgumentError("new_sr must be ≥ 1."))
    new_sr > eeg_sr(eeg) && return eeg_upsample(eeg, new_sr=new_sr)
    new_sr < eeg_sr(eeg) && return eeg_downsample(eeg, new_sr=new_sr)
    new_sr == eeg_sr(eeg) && return eeg

end

"""
    eeg_resample!(eeg; new_sr)

Resample (up- or down-sample).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate
"""
function eeg_resample!(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    eeg_tmp = eeg_resample(eeg, new_sr=new_sr)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_upsample(eeg; new_sr)

Upsample EEG.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_upsample(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    (new_sr / eeg_sr(eeg) != new_sr ÷ eeg_sr(eeg) && verbose == true) && @info "New sampling rate should be easily captured by integer fractions, e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz."
    
    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate]):eeg.eeg_time[end]
    s_upsampled, t_upsampled = s_resample(eeg.eeg_signals, t=t, new_sr=new_sr)

    t_upsampled = collect(t_upsampled)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = s_upsampled
    eeg_new.eeg_time = t_upsampled
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_upsampled, 2) * size(s_upsampled, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_upsampled, 2) * size(s_upsampled, 3)) / new_sr
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_upsampled, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_upsampled, 2) / new_sr
    eeg_new.eeg_header[:sampling_rate] = new_sr
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_upsample(EEG, new_sr=$new_sr)")

    return eeg_new
end

"""
    eeg_upsample!(eeg; new_sr)

Upsample EEG.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_upsample!(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    eeg_tmp = eeg_upsample(eeg, new_sr=new_sr)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_time = eeg_tmp.eeg_time
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_downsample(eeg; new_sr)

Downsample EEG.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_downsample(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    (new_sr < eeg_sr(eeg) && verbose == true) && @info "To prevent aliasing due to down-sampling, a low-pass filter should be applied before removing data points. The filter cutoff should be the Nyquist frequency of the new down-sampled rate, ($(new_sr / 2) Hz), not the original Nyquist frequency ($(eeg_sr(eeg) / 2) Hz)."

    (new_sr / eeg_sr(eeg) != new_sr ÷ eeg_sr(eeg) && verbose == true) && @info "New sampling rate should be easily captured by integer fractions e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz."

    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate]):eeg.eeg_time[end]
    s_downsampled, t_downsampled = s_resample(eeg.eeg_signals, t=t, new_sr=new_sr)

    t_downsampled = collect(t_downsampled)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_time = t_downsampled
    eeg_new.eeg_signals = s_downsampled
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_downsampled, 2) * size(s_downsampled, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_downsampled, 2) * size(s_downsampled, 3)) / new_sr
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_downsampled, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_downsampled, 2) / new_sr
    eeg_new.eeg_header[:sampling_rate] = new_sr
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_downsample(EEG, new_sr=$new_sr)")

    return eeg_new
end

"""
    eeg_downsample!(eeg; new_sr)

Downsample EEG.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate
"""
function eeg_downsample!(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    eeg_tmp = eeg_downsample(eeg, new_sr=new_sr)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_time = eeg_tmp.eeg_time
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_wdenoise(eeg; channel, wt)

Perform wavelet denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_wdenoise(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), wt::T) where {T <: DiscreteWavelet}

    _check_channels(eeg, channel)
    channel_n = length(channel)
    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = @views s_wdenoise(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], wt=wt)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_wdenoise(EEG, wt=$wt)")

    return eeg_new
end

"""
    eeg_wdenoise!(eeg; channel, wt)

Perform wavelet denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
"""
function eeg_wdenoise!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), wt::T) where {T <: DiscreteWavelet}

    eeg_tmp = eeg_wdenoise(eeg, channel=channel, wt=wt)
    eeg_tmp = eeg_upsample(eeg, new_sr=new_sr)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_reference_a(eeg; type, med)

Reference to auricular (A1, A2) channels. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:l`:
    - `:l`: linked - average of A1 and A2
    - `:i`: ipsilateral - A1 for left channels, A2 for right channels
    - `:c`: contraletral - A1 for right channels, A2 for left channels
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_a(eeg::NeuroAnalyzer.EEG; type::Symbol=:l, med::Bool=false)

    _check_var(type, [:l, :i, :c], "type")
    all(iszero, occursin.("a1", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain A1 channel."))
    all(iszero, occursin.("a2", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain A2 channel."))

    # keep EEG channels
    channels = eeg_signal_channels(eeg)
    signal = @view eeg.eeg_signals[channels, :, :]

    a1_idx = findfirst(isequal("A1"), eeg.eeg_header[:labels])
    a2_idx = findfirst(isequal("A2"), eeg.eeg_header[:labels])
    a1 = eeg_extract_channel(eeg, channel=a1_idx)
    a2 = eeg_extract_channel(eeg, channel=a2_idx)

    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    s_ref = similar(signal)

    if type === :l
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(mean([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            @inbounds @simd for channel_idx in 1:channel_n
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :i
        central_picks = eeg_pick(eeg, pick=:central)
        Threads.@threads for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = @views vec(mean([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            else
                reference_channel = @views vec(median([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            end
            @inbounds @simd for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg, pick=:left)
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(a1[:, :, epoch_idx])
            @inbounds @simd for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg, pick=:right)
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(a2[:, :, epoch_idx])
            @inbounds @simd for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :c
        central_picks = eeg_pick(eeg, pick=:central)
        Threads.@threads for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = @views vec(mean([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            else
                reference_channel = @views vec(median([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            end
            @inbounds @simd for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg, pick=:left)
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(a2[:, :, epoch_idx])
            @inbounds @simd for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg, pick=:right)
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(a1[:, :, epoch_idx])
            @inbounds @simd for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    end

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channels, :, :] = s_ref
    eeg_new.eeg_header[:reference] = "A ($type)"
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_reference_a(EEG, type=$type, med=$med)")

    return eeg_new
end

"""
    eeg_reference_a!(eeg; type, med)

Reference to auricular (A1, A2) channels. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:l`:
    - `:l`: linked - average of A1 and A2
    - `:i`: ipsilateral - A1 for left channels, A2 for right channels
    - `:c`: contraletral - A1 for right channels, A2 for left channels
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_a!(eeg::NeuroAnalyzer.EEG; type::Symbol=:l, med::Bool=false)

    eeg_tmp = eeg_reference_a(eeg, type=type, med=med)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_reference_m(eeg; type, med)

Reference to mastoid (M1, M2) channels. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:l`:
    - `:l`: linked - average of M1 and M2
    - `:i`: ipsilateral - M1 for left channels, M2 for right channels
    - `:c`: contraletral - M1 for right channels, M2 for left channels
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_m(eeg::NeuroAnalyzer.EEG; type::Symbol=:l, med::Bool=false)

    _check_var(type, [:l, :i, :c], "type")
    all(iszero, occursin.("m1", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain M1 channel."))
    all(iszero, occursin.("m2", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain M2 channel."))

    # keep EEG channels
    channels = eeg_signal_channels(eeg)
    signal = @view eeg.eeg_signals[channels, :, :]

    m1_idx = findfirst(isequal("M1"), eeg.eeg_header[:labels])
    m2_idx = findfirst(isequal("M2"), eeg.eeg_header[:labels])
    m1 = eeg_extract_channel(eeg, channel=m1_idx)
    m2 = eeg_extract_channel(eeg, channel=m2_idx)

    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    s_ref = similar(signal)

    if type === :l
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(mean([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            @inbounds @simd for channel_idx in 1:channel_n
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :i
        central_picks = eeg_pick(eeg, pick=:central)
        Threads.@threads for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = @views vec(mean([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            else
                reference_channel = @views vec(median([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            end
            @inbounds @simd for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg, pick=:left)
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(m1[:, :, epoch_idx])
            @inbounds @simd for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg, pick=:right)
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(m2[:, :, epoch_idx])
            @inbounds @simd for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :c
        central_picks = eeg_pick(eeg, pick=:central)
        Threads.@threads for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = @views vec(mean([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            else
                reference_channel = @views vec(median([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            end
            @inbounds @simd for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg, pick=:left)
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(m2[:, :, epoch_idx])
            @inbounds @simd for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg, pick=:right)
        Threads.@threads for epoch_idx in 1:epoch_n
            reference_channel = @views vec(m1[:, :, epoch_idx])
            for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    end

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channels, :, :] = s_ref
    eeg_new.eeg_header[:reference] = "M ($type)"
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_reference_m(EEG, type=$type, med=$med)")

    return eeg_new
end

"""
    eeg_reference_m!(eeg; type, med)

Reference to mastoid (M1, M2) channels. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:l`:
    - `:l`: linked - average of M1 and M2
    - `:i`: ipsilateral - M1 for left channels, M2 for right channels
    - `:c`: contraletral - M1 for right channels, M2 for left channels
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_m!(eeg::NeuroAnalyzer.EEG; type::Symbol=:l, med::Bool=false)

    eeg_tmp = eeg_reference_m(eeg, type=type, med=med)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_fftdenoise(eeg; channel, pad, threshold)

Perform FFT denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `threshold::Int64=100`: PSD threshold for keeping frequency components

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_fftdenoise(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, threshold::Int64=100)

    _check_channels(eeg, channel)

    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in eachindex(channel)
                eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], _ = @views s_fftdenoise(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, threshold=threshold)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_fftdenoise(EEG, channel=$channel, pad=$pad, threshold=$threshold)")

    return eeg_new
end

"""
    eeg_fftdenoise!(eeg; channel, pad, threshold)

Perform FFT denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `threshold::Int64=100`: PSD threshold for keeping frequency components
"""
function eeg_fftdenoise!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, threshold::Int64=100)

    eeg_tmp = eeg_fftdenoise(eeg, channel=channel, pad=pad, threshold=threshold)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_reference_plap(eeg; nn, weights)

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `nn::Int64=4`: use `nn` adjacent electrodes
- `weights::Bool=false`: use mean of `nn` nearest channels if false; if true, mean of `nn` nearest channels is weighted by distance to the referenced channel
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_plap(eeg::NeuroAnalyzer.EEG; nn::Int64=4, weights::Bool=false, med::Bool=false)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    # keep EEG channels
    eeg_tmp = deepcopy(eeg)
    channels = eeg_signal_channels(eeg)
    signal = eeg.eeg_signals[channels, :, :]

    channel_n = size(signal, 1)
    nn < 1 && throw(ArgumentError("nn must be ≥ 1"))
    nn > channel_n - 1 && throw(ArgumentError("nn must be < $(channel_n - 1)"))
    epoch_n = size(signal, 3)
    
    loc_x = zeros(channel_n)
    loc_y = zeros(channel_n)
    for idx in 1:channel_n
        # loc_y[idx], loc_x[idx] = pol2cart(pi / 180 * eeg.eeg_locs[!, :loc_theta][idx], eeg.eeg_locs[!, :loc_radius][idx])
        # loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        loc_x[idx] = eeg.eeg_locs[idx, :loc_x]
        loc_y[idx] = eeg.eeg_locs[idx, :loc_y]
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)

    # Euclidean distance matrix
    d = zeros(channel_n, channel_n)
    for idx1 in 1:channel_n
        for idx2 in 1:channel_n
            d[idx1, idx2] = euclidean([loc_x[idx1], loc_y[idx1]], [loc_x[idx2], loc_y[idx2]])
        end
    end
    # set weights not to reference to itself
    d[d .== 0] .= Inf

    # nn nearest neighbors index matrix
    nn_idx = zeros(Int64, channel_n, nn)
    for idx1 in 1:channel_n
        # nn_idx[idx1, :] = sortperm(d[idx1, :])[2:(nn + 1)] # 1st neighbor is the electrode itself
        nn_idx[idx1, :] = sortperm(d[idx1, :])[1:nn]
    end

    s_ref = zeros(size(signal))

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            reference_channels = @view signal[nn_idx[channel_idx, :], :, epoch_idx]
            if weights == false
                if med == false
                    reference_channel = vec(mean(reference_channels, dims=1))
                else
                    reference_channel = vec(median(reference_channels, dims=1))
                end
            else
                g = Vector{Float64}()
                for idx1 in 1:nn
                    push!(g, 1 / d[channel_idx, nn_idx[channel_idx, idx1]] / sum(1 / d[channel_idx, nn_idx[channel_idx, :]]))
                end
                if med == false
                    reference_channel = vec(mean(g .* reference_channels, dims=1))
                else
                    reference_channel = vec(median(g .* reference_channels, dims=1))
                end
            end
            s_ref[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] .- reference_channel
        end
    end

    eeg_tmp.eeg_signals[channels, :, :] = s_ref
    eeg_tmp.eeg_header[:reference] = "PLAP ($nn)"
    eeg_reset_components!(eeg_tmp)
    push!(eeg_tmp.eeg_header[:history], "eeg_reference_plap(EEG, nn=$nn, med=$med))")

    return eeg_tmp
end

"""
    eeg_reference_plap!(eeg; nn, weights)

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `nn::Int64=4`: use `nn` adjacent electrodes
- `weights::Bool=false`: use distance weights; use mean of nearest channels if false
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_plap!(eeg::NeuroAnalyzer.EEG; nn::Int64=4, weights::Bool=false, med::Bool=false)

    eeg_tmp = eeg_reference_plap(eeg, nn=nn, weights=weights, med=med)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_zero(eeg)

Zero EEG channels at the beginning and at the end.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_zero(eeg::NeuroAnalyzer.EEG)

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[:, 1, :] .= 0
    eeg_new.eeg_signals[:, end, :] .= 0

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_zero(EEG)")

    return eeg_new
end

"""
    eeg_zero!(eeg)

Zero EEG channels at the beginning and at the end.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_zero!(eeg::NeuroAnalyzer.EEG)

    eeg_tmp = eeg_zero(eeg)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_wbp(eeg; channel, pad, frq, ncyc, demean)

Perform wavelet bandpass filtering.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_wbp(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, frq::Real, ncyc::Int64=6, demean::Bool=true)

    _check_channels(eeg, channel)

    epoch_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in eachindex(channel)
            eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = @views s_wbp(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, frq=frq, fs=fs, ncyc=ncyc, demean=demean)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_wbp(EEG, channel=$channel, pad=$pad, frq=$frq, ncyc=$ncyc, demean=$demean)")

    return eeg_new
end

"""
    eeg_wbp!(eeg; channel, pad, frq, ncyc, demean)

Perform wavelet bandpass filtering.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis
"""
function eeg_wbp!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, frq::Real, ncyc::Int64=6, demean::Bool=true)

    eeg_tmp = eeg_wbp(eeg, channel=channel, pad=pad, frq=frq, ncyc=ncyc, demean=demean)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_cbp(eeg; channel, pad, frq, demean)

Perform convolution bandpass filtering.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_cbp(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, frq::Real, demean::Bool=true)

    _check_channels(eeg, channel)

    epoch_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in eachindex(channel)
            eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = @views s_cbp(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, frq=frq, fs=fs, demean=demean)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_cbp(EEG, channel=$channel, pad=$pad, frq=$frq, demean=$demean)")

    return eeg_new
end

"""
    eeg_cbp!(eeg; channel, pad, frq, demean)

Perform convolution bandpass filtering.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Tuple{Real, Real}`: filter frequency
- `demean::Bool=true`: demean signal prior to analysis
"""
function eeg_cbp!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, frq::Real, demean::Bool=true)

    eeg_tmp = eeg_cbp(eeg, channel=channel, pad=pad, frq=frq, demean=demean)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_denoise_wien(eeg; channel)

Perform Wiener deconvolution denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns
- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_denoise_wien(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channel, :, :] = s_denoise_wien(eeg.eeg_signals[channel, :, :])
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_denoise_wien(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_denoise_wien!(eeg; channel)

Perform Wiener deconvolution denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_denoise_wien!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_denoise_wien(eeg, channel=channel)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_scale(eeg; channel, factor)

Multiply EEG channel(s) by `factor`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `factor::Real`: channel signal is multiplied by `factor`

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_scale(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), factor::Real)

    _check_channels(eeg, channel)
    
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channel, :, :] = @views eeg.eeg_signals[channel, :, :] .* factor
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_scale(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_scale!(eeg; channel, factor)

Multiply EEG channel(s) by `factor`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `factor::Real`: channel signal is multiplied by `factor`
"""
function eeg_scale!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), factor::Real)

    eeg_tmp = eeg_scale(eeg, channel=channel, factor=factor)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return nothing
end