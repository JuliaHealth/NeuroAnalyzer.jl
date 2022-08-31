"""
    eeg_reference_ch(eeg; channel, med)

Reference the `eeg` to specific `channel`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_ch(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}, med::Bool=false)

    _check_channels(eeg, channel)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        if length(channel) == 1
            reference_channel = @views vec(eeg.eeg_signals[channel, :, epoch_idx])
        else
            if med == false
                reference_channel = @views vec(mean(eeg.eeg_signals[channel, :, epoch_idx], dims=1))
            else
                reference_channel = @views vec(median(eeg.eeg_signals[channel, :, epoch_idx], dims=1))
            end
        end
        Threads.@threads for channel_idx in channels
            if channel_idx != channel
                @views eeg_new.eeg_signals[channel_idx, :, epoch_idx] .-= reference_channel
            end
        end
    end

    eeg_new.eeg_header[:reference] = "channel: $channel"
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_reference_ch(EEG, channel=$channel), med=$med")

    return eeg_new
end

"""
    eeg_reference_ch!(eeg; channel, med)

Reference the `eeg` to specific channel `channel`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_ch!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}, med::Bool=false)

    eeg.eeg_signals = eeg_reference_ch(eeg, channel=channel, med=med).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_reference_ch!(EEG, channel=$channel, med=$med)")

    return nothing
end

"""
    eeg_reference_car(eeg; exclude_fpo, exclude_current, med)

Reference the `eeg` to common average reference.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `exclude_fpo::Bool=true`: exclude Fp1, Fp2, O1, O2 from CAR calculation
- `exclude_current::Bool=true`: exclude current electrode from CAR calculation
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_car(eeg::NeuroAnalyzer.EEG; exclude_fpo::Bool=false, exclude_current::Bool=false, med::Bool=false)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        reference_channels = @view eeg.eeg_signals[:, :, epoch_idx]
        if exclude_fpo == true
            l = lowercase.(eeg_labels(eeg))
            "fp1" in l && (reference_channels = reference_channels[setdiff(1:end, (findfirst(isequal("fp1"), l))), :, :])
            "fp2" in l && (reference_channels = reference_channels[setdiff(1:end, (findfirst(isequal("fp2"), l))), :, :])
            "o1" in l && (reference_channels = reference_channels[setdiff(1:end, (findfirst(isequal("o1"), l))), :, :])
            "o2" in l && (reference_channels = reference_channels[setdiff(1:end, (findfirst(isequal("o2"), l))), :, :])
        end
        if exclude_current == false
            if med == false
                reference_channel = vec(mean(reference_channels, dims=1))
            else
                reference_channel = vec(median(reference_channels, dims=1))
            end
        end
        Threads.@threads for channel_idx in 1:channel_n
            if exclude_current == false
                reference_channels = reference_channels[setdiff(1:end, (channel_idx)), :, :]
                if med == false
                    reference_channel = vec(mean(reference_channels, dims=1))
                else
                    reference_channel = vec(median(reference_channels, dims=1))
                end
            end
            if channel_idx in channels
                @views eeg_new.eeg_signals[channel_idx, :, epoch_idx] = eeg.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    end

    eeg_new.eeg_header[:reference] = "CAR"
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_reference_car(EEG, exclude_fpo=$exclude_fpo, exclude_current=$exclude_current, med=$med))")

    return eeg_new
end

"""
    eeg_reference_car!(eeg; exclude_fpo, exclude_current, med)

Reference the `eeg` to common average reference.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `exclude_fpo::Bool=true`: exclude Fp1, Fp2, O1, O2 from CAR mean calculation
- `exclude_current::Bool=true`: exclude current electrode from CAR mean calculation
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_car!(eeg::NeuroAnalyzer.EEG; exclude_fpo::Bool=false, exclude_current::Bool=false, med::Bool=false)

    eeg.eeg_signals = eeg_reference_car(eeg, exclude_fpo=exclude_fpo, exclude_current=exclude_current, med=med).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_reference_car!(EEG, exclude_fpo=$exclude_fpo, exclude_current=$exclude_current, med=$med)")

    return nothing
end

"""
    eeg_derivative(eeg; channel)

Return the derivative of the `eeg` with length same as the signal.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_derivative(eeg::NeuroAnalyzer.EEG; channel::Int64=0)

    (channel < 0 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))"))
    if channel == 0
        channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    else
        channels = channel
    end
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    
    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if channel_idx in channels
                @views eeg_new.eeg_signals[channel_idx, :, epoch_idx] = s_derivative(eeg.eeg_signals[channel_idx, :, epoch_idx])
            end
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_derivative(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_derivative!(eeg; channel)

Return the derivative of the `eeg` with length same as the signal.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel
"""
function eeg_derivative!(eeg::NeuroAnalyzer.EEG; channel::Int64=0)

    eeg.eeg_signals = eeg_derivative(eeg).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_derivative!(EEG, channel=$channel)")

    return nothing
end

"""
    eeg_detrend(eeg; channel, type, offset, order, span)

Perform piecewise detrending of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel
- `type::Symbol`, optional
    - `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:linear`: linear trend is subtracted from `signal`
    - `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` is subtracted
    - `:loess`: fit and subtract loess approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `span::Float64=0.5`: smoothing of loess

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_detrend(eeg::NeuroAnalyzer.EEG; channel::Int64=0, type::Symbol=:linear, offset::Real=0, order::Int64=1, span::Float64=0.5)

    type in [:ls, :linear, :constant, :poly, :loess, :hp] || throw(ArgumentError("type must be :ls, :linear, :constant, :poly, :loess, :hp."))
    (channel < 0 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))"))
    if channel == 0
        channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    else
        channels = channel
    end
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if channel_idx in channels
                @views eeg_new.eeg_signals[channel_idx, :, epoch_idx] = s_detrend(eeg.eeg_signals[channel_idx, :, epoch_idx], type=type, offset=offset, order=order, span=span, fs=fs)
            end
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_detrend(EEG, channel=$channel, type=$type, offset=$offset, order=$order, span=$span)")

    return eeg_new
end

"""
    eeg_detrend!(eeg; channel, type, offset, order, span)

Remove linear trend from the `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel
- `type::Symbol`, optional
    - `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:linear`: linear trend is subtracted from `signal`
    - `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` order is subtracted
    - `:loess`: fit and subtract loess approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `span::Float64`: smoothing of loess
"""
function eeg_detrend!(eeg::NeuroAnalyzer.EEG; type::Symbol=:linear, offset::Real=0, order::Int64=1, span::Float64=0.5)

    eeg.eeg_signals = eeg_detrend(eeg, channel=channel, type=type, offset=offset, order=order, span=span).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_detrend!(EEG, channel=$channel, type=$type, offset=$offset, order=$order, span=$span)")

    return nothing
end

"""
    eeg_taper(eeg; channel, taper)

Taper `eeg` with `taper`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel
- `taper::Union{Vector{Real, Vector{ComplexF64}}``

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_taper(eeg::NeuroAnalyzer.EEG; channel::Int64=0, taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    (channel < 0 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))"))
    if channel == 0
        channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    else
        channels = channel
    end
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    s_tap = zeros(eltype(taper), size(eeg.eeg_signals))

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if channel_idx in channels
                @views s_tap[channel_idx, :, epoch_idx] = s_taper(eeg.eeg_signals[channel_idx, :, epoch_idx], taper=taper)
            else
                @views s_tap[channel_idx, :, epoch_idx] = eeg.eeg_signals[channel_idx, :, epoch_idx]
            end
        end
    end

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = s_tap
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_taper(EEG, taper=$taper, channel=$channel)")

    return eeg_new
end

"""
    eeg_taper!(eeg; channel, taper)

Taper `eeg` with `taper`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel
- `taper::Union{Vector{<:Real}, Vector{ComplexF64}}``
"""
function eeg_taper!(eeg::NeuroAnalyzer.EEG; taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    eeg.eeg_signals = eeg_taper(eeg, channel=channel, taper=taper).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_taper!(EEG, taper=$taper)")

    return nothing
end

"""
    eeg_demean(eeg; channel)

Remove mean value (DC offset).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_demean(eeg::NeuroAnalyzer.EEG; channel::Int64=0)

    (channel < 0 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))"))
    if channel == 0
        channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    else
        channels = channel
    end
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if channel_idx in channels
                @views eeg_new.eeg_signals[channel_idx, :, epoch_idx] = s_demean(eeg.eeg_signals[channel_idx, :, epoch_idx])
            end
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
- `channel::Int64=0`: if specified, process only this channel
"""
function eeg_demean!(eeg::NeuroAnalyzer.EEG; channel::Int64=0)

    eeg.eeg_signals = eeg_demean(eeg, channel=channel).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_demean!(EEG, channel=$channel)")

    return nothing
end

"""
    eeg_normalize(eeg; channel, method)

Normalize each `eeg` channel.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel
- `method::Symbol`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_normalize(eeg::NeuroAnalyzer.EEG; channel::Int64=0, method::Symbol)

    (channel < 0 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))"))
    if channel == 0
        channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    else
        channels = channel
    end
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if channel_idx in channels
                @views eeg_new.eeg_signals[channel_idx, :, epoch_idx] = s_normalize(eeg.eeg_signals[channel_idx, :, epoch_idx], method=method)
            end
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_normalize(EEG, channel=$channel, method=$method)")

    return eeg_new
end

"""
    eeg_normalize!(eeg; channel, method)

Normalize each `eeg` channel.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel
- `method::Symbol`
"""
function eeg_normalize!(eeg::NeuroAnalyzer.EEG; channel::Int64=0, method::Symbol)

    eeg.eeg_signals = eeg_normalize(eeg, channel=channel, method=method).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_normalize!(EEG, channel=$channel, method=$method)")

    return nothing
end

"""
    eeg_add_noise(eeg; channel)

Add random noise to the `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_add_noise(eeg::NeuroAnalyzer.EEG; channel::Int64=0)

    (channel < 0 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))"))
    if channel == 0
        channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    else
        channels = channel
    end
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if channel_idx in channels
                @views eeg_new.eeg_signals[channel_idx, :, epoch_idx] = s_add_noise(eeg.eeg_signals[channel_idx, :, epoch_idx])
            end
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_add_noise(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_add_noise!(eeg; channel)

Add random noise to the `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, process only this channel
"""
function eeg_add_noise!(eeg::NeuroAnalyzer.EEG; channel::Int64=0)

    eeg.eeg_signals = eeg_add_noise(eeg, channel=channel).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_add_noise!(EEG, channel=$channel)")

    return nothing
end

"""
    eeg_filter(eeg; <keyword arguments>)

Apply filtering to `eeg` channels. By default it filters all signal (EEG/MEG) channels. To filter other channel type, use `channel` parameter.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, filter only this channel
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`
    - `:remez`
    - `:mavg`: moving average (with threshold)
    - `:mmed`: moving median (with threshold)
    - `:poly`: polynomial of `order` order
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64=8`: filter order, number of taps for :remez filter, k-value for :mavg and :mmed (window length = 2 × k + 1)
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
- `bw::Real=-1`: bandwidth for :iirnotch and :remez filters
- `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass), for causal filter use :onepass
- `t::Real`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{<:Real}, Nothing} - window, required for FIR filter

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_filter(eeg::NeuroAnalyzer.EEG; channel::Int64=0, fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple}=0, fs::Int64=0, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, dir::Symbol=:twopass, d::Int64=1, t::Real=0, window::Union{Vector{<:Real}, Nothing}=nothing)

    (channel < 0 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))"))
    if channel == 0
        channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    else
        channels = channel
    end
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)
    
    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if channel_idx in channels
                eeg_new.eeg_signals[channel_idx, :, epoch_idx] = @views s_filter(eeg.eeg_signals[channel_idx, :, epoch_idx], fprototype=fprototype, ftype=ftype, cutoff=cutoff, fs=fs, order=order, rp=rp, rs=rs, bw=bw, dir=dir, t=t, window=window)
            end
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_filter(EEG, fprototype=$fprototype, ftype=$ftype, cutoff=$cutoff, order=$order, rp=$rp, rs=$rs, dir=$dir, t=$t, window=$window)")

    return eeg_new
end

"""
    eeg_filter!(eeg; <keyword arguments>)

Apply filtering to `eeg` channels. By default it filters all signal (EEG/MEG) channels. To filter other channel type, use `channel` parameter.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64=0`: if specified, filter only this channel
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`
    - `:remez`
    - `:mavg`: moving average (with threshold)
    - `:mmed`: moving median (with threshold)
    - `:poly`: polynomial of `order` order
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64=8`: filter order, number of taps for :remez filter, k-value for :mavg and :mmed (window length = 2 × k + 1)
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
- `bw::Real=-1`: bandwidth for :iirnotch and :remez filters
- `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass), for causal filter use :onepass
- `t::Real`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{<:Real}, Nothing} - window, required for FIR filter
"""
function eeg_filter!(eeg::NeuroAnalyzer.EEG; channel::Int64=0, fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple}=0, fs::Int64=0, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, dir::Symbol=:twopass, t::Real=0, window::Union{Vector{<:Real}, Nothing}=nothing)

    eeg.eeg_signals = eeg_filter(eeg,
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
                                 window=window).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_filter!(EEG, channel=$channel, fprototype=$fprototype, ftype=$ftype, cutoff=$cutoff, order=$order, rp=$rp, rs=$rs, dir=$dir, window=$window)")

    return nothing
end

"""
    eeg_pca(eeg; n)

Calculate `n` first PCs for `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `n::Int64`: number of PCs

# Returns

Named tuple containing:
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
- `pc_m::PCA{Float64}`: PC mean
"""
function eeg_pca(eeg::NeuroAnalyzer.EEG; n::Int64)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]

    pc, pc_var, pc_m = s_pca(signal, n=n)

    return (pc=pc, pc_var=pc_var, pc_m=pc_m)
end

"""
    eeg_pca_reconstruct(eeg)

Reconstruct `eeg` signals using PCA components.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_pca_reconstruct(eeg::NeuroAnalyzer.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :pc in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain :pc component. Perform eeg_pca(EEG) first."))
    :pc_m in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain :pc_m component. Perform eeg_pca(EEG) first."))

    eeg_new = deepcopy(eeg)
    pc_idx = findfirst(isequal(:pc), eeg.eeg_header[:components])
    pc_m_idx = findfirst(isequal(:pc_m), eeg.eeg_header[:components])
    eeg_new.eeg_signals = s_pca_reconstruct(eeg_new.eeg_signals, pc=eeg_new.eeg_components[pc_idx], pc_m=eeg_new.eeg_components[pc_m_idx])
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_pca_reconstruct(EEG)")

    return eeg_new
end

"""
    eeg_pca_reconstruct!(eeg)

Reconstruct `eeg` signals using PCA components.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_pca_reconstruct!(eeg::NeuroAnalyzer.EEG)

    eeg.eeg_signals = eeg_pca_reconstruct(eeg).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_pca_reconstruct!(EEG)")

    return nothing
end

"""
    eeg_ica(eeg; <keyword arguments>)

Calculate `n` first ICs for `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `n::Int64`: number of ICs
- `tol::Float64=1.0e-6`: tolerance for ICA
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor: :tanh, :gaus

# Returns

Named tuple containing:
- `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
"""
function eeg_ica(eeg::NeuroAnalyzer.EEG; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    ic, ic_mw = s_ica(eeg.eeg_signals, n=n, tol=tol, iter=iter, f=f)

    return (ic=ic, ic_mw=ic_mw)
end

"""
    eeg_average(eeg)

Return the average signal of all `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_average(eeg::NeuroAnalyzer.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    eeg_new = deepcopy(eeg)
    eeg_keep_channel!(eeg_new, channel=1)
    eeg_new.eeg_signals = s_average(eeg.eeg_signals)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_average(EEG)")

    return eeg_new
end

"""
    eeg_average!(eeg)

Return the average signal of all `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_average!(eeg::NeuroAnalyzer.EEG)

    eeg.eeg_signals = eeg_average(eeg).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_average!(EEG)")

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
    s_averaged = zeros(size(eeg1.eeg_signals))

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s_averaged[channel_idx, :, epoch_idx] = @views s2_average(signal1[channel_idx, :, epoch_idx], signal2[channel_idx, :, epoch_idx])
        end
    end

    eeg_new = deepcopy(eeg1)
    eeg_new.eeg_signals = s_averaged
    eeg_reset_components!(eeg_new)
    push!(eeg.eeg_header[:history], "eeg_average(EEG1, EEG2)")

    return eeg_new
end

"""
    eeg_ica_reconstruct(eeg; ica)

Reconstruct `eeg` signals using removal of `ica` ICA components.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `ica::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_ica_reconstruct(eeg::NeuroAnalyzer.EEG; ica::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :ica in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain :ica component. Perform eeg_ica(EEG) first."))
    :ica_mw in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain :ica_mw component. Perform eeg_ica(EEG) first."))

    eeg_new = deepcopy(eeg)
    ica_a_idx = findfirst(isequal(:ica), eeg.eeg_header[:components])
    ica_mw_idx = findfirst(isequal(:ica_mw), eeg.eeg_header[:components])
    eeg_new.eeg_signals = s_ica_reconstruct(eeg_new.eeg_signals, ic=eeg_new.eeg_components[ica_a_idx], ic_mw=eeg_new.eeg_components[ica_mw_idx], ic_v=ica)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_ica_reconstruct(EEG, ica=$ica)")

    return eeg_new
end

"""
    eeg_ica_reconstruct!(eeg; ica)

Reconstruct `eeg` signals using removal of `ica` ICA components.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `ica::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove
"""
function eeg_ica_reconstruct!(eeg::NeuroAnalyzer.EEG; ica::Union{Int64, Vector{Int64}, AbstractRange})

    eeg.eeg_signals = eeg_ica_reconstruct(eeg, ica=ica).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_ica_reconstruct!(EEG, ica=$ica)")

    return nothing
end

"""
    eeg_invert_polarity(eeg; channel)

Invert polarity of `channel` of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`: channel to invert

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_invert_polarity(eeg::NeuroAnalyzer.EEG; channel::Int64)

    (channel < 1 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channel, :, :] = .- eeg_new.eeg_signals[channel, :, :]
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_invert_polarity(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_invert_polarity!(eeg; channel)

Invert polarity of `channel` of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to invert
"""
function eeg_invert_polarity!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    (channel < 1 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))

    eeg.eeg_signals[channel, :, :] = .- eeg.eeg_signals[channel, :, :]
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_invert_polarity!(EEG, channel=$channel)")

    return nothing
end

"""
    eeg_resample(eeg; new_sr)

Resample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_resample(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    new_sr > eeg_sr(eeg) && (eeg_new = eeg_upsample(eeg, new_sr=new_sr))
    new_sr < eeg_sr(eeg) && (eeg_new = eeg_downsample(eeg, new_sr=new_sr))
    new_sr == eeg_sr(eeg) && (eeg_new = eeg)

    return eeg_new
end

"""
    eeg_resample!(eeg; new_sr)

Resample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate
"""
function eeg_resample!(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    new_sr > eeg_sr(eeg) && eeg_upsample!(eeg, new_sr=new_sr)
    new_sr < eeg_sr(eeg) && eeg_downsample!(eeg, new_sr=new_sr)

    return nothing
end

"""
    eeg_upsample(eeg; new_sr)

Upsample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_upsample(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    new_sr / eeg_sr(eeg) != new_sr ÷ eeg_sr(eeg) && (@info "New sampling rate should be easily captured by integer fractions, e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.")
    
    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    s_upsampled, t_upsampled = s_resample(eeg.eeg_signals, t=t, new_sr=new_sr)

    t_upsampled = collect(t_upsampled)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = s_upsampled
    eeg_new.eeg_time = t_upsampled
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_upsampled, 2) * size(s_upsampled, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_upsampled, 2) * size(s_upsampled, 3)) / new_sr
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_upsampled, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_upsampled, 2) / new_sr
    eeg_new.eeg_header[:sampling_rate] = repeat([new_sr], eeg_channel_n(eeg_new))
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_upsample(EEG, new_sr=$new_sr)")

    return eeg_new
end

"""
    eeg_upsample!(eeg; new_sr)

Upsample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_upsample!(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    new_sr / eeg_sr(eeg) != new_sr ÷ eeg_sr(eeg) && (@info "New sampling rate should be easily captured by integer fractions e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.")

    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    eeg.eeg_signals, t_upsampled = s_resample(eeg.eeg_signals, t=t, new_sr=new_sr)

    eeg.eeg_time = collect(t_upsampled)
    eeg.eeg_header[:eeg_duration_samples] = eeg_signal_len(eeg) * eeg_epoch_n(eeg)
    eeg.eeg_header[:eeg_duration_seconds] = (eeg_signal_len(eeg) * eeg_epoch_n(eeg)) / new_sr
    eeg.eeg_header[:epoch_duration_seconds] = eeg_signal_len(eeg)
    eeg.eeg_header[:epoch_duration_seconds] = eeg_signal_len(eeg) / new_sr
    eeg.eeg_header[:sampling_rate] = repeat([new_sr], eeg.eeg_header[:channel_n])
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_upsample!(EEG, new_sr=$new_sr)")

    return nothing
end

"""
    eeg_downsample(eeg; new_sr)

Downsample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_downsample(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    new_sr < eeg_sr(eeg) && (@info "To prevent aliasing due to down-sampling, a low-pass filter should be applied before removing data points. The filter cutoff should be the Nyquist frequency of the new down-sampled rate, ($(new_sr / 2) Hz), not the original Nyquist frequency ($(eeg_sr(eeg) / 2) Hz).")

    new_sr / eeg_sr(eeg) != new_sr ÷ eeg_sr(eeg) && (@info "New sampling rate should be easily captured by integer fractions e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.")

    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    s_downsampled, t_downsampled = s_resample(eeg.eeg_signals, t=t, new_sr=new_sr)

    t_downsampled = collect(t_downsampled)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_time = t_downsampled
    eeg_new.eeg_signals = s_downsampled
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_downsampled, 2) * size(s_downsampled, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_downsampled, 2) * size(s_downsampled, 3)) / new_sr
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_downsampled, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_downsampled, 2) / new_sr
    eeg_new.eeg_header[:sampling_rate] = repeat([new_sr], eeg_channel_n(eeg_new))
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_downsample(EEG, new_sr=$new_sr)")

    return eeg_new
end

"""
    eeg_downsample!(eeg; new_sr)

Downsample all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `new_sr::Int64`: new sampling rate
"""
function eeg_downsample!(eeg::NeuroAnalyzer.EEG; new_sr::Int64)

    new_sr < eeg_sr(eeg) && (@info "To prevent aliasing due to down-sampling, a low-pass filter should be applied before removing data points. The filter cutoff should be the Nyquist frequency of the new down-sampled rate, ($(new_sr / 2) Hz), not the original Nyquist frequency ($(eeg_sr(eeg) / 2) Hz).")

    new_sr / eeg_sr(eeg) != new_sr ÷ eeg_sr(eeg) && (@info "New sampling rate should be easily captured by integer fractions e.g. 1000 Hz → 250 Hz or 256 Hz → 512 Hz.")

    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    eeg.eeg_signals, t_downsampled = s_resample(eeg.eeg_signals, t=t, new_sr=new_sr)

    eeg.eeg_time = collect(t_downsampled)
    eeg.eeg_header[:eeg_duration_samples] = eeg_signal_len(eeg) * eeg_epoch_n(eeg)
    eeg.eeg_header[:eeg_duration_seconds] = (eeg_signal_len(eeg) * eeg_epoch_n(eeg)) / new_sr
    eeg.eeg_header[:epoch_duration_samples] = eeg_signal_len(eeg)
    eeg.eeg_header[:epoch_duration_seconds] = eeg_signal_len(eeg) / new_sr
    eeg.eeg_header[:sampling_rate] = repeat([new_sr], eeg.eeg_header[:channel_n])
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_downsample!(EEG, new_sr=$new_sr)")

    return nothing
end

"""
    eeg_wdenoise(eeg; wt)

Perform wavelet denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `wt::Symbol=:db4`: wavelet type: :db2, :db4, :db8, :db10, :haar, :coif2, :coif4, :coif8

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_wdenoise(eeg::NeuroAnalyzer.EEG; wt::Symbol=:db4)

    eeg_new = deepcopy(eeg)
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    s_denoised = similar(eeg.eeg_signals)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s_denoised[channel_idx, :, epoch_idx] = @views s_wdenoise(eeg.eeg_signals[channel_idx, :, epoch_idx], wt=wt)
        end
    end

    eeg_new.eeg_signals = s_denoised
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_wdenoise(EEG, wt=$wt)")

    return eeg_new
end

"""
    eeg_wdenoise!(eeg; wt)

Perform wavelet denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `wt::Symbol=:db4`: wavelet type: db2, db4, db8, db10, haar
"""
function eeg_wdenoise!(eeg::NeuroAnalyzer.EEG; wt::Symbol=:db4)

    channel_n = eeg_channel_(eeg)
    epoch_n = eeg_epoch_n(eeg)

    s_denoised = similar(eeg.eeg_signals)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s_denoised[channel_idx, :, epoch_idx] = @views s_wdenoise(eeg.eeg_signals[channel_idx, :, epoch_idx], wt=wt)
        end
    end

    eeg.eeg_signals = s_denoised
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_wdenoise!(EEG, wt=$wt)")

    return nothing
end

"""
    eeg_reference_a(eeg; type, med)

Reference the `eeg` to auricular channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:link`: :l (linked, average of A1 and A2), :i (ipsilateral, A1 for left channels) or :c (contraletral, A1 for right channels)
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_a(eeg::NeuroAnalyzer.EEG; type::Symbol=:l, med::Bool=false)

    type in [:l, :i, :c] || throw(ArgumentError("type must be :l, :i, :c."))
    all(iszero, occursin.("a1", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain A1 channel."))
    all(iszero, occursin.("a2", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain A2 channel."))

    eeg_tmp = deepcopy(eeg)
    channel_labels = eeg_labels(eeg)
    a1_idx = findfirst(isequal("A1"), eeg_tmp.eeg_header[:labels])
    a2_idx = findfirst(isequal("A2"), eeg_tmp.eeg_header[:labels])
    a1 = eeg_extract_channel(eeg, channel=a1_idx)
    a2 = eeg_extract_channel(eeg, channel=a2_idx)
    eeg_delete_channel!(eeg_tmp, channel=a2_idx)
    eeg_delete_channel!(eeg_tmp, channel=a1_idx)
    eeg_channel_n(eeg_tmp, type=:eeg) < eeg_channel_n(eeg_tmp, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg_tmp)
    epoch_n = eeg_epoch_n(eeg_tmp)
    s_ref = zeros(size(eeg_tmp.eeg_signals))

    if type === :l
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(mean([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            Threads.@threads for channel_idx in 1:channel_n
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :i
        central_picks = eeg_pick(eeg_tmp, pick=:central)
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg_tmp, pick=:left)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(a1[:, :, epoch_idx])
            Threads.@threads for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg_tmp, pick=:right)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(a2[:, :, epoch_idx])
            Threads.@threads for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :c
        central_picks = eeg_pick(eeg_tmp, pick=:central)
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg_tmp, pick=:left)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(a2[:, :, epoch_idx])
            Threads.@threads for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg_tmp, pick=:right)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(a1[:, :, epoch_idx])
            Threads.@threads for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    end

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = vcat(s_ref, a1, a2)
    eeg_new.eeg_header[:labels] = eeg_labels(eeg_tmp)
    push!(eeg_new.eeg_header[:labels], "A1")
    push!(eeg_new.eeg_header[:labels], "A2")
    eeg_new.eeg_header[:reference] = "A ($type)"
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_reference_a(EEG, type=$type, med=$med)")

    return eeg_new
end

"""
    eeg_reference_a!(eeg; type, med)

Reference the `eeg` to auricular channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:link`: :l (linked, average of A1 and A2), :i (ipsilateral, A1 for left channels) or :c (contraletral, A1 for right channels)
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_a!(eeg::NeuroAnalyzer.EEG; type::Symbol=:l, med::Bool=false)

    type in [:l, :i, :c] || throw(ArgumentError("type must be :l, :i, :c."))
    all(iszero, occursin.("a1", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain A1 channel."))
    all(iszero, occursin.("a2", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain A2 channel."))

    eeg_tmp = deepcopy(eeg)
    channel_labels = eeg_labels(eeg)
    a1_idx = findfirst(isequal("A1"), eeg_tmp.eeg_header[:labels])
    a2_idx = findfirst(isequal("A2"), eeg_tmp.eeg_header[:labels])
    a1 = eeg_extract_channel(eeg, channel=a1_idx)
    a2 = eeg_extract_channel(eeg, channel=a2_idx)
    eeg_delete_channel!(eeg_tmp, channel=a2_idx)
    eeg_delete_channel!(eeg_tmp, channel=a1_idx)
    eeg_channel_n(eeg_tmp, type=:eeg) < eeg_channel_n(eeg_tmp, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg_tmp)
    epoch_n = eeg_epoch_n(eeg_tmp)
    s_ref = zeros(size(eeg_tmp.eeg_signals))

    if type === :l
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in 1:channel_n
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :i
        central_picks = eeg_pick(eeg_tmp, pick=:central)
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg_tmp, pick=:left)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(a1[:, :, epoch_idx])
            Threads.@threads for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg_tmp, pick=:right)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(a2[:, :, epoch_idx])
            Threads.@threads for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :c
        central_picks = eeg_pick(eeg_tmp, pick=:central)
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([a1[:, :, epoch_idx], a2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg_tmp, pick=:left)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(a2[:, :, epoch_idx])
            Threads.@threads for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg_tmp, pick=:right)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(a1[:, :, epoch_idx])
            Threads.@threads for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    end

    eeg.eeg_signals = vcat(s_ref, a1, a2)
    eeg.eeg_header[:labels] = eeg_labels(eeg_tmp)
    push!(eeg.eeg_header[:labels], "A1")
    push!(eeg.eeg_header[:labels], "A2")

    eeg.eeg_header[:reference] = "A ($type)"
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_reference_a!(EEG, type=$type, med=$med)")

    return nothing
end

"""
    eeg_reference_m(eeg; type, med)

Reference the `eeg` to mastoid channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:link`: :l (linked, average of M1 and M2), :i (ipsilateral, M1 for left channels) or :c (contraletral, M1 for right channels)
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_m(eeg::NeuroAnalyzer.EEG; type::Symbol=:l, med::Bool=false)

    type in [:l, :i, :c] || throw(ArgumentError("type must be :l, :i, :c."))
    all(iszero, occursin.("m1", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain M1 channel."))
    all(iszero, occursin.("m2", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain M2 channel."))

    eeg_tmp = deepcopy(eeg)
    channel_labels = eeg_labels(eeg)
    m1_idx = findfirst(isequal("M1"), eeg_tmp.eeg_header[:labels])
    m2_idx = findfirst(isequal("M2"), eeg_tmp.eeg_header[:labels])
    m1 = eeg_extract_channel(eeg, channel=m1_idx)
    m2 = eeg_extract_channel(eeg, channel=m2_idx)
    eeg_delete_channel!(eeg_tmp, channel=m2_idx)
    eeg_delete_channel!(eeg_tmp, channel=m1_idx)
    eeg_channel_n(eeg_tmp, type=:eeg) < eeg_channel_n(eeg_tmp, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg_tmp)
    epoch_n = eeg_epoch_n(eeg_tmp)
    s_ref = zeros(size(eeg_tmp.eeg_signals))

    if type === :l
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in 1:channel_n
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :i
        central_picks = eeg_pick(eeg_tmp, pick=:central)
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg_tmp, pick=:left)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(m1[:, :, epoch_idx])
            Threads.@threads for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg_tmp, pick=:right)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(m2[:, :, epoch_idx])
            Threads.@threads for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :c
        central_picks = eeg_pick(eeg_tmp, pick=:central)
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg_tmp, pick=:left)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(m2[:, :, epoch_idx])
            Threads.@threads for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg_tmp, pick=:right)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(m1[:, :, epoch_idx])
            Threads.@threads for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    end

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = vcat(s_ref, m1, m2)
    eeg_new.eeg_header[:labels] = eeg_labels(eeg_tmp)
    push!(eeg_new.eeg_header[:labels], "M1")
    push!(eeg_new.eeg_header[:labels], "M2")
    eeg_new.eeg_header[:reference] = "M ($type)"
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_reference_m(EEG, type=$type, med=$med)")

    return eeg_new
end

"""
    eeg_reference_m!(eeg; type, med)

Reference the `eeg` to mastoid channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:link`: :l (linked, average of M1 and M2), :i (ipsilateral, M1 for left channels) or :c (contraletral, M1 for right channels)
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_m!(eeg::NeuroAnalyzer.EEG; type::Symbol=:lm, med::Bool=false)

    type in [:l, :i, :c] || throw(ArgumentError("type must be :l, :i, :c."))
    all(iszero, occursin.("m1", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain M1 channel."))
    all(iszero, occursin.("m2", lowercase.(eeg.eeg_header[:labels]))) == false || throw(ArgumentError("EEG does not contain M2 channel."))

    eeg_tmp = deepcopy(eeg)
    channel_labels = eeg_labels(eeg)
    m1_idx = findfirst(isequal("M1"), eeg_tmp.eeg_header[:labels])
    m2_idx = findfirst(isequal("M2"), eeg_tmp.eeg_header[:labels])
    m1 = eeg_extract_channel(eeg, channel=m1_idx)
    m2 = eeg_extract_channel(eeg, channel=m2_idx)
    eeg_delete_channel!(eeg_tmp, channel=m2_idx)
    eeg_delete_channel!(eeg_tmp, channel=m1_idx)
    eeg_channel_n(eeg_tmp, type=:eeg) < eeg_channel_n(eeg_tmp, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg_tmp)
    epoch_n = eeg_epoch_n(eeg_tmp)
    s_ref = zeros(size(eeg_tmp.eeg_signals))

    if type === :l
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in 1:channel_n
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :i
        central_picks = eeg_pick(eeg_tmp, pick=:central)
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg_tmp, pick=:left)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(m1[:, :, epoch_idx])
            Threads.@threads for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg_tmp, pick=:right)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(m2[:, :, epoch_idx])
            Threads.@threads for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    elseif type === :c
        central_picks = eeg_pick(eeg_tmp, pick=:central)
        @inbounds @simd for epoch_idx in 1:epoch_n
            if med == false
                reference_channel = vec(mean([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            else
                reference_channel = vec(median([m1[:, :, epoch_idx], m2[:, :, epoch_idx]]))
            end
            Threads.@threads for channel_idx in central_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        left_picks = eeg_pick(eeg_tmp, pick=:left)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(m2[:, :, epoch_idx])
            Threads.@threads for channel_idx in left_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
        right_picks = eeg_pick(eeg_tmp, pick=:right)
        @inbounds @simd for epoch_idx in 1:epoch_n
            reference_channel = vec(m1[:, :, epoch_idx])
            Threads.@threads for channel_idx in right_picks
                s_ref[channel_idx, :, epoch_idx] = @views eeg_tmp.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
            end
        end
    end

    eeg.eeg_signals = vcat(s_ref, m1, m2)
    eeg.eeg_header[:labels] = eeg_labels(eeg_tmp)
    push!(eeg.eeg_header[:labels], "M1")
    push!(eeg.eeg_header[:labels], "M2")

    eeg.eeg_header[:reference] = "M ($type)"
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_reference_m!(EEG, type=$type)")

    return nothing
end

"""
    eeg_fftdenoise(eeg; pad, threshold)

Perform wavelet denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pad::Int64=0`: pad signal with `pad` zeros
- `threshold::Int64=100`: PSD threshold for keeping frequency components

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_fftdenoise(eeg::NeuroAnalyzer.EEG; pad::Int64=0, threshold::Int64=100)

    eeg_new = deepcopy(eeg)
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    s_denoised = similar(eeg.eeg_signals)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s_denoised[channel_idx, :, epoch_idx] = @views s_fftdenoise(eeg.eeg_signals[channel_idx, :, epoch_idx], pad=pad, threshold=threshold)
        end
    end

    eeg_new.eeg_signals = s_denoised
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_fftdenoise(EEG, pad=$pad, threshold=$threshold)")

    return eeg_new
end

"""
    eeg_fftdenoise!(eeg; pad, threshold)

Perform wavelet denoising.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pad::Int64=0`: pad signal with `pad` zeros
- `threshold::Int64=100`: PSD threshold for keeping frequency components
"""
function eeg_fftdenoise!(eeg::NeuroAnalyzer.EEG; pad::Int64=0, threshold::Int64=100)

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    s_denoised = similar(eeg.eeg_signals)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s_denoised[channel_idx, :, epoch_idx] = @views s_fftdenoise(eeg.eeg_signals[channel_idx, :, epoch_idx], pad=pad, threshold=threshold)
        end
    end

    eeg.eeg_signals = s_denoised
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_fftdenoise!(EEG, pad=$pad, threshold=$threshold)")

    return nothing
end

"""
    eeg_reference_plap(eeg; nn, weights)

Reference the `eeg` using planar Laplacian (using `nn` adjacent electrodes).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `nn::Int64=4`: number of nearest electrodes
- `weights::Bool=true`: use distance weights; use mean of nearest channels if false
- `med::Bool=false`: use median instead of mean

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reference_plap(eeg::NeuroAnalyzer.EEG; nn::Int64=4, weights::Bool=true, med::Bool=false)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))
    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    nn < channel_n - 1 || throw(ArgumentError("nn must be < $(channel_n - 1)"))
    epoch_n = eeg_epoch_n(eeg)
    
    loc_x = zeros(channel_n)
    loc_y = zeros(channel_n)
    for idx in 1:channel_n
        loc_y[idx], loc_x[idx] = pol2cart(pi / 180 * eeg.eeg_header[:loc_theta][idx],
                                          eeg.eeg_header[:loc_radius][idx])
    end
    # Euclidean distance matrix
    d = zeros(channel_n, channel_n)
    for idx1 in 1:channel_n
        for idx2 in 1:channel_n
            d[idx1, idx2] = euclidean([loc_x[idx1], loc_y[idx1]], [loc_x[idx2], loc_y[idx2]])
        end
    end
    # nn nearest neighbors index matrix
    nn_idx = zeros(Int64, channel_n, nn)
    for idx1 in 1:channel_n
        nn_idx[idx1, :] = sortperm(d[idx1, :])[2:(nn + 1)] # 1st neighbor is the electrode itself
    end

    s_ref = zeros(size(eeg.eeg_signals))

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            reference_channels = @view eeg.eeg_signals[nn_idx[channel_idx, :], :, epoch_idx]
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
                reference_channel = vec(sum(g .* reference_channels, dims=1))
            end
            s_ref[channel_idx, :, epoch_idx] = @views eeg.eeg_signals[channel_idx, :, epoch_idx] .- reference_channel
        end
    end

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = s_ref
    eeg_new.eeg_header[:reference] = "PLAP ($nn)"
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_reference_plap(EEG, nn=$nn, med=$med))")

    return eeg_new
end

"""
    eeg_reference_plap!(eeg; nn, weights)

Reference the `eeg` using planar Laplacian (using `nn` adjacent electrodes).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `nn::Int64=4`: number of nearest electrodes
- `weights::Bool=true`: use distance weights; use mean of nearest channels if false
- `med::Bool=false`: use median instead of mean
"""
function eeg_reference_plap!(eeg::NeuroAnalyzer.EEG; nn::Int64=4, weights::Bool=true, med::Bool=false)

    eeg.eeg_signals = eeg_reference_plap(eeg, nn=nn, weights=weights, med=med).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_reference_plap!(EEG, nn=$nn, med=$med)")

    return nothing
end

"""
    eeg_zero(eeg)

Zero `eeg` channels at the beginning and at the end.

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

Zero `eeg` channel at the beginning and at the end.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_zero!(eeg::NeuroAnalyzer.EEG)

    eeg.eeg_signals[:, 1, :] .= 0
    eeg.eeg_signals[:, end, :] .= 0

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_zero(EEG)")

    return nothing
end

"""
    eeg_wbp(eeg; pad, frq, ncyc, demean)

Perform wavelet bandpass filtering of the `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_wbp(eeg::NeuroAnalyzer.EEG; pad::Int64=0, frq::Real, ncyc::Int64=6, demean::Bool=true)

    eeg_new = deepcopy(eeg)

    epoch_n = eeg_epoch_n(eeg)
    channel_n = eeg_channel_n(eeg)
    fs = eeg_sr(eeg)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            eeg_new.eeg_signals[channel_idx, :, epoch_idx] = @views s_wbp(eeg.eeg_signals[channel_idx, :, epoch_idx], pad=pad, frq=frq, fs=fs, ncyc=ncyc, demean=demean)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_wbp(EEG, pad=$pad, frq=$frq, ncyc=$ncyc, demean=$demean)")

    return eeg_new
end

"""
    eeg_wbp!(eeg; pad, frq, ncyc, demean)

Perform wavelet bandpass filtering of the `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis
"""
function eeg_wbp!(eeg::NeuroAnalyzer.EEG; pad::Int64=0, frq::Real, ncyc::Int64=6, demean::Bool=true)

    eeg.eeg_signals = eeg_wbp(eeg, pad=pad, frq=frq, ncyc=ncyc, demean=demean).eeg_signals

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_wbp!(EEG, pad=$pad, frq=$frq, ncyc=$ncyc, demean=$demean)")

    return nothing
end

"""
    eeg_cbp(eeg; pad, frq, demean)

Perform convolution bandpass filtering of the `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_cbp(eeg::NeuroAnalyzer.EEG; pad::Int64=0, frq::Real, demean::Bool=true)

    eeg_new = deepcopy(eeg)

    epoch_n = eeg_epoch_n(eeg)
    channel_n = eeg_channel_n(eeg)
    fs = eeg_sr(eeg)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            eeg_new.eeg_signals[channel_idx, :, epoch_idx] = @views s_cbp(eeg.eeg_signals[channel_idx, :, epoch_idx], pad=pad, frq=frq, fs=fs, demean=demean)
        end
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_cbp(EEG, pad=$pad, frq=$frq, demean=$demean)")

    return eeg_new
end

"""
    eeg_cbp!(eeg; pad, frq, ncyc, demean)

Perform convolution bandpass filtering of the `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Tuple{Real, Real}`: filter frequency
- `demean::Bool=true`: demean signal prior to analysis
"""
function eeg_cbp!(eeg::NeuroAnalyzer.EEG; pad::Int64=0, frq::Real, demean::Bool=true)

    eeg.eeg_signals = eeg_cbp(eeg, pad=pad, frq=frq, demean=demean).eeg_signals

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_cbp!(EEG, pad=$pad, frq=$frq, demean=$demean)")

    return nothing
end

"""
    eeg_denoise_wien(eeg)

Perform Wiener deconvolution denoising of the `eeg`.

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_denoise_wien(eeg::NeuroAnalyzer.EEG)

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = s_denoise_wien(eeg.eeg_signals)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_denoise_wien(EEG)")

    return eeg_new
end

"""
    eeg_denoise_wien!(eeg)

Perform Wiener deconvolution denoising of the `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_denoise_wien!(eeg::NeuroAnalyzer.EEG)

    eeg.eeg_signals = s_denoise_wien(eeg.eeg_signals)
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_denoise_wien!(EEG)")

    return nothing
end

"""
    eeg_scale(eeg; channel, factor)

Multiply `channel` signal by `factor`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`: channel to invert
- `factor::Real`: channel signal is multiplied by factor

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_scale(eeg::NeuroAnalyzer.EEG; channel::Int64, factor::Real)

    (channel < 1 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channel, :, :] = @views eeg.eeg_signals[channel, :, :] .* factor
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_scale(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_scale!(eeg; channel)

Multiply `channel` signal by `factor`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`: channel to invert
- `factor::Real`: channel signal is multiplied by factor
"""
function eeg_scale!(eeg::NeuroAnalyzer.EEG; channel::Int64, factor::Real)

    (channel < 1 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))

    eeg.eeg_signals[channel, :, :] = @views eeg.eeg_signals[channel, :, :] .* factor
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_scale!(EEG, channel=$channel, factor=$factor)")
end