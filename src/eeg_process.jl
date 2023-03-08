"""
    eeg_derivative(eeg; channel)

Return the derivative of EEG channel(s) with length same as the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_derivative(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    ep_n = eeg_epoch_n(eeg)
    
    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in eachindex(channel)
            @views eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = s_derivative(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx])
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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_derivative!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_derivative(eeg, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_detrend(eeg; channel, type, offset, order, f)

Perform piecewise detrending of EEG channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
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

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_detrend(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)

    _check_var(type, [:ls, :linear, :constant, :poly, :loess, :hp], "type")

    ep_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in eachindex(channel)
            @views eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = s_detrend(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], type=type, offset=offset, order=order, f=f, fs=fs)
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

- `obj::NeuroAnalyzer.NEURO`
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
function eeg_detrend!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)

    eeg_tmp = eeg_detrend(eeg, channel=channel, type=type, offset=offset, order=order, f=f)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_taper(eeg; channel, taper)

Taper EEG channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `taper::Union{Vector{Real, Vector{ComplexF64}}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_taper(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    _check_channels(eeg, channel)

    ep_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `taper::Union{Vector{<:Real}, Vector{ComplexF64}}`
"""
function eeg_taper!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    eeg_tmp = eeg_taper(eeg, channel=channel, taper=taper)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_demean(eeg; channel)

Remove mean value (DC offset).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_demean(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)

    ep_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_demean!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_demean(eeg, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_normalize(eeg; channel, method)

Normalize EEG channel(s)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `method::Symbol`: method for normalization, see `normalize()` for details

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_normalize(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), method::Symbol)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            @views eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = normalize(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], method=method)
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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `method::Symbol`: method for normalization, see `normalize()` for details
"""
function eeg_normalize!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), method::Symbol)

    eeg_tmp = eeg_normalize(eeg, channel=channel, method=method)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_add_noise(eeg; channel)

Add random noise to EEG channel(s)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_add_noise(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)

    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_add_noise!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_add_noise(eeg, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_filter(eeg; <keyword arguments>)

Apply filtering to EEG channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
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
    - `:sg`: Savitzky-Golay
- `ftype::Union{Symbol, Nothing}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Real=-1`: bandwidth for `:iirnotch` and :remez filters
- `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):
    - `:twopass`
    - `:onepass`
    - `:reverse`: one pass, reverse direction
- `order::Int64=8`: filter order (6 dB/octave) for IIR filters, number of taps for `:remez` filter, attenuation (× 4 dB) for `:fir` filter, polynomial order for `:poly` filter, k-value for `:mavg` and `:mmed` (window length = 2 × k + 1)
- `t::Real`: threshold for `:mavg` and `:mmed` filters; threshold = threshold * std(signal) + mean(signal) for `:mavg` or threshold = threshold * std(signal) + median(signal) for `:mmed` filter
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: kernel for the `:conv` filter, window length for `:sg` and `:poly` filters, weighting window for `:mavg` and `:mmed`
- `preview::Bool=false`: plot filter response

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

For `:poly` filter `order` and `window` have to be set experimentally, recommended initial values are: `order=4` and `window=32`.
"""
function eeg_filter(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple{Real, Real}}=0, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, dir::Symbol=:twopass, d::Int64=1, t::Real=0, window::Union{Nothing, AbstractVector, Int64}=nothing, preview::Bool=false)

    _check_channels(eeg, channel)
    _check_var(fprototype, [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :iirnotch, :remez, :mavg, :mmed, :poly, :conv, :sg], "fprototype")
    ep_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)
    n = epoch_len(eeg)

    if preview == true
        _info("When `preview=true`, signal is not being filtered.")
        fprototype === :iirnotch && (ftype = :bs)    
        p = plot_filter_response(fs=fs, n=n, fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, window=window)
        Plots.plot(p)
        return p
    end

    eeg_new = deepcopy(eeg)

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :iirnotch, :remez]
        n = epoch_len(eeg)
        flt = s_filter_create(fprototype=fprototype, ftype=ftype, cutoff=cutoff, n=n, fs=fs, order=order, rp=rp, rs=rs, bw=bw, window=window)
    end

    # initialize progress bar
    progress_bar == true && (pb = Progress(ep_n * length(channel), 1))

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in eachindex(channel)
            if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :iirnotch, :remez]
                eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = @views s_filter_apply(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], flt=flt, dir=dir)
            else
                eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = @views s_filter(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], fprototype=fprototype, order=order, t=t, window=window)
            end
            # update progress bar
            progress_bar == true && next!(pb)
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
    - `:conv`: convolution
    - `:sg`: Savitzky-Golay
- `ftype::Union{Symbol, Nothing}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Real=-1`: bandwidth for `:iirnotch` and `:remez filters
- `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):
    - `:twopass`
    - `:onepass`
    - `:reverse`: one pass, reverse direction
- `order::Int64=8`: filter order (6 dB/octave) for IIR filters, number of taps for `:remez` filter, attenuation (× 4 dB) for `:fir` filter, polynomial order for `:poly` filter, k-value for `:mavg` and `:mmed` (window length = 2 × k + 1)
- `t::Real`: threshold for `:mavg` and `:mmed` filters; threshold = threshold * std(signal) + mean(signal) for `:mavg` or threshold = threshold * std(signal) + median(signal) for `:mmed` filter
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: kernel for the `:conv` filter, window length for `:sg` and `:poly` filters, weighting window for `:mavg` and `:mmed`
- `preview::Bool=false`: plot filter response

# Notes

For `:poly` filter `order` and `window` have to be set experimentally, recommended initial values are: `order=4` and `window=32`.
"""
function eeg_filter!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple{Real, Real}}=0, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, dir::Symbol=:twopass, t::Real=0, window::Union{Nothing, AbstractVector, Int64}=nothing, preview::Bool=false)

    if preview == true
        _info("When `preview=true`, signal is not being filtered.")
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

    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_pca(eeg; channel, n)

Perform principal component analysis (PCA).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `n::Int64`: number of PCs to calculate

# Returns

Named tuple containing:
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
- `pc_m::Vector{Float64}`: PC mean
- `pca::MultivariateStats.PCA{Float64}`: PC model
"""
function eeg_pca(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), n::Int64)

    _check_channels(eeg, channel)

    pc, pc_var, pc_m, pca = @views s_pca(obj.data[channel, :, :], n=n)

    return (pc=pc, pc_var=pc_var, pc_m=pc_m, pca=pca)
end

"""
    eeg_pca_reconstruct(eeg; channel)

Reconstruct EEG signals using embedded PCA components (`:pc`) and model (`:pca`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_pca_reconstruct(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    :pc in obj.header[:components] || throw(ArgumentError("EEG does not contain :pc component. Perform eeg_pca(EEG) first."))
    :pca in obj.header[:components] || throw(ArgumentError("EEG does not contain :pca component. Perform eeg_pca(EEG) first."))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    pc_idx = findfirst(isequal(:pc), obj.header[:components])
    pc_m_idx = findfirst(isequal(:pca), obj.header[:components])
    eeg_new.eeg_signals[channel, :, :] = @views s_pca_reconstruct(eeg_new.eeg_signals[channel, :, :], pc=eeg_new.eeg_components[pc_idx], pca=eeg_new.eeg_components[pc_m_idx])
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_pca_reconstruct(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_pca_reconstruct!(eeg; channel)

Reconstruct EEG signals using embedded PCA components (`:pc`) and model (`:pca`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_pca_reconstruct!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_pca_reconstruct(eeg, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_pca_reconstruct(eeg, pc, pca; channel)

Reconstruct EEG signals using external PCA components (`pc` and `pca`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pca::MultivariateStats.PCA{Float64}`: PC model
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_pca_reconstruct(obj::NeuroAnalyzer.NEURO, pc::Array{Float64, 3}, pca::MultivariateStats.PCA{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

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

- `obj::NeuroAnalyzer.NEURO`
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pca::MultivariateStats.PCA{Float64}`: PC model
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_pca_reconstruct!(obj::NeuroAnalyzer.NEURO, pc::Array{Float64, 3}, pca::MultivariateStats.PCA{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_pca_reconstruct(eeg, pc, pca, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_ica(eeg; <keyword arguments>)

Perform independent component analysis (ICA).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
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
function eeg_ica(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    _check_channels(eeg, channel)

    ica, ica_mw = @views s_ica(obj.data[channel, :, :], n=n, tol=tol, iter=iter, f=f)

    return (ica=ica, ica_mw=ica_mw)
end

"""
    eeg_ica_reconstruct(eeg; channel, ic)

Reconstruct EEG signals using embedded ICA components (`:ica` and `:ica_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `ic::Union{Int64, Vector{Int64}, AbstractRange}`: list of ICs to remove

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_ica_reconstruct(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), ic::Union{Int64, Vector{Int64}, AbstractRange})

    :ica in obj.header[:components] || throw(ArgumentError("EEG does not contain :ica component. Perform eeg_ica(EEG) first."))
    :ica_mw in obj.header[:components] || throw(ArgumentError("EEG does not contain :ica_mw component. Perform eeg_ica(EEG) first."))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    ica_a_idx = findfirst(isequal(:ica), obj.header[:components])
    ica_mw_idx = findfirst(isequal(:ica_mw), obj.header[:components])
    eeg_new.eeg_signals[channel, :, :] = @views s_ica_reconstruct(eeg_new.eeg_signals[channel, :, :], ica=eeg_new.eeg_components[ica_a_idx], ica_mw=eeg_new.eeg_components[ica_mw_idx], ic=ic)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_ica_reconstruct(EEG, channel=$channel, ic=$ic)")

    return eeg_new
end

"""
    eeg_ica_reconstruct!(eeg; channel, ic)

Reconstruct EEG signals using embedded ICA components (`:ica` and `:ica_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove
"""
function eeg_ica_reconstruct!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), ic::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_tmp = eeg_ica_reconstruct(eeg, channel=channel, ic=ic)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_ica_reconstruct(eeg, ica, ica_mw; channel, ic)

Reconstruct EEG signals using external ICA components (`ica` and `ica_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ica::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ica_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_ica_reconstruct(obj::NeuroAnalyzer.NEURO, ica::Array{Float64, 3}, ica_mw::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), ic::Union{Int64, Vector{Int64}, AbstractRange})

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

- `obj::NeuroAnalyzer.NEURO`
- `ica::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ica_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove
"""
function eeg_ica_reconstruct!(obj::NeuroAnalyzer.NEURO, ica::Array{Float64, 3}, ica_mw::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), ic::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_tmp = eeg_ica_reconstruct(eeg, ica, ica_mw, channel=channel, ic=ic)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_average(eeg; channel)

Return the average signal of EEG channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_average(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    eeg_keep_channel!(eeg_new, channel=1)
    eeg_new.eeg_signals = @views s_average(obj.data[channel, :, :])
    eeg_new.eeg_header[:labels]=["averaged channel"]
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_average(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_average!(eeg; channel)

Return the average signal of EEG channels.  

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_average!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_average(eeg, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_average(eeg1, eeg2)

Return the average signal of all `eeg1` and `eeg2` channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_average(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)

    size(eeg1.eeg_signals) == size(eeg2.eeg_signals) || throw(ArgumentError("Both signals must have the same size."))
    ch_n = eeg_channel_n(eeg1)
    ep_n = eeg_epoch_n(eeg1)

    eeg_new = deepcopy(eeg1)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            eeg_new.eeg_signals[channel_idx, :, epoch_idx] = @views s2_average(eeg1.eeg_signals[channel_idx, :, epoch_idx], eeg2.eeg_signals[channel_idx, :, epoch_idx])
        end
    end

    eeg_reset_components!(eeg_new)
    push!(obj.header[:history], "eeg_average(EEG1, EEG2)")

    return eeg_new
end

"""
    eeg_invert_polarity(eeg; channel)

Invert polarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_invert_polarity(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: channel(s) to invert
"""
function eeg_invert_polarity!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_invert_polarity(eeg, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_wdenoise(eeg; channel, wt)

Perform wavelet denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_wdenoise(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), wt::T) where {T <: DiscreteWavelet}

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = @views s_wdenoise(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], wt=wt)
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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
"""
function eeg_wdenoise!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), wt::T) where {T <: DiscreteWavelet}

    eeg_tmp = eeg_wdenoise(eeg, channel=channel, wt=wt)
    eeg_tmp = eeg_upsample(eeg, new_sr=new_sr)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_fftdenoise(eeg; channel, pad, threshold)

Perform FFT denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `threshold::Int64=100`: PSD threshold for keeping frequency components

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_fftdenoise(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, threshold::Int64=100)

    _check_channels(eeg, channel)

    ep_n = eeg_epoch_n(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in eachindex(channel)
                eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], _ = @views s_fftdenoise(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, threshold=threshold)
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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `threshold::Int64=100`: PSD threshold for keeping frequency components
"""
function eeg_fftdenoise!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, threshold::Int64=100)

    eeg_tmp = eeg_fftdenoise(eeg, channel=channel, pad=pad, threshold=threshold)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_zero(eeg)

Zero EEG channels at the beginning and at the end.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_zero(obj::NeuroAnalyzer.NEURO)

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

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_zero!(obj::NeuroAnalyzer.NEURO)

    eeg_tmp = eeg_zero(eeg)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_wbp(eeg; channel, pad, frq, ncyc, demean)

Perform wavelet bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_wbp(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, frq::Real, ncyc::Int64=6, demean::Bool=true)

    _check_channels(eeg, channel)

    ep_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in eachindex(channel)
            eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = @views s_wbp(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, frq=frq, fs=fs, ncyc=ncyc, demean=demean)
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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis
"""
function eeg_wbp!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, frq::Real, ncyc::Int64=6, demean::Bool=true)

    eeg_tmp = eeg_wbp(eeg, channel=channel, pad=pad, frq=frq, ncyc=ncyc, demean=demean)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_cbp(eeg; channel, pad, frq, demean)

Perform convolution bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_cbp(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, frq::Real, demean::Bool=true)

    _check_channels(eeg, channel)

    ep_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in eachindex(channel)
            eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx] = @views s_cbp(eeg_new.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, frq=frq, fs=fs, demean=demean)
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

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Tuple{Real, Real}`: filter frequency
- `demean::Bool=true`: demean signal prior to analysis
"""
function eeg_cbp!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), pad::Int64=0, frq::Real, demean::Bool=true)

    eeg_tmp = eeg_cbp(eeg, channel=channel, pad=pad, frq=frq, demean=demean)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_denoise_wien(eeg; channel)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

# Returns
- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_denoise_wien(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channel, :, :] = s_denoise_wien(eeg_new.eeg_signals[channel, :, :])
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_denoise_wien(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_denoise_wien!(eeg; channel)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
"""
function eeg_denoise_wien!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)))

    eeg_tmp = eeg_denoise_wien(eeg, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_scale(eeg; channel, factor)

Multiply EEG channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `factor::Real`: channel signal is multiplied by `factor`

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_scale(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), factor::Real)

    _check_channels(eeg, channel)
    
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals[channel, :, :] = @views eeg_new.eeg_signals[channel, :, :] .* factor
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_scale(EEG, channel=$channel)")

    return eeg_new
end

"""
    eeg_scale!(eeg; channel, factor)

Multiply EEG channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `factor::Real`: channel signal is multiplied by `factor`
"""
function eeg_scale!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), factor::Real)

    eeg_tmp = eeg_scale(eeg, channel=channel, factor=factor)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_slaplacian(eeg; m, n, s)

Transform signal channels using surface Laplacian.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `m::Int64=8`: constant positive integer for smoothness
- `n::Int64=8`: Legendre polynomial order
- `s::Float64=10^-5`: smoothing factor

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
- `G::Matrix{Float64}`: transformation matrix
- `H::Matrix{Float64}`: transformation matrix

# Source

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-7
"""
function eeg_slaplacian(obj::NeuroAnalyzer.NEURO; m::Int64=4, n::Int64=8, s::Float64=10^-5)

    obj.header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    m < 1 && throw(ArgumentError("m must be ≥ 1."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))
    s <= 0 && throw(ArgumentError("s must be > 0."))

    channels = eeg_signal_channels(eeg)
    locs = obj.locs
    ch_n = nrow(locs)
    ep_n = eeg_epoch_n(eeg)

    length(channels) > nrow(locs) && throw(ArgumentError("Some channels do not have locations."))

    G = zeros(ch_n, ch_n)
    H = zeros(ch_n, ch_n)
    cosdist = zeros(ch_n, ch_n)

    # r = locs[!, :loc_radius_sph]
    # x = locs[!, :loc_x] ./ maximum(r)
    # y = locs[!, :loc_y] ./ maximum(r)
    # z = locs[!, :loc_z] ./ maximum(r)

    x = locs[!, :loc_x]
    y = locs[!, :loc_y]
    z = locs[!, :loc_z]
    x, y, z = _locnorm(x, y, z)

    # cosine distance
    Threads.@threads for i = 1:ch_n
        @inbounds @simd  for j = 1:ch_n
            cosdist[i, j]  =  1 - (( (x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2 ) / 2 )
        end
    end

    # compute Legendre polynomial
    legpoly = zeros(n, ch_n, ch_n)
    Threads.@threads for idx1 in 1:n
        @inbounds @simd for idx2 in 1:ch_n
            legpoly[idx1, idx2, :] = @views legendre.(cosdist[idx2, :], idx1)
        end
    end

    # compute G and H
    Threads.@threads for i = 1:ch_n
        @inbounds @simd for j = 1:ch_n
            g = 0
            h = 0
            @inbounds @simd for idx_n = 1:n
                g = g + ((2 * idx_n + 1) * legpoly[idx_n, i, j] / ((idx_n * (idx_n + 1))^m))
                h = h - ((2 * (idx_n + 1) * legpoly[idx_n, i, j]) / ((idx_n * (idx_n + 1))^(m - 1)))
            end
            G[i, j] =  g / (4 * pi)
            H[i, j] = -h / (4 * pi)
        end
    end

    # add smoothing factor to the diagonal
    Gs = G + I(ch_n) * s
    GsinvS = sum(inv(Gs))

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
        data = @views obj.data[channels, :, epoch_idx]
        # dataGs = data[channels, :]' / Gs
        dataGs = Gs / data'
        # C = dataGs .- (sum(dataGs,dims=2)/sum(GsinvS))*GsinvS
        C = data .- (sum(dataGs, dims=2) / sum(GsinvS)) * GsinvS
        # compute surface Laplacian
        eeg_new.eeg_signals[channels, :, epoch_idx] = (C'*H)'
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_surface_laplacian(EEG, m=m, n=n, s=s)")

    return eeg_new, G, H
end

"""
    eeg_slaplacian!(eeg; m, n, s)

Transform signal channels using surface Laplacian.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `m::Int64=8`: constant positive integer for smoothness
- `n::Int64=8`: Legendre polynomial order
- `s::Float64=10^-5`: smoothing factor

# Returns

- `G::Matrix{Float64}`: transformation matrix
- `H::Matrix{Float64}`: transformation matrix

# Source

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-7
"""
function eeg_slaplacian!(obj::NeuroAnalyzer.NEURO; m::Int64=4, n::Int64=8, s::Float64=10^-5)

    eeg_tmp, G, H = eeg_slaplacian(eeg, m=m, n=n, s=s)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return G, H
end