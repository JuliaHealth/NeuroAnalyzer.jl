export filter_apply
export filter_create
export filter
export filter!

"""
    filter_create(signal; <keyword arguments>)

Create IIR or FIR filter.

# Arguments

- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`: second-order IIR notch filter
    - `:remez`: Remez FIR filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}=0`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
- `n::Int64`: signal length in samples
- `fs::Int64`: sampling rate
- `order::Int64=8`: filter order (6 dB/octave), number of taps for `:remez`, attenuation (× 4 dB) for `:fir` filters
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Real=-1`: bandwidth for `:iirnotch` and :remez filters
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter; default is Hamming window, number of taps is calculated using Fred Harris' rule-of-thumb

# Returns

- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`
"""
function filter_create(;fprototype::Symbol, ftype::Union{Nothing, Symbol}=nothing, cutoff::Union{Real, Tuple{Real, Real}}=0, n::Int64, fs::Int64, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, window::Union{Nothing, AbstractVector, Int64}=nothing)

    _check_var(fprototype, [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :iirnotch, :remez], "fprototype")
    if fprototype !== :iirnotch
        ftype === nothing && throw(ArgumentError("ftype must be specified."))
        _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")
    end

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    fprototype !== :iirnotch && order < 1 && throw(ArgumentError("order must be > 1."))
    order > n && throw(ArgumentError("order must be ≤ signal length ($n)."))
    ((order < 2 && fprototype !== :iirnotch && fprototype !== :remez && fprototype !== :fir) && mod(order, 2) != 0) && throw(ArgumentError("order must be even and ≥ 2."))
    window !== nothing && length(window) > n && throw(ArgumentError("window must be ≤ signal length ($n)."))

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch, :remez]
        cutoff == 0 && throw(ArgumentError("cutoff must be specified."))
    end

    if fprototype in [:iirnotch, :remez]
        bw == -1 && throw(ArgumentError("bw must be specified."))
    end

    if fprototype === :iirnotch
        if ftype !== nothing
            _info("For :iirnotch filter ftype is ignored.")
            ftype = nothing
        end
        length(cutoff) == 2 && throw(ArgumentError("For :iirnotch filter cutoff must contain only one frequency."))
    end

    if fprototype === :remez
        bw > cutoff[1] && throw(ArgumentError("For :remez filter bw must be ≤ $(cutoff[1])."))
        (length(cutoff) == 2 && bw > cutoff[2] - cutoff[1]) && throw(ArgumentError("For :remez filter bw must be ≤ $(cutoff[2] - cutoff[1])."))
    end

    if fprototype === :fir
        if window === nothing
            if ftype === :bp
                f1_stop = cutoff[1] - ((cutoff[2] - cutoff[1]) / fs) / 2
                f1_pass = cutoff[1]
                f2_stop = cutoff[2] + ((cutoff[2] - cutoff[1]) / fs) / 2
                f2_pass = cutoff[2]
                trans_bandwidth = (f2_stop - f1_stop) / fs
                n_taps = round(Int64, order * 4 / (22 * trans_bandwidth))   
            elseif ftype === :bs
                # TO DO: CHECK
                f1_pass = cutoff[1]
                f1_stop = cutoff[1] + ((cutoff[2] - cutoff[1]) / fs) / 2
                f2_stop = cutoff[2] - ((cutoff[2] - cutoff[1]) / fs) / 2
                f2_pass = cutoff[2]
                trans_bandwidth = (f2_stop - f1_stop) / fs
                n_taps = round(Int64, order * 4 / (22 * trans_bandwidth))   
            elseif ftype === :lp
                f_pass = cutoff[1]
                f_stop = cutoff[1] + minimum([maximum([0.25 * cutoff[1], 2.0]), cutoff[1]])
                trans_bandwidth = (f_stop - f_pass) / fs
                n_taps = ceil(Int64, ((order * 4) / (22 * trans_bandwidth)))
            elseif ftype === :hp
                f_pass = cutoff[1]
                f_stop = cutoff[1] - minimum([maximum([0.25 * cutoff[1], 2.0]), cutoff[1]])
                trans_bandwidth = (f_pass - f_stop) / fs
                n_taps = ceil(Int64, ((order * 4) / (22 * trans_bandwidth)))
            end

            # next power of 2
            n_taps = 2 ^ ceil(Int64, log2(n_taps))

            # filter cannot be longer than signal
            if n_taps > n
                n_taps = n
                _info("Reducing window length to $n_taps taps")
            end

            window = DSP.hamming(n_taps)

            if ftype === :hp || ftype === :bp || ftype === :bs
                mod(length(window), 2) == 0 && (window = vcat(window[1:((length(window) ÷ 2) - 1)], window[((length(window) ÷ 2) + 1):end]))
            end

            if ftype === :lp || ftype === :hp
                ftype === :lp && _info("Creating LP filter:")
                ftype === :hp && _info("Creating HP filter:")
                _info(" Using default window: hamming($n_taps)")
                _info(" Attenuation: $(order * 4) dB")
                _info(" F_pass: $(round(f_pass, digits=4)) Hz")
                _info(" F_stop: $(round(f_stop, digits=4)) Hz")
                _info(" Transition bandwidth: $(round(trans_bandwidth, digits=4)) Hz")
                _info(" Cutoff frequency: $(round((cutoff[1] - trans_bandwidth / 2), digits=4)) Hz")
            elseif ftype === :bp
                _info("Creating BP filter:")
                _info(" Using default window: hamming($n_taps)")
                _info(" Attenuation: $(order * 4) dB")
                _info(" F1_stop: $(round(f1_stop, digits=4)) Hz")
                _info(" F1_pass: $f1_pass Hz")
                _info(" F2_pass: $f2_pass Hz")
                _info(" F2_stop: $(round(f2_stop, digits=4)) Hz")
                _info(" Transition bandwidth: $(round(trans_bandwidth, digits=4)) Hz")
                _info(" Cutoff frequency: $(round((cutoff[1] - trans_bandwidth / 2), digits=4)) Hz")
            elseif ftype === :bs
                _info("Creating BS filter:")
                _info(" Using default window: hamming($n_taps)")
                _info(" Attenuation: $(order * 4) dB")
                _info(" F1_pass: $f1_pass Hz")
                _info(" F1_stop: $(round(f1_stop, digits=4)) Hz")
                _info(" F2_stop: $(round(f2_stop, digits=4)) Hz")
                _info(" F2_pass: $f2_pass Hz")
                _info(" Transition bandwidth: $(round(trans_bandwidth, digits=4)) Hz")
                _info(" Cutoff frequency: $(round((cutoff[1] - trans_bandwidth / 2), digits=4)) Hz")
            end
        else
            ftype in [:bp, :bs] && length(window) % 2 == 0 && throw(ArgumentError("For :bp and :bs filters window length must be odd."))
        end
    end

    if rp == -1
        if fprototype === :elliptic
            rp = 0.0025
        else
            rp = 2
        end
    end
    
    if rs == -1
        if fprototype === :elliptic
            rp = 40
        else
            rp = 20
        end
    end

    if ftype === :lp
        length(cutoff) != 1 && throw(ArgumentError("For :lp filter one frequency must be given."))
        responsetype = Lowpass(cutoff; fs=fs)
    elseif ftype === :hp
        length(cutoff) != 1 && throw(ArgumentError("For :hp filter one frequency must be given."))
        responsetype = Highpass(cutoff; fs=fs)
    elseif ftype === :bp
        length(cutoff) != 2 && throw(ArgumentError("For :bp filter two frequencies must be given."))
        responsetype = Bandpass(cutoff[1], cutoff[2]; fs=fs)
    elseif ftype === :bs
        length(cutoff) != 2 && throw(ArgumentError("For :bs filter two frequencies must be given."))
        responsetype = Bandstop(cutoff[1], cutoff[2]; fs=fs)
    end

    if fprototype === :butterworth
        prototype = Butterworth(order)
    elseif fprototype === :chebyshev1
        (rs < 0 || rs > fs / 2) && throw(ArgumentError("For :chebyshev1 filter rs must be ≥ 0 and ≤ $(fs / 2)."))
        prototype = Chebyshev1(order, rs)
    elseif fprototype === :chebyshev2
        (rp < 0 || rp > fs / 2) && throw(ArgumentError("For :chebyshev2 filter rp must be ≥ 0 and ≤ $(fs / 2)."))
        prototype = Chebyshev2(order, rp)
    elseif fprototype === :elliptic
        (rs < 0 || rs > fs / 2) && throw(ArgumentError("For :elliptic filter rs must be ≥ 0 and ≤ $(fs / 2)."))
        (rp < 0 || rp > fs / 2) && throw(ArgumentError("For :elliptic filter rp must be ≥ 0 and ≤ $(fs / 2)."))
        prototype = Elliptic(order, rp, rs)
    elseif fprototype === :fir
        prototype = FIRWindow(window)
    end

    if fprototype === :iirnotch
        flt = iirnotch(cutoff, bw, fs=fs)
    elseif fprototype === :remez
        ftype === :lp && (window = [(0, cutoff - bw) => 1, (cutoff + bw, fs / 2) => 0])
        ftype === :hp && (window = [(0, cutoff - bw) => 0, (cutoff + bw, fs / 2) => 1])
        ftype === :bp && (window = [(0, cutoff[1] - bw / 2) => 0, (cutoff[1] + bw / 2, cutoff[2] - bw / 2) => 1, (cutoff[2] + bw / 2, fs / 2) => 0])
        ftype === :bs && (window = [(0, cutoff[1] - bw / 2) => 1, (cutoff[1] + bw / 2, cutoff[2] - bw / 2) => 0, (cutoff[2] + bw / 2, fs / 2) => 1])
        flt = remez(order, window, Hz=fs)
    else
        flt = digitalfilter(responsetype, prototype)
    end

    return flt

end

"""
    filter_apply(s; <keyword arguments>)

Apply IIR or FIR filter.

# Arguments

- `s::AbstractVector`
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter
- `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):
    - `:twopass`
    - `:onepass`
    - `:reverse`: one pass, reverse direction

# Returns

- `s_filtered::Vector{Float64}`
"""
function filter_apply(s::AbstractVector; flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}, dir::Symbol=:twopass)

    _check_var(dir, [:twopass, :onepass, :reverse], "dir")

    dir === :twopass && (return filtfilt(flt, s))
    dir === :onepass && (return filt(flt, s))
    dir === :reverse && (return filt(flt, reverse(s)))

end

"""
    filter(obj; <keyword arguments>)

Apply filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`: second-order IIR notch filter
    - `:remez`: Remez FIR filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}=0`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Real=-1`: bandwidth for `:iirnotch` and `:remez` filters
- `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):
    - `:twopass`
    - `:onepass`
    - `:reverse`: one pass, reverse direction
- `order::Int64=8`: filter order (6 dB/octave) for IIR filters, number of taps for `:remez` filter, attenuation (× 4 dB) for `:fir` filter
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: window length for `:remez` and `:fir` filters
- `preview::Bool=false`: plot filter response

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

For `:poly` filter `order` and `window` have to be set experimentally, recommended initial values are: `order=4` and `window=32`.
"""
function filter(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), fprototype::Symbol, ftype::Union{Nothing, Symbol}=nothing, cutoff::Union{Real, Tuple{Real, Real}}=0, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, dir::Symbol=:twopass, window::Union{Nothing, AbstractVector, Int64}=nothing, preview::Bool=false)

    _check_channels(obj, ch)
    _check_var(fprototype, [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :iirnotch, :remez], "fprototype")

    ep_n = epoch_n(obj)
    fs = sr(obj)

    if preview == true
        _info("When `preview=true`, signal is not being filtered.")
        fprototype === :iirnotch && (ftype = :bs)    
        p = plot_filter_response(fs=sr(obj), n=epoch_len(obj), fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, window=window)
        Plots.plot(p)
        return p
    end

    obj_new = deepcopy(obj)

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :iirnotch, :remez]
        flt = filter_create(fprototype=fprototype, ftype=ftype, cutoff=cutoff, n=epoch_len(obj), fs=sr(obj), order=order, rp=rp, rs=rs, bw=bw, window=window)
    end

    # initialize progress bar
    progress_bar == true && (pb = Progress(ep_n * length(ch), 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:length(ch)
            obj_new.data[ch[ch_idx], :, ep_idx] = @views filter_apply(obj.data[ch[ch_idx], :, ep_idx], flt=flt, dir=dir)
            # update progress bar
            progress_bar == true && next!(pb)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.history, "filter(OBJ, ch=$ch, fprototype=$fprototype, ftype=$ftype, cutoff=$cutoff, order=$order, rp=$rp, rs=$rs, dir=$dir, window=$window)")

    return obj_new

end

"""
    filter!(obj; <keyword arguments>)

Apply filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`: second-order IIR notch filter
    - `:remez`: Remez FIR filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}=0`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Real=-1`: bandwidth for `:iirnotch` and `:remez` filters
- `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):
    - `:twopass`
    - `:onepass`
    - `:reverse`: one pass, reverse direction
- `order::Int64=8`: filter order (6 dB/octave) for IIR filters, number of taps for `:remez` filter, attenuation (× 4 dB) for `:fir` filter
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: window length for `:remez` and `:fir` filters
- `preview::Bool=false`: plot filter response

# Notes

For `:poly` filter `order` and `window` have to be set experimentally, recommended initial values are: `order=4` and `window=32`.
"""
function filter!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple{Real, Real}}=0, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, dir::Symbol=:twopass, t::Real=0, window::Union{Nothing, AbstractVector, Int64}=nothing, preview::Bool=false)

    if preview == true
        _info("When `preview=true`, signal is not being filtered.")
        fprototype === :iirnotch && (ftype = :bs)
        p = plot_filter_response(fs=sr(obj), fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, window=window)
        Plots.plot(p)
        return p
    end

    obj_new = NeuroAnalyzer.filter(obj,
                                   ch=ch,
                                   fprototype=fprototype,
                                   ftype=ftype,
                                   cutoff=cutoff,
                                   order=order,
                                   rp=rp,
                                   rs=rs,
                                   bw=bw,
                                   dir=dir,
                                   window=window)

    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
