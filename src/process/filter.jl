export filter_apply
export filter_apply!
export filter_create
export filter
export filter!

"""
    filter_create(; <keyword arguments>)

Create FIR or IIR filter.

# Arguments

- `fprototype::Symbol`: filter prototype:
    - `:fir`: FIR filter
    - `:firls`: weighted least-squares FIR filter
    - `:remez`: Remez FIR filter
    - `:butterworth`: IIR filter
    - `:chebyshev1` IIR filter
    - `:chebyshev2` IIR filter
    - `:elliptic` IIR filter
    - `:iirnotch`: second-order IIR notch filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (must be a pair of frequencies for `:bp` and `:bs`)
- `n::Int64`: signal length in samples
- `fs::Int64`: sampling rate
- `order::Union{Nothing, Int64}=nothing`: filter order (6 dB/octave) for `:butterworth`, `:chebyshev1`, `:chebyshev2`, `:elliptic` or number of taps for `:fir`, :firls` and `:remez` filters
- `rp::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Union{Nothing, Real}=nothing`: transition band width for `:iirnotch`, `:firls` and `:remez` filters
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter

# Returns

- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`
"""
function filter_create(; fprototype::Symbol, ftype::Union{Nothing, Symbol}=nothing, cutoff::Union{Real, Tuple{Real, Real}}, n::Int64, fs::Int64, order::Union{Nothing, Int64}=nothing, rp::Union{Nothing, Real}=nothing, rs::Union{Nothing, Real}=nothing, bw::Union{Nothing, Real}=nothing, w::Union{Nothing, AbstractVector}=nothing)::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}

    _check_var(fprototype, [:fir, :firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch], "fprototype")

    if fprototype !== :iirnotch
        @assert ftype !== nothing "ftype must be specified."
        _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")
    end

    @assert fs >= 1 "fs must be ≥ 1."

    fprototype in [:firls, :remez, :iirnotch] && @assert !isa(bw, Nothing) "bw must be specified."

    if !isa(bw, Nothing)
        @assert bw > 0 "bw must be > 0."
        @assert bw < cutoff[1] "bw must be < $(cutoff[1])."
        length(cutoff) == 2 && (@assert bw < fs - cutoff[2] "bw must be < $(fs - cutoff[2]).")
    end

    !isa(order, Nothing) && @assert order <= n "order must be ≤ signal length ($n)."

    !isa(w, Nothing) && @assert length(w) <= n "Length of w ($(length(w))) must be ≤ signal length ($n)."

    if fprototype in [:fir, :butterworth, :chebyshev1, :chebyshev2, :elliptic]
        if ftype === :lp
            @assert length(cutoff) == 1 "For :lp filter one frequency must be given."
            responsetype = Lowpass(cutoff)
        elseif ftype === :hp
            @assert length(cutoff) == 1 "For :hp filter one frequency must be given."
            responsetype = Highpass(cutoff)
        elseif ftype === :bp
            @assert length(cutoff) == 2 "For :bp filter two frequencies must be given."
            responsetype = Bandpass(cutoff[1], cutoff[2])
        elseif ftype === :bs
            @assert length(cutoff) == 2 "For :bs filter two frequencies must be given."
            responsetype = Bandstop(cutoff[1], cutoff[2])
        end
    end

    ## FIR filters
    if fprototype in [:fir, :firls, :remez]
        if fprototype === :fir
            w === nothing && (w = DSP.hamming(order))
            prototype = FIRWindow(w)

            ftype in [:hp, :bp, :bs] && @assert mod(order, 2) != 0 "order must be odd."
            if ftype in [:lp, :hp]
                ftype === :lp && _info("Creating LP filter:")
                ftype === :hp && _info("Creating HP filter:")
                _info(" Number of taps: $(length(w))")
            elseif ftype === :bp
                _info("Creating BP filter:")
                _info(" Number of taps: $(length(w))")
            elseif ftype === :bs
                _info("Creating BS filter:")
                _info(" Number of taps: $(length(w))")
            end

            flt = digitalfilter(responsetype, prototype, fs=fs)

            return flt

        elseif fprototype === :firls
            if ftype === :bp
                f1_stop = cutoff[1] - bw
                f1_pass = cutoff[1] + bw
                f2_pass = cutoff[2] - bw
                f2_stop = cutoff[2] + bw
                flt_shape = [0, 0, 1, 1, 0, 0]
                flt_frq = [0, f1_stop, f1_pass, f2_pass, f2_stop, fs / 2]
                if !isa(w, Nothing)
                    @assert length(w) == 6 "Length of w must be 6."
                else
                    w = ones(6)
                end
            elseif ftype === :bs
                f1_pass = cutoff[1] - bw
                f1_stop = cutoff[1] + bw
                f2_stop = cutoff[2] - bw
                f2_pass = cutoff[2] + bw
                flt_shape = [1, 1, 0, 0, 1, 1]
                flt_frq = [0, f1_pass, f1_stop, f2_stop, f2_pass, fs / 2]
                if !isa(w, Nothing)
                    @assert length(w) == 6 "Length of w must be 6."
                else
                    w = ones(6)
                end
            elseif ftype === :lp
                f_pass = cutoff[1] - bw
                f_stop = cutoff[1] + bw
                flt_shape = [1, 1, 0, 0]
                flt_frq = [0, f_pass, f_stop, fs / 2]
                if !isa(w, Nothing)
                    @assert length(w) == 4 "Length of w must be 4."
                else
                    w = ones(4)
                end
            elseif ftype === :hp
                f_pass = cutoff[1] + bw
                f_stop = cutoff[1] - bw
                flt_shape = [0, 0, 1, 1]
                flt_frq = [0, f_stop, f_pass, fs / 2]
                if !isa(w, Nothing)
                    @assert length(w) == 4 "Length of w must be 4."
                else
                    w = ones(4)
                end
            end

            if ftype in [:lp, :hp]
                ftype === :lp && _info("Creating LP filter:")
                ftype === :hp && _info("Creating HP filter:")
                _info(" Number of taps: $(order + 1)")
                _info(" Transition band width: $bw Hz")
                _info(" F_pass: $f_pass Hz")
                _info(" F_stop: $f_stop Hz")
            elseif ftype === :bp
                _info("Creating BP filter:")
                _info(" Number of taps: $(order + 1)")
                _info(" Transition band width: $bw Hz")
                _info(" F1_stop: $f1_stop Hz")
                _info(" F1_pass: $f1_pass Hz")
                _info(" F2_pass: $f2_pass Hz")
                _info(" F2_stop: $f2_stop Hz")
            elseif ftype === :bs
                _info("Creating BS filter:")
                _info(" Number of taps: $(order + 1)")
                _info(" Transition band width: $bw Hz")
                _info(" F1_pass: $f1_pass Hz")
                _info(" F1_stop: $f1_stop Hz")
                _info(" F2_stop: $f2_stop Hz")
                _info(" F2_pass: $f2_pass Hz")
            end

            flt = firls_design(order, flt_frq, flt_shape, w, true, fs=fs)

            return flt

        elseif fprototype === :remez
            if ftype === :bp
                f1_stop = cutoff[1] - bw
                f1_pass = cutoff[1] + bw
                f2_pass = cutoff[2] - bw
                f2_stop = cutoff[2] + bw
                w = [(0, f1_stop) => 0, (f1_pass, f2_pass) => 1, (f2_stop, fs / 2) => 0]
            elseif ftype === :bs
                f1_pass = cutoff[1] - bw
                f1_stop = cutoff[1] + bw
                f2_stop = cutoff[2] - bw
                f2_pass = cutoff[2] + bw
                w = [(0, f1_pass) => 1, (f1_stop, f2_stop) => 0, (f2_pass, fs / 2) => 1]
            elseif ftype === :lp
                f_pass = cutoff[1] - bw
                f_stop = cutoff[1] + bw
                w = [(0, f_pass) => 1, (f_stop, fs / 2) => 0]
            elseif ftype === :hp
                f_pass = cutoff[1] + bw
                f_stop = cutoff[1] - bw
                w = [(0, f_stop) => 0, (f_pass, fs / 2) => 1]
            end

            if ftype in [:lp, :hp]
                ftype === :lp && _info("Creating LP filter:")
                ftype === :hp && _info("Creating HP filter:")
                _info(" Number of taps: $order")
                _info(" Transition band width: $bw Hz")
                _info(" F_pass: $f_pass Hz")
                _info(" F_stop: $f_stop Hz")
            elseif ftype === :bp
                _info("Creating BP filter:")
                _info(" Number of taps: $order")
                _info(" Transition band width: $bw Hz")
                _info(" F1_stop: $f1_stop Hz")
                _info(" F1_pass: $f1_pass Hz")
                _info(" F2_pass: $f2_pass Hz")
                _info(" F2_stop: $f2_stop Hz")
            elseif ftype === :bs
                _info("Creating BS filter:")
                _info(" Number of taps: $order")
                _info(" Transition band width: $bw Hz")
                _info(" F1_pass: $f1_pass Hz")
                _info(" F1_stop: $f1_stop Hz")
                _info(" F2_stop: $f2_stop Hz")
                _info(" F2_pass: $f2_pass Hz")
            end

            flt = remez(order, w, Hz=fs)

            return flt

        end
    end

    ## IIR filters

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic]

        if rp isa Nothing
            if fprototype === :elliptic
                rp = 0.0025
            else
                rp = 2
            end
        end

        if rs isa Nothing
            if fprototype === :elliptic
                rp = 40
            else
                rp = 20
            end
        end

        if fprototype === :butterworth
            prototype = Butterworth(order)
        elseif fprototype === :chebyshev1
            _in(rs, (0, fs / 2), "rs")
            prototype = Chebyshev1(order, rs)
        elseif fprototype === :chebyshev2
            _in(rs, (0, fs / 2), "rs")
            prototype = Chebyshev2(order, rp)
        elseif fprototype === :elliptic
            _in(rs, (0, fs / 2), "rs")
            _in(rp, (0, fs / 2), "rs")
            prototype = Elliptic(order, rp, rs)
        end

        flt = digitalfilter(responsetype, prototype, fs=fs)

        return flt

    end

    if fprototype === :iirnotch
        if !isa(ftype, Nothing)
            _info("For :iirnotch filter ftype is ignored")
        end
        @assert length(cutoff) == 1 "For :iirnotch filter cutoff must contain only one frequency."

        flt = iirnotch(cutoff[1], bw, fs=fs)

        return flt

    end

end

"""
    filter_apply(s; <keyword arguments>)

Apply IIR or FIR filter.

# Arguments

- `s::AbstractVector`
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter
- `dir:Symbol=:twopass`: filter direction:
    - `:twopass`: two passes, the resulting signal has zero phase distortion, the effective filter order is doubled
    - `:onepass`: single pass, use for for causal filter
    - `:reverse`: one pass, reverse direction

# Returns

- `s_new::Vector{Float64}`
"""
function filter_apply(s::AbstractVector; flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}, dir::Symbol=:twopass)::Vector{Float64}

    _check_var(dir, [:twopass, :onepass, :reverse], "dir")

    dir === :onepass && (return filt(flt, s))
    dir === :twopass && (return filtfilt(flt, s))
    dir === :reverse && (return filt(flt, reverse(s)))

end

"""
    filter_apply(obj; <keyword arguments>)

Apply IIR or FIR filter.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter
- `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):
    - `:twopass`
    - `:onepass`
    - `:reverse`: one pass, reverse direction

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function filter_apply(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}, dir::Symbol=:twopass)::NeuroAnalyzer.NEURO

    _check_var(dir, [:twopass, :onepass, :reverse], "dir")

    ch = get_channel(obj, ch=ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    ep_n > 1 && _warn("filter() should be applied to a continuous signal.")
    _info("Signal should be tapered prior to filtering to reduce edge artifacts")

    obj_new = deepcopy(obj)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * length(ch), dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in eachindex(ch)
            obj_new.data[ch[ch_idx], :, ep_idx] = @views filter_apply(obj.data[ch[ch_idx], :, ep_idx], flt=flt, dir=dir)
            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.history, "filter(OBJ, ch=$ch, dir=$dir)")

    return obj_new

end

"""
    filter_apply!(obj; <keyword arguments>)

Apply IIR or FIR filter.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter
- `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):
    - `:twopass`
    - `:onepass`
    - `:reverse`: one pass, reverse direction

# Returns

Nothing
"""
function filter_apply!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}, dir::Symbol=:twopass)::Nothing

    _check_var(dir, [:twopass, :onepass, :reverse], "dir")

    obj_new = NeuroAnalyzer.filter_apply(obj,
                                         ch=ch,
                                         flt=flt,
                                         dir=dir)

    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end

"""
    filter(obj; <keyword arguments>)

Apply filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`: IIR filter
    - `:chebyshev1` IIR filter
    - `:chebyshev2` IIR filter
    - `:elliptic` IIR filter
    - `:iirnotch`: second-order IIR notch filter
    - `:fir`: FIR filter
    - `:firls`: weighted least-squares FIR filter
    - `:remez`: Remez FIR filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (must be a pair of frequencies for `:bp` and `:bs`)
- `rp::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Union{Nothing, Real}=nothing`: bandwidth for `:iirnotch` and `:remez` filters
- `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):
    - `:twopass`
    - `:onepass`
    - `:reverse`: one pass, reverse direction
- `order::Union{Nothing, Int64}=nothing`: default is 4, filter order (6 dB/octave) for `:butterworth`, `:chebyshev1`, `:chebyshev2`, `:elliptic` (for these filters is multiplied by two since must be even), number of taps for `:remez`, attenuation (× 4 dB) for `:fir` filters
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window, number of taps is calculated using Fred Harris' rule-of-thumb) or weights for `:firls` filter
- `preview::Bool=false`: plot filter response

# Returns

- `obj_new::NeuroAnalyzer.NEURO`

If `preview=true`, it will return `Plots.Plot{Plots.GRBackend}`.
"""
function filter(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, fprototype::Symbol, ftype::Union{Nothing, Symbol}=nothing, cutoff::Union{Real, Tuple{Real, Real}}, order::Union{Nothing, Int64}=nothing, rp::Union{Nothing, Real}=nothing, rs::Union{Nothing, Real}=nothing, bw::Union{Nothing, Real}=nothing, dir::Symbol=:twopass, w::Union{Nothing, AbstractVector}=nothing, preview::Bool=false)::Union{NeuroAnalyzer.NEURO, Plots.Plot{Plots.GRBackend}}

    _check_var(fprototype, [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :firls, :iirnotch, :remez], "fprototype")

    ch = get_channel(obj, ch=ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    ep_n > 1 && _warn("filter() should be applied to a continuous signal.")
    _info("Signal should be tapered prior to filtering to reduce edge artifacts")

    if preview
        _info("Previewing filter response, signal will not be filtered")
        if fprototype !== :firls
            order === nothing && (order = 4)
        end
        fprototype === :iirnotch && (ftype = :bs)
        p = plot_filter_response(fs=fs, n=epoch_len(obj), fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, w=w)
        Plots.plot(p)
        return p
    end

    obj_new = deepcopy(obj)

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :firls, :iirnotch, :remez]
        ftype in [:bs, :bp] && @assert length(cutoff) == 2 "For :bp and :bs filters cutoff must contain two frequencies."
        flt = filter_create(fprototype=fprototype, ftype=ftype, cutoff=cutoff, n=epoch_len(obj), fs=fs, order=order, rp=rp, rs=rs, bw=bw, w=w)
    end

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * length(ch), dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in eachindex(ch)
            obj_new.data[ch[ch_idx], :, ep_idx] = @views filter_apply(obj.data[ch[ch_idx], :, ep_idx], flt=flt, dir=dir)
            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.history, "filter(OBJ, ch=$ch, fprototype=$fprototype, ftype=$ftype, cutoff=$cutoff, order=$order, rp=$rp, rs=$rs, dir=$dir, w=$w)")

    return obj_new

end

"""
    filter!(obj; <keyword arguments>)

Apply filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`: IIR filter
    - `:chebyshev1` IIR filter
    - `:chebyshev2` IIR filter
    - `:elliptic` IIR filter
    - `:iirnotch`: second-order IIR notch filter
    - `:fir`: FIR filter
    - `:firls`: weighted least-squares FIR filter
    - `:remez`: Remez FIR filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (must be a pair of frequencies for `:bp` and `:bs`)
- `rp::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Union{Nothing, Real}=nothing`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Union{Nothing, Real}=nothing`: bandwidth for `:iirnotch` and `:remez` filters
- `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):
    - `:twopass`
    - `:onepass`
    - `:reverse`: one pass, reverse direction
- `order::Union{Nothing, Int64}=nothing`: default is 4, filter order (6 dB/octave) for `:butterworth`, `:chebyshev1`, `:chebyshev2`, `:elliptic` (for these filters is multiplied by two since must be even), number of taps for `:remez`, attenuation (× 4 dB) for `:fir` filters
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window, number of taps is calculated using Fred Harris' rule-of-thumb) or weights for `:firls` filter
- `preview::Bool=false`: plot filter response

# Returns

Nothing

If `preview=true`, it will return `Plots.Plot{Plots.GRBackend}`.
"""
function filter!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple{Real, Real}}, order::Union{Nothing, Int64}=nothing, rp::Union{Nothing, Real}=nothing, rs::Union{Nothing, Real}=nothing, bw::Union{Nothing, Real}=nothing, dir::Symbol=:twopass, t::Real=0, w::Union{Nothing, AbstractVector}=nothing, preview::Bool=false)::Union{Nothing, Plots.Plot{Plots.GRBackend}}

    if preview
        _info("Previewing filter response, signal will not be filtered")
        if fprototype !== :firls
            order === nothing && (order = 4)
        end
        fprototype === :iirnotch && (ftype = :bs)
        p = plot_filter_response(fs=sr(obj), fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, w=w)
        Plots.plot(p)
        return p
    end

    dir === :twopass && _info("Filter is applied twice, the effective filter order is doubled")

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
                                   w=w)

    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
