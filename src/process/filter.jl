export filter_create
export filter_apply
export filter_apply!
export filter
export filter!

"""
    filter_create(; <keyword arguments>)

Create FIR or IIR filter.

# Arguments

  - `fprototype::Symbol`: filter prototype:
      + `:fir`: FIR filter
      + `:firls`: weighted least-squares FIR filter
      + `:remez`: Remez FIR filter
      + `:butterworth`: IIR filter
      + `:chebyshev1` IIR filter
      + `:chebyshev2` IIR filter
      + `:elliptic` IIR filter
      + `:iirnotch`: second-order IIR notch filter
  - `ftype::Union{Nothing, Symbol}=nothing`: filter type:
      + `:lp`: low pass
      + `:hp`: high pass
      + `:bp`: band pass
      + `:bs`: band stop
  - `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (must be a pair of frequencies for `:bp` and `:bs`)
  - `fs::Int64`: signal sampling rate
  - `order::Union{Nothing, Int64}=nothing`: filter order
  - `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.5 dB
  - `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 20 dB
  - `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
  - `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter

# Returns

  - `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`
"""
function filter_create(;
        fprototype::Symbol,
        ftype::Union{Nothing, Symbol} = nothing,
        cutoff::Union{Real, Tuple{Real, Real}},
        fs::Int64,
        order::Union{Nothing, Int64} = nothing,
        rp::Union{Nothing, Real} = nothing,
        rs::Union{Nothing, Real} = nothing,
        bw::Union{Nothing, Real} = nothing,
        w::Union{Nothing, AbstractVector} = nothing,
    )::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}

    @assert fs >= 1 "fs must be ≥ 1."
    nqf = div(fs, 2)

    # check parameters

    _check_var(
        fprototype,
        [:fir, :firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch],
        "fprototype"
    )
    !isnothing(ftype) && _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")
    if fprototype === :fir
        @assert !(isnothing(order) && isnothing(w)) "Either order or w must be specified."
        if !isnothing(w)
            ftype in [:hp, :bp, :bs] && @assert mod(length(w), 2) != 0 "Length of w must be odd."
            @assert length(w) >= 1 "Length of w must be ≥ 1."
            order = length(w)
        elseif !isnothing(order)
            ftype in [:hp, :bp, :bs] && @assert mod(order, 2) != 0 "order must be odd."
            w = DSP.hamming(order)
        end
        @assert length(w) == order "Length of w ($(length(w))) and order ($order) must be equal."
    end
    if fprototype in [:firls, :remez, :iirnotch]
        @assert !isnothing(bw) "bw must be specified."
        @assert bw > 0 "bw must be > 0."
        @assert bw <= 10 "bw must be ≤ 10."
        if length(cutoff) == 1
            if bw >= cutoff
                bw = round(cutoff - 0.1, digits = 1)
                _info("bw truncated to $bw Hz")
            end
        else
            if bw >= cutoff[2]
                bw = round(cutoff[2] - 0.1, digits = 1)
                _info("bw truncated to $bw Hz")
            end
        end
    end
    if fprototype === :firls
        if ftype in [:bp, :bs]
            if !isnothing(w)
                @assert length(w) == 6 "Length of w must be 6."
            else
                w = ones(6)
            end
        elseif ftype in [:lp, :hp]
            if !isnothing(w)
                @assert length(w) == 4 "Length of w must be 4."
            else
                w = ones(4)
            end
        end
    end
    if fprototype in [:chebyshev1, :chebyshev2, :elliptic]
        if isnothing(rp)
            rp = 0.5
            _info("rp set at $rp Hz.")
        end
        if isnothing(rs)
            rs = 20
            _info("rs set at $rs Hz.")
        end
    end
    if fprototype in [:firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic]
        @assert !isnothing(order) "order must be specified."
        @assert !isnothing(ftype) "ftype must be specified."
    end
    if fprototype === :iirnotch
        !isnothing(ftype) && _info("For :iirnotch filter ftype is ignored")
        !isnothing(order) && _info("For :iirnotch filter order is ignored")
        @assert length(cutoff) == 1 "For :iirnotch filter cutoff must contain only one frequency."
    end
    if fprototype in [:fir, :butterworth, :chebyshev1, :chebyshev2, :elliptic]
        ftype in [:lp, :hp] && @assert length(cutoff) == 1 "For :$(ftype) filter, cutoff must specify only one frequency."
        ftype in [:bp, :bs] && @assert length(cutoff) == 2 "For :$(ftype) filter, cutoff must specify only one frequency."
    end
    if length(cutoff) == 1
        @assert cutoff > 0 "cutoff must be > 0 Hz."
        @assert cutoff < nqf "cutoff must be < $nqf Hz."
    else
        if cutoff[1] == cutoff[2]
            cutoff = (cutoff[1], cutoff[1] + 0.1)
        elseif cutoff[1] > cutoff[2]
            cutoff = (cutoff[2], cutoff[1])
        end
    end

    ## FIR filters

    if fprototype === :fir

        if ftype === :lp
            responsetype = Lowpass(cutoff)
        elseif ftype === :hp
            responsetype = Highpass(cutoff)
        elseif ftype === :bp
            responsetype = Bandpass(cutoff[1], cutoff[2])
        elseif ftype === :bs
            responsetype = Bandstop(cutoff[1], cutoff[2])
        end

        prototype = FIRWindow(w)

        if ftype in [:lp, :hp]
            ftype === :lp && _info("Creating LP filter:")
            ftype === :hp && _info("Creating HP filter:")
            _info(" Number of taps: $order")
        elseif ftype === :bp
            _info("Creating BP filter:")
            _info(" Number of taps: $order")
        elseif ftype === :bs
            _info("Creating BS filter:")
            _info(" Number of taps: $order")
        end

        flt = digitalfilter(responsetype, prototype, fs = fs)

        return flt

    elseif fprototype === :firls

        if ftype === :bp

            f1_stop = cutoff[1] - (bw / 2)
            f1_pass = cutoff[1] + (bw / 2)
            f2_pass = cutoff[2] - (bw / 2)
            f2_stop = cutoff[2] + (bw / 2)
            flt_shape = [0, 0, 1, 1, 0, 0]
            flt_frq = [0, f1_stop, f1_pass, f2_pass, f2_stop, nqf]

        elseif ftype === :bs

            f1_pass = cutoff[1] - (bw / 2)
            f1_stop = cutoff[1] + (bw / 2)
            f2_stop = cutoff[2] - (bw / 2)
            f2_pass = cutoff[2] + (bw / 2)
            flt_shape = [1, 1, 0, 0, 1, 1]
            flt_frq = [0, f1_pass, f1_stop, f2_stop, f2_pass, nqf]

        elseif ftype === :lp

            f_pass = cutoff[1] - (bw / 2)
            f_stop = cutoff[1] + (bw / 2)
            flt_shape = [1, 1, 0, 0]
            flt_frq = [0, f_pass, f_stop, nqf]

        elseif ftype === :hp

            f_pass = cutoff[1] + (bw / 2)
            f_stop = cutoff[1] - (bw / 2)
            flt_shape = [0, 0, 1, 1]
            flt_frq = [0, f_stop, f_pass, nqf]

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

        flt = FIRLSFilterDesign.firls_design((order - 1), flt_frq, flt_shape, w, true, fs = fs)

        return flt

    elseif fprototype === :remez

        if ftype === :bp
            f1_stop = cutoff[1] - (bw / 2)
            f1_pass = cutoff[1] + (bw / 2)
            f2_pass = cutoff[2] - (bw / 2)
            f2_stop = cutoff[2] + (bw / 2)
            w = [(0, f1_stop) => 0, (f1_pass, f2_pass) => 1, (f2_stop, nqf) => 0]
        elseif ftype === :bs
            f1_pass = cutoff[1] - (bw / 2)
            f1_stop = cutoff[1] + (bw / 2)
            f2_stop = cutoff[2] - (bw / 2)
            f2_pass = cutoff[2] + (bw / 2)
            w = [(0, f1_pass) => 1, (f1_stop, f2_stop) => 0, (f2_pass, nqf) => 1]
        elseif ftype === :lp
            f_pass = cutoff[1] - (bw / 2)
            f_stop = cutoff[1] + (bw / 2)
            w = [(0, f_pass) => 1, (f_stop, nqf) => 0]
        elseif ftype === :hp
            f_pass = cutoff[1] + (bw / 2)
            f_stop = cutoff[1] - (bw / 2)
            w = [(0, f_stop) => 0, (f_pass, nqf) => 1]
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

        flt = remez(order, w, Hz = fs, maxiter = 100)

        return flt

    end

    ## IIR filters

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic]

        if ftype === :lp
            responsetype = Lowpass(cutoff)
        elseif ftype === :hp
            responsetype = Highpass(cutoff)
        elseif ftype === :bp
            responsetype = Bandpass(cutoff[1], cutoff[2])
        elseif ftype === :bs
            responsetype = Bandstop(cutoff[1], cutoff[2])
        end

        if fprototype === :butterworth
            prototype = Butterworth(order)
        elseif fprototype === :chebyshev1
            prototype = Chebyshev1(order, rp)
        elseif fprototype === :chebyshev2
            prototype = Chebyshev2(order, rs)
        elseif fprototype === :elliptic
            prototype = Elliptic(order, rp, rs)
        end

        flt = digitalfilter(responsetype, prototype, fs = fs)

        return flt

    elseif fprototype === :iirnotch

        flt = iirnotch(cutoff[1], bw, fs = fs)

        return flt

    end

end

"""
    filter_apply(s; <keyword arguments>)

Apply IIR or FIR filter.

# Arguments

  - `s::AbstractVector`
  - `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter
  - `dir:Symbol=:twopass`: filtering direction:
      + `:twopass`: two passes, the resulting signal has zero phase distortion, the effective filter order is doubled
      + `:onepass`: single pass
      + `:reverse`: one pass, reverse direction

# Returns

  - `s_new::Vector{Float64}`
"""
function filter_apply(
        s::AbstractVector;
        flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}},
        dir::Symbol = :twopass,
    )::Vector{Float64}

    _check_var(dir, [:twopass, :onepass, :reverse], "dir")

    dir === :onepass && (return filt(flt, s))
    dir === :twopass && (return filtfilt(flt, s))
    return dir === :reverse && (return filt(flt, reverse(s)))

end

"""
    filter_apply(obj; <keyword arguments>)

Apply IIR or FIR filter.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`
  - `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter
  - `dir:Symbol=:twopass`: filtering direction:
      + `:twopass`: two passes, the resulting signal has zero phase distortion, the effective filter order is doubled
      + `:onepass`: single pass
      + `:reverse`: one pass, reverse direction

# Returns

  - `obj_new::NeuroAnalyzer.NEURO`
"""
function filter_apply(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}},
        dir::Symbol = :twopass,
    )::NeuroAnalyzer.NEURO

    _check_var(dir, [:twopass, :onepass, :reverse], "dir")

    ch = get_channel(obj, ch = ch)
    ep_n = nepochs(obj)

    ep_n > 1 && _warn("filter() should be applied to a continuous signal.")
    _info("Signal should be tapered prior to filtering to reduce edge artifacts")
    dir === :twopass && _info("Filter is applied twice, the effective filter order is doubled")

    obj_new = deepcopy(obj)

    # initialize progress bar
    progbar = Progress(ep_n * length(ch), dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(ch)
            obj_new.data[ch[ch_idx], :, ep_idx] = @views filter_apply(
                obj.data[ch[ch_idx], :, ep_idx],
                flt = flt,
                dir = dir,
            )
            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    push!(obj_new.history, "filter_apply(OBJ, ch=$ch, dir=$dir)")

    return obj_new

end

"""
    filter_apply!(obj; <keyword arguments>)

Apply IIR or FIR filter.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`
  - `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter
  - `dir:Symbol=:twopass`: filtering direction:
      + `:twopass`: two passes, the resulting signal has zero phase distortion, the effective filter order is doubled
      + `:onepass`: single pass
      + `:reverse`: one pass, reverse direction

# Returns

  - `Nothing`
"""
function filter_apply!(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}},
        dir::Symbol = :twopass,
    )::Nothing

    obj_new = filter_apply(obj, ch = ch, flt = flt, dir = dir)

    obj.data = obj_new.data
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
      + `:fir`: FIR filter
      + `:firls`: weighted least-squares FIR filter
      + `:remez`: Remez FIR filter
      + `:butterworth`: IIR filter
      + `:chebyshev1` IIR filter
      + `:chebyshev2` IIR filter
      + `:elliptic` IIR filter
      + `:iirnotch`: second-order IIR notch filter
  - `ftype::Union{Nothing, Symbol}=nothing`: filter type:
      + `:lp`: low pass
      + `:hp`: high pass
      + `:bp`: band pass
      + `:bs`: band stop
  - `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (must be a pair of frequencies for `:bp` and `:bs`)
  - `fs::Int64`: signal sampling rate
  - `order::Union{Nothing, Int64}=nothing`: filter order
  - `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.5 dB
  - `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 20 dB
  - `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
  - `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
  - `dir:Symbol=:twopass`: filtering direction:
      + `:twopass`: two passes, the resulting signal has zero phase distortion, the effective filter order is doubled
      + `:onepass`: single pass
      + `:reverse`: one pass, reverse direction
  - `preview::Bool=false`: plot filter response

# Returns

  - `obj_new::NeuroAnalyzer.NEURO`

If `preview=true`, it will return `GLMakie.Figure`.
"""
function filter(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        fprototype::Symbol,
        ftype::Union{Nothing, Symbol} = nothing,
        cutoff::Union{Real, Tuple{Real, Real}},
        order::Union{Nothing, Int64} = nothing,
        rp::Union{Nothing, Real} = nothing,
        rs::Union{Nothing, Real} = nothing,
        bw::Union{Nothing, Real} = nothing,
        w::Union{Nothing, AbstractVector} = nothing,
        dir::Symbol = :twopass,
        preview::Bool = false,
    )::Union{NeuroAnalyzer.NEURO, GLMakie.Figure}

    if preview
        _info("Previewing filter response, signal will not be filtered")
        fprototype === :iirnotch && (ftype = :bs)
        p = plot_filter(
            fs = sr(obj),
            fprototype = fprototype,
            ftype = ftype,
            cutoff = cutoff,
            order = order,
            rp = rp,
            rs = rs,
            bw = bw,
            w = w,
        )
        return p
    else
        flt = filter_create(
            fprototype = fprototype,
            ftype = ftype,
            cutoff = cutoff,
            fs = sr(obj),
            order = order,
            rp = rp,
            rs = rs,
            bw = bw,
            w = w,
        )
        obj_new = filter_apply(obj, ch = ch, flt = flt, dir = dir)

        return obj_new

    end
end

"""
    filter!(obj; <keyword arguments>)

Apply filtering.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `fprototype::Symbol`: filter prototype:
      + `:fir`: FIR filter
      + `:firls`: weighted least-squares FIR filter
      + `:remez`: Remez FIR filter
      + `:butterworth`: IIR filter
      + `:chebyshev1` IIR filter
      + `:chebyshev2` IIR filter
      + `:elliptic` IIR filter
      + `:iirnotch`: second-order IIR notch filter
  - `ftype::Union{Nothing, Symbol}=nothing`: filter type:
      + `:lp`: low pass
      + `:hp`: high pass
      + `:bp`: band pass
      + `:bs`: band stop
  - `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (must be a pair of frequencies for `:bp` and `:bs`)
  - `fs::Int64`: signal sampling rate
  - `order::Union{Nothing, Int64}=nothing`: filter order
  - `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.5 dB
  - `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 20 dB
  - `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
  - `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
  - `dir:Symbol=:twopass`: filtering direction:
      + `:twopass`: two passes, the resulting signal has zero phase distortion, the effective filter order is doubled
      + `:onepass`: single pass
      + `:reverse`: one pass, reverse direction
  - `preview::Bool=false`: plot filter response

# Returns

  - `Nothing`

If `preview=true`, it will return `GLMakie.Figure`.
"""
function filter!(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        fprototype::Symbol,
        ftype::Union{Symbol, Nothing} = nothing,
        cutoff::Union{Real, Tuple{Real, Real}},
        order::Union{Nothing, Int64} = nothing,
        rp::Union{Nothing, Real} = nothing,
        rs::Union{Nothing, Real} = nothing,
        bw::Union{Nothing, Real} = nothing,
        w::Union{Nothing, AbstractVector} = nothing,
        dir::Symbol = :twopass,
        preview::Bool = false,
    )::Union{Nothing, GLMakie.Figure}

    if preview
        _info("Previewing filter response, signal will not be filtered")
        fprototype === :iirnotch && (ftype = :bs)
        p = plot_filter(
            fs = sr(obj),
            fprototype = fprototype,
            ftype = ftype,
            cutoff = cutoff,
            order = order,
            rp = rp,
            rs = rs,
            bw = bw,
            w = w,
        )
        return p
    end

    obj_new = NeuroAnalyzer.filter(
        obj,
        ch = ch,
        fprototype = fprototype,
        ftype = ftype,
        cutoff = cutoff,
        order = order,
        rp = rp,
        rs = rs,
        bw = bw,
        dir = dir,
        w = w,
    )

    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
