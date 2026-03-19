export filter_create
export filter_apply
export filter_apply!
export filter
export filter!

"""
    filter_create(; <keyword arguments>)

Create a FIR or IIR filter object.

# Arguments

- `fprototype::Symbol`: filter prototype:
    - `:fir`: FIR filter (window method)
    - `:firls`: weighted least-squares FIR filter
    - `:remez`: Remez (Parks-McClellan) FIR filter
    - `:butterworth`: Butterworth IIR filter
    - `:chebyshev1`: Chebyshev type-I IIR filter
    - `:chebyshev2`: Chebyshev type-II IIR filter
    - `:elliptic`: elliptic IIR filter
    - `:iirnotch`: second-order IIR notch filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: cutoff frequency/ies in Hz; scalar for `:lp`/`:hp`; 2-tuple for `:bp`/`:bs`
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `order::Union{Nothing, Int64}=nothing`: filter order
- `rp::Union{Nothing, Real}=nothing`: pass-band ripple in dB (default 0.5 dB)
- `rs::Union{Nothing, Real}=nothing`: stop-band attenuation in dB (default 20 dB)
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz (required for `:firls`, `:remez`, `:iirnotch`)
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` (default: Hamming) or weight vector for `:firls`

# Returns

- `Vector{Float64}`: FIR filter coefficients (for `:fir`, `:firls`, `:remez`)
- `ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}`: IIR filter in zero-pole-gain form (for `:butterworth`, `:chebyshev1`, `:chebyshev2`, `:elliptic`)
- `Biquad{:z, Float64}`: second-order biquad filter (for `:iirnotch`)

# Throws

- `ArgumentError`: if any required argument is missing or invalid

# See also

[`filter_apply`](@ref), [`filter`](@ref)
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

    !(fs >= 1) && throw(ArgumentError("fs must be ≥ 1."))
    nqf = div(fs, 2)

    # check parameters
    _check_var(
        fprototype,
        [:fir, :firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch],
        "fprototype"
    )
    !isnothing(ftype) && _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")

    # --- :fir parameter validation ---
    if fprototype === :fir
        isnothing(bw) && throw(ArgumentError("bw must be specified for $fprototype."))
        if !isnothing(w)
            ftype in (:hp, :bp, :bs) && mod(length(w), 2) == 0 &&
                throw(ArgumentError("Length of w must be odd for :hp/:bp/:bs filters."))
            length(w) >= 1 || throw(ArgumentError("Length of w must be ≥ 1."))
            order = length(w)
        elseif !isnothing(order)
            ftype in (:hp, :bp, :bs) && mod(order, 2) == 0 &&
                throw(ArgumentError("order must be odd for :hp/:bp/:bs filters."))
            w = DSP.hamming(order)
        end
        length(w) == order || throw(ArgumentError("Length of w ($(length(w))) must equal order ($order)."))
    end

    # --- :firls / :remez / :iirnotch bw validation ---
    if fprototype in [:firls, :remez, :iirnotch]
        isnothing(bw) && throw(ArgumentError("bw must be specified for $fprototype."))
        bw > 0 || throw(ArgumentError("bw must be > 0."))
        bw <= 10 || throw(ArgumentError("bw must be ≤ 10."))
        if length(cutoff) == 1
            if bw >= cutoff
                bw = round(cutoff - 0.1, digits=1)
                _info("bw truncated to $bw Hz")
            end
        else
            if bw >= cutoff[2]
                bw = round(cutoff[2] - 0.1, digits=1)
                _info("bw truncated to $bw Hz")
            end
        end
    end

    # --- :firls weight vector defaults ---
    if fprototype === :firls
        if ftype in [:bp, :bs]
            if !isnothing(w)
                !(length(w) == 6) && throw(ArgumentError("Length of w must be 6 for :bp/:bs filter."))
            else
                w = ones(6)
            end
        elseif ftype in [:lp, :hp]
            if !isnothing(w)
                !(length(w) == 4) && throw(ArgumentError("Length of w must be 4 for :lp/:hp filter."))
            else
                w = ones(4)
            end
        end
    end

    # --- ripple defaults for equiripple IIR prototypes ---
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

    # --- order and ftype required for these prototypes ---
    if fprototype in [:firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic]
        isnothing(order) && throw(ArgumentError("order must be specified for $fprototype."))
        isnothing(ftype) && throw(ArgumentError("ftype must be specified for $fprototype."))
    end

    # --- :iirnotch specifics ---
    if fprototype === :iirnotch
        !isnothing(ftype) && _info("For :iirnotch filter ftype is ignored")
        !isnothing(order) && _info("For :iirnotch filter order is ignored")
        length(cutoff) == 1 || throw(ArgumentError("cutoff must be a scalar for :iirnotch."))
    end

    # --- cutoff arity check ---
    if fprototype in (:fir, :butterworth, :chebyshev1, :chebyshev2, :elliptic)
        if ftype in (:lp, :hp)
            length(cutoff) == 1 || throw(ArgumentError("For :$ftype, cutoff must be a scalar."))
        elseif ftype in (:bp, :bs)
            length(cutoff) == 2 || throw(ArgumentError("For :$ftype, cutoff must specify two frequencies."))
        end
    end

    # --- cutoff value checks and normalization ---
    if length(cutoff) == 1
        cutoff > 0   || throw(ArgumentError("cutoff must be > 0 Hz."))
        cutoff < nqf || throw(ArgumentError("cutoff must be < $nqf Hz (Nyquist)."))
    else
        if cutoff[1] == cutoff[2]
            cutoff = (cutoff[1], cutoff[1] + 0.1)
        elseif cutoff[1] > cutoff[2]
            cutoff = (cutoff[2], cutoff[1])
        end
    end

    # -----------------------------------------------------------------------
    # FIR filters
    # -----------------------------------------------------------------------

    if fprototype === :fir
        responsetype = if ftype === :lp; Lowpass(cutoff)
                       elseif ftype === :hp; Highpass(cutoff)
                       elseif ftype === :bp; Bandpass(cutoff[1], cutoff[2])
                       elseif ftype === :bs; Bandstop(cutoff[1], cutoff[2])
                       end
        _info("Creating $(uppercase(string(ftype))) FIR filter ($(order) taps)")
        return digitalfilter(responsetype, FIRWindow(w); fs=fs)
    end

    if fprototype === :firls
        if ftype === :bp
            f1_stop, f1_pass = cutoff[1] - bw/2, cutoff[1] + bw/2
            f2_pass, f2_stop = cutoff[2] - bw/2, cutoff[2] + bw/2
            flt_shape = [0, 0, 1, 1, 0, 0]
            flt_frq = [0, f1_stop, f1_pass, f2_pass, f2_stop, nqf]
            _info("Creating BP FIRLS filter ($order taps, bw=$bw Hz)")
            _info(" Bands: stop=[0,$f1_stop], pass=[$f1_pass,$f2_pass], stop=[$f2_stop,$nqf]")
        elseif ftype === :bs
            f1_pass, f1_stop = cutoff[1] - bw/2, cutoff[1] + bw/2
            f2_stop, f2_pass = cutoff[2] - bw/2, cutoff[2] + bw/2
            flt_shape = [1, 1, 0, 0, 1, 1]
            flt_frq = [0, f1_pass, f1_stop, f2_stop, f2_pass, nqf]
            _info("Creating BS FIRLS filter ($order taps, bw=$bw Hz)")
            _info(" Bands: stop=[0,$f1_pass], pass=[$f1_stop,$f2_stop], stop=[$f2_pass,$nqf]")
        elseif ftype === :lp
            f_pass, f_stop = cutoff - bw/2, cutoff + bw/2
            flt_shape = [1, 1, 0, 0]
            flt_frq   = [0, f_pass, f_stop, nqf]
            _info("Creating LP FIRLS filter ($order taps, bw=$bw Hz, pass=$f_pass, stop=$f_stop)")
        elseif ftype === :hp
            f_stop, f_pass = cutoff - bw/2, cutoff + bw/2
            flt_shape = [0, 0, 1, 1]
            flt_frq   = [0, f_stop, f_pass, nqf]
            _info("Creating HP FIRLS filter ($order taps, bw=$bw Hz, stop=$f_stop, pass=$f_pass)")
        end
        return FIRLSFilterDesign.firls_design(order - 1, flt_frq, flt_shape, w, true; fs=fs)
    end

    if fprototype === :remez
        if ftype === :bp
            f1_stop, f1_pass = cutoff[1] - bw/2, cutoff[1] + bw/2
            f2_pass, f2_stop = cutoff[2] - bw/2, cutoff[2] + bw/2
            w = [(0, f1_stop) => 0, (f1_pass, f2_pass) => 1, (f2_stop, nqf) => 0]
            _info("Creating BP Remez filter ($order taps, bw=$bw Hz)")
            _info(" Bands: stop=[0,$f1_stop], pass=[$f1_pass,$f2_pass], stop=[$f2_stop,$nqf]")
        elseif ftype === :bs
            f1_pass, f1_stop = cutoff[1] - bw/2, cutoff[1] + bw/2
            f2_stop, f2_pass = cutoff[2] - bw/2, cutoff[2] + bw/2
            w = [(0, f1_pass) => 1, (f1_stop, f2_stop) => 0, (f2_pass, nqf) => 1]
            _info("Creating BS Remez filter ($order taps, bw=$bw Hz)")
            _info(" Bands: stop=[0,$f1_pass], pass=[$f1_stop,$f2_stop], stop=[$f2_pass,$nqf]")
        elseif ftype === :lp
            f_pass, f_stop = cutoff - bw/2, cutoff + bw/2
            w = [(0, f_pass) => 1, (f_stop, nqf) => 0]
            _info("Creating LP Remez filter ($order taps, bw=$bw Hz, pass=$f_pass, stop=$f_stop)")
        elseif ftype === :hp
            f_stop, f_pass = cutoff - bw/2, cutoff + bw/2
            w = [(0, f_stop) => 0, (f_pass, nqf) => 1]
            _info("Creating HP Remez filter ($order taps, bw=$bw Hz, stop=$f_stop, pass=$f_pass)")
        end
        _info("Creating $(uppercase(string(ftype))) Remez filter ($order taps, bw=$bw Hz)")
        return remez(order, w; Hz=fs, maxiter=100)
    end

    # -----------------------------------------------------------------------
    # IIR filters
    # -----------------------------------------------------------------------

    if fprototype in (:butterworth, :chebyshev1, :chebyshev2, :elliptic)
        responsetype = if ftype === :lp; Lowpass(cutoff)
                       elseif ftype === :hp; Highpass(cutoff)
                       elseif ftype === :bp; Bandpass(cutoff[1], cutoff[2])
                       elseif ftype === :bs; Bandstop(cutoff[1], cutoff[2])
                       end
        prototype = if fprototype === :butterworth; Butterworth(order)
                    elseif fprototype === :chebyshev1; Chebyshev1(order, rp)
                    elseif fprototype === :chebyshev2; Chebyshev2(order, rs)
                    elseif fprototype === :elliptic; Elliptic(order, rp, rs)
                    end
        _info("Creating $(uppercase(string(ftype))) $(fprototype) filter (order=$order)")
        return digitalfilter(responsetype, prototype; fs=fs)
    end

    if fprototype === :iirnotch
        _info("Creating IIR notch filter (cutoff=$(cutoff[1]) Hz, bw=$bw Hz)")
        return iirnotch(cutoff[1], bw; fs=fs)
    end

end

"""
    filter_apply(s; <keyword arguments>)

Apply a pre-designed IIR or FIR filter to a signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter object returned by [`filter_create`](@ref)
- `dir:Symbol=:twopass`: filtering direction:
    - `:twopass`: forward pass followed by reverse pass (zero phase distortion; effective filter order is doubled)
    - `:onepass`: single forward pass (introduces phase delay)
    - `:reverse`: single reverse pass

# Returns

- `Vector{Float64}`: filtered signal of the same length as `s`

# See also

[`filter_create`](@ref), [`filter_apply(::NeuroAnalyzer.NEURO)`](@ref)
"""
function filter_apply(
        s::AbstractVector;
        flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}},
        dir::Symbol = :twopass,
    )::Vector{Float64}

    _check_var(dir, [:twopass, :onepass, :reverse], "dir")

    if dir === :onepass
        return filt(flt, s)
    elseif dir === :twopass
        return filtfilt(flt, s)
    elseif dir === :reverse
        return filt(flt, reverse(s))
    end

end

"""
    filter_apply(obj; <keyword arguments>)

Apply a pre-designed filter to selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter object returned by [`filter_create`](@ref)
- `dir:Symbol=:twopass`: filtering direction:
    - `:twopass`: forward pass followed by reverse pass (zero phase distortion; effective filter order is doubled)
    - `:onepass`: single forward pass (introduces phase delay)
    - `:reverse`: single reverse pass

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: output NEURO object with filtered data

# Notes
- For best results apply to a continuous (single-epoch) signal. A warning is issued when `nepochs(obj) > 1`.
- Taper the signal before filtering to reduce edge artifacts.

# See also

[`filter_create`](@ref), [`filter_apply!`](@ref), [`filter`](@ref)
"""
function filter_apply(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}},
    dir::Symbol = :twopass,
)::NeuroAnalyzer.NEURO

    _check_var(dir, [:twopass, :onepass, :reverse], "dir")

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)
    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    ep_n > 1 && _warn("filter_apply() should preferably be used on a continuous signal.")
    _info("Taper the signal before filtering to reduce edge artifacts.")
    dir === :twopass && _info("Two-pass filtering: effective order is doubled.")

    obj_new = deepcopy(obj)

    # initialize progress bar
    progbar = Progress(ep_n * length(ch), dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        obj_new.data[ch[ch_idx], :, ep_idx] = @views filter_apply(
            obj.data[ch[ch_idx], :, ep_idx],
            flt = flt,
            dir = dir,
        )
        # update progress bar
        progress_bar && next!(progbar)
    end

    push!(obj_new.history, "filter_apply(OBJ, ch=$ch, dir=$dir)")

    return obj_new

end

"""
    filter_apply!(obj; <keyword arguments>)

Apply a pre-designed filter in-place to selected channels of a NEURO object.

Delegates to [`filter_apply`](@ref) and copies the result back.


# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter object returned by [`filter_create`](@ref)
- `dir:Symbol=:twopass`: filtering direction:
    - `:twopass`: forward pass followed by reverse pass (zero phase distortion; effective filter order is doubled)
    - `:onepass`: single forward pass (introduces phase delay)
    - `:reverse`: single reverse pass

# Returns

- `Nothing`

# See also

[`filter_apply`](@ref), [`filter!`](@ref)
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

Design and apply a digital filter to selected channels of a NEURO object in a single call.

Combines [`filter_create`](@ref) and [`filter_apply`](@ref). When `preview=true`, the filter frequency response is plotted without modifying the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `fprototype::Symbol`: filter prototype:
    - `:fir`: FIR filter (window method)
    - `:firls`: weighted least-squares FIR filter
    - `:remez`: Remez (Parks-McClellan) FIR filter
    - `:butterworth`: Butterworth IIR filter
    - `:chebyshev1`: Chebyshev type-I IIR filter
    - `:chebyshev2`: Chebyshev type-II IIR filter
    - `:elliptic`: elliptic IIR filter
    - `:iirnotch`: second-order IIR notch filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: cutoff frequency/ies in Hz; scalar for `:lp`/`:hp`; 2-tuple for `:bp`/`:bs`
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `order::Union{Nothing, Int64}=nothing`: filter order
- `rp::Union{Nothing, Real}=nothing`: pass-band ripple in dB (default 0.5 dB)
- `rs::Union{Nothing, Real}=nothing`: stop-band attenuation in dB (default 20 dB)
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz (required for `:firls`, `:remez`, `:iirnotch`)
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` (default: Hamming) or weight vector for `:firls`
- `dir:Symbol=:twopass`: filtering direction:
    - `:twopass`: forward pass followed by reverse pass (zero phase distortion; effective filter order is doubled)
    - `:onepass`: single forward pass (introduces phase delay)
    - `:reverse`: single reverse pass
- `preview::Bool=false`: if `true`, plot the filter frequency response and return the figure without filtering the signal

# Returns

- `NeuroAnalyzer.NEURO`: filtered object (when `preview=false`)
- `GLMakie.Figure`: filter frequency-response plot (when `preview=true`)

# See also

[`filter!`](@ref), [`filter_create`](@ref), [`filter_apply`](@ref)
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

Design and apply a digital filter in-place to selected channels of a NEURO object.

When `preview=true`, the filter frequency response is plotted and returned without modifying the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `fprototype::Symbol`: filter prototype:
    - `:fir`: FIR filter (window method)
    - `:firls`: weighted least-squares FIR filter
    - `:remez`: Remez (Parks-McClellan) FIR filter
    - `:butterworth`: Butterworth IIR filter
    - `:chebyshev1`: Chebyshev type-I IIR filter
    - `:chebyshev2`: Chebyshev type-II IIR filter
    - `:elliptic`: elliptic IIR filter
    - `:iirnotch`: second-order IIR notch filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: cutoff frequency/ies in Hz; scalar for `:lp`/`:hp`; 2-tuple for `:bp`/`:bs`
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `order::Union{Nothing, Int64}=nothing`: filter order
- `rp::Union{Nothing, Real}=nothing`: pass-band ripple in dB (default 0.5 dB)
- `rs::Union{Nothing, Real}=nothing`: stop-band attenuation in dB (default 20 dB)
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz (required for `:firls`, `:remez`, `:iirnotch`)
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` (default: Hamming) or weight vector for `:firls`
- `dir:Symbol=:twopass`: filtering direction:
    - `:twopass`: forward pass followed by reverse pass (zero phase distortion; effective filter order is doubled)
    - `:onepass`: single forward pass (introduces phase delay)
    - `:reverse`: single reverse pass
- `preview::Bool=false`: if `true`, plot the filter frequency response and return the figure without filtering the signal

# Returns

- `Nothing` when `preview=false`
- `GLMakie.Figure`: filter frequency-response plot (when `preview=true`)

# See also

[`filter`](@ref), [`filter_apply!`](@ref)
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
