export plot_filter

"""
    plot_filter(<keyword arguments>)

Plot filter response.

# Arguments

  - `fs::Int64`: sampling rate
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
  - `order::Union{Nothing, Int64}=nothing`: filter order
  - `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  - `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  - `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
  - `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
  - `flim::Tuple{Real, Real}=(0, 0): frequency limit for the X-axis
  - `mono::Bool=false`: use color or gray palette
  - `gui::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_filter(;
    fs::Int64,
    fprototype::Symbol,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Real, Tuple{Real, Real}},
    order::Union{Nothing, Int64}=nothing,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    bw::Union{Nothing, Real} = nothing,
    w::Union{Nothing, AbstractVector} = nothing,
    mono::Bool = false,
    flim::Tuple{Real, Real} = (0, fs / 2),
)::GLMakie.Figure

    _check_tuple(flim, (0, fs / 2), "flim")
    @assert fs >= 1 "fs must be ≥ 1."

    nqf = div(fs, 2)
    nqf > flim[2] && (nqf = flim[2])

    # check parameters

    _check_var(
        fprototype,
        [:fir, :firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch],
        "fprototype"
    )
    !isnothing(ftype) && _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")
    @assert fs >= 1 "fs must be ≥ 1."
    if fprototype === :fir
        @assert isnothing(order) || isnothing(w) "Either order or w must be specified."
        if !isnothing(w)
            ftype in [:hp, :bp, :bs] && @assert mod(w, 2) != 0 "Length of w must be odd."
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
        if length(cutoff) == 1
            if bw >= cutoff
                bw = cutoff - 0.1
                _info("bw truncated to $bw Hz")
            end
        else
            if bw >= cutoff[2]
                bw = cutoff[2] - 0.1
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
        elseif ftype in [:lp, :hs]
            if !isnothing(w)
                @assert length(w) == 6 "Length of w must be 6."
            else
                w = ones(6)
            end
        end
    end
    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic]
        isnothing(rp) && (rp = fprototype === :elliptic ? 0.0025 : 2)
        isnothing(rs) && (rp = fprototype === :elliptic ? 40 : 20)
    end
    if fprototype !== :iirnotch
        @assert !isnothing(order) "order must be specified."
        @assert !isnothing(ftype) "ftype must be specified."
    end
    if fprototype === :iirnotch
        !isnothing(ftype) && _info("For :iirnotch filter ftype is ignored")
        !isnothing(order) && _info("For :iirnotch filter order is ignored")
        @assert length(cutoff) == 1 "For :iirnotch filter cutoff must contain only one frequency."
    end
    if fprototype in [:fir, :butterworth, :chebyshev1, :chebyshev2, :elliptic]
        if ftype === :lp
            @assert length(cutoff) == 1 "For :lp filter, cutoff must specify only one frequency."
            responsetype = Lowpass(cutoff)
        elseif ftype === :hp
            @assert length(cutoff) == 1 "For :hp filter, cutoff must specify only one frequency."
            responsetype = Highpass(cutoff)
        elseif ftype === :bp
            @assert length(cutoff) == 2 "For :bp filter, cutoff must specify two frequencies."
            responsetype = Bandpass(cutoff[1], cutoff[2])
        elseif ftype === :bs
            @assert length(cutoff) == 2 "For :bs filter, cutoff must specify two frequencies."
            responsetype = Bandstop(cutoff[1], cutoff[2])
        end
    end
    if length(cutoff) == 1
        @assert cutoff > 0 "cutoff must be > 0 Hz."
        @assert cutoff < nqf "cutoff must be < $nqf Hz."
    else
        _check_tuple(cutoff, (0, nqf), "cutoff")
    end

    # create observables

    cutoff = Observable(float.(cutoff))
    order = Observable(order)
    !isnothing(rp) && (rp = Observable(float(rp)))
    !isnothing(rs) && (rs = Observable(float(rs)))
    fprototype in [:firls, :remez, :iirnotch] && @assert !isnothing(bw) "bw must be specified."
    !isnothing(bw) && (bw = Observable(float(bw)))
    # prepare plot
    plot_size = (1400, 900)
    p = GLMakie.Figure(; size = plot_size)
    grid = p[4, 1] = GridLayout()

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic]
# cutoff
# order
# rp
# rs
    elseif fprototype in [:remez]
# cutoff
# order OR window
# bw
    elseif fprototype in [:fir]
# cutoff
# order OR window
# bw
    elseif fprototype in [:firls]
# cutoff
# order OR window
# bw
    elseif fprototype in [:iirnotch]
        _ = Label(
                grid[1, 1],
                "Cutoff [Hz]",
                fontsize = 15,
            )
        sl_cutoff = Slider(
                        grid[1, 2],
                        range = 0.5:0.1:(nqf - 0.1),
                        startvalue = cutoff[],
                        horizontal = true,
                    )
        on(sl_cutoff.value, priority=1) do val
            cutoff[] = round(val, digits=1)
            if cutoff[] > 10
                sl_bw.range = 0.1:0.1:10
            else
                if bw[] >= cutoff[]
                    bw[] = cutoff[] - 0.1
                    set_close_to!(sl_bw, bw[])
                end
                sl_bw.range = 0.1:0.1:(cutoff[] - 0.1)
            end
            notify(cutoff)
        end
        _ = Label(
                grid[2, 1],
                "Band width",
                fontsize = 15,
            )
        sl_bw = Slider(
                    grid[2, 2],
                    range = cutoff[] > 10 ? (0.1:0.1:10) : (0.1:0.1:(cutoff[] - 0.1)),
                    startvalue = bw[],
                    horizontal = true,
                )
        on(sl_bw.value) do val
            bw[] = round(val, digits=1)
            notify(bw)
        end
    end

    if fprototype in [:firls, :remez]
        _ = Label(
                grid[2, 1],
                "Band width",
                fontsize = 15,
            )
        if length(cutoff[]) == 1
            sl_bw = IntervalSlider(
                                grid[1, 2],
                                range = 0.1:0.1:(cutoff[] - 0.1),
                                startvalues = bw[],
                                horizontal = true,
                            )
        else
            sl_bw = IntervalSlider(
                                grid[1, 2],
                                range = 0.1:0.1:(cutoff[][2] - 0.1),
                                startvalues = bw[],
                                horizontal = true,
                            )
        end
        on(sl_bw.value) do val
            bw[] = val
            notify(bw)
        end
    end

    if ftype in [:hp, :lp]
        sl_cutoff = Slider(
                        grid[1, 2],
                        range = 0.1:0.1:(nqf - 0.1),
                        startvalue = cutoff,
                        horizontal = true,
                    )
        _ = Label(
                grid[2, 1],
                "Order (taps)",
                fontsize = 15,
            )
        sl_order = Slider(
                        grid[2, 2],
                        range = 1:1:500,
                        startvalue = order,
                        horizontal = true,
                    )
        on(sl_cutoff.value) do val
            cutoff[] = val
            notify(cutoff)
        end

        on(sl_order.value) do val
            order[] = val
            notify(order)
        end
    elseif ftype in [:bp, :bs]
        sl_cutoff = IntervalSlider(
                                grid[1, 2],
                                range = 0.1:0.1:(nqf - 0.1),
                                startvalues = cutoff,
                                horizontal = true,
                            )
        _ = Label(
                grid[2, 1],
                "Order (taps)",
                fontsize = 15,
            )
        sl_order = Slider(
                        grid[2, 2],
                        range = 1:2:501,
                        startvalue = order,
                        horizontal = true,
                    )
        on(sl_cutoff.interval) do val
            cutoff[] = val
            notify(cutoff)
        end

        on(sl_order.value) do val
            order[] = val
            notify(order)
        end
    end

    flt = @lift(
            filter_create(;
                fprototype = fprototype,
                ftype = ftype,
                cutoff = $cutoff,
                fs = fs,
                order = $order,
                rp = !isnothing(rp) ? $rp : nothing,
                rs = !isnothing(rs) ? $rs : nothing,
                bw = !isnothing(bw) ? $bw : nothing,
                w = w,
            )
        )

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch]

        fresp = lift(DSP.freqresp, flt)
        # convert to dB
        H = @lift(real.(20 * log10.(abs.($fresp[1]))))
        # convert rad/sample to Hz
        f = @lift(round.($fresp[2] .* fs / 2 / pi, digits=1))

        if fprototype !== :iirnotch
            fname = titlecase(String(fprototype))
            title = @lift("Filter: $(fname), type: $(uppercase(String(ftype)))\ncutoff: $(round($cutoff, digits=1)) Hz, order: $($order)\n\nFrequency response")
        else
            fname = "IIR notch"
            title = @lift("Filter: $(fname), cutoff: $(round($cutoff, digits=1)) Hz, transition band width: $(round($bw, digits=1)) Hz\n\nFrequency response")
        end

        ax1 = GLMakie.Axis(
            p[1, 1];
            xlabel = "Frequency [Hz]",
            ylabel = "Magnitude [dB]",
            title = title,
            xticks = length(flim[1]:0.1:flim[2]) > 20 ? LinearTicks(10) : LinearTicks(20),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax1, flim)
        GLMakie.ylims!(ax1, (-100, 20))
        ax1.titlesize = 18
        ax1.xlabelsize = 18
        ax1.ylabelsize = 18
        ax1.xticklabelsize = 12
        ax1.yticklabelsize = 12

        GLMakie.lines!(
                    ax1,
                    f,
                    H,
                    color = mono ? :black : :blue,
                )

        phresp = lift(DSP.phaseresp, flt)
        #phi = @lift(rad2deg.($phresp[1]))
        phi = @lift(($phresp[1]))
        # convert rad/sample to Hz
        f = @lift(round.($phresp[2] .* fs / 2 / pi, digits=1))

        ax2 = GLMakie.Axis(
            p[2, 1];
            xlabel = "Frequency [Hz]",
            ylabel = "Phase [rad]",
            title = "Phase response",
            xticks = LinearTicks(15),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax2, flim)
        ax2.titlesize = 18
        ax2.xlabelsize = 18
        ax2.ylabelsize = 18
        ax2.xticklabelsize = 12
        ax2.yticklabelsize = 12

        GLMakie.lines!(
                    ax2,
                    f,
                    phi,
                    color = mono ? :black : :blue,
                    nan_color = mono ? :black : :blue,
                )

        phresp = lift(DSP.phaseresp, flt)
        tau = @lift(-derivative(rad2deg.($phresp[1])))

        ax3 = GLMakie.Axis(
            p[3, 1];
            xlabel = "Frequency [Hz]",
            ylabel = "Group delay [samples]",
            title = "Group delay",
            xticks = LinearTicks(15),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax3, flim)
        ax3.titlesize = 18
        ax3.xlabelsize = 18
        ax3.ylabelsize = 18
        ax3.xticklabelsize = 12
        ax3.yticklabelsize = 12

        GLMakie.lines!(
                    ax3,
                    f,
                    tau,
                    color = mono ? :black : :blue,
                )

    else
        f = range(0; stop = pi, length = 1024)
        H = _fir_response(flt[], f)
        # convert to dB
        H = amp2db.(abs.(H))
        # convert rad/sample to Hz
        f = Observable(collect(f .* fs / 2 / pi))

        if fprototype === :fir
            title = "Filter: FIR, type: $(uppercase(String(ftype))), cutoff: $(round.(cutoff, digits=1)) Hz, order: $order\n\nFrequency response"
        elseif fprototype === :firls
            title = "Filter: FIR (LS), type: $(uppercase(String(ftype))), cutoff: $(round.(cutoff, digits=1)) Hz, transition band width: $bw Hz, order: $order\n\nFrequency response"
        elseif fprototype === :remez
            title = "Filter: Remez, type: $(uppercase(String(ftype))), cutoff: $(round.(cutoff, digits=1)) Hz, transition band width: $bw Hz, order: $order\n\nFrequency response"
        end

        ax1 = GLMakie.Axis(
            p[1, 1];
            xlabel = "Frequency [Hz]",
            ylabel = "Magnitude\n[dB]",
            title = title,
            xticks = LinearTicks(15),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax1, flim)
        GLMakie.ylims!(ax1, (-100, 20))
        ax1.titlesize = 18
        ax1.xlabelsize = 18
        ax1.ylabelsize = 18
        ax1.xticklabelsize = 12
        ax1.yticklabelsize = 12

        GLMakie.lines!(
            ax1,
            w[],
            H[];
            colormap = pal,
        )

        w = range(0; stop = pi, length = 1024)
        phi = _fir_response(flt[], w)
        phi = rad2deg.(-atan.(imag(phi), real(phi)))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi

        ax2 = GLMakie.Axis(
            p[2, 1];
            xlabel = "Frequency [Hz]",
            ylabel = "Phase\n[deg]",
            title = "Phase response",
            xticks = LinearTicks(15),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax2, flim)
        ax2.titlesize = 18
        ax2.xlabelsize = 18
        ax2.ylabelsize = 18
        ax2.xticklabelsize = 12
        ax2.yticklabelsize = 12

        GLMakie.lines!(ax2, w, phi; colormap = pal)

        tau = -derivative(phi)

        ax3 = GLMakie.Axis(
            p[3, 1];
            xlabel = "Frequency [Hz]",
            ylabel = "Group delay\n[samples]",
            title = "Group delay",
            xticks = LinearTicks(15),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax3, flim)
        ax3.titlesize = 18
        ax3.xlabelsize = 18
        ax3.ylabelsize = 18
        ax3.xticklabelsize = 12
        ax3.yticklabelsize = 12

        GLMakie.lines!(ax3, w, tau; colormap = pal)

    end

    if length(cutoff[]) == 1
        GLMakie.vlines!(
                    ax1,
                    cutoff,
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :red,
                )
        GLMakie.vlines!(
                    ax2,
                    cutoff,
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :red,
                )
        GLMakie.vlines!(
                    ax3,
                    cutoff,
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :red,
                )
    else
        GLMakie.vlines!(
                    ax1,
                    cutoff[1],
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :red,
                )
        GLMakie.vlines!(
                    ax1,
                    cutoff[2],
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :green,
                )
        GLMakie.vlines!(
                    ax2,
                    cutoff[1],
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :red,
                )
        GLMakie.vlines!(
                    ax2,
                    cutoff[2],
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :green,
                )
        GLMakie.vlines!(
                    ax3,
                    cutoff[1],
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :red,
                )
        GLMakie.vlines!(
                    ax3,
                    cutoff[2],
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :green,
                )
    end


    return p

end

"""
    plot_filter(obj, <keyword arguments>)

Plot filter response.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `n::Int64`: signal length in samples
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
  - `order::Union{Nothing, Int64}=nothing`: filter order
  - `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  - `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  - `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
  - `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
  - `mono::Bool=false`: use color or gray palette
  - `flim::Tuple{Real, Real}=(0, 0): frequency limit for the X-axis

# Returns

  - `p::GLMakie.Figure`
"""
function plot_filter(
    obj::NeuroAnalyzer.NEURO;
    fprototype::Symbol,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Real, Tuple{Real, Real}},
    order::Union{Nothing, Int64}=nothing,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    bw::Union{Nothing, Real} = nothing,
    w::Union{Nothing, AbstractVector} = nothing,
    mono::Bool = false,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
)::GLMakie.Figure

    p = plot_filter(;
        fs = sr(obj),
        fprototype = fprototype,
        ftype = ftype,
        cutoff = cutoff,
        order = order,
        rp = rp,
        rs = rs,
        bw = bw,
        w = w,
        mono = mono,
        flim = flim,
    )

    return p

end
