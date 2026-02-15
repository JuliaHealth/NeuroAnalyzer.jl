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
  - `order::Union{Nothing, Int64}`: filter order
  - `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  - `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  - `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
  - `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
  - `mono::Bool=false`: use color or gray palette
  - `flim::Tuple{Real, Real}=(0, 0): frequency limit for the X-axis

# Returns

  - `p::GLMakie.Figure`
"""
function plot_filter(;
    fs::Int64,
    fprototype::Symbol,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Real, Tuple{Real, Real}},
    order::Union{Nothing, Int64},
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    bw::Union{Nothing, Real} = nothing,
    w::Union{Nothing, AbstractVector} = nothing,
    mono::Bool = false,
    flim::Tuple{Real, Real} = (0, fs / 2),
)::GLMakie.Figure

    pal = mono ? :grays : :darktest
    _check_tuple(flim, "flim", (0, fs / 2))

    cutoff = Observable(float.(cutoff))
    order = Observable(order)
    !isnothing(rp) && (rp = Observable(float(rp)))
    !isnothing(rs) && (rs = Observable(float(rs)))
    fprototype in [:firls, :remez, :iirnotch] && @assert !isnothing(bw) "bw must be specified."
    !isnothing(bw) && (bw = Observable(float(bw)))

    nqf = div(fs, 2)

    # prepare plot
    plot_size = (1400, 900)
    p = GLMakie.Figure(; size = plot_size)
    grid = p[4, 1] = GridLayout()

    lab_cutoff = Label(
                    grid[1, 1],
                    "Cutoff [Hz]",
                    fontsize = 15,
                )

    if ftype in [:hp, :lp]
        sl_cutoff = Slider(
                        grid[1, 2],
                        range = 0.1:0.1:(nqf - 0.1),
                        startvalue = cutoff,
                        horizontal = true,
                    )
        lab2 = Label(
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
    else
        sl_cutoff = IntervalSlider(
                                grid[1, 2],
                                range = 0.1:0.1:(nqf - 0.1),
                                startvalues = cutoff,
                                horizontal = true,
                            )
        lab2 = Label(
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

    if fprototype in [:firls, :remez, :iirnotch]
        lab_bw = Label(
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
        w = @lift(round.($fresp[2] .* fs / 2 / pi, digits=1))

        if fprototype !== :iirnotch
            fname = titlecase(String(fprototype))
            title = @lift("Filter: $(fname), type: $(uppercase(String(ftype)))\ncutoff: $($cutoff) Hz, order: $($order)\n\nFrequency response")
        else
            fname = "IIR notch"
            title = @lift("Filter: $(fname), cutoff: $($cutoff) Hz, transition band width: $($bw) Hz\n\nFrequency response")
        end

        ax1 = GLMakie.Axis(
            p[1, 1];
            xlabel = "Frequency [Hz]",
            ylabel = "Magnitude [dB]",
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
            w,
            H,
            colormap = pal,
        )

        if length(cutoff[]) == 1
            GLMakie.vlines!(ax1, cutoff; linestyle = :dash, linewidth = 0.5, colormap = pal)
        else
            GLMakie.vlines!(ax1, cutoff[1]; linestyle = :dash, linewidth = 0.5, color = :red, colormap = pal)
            GLMakie.vlines!(ax1, cutoff[2]; linestyle = :dash, linewidth = 0.5, color = :green, colormap = pal)
        end

        phi, w = phaseresp(flt[])
        phi = rad2deg.(phi)
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi

        ax2 = GLMakie.Axis(
            p[2, 1];
            xlabel = "Frequency [Hz]",
            ylabel = "Phase [deg]",
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

        if length(cutoff[]) == 1
            GLMakie.vlines!(ax2, cutoff[]; linestyle = :dash, linewidth = 0.5, colormap = pal)
        else
            GLMakie.vlines!(ax2, cutoff[][1]; linestyle = :dash, linewidth = 0.5, color = :red, colormap = pal)
            GLMakie.vlines!(ax2, cutoff[][2]; linestyle = :dash, linewidth = 0.5, color = :green, colormap = pal)
        end

        tau = -derivative(phi)

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

        GLMakie.lines!(ax3, w, tau; colormap = pal)

        if length(cutoff[]) == 1
            GLMakie.vlines!(ax3, cutoff[]; linestyle = :dash, linewidth = 0.5, colormap = pal)
        else
            GLMakie.vlines!(ax3, cutoff[][1]; linestyle = :dash, linewidth = 0.5, color = :red, colormap = pal)
            GLMakie.vlines!(ax3, cutoff[][2]; linestyle = :dash, linewidth = 0.5, color = :green, colormap = pal)
        end

    else
        f = range(0; stop = pi, length = 1024)
        H = _fir_response(flt[], f)
        # convert to dB
        H = amp2db.(abs.(H))
        # convert rad/sample to Hz
        w = Observable(collect(f .* fs / 2 / pi))
        @show H[]
        @show w[]

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

        if length(cutoff) == 1
            GLMakie.vlines!(ax1, cutoff; linestyle = :dash, linewidth = 0.5, colormap = pal)
        else
            GLMakie.vlines!(ax1, cutoff[1]; linestyle = :dash, linewidth = 0.5, color = :red, colormap = pal)
            GLMakie.vlines!(ax1, cutoff[2]; linestyle = :dash, linewidth = 0.5, color = :green, colormap = pal)
        end

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

        if length(cutoff) == 1
            GLMakie.vlines!(ax2, cutoff; linestyle = :dash, linewidth = 0.5, colormap = pal)
        else
            GLMakie.vlines!(ax2, cutoff[1]; linestyle = :dash, linewidth = 0.5, color = :red, colormap = pal)
            GLMakie.vlines!(ax2, cutoff[2]; linestyle = :dash, linewidth = 0.5, color = :green, colormap = pal)
        end

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

        if length(cutoff) == 1
            GLMakie.vlines!(ax3, cutoff; linestyle = :dash, linewidth = 0.5, colormap = pal)
        else
            GLMakie.vlines!(ax3, cutoff[1]; linestyle = :dash, linewidth = 0.5, color = :red, colormap = pal)
            GLMakie.vlines!(ax3, cutoff[2]; linestyle = :dash, linewidth = 0.5, color = :green, colormap = pal)
        end
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
  - `order::Union{Nothing, Int64}`: filter order
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
    order::Union{Nothing, Int64},
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
