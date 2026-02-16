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
  - `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.5 dB
  - `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 20 dB
  - `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
  - `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
  - `flim::Tuple{Real, Real} = (0, fs / 2)`: frequency limit
  - `mono::Bool=false`: use color or gray palette
  - `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

  - `p::GLMakie.Figure`
  - `f::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: if `gui=true`
"""
function plot_filter(;
    fs::Int64,
    fprototype::Symbol,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Real, Tuple{Real, Real}},
    order::Union{Nothing, Int64} = nothing,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    bw::Union{Nothing, Real} = nothing,
    w::Union{Nothing, AbstractVector} = nothing,
    flim::Tuple{Real, Real} = (0, fs / 2),
    mono::Bool = false,
    gui::Bool = true,
)::Union{GLMakie.Figure, Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}

    _check_tuple(flim, (0, fs / 2), "flim")
    @assert fs >= 1 "fs must be ≥ 1."
    v = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false

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
        @assert !(isnothing(order) && isnothing(w)) "Either order or w must be specified."
        if !isnothing(w)
            ftype in [:hp, :bp, :bs] && @assert mod(length(w), 2) != 0 "Length of w must be odd."
            @assert length(w) >= 1 "Length of w must be ≥ 1."
        elseif !isnothing(order)
            ftype in [:hp, :bp, :bs] && @assert mod(order, 2) != 0 "order must be odd."
        end
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
        ftype in [:bp, :bs] && @assert length(cutoff) == 2 "For :$(ftype) filter, cutoff must specify two frequencies."
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
    fprototype in [:chebyshev1, :elliptic] && (!isnothing(rp) && (rp = Observable(float(rp))))
    fprototype in [:chebyshev2, :elliptic] && !isnothing(rs) && (rs = Observable(float(rs)))
    fprototype in [:firls, :remez, :iirnotch] && @assert !isnothing(bw) "bw must be specified."
    !isnothing(bw) && (bw = Observable(float(bw)))

    # prepare plot
    GLMakie.activate!(title = "plot_filter()")
    plot_size = (1400, 900)
    p = GLMakie.Figure(; size = plot_size)
    grid = p[4, 1] = GridLayout()

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic]
# cutoff
# order
# rp
# rs

        if ftype in [:hp, :lp]
            _ = Label(
                    grid[1, 1],
                    "Cutoff [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_cutoff = Slider(
                            grid[1, 2],
                            range = 0.5:0.1:(nqf - 0.1),
                            startvalue = cutoff[],
                            horizontal = true,
                        )
            on(sl_cutoff.value) do val
                cutoff[] = round(val, digits=1)
                notify(cutoff)
            end

            _ = Label(
                    grid[2, 1],
                    "Order [taps]",
                    fontsize = 15,
                    halign = :right,
                )
            if ftype === :lp
                sl_order = Slider(
                                grid[2, 2],
                                range = 1:1:1000,
                                startvalue = order[],
                                horizontal = true,
                            )
            elseif ftype === :hp
                sl_order = Slider(
                                grid[2, 2],
                                range = 1:2:1001,
                                startvalue = order[],
                                horizontal = true,
                            )
            end
            on(sl_order.value) do val
                order[] = val
                notify(order)
            end

            if isa(rp, Observable{Float64})
                _ = Label(
                        grid[3, 1],
                        "RP [dB]",
                        fontsize = 15,
                        halign = :right,
                    )
                sl_rp = Slider(
                            grid[3, 2],
                            # range = fprototype === :elliptic ? (0.001:0.001:0.01) : (0.5:0.5:10),
                            range = 0.1:0.1:rs[],
                            startvalue = rp[],
                            horizontal = true,
                        )
                on(sl_rp.value) do val
                    rp[] = round(val, digits=1)
                    notify(rp)
                end
            end

            if isa(rs, Observable{Float64})
                _ = Label(
                        grid[fprototype === :chebyshev2 ? 3 : 4, 1],
                        "RS [dB]",
                        fontsize = 15,
                        halign = :right,
                    )
                sl_rs = Slider(
                            grid[fprototype === :chebyshev2 ? 3 : 4, 2],
                            range = 1:1:100,
                            startvalue = rs[],
                            horizontal = true,
                        )
                on(sl_rs.value) do val
                    rs[] = round(val, digits=1)
                    sl_rp.range = 0.1:0.1:(rs[] - 0.1)
                    notify(rs)
                end
            end

        elseif ftype in [:bp, :bs]

            _ = Label(
                    grid[1, 1],
                    "Cutoff [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_cutoff = IntervalSlider(
                                    grid[1, 2],
                                    range = 0.1:0.1:(nqf - 0.1),
                                    startvalues = cutoff[],
                                    horizontal = true,
                                )
            on(sl_cutoff.interval) do val
                cutoff[] = round.(val, digits=1)
                if cutoff[][1] == cutoff[][2]
                    cutoff[] = (cutoff[][1], cutoff[][1] + 0.1)
                elseif cutoff[][1] > cutoff[][2]
                    cutoff[] = (cutoff[][2], cutoff[][1])
                end
                notify(cutoff)
            end

            _ = Label(
                    grid[2, 1],
                    "Order [taps]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_order = Slider(
                            grid[2, 2],
                            range = 1:2:1001,
                            startvalue = order[],
                            horizontal = true,
                        )
            on(sl_order.value) do val
                order[] = val
                notify(order)
            end

            if isa(rp, Observable{Float64})
                _ = Label(
                        grid[3, 1],
                        "RP [dB]",
                        fontsize = 15,
                        halign = :right,
                    )
                sl_rp = Slider(
                            grid[3, 2],
                            # range = fprototype === :elliptic ? (0.001:0.001:0.01) : (0.5:0.5:10),
                            range = 0.1:0.1:rs[],
                            startvalue = rp[],
                            horizontal = true,
                        )
                on(sl_rp.value) do val
                    rp[] = round(val, digits=1)
                    notify(rp)
                end
            end

            if isa(rs, Observable{Float64})
                _ = Label(
                        grid[fprototype === :chebyshev2 ? 3 : 4, 1],
                        "RS [dB]",
                        fontsize = 15,
                        halign = :right,
                    )
                sl_rs = Slider(
                            grid[fprototype === :chebyshev2 ? 3 : 4, 2],
                            range = 1:1:100,
                            startvalue = rs[],
                            horizontal = true,
                        )
                on(sl_rs.value) do val
                    rs[] = round(val, digits=1)
                    sl_rp.range = 0.1:0.1:(rs[] - 0.1)
                    notify(rs)
                end
            end
        end

    elseif fprototype in [:remez]

        if ftype in [:hp, :lp]
            _ = Label(
                    grid[1, 1],
                    "Cutoff [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_cutoff = Slider(
                            grid[1, 2],
                            range = 0.5:0.1:(nqf - 0.1),
                            startvalue = cutoff[],
                            horizontal = true,
                        )
            on(sl_cutoff.value) do val
                cutoff[] = round(val, digits=1)
                notify(cutoff)
            end

            _ = Label(
                    grid[2, 1],
                    "Order [taps]",
                    fontsize = 15,
                    halign = :right,
                )
            if ftype === :lp
                sl_order = Slider(
                                grid[2, 2],
                                range = 1:1:1000,
                                startvalue = order[],
                                horizontal = true,
                            )
            elseif ftype === :hp
                sl_order = Slider(
                                grid[2, 2],
                                range = 1:2:1001,
                                startvalue = order[],
                                horizontal = true,
                            )
            end
            on(sl_order.value) do val
                order[] = val
                notify(order)
            end

            _ = Label(
                    grid[3, 1],
                    "Band width [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_bw = Slider(
                        grid[3, 2],
                        range = cutoff[][1] > 10 ? (0.1:0.1:10) : (0.1:0.1:(cutoff[][1] - 0.1)),
                        startvalue = bw[],
                        horizontal = true,
                    )
            on(sl_bw.value) do val
                bw[] = round(val, digits=1)
                notify(bw)
            end

        elseif ftype in [:bp, :bs]

            _ = Label(
                    grid[1, 1],
                    "Cutoff [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_cutoff = IntervalSlider(
                                    grid[1, 2],
                                    range = 0.1:0.1:(nqf - 0.1),
                                    startvalues = cutoff[],
                                    horizontal = true,
                                )
            on(sl_cutoff.interval) do val
                cutoff[] = round.(val, digits=1)
                if cutoff[][1] == cutoff[][2]
                    cutoff[] = (cutoff[][1], cutoff[][1] + 0.1)
                elseif cutoff[][1] > cutoff[][2]
                    cutoff[] = (cutoff[][2], cutoff[][1])
                end
                notify(cutoff)
            end

            _ = Label(
                    grid[2, 1],
                    "Order [taps]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_order = Slider(
                            grid[2, 2],
                            range = 1:2:1001,
                            startvalue = order[],
                            horizontal = true,
                        )
            on(sl_order.value) do val
                order[] = val
                notify(order)
            end

            _ = Label(
                    grid[3, 1],
                    "Band width [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_bw = Slider(
                        grid[3, 2],
                        range = cutoff[][1] > 10 ? (0.1:0.1:10) : (0.1:0.1:(cutoff[][1] - 0.1)),
                        startvalue = bw[],
                        horizontal = true,
                    )
            on(sl_bw.value) do val
                bw[] = round(val, digits=1)
                notify(bw)
            end

        end

    elseif fprototype in [:fir]

        if ftype in [:hp, :lp]
            _ = Label(
                    grid[1, 1],
                    "Cutoff [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_cutoff = Slider(
                            grid[1, 2],
                            range = 0.5:0.1:(nqf - 0.1),
                            startvalue = cutoff[],
                            horizontal = true,
                        )
            on(sl_cutoff.value) do val
                cutoff[] = round(val, digits=1)
                notify(cutoff)
            end

            if isnothing(w)
                _ = Label(
                        grid[2, 1],
                        "Order [taps]",
                        fontsize = 15,
                        halign = :right,
                    )
                if ftype === :lp
                    sl_order = Slider(
                                    grid[2, 2],
                                    range = 1:1:1000,
                                    startvalue = order[],
                                    horizontal = true,
                                )
                elseif ftype === :hp
                    sl_order = Slider(
                                    grid[2, 2],
                                    range = 1:2:1001,
                                    startvalue = order[],
                                    horizontal = true,
                                )
                end
                on(sl_order.value) do val
                    order[] = val
                    notify(order)
                end
            end

        elseif ftype in [:bp, :bs]

            _ = Label(
                    grid[1, 1],
                    "Cutoff [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_cutoff = IntervalSlider(
                                    grid[1, 2],
                                    range = 0.1:0.1:(nqf - 0.1),
                                    startvalues = cutoff[],
                                    horizontal = true,
                                )
            on(sl_cutoff.interval) do val
                cutoff[] = round.(val, digits=1)
                if cutoff[][1] == cutoff[][2]
                    cutoff[] = (cutoff[][1], cutoff[][1] + 0.1)
                elseif cutoff[][1] > cutoff[][2]
                    cutoff[] = (cutoff[][2], cutoff[][1])
                end
                notify(cutoff)
            end

            if isnothing(w)
                _ = Label(
                        grid[2, 1],
                        "Order [taps]",
                        fontsize = 15,
                        halign = :right,
                    )
                sl_order = Slider(
                                grid[2, 2],
                                range = 1:2:1001,
                                startvalue = order[],
                                horizontal = true,
                            )
                on(sl_order.value) do val
                    order[] = val
                    notify(order)
                end
            end

        end

    elseif fprototype in [:firls]

        if ftype in [:hp, :lp]
            _ = Label(
                    grid[1, 1],
                    "Cutoff [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_cutoff = Slider(
                            grid[1, 2],
                            range = 0.5:0.1:(nqf - 0.1),
                            startvalue = cutoff[],
                            horizontal = true,
                        )
            on(sl_cutoff.value) do val
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
                    "Band width [Hz]",
                    fontsize = 15,
                    halign = :right,
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

            if isnothing(w)
                _ = Label(
                        grid[3, 1],
                        "Order [taps]",
                        fontsize = 15,
                        halign = :right,
                    )
                sl_order = Slider(
                                grid[3, 2],
                                range = 1:1:1000,
                                startvalue = order,
                                horizontal = true,
                            )
                on(sl_order.value) do val
                    order[] = val
                    notify(order)
                end
            end

        elseif ftype in [:bp, :bs]

            _ = Label(
                    grid[1, 1],
                    "Cutoff [Hz]",
                    fontsize = 15,
                    halign = :right,
                )
            sl_cutoff = IntervalSlider(
                                    grid[1, 2],
                                    range = 0.1:0.1:(nqf - 0.1),
                                    startvalues = cutoff,
                                    horizontal = true,
                                )
            on(sl_cutoff.values) do val
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
                    "Band width [Hz]",
                    fontsize = 15,
                    halign = :right,
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

            if isnothing(w)
                _ = Label(
                        grid[3, 1],
                        "Order [taps]",
                        fontsize = 15,
                        halign = :right,
                    )
                sl_order = Slider(
                                grid[3, 2],
                                range = 1:1:1000,
                                startvalue = order,
                                horizontal = true,
                            )
                on(sl_order.value) do val
                    order[] = val
                    notify(order)
                end
            end

        end

    elseif fprototype in [:iirnotch]

        _ = Label(
                grid[1, 1],
                "Cutoff [Hz]",
                fontsize = 15,
                halign = :right,
            )
        sl_cutoff = Slider(
                        grid[1, 2],
                        range = 0.5:0.1:(nqf - 0.1),
                        startvalue = cutoff[],
                        horizontal = true,
                    )
        on(sl_cutoff.value) do val
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
                "Band width [Hz]",
                fontsize = 15,
                halign = :right,
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
            if fprototype in [:chebyshev1, :chebyshev2, :elliptic]
                title = @lift("Filter: $(fname), type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff, digits=1)) Hz, order: $($order), RP: $($rp) dB, RS: $($rs) dB\n\nFrequency response")
            else
                title = @lift("Filter: $(fname), type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff, digits=1)) Hz, order: $($order)\n\nFrequency response")
            end
        else
            fname = "IIR notch"
            title = @lift("Filter: $(fname), cutoff: $(round.($cutoff, digits=1)) Hz, transition band width: $(round($bw, digits=1)) Hz\n\nFrequency response")
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

        fresp = lift(_fir_response, flt)
        # convert to dB
        H = @lift(amp2db.(abs.($fresp)))
        f = range(0; stop = pi, length = 1024)
        # convert rad/sample to Hz
        f = f .* fs / 2 / pi

        if fprototype === :fir
            title = @lift("Filter: FIR, type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff, digits=1)) Hz, order: $($order)\n\nFrequency response")
        elseif fprototype === :firls
            title = @lift("Filter: FIR (LS), type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff, digits=1)) Hz, transition band width: $($bw) Hz, order: $($order)\n\nFrequency response")
        elseif fprototype === :remez
            title = @lift("Filter: Remez, type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff, digits=1)) Hz, transition band width: $($bw) Hz, order: $($order)\n\nFrequency response")
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

        fresp = lift(_fir_response, flt)
        # phi = _fir_response(flt[], f)
        # convert to dB
        phi = @lift(rad2deg.(-atan.(imag($fresp), real($fresp))))
        # convert rad/sample to Hz
        f = range(0; stop = pi, length = 1024)
        f = f .* fs / 2 / pi

        ax2 = GLMakie.Axis(
                        p[2, 1];
                        xlabel = "Frequency [Hz]",
                        ylabel = "Phase\n[deg]",
                        title = "Phase response",
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
                )


        fresp = lift(_fir_response, flt)
        tau = @lift(-derivative(rad2deg.(-atan.(imag($fresp), real($fresp)))))

        ax3 = GLMakie.Axis(
                        p[3, 1];
                        xlabel = "Frequency [Hz]",
                        ylabel = "Group delay\n[samples]",
                        title = "Group delay",
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

        if isa(bw, Observable{Float64})
            if ftype === :lp
                f_pass = @lift($cutoff - ($bw / 2))
                f_stop = @lift($cutoff + ($bw / 2))
            elseif ftype === :hp
                f_pass = @lift($cutoff + ($bw / 2))
                f_stop = @lift($cutoff - ($bw / 2))
            end
            GLMakie.vlines!(
                        ax1,
                        f_pass,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax2,
                        f_pass,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax3,
                        f_pass,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax1,
                        f_stop,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax2,
                        f_stop,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax3,
                        f_stop,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
        end
    else
        c1 = lift(cutoff) do val
            val[1]
        end
        c2 = lift(cutoff) do val
            val[2]
        end
        GLMakie.vlines!(
                    ax1,
                    c1,
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :red,
                )
        GLMakie.vlines!(
                    ax1,
                    c2,
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :green,
                )
        GLMakie.vlines!(
                    ax2,
                    c1,
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :red,
                )
        GLMakie.vlines!(
                    ax2,
                    c2,
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :green,
                )
        GLMakie.vlines!(
                    ax3,
                    c1,
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :red,
                )
        GLMakie.vlines!(
                    ax3,
                    c2,
                    linestyle = :dash,
                    linewidth = 1,
                    color = mono ? :black : :green,
                )

        if isa(bw, Observable{Float64})
            if ftype === :bp
                f_pass = @lift($cutoff[2] + ($bw / 2))
                f_stop = @lift($cutoff[1] - ($bw / 2))
            elseif ftype === :bs
                f_pass = @lift($cutoff[1] - ($bw / 2))
                f_stop = @lift($cutoff[2] + ($bw / 2))
            end
            GLMakie.vlines!(
                        ax1,
                        f_pass,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax2,
                        f_pass,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax3,
                        f_pass,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax1,
                        f_stop,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax2,
                        f_stop,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
            GLMakie.vlines!(
                        ax3,
                        f_stop,
                        linestyle = :dash,
                        linewidth = 0.25,
                        color = :black,
                    )
        end

    end

    if gui
        wait(display(p))
        NeuroAnalyzer.verbose = v
        return flt[]
    else
        NeuroAnalyzer.verbose = v
        return p
    end

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
  - `flim::Tuple{Real, Real}=(0, sr(obj) / 2): frequency limit
  - `mono::Bool=false`: use color or gray palette
  - `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

  - `p::GLMakie.Figure`
  - `f::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: if `gui=true`
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
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    mono::Bool = false,
    gui::Bool = true,
)::GLMakie.Figure

    return plot_filter(;
                    fs = sr(obj),
                    fprototype = fprototype,
                    ftype = ftype,
                    cutoff = cutoff,
                    order = order,
                    rp = rp,
                    rs = rs,
                    bw = bw,
                    w = w,
                    flim = flim,
                    mono = mono,
                    gui = gui,
                )

    return p

end
