export plot_filter

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Create a standard filter-response Axis with shared lock/style settings.
function _filter_axis(fig, pos, title, ylabel, flim)
    ax = GLMakie.Axis(
        fig[pos...];
        xlabel             = "Frequency [Hz]",
        ylabel             = ylabel,
        title              = title,
        xticks             = length(flim[1]:0.1:flim[2]) > 20 ? LinearTicks(10) : LinearTicks(20),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, flim)
    _style_axis!(ax)
    return ax
end

# Add a labeled cutoff Slider (single or interval) to a grid layout row.
function _add_cutoff_slider!(grid, row, cutoff, nqf, is_interval)
    Label(grid[row, 1], "Cutoff [Hz]"; fontsize = 15, halign = :right)
    if is_interval
        sl = IntervalSlider(
            grid[row, 2];
            range       = 0.1:0.1:(nqf - 0.1),
            startvalues = cutoff[],
            horizontal  = true,
        )
        on(sl.interval) do val
            cutoff[] = round.(val; digits = 1)
            cutoff[][1] == cutoff[][2] && (cutoff[] = (cutoff[][1], cutoff[][1] + 0.1))
            cutoff[][1] > cutoff[][2] && (cutoff[] = (cutoff[][2], cutoff[][1]))
            return notify(cutoff)
        end
    else
        sl = Slider(
            grid[row, 2];
            range      = 0.5:0.1:(nqf - 0.1),
            startvalue = cutoff[],
            horizontal = true,
        )
        on(sl.value) do val
            cutoff[] = round(val; digits = 1)
            return notify(cutoff)
        end
    end
    return sl
end

# Add a labeled order Slider to a grid layout row.
function _add_order_slider!(grid, row, order, ftype)
    Label(grid[row, 1], "Order [taps]"; fontsize = 15, halign = :right)
    rng = (ftype === :lp) ? (1:1:1000) : (1:2:1001)
    sl  = Slider(
    grid[row, 2];
    range      = rng,
    startvalue = order[],   # FIX: was `order` (Observable) in :firls branch
    horizontal = true
)
    on(sl.value) do val
        order[] = val
        return notify(order)
    end
    return sl
end

# Add a labeled bandwidth Slider to a grid layout row, with dynamic range linked to cutoff.
function _add_bw_slider!(grid, row, bw, cutoff_ref)
    Label(grid[row, 1], "Band width [Hz]"; fontsize = 15, halign = :right)
    sl = Slider(
        grid[row, 2];
        range      = cutoff_ref > 10 ? (0.1:0.1:10) : (0.1:0.1:(cutoff_ref - 0.1)),
        startvalue = bw[],
        horizontal = true,
    )
    on(sl.value) do val
        bw[] = round(val; digits = 1)
        return notify(bw)
    end
    return sl
end

# Draw reactive cutoff vlines on three axes, handling single vs. two-frequency cutoff.
function _draw_cutoff_vlines!(ax1, ax2, ax3, cutoff, bw, ftype, mono)
    vl_kwargs = (linestyle = :dash, linewidth = 1)
    thin_kwargs = (linestyle = :dash, linewidth = 0.25, color = :black)

    if length(cutoff[]) == 1
        color = mono ? :black : :red
        for ax in (ax1, ax2, ax3)
            GLMakie.vlines!(ax, cutoff; vl_kwargs..., color = color)
        end
        if isa(bw, Observable{Float64})
            f_pass = ftype === :lp ?
                     @lift($cutoff - ($bw / 2)) : @lift($cutoff + ($bw / 2))
            f_stop = ftype === :lp ?
                     @lift($cutoff + ($bw / 2)) : @lift($cutoff - ($bw / 2))
            for ax in (ax1, ax2, ax3)
                GLMakie.vlines!(ax, f_pass; thin_kwargs...)
                GLMakie.vlines!(ax, f_stop; thin_kwargs...)
            end
        end
    else
        c1 = @lift($cutoff[1])
        c2 = @lift($cutoff[2])
        for ax in (ax1, ax2, ax3)
            GLMakie.vlines!(ax, c1; vl_kwargs..., color = mono ? :black : :red)
            GLMakie.vlines!(ax, c2; vl_kwargs..., color = mono ? :black : :green)
        end
        if isa(bw, Observable{Float64})
            f_pass =
                ftype === :bp ?
                @lift($cutoff[2] + ($bw / 2)) : @lift($cutoff[1] - ($bw / 2))
            f_stop =
                ftype === :bp ?
                @lift($cutoff[1] - ($bw / 2)) : @lift($cutoff[2] + ($bw / 2))
            for ax in (ax1, ax2, ax3)
                GLMakie.vlines!(ax, f_pass; thin_kwargs...)
                GLMakie.vlines!(ax, f_stop; thin_kwargs...)
            end
        end
    end
end

# ---------------------------------------------------------------------------

"""
    plot_filter(; <keyword arguments>)

Plot filter response with interactive controls for various filter types.

# Arguments

- `fs::Int64`: sampling rate in Hz; must be ≥ 1
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
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz
    - for `:lp`/`:hp`: single frequency
    - for `:bp`/`:bs`: frequency range (f1, f2)
- `order::Union{Nothing, Int64}=nothing`: filter order (number of taps for FIR, filter order for IIR)
- `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.5 dB
- `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 20 dB
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
- `flim::Tuple{Real, Real} = (0, fs / 2)`: frequency limits
- `mono::Bool=false`: if `true`, use a monochrome palette
- `gui::Bool=true`: if `true`, keep window open and interactive

# Returns

- `GLMakie.Figure`: the plotted figure, if `gui = false`
- `Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: returns the filter object, if `gui=true`

# Notes

- For IIR filters (`:butterworth`, `:chebyshev1`, etc.), default ripple values are:
    - Passband ripple (`rp`): 0.5 dB
    - Stopband attenuation (`rs`): 20 dB
- For `:elliptic` filters, defaults are 0.5 dB and 40 dB respectively.
- For FIR filters, window length must be odd.
- Bandwidth (`bw`) is required for `:firls`, `:remez`, and `:iirnotch` filters.
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
)::Union{
    GLMakie.Figure,
    Vector{Float64},
    ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64},
    Biquad{:z, Float64},
}
    # validate
    _check_tuple(flim, (0, fs / 2), "flim")
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))

    # set verbose to false during calculations
    v = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false

    # check parameters
    try
        # Nyquist frequency
        nqf = div(fs, 2)
        nqf > flim[2] && (nqf = flim[2])

        _check_var(
            fprototype,
            [
                :fir,
                :firls,
                :remez,
                :butterworth,
                :chebyshev1,
                :chebyshev2,
                :elliptic,
                :iirnotch,
            ],
            "fprototype",
        )
        !isnothing(ftype) && _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")

        if fprototype === :fir
            (isnothing(order) && isnothing(w)) &&
                throw(ArgumentError("Either order or w must be specified."))
            if !isnothing(w)
                (ftype in [:hp, :bp, :bs] && mod(length(w), 2) != 0) ||
                    throw(ArgumentError("Length of w must be odd."))
                length(w) >= 1 || throw(ArgumentError("Length of w must be ≥ 1."))
            elseif !isnothing(order)
                (ftype in [:hp, :bp, :bs] && mod(order, 2) != 0) ||
                    throw(ArgumentError("order must be odd."))
            end
        end

        if fprototype in [:firls, :remez, :iirnotch]
            isnothing(bw) && throw(ArgumentError("bw must be specified."))
            bw > 0 || throw(ArgumentError("bw must be > 0."))
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
                    length(w) == 6 || throw(ArgumentError("Length of w must be 6."))
                else
                    w = ones(6)
                end
            elseif ftype in [:lp, :hp]
                if !isnothing(w)
                    length(w) == 4 || throw(ArgumentError("Length of w must be 4."))
                else
                    w = ones(4)
                end
            end
        end

        if fprototype in [:chebyshev1, :chebyshev2, :elliptic]
            if isnothing(rp)
                rp = 0.5
                _info("rp set at $rp dB.")
            end
            if isnothing(rs)
                rs = 20
                _info("rs set at $rs dB.")
            end
        end

        if fprototype in [:firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic]
            isnothing(order) && throw(ArgumentError("order must be specified."))
            isnothing(ftype) && throw(ArgumentError("ftype must be specified."))
        end

        if fprototype === :iirnotch
            isnothing(ftype) || _info("For :iirnotch filter ftype is ignored")
            isnothing(order) || _info("For :iirnotch filter order is ignored")
            length(cutoff) == 1 || throw(
                ArgumentError(
                    "For :iirnotch filter cutoff must contain only one frequency.",
                ),
            )
        end

        if fprototype in [:fir, :butterworth, :chebyshev1, :chebyshev2, :elliptic]
            (ftype in [:lp, :hp] && length(cutoff) == 1) ||
                throw(
                    ArgumentError(
                        "For :$(ftype) filter, cutoff must specify only one frequency.",
                    ),
                )
            (ftype in [:bp, :bs] && length(cutoff) == 2) ||
                throw(
                    ArgumentError(
                        "For :$(ftype) filter, cutoff must specify two frequencies.",
                    ),
                )
        end

        if length(cutoff) == 1
            cutoff > 0 || throw(ArgumentError("cutoff must be > 0 Hz."))
            cutoff < nqf || throw(ArgumentError("cutoff must be < $nqf Hz."))
        else
            _check_tuple(cutoff, (0, nqf), "cutoff")
        end

        # wrap parameters in Observables for reactive updates
        cutoff = Observable(float.(cutoff))
        order  = Observable(order)
        fprototype in [:chebyshev1, :elliptic] &&
            !isnothing(rp) && (rp = Observable(float(rp)))
        fprototype in [:chebyshev2, :elliptic] &&
            !isnothing(rs) && (rs = Observable(float(rs)))
        !isnothing(bw) && (bw = Observable(float(bw)))

        # prepare plot
        GLMakie.activate!(; title = "plot_filter()")
        fig = GLMakie.Figure(; size = gui ? (1200, 900) : (1200, 800))

        # GUI sliders
        if gui
            grid        = fig[4, 1] = GridLayout()
            is_interval = !isnothing(ftype) && ftype in [:bp, :bs]

            if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic]
                sl_cutoff = _add_cutoff_slider!(grid, 1, cutoff, nqf, is_interval)
                sl_order  = _add_order_slider!(grid, 2, order, isnothing(ftype) ? :lp : ftype)

                if isa(rp, Observable{Float64})
                    Label(grid[3, 1], "RP [dB]"; fontsize = 15, halign = :right)
                    sl_rp = Slider(
                        grid[3, 2];
                        range      = 0.1:0.1:(isa(rs, Observable{Float64}) ? rs[] : 10.0),
                        startvalue = rp[],
                        horizontal = true,
                    )
                    on(sl_rp.value) do val
                        rp[] = round(val; digits = 1);
                        return notify(rp)
                    end
                end

                if isa(rs, Observable{Float64})
                    rs_row = fprototype === :chebyshev2 ? 3 : 4
                    Label(grid[rs_row, 1], "RS [dB]"; fontsize = 15, halign = :right)
                    sl_rs = Slider(
                        grid[rs_row, 2];
                        range      = 1:1:100,
                        startvalue = rs[],
                        horizontal = true,
                    )
                    on(sl_rs.value) do val
                        rs[] = round(val; digits = 1)
                        isa(rp, Observable{Float64}) && (sl_rp.range = 0.1:0.1:(rs[] - 0.1))
                        return notify(rs)
                    end
                end

            elseif fprototype === :remez
                sl_cutoff = _add_cutoff_slider!(grid, 1, cutoff, nqf, is_interval)
                sl_order  = _add_order_slider!(grid, 2, order, isnothing(ftype) ? :lp : ftype)
                _add_bw_slider!(grid, 3, bw, cutoff[][is_interval ? 1 : 1])

            elseif fprototype === :fir
                sl_cutoff = _add_cutoff_slider!(grid, 1, cutoff, nqf, is_interval)
                isnothing(w) &&
                    _add_order_slider!(grid, 2, order, isnothing(ftype) ? :lp : ftype)

            elseif fprototype === :firls
                sl_cutoff = _add_cutoff_slider!(grid, 1, cutoff, nqf, is_interval)
                on(sl_cutoff.value) do val  # update bw range when cutoff changes
                    c = is_interval ? cutoff[][1] : cutoff[]
                    if !isnothing(bw)
                        if c > 10
                            sl_bw.range = 0.1:0.1:10
                        else
                            bw[] >= c && (bw[] = c - 0.1; set_close_to!(sl_bw, bw[]))
                            sl_bw.range = 0.1:0.1:(c - 0.1)
                        end
                    end
                end
                sl_bw = _add_bw_slider!(grid, 2, bw, is_interval ? cutoff[][1] : cutoff[])
                isnothing(w) &&
                    _add_order_slider!(grid, 3, order, isnothing(ftype) ? :lp : ftype)

            elseif fprototype === :iirnotch
                sl_cutoff = _add_cutoff_slider!(grid, 1, cutoff, nqf, false)
                on(sl_cutoff.value) do val
                    c = cutoff[]
                    if !isnothing(bw)
                        if c > 10
                            sl_bw.range = 0.1:0.1:10
                        else
                            bw[] >= c && (bw[] = c - 0.1; set_close_to!(sl_bw, bw[]))
                            sl_bw.range = 0.1:0.1:(c - 0.1)
                        end
                    end
                end
                sl_bw = _add_bw_slider!(grid, 2, bw, cutoff[])
            end
        end

        # create filter observable
        flt = @lift(
            filter_create(
                fprototype = fprototype,
                ftype      = ftype,
                cutoff     = $cutoff,
                fs         = fs,
                order      = $order,
                rp         = !isnothing(rp) ? $rp : nothing,
                rs         = !isnothing(rs) ? $rs : nothing,
                bw         = !isnothing(bw) ? $bw : nothing,
                w          = w,
            )
        )

        # draw frequency, phase, and group-delay response plots
        if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch]
            fresp = lift(DSP.freqresp, flt)
            H     = @lift(real.(20 * log10.(abs.($fresp[1]))))
            f_hz  = @lift(round.($fresp[2] .* fs / 2 / pi; digits = 1))

            if fprototype !== :iirnotch
                fname = titlecase(String(fprototype))
                title1 = if fprototype in [:chebyshev1, :chebyshev2, :elliptic]
                    @lift(
                        "Filter: $(fname), type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff; digits=1)) Hz, order: $($order), RP: $($rp) dB, RS: $($rs) dB\n\nFrequency response"
                    )
                else
                    @lift(
                        "Filter: $(fname), type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff; digits=1)) Hz, order: $($order)\n\nFrequency response"
                    )
                end
            else
                title1 = @lift(
                    "Filter: IIR notch, cutoff: $(round.($cutoff; digits=1)) Hz, bw: $(round($bw; digits=1)) Hz\n\nFrequency response"
                )
            end

            ax1 = _filter_axis(fig, (1, 1), title1, "Magnitude [dB]", flim)
            GLMakie.ylims!(ax1, (-100, 20))
            GLMakie.lines!(ax1, f_hz, H; color = mono ? :black : :blue)

            phresp = lift(DSP.phaseresp, flt)
            phi    = @lift($phresp[1])
            f_ph   = @lift(round.($phresp[2] .* fs / 2 / pi; digits = 1))
            tau    = @lift(-derivative(rad2deg.($phresp[1])))

            ax2 = _filter_axis(fig, (2, 1), "Phase response", "Phase [rad]", flim)
            GLMakie.lines!(
                ax2,
                f_ph,
                phi;
                color = mono ? :black : :blue,
                nan_color = mono ? :black : :blue,
            )

            ax3 = _filter_axis(fig, (3, 1), "Group delay", "Group delay [samples]", flim)
            GLMakie.lines!(ax3, f_ph, tau; color = mono ? :black : :blue)

        else  # FIR family
            fresp = lift(_fir_response, flt)
            H     = @lift(amp2db.(abs.($fresp)))
            phi   = @lift(rad2deg.(-atan.(imag($fresp), real($fresp))))
            tau   = @lift(-derivative(rad2deg.(-atan.(imag($fresp), real($fresp)))))
            f_fir = range(0; stop = pi, length = 1024) .* fs / 2 / pi

            title1 = if fprototype === :fir
                @lift(
                    "Filter: FIR, type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff; digits=1)) Hz, order: $($order)\n\nFrequency response"
                )
            elseif fprototype === :firls
                @lift(
                    "Filter: FIR (LS), type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff; digits=1)) Hz, bw: $($bw) Hz, order: $($order)\n\nFrequency response"
                )
            else  # :remez
                @lift(
                    "Filter: Remez, type: $(uppercase(String(ftype))), cutoff: $(round.($cutoff; digits=1)) Hz, bw: $($bw) Hz, order: $($order)\n\nFrequency response"
                )
            end

            ax1 = _filter_axis(fig, (1, 1), title1, "Magnitude [dB]", flim)
            GLMakie.ylims!(ax1, (-100, 20))
            GLMakie.lines!(ax1, f_fir, H; color = mono ? :black : :blue)

            ax2 = _filter_axis(fig, (2, 1), "Phase response", "Phase [deg]", flim)
            GLMakie.lines!(ax2, f_fir, phi; color = mono ? :black : :blue)

            ax3 = _filter_axis(fig, (3, 1), "Group delay", "Group delay [samples]", flim)
            GLMakie.lines!(ax3, f_fir, tau; color = mono ? :black : :blue)
        end

        # draw cutoff indicator lines
        _draw_cutoff_vlines!(ax1, ax2, ax3, cutoff, bw, ftype, mono)

        if gui
            wait(display(fig))
            NeuroAnalyzer.verbose = v
            return flt[]
        else
            NeuroAnalyzer.verbose = v
            return fig
        end

    catch
        NeuroAnalyzer.verbose = v
        rethrow()
    end
end

"""
    plot_filter(obj, <keyword arguments>)

Plot the frequency response of a digital filter with customizable visualization options.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (used only for sampling rate information)
- `n::Int64`: signal length in samples for frequency response calculation
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
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz
    - for `:lp`/`:hp`: single frequency
    - for `:bp`/`:bs`: frequency range (f1, f2)
- `order::Union{Nothing, Int64}=nothing`: filter order (number of taps for FIR, filter order for IIR)
- `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2): frequency limit
- `mono::Bool=false`: if `true`, use a monochrome palette
- `gui::Bool=true`: if `true`, keep window open and interactive

# Returns

- `GLMakie.Figure`: the plotted figure, if `gui=true`
- `Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: the filter object, if `gui=false`
"""
function plot_filter(
    obj::NeuroAnalyzer.NEURO;
    fprototype::Symbol,
    ftype::Union{Nothing, Symbol} = nothing,
    cutoff::Union{Real, Tuple{Real, Real}},
    order::Union{Nothing, Int64} = nothing,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    bw::Union{Nothing, Real} = nothing,
    w::Union{Nothing, AbstractVector} = nothing,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    mono::Bool = false,
    gui::Bool = true,
)::Union{
    GLMakie.Figure,
    Vector{Float64},
    ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64},
    Biquad{:z, Float64},
}
    return plot_filter(;
        fs         = sr(obj),
        fprototype = fprototype,
        ftype      = ftype,
        cutoff     = cutoff,
        order      = order,
        rp         = rp,
        rs         = rs,
        bw         = bw,
        w          = w,
        flim       = flim,
        mono       = mono,
        gui        = gui,
    )
end
