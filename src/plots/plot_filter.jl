export plot_filter

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# complex frequency response at frequencies f (Hz)
function _filter_freqz(flt, f::AbstractVector, fs::Real)
    ω = 2π .* f ./ fs
    return flt isa Vector{Float64} ? freqresp(PolynomialRatio(flt, [1.0]), ω) :
           freqresp(flt, ω)
end

# group delay in samples: τ = -dφ/dω = -Im(H'/H); no phase unwrapping needed, NaN where |H| ≈ 0
function _group_delay(H::AbstractVector, f::AbstractVector, fs::Real)
    ω = 2π .* f ./ fs
    hmax = maximum(abs.(H))
    τ = fill(NaN, length(H))
    for i = 2:(length(H) - 1)
        abs(H[i]) > 1e-6 * hmax || continue
        τ[i] = -imag((H[i + 1] - H[i - 1]) / ((ω[i + 1] - ω[i - 1]) * H[i]))
    end
    return τ
end

# magnitude (dB), phase (deg) and group delay (samples) of the filter as applied in direction `dir`
function _filter_response(flt, f::AbstractVector, fs::Real, dir::Symbol)
    H = _filter_freqz(flt, f, fs)
    g = 20 .* log10.(max.(abs.(H), floatmin(Float64)))
    # two-pass: |H|², zero phase, zero delay
    dir === :twopass && return (2 .* g, zeros(length(f)), zeros(length(f)))
    phi = rad2deg.(DSP.unwrap(angle.(H)))
    tau = _group_delay(H, f, fs)
    # reverse pass: conj(H), i.e. negated phase and delay
    dir === :reverse && return (g, -phi, -tau)
    return (g, phi, tau)
end

# transition width the design actually has
function _bw_effective(
    flt,
    fprototype::Symbol,
    fs::Real,
    window::Symbol,
    custom_w::Bool,
    bw,
)
    if fprototype === :fir
        return custom_w ? nothing : _FIR_WINDOWS[window].k * fs / length(flt)
    elseif fprototype in (:firls, :remez, :iirnotch)
        return bw
    end
    return nothing
end

_fmt(x::Real) = string(round(x; digits = 2))
_fmt(x::Tuple) = "(" * Base.join(_fmt.(x), ", ") * ")"
_fmt(x::AbstractVector) = isempty(x) ? "not reached" : Base.join(_fmt.(x), ", ")

function _filter_title(
    flt;
    fprototype,
    ftype,
    cutoff,
    order,
    bw,
    bw_eff,
    rp,
    rs,
    fs,
    dir,
    window,
    custom_w,
    status,
)
    names = Dict(
        :fir => "FIR (window)",
        :firls => "FIR (least squares)",
        :remez => "FIR (Remez)",
        :butterworth => "Butterworth",
        :chebyshev1 => "Chebyshev I",
        :chebyshev2 => "Chebyshev II",
        :elliptic => "Elliptic",
        :iirnotch => "IIR notch",
    )
    s = "Filter: $(names[fprototype])"
    !isnothing(ftype) && (s *= ", type: $(uppercase(String(ftype)))")
    s *= ", cutoff: $(_fmt(cutoff)) Hz"
    if flt isa Vector{Float64}
        s *= ", taps: $(length(flt))"
        fprototype === :fir && (s *= ", window: $(custom_w ? "custom" : String(window))")
    elseif fprototype !== :iirnotch
        s *= ", order: $order"
    end
    !isnothing(bw) && (s *= ", bw: $(_fmt(bw)) Hz")
    !isnothing(bw_eff) && fprototype === :fir && (s *= " (effective ≈ $(_fmt(bw_eff)) Hz)")
    isnothing(bw) && !isnothing(bw_eff) && (s *= ", bw ≈ $(_fmt(bw_eff)) Hz")
    !isnothing(rp) && fprototype in (:chebyshev1, :elliptic) &&
        (s *= ", RP: $(_fmt(rp)) dB")
    !isnothing(rs) && fprototype in (:chebyshev2, :elliptic) &&
        (s *= ", RS: $(_fmt(rs)) dB")
    s *= ", fs: $fs Hz, $(dir === :twopass ? "two-pass (|H|²)" : String(dir))"

    # measured response
    r = filter_report(
        flt;
        fs = fs,
        dir = dir,
        ftype = fprototype === :iirnotch ? nothing : ftype,
        cutoff = cutoff,
        bw = bw_eff,
        n = 2^14,
        verbose = false,
    )
    m = "−3 dB: $(_fmt(r.crossings[-3.0])) Hz; −6 dB: $(_fmt(r.crossings[-6.0])) Hz; 0 Hz: $(_fmt(r.gain_dc)) dB"
    !isnothing(r.stop_att) && (m *= "; stop-band ≥ $(_fmt(r.stop_att)) dB")
    !isnothing(r.pass_dev) && (m *= "; pass-band ≤ $(_fmt(r.pass_dev)) dB")
    !isnothing(r.kernel_s) && (m *= "; kernel: $(_fmt(r.kernel_s)) s")
    s *= "\n" * m
    !isempty(status) && (s *= "\n⚠ $status — showing last valid filter")
    return s * "\n\nFrequency response"
end

# standard filter-response axis
function _filter_axis(fig, pos, title, ylabel, flim)
    ax = GLMakie.Axis(
        fig[pos...];
        xlabel = "Frequency [Hz]",
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, flim)
    _style_axis!(ax)
    return ax
end

# labeled slider writing to an Observable; current value shown on the right
function _add_slider!(grid, row, label, rng, start, obs)
    Label(grid[row, 1], label; fontsize = 15, halign = :right)
    sl = Slider(grid[row, 2]; range = rng, startvalue = start, horizontal = true)
    Label(grid[row, 3], lift(x -> string(x), sl.value); fontsize = 15, halign = :left)
    on(sl.value) do val
        return obs[] = Float64(val)
    end
    return sl
end

function _add_cutoff_slider!(grid, row, cutoff, nqf, is_interval)
    rng = 0.1:0.1:(floor(10 * nqf) / 10 - 0.1)
    if is_interval
        Label(grid[row, 1], "Cutoff [Hz]"; fontsize = 15, halign = :right)
        sl = IntervalSlider(
            grid[row, 2];
            range = rng,
            startvalues = cutoff[],
            horizontal = true,
        )
        Label(
            grid[row, 3],
            lift(x -> _fmt(Float64.(x)), sl.interval);
            fontsize = 15,
            halign = :left,
        )
        on(sl.interval) do val
            a, b = Float64.(val)
            a == b && (b = a + 0.1)
            return cutoff[] = (min(a, b), max(a, b))
        end
        return sl
    else
        return _add_slider!(grid, row, "Cutoff [Hz]", rng, cutoff[], cutoff)
    end
end

# cutoff (dashed) and transition-band edges (dotted) on all axes
function _draw_cutoff_vlines!(axes, cutoff, bw_eff, mono)
    cs = lift(c -> collect(Float64, c isa Tuple ? c : (c,)), cutoff)
    edges = lift(cutoff, bw_eff) do c, b
        isnothing(b) && return [NaN]
        cc = collect(Float64, c isa Tuple ? c : (c,))
        return vcat(cc .- b / 2, cc .+ b / 2)
    end
    for ax in axes
        GLMakie.vlines!(
            ax,
            cs;
            linestyle = :dash,
            linewidth = 1,
            color = mono ? :black : :red,
        )
        GLMakie.vlines!(ax, edges; linestyle = :dot, linewidth = 0.75, color = :black)
    end
    return nothing
end

# ---------------------------------------------------------------------------

"""
    plot_filter(; <keyword arguments>)

Plot the frequency, phase and group-delay response of a filter, with interactive controls.

The filter is designed with `filter_create` (same validation, automatic FIR order from `bw`), and the response is shown as applied in direction `dir` (for `:twopass`: |H|², zero phase, zero delay). The title reports the design parameters and the measured response (see `filter_report`).

# Arguments

- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `fprototype::Symbol`: filter prototype (`:fir`, `:firls`, `:remez`, `:butterworth`, `:chebyshev1`, `:chebyshev2`, `:elliptic`, `:iirnotch`)
- `ftype::Union{Nothing, Symbol}=nothing`: filter type (`:lp`, `:hp`, `:bp`, `:bs`); ignored for `:iirnotch`
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz
    - for `:lp`/`:hp`: single frequency
    - for `:bp`/`:bs`: frequency range (f1, f2)
- `order::Union{Nothing, Int64}=nothing`: filter order (number of taps for FIR, filter order for IIR); FIR: if `nothing`, calculated from `bw` and the GUI shows a `bw` slider instead of an order slider
- `rp::Union{Nothing, Real}=nothing`: pass-band ripple in dB (default 0.5 dB)
- `rs::Union{Nothing, Real}=nothing`: stop-band attenuation in dB (default 20 dB for IIR)
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz
- `w::Union{Nothing, AbstractVector}=nothing`: window vector for `:fir` or weight vector for `:firls`
- `window::Symbol=:hamming`: window for `:fir` (`:rect`, `:hann`, `:hamming`, `:blackman`)
- `dir::Symbol=:twopass`: filtering direction (`:twopass`, `:onepass`, `:reverse`)
- `flim::Tuple{Real, Real}=(0, fs / 2)`: frequency limits of the plot
- `n::Int64=4096`: number of frequency points within `flim`
- `mono::Bool=false`: if `true`, use a monochrome palette
- `gui::Bool=true`: if `true`, show an interactive window and wait until it is closed

# Returns

- `Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: the last valid filter, if `gui=true`
- `GLMakie.Figure`: the figure, if `gui=false`
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
    window::Symbol = :hamming,
    dir::Symbol = :twopass,
    flim::Tuple{Real, Real} = (0, fs / 2),
    n::Int64 = 4096,
    mono::Bool = false,
    gui::Bool = true,
)::Union{
    GLMakie.Figure,
    Vector{Float64},
    ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64},
    Biquad{:z, Float64},
}
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    _check_tuple(flim, (0, fs / 2), "flim")
    _check_var(dir, [:twopass, :onepass, :reverse], "dir")
    n >= 16 || throw(ArgumentError("n must be ≥ 16."))
    nqf = fs / 2    # Nyquist (not limited by flim)

    is_fir = fprototype in (:fir, :firls, :remez)
    is_iir = fprototype in (:butterworth, :chebyshev1, :chebyshev2, :elliptic)
    fprototype === :iirnotch && (ftype = nothing)
    is_interval = ftype in (:bp, :bs)
    custom_w = fprototype === :fir && !isnothing(w)

    # FIR without order/window vector: order follows bw
    order_auto = is_fir && isnothing(order) && isnothing(w)
    order_auto && isnothing(bw) &&
        throw(ArgumentError("bw or order must be specified for $fprototype."))

    # IIR ripple defaults (as in filter_create), so that sliders can start from them
    if fprototype in (:chebyshev1, :chebyshev2, :elliptic)
        isnothing(rp) && (rp = 0.5)
        isnothing(rs) && (rs = 20)
    end

    cutoff isa Tuple && cutoff[1] > cutoff[2] && (cutoff = (cutoff[2], cutoff[1]))

    design(c, o, b, p, s) = filter_create(;
        fprototype = fprototype,
        ftype = ftype,
        cutoff = c,
        fs = fs,
        order = (order_auto || fprototype === :iirnotch || custom_w) ? nothing :
                (isnothing(o) ? nothing : round(Int64, o)),
        rp = p,
        rs = s,
        bw = b,
        w = w,
        window = window,
    )

    # validate the initial design with the full filter_create checks (warnings shown)
    flt0 = design(cutoff, order, bw, rp, rs)

    v = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false

    try
        # reactive parameters
        cutoff_obs = Observable{Any}(cutoff isa Tuple ? Float64.(cutoff) : Float64(cutoff))
        order0 = flt0 isa Vector{Float64} ? length(flt0) : order
        order_obs = Observable{Union{Nothing, Float64}}(
            isnothing(order0) ? nothing : Float64(order0),
        )
        bw_obs = Observable{Union{Nothing, Float64}}(isnothing(bw) ? nothing : Float64(bw))
        rp_obs = Observable{Union{Nothing, Float64}}(isnothing(rp) ? nothing : Float64(rp))
        rs_obs = Observable{Union{Nothing, Float64}}(isnothing(rs) ? nothing : Float64(rs))
        status = Observable("")
        flt = Observable{Any}(flt0)

        # redesign on any change; invalid combinations keep the last valid filter
        onany(cutoff_obs, order_obs, bw_obs, rp_obs, rs_obs) do c, o, b, p, s
            try
                flt[] = design(c, o, b, p, s)
                status[] = ""
            catch e
                status[] = e isa ArgumentError ? e.msg : sprint(showerror, e)
            end
            return nothing
        end

        # frequency grid within flim
        f = collect(range(flim[1], flim[2]; length = n))
        resp = lift(x -> _filter_response(x, f, fs, dir), flt)
        H = lift(r -> r[1], resp)
        phi = lift(r -> r[2], resp)
        tau = lift(r -> r[3], resp)
        bw_eff = lift(
            (x, b) -> _bw_effective(x, fprototype, fs, window, custom_w, b),
            flt,
            bw_obs,
        )

        title1 = lift(
            flt,
            cutoff_obs,
            order_obs,
            bw_obs,
            bw_eff,
            rp_obs,
            rs_obs,
            status,
        ) do x, c, o, b, be, p, s, st
            return _filter_title(
                x;
                fprototype = fprototype,
                ftype = ftype,
                cutoff = c,
                order = isnothing(o) ? "" : round(Int64, o),
                bw = (fprototype === :fir && !order_auto) ? nothing : b,
                bw_eff = be,
                rp = p,
                rs = s,
                fs = fs,
                dir = dir,
                window = window,
                custom_w = custom_w,
                status = st,
            )
        end

        # figure
        GLMakie.activate!(; title = "plot_filter()")
        fig = GLMakie.Figure(; size = gui ? (1200, 950) : (1200, 850))

        ax1 = _filter_axis(fig, (1, 1), title1, "Magnitude [dB]", flim)
        GLMakie.ylims!(ax1, (-120, 10))
        GLMakie.hlines!(ax1, [-3, -6]; linestyle = :dot, linewidth = 0.75, color = :gray)
        GLMakie.lines!(ax1, f, H; color = mono ? :black : :blue)

        ax2 = _filter_axis(fig, (2, 1), "Phase response", "Phase [deg]", flim)
        GLMakie.lines!(ax2, f, phi; color = mono ? :black : :blue)

        ax3 = _filter_axis(fig, (3, 1), "Group delay", "Group delay [samples]", flim)
        GLMakie.lines!(ax3, f, tau; color = mono ? :black : :blue)
        on(tau) do t
            ft = Base.filter(isfinite, t)
            return isempty(ft) ||
                   GLMakie.ylims!(ax3, (min(0, minimum(ft)) - 1, max(0, maximum(ft)) + 1))
        end
        notify(tau)

        _draw_cutoff_vlines!((ax1, ax2, ax3), cutoff_obs, bw_eff, mono)

        # GUI sliders
        if gui
            grid = fig[4, 1] = GridLayout()
            row = 1
            _add_cutoff_slider!(grid, row, cutoff_obs, nqf, is_interval)
            row += 1

            # order: IIR always; FIR only when set manually (no window vector)
            if is_iir
                _add_slider!(grid, row, "Order", 1:1:max(20, order), order, order_obs)
                row += 1
            elseif is_fir && !order_auto && !custom_w
                odd = ftype in (:hp, :bp, :bs)
                omax = max(2001, 2 * order + 1)
                _add_slider!(
                    grid,
                    row,
                    "Order [taps]",
                    odd ? (1:2:omax) : (1:1:omax),
                    order,
                    order_obs,
                )
                row += 1
            end

            # bw: FIR with automatic order, :firls, :remez, :iirnotch
            if !isnothing(bw) && (fprototype !== :fir || order_auto)
                bmax = min(max(2 * bw, 10.0), nqf / 2)
                _add_slider!(grid, row, "Band width [Hz]", 0.05:0.05:bmax, bw, bw_obs)
                row += 1
            end

            if fprototype in (:chebyshev1, :elliptic)
                _add_slider!(grid, row, "RP [dB]", 0.1:0.1:10.0, rp, rp_obs)
                row += 1
            end
            if fprototype in (:chebyshev2, :elliptic)
                _add_slider!(grid, row, "RS [dB]", 1:1:100, rs, rs_obs)
                row += 1
            end

            wait(display(fig))
            return flt[]
        else
            return fig
        end
    finally
        NeuroAnalyzer.verbose = v
    end
end

"""
    plot_filter(obj; <keyword arguments>)

Plot the frequency response of a filter for the sampling rate of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (used only for the sampling rate)
- other arguments as in [`plot_filter`](@ref); `flim` defaults to `(0, sr(obj) / 2)`

# Returns

- `Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: the last valid filter, if `gui=true`
- `GLMakie.Figure`: the figure, if `gui=false`
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
    window::Symbol = :hamming,
    dir::Symbol = :twopass,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    n::Int64 = 4096,
    mono::Bool = false,
    gui::Bool = true,
)::Union{
    GLMakie.Figure,
    Vector{Float64},
    ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64},
    Biquad{:z, Float64},
}
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
        window = window,
        dir = dir,
        flim = flim,
        n = n,
        mono = mono,
        gui = gui,
    )
end
