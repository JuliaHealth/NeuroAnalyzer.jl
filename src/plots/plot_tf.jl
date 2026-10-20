export plot_eros
export plot_erop
export plot_imf
export plot_fi
export plot_dwc
export plot_env

"""
    plot_eros(s, f, t; <keyword arguments>)

Plot an Event-Related Oscillations (ERO) spectrogram with customizable visualization options.

# Arguments

- `sp::AbstractArray`: ERO spectrogram power values, shape (frequency, time, powers)
- `sf::AbstractVector`: vector of frequency values in Hz
- `st::AbstractVector`: vector of time values in seconds
- `db::Bool=true`: whether to display power values in decibels
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plot
- `tm::Union{Nothing, Int64, Vector{Int64}} = nothing`: time markers (in milliseconds) to be plot as vertical lines, useful for adding topoplots at these time points
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `cb::Bool=true`: if `true`, show colorbar
- `mono::Bool=false`: if `true`, use a monochrome palette
- `units::String="μV"`: power units
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_eros(
    sp::AbstractArray,
    sf::AbstractVector,
    st::AbstractVector;
    db::Bool = true,
    frq::Symbol = :lin,
    flim::Tuple{Real, Real} = (sf[1], sf[end]),
    tm::Union{Nothing, Int64, Vector{Int64}} = nothing,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    cb::Bool = true,
    mono::Bool = false,
    units::String = "μV",
    smooth::Bool = false,
    ks::Int64 = 3,
)::GLMakie.Figure
    # validate
    size(sp, 1) == length(sf) || throw(
        ArgumentError(
            "Length of sf ($(length(sf))) and number of spectrogram rows ($(size(sp, 1))) must be equal.",
        ),
    )
    size(sp, 2) == length(st) || throw(
        ArgumentError(
            "Length of st ($(length(st))) and number of spectrogram columns ($(size(sp, 2))) must be equal.",
        ),
    )
    ndims(sp) == 3 || throw(ArgumentError("sp must have 3 dimensions."))
    size(sp, 3) <= 2 || throw(ArgumentError("sp must contain ≤ 2 epochs."))
    ks > 0 || throw(ArgumentError("ks must be ≥ 1."))
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    # set color palette
    pal = mono ? :grays : :darktest

    # colorbar title
    cb_title = db ? "[dB $units^2/Hz]" : "[$units^2/Hz]"

    # set frequency limits
    if frq === :lin
        yt = flim[2] > 100 ? (flim[1]:10:flim[2]) : (flim[1]:5:flim[2])
    else
        flim = _log_flim(sf, flim, frq)
        # log axis: 10 ticks between the frequency limits
        yt = round.(logspace(flim[1], flim[2], 10); digits = 1)
    end

    # apply Gaussian filter if requested
    if smooth
        for idx in axes(sp, 3)
            sp[:, :, idx] = imfilter(@view(sp[:, :, idx]), Kernel.gaussian(ks))
        end
    end

    # resolve time markers to indices into a new local variable — do not mutate the input
    tm_indices = if !isnothing(tm)
        markers = tm isa Int64 ? [tm] : tm
        for val in markers
            val / 1000 >= st[1] || throw(
                ArgumentError(
                    "tm value ($val) is out of epoch time segment ($(st[1]):$(st[end])).",
                ),
            )
            val / 1000 <= st[end] || throw(
                ArgumentError(
                    "tm value ($val) is out of epoch time segment ($(st[1]):$(st[end])).",
                ),
            )
        end
        [vsearch(val / 1000, st) for val in markers]
    else
        Int64[]
    end

    # log frequency axis: only positive frequencies (heatmap cell edges must be > 0)
    fr = frq === :log ? findall(>(0), sf) : collect(eachindex(sf))

    function _draw_axis(fig, pos, xl, yl, tt, nticks)
        ax = GLMakie.Axis(
            fig[pos...];
            xlabel             = xl,
            ylabel             = yl,
            title              = tt,
            xticks             = LinearTicks(nticks),
            xminorticksvisible = true,
            xminorticks        = IntervalsBetween(10),
            yticks             = yt,
            yscale             = frq === :lin ? identity : log,
            _AXIS_LOCK_KWARGS...,
        )
        GLMakie.ylims!(ax, flim)
        _style_axis!(ax)
        return ax
    end

    function _draw_markers!(fig, pos, indices)
        for i in indices
            GLMakie.vlines!(fig[pos...], [st[i]]; color = :black, linewidth = 1)
        end
    end

    if size(sp, 3) == 1

        # set default values
        xl, yl, tt = _set_defaults(
            xlabel, ylabel, title,
            "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of epochs",
        )

        # prepare plot
        GLMakie.activate!(; title = "plot_eros()")
        fig = GLMakie.Figure(; size = (900, 450))

        # create axis with customizable properties
        ax = _draw_axis(fig, (1, 1), xl, yl, tt, 15)
        hm = GLMakie.heatmap!(ax, st, sf[fr], sp[fr, :, 1]'; colormap = pal)
        cb && GLMakie.Colorbar(fig[1, 2], hm; label = cb_title, labelsize = 16)
        _draw_markers!(fig, (1, 1), tm_indices)

    else

        # set default values
        xl1, yl1, tt1 = _set_defaults(
            xlabel, ylabel, title,
            "Time [ms]", "Frequency [Hz]", "ERP spectrogram",
        )
        xl2, yl2, tt2 = _set_defaults(
            xlabel, ylabel, title,
            "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of ERP epochs",
        )

        # prepare plot
        GLMakie.activate!(; title = "plot_eros()")
        fig = GLMakie.Figure(; size = (1200, 800))

        # create axis with customizable properties
        ax1 = _draw_axis(fig, (1, 1), xl1, yl1, tt1, 10)
        hm1 = GLMakie.heatmap!(ax1, st, sf[fr], sp[fr, :, 1]'; colormap = pal)
        cb && GLMakie.Colorbar(fig[1, 2], hm1; label = cb_title, labelsize = 16)
        GLMakie.vlines!(ax1, [0]; linestyle = :dash, linewidth = 0.5, color = :black)
        _draw_markers!(fig, (1, 1), tm_indices)

        # create axis with customizable properties
        ax2 = _draw_axis(fig, (2, 1), xl2, yl2, tt2, 10)
        hm2 = GLMakie.heatmap!(ax2, st, sf[fr], sp[fr, :, 2]'; colormap = pal)
        cb && GLMakie.Colorbar(fig[2, 2], hm2; label = cb_title, labelsize = 16)
        GLMakie.vlines!(ax2, [0]; linestyle = :dash, linewidth = 0.5, color = :black)
        _draw_markers!(fig, (2, 1), tm_indices)
    end

    return fig
end

"""
    plot_erop(fig, f; <keyword arguments>)

Plot the power spectrum of Event-Related Oscillations (ERO) with customizable visualization.

# Arguments

- `sp::AbstractArray`: ERO power values, shape (frequency, powers)
- `sf::AbstractVector`: vector of frequency values in Hz
- `db::Bool=true`: whether to display power values in decibels
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plot
- `cb::Bool=true`: if `true`, show colorbar
- `units::String="μV"`: power units
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_erop(
    sp::AbstractArray,
    sf::AbstractVector;
    db::Bool = true,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    flim::Tuple{Real, Real} = (sf[1], sf[end]),
    frq::Symbol = :lin,
    units::String = "μV",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    _in(flim[2], (sf[1], sf[end]), "flim")
    size(sp, 1) == length(sf) || throw(
        ArgumentError(
            "Length of sf ($(length(sf))) and number of power rows ($(size(sp, 1))) must be equal.",
        ),
    )
    ndims(sp) == 2 || throw(ArgumentError("sp must have 2 dimensions."))
    size(sp, 2) <= 2 || throw(ArgumentError("sp must contain ≤ 2 epochs."))
    _check_var(frq, [:lin, :log], "frq")

    # set frequency limits
    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        flim = (sf[2], flim[2])
    end

    # set y-axis label
    power_ylabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]"

    function _make_power_axis(fig, pos, xl, yl, tt)
        ax = GLMakie.Axis(
            fig[pos...];
            xlabel             = xl,
            ylabel             = yl,
            title              = tt,
            xminorticksvisible = true,
            xminorticks        = IntervalsBetween(10),
            xscale             = frq === :lin ? identity : log,
            _AXIS_LOCK_KWARGS...,
        )
        GLMakie.xlims!(ax, flim)
        _style_axis!(ax)
        return ax
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_erop()")

    if size(sp, 2) == 1

        # set default values
        xl, _, tt = _set_defaults(
            xlabel, ylabel, title,
            "Frequency [Hz]", power_ylabel, "Averaged power-spectra of epochs",
        )

        # prepare plot
        fig = GLMakie.Figure(; size = (900, 450))

        # create axis with customizable properties
        ax = _make_power_axis(fig, (1, 1), xl, power_ylabel, tt)
        GLMakie.lines!(ax, sf, sp[:, 1]; color = :black)

    else

        # set default values
        xl, _, tt1 = _set_defaults(
            xlabel, ylabel, title,
            "Frequency [Hz]", power_ylabel, "ERP power-spectrum",
        )
        _, _, tt2 = _set_defaults(
            xlabel, ylabel, title,
            "Frequency [Hz]", power_ylabel, "Averaged power-spectra of ERP epochs",
        )

        # prepare plot
        fig = GLMakie.Figure(; size = (1200, 800))

        # create axis with customizable properties
        ax1 = _make_power_axis(fig, (1, 1), xl, power_ylabel, tt1)
        GLMakie.lines!(ax1, sf, sp[:, 1]; color = :black)

        ax2 = _make_power_axis(fig, (2, 1), xl, power_ylabel, tt2)
        GLMakie.lines!(ax2, sf, sp[:, 2]; color = :black)
    end

    return fig
end

"""
    plot_imf(imf; <keyword arguments>)

Plot intrinsic mode functions (IMF), the residual, and the reconstructed signal from Empirical Mode Decomposition (EMD).

# Arguments

- `imf::Matrix{Float64}`: matrix where each row represents an IMF component from EMD
- `n::Int64=size(imf, 1) - 1`: number of IMF components to plot
- `t::AbstractVector`: time points corresponding to the signal values

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_imf(
    imf::Matrix{Float64};
    n::Int64 = size(imf, 1) - 1,
    t::AbstractVector,
)::GLMakie.Figure
    # validate
    n > 0 || throw(ArgumentError("n must be ≥ 1."))
    n + 1 <= size(imf, 1) || throw(ArgumentError("n must be ≤ $(size(imf, 1) - 1)."))
    size(imf, 2) == length(t) || throw(
        ArgumentError(
            "Length of t ($(length(t))) and number of imf columns ($(size(imf, 2))) must be equal.",
        ),
    )

    # the last row of imf is the residual; the reconstruction is the sum of all rows
    s_restored = sum(imf; dims = 1)[:]
    imf_plot   = vcat(imf, s_restored')

    # set y-axis limits
    ylim   = (floor(minimum(imf_plot); digits = 0), ceil(maximum(imf_plot); digits = 0))
    ylim   = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # prepare plot
    GLMakie.activate!(; title = "plot_imf()")
    fig = GLMakie.Figure(; size = (1200, 800))
    nr  = ceil(Int64, (n + 1) / 2)

    idx  = 1
    cidx = 1
    for idx1 = 1:nr
        cidx = 1
        for idx2 = 1:2
            if idx <= n + 1
                label = idx == n + 1 ? "Residual" : "IMF: $idx"
                # create axis with customizable properties
                ax = GLMakie.Axis(
                    fig[idx1, idx2];
                    xlabel             = "Time [s]",
                    ylabel             = "",
                    title              = label,
                    xticks             = LinearTicks(10),
                    xminorticksvisible = true,
                    xminorticks        = IntervalsBetween(10),
                    yticks             = yticks,
                    _AXIS_LOCK_KWARGS...,
                )
                GLMakie.ylims!(ax, ylim)
                _style_axis!(ax)
                GLMakie.lines!(ax, t, imf[idx, :]; color = :black)
                idx  += 1
                cidx += 1
            end
        end
    end

    # place the reconstructed signal spanning both columns
    row = cidx == 1 ? nr : nr + 1
    ax  = GLMakie.Axis(
    fig[row, 1:2];
    xlabel             = "Time [s]",
    ylabel             = "",
    title              = "Reconstructed signal",
    xticks             = LinearTicks(10),
    xminorticksvisible = true,
    xminorticks        = IntervalsBetween(10),
    yticks             = yticks,
    _AXIS_LOCK_KWARGS...
)
    GLMakie.ylims!(ax, ylim)
    _style_axis!(ax)
    GLMakie.lines!(ax, t, s_restored; color = :black)

    return fig
end

"""
    plot_fi(fi, st; <keyword arguments>)

Plot instantaneous frequency over time from time-frequency analysis (e.g., Hilbert-Huang Transform).

# Arguments

- `fi::Vector{Float64}`: vector of instantaneous frequency values in Hz
- `st::Vector{Float64}`: vector of time points in seconds corresponding to frequency values
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Power [μV^2/Hz]
- `title::String="default"`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_fi(
    fi::Vector{Float64},
    st::Vector{Float64};
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
)::GLMakie.Figure
    # validate
    length(fi) == length(st) || throw(
        ArgumentError(
            "Length of frequencies ($(length(fi))) and time points ($(length(st))) must be equal.",
        ),
    )

    # set default values
    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Frequency [Hz]", "")

    # prepare plot
    GLMakie.activate!(; title = "plot_fi()")
    fig = GLMakie.Figure(; size = (900, 450))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xl,
        ylabel             = yl,
        title              = tt,
        xticks             = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, _xlims(st))
    _style_axis!(ax)

    GLMakie.lines!(st, fi; linewidth = 1, color = :black)

    return fig
end

"""
    plot_dwc(dc; <keyword arguments>)

Plot discrete wavelet decomposition (DWC) coefficients showing signal decomposition across different scales.

# Arguments

- `dc::Matrix{Float64}`: matrix where each row represents wavelet coefficients at a specific decomposition level
- `n::Int64=size(dc, 1) - 1`: number of decomposition levels to plot
- `t::AbstractVector`: time points corresponding to the signal values

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_dwc(
    dc::Matrix{Float64};
    n::Int64 = size(dc, 1) - 1,
    t::AbstractVector,
)::GLMakie.Figure
    # validate
    n > 1 || throw(ArgumentError("n must be > 1."))
    n <= size(dc, 1) - 1 || throw(ArgumentError("n must be ≤ $(size(dc, 1) - 1)."))
    size(dc, 2) == length(t) || throw(
        ArgumentError(
            "Length of t ($(length(t))) and number of dc columns ($(size(dc, 2))) must be equal.",
        ),
    )

    # set y-axis limits
    ylim   = (floor(minimum(dc); digits = 0), ceil(maximum(dc); digits = 1))
    ylim   = _tuple_max(ylim)
    yticks = unique([ylim[1], 0, ylim[2]])

    # prepare plot
    GLMakie.activate!(; title = "plot_dwc()")
    fig = GLMakie.Figure(; size = (1200, 800))
    nr  = ceil(Int64, (n + 1) / 2)

    # dc[1, :] is the original signal; coefficients start at dc[2, :]
    idx  = 2
    cidx = 1
    for idx1 = 1:nr
        cidx = 1
        for idx2 = 1:2
            if idx < n + 2
                ax = GLMakie.Axis(
                    fig[idx1, idx2];
                    xlabel             = "Time [s]",
                    ylabel             = "",
                    title              = "Coefficient #$(idx - 1)",
                    xticks             = LinearTicks(10),
                    xminorticksvisible = true,
                    xminorticks        = IntervalsBetween(10),
                    yticks             = yticks,
                    _AXIS_LOCK_KWARGS...,
                )
                GLMakie.ylims!(ax, ylim)
                _style_axis!(ax)
                GLMakie.lines!(ax, t, dc[idx, :]; color = :black)
                idx  += 1
                cidx += 1
            end
        end
    end

    # place the original signal spanning both columns
    row = cidx == 1 ? nr : nr + 1
    ax  = GLMakie.Axis(
    fig[row, 1:2];
    xlabel             = "Time [s]",
    ylabel             = "",
    title              = "Original signal",
    xticks             = LinearTicks(10),
    xminorticksvisible = true,
    xminorticks        = IntervalsBetween(10),
    yticks             = yticks,
    _AXIS_LOCK_KWARGS...
)
    GLMakie.xlims!(ax, _xlims(t))
    GLMakie.ylims!(ax, ylim)
    _style_axis!(ax)
    GLMakie.lines!(ax, t, dc[1, :]; color = :black)

    return fig
end

"""
    plot_env(x, e; <keyword arguments>)

Plot an envelope (e.g. returned by `tenv()`, `henv()`, `penv()`, `senv()` or their `_mean` / `_median` variants), optionally with its 95% confidence interval.

# Arguments

- `x::AbstractVector`: time points (amplitude, Hilbert and spectrogram envelopes) or frequencies (power envelope)
- `e::AbstractVector`: envelope
- `el::Union{Nothing, AbstractVector}=nothing`: lower bound of the 95% CI
- `eu::Union{Nothing, AbstractVector}=nothing`: upper bound of the 95% CI
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Amplitude"`: y-axis label
- `title::String=""`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_env(
    x::AbstractVector,
    e::AbstractVector;
    el::Union{Nothing, AbstractVector} = nothing,
    eu::Union{Nothing, AbstractVector} = nothing,
    xlabel::String = "Time [s]",
    ylabel::String = "Amplitude",
    title::String = "",
)::GLMakie.Figure
    # validate
    length(x) == length(e) ||
        throw(ArgumentError("Lengths of x ($(length(x))) and e ($(length(e))) must be equal."))
    isnothing(el) == isnothing(eu) || throw(ArgumentError("Both el and eu must be provided."))
    if !isnothing(el)
        length(el) == length(eu) == length(e) ||
            throw(ArgumentError("el and eu must have the same length as e."))
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_env()")
    fig = GLMakie.Figure(; size = (900, 450))
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0.05, 0.05),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)

    # 95% CI
    !isnothing(el) && GLMakie.band!(ax, x, el, eu; alpha = 0.25, color = :grey, strokewidth = 0.5)
    GLMakie.lines!(ax, x, e; linewidth = 2, color = :black)

    return fig
end
