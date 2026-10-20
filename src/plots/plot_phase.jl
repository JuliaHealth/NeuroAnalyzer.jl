export plot_polar
export plot_phase
export plot_ispc
export plot_pli
export plot_itpc
export plot_itpc_spec
export plot_itpc_frq

"""
    plot_polar(s; <keyword arguments>)

Create a polar plot from signal data with customizable styling.

# Arguments

- `s::Union{AbstractVector, AbstractMatrix}`: input data to plot; vector must contain 2 values: phases and lengths; matrix must contain 2 columns: phases and lengths
- `m::Tuple{Real, Real}=(0, 0)`: major value to plot
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `ticks::Bool=false`: draw x- and y-axis ticks

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_polar(
    s::Union{AbstractVector, AbstractMatrix};
    m::Tuple{Real, Real} = (0, 0),
    title::String = "",
    mono::Bool = false,
    ticks::Bool = true,
)::GLMakie.Figure
    # validate
    size(s, 1) == 2 && (s = s')
    length(m) == 2 ||
        throw(ArgumentError("m must contain 2 values: phases and lengths."))
    ndims(s) > 1 && size(s, 2) != 2 &&
        throw(ArgumentError("signal must contain 2 columns: phases and lengths."))

    # prepare plot
    GLMakie.activate!(; title = "plot_polar()")
    fig = GLMakie.Figure(; size = (800, 800))

    # create axis with customizable properties
    ax = GLMakie.PolarAxis(
        fig[1, 1];
        title         = title,
        thetazoomlock = true,
        rzoomlock     = true,
    )
    !ticks && hidespines!(ax)

    if ndims(s) == 1
        for idx in eachindex(s)
            GLMakie.lines!([0, s[idx]], [0, 1]; linewidth = 2, color = :black)
        end
    else
        for idx in axes(s, 1)
            GLMakie.lines!([0, s[idx, 1]], [0, s[idx, 2]]; linewidth = 2, color = :black)
        end
    end

    if m != (0, 0)
        GLMakie.lines!(
            [0, m[1]], [0, m[2]];
            linewidth = 2,
            color     = mono ? :darkgray : :red,
        )
    end

    return fig
end

"""
    plot_phase(ph, sf; <keyword arguments>)

Plot phase values from time-frequency analysis with customizable styling and units.

# Arguments

- `ph::Vector{Float64}`: vector of phase values (radians or degrees)
- `sf::Vector{Float64}`: vector of corresponding frequencies or time points
- `unit::Symbol=:rad`: phase unit specification (`:rad` radians or `:deg` degrees)
- `type::Symbol=:line`: plot type (`:line`: line plot, `:stem`: stem plot with markers)
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_phase(
    ph::Vector{Float64},
    sf::Vector{Float64};
    unit::Symbol = :rad,
    type::Symbol = :line,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
)::GLMakie.Figure
    # validate
    _check_var(unit, [:rad, :deg], "unit")
    _check_var(type, [:line, :stem], "type")
    length(ph) == length(sf) || throw(
        ArgumentError(
            "Length of phases ($(length(ph))) and frequencies ($(length(sf))) must be equal.",
        ),
    )

    # set default values
    xl, yl, tt = _set_defaults(
        xlabel, ylabel, title,
        "Frequency [Hz]",
        unit === :rad ? "Phase [rad]" : "Phase [°]",
        "",
    )

    # prepare plot
    GLMakie.activate!(; title = "plot_phase()")
    fig = GLMakie.Figure(; size = (900, 450))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xl,
        ylabel             = yl,
        title              = tt,
        _frq_axis_kwargs(:lin)...,
        xautolimitmargin   = (0, 0),
        yautolimitmargin   = (0.05, 0.05),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, _xlims(sf))
    _style_axis!(ax)

    if type === :line
        GLMakie.lines!(sf, ph; linewidth = 1, color = :black)
    else
        GLMakie.stem!(sf, ph; markersize = 10, color = :black)
    end

    return fig
end

# signals, signal difference, phases, phase difference and their polar views of a phase coupling
# measure (ISPC, PLI); `value` / `angle`: the length / angle of the resultant vector
function _plot_phase_coupling(
    t::AbstractVector,
    s1::AbstractVector,
    s2::AbstractVector,
    sd::AbstractVector,
    phd::AbstractVector,
    s1ph::AbstractVector,
    s2ph::AbstractVector,
    value::Real,
    angle::Real,
    name::String,
    unit::String,
)::GLMakie.Figure
    # validate
    n = length(t)
    all(==(n), length.((s1, s2, sd, phd, s1ph, s2ph))) ||
        throw(ArgumentError("t, s1, s2, sd, phd, s1ph and s2ph must have the same length."))

    GLMakie.activate!(; title = "plot_$(lowercase(name))()")
    fig = GLMakie.Figure(; size = (1200, 900))
    yl = isempty(unit) ? "Amplitude" : "Amplitude [$unit]"
    function line_axis(pos, ylabel, title)
        ax = GLMakie.Axis(
            fig[pos...];
            xlabel = "Time [s]",
            ylabel = ylabel,
            title = title,
            _TIME_AXIS_KWARGS...,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0.05, 0.05),
            _AXIS_LOCK_KWARGS...,
        )
        _style_axis!(ax)
        return ax
    end
    ax1 = line_axis((1, 1), yl, "Signals")
    GLMakie.lines!(ax1, t, s1; color = :black)
    GLMakie.lines!(ax1, t, s2; color = :gray)
    ax2 = line_axis((1, 2), yl, "Signal difference")
    GLMakie.linkyaxes!(ax1, ax2)
    GLMakie.lines!(ax2, t, sd; color = :black, linewidth = 1)
    ax3 = line_axis((2, 1), "Phase [rad]", "Phases")
    GLMakie.lines!(ax3, t, s1ph; color = :black)
    GLMakie.lines!(ax3, t, s2ph; color = :gray)
    ax4 = line_axis((2, 2), "Phase [rad]", "Phase difference")
    GLMakie.linkyaxes!(ax3, ax4)
    GLMakie.lines!(ax4, t, phd; color = :black, linewidth = 1)

    # polar views: unit vectors at every phase (one line segment each, drawn as a single plot)
    rays(ph) = (vec(vcat(zeros(length(ph))', ph', fill(NaN, length(ph))')),
                vec(vcat(zeros(length(ph))', ones(length(ph))', fill(NaN, length(ph))')))
    ax5 = GLMakie.PolarAxis(fig[3, 1]; title = "Phases", thetazoomlock = true, rzoomlock = true)
    GLMakie.hiderdecorations!(ax5)
    GLMakie.lines!(ax5, rays(s1ph)...; color = :black, linewidth = 0.5)
    GLMakie.lines!(ax5, rays(s2ph)...; color = :gray, linewidth = 0.5)
    ax6 = GLMakie.PolarAxis(
        fig[3, 2];
        title = "Phase difference and $name = $(round(value; digits = 2))",
        thetazoomlock = true,
        rzoomlock = true,
    )
    GLMakie.hiderdecorations!(ax6)
    GLMakie.lines!(ax6, rays(phd)...; color = :black, linewidth = 0.5)
    GLMakie.lines!(ax6, [0, angle], [0, value]; color = :red, linewidth = 2)

    return fig
end

"""
    plot_ispc(t, s1, s2; <keyword arguments>)

Plot ISPC (Inter-Subject-Phase Clustering) of two signals, as calculated by `ispc()`.

# Arguments

- `t::AbstractVector`: time points
- `s1::AbstractVector`: signal 1
- `s2::AbstractVector`: signal 2
- `ispcv::Real`: ISPC value
- `ispca::Real`: ISPC angle
- `sd::AbstractVector`: signal difference (s2 - s1)
- `phd::AbstractVector`: phase difference (s2 - s1)
- `s1ph::AbstractVector`: signal 1 phases
- `s2ph::AbstractVector`: signal 2 phases
- `unit::String=""`: signal unit

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_ispc(
    t::AbstractVector,
    s1::AbstractVector,
    s2::AbstractVector;
    ispcv::Real,
    ispca::Real,
    sd::AbstractVector,
    phd::AbstractVector,
    s1ph::AbstractVector,
    s2ph::AbstractVector,
    unit::String = "",
)::GLMakie.Figure
    return _plot_phase_coupling(t, s1, s2, sd, phd, s1ph, s2ph, ispcv, ispca, "ISPC", unit)
end

"""
    plot_pli(t, s1, s2; <keyword arguments>)

Plot PLI (Phase Lag Index) of two signals, as calculated by `pli()`.

# Arguments

- `t::AbstractVector`: time points
- `s1::AbstractVector`: signal 1
- `s2::AbstractVector`: signal 2
- `pv::Real`: PLI value
- `sd::AbstractVector`: signal difference (s2 - s1)
- `phd::AbstractVector`: phase difference (s2 - s1)
- `s1ph::AbstractVector`: signal 1 phases
- `s2ph::AbstractVector`: signal 2 phases
- `unit::String=""`: signal unit

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_pli(
    t::AbstractVector,
    s1::AbstractVector,
    s2::AbstractVector;
    pv::Real,
    sd::AbstractVector,
    phd::AbstractVector,
    s1ph::AbstractVector,
    s2ph::AbstractVector,
    unit::String = "",
)::GLMakie.Figure
    return _plot_phase_coupling(t, s1, s2, sd, phd, s1ph, s2ph, pv, 0.0, "PLI", unit)
end

"""
    plot_itpc(itpcph; <keyword arguments>)

Plot ITPC (Inter-Trial-Phase Clustering) at a time point, as calculated by `itpc()`: the distribution of the phases across trials and their resultant vector.

# Arguments

- `itpcph::AbstractVector`: phases of the trials at the time point
- `itpcv::Real`: ITPC (or ITPCz / weighted ITPC) value
- `itpca::Real`: ITPC angle
- `label::String="ITPC"`: name of the value (e.g. `"ITPCz"`, `"wITPC"`)
- `title::String=""`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_itpc(
    itpcph::AbstractVector;
    itpcv::Real,
    itpca::Real,
    label::String = "ITPC",
    title::String = "",
)::GLMakie.Figure
    isempty(itpcph) && throw(ArgumentError("itpcph must not be empty."))

    GLMakie.activate!(; title = "plot_itpc()")
    fig = GLMakie.Figure(; size = (1200, 600))
    !isempty(title) && GLMakie.Label(fig[0, 1:2], title; fontsize = 20)
    ax1 = GLMakie.Axis(
        fig[1, 1];
        xlabel = "Phase angle [rad]",
        ylabel = "Count/bin",
        title = "Phase angles across trials",
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax1)
    GLMakie.hist!(ax1, itpcph; bins = 15, color = :grey)
    ax2 = GLMakie.PolarAxis(
        fig[1, 2];
        title = "Phases, $label = $(round(itpcv; digits = 2))",
        thetazoomlock = true,
        rzoomlock = true,
    )
    GLMakie.hiderdecorations!(ax2)
    n = length(itpcph)
    GLMakie.lines!(ax2, vec(vcat(zeros(n)', itpcph', fill(NaN, n)')),
                   vec(vcat(zeros(n)', ones(n)', fill(NaN, n)')); color = :black, linewidth = 0.5)
    GLMakie.lines!(ax2, [0, itpca], [0, itpcv]; color = :red, linewidth = 2)

    return fig
end

"""
    plot_itpc_spec(t, f, itpcs; <keyword arguments>)

Plot the time-frequency map of ITPC (Inter-Trial-Phase Clustering), as calculated by `itpc_spec()`.

# Arguments

- `t::AbstractVector`: time points
- `f::AbstractVector`: frequencies
- `itpcs::AbstractMatrix`: ITPC (or ITPCz), frequencies × time points
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits
- `frq::Symbol=:lin`: frequency axis scaling (`:lin` or `:log`)
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Frequency [Hz]"`: y-axis label
- `title::String="ITPC spectrogram"`: plot title
- `cb_title::String="ITPC"`: color bar title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_itpc_spec(
    t::AbstractVector,
    f::AbstractVector,
    itpcs::AbstractMatrix;
    flim::Tuple{Real, Real} = (f[1], f[end]),
    frq::Symbol = :lin,
    xlabel::String = "Time [s]",
    ylabel::String = "Frequency [Hz]",
    title::String = "ITPC spectrogram",
    cb_title::String = "ITPC",
    mono::Bool = false,
)::GLMakie.Figure
    return plot_spectrogram(
        Vector{Float64}(t),
        Vector{Float64}(f),
        Matrix{Float64}(itpcs);
        db = false,
        frq = frq,
        flim = flim,
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        cb_title = cb_title,
        mono = mono,
    )
end

"""
    plot_itpc_frq(t, s; <keyword arguments>)

Plot ITPC (Inter-Trial-Phase Clustering) over time at one frequency (a row of the `itpc_spec()` result).

# Arguments

- `t::AbstractVector`: time points
- `s::AbstractVector`: ITPC (or ITPCz) at the frequency
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="ITPC"`: y-axis label
- `title::String=""`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_itpc_frq(
    t::AbstractVector,
    s::AbstractVector;
    xlabel::String = "Time [s]",
    ylabel::String = "ITPC",
    title::String = "",
)::GLMakie.Figure
    length(t) == length(s) ||
        throw(ArgumentError("Lengths of t ($(length(t))) and s ($(length(s))) must be equal."))

    GLMakie.activate!(; title = "plot_itpc_frq()")
    fig = GLMakie.Figure(; size = (900, 450))
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        _TIME_AXIS_KWARGS...,
        xautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)
    # ITPC is in [0, 1]; ITPCz is not bounded
    maximum(s) <= 1 && GLMakie.ylims!(ax, 0, 1)
    GLMakie.lines!(ax, t, s; color = :black)

    return fig
end
