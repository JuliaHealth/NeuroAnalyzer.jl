export plot_acf

"""
    plot_acf(obj; <keyword arguments>)

Plot the auto-correlation function (ACF, see `acor()`) of channels, averaged across epochs, with auto-correlation windows (ACW, see `acw()`) marked on each curve.

ACW markers are the ACW values averaged across epochs (not the ACW of the averaged ACF); dashed lines show their levels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `ep::Union{Nothing, Int64, AbstractRange{<:Integer}, Vector{<:Integer}}=nothing`: epoch(s) to average; default: all epochs (for ERP: the averaged epoch)
- `l::Real=1.0`: maximum lag in seconds (limited to the epoch length)
- `xlim::Union{Nothing, Real}=nothing`: upper limit of the lag axis in seconds; default: `l`
- `markers::Vector{Symbol}=[:acw50, :acw0]`: ACWs to mark: `:acw50`, `:acwe`, `:acw0`
- `method::Symbol=:sum`: method of calculating auto-correlation (see `acor()`)
- `title::String="Auto-correlation function"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `leg::Bool=true`: if `true`, show the channel legend

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_acf(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ep::Union{Nothing, Int64, AbstractRange{<:Integer}, Vector{<:Integer}} = nothing,
    l::Real = 1.0,
    xlim::Union{Nothing, Real} = nothing,
    markers::Vector{Symbol} = [:acw50, :acw0],
    method::Symbol = :sum,
    title::String = "Auto-correlation function",
    mono::Bool = false,
    leg::Bool = true,
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)::GLMakie.Figure
    # validate
    l > 0 || throw(ArgumentError("l must be > 0."))
    isnothing(xlim) || xlim > 0 || throw(ArgumentError("xlim must be > 0."))
    for m in markers
        _check_var(m, [:acw50, :acwe, :acw0], "markers")
    end
    ch_idx = _get_ch(obj, ch, exclude_bad_channels)
    isempty(ch_idx) && throw(ArgumentError("No channels selected."))

    # maximum lag in samples, limited to the epoch length
    l_smp = min(round(Int64, l * sr(obj)), epoch_len(obj) - 1)

    ac = acor(
        obj;
        ch = ch,
        l = l_smp,
        method = method,
        exclude_bad_channels = exclude_bad_channels,
    )

    # lags ≥ 0, in ms
    lags = ac.lags[(l_smp + 1):end]
    lags_ms = lags .* 1000
    r = ac.autocor[:, (l_smp + 1):end, :]

    # epochs to average: for ERP the averaged epoch (epoch 1 of the acor() output)
    ep_idx = if isnothing(ep)
        datatype(obj) in ("erp", "erf") ? [1] : collect(1:size(r, 3))
    else
        collect(ep)
    end
    all(in(1:size(r, 3)), ep_idx) ||
        throw(ArgumentError("ep must be in [1, $(size(r, 3))]."))

    ch_n = size(r, 1)
    clabels = labels(obj)[ch_idx]
    levels = Dict(:acw50 => 0.5, :acwe => exp(-1), :acw0 => 0.0)
    shapes = Dict(:acw50 => :circle, :acwe => :diamond, :acw0 => :rect)

    # prepare plot
    GLMakie.activate!(; title = "plot_acf()")
    fig = GLMakie.Figure(; size = (800, 400))
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = "Lag [ms]",
        ylabel             = "Auto-correlation",
        title              = title,
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)

    # reference lines at the marked levels
    isempty(markers) ||
        GLMakie.hlines!(ax, [levels[m] for m in markers]; color = :gray, linestyle = :dash)

    cmap = _trace_colors(ch_n; mono = mono)
    for idx in 1:ch_n
        m = vec(mean(@view(r[idx, :, ep_idx]); dims = 2))
        GLMakie.lines!(ax, lags_ms, m; color = cmap[idx], linewidth = 1.5, label = clabels[idx])
        for mk in markers
            # ACW per epoch, averaged across the epochs where the level is reached
            a = [_acw(@view(r[idx, :, e]), lags, levels[mk]) for e in ep_idx]
            a = filter(!isnan, a)
            isempty(a) && continue
            GLMakie.scatter!(
                ax,
                [mean(a) * 1000],
                [levels[mk]];
                color      = cmap[idx],
                marker     = shapes[mk],
                markersize = 12,
            )
        end
    end

    GLMakie.xlims!(ax, 0, 1000 * (isnothing(xlim) ? lags[end] : min(xlim, lags[end])))

    # legend outside the plot (not drawn if it does not fit)
    leg && _side_legend!(fig, ax)

    return fig
end
