export plot_phsd
export plot_phsd_3d
export plot_phsd_topo

"""
    plot_phsd(f, ph; <keyword arguments>)

Plot PHSD (phase spectral density).

# Arguments

- `f::Vector{Float64}`: vector of frequency values in Hz
- `ph::Vector{Float64}`: vector of phase values in radians
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plots
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_phsd(
    f::Vector{Float64},
    ph::Vector{Float64};
    flim::Tuple{Real, Real} = (f[1], f[end]),
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    frq::Symbol = :lin,
)::GLMakie.Figure
    # validate
    length(ph) == length(f) || throw(
        ArgumentError("Length of powers vector must equal length of frequencies vector."),
    )
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(f), "flim")

    # prepare log-scaled frequencies axis
    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(f[2]) Hz.")
        flim = (f[2], flim[2])
    end

    # frequency limits
    f1 = vsearch(flim[1], f)
    f2 = vsearch(flim[2], f)

    # prepare plot
    GLMakie.activate!(; title = "plot_phsd()")
    plot_size = (900, 450)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(5),
        xscale = frq === :lin ? identity : log,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0.1, 0.1),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.autolimits!(ax)
    _style_axis!(ax)

    # draw phases
    GLMakie.lines!(
        ax,
        f[f1:f2],
        ph[f1:f2];
        linewidth = 2,
        color = :black,
    )

    return fig
end

"""
    plot_phsd(f, ph; <keyword arguments>)

Plot multi-channel Phase Spectral Density (PHSD) with customizable visualization options.

# Arguments

- `f::Vector{Float64}`: vector of frequency values in Hz
- `ph::Matrix{Float64}`: matrix of phase values, shape (channels, frequencies)
- `clabels::Vector{String}=string.(1:size(sp, 1))`: channel labels
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plots
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `avg::Bool=false`: if `true`, plot averaged PHSD
- `ci95::Bool=false`: if `true`, plot mean and ±95% CI of averaged PHSDs
- `leg::Bool=true`: if `true`, add legend with channel labels

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_phsd(
    f::Vector{Float64},
    ph::Matrix{Float64};
    clabels::Vector{String} = string.(1:size(ph, 1)),
    flim::Tuple{Real, Real} = (f[1], f[end]),
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    mono::Bool = false,
    frq::Symbol = :lin,
    avg::Bool = false,
    ci95::Bool = false,
    leg::Bool = true,
)::GLMakie.Figure
    # validate
    size(ph, 2) == length(f) || throw(
        ArgumentError("Length of powers vector must equal length of frequencies vector."),
    )
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(f), "flim")
    avg && ci95 && throw(ArgumentError("avg and ci95 cannot both be true."))

    # set color palette
    pal = mono ? :grays : :darktest

    # number of channels
    ch_n = size(p, 1)

    # frequency limits
    f1 = vsearch(flim[1], f)
    f2 = vsearch(flim[2], f)

    # get mean and 95%CI
    if ci95
        msci95_data = NeuroAnalyzer.msci95(ph[:, f1:f2])
        s_m = msci95_data.sm
        s_l = msci95_data.ll
        s_u = msci95_data.ul
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_phsd()")
    plot_size = (900, 450)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax  = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xlabel,
        ylabel             = ylabel,
        title              = title,
        xticks             = LinearTicks(15),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        xscale             = frq === :lin ? identity : log,
        xautolimitmargin   = (0, 0),
        yautolimitmargin   = (0.1, 0.1),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.autolimits!(ax)
    _style_axis!(ax)

    if ci95
        # draw 95% CI
        GLMakie.band!(ax, f[f1:f2], s_u, s_l; alpha = 0.25, color = :grey, strokewidth = 0.5)

        # draw mean
        GLMakie.lines!(ax, f[f1:f2], s_m; color = :black, linewidth = 2)
    else
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            Makie.lines!(
                ax,
                f[f1:f2],
                p[idx, f1:f2];
                color = cmap[idx],
                colormap = pal,
                colorrange = 1:ch_n,
                linewidth = 2,
                label = clabels[idx],
            )
        end

        # draw averaged channels
        if avg
            s = mean(p[f1:f2]; dims = 1)[:]
            GLMakie.lines!(ax, f[f1:f2], s; linewidth = 4, color = :black)

        end

        # add legend if requested
        (leg && ch_n < 30) && axislegend(; position = :rt, colormap = pal)
    end

    return fig
end

"""
    plot_phsd_3d(f, ph; <keyword arguments>)

Plot a 3D representation of multi-channel Phase Spectral Density (PHSD).

# Arguments

- `f::Vector{Float64}`: vector of frequency values in Hz
- `ph::Matrix{Float64}`: matrix of phase values, shape (channels, frequencies)
- `clabels::Vector{String}=string.(1:size(ph, 1))`: channel labels
- `db::Bool=true`: whether powers are normalized to dB
- `flim::Tuple{Real, Real}=(f[1], f[end]): frequency limits
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `zlabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `variant::Symbol=:w`: 3D visualization type:
    - `:w`: waterfall plot
    - `:s`: surface plot

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_phsd_3d(
    f::Vector{Float64},
    ph::Matrix{Float64};
    clabels::Vector{String} = string.(1:size(ph, 1)),
    db::Bool = true,
    flim::Tuple{Real, Real} = (f[1], f[end]),
    xlabel::String = "",
    ylabel::String = "",
    zlabel::String = "",
    title::String = "",
    mono::Bool = false,
    frq::Symbol = :lin,
    variant::Symbol,
)::GLMakie.Figure
    # validate
    _check_var(variant, [:w, :s], "variant")
    size(ph, 2) == length(f) || throw(
        ArgumentError("Length of powers vector must equal length of frequencies vector."),
    )
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(f), "flim")

    # set color palette
    pal = mono ? :grays : :darktest

    # number of channels
    ch_n = size(ph, 1)

    # prepare log-scaled frequencies axis
    if frq === :log && flim[1] == 0
        _warn("Currently log scale is not supported by Makie.")
        _warn("Lower frequency bound truncated to $(f[2]) Hz.")
        flim = (f[2], flim[2])
    end

    # frequency limits
    f1 = vsearch(flim[1], f)
    f2 = vsearch(flim[2], f)

    yts = ch_n > 64 ? 5 : ch_n > 32 ? 2 : 1

    # prepare plot
    GLMakie.activate!(; title = "plot_psd_3d()")
    plot_size = (900, 450)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis3(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        zlabel = zlabel,
        title = title,
        xticks = LinearTicks(15),
        # xminorticksvisible=true,
        # xminorticks=IntervalsBetween(10),
        # xscale=frq === :lin ? identity : log,
        yticks = (1:yts:ch_n, clabels[1:yts:end]),
        zoommode = :disable,
        xtranslationlock = true,
        ytranslationlock = true,
        ztranslationlock = true,
        aspect = (1, 1, 0.5),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0.1, 0.1),
        zautolimitmargin = (0, 0),
    )
    GLMakie.autolimits!(ax)
    _style_axis!(ax)

    # plot powers
    if variant === :w

        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            GLMakie.lines!(
                f,
                ones(length(f)) .* idx,
                ph[idx, f1:f2];
                linewidth = 2,
                color = mono ? :black : cmap[idx],
                colormap = pal,
                colorrange = 1:ch_n,
            )
        end

    elseif variant === :s

        # plot powers
        cmap = GLMakie.resample_cmap(pal, ch_n)
        GLMakie.surface!(f, eachindex(clabels), ph[:, f1:f2]'; colormap = pal)

    end

    return fig
end

"""
    plot_phsd_topo(locs, f, ph; <keyword arguments>)

Plot a topographical map of Phase Spectral Density (PHSD) across channel locations.

# Arguments

- `locs::DataFrame`: channel location data
- `f::Vector{Float64}`: vector of frequency values in Hz
- `ph::Matrix{Float64}`: matrix of phase values, shape (channels, frequencies)
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plots
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: if `true`, draw head outline

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_phsd_topo(
    locs::DataFrame,
    f::Vector{Float64},
    ph::Matrix{Float64};
    flim::Tuple{Real, Real} = (f[1], f[end]),
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    frq::Symbol = :lin,
    cart::Bool = false,
    head::Bool = true,
)::GLMakie.Figure
    # validate
    size(ph, 2) == length(f) || throw(
        ArgumentError("Length of powers vector must equal length of frequencies vector."),
    )
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(f), "flim")

    # prepare log-scaled frequencies axis
    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(f[2]) Hz.")
        flim = (f[2], flim[2])
    end

    # frequency limits
    f1 = vsearch(flim[1], f)
    f2 = vsearch(flim[2], f)

    # plot parameters
    ch_n = size(ph, 1)
    if ch_n <= 64
        plot_size   = (1000, 1000)
        marker_size = (150, 75)
        xl = 1.2
        yl = 1.2
    elseif ch_n <= 100
        plot_size   = (1200, 1200)
        marker_size = (110, 55)
        xl = 1.5
        yl = 1.5
    else
        plot_size   = (1400, 1400)
        marker_size = (90, 45)
        xl = 1.5
        yl = 1.5
    end

    # get locations
    if !cart
        loc_x = zeros(DataFrames.nrow(locs))
        loc_y = zeros(DataFrames.nrow(locs))
        for idx in axes(locs, 1)
            loc_x[idx], loc_y[idx] =
                pol2cart(locs.loc_radius[idx], locs.loc_theta[idx])
        end
    else
        loc_x = locs.loc_x
        loc_y = locs.loc_y
    end

    # prepare PHSD plots
    fig_vec      = GLMakie.Figure[]
    fig_full_vec = GLMakie.Figure[]
    for idx in axes(p, 1)
        fig_mini = GLMakie.Figure(; size = marker_size, figure_padding = 0)
        ax = GLMakie.Axis(
            fig_mini[1, 1];
            xlabel           = "",
            ylabel           = "",
            title            = locs[idx, :label],
            xscale           = frq === :lin ? identity : log,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0.1, 0.1),
        )
        hidedecorations!(ax)
        GLMakie.autolimits!(ax)
        ax.titlesize = 8
        GLMakie.lines!(ax, f, ph[idx, f1:f2]; linewidth = 1, color = :black)
        push!(fig_vec, fig_mini)

        fig_full = plot_psd(
            f, ph[idx, f1:f2];
            xlabel = xlabel,
            ylabel = ylabel,
            title  = locs[idx, :label] * ": " * title,
            flim   = flim,
            frq    = frq,
        )
        push!(fig_full_vec, fig_full)
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_phsd_topo()")
    fig = GLMakie.Figure(; size = plot_size, figure_padding = 0)

    # create axis with customizable properties
    ax  = GLMakie.Axis(
        fig[1, 1];
        xlabel = "",
        ylabel = "",
        title  = title,
        aspect = 1,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, (-xl, xl))
    GLMakie.ylims!(ax, (-yl, yl))
    hidespines!(ax)
    hidedecorations!(ax)
    ax.titlesize = 18

    # draw head outline
    head && _draw_head_outline!(ax; lw = 3)

    for idx in axes(p, 1)
        io = IOBuffer()
        show(io, MIME"image/png"(), fig_vec[idx])
        pp = FileIO.load(io)
        GLMakie.scatter!(
            loc_x[idx],
            loc_y[idx];
            marker = pp,
            markersize = marker_size,
            markerspace = :pixel,
        )
    end

    # PHSD positions
    loc_x_range = [(loc_x[i] - 0.15, loc_x[i] + 0.15) for i in eachindex(loc_x)]
    loc_y_range = [(loc_y[i] - 0.1,  loc_y[i] + 0.1)  for i in eachindex(loc_y)]

    # mouse events
    on(events(fig).mousebutton) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                ax_x = mouseposition(ax)[1]
                ax_y = mouseposition(ax)[2]
                for idx in eachindex(loc_x)
                    if ax_x >= loc_x_range[idx][1] &&
                       ax_x <= loc_x_range[idx][2] &&
                       ax_y >= loc_y_range[idx][1] &&
                       ax_y <= loc_y_range[idx][2]
                        display(GLMakie.Screen(), fig_full_vec[idx])
                        break
                    end
                end
            end
        end
    end

    return fig
end


"""
    plot_phsd(obj; <keyword arguments>)

Plot Phase Spectral Density (PHSD) using various estimation methods with customizable visualization.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `seg::Tuple{Real, Real}=(0, 10)`: time segment to analyze (from, to) in seconds; default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits for the plots
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is Phase [rad]
- `zlabel::String="default"`: z-axis label for 3-d plots, default is Phase [rad]
- `title::String="default"`: plot title, default is PHSD [frequency limit: 0-128 Hz] [epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: if `true`, use a monochrome palette
- `type::Symbol=:normal`: plot type:
    - `:normal` single channel or butterfly for multichannel
    - `:w3d`: 3-d waterfall
    - `:s3d`: 3-d surface
    - `:topo`: topographical
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: if `true`, draw head outline
- `leg::Bool=true`: if `true`, add legend with channel labels
- `avg::Bool=false`: if `true`, plot averaged PSD
- `ci95::Bool=false`: if `true`, plot mean and ±95% CI of averaged PSDs

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_phsd(
    obj::NeuroAnalyzer.NEURO;
    seg::Tuple{Real, Real} = (0, 10),
    ep::Int64 = 0,
    ch::Union{String, Vector{String}, Regex} = "all",
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    frq::Symbol = :lin,
    xlabel::String = "default",
    ylabel::String = "default",
    zlabel::String = "default",
    title::String = "default",
    mono::Bool = false,
    type::Symbol = :normal,
    cart::Bool = false,
    head::Bool = true,
    leg::Bool = true,
    avg::Bool = false,
    ci95::Bool = false,
)::GLMakie.Figure
    # validate
    _check_var(type, [:normal, :w3d, :s3d, :topo], "type")
    _check_var(frq, [:lin, :log], "frq")
    avg && ci95 && throw(ArgumentError("avg and ci95 cannot both be true."))

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    length(ch) == 1 && (ch = ch[1])

    # get signal for specified epochs
    if nepochs(obj) == 1
        ep == 0 || throw(ArgumentError("For continuous object, ep must not be specified."))
        if obj.time_pts[end] < 10 && seg == (0, 10)
            seg = (0, obj.time_pts[end])
        else
            _check_segment(obj, seg)
        end
        seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))
        signal = @view(obj.data[ch, seg[1]:seg[2], 1])
        t = obj.time_pts[seg[1]:seg[2]]
        _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    else
        ep != 0 || throw(ArgumentError("For epoched object, ep must be specified."))
        t = obj.epoch_time
        _check_epochs(obj, ep)
        signal = @view(obj.data[ch, :, ep])
    end

    # channel labels
    clabels = labels(obj)[ch]

    # frequency limits
    fs = sr(obj)
    _check_tuple(flim, (0, sr(obj) / 2), "flim")

    # calculate PHSD
    phsd_data = phsd(signal; fs = fs)
    sp, sf = phsd_data.sp, phsd_data.sf

    # factor out repeated title suffix pattern
    ep_suffix = ep != 0 ? "\n[epoch: $ep]" : "\n[time window: $t_s1:$t_s2]"
    title == "default" && (title = "PHSD$ep_suffix")

    if type === :normal
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Phase [rad]")
        if length(ch) == 1
            fig = plot_phsd(
                sf, sp;
                xlabel = xlabel,
                ylabel = ylabel,
                title  = title,
                flim   = flim,
                frq    = frq,
            )
        else
            fig = plot_phsd(
                sf, sp;
                xlabel   = xlabel,
                ylabel   = "",
                clabels  = clabels,
                title    = title,
                flim     = flim,
                frq      = frq,
                avg      = avg,
                ci95     = ci95,
                leg      = leg,
                mono     = mono,
            )
        end

    elseif type in [:w3d, :s3d]
        ndims(sp) >= 2 ||
            throw(ArgumentError("For type=:$type plot the signal must contain ≥ 2 channels."))
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "")
        zlabel == "default" && (zlabel = "Phase [rad]")
        fig = plot_phsd_3d(
            sf, sp;
            clabels = clabels,
            xlabel  = xlabel,
            ylabel  = ylabel,
            zlabel  = zlabel,
            title   = title,
            flim    = flim,
            frq     = frq,
            mono    = mono,
            variant = type === :w3d ? :w : :s,
        )

    elseif type === :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Phase [rad]")
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        length(unique(obj.header.recording[:channel_type][ch])) == 1 ||
            throw(ArgumentError(
                "For multi-channel topo plot all channels must be of the same type.",
            ))
        _has_locs(obj)
        chs  = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        ndims(sp) == 1 && (sp = reshape(sp, 1, length(sp)))
        fig = plot_phsd_topo(
            locs, sf, sp;
            xlabel = xlabel,
            ylabel = ylabel,
            title  = title,
            flim   = flim,
            frq    = frq,
            cart   = cart,
            head   = head,
        )

    end

    return fig
end