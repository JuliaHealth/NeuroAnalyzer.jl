export plot_psd
export plot_psd_3d
export plot_psd_topo

"""
    plot_psd(f, p; <keyword arguments>)

Plot the Power Spectral Density (PSD) with customizable visualization options.

# Arguments

- `f::Vector{Float64}`: vector of frequency values in Hz
- `p::Vector{Float64}`: vector of power spectral density values
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plots
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_psd(
    f::Vector{Float64},
    p::Vector{Float64};
    flim::Tuple{Real, Real} = (f[1], f[end]),
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    frq::Symbol = :lin,
)::GLMakie.Figure
    # validate
    length(p) == length(f) || throw(
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
    GLMakie.activate!(; title = "plot_psd()")
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

    GLMakie.lines!(
        ax,
        f[f1:f2],
        p[f1:f2];
        linewidth = 2,
        color = :black,
    )

    return fig
end

"""
    plot_psd(f, p; <keyword arguments>)

Plot multi-channel Power Spectral Density (PSD) with customizable visualization options.
# Arguments

- `f::Vector{Float64}`: vector of frequency values in Hz
- `p::Matrix{Float64}`: matrix of power spectral density values, shape (channels, frequencies)
- `clabels::Vector{String}=string.(1:size(p, 1))`: channel labels
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plots
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `avg::Bool=false`: if `true`, plot averaged PSD across channels
- `ci95::Bool=false`: if `true`, plot mean and ±95% confidence interval of averaged PSDs
- `leg::Bool=true`: if `true`, add legend with channel labels

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_psd(
    f::Vector{Float64},
    p::Matrix{Float64};
    clabels::Vector{String} = string.(1:size(p, 1)),
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
    size(p, 2) == length(f) || throw(
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
        msci95_data = NeuroAnalyzer.msci95(p[:, f1:f2])
        s_m = msci95_data.sm
        s_l = msci95_data.ll
        s_u = msci95_data.ul
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_psd()")
    plot_size = (900, 450)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis(
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
        GLMakie.band!(
            ax,
            f[f1:f2],
            s_u,
            s_l;
            alpha = 0.25,
            color = :grey,
            strokewidth = 0.5,
        )
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
    plot_psd_3d(f, p; <keyword arguments>)

Plot a 3D representation of multi-channel Power Spectral Density (PSD).

# Arguments

- `f::Vector{Float64}`: vector of frequency values in Hz
- `p::Matrix{Float64}`: matrix of power values, shape (channels, frequencies)
- `clabels::Vector{String}=string.(1:size(p, 1))`: channel labels
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
function plot_psd_3d(
    f::Vector{Float64},
    p::Matrix{Float64};
    clabels::Vector{String} = string.(1:size(p, 1)),
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
    size(p, 2) == length(f) || throw(
        ArgumentError("Length of powers vector must equal length of frequencies vector."),
    )
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(f), "flim")

    # set color palette
    pal = mono ? :grays : :darktest

    # number of channels
    ch_n = size(p, 1)

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
                p[idx, f1:f2];
                linewidth = 2,
                color = mono ? :black : cmap[idx],
                colormap = pal,
                colorrange = 1:ch_n,
            )
        end

    elseif variant === :s

        # plot powers
        cmap = GLMakie.resample_cmap(pal, ch_n)
        GLMakie.surface!(f, eachindex(clabels), p[:, f1f2]'; colormap = pal)
    end

    return fig
end

"""
    plot_psd_topo(locs, f, p; <keyword arguments>)

Plot a topographical map of Power Spectral Density (PSD) across channel locations.

# Arguments

- `locs::DataFrame`: channel location data
- `f::Vector{Float64}`: vector of frequency values in Hz
- `p::Matrix{Float64}`: matrix of power spectral density values, shape (channels, frequencies)
- `flim::Tuple{Real, Real}=(f[1], f[end]): frequency limits for the plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: if `true`, draw head outline

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_psd_topo(
    locs::DataFrame,
    f::Vector{Float64},
    p::Matrix{Float64};
    flim::Tuple{Real, Real} = (f[1], f[end]),
    title::String = "",
    xlabel::String = "",
    ylabel::String = "",
    frq::Symbol = :lin,
    cart::Bool = false,
    head::Bool = true,
)::GLMakie.Figure
    # validate
    size(p, 2) == length(f) || throw(
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
    ch_n = size(p, 1)
    if ch_n <= 64
        plot_size   = (1000, 1000)
        marker_size = (150, 75)
        xl          = 1.2
        yl          = 1.2
    elseif ch_n <= 100
        plot_size   = (1200, 1200)
        marker_size = (110, 55)
        xl          = 1.5
        yl          = 1.5
    else
        plot_size   = (1400, 1400)
        marker_size = (90, 45)
        xl          = 1.5
        yl          = 1.5
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

    # prepare PSD plots
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
        GLMakie.lines!(ax, f, p[idx, f1:f2]; linewidth = 1, color = :black)
        push!(fig_vec, fig_mini)

        fig_full = plot_psd(
            f, p[idx, f1:f2];
            xlabel = xlabel,
            ylabel = ylabel,
            title  = locs[idx, :label] * ": " * title,
            flim   = flim,
            frq    = frq,
        )
        push!(fig_full_vec, fig_full)
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_psd_topo()")
    fig = GLMakie.Figure(; size = plot_size, figure_padding = 0)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = "",
        ylabel = "",
        title = title,
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

    # PSD positions
    loc_x_range = [(loc_x[i] - 0.15, loc_x[i] + 0.15) for i in eachindex(loc_x)]
    loc_y_range = [(loc_y[i] - 0.1, loc_y[i] + 0.1) for i in eachindex(loc_y)]

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
    plot_psd(obj; <keyword arguments>)

Plot Power Spectral Density (PSD) using various estimation methods with customizable visualization.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `seg::Tuple{Real, Real}=(0, 10)`: time segment to analyze (from, to) in seconds; default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}=datatype(obj)`: channel name(s)
- `db::Bool=true`: if `true`, normalize powers to dB
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-taper periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits for the plots
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `ref::Symbol=:abs`: PSD reference type:
    - `:abs`: absolute power (no reference)
    - `:total`: relative to total power
    - `:delta`: relative to delta band power
    - `:theta`: relative to theta band power
    - `:alpha`: relative to alpha band power
    - `:beta`: relative to beta band power
    - `:beta_high`: relative to high beta band power
    - `:gamma`: relative to gamma band power
    - `:gamma_1`: relative to gamma-1 band power
    - `:gamma_2`: relative to gamma-2 band power
    - `:gamma_lower`: relative to lower gamma band power
    - `:gamma_higher`: relative to higher gamma band power
- `demean::Bool=true`: subtract DC component before estimating PSD
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `zlabel::String="default"`: z-axis label for 3-d plots
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `type::Symbol=:normal`: plot type:
    - `:normal` single channel or butterfly plot for multichannel
    - `:w3d`: 3D waterfall plot
    - `:s3d`: 3D surface plot
    - `:topo`: topographical plot
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: if `true`, draw head outline
- `leg::Bool=true`: if `true`, add legend with channel labels
- `avg::Bool=false`: if `true`, plot averaged PSD across channels
- `ci95::Bool=false`: if `true`, plot mean and ±95% confidence interval of averaged PSDs

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_psd(
    obj::NeuroAnalyzer.NEURO;
    seg::Tuple{Real, Real} = (0, 10),
    ep::Int64 = 0,
    ch::Union{String, Vector{String}, Regex} = datatype(obj),
    db::Bool = true,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    ref::Symbol = :abs,
    demean::Bool = true,
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
    _check_var(method, [:welch, :fft, :stft, :mt, :mw, :gh], "method")
    _check_var(
        ref,
        [
            :abs,
            :total,
            :delta,
            :theta,
            :alpha,
            :alpha_lower,
            :alpha_higher,
            :beta,
            :beta_lower,
            :beta_higher,
            :gamma,
            :gamma_1,
            :gamma_2,
            :gamma_lower,
            :gamma_higher,
        ],
        "ref",
    )
    _check_var(frq, [:lin, :log], "frq")

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
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

    # set units
    units = _ch_units(obj, labels(obj)[ch[1]])

    # frequency limits
    fs = sr(obj)
    ref !== :abs && (flim = band_frq(obj; band = ref))
    _check_tuple(flim, (0, sr(obj) / 2), "flim")

    # factor out repeated title fragments
    _method_label = Dict(
        :welch => "Welch's periodogram",
        :fft   => "fast Fourier transform",
        :stft  => "short-time Fourier transform",
        :mt    => "multi-taper",
        :mw    => "Morlet wavelet convolution",
        :gh    => "Gaussian and Hilbert transform",
    )
    method_label = _method_label[method]
    ep_suffix = ep != 0 ? "\n[epoch: $ep]" : "\n[time window: $t_s1:$t_s2]"
    ref_prefix = if ref === :abs
        "Absolute"
    elseif ref === :total
        "relative to total power"
    else
        "relative to $(replace(string(ref), "_" => " ")) power"
    end
    default_title =
        ref === :abs ?
        "Absolute PSD ($method_label)$ep_suffix" :
        "PSD ($method_label) $ref_prefix$ep_suffix"

    # common psd / psd_rel call arguments, varying only by method
    _common_kwargs = (
        fs = fs, db = db, method = method, w = w, demean = demean,
        wlen = wlen, woverlap = woverlap, nt = nt, ncyc = ncyc, gw = gw,
    )

    if ref === :abs
        p, f = psd(signal; _common_kwargs...)
    elseif ref === :total
        p, f = psd_rel(signal; _common_kwargs...)
    else
        p, f = psd_rel(signal; flim = flim, _common_kwargs...)
    end

    title == "default" && (title = default_title)

    if type === :normal
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (
            ylabel =
                ref !== :abs ?
                "Power ratio" :
                (db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        )

        if length(ch) == 1
            fig = plot_psd(
                f,
                p;
                xlabel = xlabel,
                ylabel = ylabel,
                title = title,
                flim = flim,
                frq = frq,
            )
        else
            fig = plot_psd(
                f,
                p;
                xlabel = xlabel,
                ylabel = "",
                clabels = clabels,
                title = title,
                flim = flim,
                frq = frq,
                avg = avg,
                ci95 = ci95,
                leg = leg,
                mono = mono,
            )
        end

    elseif type in [:w3d, :s3d]
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "")
        zlabel == "default" &&
            (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        return plot_psd_3d(
            f,
            p;
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
        ylabel == "default" &&
            (ylabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        length(unique(obj.header.recording[:channel_type][ch])) != 1 && throw(
            ArgumentError(
                "For multi-channel topo plot all channels must be of the same type.",
            ),
        )
        _has_locs(obj)
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        ndims(p) == 1 && (p = reshape(p, 1, length(p)))
        fig = plot_psd_topo(
            locs,
            f,
            p;
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
