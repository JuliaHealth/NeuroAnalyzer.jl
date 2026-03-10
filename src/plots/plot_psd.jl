export plot_psd
export plot_psd_3d
export plot_psd_topo

"""
    plot_psd(sf, sp; <keyword arguments>)

Plot PSD (power spectrum density).

# Arguments

  - `sf::Vector{Float64}`: frequencies
  - `sp::Vector{Float64}`: powers
  - `flim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limits
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `frq::Symbol=:lin`: frequency scaling - `:lin` or `:log`

# Returns

  - `p::GLMakie.Figure`
"""
function plot_psd(
        sf::Vector{Float64},
        sp::Vector{Float64};
        flim::Tuple{Real, Real} = (sf[1], sf[end]),
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        frq::Symbol = :lin,
    )::GLMakie.Figure

    @assert length(sp) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        flim = (sf[2], flim[2])
    end

    f1 = vsearch(flim[1], sf)
    f2 = vsearch(flim[2], sf)

    # prepare plot
    GLMakie.activate!(title = "plot_psd()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(5),
        xscale = frq === :lin ? identity : log,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0.1, 0.1),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.autolimits!(ax)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # draw powers
    Makie.lines!(ax,
                 sf[f1:f2],
                 sp[f1:f2],
                 linewidth = 2,
                 color = :black)

    return p

end

"""
    plot_psd(sf, sp; <keyword arguments>)

Plot multi-channel PSD (power spectrum density).

# Arguments

  - `sf::Vector{Float64}`: frequencies
  - `sp::Matrix{Float64}`: powers
  - `clabels::Vector{String}=string.(1:size(sp, 1))`: channel labels
  - `flim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limits
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette
  - `frq::Symbol=:lin`: frequency scaling - `:lin` or `:log`
  - `avg::Bool=false`: if true, plot averaged PSD
  - `ci95::Bool=false`: if true, plot mean and ±95% CI of averaged PSDs
  - `leg::Bool=true`: if true, add legend with channel labels

# Returns

  - `p::GLMakie.Figure`
"""
function plot_psd(
        sf::Vector{Float64},
        sp::Matrix{Float64};
        clabels::Vector{String} = string.(1:size(sp, 1)),
        flim::Tuple{Real, Real} = (sf[1], sf[end]),
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        mono::Bool = false,
        frq::Symbol = :lin,
        avg::Bool = false,
        ci95::Bool = false,
        leg::Bool = true,
    )::GLMakie.Figure

    ch_n = size(sp, 1)

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    pal = mono ? :grays : :darktest

    f1 = vsearch(flim[1], sf)
    f2 = vsearch(flim[2], sf)

    # get mean and 95%CI
    if ci95
        s_m, _, s_u, s_l = NeuroAnalyzer.msci95(sp[f1:f2])
    end

    # prepare plot
    GLMakie.activate!(title = "plot_psd()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(15),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xscale = frq === :lin ? identity : log,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0.1, 0.1),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.autolimits!(ax)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    if ci95
        # draw 95% CI
        Makie.band!(ax, sf[f1:f2], s_u, s_l; alpha = 0.25, color = :grey, strokewidth = 0.5)

        # draw mean
        Makie.lines!(ax, sf[f1:f2], s_m; color = :black, linewidth = 2)
    else
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            Makie.lines!(
                ax,
                sf[f1:f2],
                sp[idx, f1:f2];
                color = cmap[idx],
                colormap = pal,
                colorrange = 1:ch_n,
                linewidth = 2,
                label = clabels[idx],
            )
        end

        # draw averaged channels
        if avg
            s = mean(sp[f1:f2], dims = 1)[:]
            Makie.lines!(
                ax,
                sf[f1:f2],
                s,
                colormap = pal,
                linewidth = 4,
                color = :black,
            )
        end

        (leg && ch_n < 30) && axislegend(; position = :rt, colormap = pal)

    end

    return p

end

"""
    plot_psd_3d(sf, sp; <keyword arguments>)

Plot 3-d PSD (power spectrum density).

# Arguments

  - `sf::Vector{Float64}`: frequencies
  - `sp::Array{Float64, 3}`: powers
  - `clabels::Vector{String}=string.(1:size(sp, 1))`: channel labels
  - `db::Bool=true`: whether powers are normalized to dB
  - `flim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limits
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `zlabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette
  - `frq::Symbol=:lin`: frequency scaling - `:lin` or `:log`
  - `variant::Symbol`: waterfall (`:w`) or surface (`:s`)

# Returns

  - `p::GLMakie.Figure`
"""
function plot_psd_3d(
        sf::Vector{Float64},
        sp::Matrix{Float64};
        clabels::Vector{String} = string.(1:size(sp, 1)),
        db::Bool = true,
        flim::Tuple{Real, Real} = (sf[1], sf[end]),
        xlabel::String = "",
        ylabel::String = "",
        zlabel::String = "",
        title::String = "",
        mono::Bool = false,
        frq::Symbol = :lin,
        variant::Symbol,
    )::GLMakie.Figure

    _check_var(variant, [:w, :s], "variant")
    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    ch_n = size(sp, 1)

    pal = mono ? :grays : :darktest

    if frq === :log && flim[1] == 0
        _warn("Currently log scale is not supported by Makie.")
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        flim = (sf[2], flim[2])
    end

    if ch_n > 64
        yts = 5
    elseif ch_n > 32
        yts = 2
    else
        yts = 1
    end

    # prepare plot
    GLMakie.activate!(title = "plot_psd()")
    if variant === :w
        plot_size = (900, 450)
        p = GLMakie.Figure(size = plot_size)
        ax = GLMakie.Axis3(
            p[1, 1];
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
        GLMakie.xlims!(ax, flim)
        ax.titlesize = 18
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        # plot powers
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            Makie.lines!(
                sf,
                ones(length(sf)) .* idx,
                sp[idx, :];
                linewidth = 2,
                color = mono ? :black : cmap[idx],
                colormap = pal,
                colorrange = 1:ch_n,
            )
        end
    else
        f1 = vsearch(flim[1], sf)
        f2 = vsearch(flim[2], sf)
        plot_size = (900, 450)
        p = GLMakie.Figure(size = plot_size)
        ax = GLMakie.Axis3(
            p[1, 1];
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
            yautolimitmargin = (0, 0),
            zautolimitmargin = (0, 0),
        )
        ax.titlesize = 18
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        # plot powers
        cmap = GLMakie.resample_cmap(pal, ch_n)
        Makie.surface!(sf[f1:f2], eachindex(clabels), sp[:, f1:f2]'; colormap = pal)
    end

    return p

end

"""
    plot_psd_topo(locs, sf, sp; <keyword arguments>)

Plot topographical map of PSDs (power spectrum density).

# Arguments

  - `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
  - `sf::Vector{Float64}`: frequencies
  - `sp::Matrix{Float64}`: powers
  - `flim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limits
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `frq::Symbol=:lin`: frequency scaling - `:lin` or `:log`
  - `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
  - `head::Bool=true`: plot head shape

# Returns

  - `p::GLMakie.Figure`
"""
function plot_psd_topo(
        locs::DataFrame,
        sf::Vector{Float64},
        sp::Matrix{Float64};
        flim::Tuple{Real, Real} = (sf[1], sf[end]),
        title::String = "",
        xlabel::String = "",
        ylabel::String = "",
        frq::Symbol = :lin,
        cart::Bool = false,
        head::Bool = true,
    )::GLMakie.Figure

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        flim = (sf[2], flim[2])
    end

    pos = collect(1:DataFrames.nrow(locs))

    # plot parameters
    if size(sp, 1) <= 64
        plot_size = (1000, 1000)
        marker_size = (150, 75)
        xl = 1.2
        yl = 1.2
    elseif _in(size(sp, 1), (64, 100))
        plot_size = (1200, 1200)
        marker_size = (110, 55)
        xl = 1.5
        yl = 1.5
    else
        plot_size = (1400, 1400)
        marker_size = (90, 45)
        xl = 1.5
        yl = 1.5
    end

    # get locations
    if !cart
        loc_x = zeros(size(locs, 1))
        loc_y = zeros(size(locs, 1))
        for idx in axes(locs, 1)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
    else
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
    end

    # prepare PSD plots
    pp_vec = GLMakie.Figure[]
    pp_full_vec = GLMakie.Figure[]
    for idx in axes(sp, 1)
        pp = GLMakie.Figure(
            size = marker_size,
            figure_padding = 0,
        )
        ax = GLMakie.Axis(
            pp[1, 1];
            xlabel = "",
            ylabel = "",
            title = locs[idx, :label],
            xscale = frq === :lin ? identity : log,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0.1, 0.1),
        )
        hidedecorations!(ax)
        GLMakie.xlims!(ax, flim)
        ax.titlesize = 8
        # plot powers
        GLMakie.lines!(ax, sf, sp[idx, :]; linewidth = 1, color = :black)
        push!(pp_vec, pp)
        pp_full = plot_psd(
            sf,
            sp[idx, :];
            xlabel = xlabel,
            ylabel = ylabel,
            title = locs[idx, :label] * ": " * title,
            flim = flim,
            frq = frq,
        )
        push!(pp_full_vec, pp_full)
    end

    # prepare plot
    GLMakie.activate!(title = "plot_psd()")
    p = GLMakie.Figure(
        size = plot_size,
        figure_padding = 0,
    )
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = "",
        ylabel = "",
        title = title,
        aspect = 1,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.xlims!(ax, (-xl, xl))
    GLMakie.ylims!(ax, (-yl, yl))
    hidespines!(ax)
    hidedecorations!(ax)
    ax.titlesize = 18

    if head
        # nose
        GLMakie.lines!(ax, [-0.2, 0], [0.98, 1.08]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [0.2, 0], [0.98, 1.08]; linewidth = 3, color = :black)

        # ears
        # left
        GLMakie.lines!(ax, [-0.995, -1.03], [0.1, 0.15]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.03, -1.06], [0.15, 0.16]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.06, -1.1], [0.16, 0.14]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.1, -1.12], [0.14, 0.05]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.12, -1.1], [0.05, -0.1]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.1, -1.13], [-0.1, -0.3]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.13, -1.09], [-0.3, -0.37]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.09, -1.02], [-0.37, -0.39]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.02, -0.98], [-0.39, -0.33]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-0.98, -0.975], [-0.33, -0.22]; linewidth = 3, color = :black)
        # right
        GLMakie.lines!(ax, [0.995, 1.03], [0.1, 0.15]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.03, 1.06], [0.15, 0.16]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.06, 1.1], [0.16, 0.14]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.1, 1.12], [0.14, 0.05]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.12, 1.1], [0.05, -0.1]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.1, 1.13], [-0.1, -0.3]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.13, 1.09], [-0.3, -0.37]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.09, 1.02], [-0.37, -0.39]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.02, 0.98], [-0.39, -0.33]; linewidth = 3, color = :black)
        GLMakie.lines!(ax, [0.98, 0.975], [-0.33, -0.22]; linewidth = 3, color = :black)

        # head
        GLMakie.arc!(ax, (0, 0), 1, 0, 2pi; linewidth = 3, color = :black)
    end

    for idx in axes(sp, 1)
        io = IOBuffer()
        show(io, MIME"image/png"(), pp_vec[idx])
        pp = FileIO.load(io)
        GLMakie.scatter!(loc_x[idx], loc_y[idx]; marker = pp, markersize = marker_size, markerspace = :pixel)
    end

    loc_x_range = Tuple{Float64, Float64}[]
    loc_y_range = Tuple{Float64, Float64}[]
    for idx in eachindex(loc_x)
        push!(loc_x_range, (loc_x[idx] - 0.15, loc_x[idx] + 0.15))
        push!(loc_y_range, (loc_y[idx] - 0.1, loc_y[idx] + 0.1))
    end
    on(events(p).mousebutton) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                ax_x = mouseposition(ax)[1]
                ax_y = mouseposition(ax)[2]
                for idx in eachindex(loc_x)
                    if ax_x >= loc_x_range[idx][1] &&
                            ax_x <= loc_x_range[idx][2] &&
                            ax_y >= loc_y_range[idx][1] &&
                            ax_y <= loc_y_range[idx][2]
                        display(GLMakie.Screen(), pp_full_vec[idx])
                        break
                    end
                end
            end
        end
    end

    return p

end

"""
    plot_psd(obj; <keyword arguments>)

Plot PSD (power spectrum density).

# Arguments

  - `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  - `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  - `ep::Int64=0`: epoch to display
  - `ch::Union{String, Vector{String}, Regex}=datatype(obj)`: channel name or list of channel names
  - `db::Bool=true`: normalize powers to dB
  - `method::Symbol=:welch`: PSD method:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-taper periodogram
      + `:stft`: short-time Fourier transform
      + `:mw`: Morlet wavelet convolution
      + `:gh`: Gaussian and Hilbert transform
  - `nt::Int64=7`: number of Slepian tapers
  - `wlen::Int64=fs`: window length in samples, default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
  - `w::Bool=true`: if true, apply Hanning window
  - `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
  - `gw::Real=5`: Gaussian width in Hz
  - `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher`
  - `demean::Bool=true`: subtract DC before calculating PSD
  - `frq::Symbol=:lin`: frequency scaling - `:lin` or `:log`
  - `xlabel::String="default"`: x-axis label
  - `ylabel::String="default"`: y-axis label
  - `zlabel::String="default"`: z-axis label for 3-d plots
  - `title::String="default"`: plot title
  - `mono::Bool=false`: use color or gray palette
  - `type::Symbol=:normal`: plot type:
      + `:normal` single channel or butterfly for multichannel
      + `:w3d`: 3-d waterfall
      + `:s3d`: 3-d surface
      + `:topo`: topographical
  - `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
  - `head::Bool=true`: plot head shape
  - `leg::Bool=true`: if true, add legend with channel labels
  - `avg::Bool=false`: if true, plot averaged PSD
  - `ci95::Bool=false`: if true, plot mean and ±95% CI of averaged PSDs

# Returns

  - `p::GLMakie.Figure`
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

    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    length(ch) == 1 && (ch = ch[1])

    if nepochs(obj) == 1
        @assert ep == 0 "For continuous object, ep must not be specified."
        if obj.time_pts[end] < 10 && seg == (0, 10)
            seg = (0, obj.time_pts[end])
        else
            _check_segment(obj, seg)
        end
        seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))
        signal = @views obj.data[ch, seg[1]:seg[2], 1]
        t = obj.time_pts[seg[1]:seg[2]]
        _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    else
        @assert ep != 0 "For epoched object, ep must be specified."
        t = obj.epoch_time
        _check_epochs(obj, ep)
        signal = @views obj.data[ch, :, ep]
    end

    # channel labels
    clabels = labels(obj)[ch]

    # set units
    units = _ch_units(obj, labels(obj)[ch[1]])

    # frequency limits
    fs = sr(obj)
    ref !== :abs && (flim = band_frq(obj, band = ref))
    _check_tuple(flim, (0, sr(obj) / 2), "flim")

    # calculate PSD
    if ref === :abs
        if method === :welch
            sp, sf = psd(
                        signal,
                        fs = fs,
                        db = db,
                        method = :welch,
                        wlen = wlen,
                        woverlap = woverlap,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "Absolute PSD (Welch's periodogram)\n[epoch: $ep]")
            else
                title == "default" && (title = "Absolute PSD (Welch's periodogram)\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :fft
            sp, sf = psd(
                        signal,
                        fs = fs,
                        db = db,
                        method = :fft,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "Absolute PSD (fast Fourier transform)\n[epoch: $ep]")
            else
                title == "default" && (title = "Absolute PSD (fast Fourier transform)\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :stft
            sp, sf = psd(
                        signal,
                        fs = fs,
                        db = db,
                        method = :stft,
                        wlen = wlen,
                        woverlap = woverlap,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "Absolute PSD (short-time Fourier transform)\n[epoch: $ep]")
            else
                title == "default" &&
                    (title = "Absolute PSD (short-time Fourier transform)\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :mt
            sp, sf = psd(
                        signal,
                        fs = fs,
                        db = db,
                        method = :mt,
                        nt = nt,
                        wlen = wlen,
                        woverlap = woverlap,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "Absolute PSD (multi-taper)\n[epoch: $ep]")
            else
                title == "default" && (title = "Absolute PSD (multi-taper)\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :mw
            sp, sf = psd(
                        signal,
                        fs = fs,
                        db = db,
                        method = :mw,
                        ncyc = ncyc,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "Absolute PSD (Morlet wavelet convolution)\n[epoch: $ep]")
            else
                title == "default" && (title = "Absolute PSD (Morlet wavelet convolution)\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :gh
            sp, sf = psd(
                        signal,
                        fs = fs,
                        db = db,
                        method = :gh,
                        gw = gw,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "Absolute PSD (Gaussian and Hilbert transform)\n[epoch: $ep]")
            else
                title == "default" &&
                    (title = "Absolute PSD (Gaussian and Hilbert transform)\n[time window: $t_s1:$t_s2]")
            end
        end
    elseif ref === :total
        if method === :welch
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :welch,
                        wlen = wlen,
                        woverlap = woverlap,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "PSD (Welch's periodogram) relative to total power\n[epoch: $ep]")
            else
                title == "default" &&
                    (title = "PSD (Welch's periodogram) relative to total power\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :fft
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :fft,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "PSD (fast Fourier transform) relative to total power\n[epoch: $ep]")
            else
                title == "default" &&
                    (title = "PSD (fast Fourier transform) relative to total power\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :stft
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :stft,
                        wlen = wlen,
                        woverlap = woverlap,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" &&
                    (title = "PSD (short-time Fourier transform) relative to total power\n[epoch: $ep]")
            else
                title == "default" &&
                    (title = "PSD (short-time Fourier transform) relative to total power\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :mt
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :mt,
                        nt = nt,
                        wlen = wlen,
                        woverlap = woverlap,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "PSD (multi-taper) relative to total power\n[epoch: $ep]")
            else
                title == "default" && (title = "PSD (multi-taper) relative to total power\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :mw
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :mw,
                        ncyc = ncyc,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (title = "PSD (Morlet wavelet convolution) relative to total power\n[epoch: $ep]")
            else
                title == "default" &&
                    (title = "PSD (Morlet wavelet convolution) relative to total power\n[time window: $t_s1:$t_s2]")
            end
        elseif method === :gh
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :gh,
                        gw = gw,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" &&
                    (title = "PSD (Gaussian and Hilbert transform) relative to total power\n[epoch: $ep]")
            else
                title == "default" &&
                    (title = "PSD (Gaussian and Hilbert transform) relative to total power\n[time window: $t_s1:$t_s2]")
            end
        end
    else
        if method === :welch
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :welch,
                        flim = flim,
                        wlen = wlen,
                        woverlap = woverlap,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (
                    title = "PSD (Welch's periodogram) relative to $(replace(string(ref), "_" => " ")) power\n[epoch: $ep]"
                )
            else
                title == "default" && (
                    title = "PSD (Welch's periodogram) relative to $(replace(string(ref), "_" => " ")) power\n[time window: $t_s1:$t_s2]"
                )
            end
        elseif method === :fft
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :fft,
                        flim = flim,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (
                    title = "PSD (fast Fourier transform) relative to $(replace(string(ref), "_" => " ")) power\n[epoch: $ep]"
                )
            else
                title == "default" && (
                    title = "PSD (fast Fourier transform) relative to $(replace(string(ref), "_" => " ")) power\n[time window: $t_s1:$t_s2]"
                )
            end
        elseif method === :stft
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :stft,
                        flim = flim,
                        wlen = wlen,
                        woverlap = woverlap,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (
                    title = "PSD (short-time Fourier transform) relative to $(replace(string(ref), "_" => " ")) power\n[epoch: $ep]"
                )
            else
                title == "default" && (
                    title = "PSD (short-time Fourier transform) relative to $(replace(string(ref), "_" => " ")) power\n[time window: $t_s1:$t_s2]"
                )
            end
        elseif method === :mt
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :mt,
                        flim = flim,
                        nt = nt,
                        wlen = wlen,
                        woverlap = woverlap,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" &&
                    (title = "PSD (multi-taper) relative to $(replace(string(ref), "_" => " ")) power\n[epoch: $ep]")
            else
                title == "default" && (
                    title = "PSD (multi-taper) relative to $(replace(string(ref), "_" => " ")) power\n[time window: $t_s1:$t_s2]"
                )
            end
        elseif method === :mw
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :mw,
                        flim = flim,
                        ncyc = ncyc,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (
                    title = "PSD (Morlet wavelet convolution) relative to $(replace(string(ref), "_" => " ")) power\n[epoch: $ep]"
                )
            else
                title == "default" && (
                    title = "PSD (Morlet wavelet convolution) relative to $(replace(string(ref), "_" => " ")) power\n[time window: $t_s1:$t_s2]"
                )
            end
        elseif method === :gh
            sp, sf = psd_rel(
                        signal,
                        fs = fs,
                        db = db,
                        method = :gh,
                        flim = flim,
                        gw = gw,
                        w = w,
                        demean = demean,
                    )
            if ep != 0
                title == "default" && (
                    title = "PSD (Gaussian and Hilbert transform) relative to $(replace(string(ref), "_" => " ")) power\n[epoch: $ep]"
                )
            else
                title == "default" && (
                    title = "PSD (Gaussian and Hilbert transform) relative to $(replace(string(ref), "_" => " ")) power\n[time window: $t_s1:$t_s2]"
                )
            end
        end
    end

    if type === :normal
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if ref !== :abs
            ylabel == "default" && (ylabel = "Power ratio")
        else
            ylabel == "default" && (ylabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        end
        if length(ch) == 1
            p = plot_psd(sf, sp; xlabel = xlabel, ylabel = ylabel, title = title, flim = flim, frq = frq)
        else
            p = plot_psd(
                sf,
                sp,
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
    elseif type === :w3d || type === :s3d
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "")
        zlabel == "default" && (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        ch_t = obj.header.recording[:channel_type]
        p = plot_psd_3d(
            sf,
            sp,
            clabels = clabels,
            xlabel = xlabel,
            ylabel = ylabel,
            zlabel = zlabel,
            title = title,
            flim = flim,
            frq = frq,
            mono = mono,
            variant = type === :w3d ? :w : :s,
        )
    elseif type === :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        @assert length(unique(obj.header.recording[:channel_type][ch])) == 1 "For multi-channel topo plot all channels must be of the same type."
        _has_locs(obj)
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        ndims(sp) == 1 && (sp = reshape(sp, 1, length(sp)))
        p = plot_psd_topo(
            locs,
            sf,
            sp,
            xlabel = xlabel,
            ylabel = ylabel,
            title = title,
            flim = flim,
            frq = frq,
            cart = cart,
            head = head,
        )
    end

    return p

end
