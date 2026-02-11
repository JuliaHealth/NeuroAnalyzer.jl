export plot_phsd
export plot_phsd_3d
export plot_phsd_topo

"""
    plot_phsd(f, ph; <keyword arguments>)

Plot PHSD (phase spectral density).

# Arguments

- `f::Vector{Float64}`: frequencies
- `ph::Vector{Float64}`: phases
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling

# Returns

- `p::GLMakie.Figure`
"""
function plot_phsd(f::Vector{Float64}, ph::Vector{Float64}; flim::Tuple{Real, Real}=(f[1], f[end]), xlabel::String="", ylabel::String="", title::String="", frq::Symbol=:lin)::GLMakie.Figure

    @assert length(ph) == length(f) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, "flim")

    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(f[2]) Hz.")
        flim = (f[2], flim[2])
    end

    # prepare plot
    plot_size = (900, 450)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(15),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      xscale=frq === :lin ? identity : log,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0.1, 0.1),
                      xzoomlock=true,
                      yzoomlock=true,
                      xpanlock=true,
                      ypanlock=true,
                      xrectzoom=false,
                      yrectzoom=false)
    GLMakie.xlims!(ax, flim)
    GLMakie.ylims!(ax, extrema(ph))
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # draw powers
    Makie.lines!(f,
                 ph,
                 linewidth=2,
                 color=:black)

    return p

end

"""
    plot_phsd(f, ph; <keyword arguments>)

Plot multi-channel PHSD (phase spectral density).

# Arguments

- `f::Vector{Float64}`: frequencies
- `ph::Matrix{Float64}`: phases
- `clabels::Vector{String}=string.(1:size(sp, 1))`: channel labels
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `avg::Bool=false`: if true, plot averaged PHSD
- `ci95::Bool=false`: if true, plot mean and ±95% CI of averaged PHSDs
- `leg::Bool=true`: if true, add legend with channel labels

# Returns

- `p::GLMakie.Figure`
"""
function plot_phsd(f::Vector{Float64}, ph::Matrix{Float64}; clabels::Vector{String}=string.(1:size(ph, 1)), flim::Tuple{Real, Real}=(f[1], f[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, frq::Symbol=:lin, avg::Bool=false, ci95::Bool=false, leg::Bool=true)::GLMakie.Figure

    ch_n = size(ph, 1)

    @assert size(ph, 2) == length(f) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, "flim")

    pal = mono ? :grays : :darktest

    # get mean and 95%CI
    if ci95
        s_m, _, s_u, s_l = NeuroAnalyzer.msci95(ph)
    end

    # prepare plot
    plot_size = (900, 450)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(15),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      xscale=frq === :lin ? identity : log,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0.1, 0.1),
                      xzoomlock=true,
                      yzoomlock=true,
                      xpanlock=true,
                      ypanlock=true,
                      xrectzoom=false,
                      yrectzoom=false)
    GLMakie.xlims!(ax, flim)
    if ci95
        GLMakie.ylims!(ax, minimum(s_l), maximum(s_u))
    else
        GLMakie.ylims!(ax, extrema(ph))
    end
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    if ci95
        # draw 95% CI
        Makie.band!(f,
                    s_u,
                    s_l,
                    alpha=0.25,
                    color=:grey,
                    strokewidth=0.5)

        # draw mean
        Makie.lines!(f,
                     s_m,
                     color=:black,
                     linewidth=2)
    else
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            Makie.lines!(f,
                         ph[idx, :],
                         color=cmap[idx],
                         colormap=pal,
                         colorrange=1:ch_n,
                         linewidth=2,
                         label=clabels[idx])
        end

        # draw averaged channels
        if avg
            s = mean(ph, dims=1)[:]
            Makie.lines!(f,
                         s,
                         colormap=pal,
                         linewidth=4,
                         color=:black)
        end

        (leg && ch_n < 30) && axislegend(position=:rt,
                                         colormap=pal)

    end

    return p

end

"""
    plot_phsd_3d(f, ph; <keyword arguments>)

Plot 3-d PHSD (phase phectral density).

# Arguments

- `f::Vector{Float64}`: frequencies
- `ph::Array{Float64, 3}`: phases
- `clabels::Vector{String}=string.(1:size(ph, 1))`: channel labels
- `db::Bool=true`: whether powers are normalized to dB
- `flim::Tuple{Real, Real}=(f[1], f[end]): frequency limits
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `zlabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `variant::Symbol`: waterfall (`:w`) or surface (`:s`)

# Returns

- `p::GLMakie.Figure`
"""
function plot_phsd_3d(f::Vector{Float64}, ph::Matrix{Float64}; clabels::Vector{String}=string.(1:size(ph, 1)), db::Bool=true, flim::Tuple{Real, Real}=(f[1], f[end]), xlabel::String="", ylabel::String="", zlabel::String="", title::String="", mono::Bool=false, frq::Symbol=:lin, variant::Symbol)::GLMakie.Figure

    _check_var(variant, [:w, :s], "variant")
    @assert size(ph, 2) == length(f) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, "flim")

    ch_n = size(ph, 1)

    pal = mono ? :grays : :darktest

    if frq === :log && flim[1] == 0
        _warn("Currently log scale is not supported by Makie.")
        _warn("Lower frequency bound truncated to $(f[2]) Hz.")
        flim = (f[2], flim[2])
    end

    if ch_n > 64
        yts = 5
    elseif ch_n > 32
        yts = 2
    else
        yts = 1
    end

    # prepare plot
    if variant === :w
        plot_size = (900, 450)
        p = GLMakie.Figure(size=plot_size)
        ax = GLMakie.Axis3(p[1, 1],
                           xlabel=xlabel,
                           ylabel=ylabel,
                           zlabel=zlabel,
                           title=title,
                           xticks=LinearTicks(15),
                           # xminorticksvisible=true,
                           # xminorticks=IntervalsBetween(10),
                           # xscale=frq === :lin ? identity : log,
                           yticks=(1:yts:ch_n, clabels[1:yts:end]),
                           zoommode=:disable,
                           xtranslationlock=true,
                           ytranslationlock=true,
                           ztranslationlock=true,
                           aspect=(1, 1, 0.5),
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0),
                           zautolimitmargin=(0, 0))
        GLMakie.xlims!(ax, flim)
        ax.titlesize = 20
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        # plot powers
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            Makie.lines!(f,
                         ones(length(f)) .* idx,
                         ph[idx, :],
                         linewidth=2,
                         color=mono ? :black : cmap[idx],
                         colormap=pal,
                         colorrange=1:ch_n)
        end
    else
        f1 = vsearch(flim[1], f)
        f2 = vsearch(flim[2], f)
        plot_size = (900, 450)
        p = GLMakie.Figure(size=plot_size)
        ax = GLMakie.Axis3(p[1, 1],
                           xlabel=xlabel,
                           ylabel=ylabel,
                           zlabel=zlabel,
                           title=title,
                           xticks=LinearTicks(15),
                           # xminorticksvisible=true,
                           # xminorticks=IntervalsBetween(10),
                           # xscale=frq === :lin ? identity : log,
                           yticks=(1:yts:ch_n, clabels[1:yts:end]),
                           zoommode=:disable,
                           xtranslationlock=true,
                           ytranslationlock=true,
                           ztranslationlock=true,
                           aspect=(1, 1, 0.5),
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0),
                           zautolimitmargin=(0, 0))
        ax.titlesize = 20
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        # plot powers
        cmap = GLMakie.resample_cmap(pal, ch_n)
        Makie.surface!(f[f1:f2],
                       eachindex(clabels),
                       ph[:, f1:f2]',
                       colormap=pal)
    end

    return p

end

"""
    plot_phsd_topo(locs, f, ph; <keyword arguments>)

Plot topographical map of PHSDs (phase spectral density).

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `f::Vector{Float64}`: frequencies
- `ph::Array{Float64, 3}`: phases
- `flim::Tuple{Real, Real}=(f[1], f[end]): frequency limits
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: plot head shape

# Returns

- `p::GLMakie.Figure`
"""
function plot_phsd_topo(locs::DataFrame, f::Vector{Float64}, ph::Matrix{Float64}; flim::Tuple{Real, Real}=(f[1], f[end]), xlabel::String="", ylabel::String="", title::String="", frq::Symbol=:lin, cart::Bool=false, head::Bool=true)::GLMakie.Figure

    @assert size(ph, 2) == length(f) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, "flim")

    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(f[2]) Hz.")
        flim = (f[2], flim[2])
    end

    # plot parameters
    if size(ph, 1) <= 64
        plot_size = (1000, 1000)
        marker_size = (120, 80)
        xl = 1.2
        yl = 1.2
    elseif _in(size(ph, 1), (64, 100))
        plot_size = (1200, 1200)
        marker_size = (120, 80)
        xl = 1.5
        yl = 1.5
    else
        plot_size = (1500, 1500)
        marker_size = (85, 50)
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

    # prepare PHSD plots
    pp_vec = GLMakie.Figure[]
    pp_full_vec = GLMakie.Figure[]
    for idx in axes(ph, 1)
        pp = GLMakie.Figure(size=marker_size,
                            figure_padding=0)
        ax = GLMakie.Axis(pp[1, 1],
                          xlabel="",
                          ylabel="",
                          title=locs[idx, :label],
                          xscale=frq === :lin ? identity : log,
                          xautolimitmargin=(0, 0),
                          yautolimitmargin=(0.1, 0.1))
        hidedecorations!(ax)
        GLMakie.xlims!(ax, flim)
        ax.titlesize = 8
        # plot phases
        GLMakie.lines!(ax,
                       f,
                       ph[idx, :],
                       linewidth=1,
                       color=:black)
        push!(pp_vec, pp)
        pp_full = plot_phsd(f,
                            ph[idx, :],
                            xlabel=xlabel,
                            ylabel=ylabel,
                            title=locs[idx, :label] * ": " * title,
                            flim=flim,
                            frq=frq)
        push!(pp_full_vec, pp_full)
    end

    # prepare plot
    plot_size = (mplot_size, mplot_size)
    p = GLMakie.Figure(size=plot_size,
                       figure_padding=0)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel="",
                      ylabel="",
                      title=title,
                      aspect=1,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0.1, 0.1))
    GLMakie.xlims!(ax, (-xl, xl))
    GLMakie.ylims!(ax, (-yl, yl))
    hidespines!(ax)
    hidedecorations!(ax)
    ax.titlesize = 20

    if head
        # nose
        GLMakie.lines!(ax, [-0.1, 0], [0.995, 1.1], linewidth=3, color=:black)
        GLMakie.lines!(ax, [0, 0.1], [1.1, 0.995], linewidth=3, color=:black)

        # ears
        # left
        GLMakie.lines!(ax, [-0.995, -1.03], [0.1, 0.15], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.03, -1.06], [0.15, 0.16], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.06, -1.1], [0.16, 0.14], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.1, -1.12], [0.14, 0.05], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.12, -1.10], [0.05, -0.1], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.10, -1.13], [-0.1, -0.3], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.13, -1.09], [-0.3, -0.37], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.09, -1.02], [-0.37, -0.39], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.02, -0.98], [-0.39, -0.33], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-0.98, -0.975], [-0.33, -0.22], linewidth=3, color=:black)
        # right
        GLMakie.lines!(ax, [0.995, 1.03], [0.1, 0.15], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.03, 1.06], [0.15, 0.16], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.06, 1.1], [0.16, 0.14], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.1, 1.12], [0.14, 0.05], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.12, 1.10], [0.05, -0.1], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.10, 1.13], [-0.1, -0.3], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.13, 1.09], [-0.3, -0.37], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.09, 1.02], [-0.37, -0.39], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.02, 0.98], [-0.39, -0.33], linewidth=3, color=:black)
        GLMakie.lines!(ax, [0.98, 0.975], [-0.33, -0.22], linewidth=3, color=:black)

        # head
        GLMakie.arc!(ax,(0, 0), 1, 0, 2pi, linewidth=3, color=:black)
    end

    for idx in axes(ph, 1)
        io = IOBuffer()
        show(io, MIME"image/png"(), pp_vec[idx])
        pp = FileIO.load(io)
        GLMakie.scatter!(loc_x[idx],
                         loc_y[idx],
                         marker=pp,
                         markersize=marker_size,
                         markerspace=:pixel)
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
                    if ax_x >= loc_x_range[idx][1] && ax_x <= loc_x_range[idx][2] && ax_y >= loc_y_range[idx][1] && ax_y <= loc_y_range[idx][2]
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
    plot_phsd(obj; <keyword arguments>)

Plot PHSD (phase spectral density).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is Phase [rad]
- `zlabel::String="default"`: z-axis label for 3-d plots, default is Phase [rad]
- `title::String="default"`: plot title, default is PHSD [frequency limit: 0-128 Hz] [epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or gray palette
- `type::Symbol=:normal`: plot type:
    - `:normal` single channel or butterfly for multichannel
    - `:w3d`: 3-d waterfall
    - `:s3d`: 3-d surface
    - `:topo`: topographical
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: plot head shape
- `leg::Bool=true`: if true, add legend with channel labels
- `avg::Bool=false`: if true, plot averaged PSD
- `ci95::Bool=false`: if true, plot mean and ±95% CI of averaged PSDs

# Returns

- `p::GLMakie.Figure`
"""
function plot_phsd(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, ch::Union{String, Vector{String}, Regex}="all", flim::Tuple{Real, Real}=(0, sr(obj) / 2), frq::Symbol=:lin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, cart::Bool=false, head::Bool=true, leg::Bool=true, avg::Bool=false, ci95::Bool=false)::GLMakie.Figure

    _check_var(type, [:normal, :w3d, :s3d, :topo], "type")
    _check_var(frq, [:lin, :log], "frq")

    ch = get_channel(obj, ch=ch)
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

    # frequency limits
    fs = sr(obj)
    _check_tuple(flim, "flim", (0, sr(obj) / 2))

    # calculate PHSD
    sp, sf = phsd(signal, fs=fs)
    if ep != 0
        title == "default" && (title = "PHSD\n[epoch: $ep]")
    else
        title == "default" && (title = "PHSD\n[time window: $t_s1:$t_s2]")
    end

    if type === :normal
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Phase [rad]")
        if length(ch) == 1
            p = plot_phsd(sf,
                          sp,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          flim=flim,
                          frq=frq)
        else
            p = plot_phsd(sf,
                          sp,
                          xlabel=xlabel,
                          ylabel="",
                          clabels=clabels,
                          title=title,
                          flim=flim,
                          frq=frq,
                          avg=avg,
                          ci95=ci95,
                          leg=leg,
                          mono=mono)
        end
    elseif type === :w3d || type === :s3d
        ch_t = obj.header.recording[:channel_type]
        @assert ndims(sp) >= 2 "For type=:$type plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "")
        zlabel == "default" && (zlabel = "Phase [rad]")
        p = plot_phsd_3d(sf,
                         sp,
                         clabels=clabels,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         zlabel=zlabel,
                         title=title,
                         flim=flim,
                         frq=frq,
                         mono=mono,
                         variant=type === :w3d ? :w : :s)
    elseif type === :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Phase [rad]")
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        @assert length(unique(obj.header.recording[:channel_type][ch])) == 1 "For multi-channel topo plot all channels must be of the same type."
        _has_locs(obj)
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        ndims(sp) == 1 && (sp = reshape(sp, 1, length(sp)))
        p = plot_phsd_topo(locs,
                           sf,
                           sp,
                           xlabel=xlabel,
                           ylabel=ylabel,
                           title=title,
                           flim=flim,
                           frq=frq,
                           cart=cart,
                           head=head)
    end

    return p

end
