export mplot_psd
export mplot_psd_avg
export mplot_psd_butterfly
export mplot_psd_3d
export mplot_psd_topo

"""
    mplot_psd(sf, sp; <keyword arguments>)

Plot PSD (power spectrum density).

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Vector{Float64}`: powers
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_psd(sf::Vector{Float64}, sp::Vector{Float64}; frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, frq::Symbol=:lin, kwargs...)::GLMakie.Figure

    @assert length(sp) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    if frq === :log && frq_lim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        frq_lim = (sf[2], frq_lim[2])
    end

    # prepare plot
    plot_size = (1200, 600)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(15),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      xscale=frq === :lin ? identity : log10,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.xlims!(ax, frq_lim)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot powers
    Makie.lines!(sf,
                 sp,
                 color=:black)

    return p

end

"""
    mplot_psd(sf, sp; <keyword arguments>)

Plot multi-channel PSD (power spectrum density).

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Matrix{Float64}`: powers
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_psd(sf::Vector{Float64}, sp::Matrix{Float64}; clabels::Vector{String}=[""], frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, frq::Symbol=:lin, kwargs...)::GLMakie.Figure

    ch_n = size(sp, 1)
    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    # reverse so 1st channel is on top
    sp = @views reverse(sp[:, eachindex(sf)], dims = 1)
    # also, reverse colors if palette is not mono
    if mono
        pal = :grays
        channel_color = repeat([:black], ch_n)
    else
        pal = :darktest
        channel_color = ch_n:-1:1
    end

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    # normalize and shift so all channels are visible
    # each channel is between -1.0 and +1.0
    for idx in 1:ch_n
        # scale by 0.5 so maxima do not overlap
        sp[idx, :] = @views normalize_minmax(sp[idx, :]) .* 0.5 .+ (idx - 1)
    end

    if frq === :log && frq_lim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        frq_lim = (sf[2], frq_lim[2])
    end

    # prepare plot
    plot_size = 100 + 40 * ch_n <= 800 ? (1200, 600) : (1200, 100 + 40 * ch_n)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=100 + 40 * ch_n < 800 ? xlabel : "",
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(15),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      yticks=((ch_n - 1):-1:0, clabels),
                      yticksvisible=false,
                      xscale=frq === :lin ? identity : log10,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.xlims!(ax, frq_lim)
    GLMakie.ylims!(ax, -0.5, ch_n)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = ch_n <= 64 ? 12 : 10;

    # plot channels
    cmap = reverse(GLMakie.resample_cmap(pal, ch_n))
    for idx in 1:ch_n
        Makie.lines!(sf,
                     sp[idx, :],
                     linewidth=1,
                     color=mono ? :black : cmap[idx],
                     colormap=pal,
                     colorrange=1:ch_n)
    end

    return p

end

"""
    mplot_psd_avg(sf, sp; <keyword arguments>)

Plot PSD mean and ±95% CI of averaged channels.

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Matrix{Float64}`: powers
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_psd_avg(sf::Vector{Float64}, sp::Matrix{Float64}; frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, frq::Symbol=:lin, kwargs...)::GLMakie.Figure

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = NeuroAnalyzer.msci95(sp)

    if frq === :log && frq_lim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        frq_lim = (sf[2], frq_lim[2])
    end

    # prepare plot
    plot_size = (1200, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(15),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      xscale=frq === :lin ? identity : log10,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.xlims!(ax, frq_lim)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot upper 95% CI
    Makie.band!(sf,
                s_u,
                s_l,
                alpha=0.25,
                color=:grey,
                strokewidth=0.5)
    # plot mean
    Makie.lines!(sf,
                 s_m,
                 color=:black,
                 linewidth=1)

    return p

end

"""
    mplot_psd_butterfly(sf, sp; <keyword arguments>)

Butterfly PSD plot.

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `avg::Bool=false`: plot average channels
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_psd_butterfly(sf::Vector{Float64}, sp::Matrix{Float64}; clabels::Vector{String}=[""], frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, frq::Symbol=:lin, avg::Bool=false, kwargs...)::GLMakie.Figure

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    ch_n = size(sp, 1)
    pal = mono ? :grays : :darktest

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    if frq === :log && frq_lim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        frq_lim = (sf[2], frq_lim[2])
    end

    # prepare plot
    plot_size = (1200, 600)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(15),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      xscale=frq === :lin ? identity : log10,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.xlims!(ax, frq_lim)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    cmap = GLMakie.resample_cmap(pal, ch_n)
    for idx in 1:ch_n
        Makie.lines!(sf,
                     sp[idx, :],
                     color=cmap[idx],
                     colormap=pal,
                     colorrange=1:ch_n,
                     linewidth=0.5,
                     label=clabels[idx])
    end
    ch_n < 40 && axislegend(position = :rb)

    # plot averaged channels
    if avg
        s = mean(sp, dims=1)[:]
        Makie.lines!(sf,
                     s,
                     linewidth=2,
                     color=:black)
    end

    return p

end

"""
    mplot_psd_w3d(sf, sp; <keyword arguments>)

Plot 3-d waterfall PSD plot.

# Arguments

- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `db::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `zlabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `variant::Symbol`: waterfall (`:w`) or surface (`:s`)
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_psd_3d(sf::Vector{Float64}, sp::Matrix{Float64}; clabels::Vector{String}=[""], db::Bool=true, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", zlabel::String="", title::String="", mono::Bool=false, frq::Symbol=:lin, variant::Symbol, kwargs...)::GLMakie.Figure

    _check_var(variant, [:w, :s], "variant")
    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    ch_n = size(sp, 1)

    pal = mono ? :grays : :darktest

    # channel labels
    clabels == [""] && (clabels = repeat([""], ch_n))

    if frq === :log && frq_lim[1] == 0
        _warn("Currently log10 scale is not supported by Makie.")
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        frq_lim = (sf[2], frq_lim[2])
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
        plot_size = (1200, 600)
        p = GLMakie.Figure(size=plot_size)
        ax = GLMakie.Axis3(p[1, 1],
                           xlabel=xlabel,
                           ylabel=ylabel,
                           zlabel=zlabel,
                           title=title,
                           xticks=LinearTicks(15),
                           # xminorticksvisible=true,
                           # xminorticks=IntervalsBetween(10),
                           # xscale=frq === :lin ? identity : log10,
                           yticks=(1:yts:ch_n, clabels[1:yts:end]),
                           zoommode=:disable,
                           xtranslationlock=true,
                           ytranslationlock=true,
                           ztranslationlock=true,
                           aspect=(1, 1, 0.5),
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0),
                           zautolimitmargin=(0, 0);
                           kwargs...)
        GLMakie.xlims!(ax, frq_lim)
        ax.titlesize = 20
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        # plot powers
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            Makie.lines!(sf,
                         ones(length(sf)) .* idx,
                         sp[idx, :],
                         linewidth=1,
                         color=mono ? :black : cmap[idx],
                         colormap=pal,
                         colorrange=1:ch_n)
        end
    else
        f1 = vsearch(frq_lim[1], sf)
        f2 = vsearch(frq_lim[2], sf)
        plot_size = (1200, 600)
        p = GLMakie.Figure(size=plot_size)
        ax = GLMakie.Axis3(p[1, 1],
                           xlabel=xlabel,
                           ylabel=ylabel,
                           zlabel=zlabel,
                           title=title,
                           xticks=LinearTicks(15),
                           # xminorticksvisible=true,
                           # xminorticks=IntervalsBetween(10),
                           # xscale=frq === :lin ? identity : log10,
                           yticks=(1:yts:ch_n, clabels[1:yts:end]),
                           zoommode=:disable,
                           xtranslationlock=true,
                           ytranslationlock=true,
                           ztranslationlock=true,
                           aspect=(1, 1, 0.5),
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0),
                           zautolimitmargin=(0, 0);
                           kwargs...)
        ax.titlesize = 20
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        # plot powers
        cmap = GLMakie.resample_cmap(pal, ch_n)
        Makie.surface!(sf[f1:f2],
                       eachindex(clabels),
                       sp[:, f1:f2]',
                       colormap=pal)
    end

    return p

end

"""
    mplot_psd_topo(locs, sf, sp; <keyword arguments>)

Plot topographical map of PSDs.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: unused, for compatibility only
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: plot head shape
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_psd_topo(locs::DataFrame, sf::Vector{Float64}, sp::Matrix{Float64}; frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), title::String="", mono::Bool=true, frq::Symbol=:lin, cart::Bool=false, head::Bool=true, kwargs...)::GLMakie.Figure

    @assert size(sp, 2) == length(sf) "Length of powers vector must equal length of frequencies vector."
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    if frq === :log && frq_lim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz.")
        frq_lim = (sf[2], frq_lim[2])
    end

    # plot parameters
    if size(sp, 1) <= 64
        mplot_size = 1000
        marker_size = (120, 80)
        xl = 1.2
        yl = 1.2
    elseif _in(size(sp, 1), (64, 100))
        mplot_size = 1200
        marker_size = (120, 80)
        xl = 1.5
        yl = 1.5
    else
        mplot_size = 1500
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

    # prepare PSD plots
    pp_vec = GLMakie.Figure[]
    for idx in axes(sp, 1)
        pp = GLMakie.Figure(size=marker_size,
                            figure_padding=0)
        ax = GLMakie.Axis(pp[1, 1],
                          xlabel="",
                          ylabel="",
                          title=locs[idx, :label],
                          xscale=frq === :lin ? identity : log10,
                          xautolimitmargin=(0, 0),
                          yautolimitmargin=(0, 0);
                          kwargs...)
        hidedecorations!(ax)
        GLMakie.xlims!(ax, frq_lim)
        ax.titlesize = 8
        # plot powers
        GLMakie.lines!(sf,
                       sp[idx, :],
                       linewidth=1,
                       color=:black)
        push!(pp_vec, pp)
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
                      yautolimitmargin=(0, 0);
                      kwargs...)
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

    for idx in axes(sp, 1)
        io = IOBuffer()
        show(io, MIME"image/png"(), pp_vec[idx])
        pp = FileIO.load(io)
        GLMakie.scatter!(loc_x[idx],
                         loc_y[idx],
                         marker=pp,
                         markersize=marker_size,
                         markerspace=:pixel)
    end

    return p

end

"""
    mplot_psd(obj; <keyword arguments>)

Plot power spectrum density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `db::Bool=true`: normalize powers to dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-taper periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher`
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is `Power [dB units^2/Hz] or Power [units^2/Hz]`
- `zlabel::String="default"`: z-axis label for 3-d plots, default is `Power [dB units^2/Hz] or Power [units^2/Hz]`
- `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or gray palette
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:butterfly`
    - `:mean`
    - `:w3d`: 3-d waterfall
    - `:s3d`: 3-d surface
    - `:topo`: topographical
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: plot head shape
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_psd(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, ch::Union{String, Vector{String}, Regex}, db::Bool=true, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, ref::Symbol=:abs, frq::Symbol=:lin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, cart::Bool=false, head::Bool=true, kwargs...)::GLMakie.Figure

    _check_var(type, [:normal, :butterfly, :mean, :w3d, :s3d, :topo], "type")
    _check_var(method, [:welch, :fft, :stft, :mt, :mw, :gh], "method")
    _check_var(ref, [:abs, :total, :delta, :theta, :alpha, :alpha_lower, :alpha_higher, :beta, :beta_lower, :beta_higher, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "ref")
    _check_var(frq, [:lin, :log], "frq")

    ch = get_channel(obj, ch=ch)
    _check_tuple(frq_lim, "frq_lim", (0, sr(obj) / 2))

    @assert seg[1] != seg[2] "Signal is too short for analysis."

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    if ep != 0
        _check_epochs(obj, ep)
        if nepochs(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    # get time vector
    if seg[2] <= epoch_len(obj)
        signal = obj.data[ch, seg[1]:seg[2], 1]
    else
        signal = epoch(obj, ep_n=1).data[ch, seg[1]:seg[2], 1]
    end
    # t = _get_t(seg[1], seg[2], sr(obj))
    t = obj.time_pts[seg[1]:seg[2]]
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    clabels = labels(obj)

    # set units
    units = _ch_units(obj, labels(obj)[ch[1]])

    ref !== :abs && (frq_lim = band_frq(obj, band=ref))

    # get frequency range
    fs = sr(obj)
    _check_tuple(frq_lim, "frq_lim", (0, sr(obj) / 2))

    # # get time vector
    # _, t_s1, _, t_s2 = _convert_t(obj.epoch_time[1], obj.epoch_time[end])

    if ref === :abs
        if method === :welch
            sp, sf = psd(signal, fs=fs, db=db, method=:welch, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (Welch's periodogram)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd(signal, fs=fs, db=db, method=:fft, w=w)
            title == "default" && (title = "Absolute PSD (fast Fourier transform)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (short-time Fourier transform)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (multi-taper)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd(signal, fs=fs, db=db, method=:mw, ncyc=ncyc, w=w)
            title == "default" && (title = "Absolute PSD (Morlet wavelet convolution)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "Absolute PSD (Gaussian and Hilbert transform)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        end
    elseif ref === :total
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, db=db, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (Welch's periodogram) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:fft, w=w)
            title == "default" && (title = "PSD (fast Fourier transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (short-time Fourier transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (multi-taper) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mw, ncyc=ncyc, w=w)
            title == "default" && (title = "PSD (Morlet wavelet convolution) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "PSD (Gaussian and Hilbert transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        end
    else
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, db=db, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (Welch's periodogram) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:fft, frq_lim=frq_lim, w=w)
            title == "default" && (title = "PSD (fast Fourier transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:stft, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (short-time Fourier transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mt, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (multi-taper) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mw, frq_lim=frq_lim, ncyc=ncyc, w=w)
            title == "default" && (title = "PSD (Morlet wavelet convolution) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "PSD (Gaussian and Hilbert transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        end
    end

    # set labels
    if type !== :w3d && type !== :s3d && type !== :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if ref !== :abs
            ylabel == "default" && (ylabel = "Power ratio")
        end
        ylabel == "default" && (ylabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
    end

    if type === :normal
        ch_t = obj.header.recording[:channel_type]
        if length(ch) > 1
            ch_t_uni = unique(ch_t[ch])
            @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels must be of the same type."
        end
        if size(sp, 1) == 1
            p = mplot_psd(sf,
                         sp[1, :],
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         frq_lim=frq_lim,
                         frq=frq,
                         mono=mono;
                         kwargs...)
        else
            p = mplot_psd(sf,
                         sp,
                         xlabel=xlabel,
                         ylabel="",
                         clabels=clabels[ch],
                         title=title,
                         frq_lim=frq_lim,
                         frq=frq,
                         mono=mono;
                         kwargs...)
        end
    elseif type === :butterfly
        pl = mplot_locs(obj, ch=labels(obj)[ch], selected=ch=labels(obj)[ch], ps=:s)
        ch_t = obj.header.recording[:channel_type]
        ch_t_uni = unique(ch_t[ch])
        @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels must be of the same type."
        @assert ndims(sp) >= 2 "For type=:butterfly plot the signal must contain ≥ 2 channels."
        title = replace(title, "channel" => "channels")
        p = mplot_psd_butterfly(sf,
                               sp,
                               clabels=clabels[ch],
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               frq_lim=frq_lim,
                               frq=frq,
                               mono=mono;
                               kwargs...)
        GLMakie.scatter!(p[1, 1],
                         pl)
    elseif type === :mean
        ch_t = obj.header.recording[:channel_type]
        ch_t_uni = unique(ch_t[ch])
        @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels must be of the same type."
        @assert ndims(sp) >= 2 "For type=:mean plot the signal must contain ≥ 2 channels."
        title = replace(title, "PSD" => "PSD [mean ± 95%CI]")
        title = replace(title, "channel" => "averaged channels")
        p = mplot_psd_avg(sf,
                         sp,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         frq_lim=frq_lim,
                         frq=frq,
                         mono=mono;
                         kwargs...)
    elseif type === :w3d
        ch_t = obj.header.recording[:channel_type]
        ch_t_uni = unique(ch_t[ch])
        @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels must be of the same type."
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "")
        zlabel == "default" && (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = mplot_psd_3d(sf,
                        sp,
                        clabels=clabels[ch],
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        frq_lim=frq_lim,
                        frq=frq,
                        mono=mono,
                        variant=:w;
                        kwargs...)
    elseif type === :s3d
        ch_t = obj.header.recording[:channel_type]
        ch_t_uni = unique(ch_t[ch])
        @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels must be of the same type."
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "")
        zlabel == "default" && (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = mplot_psd_3d(sf,
                        sp,
                        clabels=clabels[ch],
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        frq_lim=frq_lim,
                        frq=frq,
                        mono=mono,
                        variant=:s;
                        kwargs...)
    elseif type === :topo
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        @assert length(unique(obj.header.recording[:channel_type][ch])) == 1 "For multi-channel PSD plots all channels must be of the same type."
        _has_locs(obj)
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        ndims(sp) == 1 && (sp = reshape(sp, 1, length(sp)))
        xlabel == "default" && (xlabel = "")
        ylabel == "default" && (ylabel = "")
        title = replace(title, "channel" => "channels")
        p = mplot_psd_topo(locs,
                          sf,
                          sp,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          frq_lim=frq_lim,
                          frq=frq,
                          cart=cart,
                          head=head;
                          kwargs...)
    end

    return p

end

"""
    mplot_psd(obj; <keyword arguments>)

Plot power spectrum density of embedded or external component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `db::Bool=true`: normalize powers to dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-taper periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher`
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is `Power [dB units^2/Hz] or Power [units^2/Hz]`
- `zlabel::String="default"`: z-axis label for 3-d plots, default is `Power [dB units^2/Hz] or Power [units^2/Hz]`
- `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or gray palette
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:butterfly`
    - `:mean`
    - `:w3d`: 3-d waterfall
    - `:s3d`: 3-d surface
    - `:topo`: topographical
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: plot head shape
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_psd(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, db::Bool=true, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, ref::Symbol=:abs, frq::Symbol=:lin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, cart::Bool=false, head::Bool=true, kwargs...)::GLMakie.Figure

    _check_var(type, [:normal, :butterfly, :mean, :w3d, :s3d, :topo], "type")
    _check_var(method, [:welch, :fft, :stft, :mt, :mw, :gh], "method")
    _check_var(ref, [:abs, :total, :delta, :theta, :alpha, :alpha_lower, :alpha_higher, :beta, :beta_lower, :beta_higher, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "ref")
    _check_var(frq, [:lin, :log], "frq")

    units = "A.U."

    @assert seg[1] != seg[2] "Signal is too short for analysis."

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    if ep != 0
        _check_epochs(obj, ep)
        if nepochs(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    # select component channel, default is all channels
    c isa Symbol && (c = _get_component(obj, c))
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    clabels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (clabels = [clabels])

    ref !== :abs && (frq_lim = band_frq(obj, band=ref))

    # get frequency range
    fs = sr(obj)
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    # get time vector
    if seg[2] <= epoch_len(obj)
        signal = c[c_idx, seg[1]:seg[2], 1]
    else
        signal = reshape(c, size(c, 1), :, 1)[c_idx, seg[1]:seg[2], 1]
    end
    # t = _get_t(seg[1], seg[2], sr(obj))
    t = obj.time_pts[seg[1]:seg[2]]
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    if ref === :abs
        if method === :welch
            sp, sf = psd(signal, fs=fs, db=db, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (Welch's periodogram)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd(signal, fs=fs, db=db, method=:fft, w=w)
            title == "default" && (title = "Absolute PSD (fast Fourier transform)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (short-time Fourier transform)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (multi-taper)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd(signal, fs=fs, db=db, method=:mw, ncyc=ncyc, w=w)
            title == "default" && (title = "Absolute PSD (Morlet wavelet convolution)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "Absolute PSD (Gaussian and Hilbert transform)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end
    elseif ref === :total
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, db=db, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (Welch's periodogram) relative to total power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:fft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (fast Fourier transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (short-time Fourier transform) relative to total power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (multi-taper) relative to total power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mw, ncyc=ncyc, w=w)
            title == "default" && (title = "PSD (Morlet wavelet convolution) relative to total power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "PSD (Gaussian and Hilbert transform) relative to total power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end
    else
        if method === :welch
            sp, sf = psd_rel(signal, fs=fs, db=db, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) relative to $ref power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :fft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:fft, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (fast Fourier transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :stft
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:stft, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "PSD (short-time Fourier transform) relative to $(replace(string(ref), "_"=>" ")) power\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mt, frq_lim=frq_lim, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Absolute PSD (multi-taper) relative to $(replace(string(ref), "_"=>" ")) power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:mw, frq_lim=frq_lim, ncyc=ncyc, w=w)
            title == "default" && (title = "PSD (Morlet wavelet convolution) relative to $(replace(string(ref), "_"=>" ")) power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd_rel(signal, fs=fs, db=db, method=:gh, gw=gw, w=w)
            title == "default" && (title = "PSD (Gaussian and Hilbert transform) relative to $(replace(string(ref), "_"=>" ")) power\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end
    end

    # set labels
    if type !== :w3d && type !== :s3d && type !== :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if ref !== :abs
            ylabel == "default" && (ylabel = "Power ratio")
        end
        ylabel == "default" && (ylabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
    end

    if type === :normal
        @assert ndims(sp) == 1 "For type=:normal the signal must contain 1 c_idx."
        p = mplot_psd(sf,
                     sp,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     title=title,
                     frq_lim=frq_lim,
                     frq=frq,
                     mono=mono;
                     kwargs...)
    elseif type === :butterfly
        @assert ndims(sp) >= 2 "For type=:butterfly plot the signal must contain ≥ 2 c_idxs."
        title = replace(title, "component" => "components")
        p = mplot_psd_butterfly(sf,
                               sp,
                               clabels=clabels,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               frq_lim=frq_lim,
                               frq=frq,
                               mono=mono;
                               kwargs...)
    elseif type === :mean
        @assert ndims(sp) >= 2 "For type=:mean plot the signal must contain ≥ 2 c_idxs."
        title = replace(title, "PSD" => "PSD [mean ± 95%CI]")
        title = replace(title, "component" => "averaged components")
        p = mplot_psd_avg(sf,
                         sp,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         frq_lim=frq_lim,
                         frq=frq,
                         mono=mono;
                         kwargs...)
    elseif type === :w3d
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "")
        zlabel == "default" && (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = mplot_psd_3d(sf,
                        sp,
                        clabels=clabels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        frq_lim=frq_lim,
                        frq=frq,
                        mono=mono,
                        variant=:w;
                        kwargs...)
    elseif type === :s3d
        @assert ndims(sp) >= 2 "For type=:w3d plot the signal must contain ≥ 2 channels."
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "")
        zlabel == "default" && (zlabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = mplot_psd_3d(sf,
                        sp,
                        clabels=clabels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        frq_lim=frq_lim,
                        frq=frq,
                        mono=mono,
                        variant=:s;
                        kwargs...)
    elseif type === :topo
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        @assert length(unique(obj.header.recording[:channel_type][ch])) == 1 "For multi-channel PSD plots all channels must be of the same type."
        _has_locs(obj)
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        ndims(sp) == 1 && (sp = reshape(sp, 1, length(sp)))
        xlabel == "default" && (xlabel = "")
        ylabel == "default" && (ylabel = "")
        title = replace(title, "channel" => "channels")
        p = mplot_psd_topo(locs,
                          sf,
                          sp,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          frq_lim=frq_lim,
                          frq=frq,
                          cart=cart,
                          head=head;
                          kwargs...)
    end

    return p

end
