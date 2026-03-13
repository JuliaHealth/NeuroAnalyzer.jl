export plot_spectrogram
export plot_spectrogram_topo

"""
    plot_spectrogram(st, sf, sp; <keyword arguments>)

Plot single-channel spectrogram.

# Arguments

- `st::Vector{Float64}`: time
- `sf::Vector{<:Real}`: frequencies
- `sp::Matrix{Float64}`: powers
- `db::Bool=true`: whether powers are normalized to dB
- `frq::Symbol=:lin`: frequency scaling - `:lin` or `:log`
- `flim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `units::String=""`
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `cb::Bool=true`: plot color bar
- `cb_title::String=""`: color bar label
- `threshold::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
    - `:in`: draw region is values are in the threshold values, including threshold boundaries
    - `:bin`: draw region is values are between the threshold values, excluding threshold boundaries

# Returns

- `GLMakie.Figure`
"""
function plot_spectrogram(
        st::Vector{Float64},
        sf::Vector{<:Real},
        sp::Matrix{Float64};
        db::Bool = true,
        frq::Symbol = :lin,
        flim::Tuple{Real, Real} = (sf[1], sf[end]),
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        mono::Bool = false,
        units::String = "",
        smooth::Bool = false,
        n::Int64 = 3,
        cb::Bool = true,
        cb_title::String = "",
        threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
        threshold_type::Symbol = :neq,
    )::GLMakie.Figure

    @assert size(sp, 2) == length(st) "Size of powers ($(size(sp, 2))) and time vector ($(length(st))) do not match."
    @assert size(sp, 1) == length(sf) "Size of powers ($(size(sp, 1))) and frequencies vector ($(length(sf))) do not match."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    pal = mono ? :grays : :darktest

    if smooth
        sp = imfilter(sp, Kernel.gaussian(n))
    end

    sp = sp'

    if !isnothing(threshold)
        sp_threshold = deepcopy(sp)
        _, bm = seg_extract(sp, threshold = threshold, threshold_type = threshold_type)
        sp_threshold[.!bm] .= NaN
    end

    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        flim = (sf[2], flim[2])
    end

    # prepare plot
    GLMakie.activate!(title = "plot_spectrogram()")
    plot_size = (1200, 800)
    fig = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(10),
        yscale = frq === :lin ? identity : log,
        xgridvisible = false,
        ygridvisible = false,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.xlims!(ax, (st[1], st[end]))
    GLMakie.ylims!(ax, flim)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # draw spectrogram
    if !isnothing(threshold)
        hm = GLMakie.heatmap!(ax, st, sf, sp_threshold, colorrange = extrema(sp[.!isnan.(sp)]), colormap = pal)
    else
        hm = GLMakie.heatmap!(ax, st, sf, sp, colormap = pal)
    end

    # draw colorbar
    if cb
        Colorbar(fig[1, 2], hm, label = cb_title, labelsize = 16)
    end

    return fig

end

"""
    plot_spectrogram(sf, sp; <keyword arguments>)

Plot multiple-channel spectrogram.

# Arguments

- `sf::Vector{<:Real}`: frequencies
- `sp::Matrix{Float64}`: powers
- `clabels::Vector{String}=string.(1:size(sp, 1))`: channel labels
- `db::Bool=true`: whether powers are normalized to dB
- `frq::Symbol=:lin`: frequency scaling - `:lin` or `:log`
- `flim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `units::String=""`
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `cb::Bool=true`: plot color bar
- `cb_title::String=""`: color bar label
- `threshold::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
    - `:in`: draw region is values are in the threshold values, including threshold boundaries
    - `:bin`: draw region is values are between the threshold values, excluding threshold boundaries

# Returns

- `GLMakie.Figure`
"""
function plot_spectrogram(
        sf::Vector{<:Real},
        sp::Matrix{Float64};
        clabels::Vector{String} = string.(1:size(sp, 1)),
        db::Bool = true,
        frq::Symbol = :lin,
        flim::Tuple{Real, Real} = (sf[1], sf[end]),
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        mono::Bool = false,
        units::String = "",
        smooth::Bool = false,
        n::Int64 = 3,
        cb::Bool = true,
        cb_title::String = "",
        threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
        threshold_type::Symbol = :neq,
    )::GLMakie.Figure

    @assert size(sp, 1) == length(clabels) "Size of powers ($(size(sp, 1))) and channels vector ($(length(clabels))) do not match."
    @assert size(sp, 2) == length(sf) "Size of powers ($(size(sp, 2))) and frequencies vector ($(length(sf))) do not match."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    pal = mono ? :grays : :darktest

    if smooth
        sp = imfilter(sp, Kernel.gaussian(n))
    end

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        flim = (sf[2], flim[2])
    end

    ch = collect(eachindex(clabels)) .- 0.5
    ch_n = length(ch)
    reverse!(sp, dims = 1)

    # prepare plot
    GLMakie.activate!(title = "plot_spectrogram()")
    plot_size = (1200, 800)
    fig = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(15),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        yticks = (0.5:1:ch_n, reverse(clabels)),
        yticksvisible = false,
        xscale = frq === :lin ? identity : log,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.xlims!(ax, flim)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    hm = GLMakie.heatmap!(ax, sf, ch, sp', colormap = pal)

    # draw thresholded region
    if !isnothing(threshold)
        _, bm = seg_extract(sp, threshold = threshold, threshold_type = threshold_type)
        reg = ones(size(sp)) .* minimum(sp)
        reg[bm] .= maximum(sp)
        GLMakie.contour!(ax, sf, ch, reg', levels = 1, color = :black, linewidth = 2)
    end

    # draw colorbar
    if cb
        Colorbar(fig[1, 2], hm, label = cb_title, labelsize = 16)
    end

    return fig

end


"""
    plot_spectrogram_topo(locs, st, sf, sp; <keyword arguments>)

Plot topographical map of spectrograms.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `st::Vector{Float64}`: time
- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `db::Bool=true`: whether powers are normalized to dB
- `flim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `units::String=""`
- `cb::Bool=true`: plot color bar
- `cb_title::String=""`: color bar label
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `mono::Bool=false`: unused, for compatibility only
- `frq::Symbol=:lin`: frequency scaling - `:lin` or `:log`
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: plot head shape

# Returns

- `GLMakie.Figure`
"""
function plot_spectrogram_topo(
        locs::DataFrame,
        st::Vector{Float64},
        sf::Vector{Float64},
        sp::Array{Float64, 3};
        db::Bool = true,
        flim::Tuple{Real, Real} = (sf[1], sf[end]),
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        units::String = "",
        cb::Bool = true,
        cb_title::String = "",
        smooth::Bool = false,
        n::Int64 = 3,
        mono::Bool = true,
        frq::Symbol = :lin,
        cart::Bool = false,
        head::Bool = true,
    )::GLMakie.Figure

    @assert size(sp, 3) == DataFrames.nrow(locs) "Size of powers ($(size(sp, 3))) and number of locs ($(DataFrames.nrow(locs))) do not match."
    @assert size(sp, 2) == length(st) "Size of powers ($(size(sp, 2))) and time vector ($(length(st))) do not match."
    @assert size(sp, 1) == length(sf) "Size of powers ($(size(sp, 1))) and frequencies vector ($(length(sf))) do not match."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    pal = mono ? :grays : :darktest

    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        flim = (sf[2], flim[2])
    end

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

    if smooth
        sp = imfilter(sp, Kernel.gaussian(n))
    end

    # prepare spectrogram plots
    pp_vec = GLMakie.Figure[]
    pp_full_vec = GLMakie.Figure[]
    for idx in axes(sp, 3)
        pp = GLMakie.Figure(size = marker_size, figure_padding = 0)
        ax = GLMakie.Axis(
            pp[1, 1];
            xlabel = "",
            ylabel = "",
            aspect = nothing,
            title = locs[idx, :label],
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
        )
        hidedecorations!(ax)
        GLMakie.xlims!(ax, flim)
        ax.titlesize = 8
        # plot powers
        GLMakie.heatmap!(ax, sf, st, sp[:, :, idx]', colormap = pal)
        push!(pp_vec, pp)
        pp_full = plot_spectrogram(
            st,
            sf,
            sp[:, :, idx];
            db = db,
            frq = frq,
            flim = flim,
            xlabel = xlabel,
            ylabel = ylabel,
            title = locs[idx, :label] * ": " * title,
            mono = mono,
            units = units,
            smooth = smooth,
            n = n,
            cb = cb,
            cb_title = cb_title,
        )
        push!(pp_full_vec, pp_full)
    end

    # prepare plot
    GLMakie.activate!(title = "plot_spectrogram()")
    fig = GLMakie.Figure(
        size = plot_size,
        figure_padding = 0,
    )
    ax = GLMakie.Axis(
        fig[1, 1];
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
        GLMakie.lines!(ax, [-0.2, 0], [0.98, 1.08], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [0.2, 0], [0.98, 1.08], linewidth = 3, color = :black)

        # ears
        # left
        GLMakie.lines!(ax, [-0.995, -1.03], [0.1, 0.15], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.03, -1.06], [0.15, 0.16], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.06, -1.1], [0.16, 0.14], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.1, -1.12], [0.14, 0.05], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.12, -1.1], [0.05, -0.1], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.1, -1.13], [-0.1, -0.3], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.13, -1.09], [-0.3, -0.37], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.09, -1.02], [-0.37, -0.39], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.02, -0.98], [-0.39, -0.33], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-0.98, -0.975], [-0.33, -0.22], linewidth = 3, color = :black)
        # right
        GLMakie.lines!(ax, [0.995, 1.03], [0.1, 0.15], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.03, 1.06], [0.15, 0.16], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.06, 1.1], [0.16, 0.14], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.1, 1.12], [0.14, 0.05], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.12, 1.1], [0.05, -0.1], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.1, 1.13], [-0.1, -0.3], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.13, 1.09], [-0.3, -0.37], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.09, 1.02], [-0.37, -0.39], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.02, 0.98], [-0.39, -0.33], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [0.98, 0.975], [-0.33, -0.22], linewidth = 3, color = :black)

        # head
        GLMakie.arc!(ax, (0, 0), 1, 0, 2pi, linewidth = 3, color = :black)
    end

    for idx in axes(sp, 3)
        io = IOBuffer()
        show(io, MIME"image/png"(), pp_vec[idx])
        pp = FileIO.load(io)
        GLMakie.scatter!(loc_x[idx], loc_y[idx], marker = pp, markersize = marker_size, markerspace = :pixel)
    end

    loc_x_range = Tuple{Float64, Float64}[]
    loc_y_range = Tuple{Float64, Float64}[]
    for idx in eachindex(loc_x)
        push!(loc_x_range, (loc_x[idx] - 0.15, loc_x[idx] + 0.15))
        push!(loc_y_range, (loc_y[idx] - 0.1, loc_y[idx] + 0.1))
    end
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
                        display(GLMakie.Screen(), pp_full_vec[idx])
                        break
                    end
                end
            end
        end
    end

    return fig

end

"""
    plot_spectrogram(obj; <keyword arguments>)

Plots spectrogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}=datatype(obj)`: channel name or list of channel names
- `db::Bool=true`: normalize powers to dB; for CWT scaleogram: normalize to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal
- `method::Symbol=:stft`: spectrogram method:
    - `:stft`: short-time Fourier
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
    - `:hht`: Hilbert-Huang transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `gw::Real=10`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `frq::Symbol=:lin`: frequency scaling - `:lin` or `:log`
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: y-axis limits
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `mono::Bool=false`: use color or gray palette
- `markers::Bool`: draw markers if available
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `cb::Bool=true`: plot color bar
- `threshold::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
    - `:in`: draw region is values are in the threshold values, including threshold boundaries
    - `:bin`: draw region is values are between the threshold values, excluding threshold boundaries
- `type::Symbol=:normal`:
    - `:normal`
    - `:topo`
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: plot head shape

# Returns

- `GLMakie.Figure`
"""
function plot_spectrogram(
    obj::NeuroAnalyzer.NEURO;
    seg::Tuple{Real, Real} = (0, 10),
    ep::Int64 = 0,
    ch::Union{String, Vector{String}, Regex} = datatype(obj),
    db::Bool = true,
    method::Symbol = :stft,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    gw::Real = 10,
    wt::T = wavelet(Morlet(2π), β = 2),
    frq::Symbol = :lin,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    mono::Bool = false,
    markers::Bool = true,
    smooth::Bool = false,
    n::Int64 = 3,
    cb::Bool = true,
    threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    threshold_type::Symbol = :neq,
    type::Symbol = :normal,
    cart::Bool = false,
    head::Bool = true,
)::GLMakie.Figure where {T <: CWT}

    _check_var(type, [:normal, :topo], "type")
    _check_var(method, [:stft, :mt, :mw, :gh, :cwt, :hht], "method")
    @assert n > 0 "n must be ≥ 1."

    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    if method === :cwt
        if type === :normal
            @assert length(ch) == 1 "For :cwt method only one channel must be selected."
        end
    end
    if type === :topo
        @assert method !== :hht "For :hht method topographical map is not available."
    end
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
    _check_tuple(flim, (0, sr(obj) / 2), "flim")

    # calculate spectrogram
    if length(ch) == 1 || type === :topo
        if method === :stft
            sp, sf, st = NeuroAnalyzer.spectrogram(
                signal, fs = fs, db = false, method = :stft, wlen = wlen, woverlap = woverlap, w = w
            )
            if ep != 0
                title == "default" && (title = "Spectrogram (short-time Fourier)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (short-time Fourier)")
            end
        elseif method === :mt
            sp, sf, st = NeuroAnalyzer.spectrogram(
                signal, fs = fs, db = false, method = :mt, nt = nt, wlen = wlen, woverlap = woverlap, w = w
            )
            if ep != 0
                title == "default" && (title = "Spectrogram (multi-tapered)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (multi-tapered)")
            end
        elseif method === :mw
            _, sp, _, sf, st = NeuroAnalyzer.mwspectrogram(signal, fs = fs, ncyc = ncyc, db = false, w = w)
            if ep != 0
                title == "default" && (title = "Spectrogram (Morlet wavelet)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (Morlet wavelet)")
            end
        elseif method === :gh
            sp, _, sf, st = NeuroAnalyzer.ghtspectrogram(signal, fs = fs, db = false, gw = gw, w = w)
            if ep != 0
                title == "default" && (title = "Spectrogram (Gaussian-Hilbert)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (Gaussian-Hilbert)")
            end
        elseif method === :cwt
            _log_off()
            sp, sf, st = NeuroAnalyzer.cwtspectrogram(signal, fs = fs, wt = wt)
            _log_on()
            sf[1] > flim[1] && (flim = (sf[1], flim[2]))
            sf[end] < flim[2] && (flim = (flim[1], sf[end]))
            if ep != 0
                title == "default" && (title = "CWT Scaleogram\n[epoch: $ep]")
            else
                title == "default" && (title = "CWT Scaleogram")
            end
        elseif method === :hht
            imf = emd(signal, t)
            sp, _, sf, st = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs = fs, db = false)
            if ep != 0
                title == "default" && (title = "Spectrogram (Hilbert-Huang)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (Hilbert-Huang)")
            end
        end
    elseif length(ch) > 1 && type === :normal
        if method === :stft
            sp, sf = psd(signal, fs = fs, db = db, method = :stft, nt = nt, wlen = wlen, woverlap = woverlap, w = w)
            if ep != 0
                title == "default" && (title = "Spectrogram (short-time Fourier)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (short-time Fourier)")
            end
        elseif method === :mt
            sp, sf = psd(signal, fs = fs, db = db, method = :mt, nt = nt, wlen = wlen, woverlap = woverlap, w = w)
            if ep != 0
                title == "default" && (title = "Spectrogram (multi-tapered)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (multi-tapered)")
            end
        elseif method === :mw
            sp, sf = psd(signal, fs = fs, db = db, method = :mw, w = w, ncyc = ncyc)
            if ep != 0
                title == "default" && (title = "Spectrogram (Morlet wavelet)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (Morlet wavelet)")
            end
        elseif method === :gh
            sp, sf = psd(signal, fs = fs, db = db, method = :gh, w = w, gw = gw)
            if ep != 0
                title == "default" && (title = "Spectrogram (Gaussian-Hilbert)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (Gaussian-Hilbert)")
            end
        elseif method === :cwt
            _log_off()
            sp, sf = psd(signal, fs = fs, method = :cwt, wt = wt)
            _log_on()
            sf[1] > flim[1] && (flim = (sf[1], flim[2]))
            sf[end] < flim[2] && (flim = (flim[1], sf[end]))
            if ep != 0
                title == "default" && (title = "CWT Scaleogram\n[epoch: $ep]")
            else
                title == "default" && (title = "CWT Scaleogram")
            end
        elseif method === :hht
            imf = emd(signal[1, :], t)
            sp_tmp, _, sf, st = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs = fs, db = db)
            sp = zeros(size(signal, 1), length(sp_tmp))
            sp[1, :] = sp_tmp
            for idx in axes(signal, 1)[(begin + 1):end]
                imf = emd(signal[idx, :], t)
                sp[idx, :], _, _, _ = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs = fs, db = db)
            end
            if ep != 0
                title == "default" && (title = "Spectrogram (Hilbert-Huang)\n[epoch: $ep]")
            else
                title == "default" && (title = "Spectrogram (Hilbert-Huang)")
            end
        end
    end

    f1 = vsearch(flim[1], sf)
    f2 = vsearch(flim[2], sf)
    sf = sf[f1:f2]
    if length(ch) == 1 && type === :normal
        sp = sp[f1:f2, :]
        st .+= t[1]
        method !== :cwt && db && (sp = pow2db.(sp))
    elseif length(ch) >= 1 && type === :topo
        sp = sp[f1:f2, :, :]
        st .+= t[1]
        method !== :cwt && db && (sp = pow2db.(sp))
    else
        sp = sp[:, f1:f2]
    end

    cb_title = method === :cwt ? "Magnitude" : "Power"
    method !== :cwt && (cb_title *= db ? " [dB $units^2/Hz]" : " [$units^2/Hz]")

    if length(ch) == 1 && type === :normal
        xlabel == "default" && (xlabel = "Time [s]")
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        fig = plot_spectrogram(
            st,
            sf,
            sp;
            db = db,
            frq = frq,
            flim = flim,
            xlabel = xlabel,
            ylabel = ylabel,
            title = title,
            mono = mono,
            units = units,
            smooth = smooth,
            n = n,
            cb = cb,
            cb_title = cb_title,
            threshold = threshold,
            threshold_type = threshold_type,
        )
    elseif length(ch) > 1 && type === :normal
        ylabel == "default" && (ylabel = "")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        fig = plot_spectrogram(
            sf,
            sp;
            clabels = clabels,
            db = db,
            frq = frq,
            flim = flim,
            xlabel = xlabel,
            ylabel = ylabel,
            title = title,
            mono = mono,
            units = units,
            smooth = smooth,
            n = n,
            cb = cb,
            cb_title = cb_title,
            threshold = threshold,
            threshold_type = threshold_type,
        )
    elseif type === :topo
        xlabel == "default" && (xlabel = "Time [s]")
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        @assert length(unique(obj.header.recording[:channel_type][ch])) == 1 "For multi-channel topo plot all channels must be of the same type."
        _has_locs(obj)
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        ndims(sp) == 1 && (sp = reshape(sp, 1, length(sp)))
        fig = plot_spectrogram_topo(
            locs,
            st,
            sf,
            sp;
            frq = frq,
            flim = flim,
            xlabel = xlabel,
            ylabel = ylabel,
            title = title,
            mono = mono,
            units = units,
            cart = cart,
            smooth = smooth,
            n = n,
            cb = cb,
            cb_title = cb_title,
            head = head,
        )
    end

    # plot markers if available
    if length(ch) == 1 && markers && _has_markers(obj)
        markers_pos = obj.markers[!, :start]
        markers_id = obj.markers[!, :id]
        markers_desc = obj.markers[!, :value]
        for idx in eachindex(markers_pos)
            if _in(markers_pos[idx], (obj.time_pts[1], obj.time_pts[end]))
                GLMakie.vlines!(fig[1, 1], markers_pos[idx], linestyle = :dash, linewidth = 1, color = :black)
                if length(ch) > 1
                    GLMakie.textlabel!(
                        fig[1, 1],
                        (markers_pos[idx] + 0.07, ch_n > 20 ? 20.4 : ch_n + 0.4),
                        text = "$(markers_id[idx]) / $(markers_desc[idx])",
                        text_align = (:left, :center),
                        fontsize = 8,
                        text_rotation = pi / 2,
                    )
                else
                    GLMakie.textlabel!(
                        fig[1, 1],
                        (markers_pos[idx] + 0.07, 0.97 * minimum(obj.data[ch, :, :])),
                        text = "$(markers_id[idx]) / $(markers_desc[idx])",
                        text_align = (:left, :center),
                        fontsize = 8,
                        text_rotation = pi / 2,
                    )
                end
            end
        end
    end

    return fig

end
