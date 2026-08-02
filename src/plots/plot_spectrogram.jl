export plot_spectrogram
export plot_spectrogram_topo

"""
    plot_spectrogram(st, sf, sp; <keyword arguments>)

Plot a single-channel spectrogram (time vs. frequency).

# Arguments

- `st::Vector{Float64}`: vector of time values in seconds
- `sf::Vector{Float64}`: vector of frequency values in Hz
- `sp::Matrix{Float64}`: spectrogram power values
- `db::Bool=true`: whether to display power values in decibels
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `units::String=""`: power units
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing
- `cb::Bool=true`: if `true`, show color bar
- `cb_title::String=""`: colorbar title
- `threshold::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: threshold for marking regions
    - if `Real`, use a single threshold value
    - if `Tuple{Real, Real}`, use a range for `:in` or `:bin` thresholding
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: values equal to threshold
    - `:neq`: values not equal to threshold
    - `:geq`: values ≥ threshold
    - `:leq`: values ≤ threshold
    - `:g`: values > threshold
    - `:l`: values < threshold
    - `:in`: values in the threshold range (inclusive)
    - `:bin`: values in the threshold range (exclusive)

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_spectrogram(
    st::Vector{Float64},
    sf::Vector{Float64},
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
    ks::Int64 = 3,
    cb::Bool = true,
    cb_title::String = "",
    threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    threshold_type::Symbol = :neq,
)::GLMakie.Figure
    # validate
    size(sp, 2) == length(st) || throw(
        ArgumentError(
            "Size of powers ($(size(sp, 2))) and time vector ($(length(st))) do not match.",
        ),
    )
    size(sp, 1) == length(sf) || throw(
        ArgumentError(
            "Size of powers ($(size(sp, 1))) and frequencies vector ($(length(sf))) do not match.",
        ),
    )
    ks > 0 || throw(ArgumentError("ks must be ≥ 1."))
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    # set color palette
    pal = mono ? :grays : :darktest

    # apply Gaussian filter if requested
    smooth && (sp = imfilter(sp, Kernel.gaussian(ks)))

    # transpose for GLMakie heatmap (expects x columns, y rows)
    sp = sp'

    # prepare and apply thresholding mask
    if !isnothing(threshold)
        sp_threshold = deepcopy(sp)
        _, bm = seg_extract(sp; threshold = threshold, threshold_type = threshold_type)
        sp_threshold[.!bm] .= NaN
    end

    # prepare log-scaled frequencies axis
    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        flim = (sf[2], flim[2])
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_spectrogram()")
    plot_size = (1200, 800)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
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
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, (st[1], st[end]))
    GLMakie.ylims!(ax, flim)
    _style_axis!(ax)

    # draw spectrogram
    if !isnothing(threshold)
        hm = GLMakie.heatmap!(
            ax,
            st,
            sf,
            sp_threshold;
            colorrange = extrema(sp[.!isnan.(sp)]),
            colormap = pal,
        )
    else
        hm = GLMakie.heatmap!(ax, st, sf, sp; colormap = pal)
    end

    # draw colorbar if requested
    cb && GLMakie.Colorbar(fig[1, 2], hm; label = cb_title, labelsize = 16)

    return fig
end

"""
    plot_spectrogram(sf, sp; <keyword arguments>)

Plot multiple-channel spectrogram.

# Arguments

- `sf::Vector{Float64}`: vector of frequency values in Hz
- `sp::Matrix{Float64}`: spectrogram power values
- `clabels::Vector{String}=string.(1:size(sp, 1))`: channel labels
- `db::Bool=true`: whether powers are normalized to dB
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `units::String=""`: power units
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing
- `cb::Bool=true`: if `true`, show color bar
- `cb_title::String=""`: colorbar title
- `threshold::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: threshold for marking regions
    - if `Real`, use a single threshold value
    - if `Tuple{Real, Real}`, use a range for `:in` or `:bin` thresholding
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: values equal to threshold
    - `:neq`: values not equal to threshold
    - `:geq`: values ≥ threshold
    - `:leq`: values ≤ threshold
    - `:g`: values > threshold
    - `:l`: values < threshold
    - `:in`: values in the threshold range (inclusive)
    - `:bin`: values in the threshold range (exclusive)

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_spectrogram(
    sf::Vector{Float64},
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
    ks::Int64 = 3,
    cb::Bool = true,
    cb_title::String = "",
    threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    threshold_type::Symbol = :neq,
)::GLMakie.Figure
    # validate
    size(sp, 1) == length(clabels) || throw(
        ArgumentError(
            "Size of powers ($(size(sp, 1))) and channels vector ($(length(clabels))) do not match.",
        ),
    )
    size(sp, 2) == length(sf) || throw(
        ArgumentError(
            "Size of powers ($(size(sp, 2))) and frequencies vector ($(length(sf))) do not match.",
        ),
    )
    ks > 0 || throw(ArgumentError("ks must be ≥ 1."))
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    # set color palette
    pal = mono ? :grays : :darktest

    # apply Gaussian filter if requested
    smooth && (sp = imfilter(sp, Kernel.gaussian(ks)))

    # prepare log-scaled frequencies axis
    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        flim = (sf[2], flim[2])
    end

    # channel labels
    ch = collect(eachindex(clabels)) .- 0.5
    ch_n = length(ch)
    reverse!(sp; dims = 1)

    # prepare plot
    GLMakie.activate!(; title = "plot_spectrogram()")
    plot_size = (1200, 800)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
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
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, flim)
    _style_axis!(ax)

    hm = GLMakie.heatmap!(ax, sf, ch, sp'; colormap = pal)

    # draw thresholded region
    if !isnothing(threshold)
        _, bm = seg_extract(sp; threshold = threshold, threshold_type = threshold_type)
        reg = ones(size(sp)) .* minimum(sp)
        reg[bm] .= maximum(sp)
        GLMakie.contour!(ax, sf, ch, reg'; levels = 1, color = :black, linewidth = 2)
    end

    # draw colorbar
    cb && GLMakie.Colorbar(fig[1, 2], hm; label = cb_title, labelsize = 16)

    return fig
end

"""
    plot_spectrogram_topo(locs, st, sf, sp; <keyword arguments>)

Plot a topographical map of spectrogram data across channel locations with customizable visualization.

# Arguments

- `locs::DataFrame`: channel location data
- `st::Vector{Float64}`: time points in seconds corresponding to spectrogram columns
- `sf::Vector{Float64}`: frequency values in Hz corresponding to spectrogram rows
- `sp::Array{Float64, 3}`: spectrogram power values
- `db::Bool=true`: whether powers are normalized to dB
- `flim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limits for the plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `units::String=""`: power units
- `cb::Bool=true`: if `true`, show color bar
- `cb_title::String=""`: colorbar title
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing
- `mono::Bool=false`: unused, for compatibility only
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: if `true`, draw head outline

# Returns

- `GLMakie.Figure`: the plotted figure
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
    ks::Int64 = 3,
    mono::Bool = true,
    frq::Symbol = :lin,
    cart::Bool = false,
    head::Bool = true,
)::GLMakie.Figure
    # validate
    size(sp, 3) == DataFrames.nrow(locs) || throw(
        ArgumentError(
            "Size of powers ($(size(sp, 3))) and number of locs ($(DataFrames.nrow(locs))) do not match.",
        ),
    )
    size(sp, 2) == length(st) || throw(
        ArgumentError(
            "Size of powers ($(size(sp, 2))) and time vector ($(length(st))) do not match.",
        ),
    )
    size(sp, 1) == length(sf) || throw(
        ArgumentError(
            "Size of powers ($(size(sp, 1))) and frequencies vector ($(length(sf))) do not match.",
        ),
    )
    ks > 0 || throw(ArgumentError("ks must be ≥ 1."))
    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    # set color palette
    pal = mono ? :grays : :darktest

    # prepare log-scaled frequencies axis
    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        flim = (sf[2], flim[2])
    end

    # number of channels
    ch_n = size(sp, 1)

    # plot parameters
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
    if cart
        loc_x = locs.loc_x
        loc_y = locs.loc_y
    else
        loc_x = zeros(size(locs, 1))
        loc_y = zeros(size(locs, 1))
        for idx in axes(locs, 1)
            loc_x[idx], loc_y[idx] =
                pol2cart(locs.loc_radius[idx], locs.loc_theta[idx])
        end
    end

    # apply Gaussian filter if requested
    smooth && (sp = imfilter(sp, Kernel.gaussian(ks)))

    # prepare spectrogram plots
    pp_vec      = GLMakie.Figure[]
    pp_full_vec = GLMakie.Figure[]
    for idx in axes(sp, 3)
        pp = GLMakie.Figure(; size = marker_size, figure_padding = 0)
        ax = GLMakie.Axis(
            pp[1, 1];
            xlabel           = "",
            ylabel           = "",
            aspect           = nothing,
            title            = locs[idx, :label],
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
        )
        hidedecorations!(ax)
        GLMakie.xlims!(ax, flim)
        ax.titlesize = 8
        GLMakie.heatmap!(ax, sf, st, sp[:, :, idx]'; colormap = pal)
        push!(pp_vec, pp)

        pp_full = plot_spectrogram(
            st, sf, sp[:, :, idx];
            db       = db,
            frq      = frq,
            flim     = flim,
            xlabel   = xlabel,
            ylabel   = ylabel,
            title    = locs[idx, :label] * ": " * title,
            mono     = mono,
            units    = units,
            smooth   = smooth,
            ks       = ks,
            cb       = cb,
            cb_title = cb_title,
        )
        push!(pp_full_vec, pp_full)
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_spectrogram_topo()")
    fig = GLMakie.Figure(; size = plot_size, figure_padding = 0)
    ax  = GLMakie.Axis(
    fig[1, 1];
    xlabel = "",
    ylabel = "",
    title  = title,
    aspect = 1,
    _AXIS_LOCK_KWARGS...
)
    GLMakie.xlims!(ax, (-xl, xl))
    GLMakie.ylims!(ax, (-yl, yl))
    hidespines!(ax)
    hidedecorations!(ax)
    ax.titlesize = 18

    # draw head outline
    head && _draw_head_outline!(ax; lw = 3)

    for idx in axes(sp, 3)
        io = IOBuffer()
        show(io, MIME"image/png"(), pp_vec[idx])
        pp = FileIO.load(io)
        GLMakie.scatter!(
            loc_x[idx],
            loc_y[idx];
            marker      = pp,
            markersize  = marker_size,
            markerspace = :pixel,
        )
    end

    # spectrogram positions
    loc_x_range = [(loc_x[idx] - 0.15, loc_x[idx] + 0.15) for idx in eachindex(loc_x)]
    loc_y_range = [(loc_y[idx] - 0.1, loc_y[idx] + 0.1) for idx in eachindex(loc_y)]

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

Plot a spectrogram or scalogram using various time-frequency analysis methods.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `seg::Tuple{Real, Real}=(0, 10)`: time segment to analyze (from, to) in seconds; default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}=datatype(obj)`: channel name(s)
- `db::Bool=true`: if `true`, normalize powers to dB; for CWT scaleogram: normalize to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal
- `method::Symbol=:stft`: spectrogram estimation method:
    - `:stft`: short-time Fourier
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
    - `:hht`: Hilbert-Huang transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `gw::Real=10`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: y-axis frequency limits (min, max) in Hz
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `markers::Bool=true`: if `true`, draw markers if available
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing
- `cb::Bool=true`: if `true`, show color bar
- `threshold::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: threshold for marking regions
    - if `Real`, use a single threshold value
    - if `Tuple{Real, Real}`, use a range for `:in` or `:bin` thresholding
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: values equal to threshold
    - `:neq`: values not equal to threshold
    - `:geq`: values ≥ threshold
    - `:leq`: values ≤ threshold
    - `:g`: values > threshold
    - `:l`: values < threshold
    - `:in`: values in the threshold range (inclusive)
    - `:bin`: values in the threshold range (exclusive)
- `type::Symbol=:normal`: plot type:
    - `:normal`: standard spectrogram
    - `:topo`: topographical map of spectrograms
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: if `true`, draw head outline

# Returns

- `GLMakie.Figure`: the plotted figure
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
    wt::T = wavelet(Morlet(2π); β = 2),
    frq::Symbol = :lin,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    mono::Bool = false,
    markers::Bool = true,
    smooth::Bool = false,
    ks::Int64 = 3,
    cb::Bool = true,
    threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    threshold_type::Symbol = :neq,
    type::Symbol = :normal,
    cart::Bool = false,
    head::Bool = true,
)::GLMakie.Figure where {T <: CWT}
    # validate
    _check_var(type, [:normal, :topo], "type")
    _check_var(method, [:stft, :mt, :mw, :gh, :cwt, :hht], "method")
    ks > 0 || throw(ArgumentError("ks must be ≥ 1."))

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    if method === :cwt
        if type === :normal
            length(ch) == 1 ||
                throw(ArgumentError("For :cwt method only one channel must be selected."))
        end
    end
    if type === :topo
        method !== :hht ||
            throw(ArgumentError("For :hht method topographical map is not available."))
    end
    length(ch) == 1 && (ch = ch[1])

    # number of channels
    ch_n = length(ch)

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
    _check_tuple(flim, (0, sr(obj) / 2), "flim")

    # factor out the repeated ep suffix pattern
    ep_suffix = ep != 0 ? "\n[epoch: $ep]" : ""

    # calculate spectrogram
    if length(ch) == 1 || type === :topo
        if method === :stft
            spec_data = NeuroAnalyzer.spectrogram(
                signal; fs = fs, db = false, method = :stft,
                wlen = wlen, woverlap = woverlap, w = w,
            )
            sp, sf, st = spec_data.p, spec_data.f, spec_data.t
            title == "default" && (title = "Spectrogram (short-time Fourier)$ep_suffix")

        elseif method === :mt
            spec_data = NeuroAnalyzer.spectrogram(
                signal; fs = fs, db = false, method = :mt,
                nt = nt, wlen = wlen, woverlap = woverlap, w = w,
            )
            sp, sf, st = spec_data.p, spec_data.f, spec_data.t
            title == "default" && (title = "Spectrogram (multi-tapered)$ep_suffix")

        elseif method === :mw
            spec_data =
                NeuroAnalyzer.mwspectrogram(signal; fs = fs, ncyc = ncyc, db = false, w = w)
            sp, sf, st = spec_data.p, spec_data.f, spec_data.t
            title == "default" && (title = "Spectrogram (Morlet wavelet)$ep_suffix")

        elseif method === :gh
            spec_data =
                NeuroAnalyzer.ghtspectrogram(signal; fs = fs, db = false, gw = gw, w = w)
            sp, sf, st = spec_data.p, spec_data.f, spec_data.t
            title == "default" && (title = "Spectrogram (Gaussian-Hilbert)$ep_suffix")

        elseif method === :cwt
            spec_data = NeuroAnalyzer.cwtspectrogram(signal; fs = fs, wt = wt)
            sp, sf, st = spec_data.m, spec_data.f, spec_data.t
            sf[1] > flim[1] && (flim = (sf[1], flim[2]))
            sf[end] < flim[2] && (flim = (flim[1], sf[end]))
            title == "default" && (title = "CWT Scaleogram$ep_suffix")

        elseif method === :hht
            imf = emd(signal, t)
            spec_data =
                NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], t; fs = fs, db = false)
            sp, sf, st = spec_data.p, spec_data.f, spec_data.t
            title == "default" && (title = "Spectrogram (Hilbert-Huang)$ep_suffix")
        end
    elseif length(ch) > 1 && type === :normal
        if method === :stft
            psd_data = psd(
                signal; fs = fs, db = db, method = :stft,
                nt = nt, wlen = wlen, woverlap = woverlap, w = w,
            )
            sp, sf = psd_data.p, psd_data.f
            title == "default" && (title = "Spectrogram (short-time Fourier)$ep_suffix")

        elseif method === :mt
            psd_data = psd(
                signal; fs = fs, db = db, method = :mt,
                nt = nt, wlen = wlen, woverlap = woverlap, w = w,
            )
            sp, sf = psd_data.p, psd_data.f
            title == "default" && (title = "Spectrogram (multi-tapered)$ep_suffix")

        elseif method === :mw
            psd_data = psd(signal; fs = fs, db = db, method = :mw, w = w, ncyc = ncyc)
            sp, sf = psd_data.p, psd_data.f
            title == "default" && (title = "Spectrogram (Morlet wavelet)$ep_suffix")

        elseif method === :gh
            psd_data = psd(signal; fs = fs, db = db, method = :gh, w = w, gw = gw)
            sp, sf = psd_data.p, psd_data.f
            title == "default" && (title = "Spectrogram (Gaussian-Hilbert)$ep_suffix")

        elseif method === :cwt
            psd_data = psd(signal; fs = fs, method = :cwt, wt = wt)
            sp, sf = psd_data.p, psd_data.f
            sf[1] > flim[1] && (flim = (sf[1], flim[2]))
            sf[end] < flim[2] && (flim = (flim[1], sf[end]))
            title == "default" && (title = "CWT Scaleogram$ep_suffix")

        elseif method === :hht
            imf = emd(signal[1, :], t)
            spec_data =
                NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], t; fs = fs, db = db)
            sp_tmp, sf = spec_data.p, spec_data.f
            sp = zeros(size(signal, 1), length(sp_tmp))
            sp[1, :] = sp_tmp
            for idx in axes(signal, 1)[(begin + 1):end]
                imf = emd(signal[idx, :], t)
                sp[idx, :] =
                    NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], t; fs = fs, db = db).p
            end
            title == "default" && (title = "Spectrogram (Hilbert-Huang)$ep_suffix")
        end
    end

    # frequency limit
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

    # set plot labels
    cb_title = method === :cwt ? "Magnitude" : "Power"
    method !== :cwt && (cb_title *= db ? " [dB $units^2/Hz]" : " [$units^2/Hz]")

    if length(ch) == 1 && type === :normal
        xlabel == "default" && (xlabel = "Time [s]")
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        fig = plot_spectrogram(
            st, sf, sp;
            db             = db,
            frq            = frq,
            flim           = flim,
            xlabel         = xlabel,
            ylabel         = ylabel,
            title          = title,
            mono           = mono,
            units          = units,
            smooth         = smooth,
            ks             = ks,
            cb             = cb,
            cb_title       = cb_title,
            threshold      = threshold,
            threshold_type = threshold_type,
        )

    elseif length(ch) > 1 && type === :normal
        ylabel == "default" && (ylabel = "")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        fig = plot_spectrogram(
            sf, sp;
            clabels        = clabels,
            db             = db,
            frq            = frq,
            flim           = flim,
            xlabel         = xlabel,
            ylabel         = ylabel,
            title          = title,
            mono           = mono,
            units          = units,
            smooth         = smooth,
            ks             = ks,
            cb             = cb,
            cb_title       = cb_title,
            threshold      = threshold,
            threshold_type = threshold_type,
        )

    elseif type === :topo
        xlabel == "default" && (xlabel = "Time [s]")
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        length(unique(obj.header.recording[:channel_type][ch])) == 1 || throw(
            ArgumentError(
                "For multi-channel topo plot all channels must be of the same type.",
            ),
        )
        _has_locs(obj)
        chs  = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        ndims(sp) == 1 && (sp = reshape(sp, 1, 1, length(sp)))
        fig = plot_spectrogram_topo(
            locs, st, sf, sp;
            frq      = frq,
            flim     = flim,
            xlabel   = xlabel,
            ylabel   = ylabel,
            title    = title,
            mono     = mono,
            units    = units,
            cart     = cart,
            smooth   = smooth,
            ks       = ks,
            cb       = cb,
            cb_title = cb_title,
            head     = head,
        )
    end

    # draw markers if available — only supported for single-channel plots
    if ch_n == 1 && markers && _has_markers(obj)
        for idx in eachindex(obj.markers[!, :start])
            mpos = obj.markers[idx, :start]
            if _in(mpos, (obj.time_pts[1], obj.time_pts[end]))
                GLMakie.vlines!(
                    fig[1, 1], mpos;
                    linestyle = :dash,
                    linewidth = 1,
                    color     = :black,
                )
                GLMakie.textlabel!(
                    fig[1, 1],
                    (mpos + 0.07, 0.97 * minimum(obj.data[ch, :, :]));
                    text = "$(obj.markers[idx, :id]) / $(obj.markers[idx, :value])",
                    text_align = (:left, :center),
                    fontsize = 8,
                    text_rotation = pi / 2,
                )
            end
        end
    end

    return fig
end
