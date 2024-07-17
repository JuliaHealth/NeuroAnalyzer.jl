export plot_spectrogram

"""
    plot_spectrogram(st, sf, sp; <keyword arguments>)

Plot single-channel spectrogram.

# Arguments

- `st::Vector{Float64}`: time
- `sf::Vector{<:Real}`: frequencies
- `sp::Array{Float64, 2}`: powers
- `db::Bool=true`: whether powers are normalized to dB
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `units::String=""`
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(st::Vector{Float64}, sf::Vector{<:Real}, sp::Array{Float64, 2}; db::Bool=true, frq::Symbol=:lin, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, units::String="", smooth::Bool=false, n::Int64=3, kwargs...)

    @assert size(sp, 2) == length(st) "Size of powers $(size(sp, 2)) and time vector $(length(st)) do not match."
    @assert size(sp, 1) == length(sf) "Size of powers $(size(sp, 1)) and frequencies vector $(length(sf)) do not match."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest
    cb_title = db ? "[dB $units^2/Hz]" : "[$units^2/Hz]"

    if smooth
        sp = imfilter(sp, Kernel.gaussian(n))
    end

    if frq === :lin
        ysc = :identity
        yt = _ticks(frq_lim)
    else
        if frq_lim[1] == 0
            frq_lim = (0.001, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.001 Hz")
            sf[1] == 0 && (sf[1] = 0.001)
            yt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        else
            yt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        end
        ysc = :log10
    end

    p = Plots.heatmap(st,
                      sf,
                      sp,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      ylims=frq_lim,
                      xticks=_ticks(st),
                      yticks=yt,
                      yscale=ysc,
                      title=title,
                      size=(1200, 800),
                      margins=20Plots.px,
                      seriescolor=pal,
                      background_color=:black,
                      background_color_outside=:white,
                      foreground_color=:black,
                      colorbar_title=cb_title,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      xtickfontsize=6,
                      ytickfontsize=6;
                      kwargs...)

    return p

end

"""
    plot_spectrogram(sch, sf, sp; <keyword arguments>)

Plot multiple-channel spectrogram.

# Arguments

- `sch::Vector{String}`: channel labels
- `sf::Vector{<:Real}`: frequencies
- `sp::Array{Float64, 2}`: powers
- `db::Bool=true`: whether powers are normalized to dB
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `units::String=""`
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(sch::Vector{String}, sf::Vector{<:Real}, sp::Array{Float64, 2}; db::Bool=true, frq::Symbol=:lin, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, units::String="", smooth::Bool=false, n::Int64=3, kwargs...)

    @assert size(sp, 1) == length(sch) "Size of powers $(size(sp, 1)) and channels vector $(length(sch)) do not match."
    @assert size(sp, 2) == length(sf) "Size of powers $(size(sp, 2)) and frequencies vector $(length(sf)) do not match."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest
    cb_title = db ? "[dB $units^2/Hz]" : "[$units^2/Hz]"

    if smooth
        sp = imfilter(sp, Kernel.gaussian(n))
    end

    if frq === :lin
        xsc = :identity
        xt = _ticks(frq_lim)
    else
        if frq_lim[1] == 0
            frq_lim = (0.001, frq_lim[2])
            _warn("Lower frequency bound truncated to 0.001 Hz")
            sf[1] == 0 && (sf[1] = 0.001)
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        else
            xt = (round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3), string.(round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=3)))
        end
        xsc = :log10
    end

    ch = collect(eachindex(sch)) .- 0.5
    p = Plots.heatmap(sf,
                      ch,
                      sp,
                      xlabel=xlabel,
                      xticks=xt,
                      xlims=frq_lim,
                      xscale=xsc,
                      ylabel=ylabel,
                      yticks=(ch, sch),
                      title=title,
                      size=(1200, 800),
                      margins=20Plots.px,
                      seriescolor=pal,
                      background_color=:black,
                      background_color_outside=:white,
                      foreground_color=:black,
                      colorbar_title=cb_title,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      xtickfontsize=6,
                      ytickfontsize=6;
                      kwargs...)

    return p

end

"""
    plot_spectrogram(obj; <keyword arguments>)

Plots spectrogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
- `db::Bool=true`: normalize powers to dB; for CWT scaleogram: normalize to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal
- `method::Symbol=:stft`: method of calculating spectrogram:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `wt<:CWT=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: y-axis limits
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or gray palette
- `markers::Bool`: draw markers if available
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, ch::Union{Int64, Vector{Int64}, <:AbstractRange}, db::Bool=true, method::Symbol=:stft, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128), frq::Symbol=:lin, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), ncyc::Union{Int64, Tuple{Int64, Int64}}=32, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, smooth::Bool=false, n::Int64=3, kwargs...) where {T <: CWT}

    @assert seg[1] != seg[2] "Signal is too short for analysis."

    _check_var(method, [:stft, :mt, :mw, :gh, :cwt], "method")
    _check_channels(obj, ch)

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

    # set units
    units = _ch_units(obj, ch[1])

    clabels = labels(obj)[ch]
    ch isa Int64 && (clabels = [clabels])

    # get frequency range
    fs = sr(obj)
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    # calculate spectrogram
    length(ch) > 1 && @assert length(signal) / length(ch) >= 4 * sr(obj) "For multi-channel plot, signal length must be ≥ 4 × sampling rate (4 × $(sr(obj)) samples)."

    if ch isa Int64
        xlabel == "default" && (xlabel = "Time [s]")
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        if method === :stft
            sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=false, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Spectrogram (short-time Fourier transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=false, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Spectrogram (multi-tapered periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            _, sp, _, sf, st = NeuroAnalyzer.mwspectrogram(signal, fs=fs, ncyc=ncyc, db=false, w=w)
            title == "default" && (title = "Spectrogram (Morlet-wavelet transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, _, sf, st = NeuroAnalyzer.ghspectrogram(signal, fs=fs, db=false, gw=gw, w=w)
            title == "default" && (title = "Spectrogram (Gaussian and Hilbert transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf, st = NeuroAnalyzer.cwtspectrogram(signal, fs=fs, wt=wt, norm=db)
            sf[1] > frq_lim[1] && (frq_lim = (sf[1], frq_lim[2]))
            sf[end] < frq_lim[2] && (frq_lim = (frq_lim[1], sf[end]))
            title == "default" && (title = "CWT Scaleogram [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        end

        f1 = vsearch(frq_lim[1], sf)
        f2 = vsearch(frq_lim[2], sf)
        sf = sf[f1:f2]
        sp = sp[f1:f2, :]

        st .+= t[1]

        if method === :cwt
            p = plot_spectrogram(st, sf, sp, db=db, frq=frq, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, units=units, smooth=smooth, n=n, cb_title="Magnitude")
        else
            db && (sp = pow2db.(sp))
            p = plot_spectrogram(st, sf, sp, db=db, frq=frq, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, units=units, smooth=smooth, n=n, kwargs=kwargs)
        end

        # plot markers if available
        # TODO: draw markers length
        if markers && _has_markers(obj)
            markerspos = obj.markers[!, :start] ./ sr(obj)
            markers_desc = obj.markers[!, :description]
            p = Plots.vline!(markerspos,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:black,
                             label=false)
            for idx in eachindex(markers_desc)
                p = Plots.plot!(annotations=(markerspos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
            end
        end

    else
        ch_t = obj.header.recording[:channel_type]
        if length(ch) > 1
            ch_t_uni = unique(ch_t[ch])
            @assert length(ch_t_uni) == 1 "For multi-channel spectrogram plot, all channels should be of the same type."
        end
        ylabel == "default" && (ylabel = "Channel")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if method === :stft
            sp, sf = psd(signal, fs=fs, db=db, method=:stft, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Spectrogram (short-time Fourier transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channels: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Spectrogram (multi-tapered periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channels: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = psd(signal, fs=fs, db=db, method=:mw, w=w, ncyc=ncyc)
            title == "default" && (title = "Spectrogram (Morlet-wavelet transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channels: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, sf = psd(signal, fs=fs, db=db, method=:gh, w=w, gw=gw)
            title == "default" && (title = "Spectrogram (Gaussian and Hilbert transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channels: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf = psd(signal, fs=fs, db=db, method=:cwt, wt=wt)
            sf[1] > frq_lim[1] && (frq_lim = (sf[1], frq_lim[2]))
            sf[end] < frq_lim[2] && (frq_lim = (frq_lim[1], sf[end]))
            title == "default" && (title = "CWT Scaleogram [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channels: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")
        end

        f1 = vsearch(frq_lim[1], sf)
        f2 = vsearch(frq_lim[2], sf)
        sf = sf[f1:f2]
        sp = sp[:, f1:f2]

        p = plot_spectrogram(clabels, sf, sp, db=db, frq=frq, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, units=units, kwargs=kwargs)
    end

    Plots.plot(p)

    return p

end

"""
    plot_spectrogram(obj, c; <keyword arguments>)

Plots spectrogram of embedded or external component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Union{Symbol, AbstractArray}`: component to plot
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
- `db::Bool=true`: normalize powers to dB; for CWT scaleogram: normalize to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal
- `method::Symbol=:stft`: method of calculating spectrogram:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `wt<:CWT=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: y-axis limits
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[component: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or gray palette
- `markers::Bool`: draw markers if available
- `units::String=""`
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; seg::Tuple{Real, Real}=(0, 10), ep::Union{Int64, AbstractRange}=1, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, db::Bool=true, method::Symbol=:stft, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128), ncyc::Union{Int64, Tuple{Int64, Int64}}=32, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, units::String="", smooth::Bool=false, n::Int64=3, kwargs...) where {T <: CWT}

    @assert seg[1] != seg[2] "Signal is too short for analysis."
    @assert n > 0 "n must be ≥ 1."

    _check_var(method, [:stft, :mt, :mw, :gh, :cwt], "method")

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
    c_idx isa Int64 && (clabels = [clabels])

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

    # get frequency range
    fs = sr(obj)
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    # calculate spectrogram
    length(c_idx) > 1 && @assert length(signal) / length(c_idx) >= 4 * sr(obj) "For multi-channel plot, signal length must be ≥ 4 × sampling rate (4 × $(sr(obj)) samples)."

    if length(c_idx) == 1
        xlabel == "default" && (xlabel = "Time [s]")
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        if method === :stft
            sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[f1:f2, :]
            title == "default" && (title = "Spectrogram (short-time Fourier transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[f1:f2, :]
            title == "default" && (title = "Spectrogram (multi-tapered periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            _, sp, _, sf = NeuroAnalyzer.mwspectrogram(signal, fs=fs, ncyc=ncyc, db=db, w=w)
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            title == "default" && (title = "Spectrogram (Morlet-wavelet transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, _, sf = NeuroAnalyzer.ghspectrogram(signal, fs=fs, db=db, gw=gw, w=w)
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            title == "default" && (title = "Spectrogram (Gaussian and Hilbert transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf = NeuroAnalyzer.cwtspectrogram(signal, fs=fs, wt=wt, norm=db)
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            sf[1] > frq_lim[1] && (frq_lim = (sf[1], frq_lim[2]))
            sf[end] < frq_lim[2] && (frq_lim = (frq_lim[1], sf[end]))
            title == "default" && (title = "CWT Scaleogram [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end

        st .+= t[1]

        p = plot_spectrogram(st, sf, sp, db=db, frq=frq, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, units=units, smooth=smooth, n=n, kwargs=kwargs)

        # plot markers if available
        # TODO: draw markers length
        if markers && _has_markers(obj)
            markerspos = obj.markers[!, :start] ./ sr(obj)
            markers_desc = obj.markers[!, :description]
            p = Plots.vline!(markerspos,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:black,
                             label=false)
            for idx in eachindex(markers_desc)
                p = Plots.plot!(annotations=(markerspos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
            end
        end

    else
        ylabel == "default" && (ylabel = "Component")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if method === :stft
            sp, sf = psd(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[:, f1:f2]
            title == "default" && (title = "Spectrogram (short-time Fourier transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[:, f1:f2]
            title == "default" && (title = "Spectrogram (multi-tapered periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = mwpsd(signal, fs=fs, ncyc=ncyc, db=db, w=w)
            sf = linspace(0, frq_lim[2], size(sp, 2))
            title == "default" && (title = "Spectrogram (Morlet-wavelet transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, _, sf = NeuroAnalyzer.ghspectrogram(signal, fs=fs, db=db, gw=gw, w=w)
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            title == "default" && (title = "Spectrogram (Gaussian and Hilbert transform) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            sp, sf = NeuroAnalyzer.cwtspectrogram(signal, fs=fs, wt=wt, norm=db)
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            sf[1] > frq_lim[1] && (frq_lim = (sf[1], frq_lim[2]))
            sf[end] < frq_lim[2] && (frq_lim = (frq_lim[1], sf[end]))
            title == "default" && (title = "CWT Scaleogram [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end

        p = plot_spectrogram(clabels, sf, sp, db=db, frq=frq, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, units=units, smooth=smooth, n=n, kwargs=kwargs)
    end

    Plots.plot(p)

    return p

end
