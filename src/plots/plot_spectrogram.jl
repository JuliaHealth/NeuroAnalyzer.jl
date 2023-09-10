export plot_spectrogram

"""
    plot_spectrogram(st, sf, sp; <keyword arguments>)

Plot single-channel spectrogram.

# Arguments

- `st::Vector{Float64}`: time
- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 2}`: powers
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `units::String=""`
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(st::Vector{Float64}, sf::Vector{Float64}, sp::Array{Float64, 2}; norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, units::String="", kwargs...)

    @assert size(sp, 2) == length(st) "Size of powers $(size(sp, 2)) and time vector $(length(st)) do not match."
    @assert size(sp, 1) == length(sf) "Size of powers $(size(sp, 1)) and frequencies vector $(length(sf)) do not match."

    pal = mono == true ? :grays : :darktest
    cb_title = norm == true ? "[dB/Hz]" : "[$units^2/Hz]"

    p = Plots.heatmap(st,
                      sf,
                      sp,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      ylims=frq_lim,
                      xticks=_ticks(st),
                      yticks=_ticks(frq_lim),
                      title=title,
                      size=(1200, 800),
                      margins=20Plots.px,
                      seriescolor=pal,
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
- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 2}`: powers
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `units::String=""`
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(sch::Vector{String}, sf::Vector{Float64}, sp::Array{Float64, 2}; norm::Bool=true, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, units::String="", kwargs...)

    @assert size(sp, 1) == length(sch) "Size of powers $(size(sp, 1)) and channels vector $(length(sch)) do not match."
    @assert size(sp, 2) == length(sf) "Size of powers $(size(sp, 2)) and frequencies vector $(length(sf)) do not match."

    pal = mono == true ? :grays : :darktest
    cb_title = norm == true ? "[dB/Hz]" : "[$units^2/Hz]"
    
    ch = collect(eachindex(sch)) .- 0.5
    p = Plots.heatmap(sf,
                      ch,
                      sp,
                      xlabel=xlabel,
                      xticks=_ticks(sf),
                      ylabel=ylabel,
                      yticks=(ch, sch),
                      title=title,
                      size=(1200, 800),
                      margins=20Plots.px,
                      seriescolor=pal,
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
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: Use color or gray palette
- `markers::Bool`: draw markers if available
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, ch::Union{Int64, Vector{Int64}, <:AbstractRange}, norm::Bool=true, method::Symbol=:standard, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, gw::Real=5, wt::T=wavelet(Morlet(2π), β=2), frq::Symbol=:log, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, kwargs...) where {T <: CWT}

    _check_var(method, [:standard, :stft, :mt, :mw, :gh, :cwt], "method")

    _check_channels(obj, ch)

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

    # set units
    units = _set_units(obj, ch[1])

    clabels = labels(obj)[ch]
    length(ch) == 1 && (clabels = [clabels])

    # get frequency range
    fs = sr(obj)
    frq_lim = tuple_order(frq_lim)
    @assert !(frq_lim[1] < 0 || frq_lim[2] < 0 || frq_lim[1] > sr(obj) / 2 || frq_lim[2] > sr(obj) / 2) "frq_lim must be in [0, $(sr(obj) / 2)]."

    # calculate spectrogram
    length(ch) > 1 && @assert length(signal) / length(ch) >= 4 * sr(obj) "For multi-channel plot, signal length must be ≥ 4 × sampling rate (4 × $(sr(obj)) samples)."

    if length(ch) == 1
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        xlabel == "default" && (xlabel = "Time [s]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f, st = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=false, st=false, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f, st = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=true, st=false, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            s_p, s_f, st = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=false, st=true, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            title = replace(title, "method" => "(short-time Fourier transform)")
        elseif method === :mw
            _, s_p, _, s_f = NeuroAnalyzer.wspectrogram(signal, fs=fs, frq_lim=frq_lim, frq=frq, frq_n=frq_n, ncyc=ncyc, norm=false)
            st = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        elseif method === :gh
            s_p, _, s_f = NeuroAnalyzer.ghspectrogram(signal, fs=fs, frq_lim=frq_lim, norm=false, frq=frq, frq_n=frq_n, gw=gw)
            st = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(Gaussian and Hilbert transform)")
        elseif method === :cwt
            s_p, s_f = NeuroAnalyzer.cwtspectrogram(signal, fs=fs, frq_lim=frq_lim, norm=false, wt=wt)
            st = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(continuous wavelet transformation)")
        end

        st .+= t[1]

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(st, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, units=units, kwargs=kwargs)

        # plot markers if available
        # TODO: draw markers length
        if markers == true && _has_markers(obj) == true
            markers_pos = obj.markers[!, :start] ./ sr(obj)
            markers_desc = obj.markers[!, :description]
            p = Plots.vline!(markers_pos,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:black,
                             label=false)
            for idx in eachindex(markers_desc)
                p = Plots.plot!(annotations=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
            end
        end

    else
        ch_t = obj.header.recording[:channel_type]
        if length(ch) > 1
            ch_t_uni = unique(ch_t[ch])
            @assert length(ch_t_uni) == 1 "For multi-channel PSD plots all channels should be of the same type."
        end
        ylabel == "default" && (ylabel = "Channel")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channels: $(_channel2channel_name(ch)), epoch: $ep, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=false, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=true, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            s_p, s_f = psd(signal, fs=fs, norm=false, st=true, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(short-time Fourier transform)")
        elseif method === :mw
            s_p, s_f = psd_mw(signal, fs=fs, frq_lim=frq_lim, ncyc=ncyc, norm=false, frq=frq, frq_n=frq_n)
            s_f = linspace(0, frq_lim[2], size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        elseif method === :gh
            s_p, _, s_f = NeuroAnalyzer.ghspectrogram(signal, fs=fs, frq_lim=frq_lim, norm=false, gw=gw, frq=frq, frq_n=frq_n)
            st = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(Gaussian and Hilbert transform)")
        elseif method === :cwt
            s_p, s_f = NeuroAnalyzer.cwtspectrogram(signal, fs=fs, frq_lim=frq_lim, norm=false, wt=wt)
            st = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(continuous wavelet transformation)")
        end
        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(clabels, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, units=units, kwargs=kwargs)
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
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: y-axis limits
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[component: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: Use color or gray palette
- `markers::Bool`: draw markers if available
- `units::String=""`
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; seg::Tuple{Real, Real}=(0, 10), ep::Union{Int64, AbstractRange}=1, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, norm::Bool=true, method::Symbol=:standard, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), gw::Real=5, wt::T=wavelet(Morlet(2π), β=2), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, units::String="", kwargs...) where {T <: CWT}

    _check_var(method, [:standard, :stft, :mt, :mw, :gh, :cwt], "method")

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

    # select component c_idxs, default is all c_idxs
    c isa Symbol && (c = _get_component(obj, c))
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    clabels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (clabels = [clabels])

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
    frq_lim = tuple_order(frq_lim)
    @assert !(frq_lim[1] < 0 || frq_lim[2] < 0 || frq_lim[1] > sr(obj) / 2 || frq_lim[2] > sr(obj) / 2) "frq_lim must be in [0, $(sr(obj) / 2)]."

    # calculate spectrogram
    length(c_idx) > 1 && @assert length(signal) / length(c_idx) >= 4 * sr(obj) "For multi-channel plot, signal length must be ≥ 4 × sampling rate (4 × $(sr(obj)) samples)."

    if length(c_idx) == 1
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        xlabel == "default" && (xlabel = "Time [s]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f, st = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=false, st=false, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # st = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f, st = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=true, st=false, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # st = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            s_p, s_f, st = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=false, st=true, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # st = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(short-time Fourier transform)")
        elseif method === :mw
            _, s_p, _, s_f = NeuroAnalyzer.wspectrogram(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            st = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        st .+= t[1]

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(st, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, units=units, kwargs=kwargs)

        # plot markers if available
        # TODO: draw markers length
        if markers == true && _has_markers(obj) == true
            markers_pos = obj.markers[!, :start] ./ sr(obj)
            markers_desc = obj.markers[!, :description]
            p = Plots.vline!(markers_pos,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:black,
                             label=false)
            for idx in eachindex(markers_desc)
                p = Plots.plot!(annotations=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
            end
        end

    else
        ylabel == "default" && (ylabel = "Component")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=false, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=true, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            s_p, s_f = psd(signal, fs=fs, norm=false, st=true, mt=false, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(short-time Fourier transform)")
        elseif method === :mw
            s_p, s_f = psd_mw(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_f = linspace(0, frq_lim[2], size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(clabels, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, units=units, kwargs=kwargs)
    end

    Plots.plot(p)

    return p

end
