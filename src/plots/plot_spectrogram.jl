export plot_spectrogram

"""
    plot_spectrogram(s_t, s_frq, s_pow; <keyword arguments>)

Plot single-channel spectrogram.

# Arguments

- `s_t::Vector{Float64}`: time
- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Array{Float64, 2}`: powers
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(s_t::Vector{Float64}, s_frq::Vector{Float64}, s_pow::Array{Float64, 2}; norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    size(s_pow, 2) == length(s_t) || throw(ArgumentError("Size of powers $(size(s_pow, 2)) and time vector $(length(s_t)) do not match."))
    size(s_pow, 1) == length(s_frq) || throw(ArgumentError("Size of powers $(size(s_pow, 1)) and frequencies vector $(length(s_frq)) do not match."))

    pal = mono == true ? :grays : :darktest
    cb_title = norm == true ? "[dB/Hz]" : "[μV^2/Hz]"

    p = Plots.heatmap(s_t,
                      s_frq,
                      s_pow,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      ylims=frq_lim,
                      xticks=_ticks(s_t),
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
    plot_spectrogram(s_ch, s_frq, s_pow; <keyword arguments>)

Plot multiple-channel spectrogram.

# Arguments

- `s_ch::Vector{String}`: channel labels
- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Array{Float64, 2}`: powers
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(s_ch::Vector{String}, s_frq::Vector{Float64}, s_pow::Array{Float64, 2}; norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    size(s_pow, 1) == length(s_ch) || throw(ArgumentError("Size of powers $(size(s_pow, 1)) and channels vector $(length(s_ch)) do not match."))
    size(s_pow, 2) == length(s_frq) || throw(ArgumentError("Size of powers $(size(s_pow, 2)) and frequencies vector $(length(s_frq)) do not match."))

    pal = mono == true ? :grays : :darktest
    cb_title = norm == true ? "[dB/Hz]" : "[μV^2/Hz]"
    
    ch = collect(1:length(s_ch)) .- 0.5
    p = Plots.heatmap(s_frq,
                      ch,
                      s_pow,
                      xlabel=xlabel,
                      xticks=_ticks(s_frq),
                      ylabel=ylabel,
                      yticks=(ch, s_ch),
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
- `epoch::Int64`: epoch to display
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or grey palette
- `markers::Bool`: draw markers if available
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(obj::NeuroAnalyzer.NEURO; epoch::Union{Int64, AbstractRange}, channel::Union{Int64, Vector{Int64}, <:AbstractRange}, norm::Bool=true, method::Symbol=:standard, nt::Int64=8, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, kwargs...)

    _check_var(method, [:standard, :stft, :mt, :mw], "method")

    _check_epochs(obj, epoch)
    _check_channels(obj, channel)

    clabels = labels(obj)[channel]
    length(channel) == 1 && (clabels = [clabels])

    # get frequency range
    fs = sr(obj)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_lim = tuple_order(frq_lim)
    (frq_lim[1] < 0 || frq_lim[2] > fs / 2) && throw(ArgumentError("frq_lim must be ≥ 0 and ≤ $(fs / 2)."))

    # calculate spectrogram
    signal = obj.data[channel, :, epoch]
    length(channel) > 1 && length(signal) / length(channel) < 4 * sr(obj) && throw(ArgumentError("For multi-channel plot, signal length must be ≥ 4 × sampling rate (4 × $(sr(obj)) samples)."))

    # get time vector
    _, t_s1, _, t_s2 = _convert_t(obj.epoch_time[1], obj.epoch_time[end])

    if length(channel) == 1
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        xlabel == "default" && (xlabel = "Time [s]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $(_channel2channel_name(channel)), epoch: $epoch, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f, s_t = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=false, st=false, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f, s_t = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=true, st=false, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            s_p, s_f, s_t = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=false, st=true, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(short-time Fourier transform)")
        elseif method === :mw
            _, s_p, _, s_f = NeuroAnalyzer.wspectrogram(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(s_t, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, kwargs=kwargs)

        # plot markers if available
        # TODO: draw markers length
        if markers == true && obj.header.has_markers == true
            markers_pos = obj.markers[!, :start] ./ sr(obj)
            markers_desc = obj.markers[!, :description]
            p = Plots.vline!(markers_pos,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:black,
                             label=false)
            for idx in eachindex(markers_desc)
                p = Plots.plot!(annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
            end
        end

    else
        ylabel == "default" && (ylabel = "Channel")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channels: $(_channel2channel_name(channel)), epoch: $epoch, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=true, nt=nt)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            _info("Method :stft is not available for multi-channel spectrogram, using standard periodogram.")
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mw
            s_p, s_f = mwpsd(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_f = linspace(0, frq_lim[2], size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(clabels, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, kwargs=kwargs)
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
- `epoch::Int64`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[component: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or grey palette
- `markers::Bool`: draw markers if available
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; epoch::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, norm::Bool=true, method::Symbol=:standard, nt::Int64=8, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, kwargs...)

    _check_var(method, [:standard, :stft, :mt, :mw], "method")

    _check_epochs(obj, epoch)

    # select component c_idxs, default is all c_idxs
    typeof(c) == Symbol && (c = _get_component(obj, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    clabels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (clabels = [clabels])

    # get frequency range
    fs = sr(obj)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_lim = tuple_order(frq_lim)
    (frq_lim[1] < 0 || frq_lim[2] > fs / 2) && throw(ArgumentError("frq_lim must be ≥ 0 and ≤ $(fs / 2)."))

    # calculate spectrogram
    signal = c[c_idx, :, epoch]
    length(c_idx) > 1 && length(signal) / length(c_idx) < 4 * sr(obj) && throw(ArgumentError("For multi-channel plot, signal length must be ≥ 4 × sampling rate (4 × $(sr(obj)) samples)."))

    # get time vector
    _, t_s1, _, t_s2 = _convert_t(obj.epoch_time[1], obj.epoch_time[end])

    if length(c_idx) == 1
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        xlabel == "default" && (xlabel = "Time [s]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f, s_t = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=false, st=false, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f, s_t = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=true, st=false, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            s_p, s_f, s_t = NeuroAnalyzer.spectrogram(signal, fs=fs, norm=false, mt=false, st=true, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(short-time Fourier transform)")
        elseif method === :mw
            _, s_p, _, s_f = NeuroAnalyzer.wspectrogram(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(s_t, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, kwargs=kwargs)

        # plot markers if available
        # TODO: draw markers length
        if markers == true && obj.header.has_markers == true
            markers_pos = obj.markers[!, :start] ./ sr(obj)
            markers_desc = obj.markers[!, :description]
            p = Plots.vline!(markers_pos,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:black,
                             label=false)
            for idx in eachindex(markers_desc)
                p = Plots.plot!(annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
            end
        end

    else
        ylabel == "default" && (ylabel = "Component")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[components: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=true, nt=nt)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            _info("Method :stft is not available for multi-channel spectrogram, using standard periodogram.")
            s_p, s_f = psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mw
            s_p, s_f = mwpsd(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_f = linspace(0, frq_lim[2], size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(clabels, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, kwargs=kwargs)
    end

    Plots.plot(p)

    return p
end
