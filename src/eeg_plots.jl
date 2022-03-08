"""
    signal_plot(t, signal; <keyword arguments>)

Plot multi-channel `signal`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
- `signal::AbstractArray`
- `labels::Vector{String}=[""]`: labels vector
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::AbstractArray; labels::Vector{String}=[""], xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    channel_n = size(signal, 1)

    # reverse so 1st channel is on top
    channel_color = channel_n:-1:1
    signal = reverse(signal[:, :], dims = 1)
    s_normalized = zeros(size(signal))

    # normalize and shift so all channels are visible
    variances = var(signal, dims=2)
    mean_variance = mean(variances)
    for idx in 1:channel_n
        s = @view signal[idx, :]
        s_normalized[idx, :] = (s .- mean(s)) ./ mean_variance .+ (idx - 1)
    end

    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(floor(t[1], digits=2), ceil(t[end], digits=2)),
             xticks=floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2),
             ylims=(-0.5, channel_n-0.5),
             title=title,
             palette=:darktest,
             titlefontsize=10,
             xlabelfontsize=8,
             ylabelfontsize=8,
             xtickfontsize=8,
             ytickfontsize=8;
             kwargs...)
    for idx in 1:channel_n
        p = plot!(t,
                  s_normalized[idx, 1:length(t)],
                  label="",
                  color=channel_color[idx])
    end

    p = plot!(yticks=((channel_n - 1):-1:0, labels))

    return p
end

"""
    signal_plot(t, signal; <keyword arguments>)

Plot single-channel `signal`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
- `signal::Vector{Float64}`
- `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Amplitude [μV]"`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Vector{Float64}; ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    ylim == (0, 0) && (ylim = (floor(minimum(signal), digits=0), ceil(maximum(signal), digits=0)))
    abs(ylim[1]) > abs(ylim[2]) && (ylim = (-abs(ylim[1]), abs(ylim[1])))
    abs(ylim[1]) < abs(ylim[2]) && (ylim = (-abs(ylim[2]), abs(ylim[2])))
    ylim = tuple_order(ylim)

    hl = plot((size(signal, 2), 0), seriestype=:hline, linewidth=0.5, linealpha=0.5, linecolor=:gray, label="")
    p = plot!(t,
              signal[1:length(t)],
              color=1,
              label="",
              legend=false,
              title=title,
              xlabel=xlabel,                  
              xlims=(floor(t[1], digits=2), ceil(t[end], digits=2)),
              xticks=floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2),
              ylabel=ylabel,
              ylims=ylim,
              yguidefontrotation=0,
              yticks=[ylim[1], 0, ylim[2]],
              palette=:darktest,
              grid=false,
              titlefontsize=10,
              xlabelfontsize=8,
              ylabelfontsize=8,
              xtickfontsize=4,
              ytickfontsize=4,
              margins=10Plots.px;                  
              kwargs...)

    plot(p)

    return p
end

"""
    eeg_plot(eeg; <keyword arguments>)

Plot `eeg` channels. If signal is multi-channel, only channel amplitudes are plotted. For single-channel signal, the histogram, amplitude, power density and spectrogram are plotted.

# Arguments

- `eeg::NeuroJ.EEG`: EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epochs to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
- `offset::Int64=0`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
- `labels::Vector{String}=[""]`: channel labels vector
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `head::Bool=true`: add head with electrodes
- `hist::Symbol[:hist, :kd]=:hist`: histogram type
- `norm::Bool=true`: convert power to dB
- `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: frequency limit for PSD and spectrogram
- `kwargs`: other arguments for plot() function; <keyword arguments>

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], xlabel::String="Time [s]", ylabel::String="", title::String="", head::Bool=true, hist::Symbol=:hist, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

    hist in [:hist, :kd] || throw(ArgumentError("hist must be :hist or :kd."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))

    (frq_lim[1] < 0 || frq_lim[1] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    (frq_lim[2] < 0 || frq_lim[2] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    frq_lim == (0, 0) && (frq_lim = (0, eeg_sr(eeg) / 2))
    frq_lim = tuple_order(frq_lim)

    (epoch != 1 && (offset != 0 || len != 0)) && throw(ArgumentError("For epoch ≠ 1, offset and len must not be specified."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    (length(epoch) == 1 && (epoch < 1 || epoch > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    (length(epoch) > 1 && (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    if length(epoch) > 1
        sort!(epoch)
        (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
        len = eeg_epoch_len(eeg) * length(epoch)
        offset = eeg_epoch_len(eeg) * (epoch[1] - 1)
        epoch = 1
    end

    # default length is one epoch or 20 seconds
    if len == 0
        if eeg_epoch_len(eeg) > 20 * eeg_sr(eeg)
            len = 20 * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end
    
    # select channels, default is 1:20 or all channels
    if channel == 0
        if eeg_channel_n(eeg) >= 20
            channel = 1:20
        else
            channel = 1:eeg_channel_n(eeg)
        end
    end
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    length(channel) > 1 && sort!(channel)
    for idx in 1:length(channel)
        (channel[idx] < 1 || channel[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
    length(channel) == 1 && eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    # get epochs markers for len > epoch_len
    epoch_markers = Vector{Int64}[]
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_epoch_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> floor(Int64, offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= ceil(Int64, (offset + len) / eeg_sr(eeg))]
    else
        eeg_tmp = eeg
    end

    labels = eeg_labels(eeg_tmp)

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]
    t[1] = floor(t[1], digits=2)
    t[end] = ceil(t[end], digits=2)

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
    (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    if length(channel) == 1
        title = ""
        ylabel = "Amplitude [μV]"
        channel_name = labels[1]
        labels = [""]
        signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + length(t)), epoch]
        signal = vec(signal)
    else
        signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + length(t)), epoch]
    end

    t_1 = floor(t[1], digits=2)
    t_2 = ceil(t[end], digits=2)
    t_1 < 1.0 && (t_s1 = string(floor(t_1 * 1000, digits=2)) * " ms")
    t_1 >= 1.0 && (t_s1 = string(floor(t_1, digits=2)) * " s")
    t_2 < 1.0 && (t_s2 = string(ceil(t_2 * 1000, digits=2)) * " ms")
    t_2 >= 1.0 && (t_s2 = string(ceil(t_2, digits=2)) * " s")

    epoch_tmp = epoch
    offset > eeg_epoch_len(eeg) && (epoch_tmp = floor(Int64, offset / eeg_epoch_len(eeg)) + 1)
    title == "" && (title = "Signal\n[epoch: $(string(epoch_tmp)), time window: $t_s1:$t_s2]")

    p = signal_plot(t,
                    signal,
                    labels=labels,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    title=title;
                    kwargs...)

    # add epochs markers
    if length(epoch_markers) > 0 && len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                   linestyle=:dash,
                   linewidth=0.2,
                   linecolor=:black,
                   label="")
        if typeof(signal) == Vector{Float64}
            for idx in 1:length(epoch_markers)
                p = plot!(annotation=((epoch_markers[idx] - (length(signal) / eeg_sr(eeg) / 40)), maximum(ceil.(abs.(signal))), text("E$(string(floor(Int64, epoch_markers[idx] / (eeg_epoch_len(eeg) / eeg_sr(eeg)))))", pointsize=4, halign=:left, valign=:top)))
            end
        else
            for idx in 1:length(epoch_markers)
                p = plot!(annotation=((epoch_markers[idx] - (size(signal, 2) / eeg_sr(eeg) / 40)), length(channel)-0.5, text("E$(string(floor(Int64, epoch_markers[idx] / (eeg_epoch_len(eeg) / eeg_sr(eeg)))))", pointsize=6, halign=:left, valign=:top)))
            end
        end
    end

    if typeof(signal) == Vector{Float64}
        # cannot plot electrodes without locations
        eeg.eeg_header[:channel_locations] == false && (head = false)
        psd = eeg_plot_psd(eeg, epoch=epoch, channel=channel, len=len, offset=offset, frq_lim=frq_lim, title="PSD [dB]\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", norm=norm, legend=false)
        frq_lim == (0, 0) && (frq_lim = (0, div(eeg_sr(eeg), 2)))
        s = eeg_plot_spectrogram(eeg, epoch=epoch, channel=channel, len=len, offset=offset, frq_lim=frq_lim, title="Spectrogram [dB]\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", legend=false)
        ht_a = eeg_plot_histogram(eeg, epoch=epoch, channel=channel, len=len, offset=offset, type=hist, labels=[""], legend=false, title="Signal\nhistogram")
        _, _, _, s_phase = signal_spectrum(signal)
        ht_p = signal_plot_histogram(rad2deg.(s_phase), offset=offset, len=len, type=:kd, labels=[""], legend=false, title="Phase\nhistogram", xticks=[-180, 0, 180], linecolor=:black)
        if head == true
            hd = eeg_plot_electrodes(eeg, labels=false, selected=channel, small=true, title="Channel: $channel\nLabel: $channel_name")
            l = @layout [a{0.33h} b{0.2w}; c{0.33h} d{0.2w}; e{0.33h} f{0.2w}]
            p = plot(p, ht_a, psd, ht_p, s, hd, layout=l)
        else
            l = @layout [a{0.33h} b{0.2w}; c{0.33h} d{0.2w}; e{0.33h} _]
            p = plot(p, ht_a, psd, ht_p, s, layout=l)
        end
    end

    plot(p)

    return p
end

"""
    eeg_draw_head(p, loc_x, loc_y; head_labels, kwargs)

Draw head over a topographical plot `p`.

# Arguments

- `p::Plots.Plot{Plots.GRBackend}`: electrodes plot
- `loc_x::Vector{Float64}`: vector of x electrode position
- `loc_y::Vector{Float64}`: vector of y electrode position
- `head_labels::Bool=true`: add text labels to the plot
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_draw_head(p::Plots.Plot{Plots.GRBackend}, loc_x::Vector{Float64}, loc_y::Vector{Float64}; head_labels::Bool=true, kwargs...)

    pts = Plots.partialcircle(0, 2π, 100, maximum(loc_x))
    x, y = Plots.unzip(pts)
    x = x .* 4
    y = y .* 4
    head = Shape(x, y)
    nose = Shape([(-0.05, maximum(y)), (0, maximum(y) + 0.1 * maximum(y)), (0.05, maximum(y))])
    ear_l = Shape([(minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), -0.08), (minimum(x) + 0.1 * minimum(x), 0.08), (minimum(x), 0.1)])
    ear_r = Shape([(maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), -0.08), (maximum(x) + 0.1 * maximum(x), 0.08), (maximum(x), 0.1)])
    p = plot!(p, head, fill=nothing, label="")
    p = plot!(nose, fill=nothing, label="")
    p = plot!(ear_l, fill=nothing, label="")
    p = plot!(ear_r, fill=nothing, label="")
    if head_labels == true
        p = plot!(annotation=(0, 1 - maximum(y) / 5, text("Inion", pointsize=12, halign=:center, valign=:center)))
        p = plot!(annotation=(0, -1 - minimum(y) / 5, text("Nasion", pointsize=12, halign=:center, valign=:center)))
        p = plot!(annotation=(-1 - minimum(x) / 5, 0, text("Left", pointsize=12, halign=:center, valign=:center, rotation=90)))
        p = plot!(annotation=(1 - maximum(x) / 5, 0, text("Right", pointsize=12, halign=:center, valign=:center, rotation=-90)))
    end
    p = plot!(; kwargs...)

    return p
end

"""
    eeg_plot_filter_response(eeg; <keyword arguments>)

Plot filter response.

# Arguments

- `eeg::NeuroJ.EEG`
- `fprototype::Symbol`: filter class: :butterworth, :chebyshev1, :chebyshev2, :elliptic
- `ftype::Symbol`: filter type: :lp, :hp, :bp, :bs
- `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64`: filter order
- `rp::Union{Int64, Float64}`: dB ripple in the passband
- `rs::Union{Int64, Float64}`: dB attenuation in the stopband
- `window::window::Union{Vector{Float64}, Nothing}`: window, required for FIR filter
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_filter_response(eeg::NeuroJ.EEG; fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Tuple}, order::Int64, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, window::Union{Vector{Float64}, Nothing}=nothing, kwargs...)

    fs = eeg_sr(eeg)

    ftype in [:lp, :hp, :bp, :bs] || throw(ArgumentError("ftype must be :bp, :hp, :bp or :bs."))
    fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic] || throw(ArgumentError("fprototype must be :butterworth, :chebyshev1:, :chebyshev2 or :elliptic."))

    if ftype === :lp
        length(cutoff) != 1 && throw(ArgumentError("For :lp filter one frequency must be given."))
        responsetype = Lowpass(cutoff; fs=fs)
    elseif ftype === :hp
        length(cutoff) != 1 && throw(ArgumentError("For :hp filter one frequency must be given."))
        responsetype = Highpass(cutoff; fs=fs)
    elseif ftype === :bp
        length(cutoff) != 2 && throw(ArgumentError("For :bp filter two frequencies must be given."))
        responsetype = Bandpass(cutoff[1], cutoff[2]; fs=fs)
    elseif ftype === :bs
        length(cutoff) != 2 && throw(ArgumentError("For :bs filter two frequencies must be given."))
        cutoff = tuple_order(cutoff)
        responsetype = Bandstop(cutoff[1], cutoff[2]; fs=fs)
    end

    fprototype === :butterworth && (prototype = Butterworth(order))
    if fprototype === :chebyshev1
        (rs < 0 || rs > eeg_sr(eeg) / 2) && throw(ArgumentError("For :chebyshev1 filter rs must be > 0 and ≤ $(eeg_sr(eeg) / 2)."))
        prototype = Chebyshev1(order, rs)
    end
    if fprototype === :chebyshev2
        (rp < 0 || rp > eeg_sr(eeg) / 2) && throw(ArgumentError("For :chebyshev2 filter rp must be > 0 and ≤ $(eeg_sr(eeg) / 2)."))
        prototype = Chebyshev2(order, rp)
    end
    if fprototype === :elliptic
        (rs < 0 || rs > eeg_sr(eeg) / 2) && throw(ArgumentError("For :elliptic filter rs must be > 0 and ≤ $(eeg_sr(eeg) / 2)."))
        (rp < 0 || rp > eeg_sr(eeg) / 2) && throw(ArgumentError("For :elliptic filter rp must be > 0 and ≤ $(eeg_sr(eeg) / 2)."))
        prototype = Elliptic(order, rp, rs)
    end

    ffilter = digitalfilter(responsetype, prototype)

    H, w = freqresp(ffilter)
    # convert to dB
    H = 20 * log10.(abs.(H))
    # convert rad/sample to Hz
    w = w .* fs / 2 / pi
    x_max = w[end]
    ftype === :hp && (x_max = cutoff * 10)
    p1 = plot(w,
              H,
              title="Filter: $(titlecase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order\nFrequency response",
              xlims=(0, x_max),
              ylabel="Magnitude\n[dB]",
              xlabel="Frequency [Hz]",
              label="",
              titlefontsize=10,
              xlabelfontsize=8,
              ylabelfontsize=8,
              xtickfontsize=4,
              ytickfontsize=4)
    if length(cutoff) == 1
        p1 = plot!((0, cutoff),
                   seriestype=:vline,
                   linestyle=:dash,
                   label="")
    else
        p1 = plot!((0, cutoff[1]),
                   seriestype=:vline,
                   linestyle=:dash,
                   label="")
        p1 = plot!((0, cutoff[2]),
                   seriestype=:vline,
                   linestyle=:dash,
                   label="")
    end

    phi, w = phaseresp(ffilter)
    phi = rad2deg.(angle.(phi))
    # convert rad/sample to Hz
    w = w .* fs / 2 / pi
    x_max = w[end]
    ftype === :hp && (x_max = cutoff * 10)
    p2 = plot(w,
              phi,
              title="Phase response",
              ylims=(-180, 180),
              xlims=(0, x_max),
              ylabel="Phase\n[°]",
              xlabel="Frequency [Hz]",
              label="",
              titlefontsize=10,
              xlabelfontsize=8,
              ylabelfontsize=8,
              xtickfontsize=4,
              ytickfontsize=4)
    if length(cutoff) == 1
        p2 = plot!((0, cutoff),
                   seriestype=:vline,
                   linestyle=:dash,
                   label="")
    else
        p2 = plot!((0, cutoff[1]),
                   seriestype=:vline,
                   linestyle=:dash,
                   label="")
        p2 = plot!((0, cutoff[2]),
                   seriestype=:vline,
                   linestyle=:dash,
                   label="")
    end

    tau, w = grpdelay(ffilter)
    tau = abs.(tau)
    # convert rad/sample to Hz
    w = w .* fs / 2 / pi
    x_max = w[end]
    ftype === :hp && (x_max = cutoff * 10)
    p3 = plot(w,
              tau,
              title="Group delay",
              xlims=(0, x_max),
              ylabel="Group delay\n[samples]",
              xlabel="Frequency [Hz]",
              label="",
              titlefontsize=10,
              xlabelfontsize=8,
              ylabelfontsize=8,
              xtickfontsize=4,
              ytickfontsize=4)
    if length(cutoff) == 1
        p3 = plot!((0, cutoff),
                   seriestype=:vline,
                   linestyle=:dash,
                   label="")
    else
        p3 = plot!((0, cutoff[1]),
                   seriestype=:vline,
                   linestyle=:dash,
                   label="")
        p3 = plot!((0, cutoff[2]),
                   seriestype=:vline,
                   linestyle=:dash,
                   label="")
    end

    p = plot(p1, p2, p3, layout=(3, 1), palette=:darktest; kwargs...)

    return p
end

"""
    signal_plot_avg(t, signal; <keyword arguments>)

Plot averaged `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange`
- `signal::Matrix{Float64}`
- `norm::Bool=true`: normalize the `signal` prior to calculations
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Amplitude [μV]"`: y-axis label
- `title::String=""`: plot title
- `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_avg(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Matrix{Float64}; norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    if norm == true
        s_normalized = signal_normalize_zscore(signal)
    else
        s_normalized = signal
    end
    
    s_normalized_m, s_normalized_s, s_normalized_u, s_normalized_l = signal_ci95(s_normalized)

    ylim == (0, 0) && (ylim = (floor(minimum(s_normalized_l), digits=0), ceil(maximum(s_normalized_u), digits=0)))
    abs(ylim[1]) > abs(ylim[2]) && (ylim = (-abs(ylim[1]), abs(ylim[1])))
    abs(ylim[1]) < abs(ylim[2]) && (ylim = (-abs(ylim[2]), abs(ylim[2])))
    ylim[1] > ylim[2] && (ylim = (ylim[2], ylim[1]))

    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(floor(t[1], digits=2), ceil(t[end], digits=2)),
             xticks=floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2),
             ylims=ylim,
             title=title,
             palette=:darktest,
             titlefontsize=10,
             xlabelfontsize=8,
             ylabelfontsize=8,
             xtickfontsize=4,
             ytickfontsize=4;
             kwargs...)
    p = plot!(t,
              s_normalized_u[1:length(t)],
              fillrange = s_normalized_l,
              fillalpha = 0.35, 
              label=false,
              t=:line,
              c=:grey,
              lw=0.5)
    p = plot!(t,
              s_normalized_l[1:length(t)],
              label=false,
              t=:line,
              c=:grey,
              lw=0.5)
    p = plot!(t,
              s_normalized_m[1:length(t)],
              label=false,
              t=:line,
              c=:black)

    return p
end

"""
    eeg_plot_avg(eeg; <keyword arguments>)

Plot averaged `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`: EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
- `offset::Int64=0`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
- `norm::Bool=true`: normalize the `signal` prior to calculations
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Amplitude [μV]"`: y-axis label
- `title::String=""`: plot title
- `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
- `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: frequency limit for PSD and spectrogram
- `hist::Symbol=:hist`: histogram type: :hist, :kd
- `head::Bool=true`: add head plot
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_avg(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), hist::Symbol=:hist, head::Bool=true, kwargs...)

    (frq_lim[1] < 0 || frq_lim[1] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    (frq_lim[2] < 0 || frq_lim[2] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))
    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    (epoch != 1 && (offset != 0 || len != 0)) && throw(ArgumentError("For epoch ≠ 1, offset and len must not be specified."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    (length(epoch) == 1 && (epoch < 1 || epoch > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    (length(epoch) > 1 && (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    if length(epoch) > 1
        sort!(epoch)
        (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
        len = eeg_epoch_len(eeg) * length(epoch)
        offset = eeg_epoch_len(eeg) * (epoch[1] - 1)
        epoch = 1
    end

    # default length is one epoch or 20 seconds
    if len == 0
        if eeg_epoch_len(eeg) > 20 * eeg_sr(eeg)
            len = 20 * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end

    # select channels, default is all channels
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    channel == 0 && (channel = 1:eeg_channel_n(eeg))
    length(channel) == 1 && throw(ArgumentError("At least 2 channels are required."))
    length(channel) > 1 && sort!(channel)
    for idx in 1:length(channel)
        (channel[idx] < 1 || channel[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end

    # get epochs markers for len > epoch_len
    epoch_markers = Vector{Int64}[]
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_epoch_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> floor(Int64, offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= ceil(Int64, (offset + len) / eeg_sr(eeg))]
    else
        eeg_tmp = eeg
    end

    labels = eeg_labels(eeg_tmp)

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]
    t[1] = floor(t[1], digits=2)
    t[end] = ceil(t[end], digits=2)

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
    (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + length(t)), epoch]

    t_1 = floor(t[1], digits=2)
    t_2 = ceil(t[end], digits=2)
    t_1 < 1.0 && (t_s1 = string(floor(t_1 * 1000, digits=2)) * " ms")
    t_1 >= 1.0 && (t_s1 = string(floor(t_1, digits=2)) * " s")
    t_2 < 1.0 && (t_s2 = string(ceil(t_2 * 1000, digits=2)) * " ms")
    t_2 >= 1.0 && (t_s2 = string(ceil(t_2, digits=2)) * " s")
    epoch_tmp = epoch
    offset > eeg_epoch_len(eeg) && (epoch_tmp = floor(Int64, offset / eeg_epoch_len(eeg)) + 1)
    title == "" && (title = "Signal averaged\n[epoch: $(string(epoch_tmp)), time window: $t_s1:$t_s2]")

    p = signal_plot_avg(t,
                        signal,
                        offset=offset,
                        norm=norm,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        title=title,
                        ylim=ylim;
                        kwargs...)

    # add epochs markers
    if norm == true
        s_normalized = signal_normalize_zscore(signal)
    else
        s_normalized = signal
    end
    s_normalized_m, s_normalized_s, s_normalized_u, s_normalized_l = signal_ci95(s_normalized)
    ylim = (floor(minimum(s_normalized_l), digits=0), ceil(maximum(s_normalized_u), digits=0))
    abs(ylim[1]) > abs(ylim[2]) && (ylim = (-abs(ylim[1]), abs(ylim[1])))
    abs(ylim[1]) < abs(ylim[2]) && (ylim = (-abs(ylim[2]), abs(ylim[2])))
    if length(epoch_markers) > 0 && len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                   linestyle=:dash,
                   linewidth=0.2,
                   linecolor=:black,
                   label="")
        for idx in 1:length(epoch_markers)
            p = plot!(annotation=((epoch_markers[idx] - 1), ylim[2], text("E$(string(floor(Int64, epoch_markers[idx] / (eeg_epoch_len(eeg) / eeg_sr(eeg)))))", pointsize=4, halign=:center, valign=:top)))
        end
    end

    # cannot plot electrodes without locations
    eeg.eeg_header[:channel_locations] == false && (head = false)
    frq_lim == (0, 0) && (frq_lim = (0, div(eeg_sr(eeg), 2)))
    eeg_avg = eeg_average(eeg)
    psd = eeg_plot_psd(eeg_tmp, epoch=epoch, channel=channel, len=len, offset=offset, average=true, norm=true, frq_lim=frq_lim, title="PSD averaged with 95%CI [dB]\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", legend=false)
    s = eeg_plot_spectrogram(eeg_avg, epoch=epoch, channel=1, len=len, offset=offset, frq_lim=frq_lim, title="Spectrogram averaged\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", legend=false)
    ht_a = eeg_plot_histogram(eeg_avg, epoch=epoch, channel=1, len=len, offset=offset, type=hist, labels=[""], legend=false, title="Signal\nhistogram")
    _, _, _, s_phase = signal_spectrum(s_normalized_m)
    ht_p = signal_plot_histogram(rad2deg.(s_phase), offset=offset, len=len, type=:kd, labels=[""], legend=false, title="Phase\nhistogram", xticks=[-180, 0, 180], linecolor=:black)
    if head == true
        if collect(channel[1]:channel[end]) == channel
            channel_list = string(channel[1]) * ":" * string(channel[end])
        else
            channel_list = "" 
            for idx in 1:(length(channel) - 1)
                channel_list *= string(channel[idx])
                channel_list *= ", "
            end
            channel_list *= string(channel[end])
        end
        hd = eeg_plot_electrodes(eeg, labels=false, selected=channel, small=true, title="Channels\n$channel_list")
        l = @layout [a{0.33h} b{0.2w}; c{0.33h} d{0.2w}; e{0.33h} f{0.2w}]
        p = plot(p, ht_a, psd, ht_p, s, hd, layout=l)
    else
        l = @layout [a{0.33h} b{0.2w}; c{0.33h} d{0.2w}; e{0.33h} _]
        p = plot(p, ht_a, psd, ht_p, s, layout=l)
    end

    plot(p)

    return p
end

"""
    signal_plot_butterfly(t, signal; <keyword arguments>)

Butterfly plot of `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange`
- `signal::Matrix{Float64}`
- `labels::Vector{String}=[""]`: channel labels vector
- `norm::Bool=true`: normalize the `signal` prior to calculations
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Amplitude [μV]"`: y-axis label
- `title::String=""`: plot title
- `ylim::Tuple`: y-axis limits, default (0, 0)
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_butterfly(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    channel_n = size(signal, 1)

    if norm == true
        s_normalized = signal_normalize_zscore(reshape(signal, size(signal, 1), size(signal, 2), 1))
    else
        s_normalized = signal
    end

    ylim == (0, 0) && (ylim = (floor(minimum(s_normalized), digits=0), ceil(maximum(s_normalized), digits=0)))
    abs(ylim[1]) > abs(ylim[2]) && (ylim = (-abs(ylim[1]), abs(ylim[1])))
    abs(ylim[1]) < abs(ylim[2]) && (ylim = (-abs(ylim[2]), abs(ylim[2])))
    ylim[1] > ylim[2] && (ylim = (ylim[2], ylim[1]))

    if labels == [""]
        labels = Vector{String}(undef, channel_n)
        for idx in 1:channel_n
            labels[idx] = ""
        end
    end
    
    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(floor(t[1], digits=2), ceil(t[end], digits=2)),
             xticks=floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2),
             ylims=ylim,
             title=title,
             palette=:darktest,
             titlefontsize=10,
             xlabelfontsize=8,
             ylabelfontsize=8,
             xtickfontsize=4,
             ytickfontsize=4;
             kwargs...)
    for idx in 1:channel_n
        p = plot!(t,
                  s_normalized[idx, 1:length(t)],
                  t=:line,
                  color=idx,
                  label=labels[idx])
    end

    return p
end

"""
    eeg_plot_butterfly(eeg; <keyword arguments>)

Butterfly plot of `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`: EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
- `offset::Int64=0`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
- `labels::Vector{String}=[""]`: channel labels vector
- `norm::Bool=false`: normalize the `signal` prior to calculations
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Amplitude [μV]"`: y-axis label
- `title::String=""`: plot title
- `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
- `head::Bool=true`: add head with electrodes
- `hist::Bool=true`: add histograms
- `average::Bool=true`: plot averaged signal
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_butterfly(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), head::Bool=true, hist::Bool=true, average::Bool=false, kwargs...)

    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))
    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    (epoch != 1 && (offset != 0 || len != 0)) && throw(ArgumentError("For epoch ≠ 1, offset and len must not be specified."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    (length(epoch) == 1 && (epoch < 1 || epoch > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    (length(epoch) > 1 && (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    if length(epoch) > 1
        sort!(epoch)
        (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
        len = eeg_epoch_len(eeg) * length(epoch)
        offset = eeg_epoch_len(eeg) * (epoch[1] - 1)
        epoch = 1
    end

    # default length is one epoch or 20 seconds
    if len == 0
        if eeg_epoch_len(eeg) > 20 * eeg_sr(eeg)
            len = 20 * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end

    # select channels, default is all channels
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    channel == 0 && (channel = 1:eeg_channel_n(eeg))
    length(channel) == 1 && throw(ArgumentError("number of channels must be ≥ 2."))
    length(channel) > 1 && sort!(channel)
    for idx in 1:length(channel)
        (channel[idx] < 1 || channel[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end

    # get epochs markers for len > epoch_len
    epoch_markers = Vector{Int64}[]
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_epoch_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> floor(Int64, offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= ceil(Int64, (offset + len) / eeg_sr(eeg))]
    else
        eeg_tmp = eeg
    end

    labels = eeg_labels(eeg_tmp)

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]
    t[1] = floor(t[1], digits=2)
    t[end] = ceil(t[end], digits=2)

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
   (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + length(t)), epoch]
    ndims(signal) == 1 && (signal = vec(signal))

    if average == false
        p = signal_plot_butterfly(t,
                                  signal,
                                  offset=offset,
                                  labels=labels,
                                  norm=norm,
                                  xlabel=xlabel,
                                  ylabel=ylabel,
                                  title=title,
                                  ylim=ylim;
                                  kwargs...)
    else
        p = signal_plot_avg(t,
                          signal,
                          offset=offset,
                          labels=labels,
                          norm=norm,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          ylim=ylim;
                          kwargs...)
    end

    # add epochs markers
    if norm == true
        s_normalized = signal_normalize_zscore(reshape(signal, eeg_channel_n(eeg), size(signal, 2), 1))
    else
        s_normalized = signal
    end
    ylim == (0, 0) && (ylim = (floor(minimum(s_normalized), digits=0), ceil(maximum(s_normalized), digits=0)))
    abs(ylim[1]) > abs(ylim[2]) && (ylim = (-abs(ylim[1]), abs(ylim[1])))
    abs(ylim[1]) < abs(ylim[2]) && (ylim = (-abs(ylim[2]), abs(ylim[2])))
    if length(epoch_markers) > 0 && len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                   linestyle=:dash,
                   linewidth=0.2,
                   linecolor=:black,
                   label="")
        for idx in 1:length(epoch_markers)
            p = plot!(annotation=((epoch_markers[idx] - 1), maximum(ceil.(abs.(signal))), text("E$(string(floor(Int64, epoch_markers[idx] / (eeg_epoch_len(eeg) / eeg_sr(eeg)))))", pointsize=6, halign=:center, valign=:top)))
        end
    end

    return p
end

"""
    signal_plot_psd(s_powers, s_freqs; <keyword arguments>)

Plot power spectrum density.

# Arguments

- `s_powers::Vector{Float64}`: signal powers
- `s_freqs::Vector{Float64}`: signal frequencies
- `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: x-axis limit
- `xlabel::String="Frequency [Hz]"`: x-axis label
- `ylabel::String="Power [μV^2/Hz]"`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_psd(s_powers::Vector{Float64}, s_freqs::Vector{Float64}; frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="", kwargs...)

    (frq_lim[1] < 0 || frq_lim[2] < 0) && throw(ArgumentError("frq_lim must be ≥ 0."))
    frq_lim == (0, 0) && (frq_lim = (0, s_freqs[end]))
    frq_lim[1] > frq_lim[2] && (frq_lim = (frq_lim[2], frq_lim[1]))

    p = plot(s_freqs,
             s_powers,
             xlabel=xlabel,
             ylabel=ylabel,
             xlims=frq_lim,
             legend=false,
             t=:line,
             c=:black,
             title=title,
             palette=:darktest,
             titlefontsize=10,
             xlabelfontsize=8,
             ylabelfontsize=8,
             xtickfontsize=4,
             ytickfontsize=4;
             kwargs...)

    plot(p)

    return p
end

"""
    signal_plot_psd(signal; <keyword arguments>)

Plot `signal` channel power spectrum density.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64`: sampling frequency
- `norm::Bool=false`: converts power to dB
- `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: x-axis limit
- `xlabel::String="Frequency [Hz]"`: x-axis label
- `ylabel::String="Power [μV^2/Hz]"`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_psd(signal::Vector{Float64}; fs::Int64, norm::Bool=false, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="", kwargs...)

    fs < 0 && throw(ArgumentError("fs must be ≥ 0."))
    (frq_lim[1] < 0 || frq_lim[2]) < 0 && throw(ArgumentError("frq_lim must be ≥ 0."))
    frq_lim == (0, 0) && (frq_lim = (0, fs / 2))
    frq_lim[1] > frq_lim[2] && (frq_lim = (frq_lim[2], frq_lim[1]))

    s_freqs, s_powers = signal_psd(signal, fs=fs, norm=norm)
    norm == true && (ylabel::String="Power [dB]")

    p = signal_plot_psd(s_powers,
                        s_freqs,
                        title=title,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        xlims=frq_lim;
                        kwargs...)

    plot(p)

    return p
end

"""
    signal_plot_psd(signal; <keyword arguments>)

Plot `signal` channels power spectrum density.

# Arguments

- `signal::Matrix{Float64}`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: power in dB
- `average::Bool=false`: plots average power and 95%CI for all channels
- `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: x-axis limit
- `labels::Vector{String}=[""]`: channel labels vector
- `xlabel::String="Frequency [Hz]"`: x-axis label
- `ylabel::String="Power [μV^2/Hz]"`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false, average::Bool=false, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), labels::Vector{String}=[""], xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="", kwargs...)

    fs < 0 && throw(ArgumentError("fs must be ≥ 0."))
    (frq_lim[1] < 0 || frq_lim[2]) < 0 && throw(ArgumentError("frq_lim must be ≥ 0."))
    frq_lim[1] > frq_lim[2] && (frq_lim = (frq_lim[2], frq_lim[1]))

    norm == true && (ylabel="Power [dB]")

    channel_n = size(signal, 1)
    signal = reshape(signal, size(signal, 1), size(signal, 2), 1)
    s_powers, s_freqs = signal_psd(signal, fs=fs, norm=norm)
    s_powers = s_powers[:, :, 1]
    s_freqs = s_freqs[:, :, 1]
    frq_lim == (0, 0) && (frq_lim = (0, s_freqs[1, end]))

    if average == true
        s_powers_m, s_powers_s, s_powers_u, s_powers_l = signal_ci95(s_powers)
        s_freqs = s_freqs[1, :]
        channel_n = 1
        labels == [""]
    end

    if labels == [""]
        labels = Vector{String}(undef, channel_n)
        for idx in 1:channel_n
            labels[idx] = ""
        end
    end

    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=frq_lim,
             title=title,
             palette=:darktest,
             titlefontsize=10,
             xlabelfontsize=8,
             ylabelfontsize=8,
             xtickfontsize=4,
             ytickfontsize=4;
             kwargs...)
    if average == true
        p = plot!(s_freqs,
                  s_powers_u,
                  fillrange=s_powers_l,
                  fillalpha=0.35,
                  label=false,
                  t=:line,
                  c=:grey,
                  lw=0.5)
        p = plot!(s_freqs,
                  s_powers_l,
                  label=false,
                  t=:line,
                  c=:grey,
                  lw=0.5)
        p = plot!(s_freqs,
                  s_powers_m,
                  label=false,
                  t=:line,
                  c=:black)
    else
        for idx in 1:channel_n
            p = plot!(s_freqs[idx, :],
                      s_powers[idx, :],
                      label=labels[idx],
                      t=:line)
        end
    end

    return p
end

"""
    eeg_plot_psd(eeg; <keyword arguments>)

Plot `eeg` channels power spectrum density.

# Arguments

- `eeg::NeuroJ.EEG`: EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
- `offset::Int64=0`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
- `labels::Vector{String}=[""]`: channel labels vector
- `norm::Bool=false`: power in dB
- `average::Bool=false`: plots average power and 95%CI for all channels
- `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: x-axis limit
- `xlabel::String="Frequency [Hz]`: x-axis label
- `ylabel::String="Power [μV^2/Hz]"`: y-axis label
- `title::String=""`: plot title
- `head::Bool=false`: add head with electrodes
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_psd(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], norm::Bool=false, average::Bool=false, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="", head::Bool=false, kwargs...)

    (frq_lim[1] < 0 || frq_lim[1] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    (frq_lim[2] < 0 || frq_lim[2] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))
    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    (epoch != 1 && (offset != 0 || len != 0)) && throw(ArgumentError("For epoch ≠ 1, offset and len must not be specified."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    (length(epoch) == 1 && (epoch < 1 || epoch > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    (length(epoch) > 1 && (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    if length(epoch) > 1
        sort!(epoch)
        (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
        len = eeg_epoch_len(eeg) * length(epoch)
        offset = eeg_epoch_len(eeg) * (epoch[1] - 1)
        epoch = 1
    end

    # default length is one epoch or 20 seconds
    if len == 0
        if eeg_epoch_len(eeg) > 20 * eeg_sr(eeg)
            len = 20 * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end

    # select channels, default is 1:20 or all channels
    if channel == 0
        if eeg_channel_n(eeg) >= 20
            channel = 1:20
        else
            channel = 1:eeg_channel_n(eeg)
        end
    end
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    length(channel) > 1 && sort!(channel)
    for idx in 1:length(channel)
        (channel[idx] < 1 || channel[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
    (length(channel) == 1 && average == true) && throw(ArgumentError("channel must contain ≥ 2 channels if average=true"))

    # get epochs markers for len > epoch_len
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
    else
        eeg_tmp = eeg
    end

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
    (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    fs = eeg_sr(eeg_tmp)
    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + len), epoch]
    ndims(signal) == 1 && (signal = vec(signal))

    labels = eeg_labels(eeg_tmp)

    p = signal_plot_psd(signal,
                        fs=fs,
                        labels=labels,
                        norm=norm,
                        average=average,
                        frq_lim=frq_lim,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        title=title;
                        kwargs...)

    return p
end

"""
    eeg_plot_electrodes(eeg; <keyword arguments>)

Plot `eeg` electrodes.

# Arguments

- `eeg:EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: which channel should be highlighted, default is all channels
- `labels::Bool=true`: plot electrode labels
- `head::Bool`=true: plot head
- `head_labels::Bool=false`: plot head labels
- `small::Bool=false`: draws small plot
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_electrodes(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head::Bool=true, head_labels::Bool=false, small::Bool=false, kwargs...)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() first."))
    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    # select channels, default is all channels
    channel == 0 && (channel = 1:eeg_channel_n(eeg))
    length(channel) > 1 && sort!(channel)
    for idx in 1:length(channel)
        (channel[idx] < 1 || channel[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
    typeof(channel) <: AbstractRange && (channel = collect(channel))

    # select channels, default is all channels
    if selected != 0
        length(selected) > 1 && sort!(selected)
        for idx in 1:length(selected)
            (selected[idx] < 1 || selected[idx] > eeg.eeg_header[:channel_n]) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
        end
        typeof(selected) <: AbstractRange && (selected = collect(selected))
        length(selected) > 1 && (intersect(selected, channel) == selected || throw(ArgumentError("channel must include selected.")))
        length(selected) == 1 && (intersect(selected, channel) == [selected] || throw(ArgumentError("channel must include selected.")))
    end
    
    eeg_tmp = eeg_keep_channel(eeg, channel=channel)

    # look for location data
    loc_x = zeros(eeg_channel_n(eeg_tmp, type=:eeg))
    loc_y = zeros(eeg_channel_n(eeg_tmp, type=:eeg))
    for idx in 1:eeg_channel_n(eeg_tmp, type=:eeg)
        loc_y[idx], loc_x[idx] = pol2cart(pi / 180 * eeg_tmp.eeg_header[:loc_theta][idx],
                                          eeg_tmp.eeg_header[:loc_radius][idx])
    end
    x_lim = (findmin(loc_x)[1] * 1.8, findmax(loc_x)[1] * 1.8)
    y_lim = (findmin(loc_y)[1] * 1.8, findmax(loc_y)[1] * 1.8)

    if small == true
        plot_size = (800, 800)
        marker_size = 4
        labels = false
    else
        plot_size = (800, 800)
        marker_size = 6
        font_size = 4
    end

    p = plot(grid=true,
             framestyle=:none,
             palette=:darktest,
             size=plot_size,
             markerstrokewidth=0,
             border=:none,
             aspect_ratio=1,
             margins=-20Plots.px,
             titlefontsize=10;
             kwargs...)
    if length(selected) == eeg_tmp.eeg_header[:channel_n]
        for idx in 1:eeg_tmp.eeg_header[:channel_n]
            p = plot!((loc_x[idx], loc_y[idx]),
                      color=idx,
                      seriestype=:scatter,
                      xlims=x_lim,
                      ylims=x_lim,
                      grid=true,
                      label="",
                      markersize=marker_size,
                      markerstrokewidth=0,
                      markerstrokealpha=0;
                      kwargs...)
        end
    else
        p = plot!(loc_x,
                  loc_y,
                  seriestype=:scatter,
                  color=:black,
                  alpha=0.2,
                  xlims=x_lim,
                  ylims=y_lim,
                  grid=true,
                  label="",
                  markersize=marker_size,
                  markerstrokewidth=0,
                  markerstrokealpha=0;
                  kwargs...)
        if selected != 0
            eeg_tmp = eeg_keep_channel(eeg, channel=selected)
            loc_x = zeros(eeg_channel_n(eeg_tmp, type=:eeg))
            loc_y = zeros(eeg_channel_n(eeg_tmp, type=:eeg))
            for idx in 1:eeg_channel_n(eeg_tmp, type=:eeg)
                loc_y[idx], loc_x[idx] = pol2cart(pi / 180 * eeg_tmp.eeg_header[:loc_theta][idx],
                                                  eeg_tmp.eeg_header[:loc_radius][idx])
            end
            for idx in 1:eeg_tmp.eeg_header[:channel_n]
                p = plot!((loc_x[idx], loc_y[idx]), color=idx, seriestype=:scatter, xlims=x_lim, ylims=x_lim, grid=true, label="", markersize=marker_size, markerstrokewidth=0, markerstrokealpha=0)
            end
        end
    end
    if labels == true
        for idx in 1:length(eeg_labels(eeg_tmp))
            plot!(annotation=(loc_x[idx], loc_y[idx] + 0.05, text(eeg_labels(eeg_tmp)[idx], pointsize=font_size)))
        end
        p = plot!()
    end
    if head == true
        # for some reason head is enlarged for channel > 1
        eeg_tmp = eeg_keep_channel(eeg, channel=1)
        loc_x = zeros(eeg_channel_n(eeg_tmp, type=:eeg))
        loc_y = zeros(eeg_channel_n(eeg_tmp, type=:eeg))
        for idx in 1:eeg_channel_n(eeg_tmp, type=:eeg)
            loc_y[idx], loc_x[idx] = pol2cart(pi / 180 * eeg_tmp.eeg_header[:loc_theta][idx],
                                              eeg_tmp.eeg_header[:loc_radius][idx])
        end
        hd = eeg_draw_head(p, loc_x, loc_x, head_labels=head_labels)
        plot!(hd)
    end

    plot(p)

    return p
end

"""
    eeg_plot_matrix(eeg, m; <keyword arguments>)

Plot matrix `m` of `eeg` channels.

# Arguments

- `eeg:EEG`
- `m::Union{Matrix{Float64}, Array{Float64, 3}}`: channels by channels matrix
- `epoch::Int64=1`: epoch number to display
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_matrix(eeg::NeuroJ.EEG, m::Union{Matrix{Float64}, Array{Float64, 3}}; epoch::Int64=1, kwargs...)
    (epoch < 1 || epoch > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    labels = eeg_labels(eeg)
    channel_n = size(m, 1)
    ndims(m) == 3 && (m = m[:, :, epoch])

    p = heatmap(m,
                xticks=(1:channel_n, labels),
                yticks=(1:channel_n, eeg_labels(eeg)),
                titlefontsize=10,
                xlabelfontsize=8,
                ylabelfontsize=8,
                xtickfontsize=4,
                ytickfontsize=4;
                kwargs...)
    plot(p)

    return p
end

"""
    eeg_plot_matrix(eeg, cov_m, lags; <keyword arguments>)

Plot covariance matrix `m` of `eeg` channels.

# Arguments

- `eeg:EEG`
- `cov_m::Union{Matrix{Float64}, Array{Float64, 3}}`: covariance matrix
- `lags::Union{Vector{Int64}, Vector{Float64}}`: covariance lags
- `channel::Union{Int64, Vector{Int64}, UnitRange{Int64}, Nothing}`: channel to display
- `epoch::Int64=1`: epoch number to display
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_covmatrix(eeg::NeuroJ.EEG, cov_m::Union{Matrix{Float64}, Array{Float64, 3}}, lags::Union{Vector{Int64}, Vector{Float64}}; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch::Int64=1, kwargs...)

    (epoch < 1 || epoch > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    # select channels, default is 1:20 or all channels
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    channel == 0 && (channel = 1:eeg_channel_n(eeg))
    length(channel) > 1 && sort!(channel)
    for idx in 1:length(channel)
        (channel[idx] < 1 || channel[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end

    labels = eeg_labels(eeg)
    ndims(cov_m) == 3 && (cov_m = cov_m[:, :, epoch])
    p = []
    for idx in channel
        push!(p, plot(lags, cov_m[idx, :], title="ch: $(labels[idx])", label="", titlefontsize=6, xlabelfontsize=8, ylabelfontsize=8, xtickfontsize=4, ytickfontsize=4, lw=0.5))
    end
    p = plot(p...; kwargs...)

    return p
end

"""
    signal_plot_spectrogram(signal; <keyword arguments>)

Plot spectrogram of `signal`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64`: sampling frequency
- `offset::Int64=0`: displayed segment offset in samples
- `norm::Bool=true`: normalize powers to dB
- `demean::Bool=true`: demean signal prior to analysis
- `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Frequency [Hz]"`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_spectrogram(signal::Vector{Float64}; fs::Int64, offset::Int64=0, norm::Bool=true, demean::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel="Time [s]", ylabel="Frequency [Hz]", title="", kwargs...)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1 Hz."))
    (frq_lim[1] < 0 || frq_lim[1] > fs / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(fs / 2)."))
    (frq_lim[2] < 0 || frq_lim[2] > fs / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(fs / 2)."))
    frq_lim = tuple_order(frq_lim)
    frq_lim == (0, 0) && (frq_lim = (0, fs / 2))

    demean == true && (signal = signal_demean(signal))

    nfft = length(signal)
    interval = fs
    overlap = round(Int64, fs * 0.85)

    frq_lim == (0, 0) && (frq_lim = (0, fs/2))
    spec = spectrogram(signal, interval, overlap, nfft=nfft, fs=fs, window=hanning)
    t = collect(spec.time) .+ (offset / fs)

    if norm == false
        p = heatmap(t,
                    spec.freq,
                    spec.power,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    ylims=frq_lim,
                    xticks=floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2),
                    title=title,
                    colorbar_title="Power/frequency [μV^2/Hz]",
                    titlefontsize=10,
                    xlabelfontsize=8,
                    ylabelfontsize=8,
                    xtickfontsize=4,
                    ytickfontsize=4;
                    kwargs...)
    else
        # in dB
        p = heatmap(t,
                    spec.freq,
                    pow2db.(spec.power),
                    xlabel=xlabel,
                    ylabel=ylabel,
                    ylims=frq_lim,
                    xticks=floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2),
                    title=title,
                    colorbar_title="Power/frequency [dB/Hz]",
                    titlefontsize=10,
                    xlabelfontsize=8,
                    ylabelfontsize=8,
                    xtickfontsize=4,
                    ytickfontsize=4;
                    kwargs...)
    end

    return p
end

"""
    eeg_plot_spectrogram(eeg; <keyword arguments>)

Plots spectrogram of `eeg` channel.

# Arguments

- `eeg:EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epoch to plot
- `channel::Int64`: channel to plot
- `offset::Int64=0`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
- `norm::Bool=true`: normalize powers to dB
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Frequency [Hz]"`: y-axis label
- `title::String=""`: plot title
- `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_spectrogram(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Int64, offset::Int64=0, len::Int64=0, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel::String="Time [s]", ylabel::String="Frequency [Hz]", title::String="", kwargs...)

    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))
    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    (frq_lim[1] < 0 || frq_lim[1] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    (frq_lim[2] < 0 || frq_lim[2] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    frq_lim == (0, 0) && (frq_lim = (0, eeg_sr(eeg) / 2))
    frq_lim = tuple_order(frq_lim)
    fs = eeg_sr(eeg)

    (epoch != 1 && (offset != 0 || len != 0)) && throw(ArgumentError("For epoch ≠ 1, offset and len must not be specified."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    (length(epoch) == 1 && (epoch < 1 || epoch > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    (length(epoch) > 1 && (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    if length(epoch) > 1
        sort!(epoch)
        (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
        len = eeg_epoch_len(eeg) * length(epoch)
        offset = eeg_epoch_len(eeg) * (epoch[1] - 1)
        epoch = 1
    end

    # default length is one epoch or 20 seconds
    if len == 0
        if eeg_epoch_len(eeg) > 20 * eeg_sr(eeg)
            len = 20 * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end

    # get epochs markers for len > epoch_len
    epoch_markers = Vector{Int64}[]
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_epoch_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> floor(Int64, offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= ceil(Int64, (offset + len) / eeg_sr(eeg))]
    else
        eeg_tmp = eeg
    end

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
    (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + len), epoch]

    p = signal_plot_spectrogram(signal,
                                fs=fs,
                                offset=offset,
                                norm=norm,
                                xlabel=xlabel,
                                ylabel=ylabel,
                                frq_lim=frq_lim,
                                title=title;
                                kwargs...)

    # add epochs markers
    if length(epoch_markers) > 0 && len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                   linestyle=:dash,
                   linewidth=0.2,
                   linecolor=:black,
                   label="")
        for idx in 1:length(epoch_markers)
            p = plot!(annotation=((epoch_markers[idx] - 1), frq_lim[2], text("E$(string(floor(Int64, epoch_markers[idx] / (eeg_epoch_len(eeg) / eeg_sr(eeg)))))", pointsize=4, halign=:center, valign=:top)))
        end
    end

    plot(p)

    return p
end

"""
    signal_plot_histogram(signal; <keyword arguments>)

Plot histogram of `signal`.

# Arguments

- `signal::Vector{Float64}`
- `type::Symbol`: type of histogram: regular `:hist` or kernel density `:kd`
- `label::String=""`: channel label
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_histogram(signal::Vector{Float64}; type::Symbol=:hist, label::String="", xlabel::String="", ylabel::String="", title::String="", kwargs...)

    type === :kd && (type = :density)

    p = plot(signal,
             seriestype=type,
             xlabel=xlabel,
             ylabel=ylabel,
             label=label,
             title=title,
             palette=:darktest,
             grid=false,
             linecolor=1,
             fillcolor=1,
             linewidth=0.5,
             margins=0Plots.px,
             yticks=false,
             titlefontsize=10,
             xlabelfontsize=8,
             ylabelfontsize=8,
             xtickfontsize=4,
             ytickfontsize=4,
             xticks=[floor(minimum(signal), digits=1), 0, ceil(maximum(signal), digits=1)];
             kwargs...)

    plot(p)

    return p
end


"""
    signal_plot_histogram(signal; <keyword arguments>)

Plot histogram of `signal`.

# Arguments

- `signal::Matrix{Float64}`
- `type::Symbol`: type of histogram: :hist or :kd
- `labels::Vector{String}=[""]`
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_histogram(signal::Union{Vector{Float64}, Matrix{Float64}}; type::Symbol=:hist, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", kwargs...)

    channel_n = size(signal, 1)

    # reverse so 1st channel is on top
    signal = reverse(signal, dims = 1)

    labels == [""] && (labels = repeat([""], channel_n))

    # plot channels
    if channel_n == 1
            p = signal_plot_histogram(signal,
                                      type=type,
                                      labels="",
                                      xlabel=xlabel,
                                      ylabel=ylabel,
                                      title=title,
                                      legend=false,
                                      linecolor=1,
                                      fillcolor=1,
                                      linewidth=0.5,
                                      margins=0Plots.px,
                                      yticks=true,
                                      titlefontsize=10,
                                      xlabelfontsize=8,
                                      ylabelfontsize=8,
                                      xtickfontsize=4,
                                      ytickfontsize=4,
                                      xticks=[floor(minimum(signal), digits=1), 0, ceil(maximum(signal), digits=1)];
                                      kwargs...
)    else
        p = []
        for idx in 1:channel_n
            push!(p, signal_plot_histogram(signal[idx, :],
                                           type=type,
                                           label="",
                                           xlabel=xlabel,
                                           ylabel=labels[idx],
                                           title=title,
                                           linecolor=idx,
                                           fillcolor=idx,
                                           linewidth=0.5,
                                           left_margin=30Plots.px,
                                           yticks=true,
                                           titlefontsize=10,
                                           xlabelfontsize=8,
                                           ylabelfontsize=8,
                                           xtickfontsize=4,
                                           ytickfontsize=4,
                                           xticks=[floor(minimum(signal), digits=1), 0, ceil(maximum(signal), digits=1)];
                                           kwargs...))
        end
        p = plot(p...,
                 palette=:darktest,
                 layout=(channel_n, 1);
                 kwargs...)
    end

    plot(p)

    return p
end

"""
    eeg_plot_histogram(eeg; <keyword arguments>)

Plot `eeg` channel histograms.

# Arguments

- `eeg::NeuroJ.EEG`: EEG object
- `type::Symbol: type of histogram: :hist or :kd
- `epoch::Int64=1`: epoch number to display
- `channel::Int64`: channel to display
- `offset::Int64=0`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default 1 epoch or 20 seconds
- `label::String=""`: channel label
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_histogram(eeg::NeuroJ.EEG; type::Symbol=:hist, epoch::Int64=1, channel::Int64, offset::Int64=0, len::Int64=0, label::String="", xlabel::String="", ylabel::String="", title::String="", kwargs...)

    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))
    (epoch < 1 || epoch > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    (channel < 1 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))

    # default length is one epoch or 20 seconds
    if len == 0
        if eeg_epoch_len(eeg) > 20 * eeg_sr(eeg)
            len = 20 * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end

    # get epochs markers for len > epoch_len
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
    else
        eeg_tmp = eeg
    end

    label == "" && (label = eeg_labels(eeg_tmp)[channel])

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
    (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    signal = vec(eeg_tmp.eeg_signals[channel, (1 + offset):(offset + len), epoch])

    p = signal_plot_histogram(signal,
                              type=type,
                              labels=label,
                              xlabel=xlabel,
                              ylabel=ylabel,
                              title=title;
                              kwargs...)

    return p
end

"""
    signal_plot_ica(t, ica; <keyword arguments>)

Plot `ica` components against time vector `t`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`: the time vector
- `ica::Vector{Float64}`
- `label::String=""`: channel label
- `norm::Bool=true`: normalize the `ica` prior to calculations
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Amplitude [μV]"`: y-axis label
- `title::String=""`: plot title
- `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits (-ylim:ylim)
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_ica(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, ica::Vector{Float64}; label::String="", norm::Bool=true, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    if ylim === nothing
        ylim = maximum(ica) * 1.5
        ylim = ceil(Int64, ylim)
    end

    p = plot(t,
             ica[1:length(t)],
             xlabel=xlabel,
             ylabel=ylabel,
             label="",
             xlims=(floor(t[1], digits=2), ceil(t[end], digits=2)),
             xticks=floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2),
             ylims=(-ylim, ylim),
             title=title,
             titlefontsize=10,
             xlabelfontsize=8,
             ylabelfontsize=8,
             xtickfontsize=4,
             ytickfontsize=4,
             palette=:darktest;
             kwargs...)

    plot(p)

    return p
end

"""
    signal_plot_ica(t, ica; <keyword arguments>)

Plots `ica` components.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
- `ica::Matrix{Float64}`
- `labels::Vector{String}=[""]`: labels vector
- `norm::Bool=true`: normalize the ICs prior to calculations
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_ica(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange,}, ica::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    ica_n = size(ica, 1)

    # reverse so 1st channel is on top
    ica_color = ica_n:-1:1
    ica = reverse(ica[:, :], dims = 1)
    ica_normalized = zeros(size(ica))

    if labels == [""]
        labels = Vector{String}(undef, size(ica, 1))
        for idx in 1:size(ica, 1)
            labels[idx] = "IC $idx"
        end
    end

    if norm == true
        # normalize and shift so all channels are visible
        variances = var(ica, dims=2)
        mean_variance = 10 * mean(variances)
        for idx in 1:ica_n
            i = @view ica[idx, :]
            ica_normalized[idx, :] = (i .- mean(i)) ./ mean_variance .+ (idx - 1)
        end
    else
        ica_normalized = ica
    end

    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(floor(t[1], digits=2), ceil(t[end], digits=2)),
             xticks=floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2),
             ylims=(-0.5, ica_n-0.5),
             title=title,
             titlefontsize=10,
             xlabelfontsize=8,
             ylabelfontsize=8,
             xtickfontsize=4,
             ytickfontsize=4,
             palette=:darktest;
             kwargs...)
    for idx in 1:ica_n
        p = plot!(t,
                  ica_normalized[idx, 1:length(t)],
                  label="",
                  color=ica_color[idx])
    end
    p = plot!(p, yticks=((ica_n - 1):-1:0, labels))

    return p
end

"""
    eeg_plot_ica(eeg; <keyword arguments>)

Plots embedded ICs components.

# Arguments

- `eeg::NeuroJ.EEG`
- `epoch::Int64=1`: epoch number to display
- `offset::Int64=0`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default 1 epoch or 20 seconds
- `ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing`: which IC to plot, default all
- `norm::Bool=true`: normalize the ICs prior to calculations
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_ica(eeg::NeuroJ.EEG; epoch::Int64=1, offset::Int64=0, len::Int64=0, ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    :ica in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain ICA. Perform eeg_ica!(EEG) first."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))

    len > eeg_signal_len(eeg) && throw(ArgumentError("len must be < $(eeg_signal_len(eeg))."))
    (epoch < 1 || epoch > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))

    # default length is one epoch or 20 seconds
    if len == 0
        if eeg_epoch_len(eeg) > 20 * eeg_sr(eeg)
            len = 20 * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end

    # select ICs, default is all
    ica_idx = findfirst(isequal(:ica), eeg.eeg_header[:components])
    ic === nothing && (ic = 1:size(eeg.eeg_components[ica_idx], 1))
    typeof(ic) <: AbstractRange && (ic = collect(ic))
    if typeof(ic) == Vector{Int64}
        sort(ic)
        for idx in 1:length(ic)
            ic[idx] > size(eeg.eeg_components[ica_idx], 1) && throw(ArgumentError("ic must be ≥ 1 and ≤ $(size(eeg.eeg_components[ica_idx], 1))."))
        end
    else
        (ic <= 0 || ic > size(eeg.eeg_components[ica_idx], 1)) && throw(ArgumentError("ic must be > 0 and ≤ $(size(eeg.eeg_components[ica_idx], 1))."))
    end
    (typeof(ic) == Int64 && ic > size(eeg.eeg_components[ica_idx], 1)) && throw(ArgumentError("ic must be ≥ 1 and ≤ $(size(eeg.eeg_components[ica_idx], 1))."))

    # get epochs markers for len > epoch_len
    epoch_markers = Vector{Int64}[]   
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_epoch_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        eeg_h = eeg_history(eeg)
        ica_par = ""
        for idx in 1:length(eeg_h)
            occursin("eeg_ica!", eeg_h[idx]) && (ica_par = eeg_h[idx])
        end
        ica_n = r"n\=([0-9]+)"
        m = match(ica_n, ica_par)
        ica_n = replace(m.match, "n=" => "")
        ica_n = replace(ica_n, "))" => "")
        ica_n = parse(Int64, ica_n)
        ica_tol = r"tol\=([0-9].[0-9]+)"
        m = match(ica_tol, ica_par)
        ica_tol = replace(m.match, "tol=" => "")
        ica_tol = replace(ica_tol, "))" => "")
        ica_tol = parse(Float64, ica_tol)
        ica_iter = r"iter\=([0-9]+)"
        m = match(ica_iter, ica_par)
        ica_iter = replace(m.match, "iter=" => "")
        ica_iter = replace(ica_iter, "))" => "")
        ica_iter = parse(Int64, ica_iter)
        ica_f = r"f\=.+"
        m = match(ica_f, ica_par)
        ica_f = replace(m.match, "f=" => "")
        ica_f = replace(ica_f, "))" => "")
        ica_f = Symbol(ica_f)
        eeg_ica!(eeg_tmp, n=ica_n, iter=ica_iter, tol=ica_tol, f=ica_f)
    else
        eeg_tmp = eeg
    end

    labels = Vector{String}(undef, size(eeg_tmp.eeg_components[ica_idx], 1))
    for idx in 1:size(eeg_tmp.eeg_components[ica_idx], 1)
        labels[idx] = "IC $idx"
    end
    labels = labels[ic]

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]
    t[1] = floor(t[1], digits=2)
    t[end] = ceil(t[end], digits=2)

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
    (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    ica_m = eeg_tmp.eeg_components[ica_idx][ic, (1 + offset):(offset + Int(len)), epoch]

    p = signal_plot_ica(t,
                        ica_m,
                        norm=norm,
                        label=labels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        title=(title * "\nIC #$(ic)"),
                        kwargs...)
    # add epochs markers
    if length(epoch_markers) > 0 && len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                  linestyle=:dash,
                  linewidth=0.5,
                  linecolor=:black,
                  label="")
        for idx in 1:(floor(Int64, (offset + len) / eeg_epoch_len(eeg)))
            p = plot!(annotation=((epoch_markers[idx] - 1), ((length(ic) - 1) * 1.02), text("E$idx", pointsize=6, halign=:center, valign=:center)))
        end
    end

    if length(ic) == 1 || typeof(ic) == Int64
        p = plot!(p, title="IC #$(lpad(string(ic), 3, "0"))")
        ica_weights_idx = findfirst(isequal(:ica_mw), eeg.eeg_header[:components])
        ica_weights = round.(eeg.eeg_components[ica_weights_idx][:, ic, epoch], digits=1)
        eeg_tmp = deepcopy(eeg)
        eeg_tmp.eeg_header[:labels] = string.(ica_weights)
        hd = eeg_plot_electrodes(eeg_tmp, labels=true, selected=0, head=true, title="Weights", alpha=1)
        p = plot!(p, hd, layout=(2, 1))
    end

    plot(p)

    return p
end

"""
    eeg_plot_topo(eeg; <keyword arguments>)

Plot topographical view of `eeg` component.

# Arguments

- `eeg::NeuroJ.EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epochs to display
- `offset::Int64=1`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default is 1 second
- `m::Symbol=:shepard`: interpolation method :shepard (Shepard), :mq (Multiquadratic), :tp (ThinPlate)
- `c::Symbol=:amp`: component name (:ica, :pca, :amp, :power)
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange, Tuple, Nothing}=nothing`: component index, e.g. ICA number or frequency range
- `norm::Bool=true`: convert power as dB
- `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: frequency limit for PSD and spectrogram
- `head_labels::Bool=false`: plot head labels
- `cb::Bool=false`: add color bars to plots
- `cb_label::String=""`: color bar label
- `average::Bool=true`: plot averaged signal and PSD
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_topo(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, offset::Int64=0, len::Int64=0, m::Symbol=:shepard, c::Symbol=:amp, c_idx::Union{Int64, Vector{Int64}, AbstractRange, Tuple, Nothing}=nothing, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0,0), head_labels::Bool=false, cb::Bool=false, cb_label::String="", average::Bool=true, title::String="", kwargs...)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))

    m in [:shepard, :mq, :tp] || throw(ArgumentError("m must be :shepard, :mq or :tp."))
    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() first."))
    (c === :amp  || c === :power || c in eeg.eeg_header[:components]) || throw(ArgumentError("Component $c not found."))
    frq_lim = tuple_order(frq_lim)

    len < 0 && throw(ArgumentError("len must be > 0."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    (epoch < 1 || epoch > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    # default length is 100 ms
    len == 0 && (len = floor(Int64, eeg_sr(eeg) / 10))
    # offset == 0 && (offset = (epoch - 1) * eeg_epoch_len(eeg) + 1)
    epoch == 0 && (epoch = floor(Int64, offset / eeg_epoch_len(eeg)) + 1)
    offset + len > epoch * eeg_epoch_len(eeg) && throw(ArgumentError("offset + len must be ≤ $(epoch * eeg_epoch_len(eeg))."))
    # offset < ((epoch - 1) * eeg_epoch_len(eeg) + 1 + offset) && throw(ArgumentError("offset must be ≥ $((epoch - 1) * eeg_epoch_len(eeg) + 1)."))

#=
    if offset > (eeg_epoch_len(eeg) - len) && eeg_epoch_n(eeg) > 1
        epoch = 0
        while offset > eeg_epoch_len(eeg)
            epoch += 1
            offset -= eeg_epoch_len(eeg)
        end
    else
        epoch = 1
    end
=#

    typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
    # ignore c_idx for components other than ICA
    if typeof(c_idx) == Vector{Int64} && c === :ica
        component_idx = findfirst(isequal(c), eeg.eeg_header[:components])
        for idx in length(c_idx):-1:1
            c_idx[idx] > size(eeg.eeg_components[component_idx], 1) && throw(ArgumentError("For component $c range must be 1:$(size(eeg.eeg_components[component_idx], 1))"))
        end
        p = []
        length(c_idx) <= 10 && (l_row = 2)
        length(c_idx) > 10 && (l_row = 4)
        l = (l_row, ceil(Int64, length(c_idx) / l_row))
        for idx in 1:length(c_idx)
            push!(p, eeg_plot_topo(eeg; offset=offset, len=len, m=m, c=c, c_idx=c_idx[idx], head_labels=head_labels, title="", colorbar=cb))
        end
        for idx in (length(c_idx) + 1):l[1]*l[2]
            push!(p, plot(border=:none, title=""))
        end
        p = plot!(p..., layout=l, size=(l[2] * 800, l[1] * 800))
        plot(p)

        if figure !== ""
            isfile(figure) && @warn "File $figure will be overwritten."
            savefig(p, figure)
        end

        return p
    end

    t_1 = floor(eeg.eeg_time[1 + offset], digits=2)
    t_2 = ceil(eeg.eeg_time[1 + offset + len], digits=2)
    t_1 < 1.0 && (t_s1 = string(floor(t_1 * 1000, digits=2)) * " ms")
    t_1 >= 1.0 && (t_s1 = string(floor(t_1, digits=2)) * " s")
    t_2 < 1.0 && (t_s2 = string(ceil(t_2 * 1000, digits=2)) * " ms")
    t_2 >= 1.0 && (t_s2 = string(ceil(t_2, digits=2)) * " s")

    if c === :amp
        s_non_interpolated = mean(eeg.eeg_signals[:, 1 + offset:(offset + len), epoch], dims=2)
        title = "Unweighted amplitude"
    elseif c === :ica
        s = eeg.eeg_signals[:, :, epoch]
        s = reshape(s, size(s, 1), size(s, 2), 1)
        component_idx = findfirst(isequal(c), eeg.eeg_header[:components])
        a = eeg.eeg_components[component_idx][c_idx, :, epoch]
        component_idx = findfirst(isequal(:ica_mw), eeg.eeg_header[:components])
        m_v = eeg.eeg_components[component_idx][:, c_idx, epoch]
        s_w = m_v * a'
        s_non_interpolated = mean(s_w[:, 1 + offset:(offset + len), epoch], dims=2)
        title = "$(uppercase(string(c))) #$(lpad(string(c_idx), 3, "0"))"
    elseif c === :pca
        s = eeg.eeg_signals[:, :, epoch]
        s = reshape(s, size(s, 1), size(s, 2), 1)
        pca_idx = findfirst(isequal(:pca), eeg.eeg_header[:components])
        pca = eeg.eeg_components[pca_idx][:, :, epoch]
        pca_m_idx = findfirst(isequal(:pca_m), eeg.eeg_header[:components])
        pca_m = eeg.eeg_components[pca_m_idx]
        s_reconstructed = reconstruct(pca_m, pca)
        s_non_interpolated = mean(s_reconstructed[:, 1 + offset:(offset + len), epoch], dims=2)
        title = "$(uppercase(string(c))) #1:$(size(pca, 1)) reconstruction [A.U.]"
    elseif c === :power
        if typeof(c_idx) <: Tuple
            s = eeg.eeg_signals[:, 1 + 1 + offset:(offset + len), epoch]
            s_non_interpolated = zeros(size(s, 1))
            for idx in 1:size(s_non_interpolated, 1)
                s_non_interpolated[idx] = signal_band_power(s[idx, :], fs=eeg_sr(eeg), f=c_idx)
            end
            norm == true && s_non_interpolated[s_non_interpolated .<= 0] .= eps()
            norm == true && (s_non_interpolated = pow2db.(s_non_interpolated))
            title = "Power [A.U.]\n[frequency range: $(c_idx[1]):$(c_idx[end]) Hz]"
        elseif typeof(c_idx) == Int64 || typeof(c_idx) == Float64
            s = eeg.eeg_signals[:, 1 + offset:(offset + len), epoch]
            s = reshape(s, size(s, 1), size(s, 2), 1)
            s_psd, s_frq = signal_psd(s, fs=eeg_sr(eeg), norm=norm)
            # _, _, s_p, _ = signal_spectrum(s)
            # norm == true && (s_p = pow2db.(s_p[:, :, 1]))
            # s_f = freqs(s, eeg_sr(eeg))
            frq_idx = vsearch(c_idx, s_frq[1, :])
            s_non_interpolated = zeros(size(s, 1))
            for idx in 1:size(s_non_interpolated, 1)
                s_non_interpolated[idx] = s_psd[idx, frq_idx]
            end
            title = "Power [A.U.]\n[frequency: $c_idx Hz]"
        end
    else
        throw(ArgumentError("Component $c not found."))
    end

    # plot signal at electrodes at time
    loc_x = zeros(eeg_channel_n(eeg))
    loc_y = zeros(eeg_channel_n(eeg))
    for idx in 1:eeg_channel_n(eeg)
        loc_y[idx], loc_x[idx] = pol2cart(pi / 180 * eeg.eeg_header[:loc_theta][idx],
                                          eeg.eeg_header[:loc_radius][idx])
    end
    x_lim = (findmin(loc_x)[1] * 1.8, findmax(loc_x)[1] * 1.8)
    y_lim = (findmin(loc_y)[1] * 1.8, findmax(loc_y)[1] * 1.8)

    # interpolate
    x_lim_int = (findmin(loc_x)[1] * 1.4, findmax(loc_x)[1] * 1.4)
    y_lim_int = (findmin(loc_y)[1] * 1.4, findmax(loc_y)[1] * 1.4)
    interpolation_factor = 100
    interpolated_x = linspace(x_lim_int[1], x_lim_int[2], interpolation_factor)
    interpolated_y = linspace(y_lim_int[1], y_lim_int[2], interpolation_factor)
    interpolation_m = Matrix{Tuple{Float64, Float64}}(undef, interpolation_factor, interpolation_factor)
    for idx1 in 1:interpolation_factor
        for idx2 in 1:interpolation_factor
            interpolation_m[idx1, idx2] = (interpolated_x[idx1], interpolated_y[idx2])
        end
    end
    s_interpolated = zeros(interpolation_factor, interpolation_factor)
    electrode_locations = [loc_x loc_y]'
    m === :shepard && (itp = ScatteredInterpolation.interpolate(Shepard(), electrode_locations, s_non_interpolated))
    m === :mq && (itp = ScatteredInterpolation.interpolate(Multiquadratic(), electrode_locations, s_non_interpolated))
    m === :tp && (itp = ScatteredInterpolation.interpolate(ThinPlate(), electrode_locations, s_non_interpolated))
    for idx1 in 1:interpolation_factor
        for idx2 in 1:interpolation_factor
            s_interpolated[idx1, idx2] = ScatteredInterpolation.evaluate(itp, [interpolation_m[idx1, idx2][1]; interpolation_m[idx1, idx2][2]])[1]
        end
    end

    s_interpolated = signal_normalize_minmax(s_interpolated)

    p = plot(grid=false,
             framestyle=:none,
             border=:none,
             margins=0Plots.px,
             aspect_ratio=1,
             titlefontsize=10,
             xlabelfontsize=8,
             ylabelfontsize=8,
             xtickfontsize=4,
             ytickfontsize=4,
             title=title,
             size=(800, 800);
             kwargs...)
    p = plot!(interpolated_x, interpolated_y, s_interpolated, fill=:darktest, seriestype=:heatmap,
             colorbar_title="[A.U.]", clims=(-1, 1), levels=10, linewidth=0)
    p = plot!(interpolated_x, interpolated_y, s_interpolated, fill=:darktest, seriestype=:contour,
             colorbar_title="[A.U.]", clims=(-1, 1), levels=5, linecolor=:black, linewidth=0.5)
    p = plot!((loc_x, loc_y), color=:black, seriestype=:scatter, xlims=x_lim, ylims=x_lim, grid=true, label="", markersize=2, markerstrokewidth=0, markerstrokealpha=0)
    # draw head
    pts = Plots.partialcircle(0, 2π, 100, maximum(loc_x))
    x, y = Plots.unzip(pts)
    x = x .* 1.2
    y = y .* 1.2
    head = Shape(x, y)
    nose = Shape([(-0.05, maximum(y)), (0, maximum(y) + 0.1 * maximum(y)), (0.05, maximum(y))])
    ear_l = Shape([(minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), -0.08), (minimum(x) + 0.1 * minimum(x), 0.08), (minimum(x), 0.1)])
    ear_r = Shape([(maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), -0.08), (maximum(x) + 0.1 * maximum(x), 0.08), (maximum(x), 0.1)])
    for idx = 0:0.001:1
        peripheral = Shape(x .* (1 + idx), y .* (1 + idx))
        p = plot!(p, peripheral, fill=nothing, label="", linecolor=:white, linewidth=1)
    end
    p = plot!(p, head, fill=nothing, label="", linewidth=1)
    p = plot!(p, nose, fill=nothing, label="", linewidth=1)
    p = plot!(p, ear_l, fill=nothing, label="", linewidth=1)
    p = plot!(p, ear_r, fill=nothing, label="", linewidth=1)
    if head_labels == true
        p = plot!(annotation=(0, 1 - maximum(y) / 5, text("Inion", pointsize=12, halign=:center, valign=:center)))
        p = plot!(annotation=(0, -1 - minimum(y) / 5, text("Nasion", pointsize=12, halign=:center, valign=:center)))
        p = plot!(annotation=(-1 - minimum(x) / 5, 0, text("Left", pointsize=12, halign=:center, valign=:center, rotation=90)))
        p = plot!(annotation=(1 - maximum(x) / 5, 0, text("Right", pointsize=12, halign=:center, valign=:center, rotation=-90)))
    end
    p = plot!(p, xlims=(x_lim_int), ylims=(y_lim_int))

    if c !== :ica
        offset += (epoch - 1) * eeg_epoch_len(eeg) + 1
        frq_lim == (0, 0) && (frq_lim = (0, div(eeg_sr(eeg), 2)))
        average == true && (s = eeg_plot_butterfly(eeg, channel=0, len=len, offset=offset, title="Averaged signal\n[epoch: $epoch, time window: $t_s1:$t_s2]", hist=false, head=false, legend=false, average=average))
        average == false && (s = eeg_plot_butterfly(eeg, channel=0, len=len, offset=offset, title="Butterfly plot\n[epoch: $epoch, time window: $t_s1:$t_s2]", hist=false, head=false, legend=false, average=average))
        average == true && (ps = eeg_plot_psd(eeg, channel=0, len=len, offset=offset, frq_lim=frq_lim, title="Averaged PSD\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", norm=norm, legend=false, average=average))
        average == false && (ps = eeg_plot_psd(eeg, channel=0, len=len, offset=offset, frq_lim=frq_lim, title="PSD\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", norm=norm, legend=false, average=average))
        average == false && (h = eeg_plot_electrodes(eeg, channel=0, selected=1:eeg_channel_n(eeg), labels=true, head_labels=false, title="Channels\n[1:$(length(eeg_labels(eeg)))]", marksersize=1))
        average == true && (h = eeg_plot_electrodes(eeg, channel=0, selected=0, labels=true, head_labels=false, title="Channels\n1:$(eeg_channel_n(eeg, type=:eeg))", alpha=1, marksersize=1))
        l = @layout [a{0.7w} b{0.3w}; c{0.5h}; d{0.3w}]
        len >= 4 * eeg_sr(eeg) && (p = plot(ps, p, s, h, layout=(2, 2)))
        len < 4 * eeg_sr(eeg) && (p = plot(p, s, layout=(2, 1)))
    end

    plot(p)

    return p
end

"""
    signal_plot_band(signal; <keyword arguments>)

Plot absolute/relative bands powers of a single-channel `signal`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64`: sampling rate
- `band::Vector{Symbol}=[:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]`: band names, e.g. [:delta, alpha] (see `eeg_band()`)
- `type::Symbol`: plots absolute (:abs) or relative power (:rel)
- `norm::Bool=true`: convert power to dB if true
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_band(signal::Vector{Float64}; fs::Int64, band::Vector{Symbol}=[:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], type::Symbol, norm::Bool=true, xlabel::String="", ylabel::String="", title::String="", kwargs...)

    fs < 0 && throw(ArgumentError("fs must be ≥ 0."))
    type in [:abs, :rel] || throw(ArgumentError("type must be :abs or :rel."))
    band_frq = Array{Tuple{Float64, Float64}}(undef, length(band))
    for idx in 1:length(band)
        band[idx] in [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher] || throw(ArgumentError("band must be: :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower or :gamma_higher."))
        band_frq[idx] = signal_band(fs, band[idx])
    end
    for idx in 1:length(band_frq)
        band_frq[idx][1] > fs / 2 && (band_frq[idx] = (fs / 2, band_frq[idx][2]))
        band_frq[idx][2] > fs / 2 && (band_frq[idx] = (band_frq[idx][1], fs / 2))
    end

    total_pow = round(signal_total_power(signal, fs=fs), digits=2)
    abs_band_pow = zeros(length(band))
    for idx in 1:length(band)
        abs_band_pow[idx] = round(signal_band_power(signal, fs=fs, f=band_frq[idx]), digits=2)
    end
    for idx in 1:length(band)
        abs_band_pow[idx] = round(signal_band_power(signal, fs=fs, f=band_frq[idx]), digits=2)
    end
    prepend!(abs_band_pow, total_pow)
    norm == true && (abs_band_pow = pow2db.(abs_band_pow))
    rel_band_pow = abs_band_pow ./ total_pow
    labels = Array{String}(undef, (length(band) + 1))
    labels[1] = "total"
    for idx in 2:(length(band) + 1)
        labels[idx] = String(band[idx - 1])
        labels[idx] = replace(labels[idx], "_" => " ")
        labels[idx] = replace(labels[idx], "alpha" => "α")
        labels[idx] = replace(labels[idx], "beta" => "β")
        labels[idx] = replace(labels[idx], "delta" => "δ")
        labels[idx] = replace(labels[idx], "theta" => "θ")
        labels[idx] = replace(labels[idx], "gamma" => "γ")
        labels[idx] = replace(labels[idx], " high" => "\nhigh")
        labels[idx] = replace(labels[idx], " lower" => "\nlower")
        labels[idx] = replace(labels[idx], " higher" => "\nhigher")
    end
    if type === :abs
        ylabel == "" && (ylabel = "Absolute power")
        norm == true && (ylabel *= " [dB]")
        norm == false && (ylabel *= " [μV^2/Hz]")
        p = plot(labels,
                 abs_band_pow,
                 seriestype=:bar,
                 label="",
                 xlabel=xlabel,
                 ylabel=ylabel,
                 title=title,
                 palette=:darktest,
                 titlefontsize=10,
                 xlabelfontsize=8,
                 ylabelfontsize=8,
                 xtickfontsize=8,
                 ytickfontsize=8;
                 kwargs...)
    else
        ylabel == "" && (ylabel = "Relative power")
        p = plot(labels,
                 rel_band_pow,
                 seriestype=:bar,
                 label="",
                 xlabel=xlabel,
                 ylabel=ylabel,
                 title=title,
                 palette=:darktest,
                 titlefontsize=10,
                 xlabelfontsize=8,
                 ylabelfontsize=8,
                 xtickfontsize=8,
                 ytickfontsize=8;
                 kwargs...)
    end

    return p
end

"""
    signal_plot_band(signal; <keyword arguments>)

Plot absolute/relative band power of `signal` channels.

# Arguments

- `signal::Matrix{Float64}`
- `fs::Int64`: sampling rate
- `band:Symbols=:total`: band name, e.g. :delta (see `eeg_band()`)
- `type::Symbol`: plots absolute (:abs) or relative power (:rel)
- `norm::Bool=true`: convert power to dB if true
- `labels::Vector{String}=[""]`: x-axis labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function signal_plot_band(signal::Matrix{Float64}; fs::Int64, band::Symbol, type::Symbol, norm::Bool=true, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", kwargs...)

    fs < 0 && throw(ArgumentError("fs must be ≥ 0."))
    (band === :total && type === :rel) && throw(ArgumentError("for band :total, type must be :abs."))
    type in [:abs, :rel] || throw(ArgumentError("type must be :abs or :rel."))
    band in [:total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher] || throw(ArgumentError("band must be: :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower or :gamma_higher."))
    
    band_frq = signal_band(fs, band)
    channel_n = size(signal, 1)

    total_pow = zeros(channel_n)
    for idx in 1:channel_n
        total_pow[idx] = round(signal_total_power(signal[idx, :], fs=fs), digits=2)
    end

    if band !== :total
        abs_band_pow = zeros(channel_n)
        for idx in 1:channel_n
            abs_band_pow[idx] = round(signal_band_power(signal[idx, :], fs=fs, f=band_frq), digits=2)
        end
    else
        abs_band_pow = total_pow
    end

    norm == true && (total_pow = pow2db.(total_pow))
    norm == true && (abs_band_pow = pow2db.(abs_band_pow))
    type === :rel && (rel_band_pow = abs_band_pow ./ total_pow)
    if labels == [""]
        labels = Vector{String}(undef, channel_n)
        for idx in 1:channel_n
            labels[idx] = "ch " * string(idx)
        end
    end
    if type === :abs
        ylabel == "" && (ylabel = "Absolute power")
        norm == true && (ylabel *= " [dB]")
        norm == false && (ylabel *= " [μV^2/Hz]")
        p = plot(labels,
                 abs_band_pow,
                 seriestype=:bar,
                 label="",
                 xticks=(0.5:(length(labels) - 0.5), labels),
                 xlabel=xlabel,
                 ylabel=ylabel,
                 title=title,
                 palette=:darktest,
                 titlefontsize=10,
                 xlabelfontsize=8,
                 ylabelfontsize=8,
                 xtickfontsize=8,
                 ytickfontsize=8;
                 kwargs...)
    else
        ylabel == "" && (ylabel = "Relative power")
        p = plot(labels,
                 rel_band_pow,
                 seriestype=:bar,
                 label="",
                 xticks=(0.5:(length(labels) - 0.5), labels),
                 xlabel=xlabel,
                 ylabel=ylabel,
                 title=title,
                 palette=:darktest,
                 titlefontsize=10,
                 xlabelfontsize=8,
                 ylabelfontsize=8,
                 xtickfontsize=8,
                 ytickfontsize=8;
                 kwargs...)
    end

    return p
end

"""
    eeg_plot_band(eeg; <keyword arguments>)

Plots `eeg` channels. If signal is multichannel, only channel amplitudes are plotted. For single-channel signal, the histogram, amplitude, power density and spectrogram are plotted.

# Arguments

- `eeg::NeuroJ.EEG`: EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epochs to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channels to display
- `offset::Int64=0`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
- `band:Vector{Symbols}=:all`: band name, e.g. :delta (see `eeg_band()`)
- `type::Symbol`: plots absolute (:abs) or relative power (:rel)
- `norm::Bool=true`: convert power to dB if true
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_band(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}, offset::Int64=0, len::Int64=0, band::Union{Symbol, Vector{Symbol}}=:all, type::Symbol, norm::Bool=true, xlabel::String="", ylabel::String="", title::String="", kwargs...)

    (band === :all && (typeof(channel) != Int64 || length(channel) != 1)) && throw(ArgumentError("For band :all only one channel may be specified."))
    band === :all && (band = [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher])
    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before plotting."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))
    (typeof(channel) == Int64 && (channel < 1 || channel > eeg_channel_n(eeg))) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))

    (epoch != 1 && (offset != 0 || len != 0)) && throw(ArgumentError("For epoch ≠ 1, offset and len must not be specified."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    (length(epoch) == 1 && (epoch < 1 || epoch > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    (length(epoch) > 1 && (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg))) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    if length(epoch) > 1
        sort!(epoch)
        (epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
        len = eeg_epoch_len(eeg) * length(epoch)
        offset = eeg_epoch_len(eeg) * (epoch[1] - 1)
        epoch = 1
    end

    # default length is one epoch or 20 seconds
    if len == 0
        if eeg_epoch_len(eeg) > 20 * eeg_sr(eeg)
            len = 20 * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end

    # get epochs markers for len > epoch_len
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
    else
        eeg_tmp = eeg
    end

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]
    t[1] = floor(t[1], digits=2)
    t[end] = ceil(t[end], digits=2)

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
    (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + length(t)), epoch]

    if typeof(band) == Symbol
        typeof(channel) <: AbstractRange && (channel = collect(channel))
        if collect(channel[1]:channel[end]) == channel
            channel_list = string(channel[1]) * ":" * string(channel[end])
        elseif typeof(channel) != Int64
            channel_list = "" 
            for idx in 1:(length(channel) - 1)
                channel_list *= string(channel[idx])
                channel_list *= ", "
            end
            channel_list *= string(channel[end])
        else
            channel_list = string(channel)
        end
        typeof(signal) == Vector{Float64} && (signal = reshape(signal, 1, length(signal)))
        epoch_tmp = epoch
        offset > eeg_epoch_len(eeg) && (epoch_tmp = floor(Int64, offset / eeg_epoch_len(eeg)) + 1)
        title == "" && (title = "$(titlecase(string(band))) power\n[epoch: $(string(epoch_tmp)), channel(s): $channel_list]\n[offset: $offset samples, length: $len samples]")
        labels = eeg_labels(eeg)[channel]
        typeof(labels) == String && (labels = [labels])
        p = signal_plot_band(signal,
                             band=band,
                             fs=eeg_sr(eeg),
                             type=type,
                             labels=labels,
                             norm=norm,
                             xlabel=xlabel,
                             ylabel=ylabel,
                             title=title;
                             kwargs...)
    else
        signal = vec(signal)
        epoch_tmp = epoch
        offset > eeg_epoch_len(eeg) && (epoch_tmp = floor(Int64, offset / eeg_epoch_len(eeg)) + 1)
        title == "" && (title = "Band powers\n[epoch: $(string(epoch_tmp)), channel: $channel ($(eeg_labels(eeg)[channel])), offset: $offset samples, length: $len samples]")

        p = signal_plot_band(signal,
                             band=band,
                             fs=eeg_sr(eeg),
                             type=type,
                             norm=norm,
                             xlabel=xlabel,
                             ylabel=ylabel,
                             title=title;
                             kwargs...)
    end

    plot(p)

    return p
end

"""
    eeg_plot_save(p; file_name::String)

Saves plot as file (PDF/PNG/TIFF). File format is determined using `file_name` extension.

# Arguments

- `p::Plots.Plot{Plots.GRBackend}`
- `file_name::String`
"""
function eeg_plot_save(p::Plots.Plot{Plots.GRBackend}; file_name::String)

    ext = splitext(file_name)[2]
    ext in [".png", ".pdf", ".jpg", ".tiff"] || throw(ArgumentError("Fiel format must be: .png, .pdf, .tiff or .jpg"))
    isfile(file_name) && @warn "File $file_name will be overwritten."
    savefig(p, file_name)

    return
end