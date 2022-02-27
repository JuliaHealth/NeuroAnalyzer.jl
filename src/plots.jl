"""
    signal_plot(t, signal, labels="", xlabel="Time [s]", ylabel="", title="", kwargs...)

Plots `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
- `signal::AbstractArray`
- `labels::Vector{String}` - labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::AbstractArray; labels::Vector{String}=[""], xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    channel_n = eeg_channel_n(eeg)

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
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
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
    signal_plot(t, signal; ylim=(0, 0), xlabel="Time [s]", ylabel="Amplitude [μV]", title="", kwargs...)

Plots `signal` channel.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
- `signal::Vector{Float64}`
- `ylim::Tuple` - y-axis limit (-ylim:ylim)
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Vector{Float64}; ylim::Tuple=(0, 0), xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", kwargs...)

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
              xlims=(floor(t[1]), ceil(t[end])),
              xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
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
    eeg_plot(eeg; epoch=1, channel=0, offset=0, len=0, labels=[""], xlabel="Time  ylabel="Channels", title="", head=false, hist=:hist, frq_lim=(0, 0), figure="", kwargs...)

Plots `eeg` channels. If signal is multichannel, only channel amplitudes are plotted. For single-channel signal, the histogram, amplitude, power density and spectrogram are plotted.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epochs to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channels to display
- `offset::Int64` - displayed segment offset in samples
- `len::Int64` - displayed segment length in samples, default 1 epoch or 20 seconds
- `labels::Vector{String}` - channel labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `head::Bool` - add head with electrodes
- `hist::Symbol` - histogram type
- `frq_lim::Tuple` - frequency limit for PSD and spectrogram
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], xlabel::String="Time [s]", ylabel::String="", title::String="", head::Bool=true, hist::Symbol=:hist, frq_lim::Tuple=(0, 0), figure::String="", kwargs...)

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

    # get epochs markers for len > epoch_len
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_signal_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> Int(offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= Int((offset + len) / eeg_sr(eeg))]
    else
        eeg_tmp = eeg
    end

    labels = eeg_labels(eeg_tmp)

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]

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

    t_1 = t[1]
    t_2 = t[end]
    t_1 < 1.0 && (t_s1 = string(round(t_1 * 1000, digits=2)) * " ms")
    t_1 >= 1.0 && (t_s1 = string(round(t_1, digits=2)) * " s")
    t_2 < 1.0 && (t_s2 = string(round(t_2 * 1000, digits=2)) * " ms")
    t_2 >= 1.0 && (t_s2 = string(round(t_2, digits=2)) * " s")

    title == "" && (title = "Signal\n[epoch: $epoch, channel: $channel ($(eeg_labels(eeg)[channel])), time window: $t_s1:$t_s2]")

    p = signal_plot(t,
                    signal,
                    labels=labels,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    title=title;
                    kwargs...)

    # add epochs markers
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                   linestyle=:dash,
                   linewidth=0.2,
                   linecolor=:black,
                   label="")
        if typeof(signal) == Vector{Float64}
            for idx in 1:length(epoch_markers)
                p = plot!(annotation=((epoch_markers[idx] - 1), (maximum(ceil.(abs.(signal))) * 1.02), text("E$idx", pointsize=6, halign=:center, valign=:center)))
            end
        else
            for idx in 1:length(epoch_markers)
                p = plot!(annotation=((epoch_markers[idx] - 1), ((channel[end] - 1) * 1.02), text("E$idx", pointsize=6, halign=:center, valign=:center)))
            end
        end
    end

    if typeof(signal) == Vector{Float64}
        # cannot plot electrodes without locations
        psd = eeg_plot_psd(eeg, epoch=epoch, channel=channel, len=len, offset=offset, frq_lim=frq_lim, title="PSD\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", legend=false)
        eeg.eeg_header[:channel_locations] == false && (head = false)
        frq_lim == (0, 0) && (frq_lim = (0, div(eeg_sr(eeg), 2)))
        s = eeg_plot_spectrogram(eeg, epoch=epoch, channel=channel, len=len, offset=offset, frq_lim=frq_lim, title="Spectrogram\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", legend=false)
        ht_a = eeg_plot_histogram(eeg, epoch=epoch, channel=channel, len=len, offset=offset, type=hist, labels=[""], legend=false, title="Signal")
        _, _, _, s_phase = signal_spectrum(signal)
        ht_p = signal_plot_histogram(rad2deg.(s_phase), offset=offset, len=len, type=:kd, labels=[""], legend=false, title="Phase", xticks=[-180, 0, 180], linecolor=:black)
        if head == true
            hd = eeg_plot_electrodes(eeg, labels=false, selected=channel, small=true, title=channel_name)
            l = @layout [a{0.33h} b{0.2w}; c{0.33h} d{0.2w}; e{0.33h} f{0.2w}]
            p = plot(p, ht_a, psd, ht_p, s, hd, layout=l)
        else
            l = @layout [a{0.33h} b{0.2w}; c{0.33h} d{0.2w}; e{0.33h} _]
            p = plot(p, ht_a, psd, ht_p, s, layout=l)
        end
    end

    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    eeg_draw_head(p, loc_x, loc_y; head_labels=true)

Draws head over a topographical plot `p`.

# Arguments

- `p::Plot` - electrodes plot
- `loc_x::Vector{Float64}` - vector of x electrode position
- `loc_y::Vector{Float64}` - vector of y electrode position
- `head_labels::Bool` - add text labels to the plot

# Returns

- `p::Plot`
"""
function eeg_draw_head(p, loc_x::Vector{Float64}, loc_y::Vector{Float64}; head_labels::Bool=true, kwargs...)

    pts = Plots.partialcircle(0, 2π, 100, maximum(loc_x))
    x, y = Plots.unzip(pts)
    x = x .* 4
    y = y .* 4
    head = Shape(x, y)
    nose = Shape([(-0.1, maximum(y)), (0, maximum(y) + 0.1 * maximum(y)), (0.1, maximum(y))])
    ear_l = Shape([(minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), 0.1), (minimum(x), 0.1)])
    ear_r = Shape([(maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), 0.1), (maximum(x), 0.1)])
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
    eeg_plot_filter_response(eeg; fprototype, ftype, cutoff, order, rp, rs, window=nothing, figure, kwargs...)

Returns zero phase distortion filter response.

# Arguments

- `fprototype::Symbol[:butterworth, :chebyshev1, :chebyshev2, :elliptic]
- `ftype::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Union{Int64, Float64, Tuple}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64` - filter order
- `rp::Union{Int64, Float64}` - dB ripple in the passband
- `rs::Union{Int64, Float64}` - dB attenuation in the stopband
- `window::window::Union{Vector{Float64}, Nothing} - window, required for FIR filter
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_filter_response(eeg::EEG; fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Tuple}, order::Int64, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, window::Union{Vector{Float64}, Nothing}=nothing, figure::String="", kwargs...)

    fs = eeg_sr(eeg)

    ftype in [:lp, :hp, :bp, :bs] || throw(ArgumentError("Filter type must be :bp, :hp, :bp or :bs."))
    fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic] || throw(ArgumentError("Filter prototype must be :butterworth, :chebyshev1:, :chebyshev2 or :elliptic."))

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
              title="Filter: $(titlecase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order\n\nFrequency response",
              xlims=(0, x_max),
              ylabel="Magnitude [dB]",
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
              ylabel="Phase [°]",
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
              ylabel="Group delay [samples]",
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

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_avg(t, signal; norm=true, xlabel"Time [s]", ylabel="Amplitude [μV]" title="", ylim=(0, 0))

Plots averaged `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange`
- `signal::Matrix{Float64}`
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `ylim::Tuple` - y-axis limits
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
- `ylim::Tuple`
"""
function signal_plot_avg(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Matrix{Float64}; norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple=(0, 0), kwargs...)

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
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
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
    eeg_plot_avg(eeg; epoch=1, channel=0, offset=0, len=0, labels=[""], norm=true, xlabel="Time  ylabel="Channels", title="", ylim=(0, 0), frq_lim=(0, 0), head=true, figure="", kwargs...)

Plots averaged `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channel to display
- `offset::Int64` - displayed segment offset in samples
- `len::Int64` - displayed segment length in samples, default 1 epoch or 20 seconds
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `ylim::Tuple` - y-axis limits
- `frq_lim::Tuple` - frequency limit for PSD and spectrogram
- `hist::Symbol[:hist, :kd]` - histogram type
- `head::Bool` - add head plot
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_avg(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple=(0, 0), frq_lim::Tuple=(0, 0), hist::Symbol=:hist, head::Bool=true, figure::String="", kwargs...)

    (frq_lim[1] < 0 || frq_lim[1] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    (frq_lim[2] < 0 || frq_lim[2] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))

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
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_signal_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> Int(offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= Int((offset + len) / eeg_sr(eeg))]
    else
        eeg_tmp = eeg
    end

    labels = eeg_labels(eeg_temp)

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
    (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + length(t)), epoch]

    title == "" && (title = "Signal averaged\n[epoch: $epoch, channel: $channel ($(eeg_labels(eeg)[channel])), offset: $offset samples, length: $len samples]")

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
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                   linestyle=:dash,
                   linewidth=0.2,
                   linecolor=:black,
                   label="")
        for idx in 1:length(epoch_markers)
            p = plot!(annotation=((epoch_markers[idx] - 1), (ylim[2] * 1.02), text("E$idx", pointsize=6, halign=:center, valign=:center)))
        end
    end

    # cannot plot electrodes without locations
    eeg.eeg_header[:channel_locations] == false && (head = false)
    frq_lim == (0, 0) && (frq_lim = (0, div(eeg_sr(eeg), 2)))
    eeg_avg = eeg_average(eeg)
    psd = eeg_plot_psd(eeg_tmp, epoch=epoch, channel=channel, len=len, offset=offset, average=true, frq_lim=frq_lim, title="PSD averaged with 95%CI\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", legend=false)
    s = eeg_plot_spectrogram(eeg_avg, epoch=epoch, channel=1, len=len, offset=offset, frq_lim=frq_lim, title="Spectrogram averaged\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", legend=false)
    ht_a = eeg_plot_histogram(eeg_avg, epoch=epoch, channel=1, len=len, offset=offset, type=hist, labels=[""], legend=false, title="Signal")
    _, _, _, s_phase = signal_spectrum(s_normalized_m)
    ht_p = signal_plot_histogram(rad2deg.(s_phase), offset=offset, len=len, type=:kd, labels=[""], legend=false, title="Phase", xticks=[-180, 0, 180], linecolor=:black)
    if head == true
        hd = eeg_plot_electrodes(eeg, labels=false, selected=channel, small=true, title="Channels")
        l = @layout [a{0.33h} b{0.2w}; c{0.33h} d{0.2w}; e{0.33h} f{0.2w}]
        p = plot(p, ht_a, psd, ht_p, s, hd, layout=l)
    else
        l = @layout [a{0.33h} b{0.2w}; c{0.33h} d{0.2w}; e{0.33h} _]
        p = plot(p, ht_a, psd, ht_p, s, layout=l)
    end

    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_butterfly(t, signal; labels=[""], norm=true, xlabel"Time [s]", ylabel="Amplitude [μV]", title="Butterfly plot", ylim=(0, 0), kwargs...)

Butterfly plot of `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange`
- `signal::Matrix{Float64}`
- `labels::Vector{String}` - channel labels vector
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `ylim::Tuple` - y-axis limits
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_butterfly(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Butterfly plot", ylim::Tuple=(0, 0), kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    channel_n = eeg_channel_n(eeg)

    if norm == true
        s_normalized = signal_normalize_zscore(reshape(signal, eeg_channel_n(eeg), size(signal, 2), 1))
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
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
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
    eeg_plot_butterfly(eeg; epoch=1, channel=0, offset=0, len=0, labels=[""], norm=true, xlabel="Time  ylabel="Channels", title="Butterfly plot", ylim=nothing, head=true, figure="", kwargs...)

Butterfly plot of `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channel to display
- `offset::Int64` - displayed segment offset in samples
- `len::Int64` - displayed segment length in samples, default 1 epoch or 20 seconds
- `labels::Vector{String}` - channel labels vector
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `ylim::Union{Int64, Float64, Nothing}` - y-axis limits (-ylim:ylim)
- `head::Bool` - add head with electrodes
- `hist::Bool` - add histograms
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_butterfly(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Butterfly plot", ylim::Tuple=(0, 0), head::Bool=true, hist::Bool=true, figure::String="", kwargs...)

    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))

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
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_signal_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> Int(offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= Int((offset + len) / eeg_sr(eeg))]
    else
        eeg_tmp = eeg
    end

    labels = eeg_labels(eeg_temp)

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
   (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + length(t)), epoch]
    ndims(signal) == 1 && (signal = vec(signal))

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

    # add epochs markers
    if norm == true
        s_normalized = signal_normalize_zscore(reshape(signal, eeg_channel_n(eeg), size(signal, 2), 1))
    else
        s_normalized = signal
    end
    ylim == (0, 0) && (ylim = (floor(minimum(s_normalized), digits=0), ceil(maximum(s_normalized), digits=0)))
    abs(ylim[1]) > abs(ylim[2]) && (ylim = (-abs(ylim[1]), abs(ylim[1])))
    abs(ylim[1]) < abs(ylim[2]) && (ylim = (-abs(ylim[2]), abs(ylim[2])))
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                   linestyle=:dash,
                   linewidth=0.2,
                   linecolor=:black,
                   label="")
        for idx in 1:length(epoch_markers)
            p = plot!(annotation=((epoch_markers[idx] - 1), (maximum(ceil.(abs.(signal))) * 1.02), text("E$idx", pointsize=6, halign=:center, valign=:center)))
        end
    end

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_psd(s_powers, s_freqs; frq_lim=(0, 0), xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="", kwargs...)

Plots power spectrum density.

# Arguments

- `s_powers::Vector{Float64}`
- `s_freqs::Vector{Float64}`
- `frq_lim::Tuple` - x-axis limit
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_psd(s_powers::Vector{Float64}, s_freqs::Vector{Float64}; frq_lim::Tuple=(0, 0), xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="", kwargs...)

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
    signal_plot_psd(signal; fs, norm=false, frq_lim=nothing, xlabel="Frequency [Hz]", ylabel="Power [μV^2/ title="", kwargs...)

Plots power spectrum density.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64` - sampling frequency
- `norm::Bool` - converts power to dB
- `frq_lim::Tuple` - x-axis limit
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_psd(signal::Vector{Float64}; fs::Int64, norm::Bool=false, frq_lim::Tuple=(0, 0), xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="", kwargs...)

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
    signal_plot_psd(signal; fs, norm=false, average=false, frq_lim=nothing, labels=[""], xlabel="Frequency [ ylabel="Power [μV^2/Hz]", title="", kwargs...)

Plots power spectrum density.

# Arguments

- `signal::Matrix{Float64}`
- `fs::Int64` - sampling rate
- `norm::Bool` - power in dB
- `average::Bool` - plots average power and 95%CI for all channels
- `frq_lim::Tuple` - x-axis limit
- `labels::Vector{String}` - channel labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false, average::Bool=false, frq_lim::Tuple=(0, 0), labels::Vector{String}=[""], xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="", kwargs...)

    fs < 0 && throw(ArgumentError("fs must be ≥ 0."))
    (frq_lim[1] < 0 || frq_lim[2]) < 0 && throw(ArgumentError("frq_lim must be ≥ 0."))
    frq_lim[1] > frq_lim[2] && (frq_lim = (frq_lim[2], frq_lim[1]))

    norm == true && (ylabel="Power [dB]")

    channel_n = eeg_channel_n(eeg)
    signal = reshape(signal, eeg_channel_n(eeg), size(signal, 2), 1)
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
    eeg_plot_psd(eeg; epoch=1, channel=0, offset=0, len=0, labels=[""], norm=false, average=false, frq_lim=not xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="", head=false, figure="", kwargs...)

Plots power spectrum density.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channel to display
- `offset::Int64` - displayed segment offset in samples
- `len::Int64` - displayed segment length in samples, default 1 epoch or 20 seconds
- `labels::Vector{String}` - channel labels vector
- `norm::Bool` - power in dB
- `average::Bool` - plots average power and 95%CI for all channels
- `frq_lim::Tuple` - x-axis limit
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `head::Bool` - add head with electrodes
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_psd(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], norm::Bool=false, average::Bool=false, frq_lim::Tuple=(0, 0), xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="", head::Bool=false, figure::String="", kwargs...)

    (frq_lim[1] < 0 || frq_lim[1] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    (frq_lim[2] < 0 || frq_lim[2] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))

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

    labels = eeg_labels(eeg_temp)

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

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    eeg_plot_electrodes(eeg; channel=0, selected=0, labels=true, head=true, head_labels=false, small=false)

Plots electrodes.

# Arguments

- `eeg:EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channel to display
- `selected::Union{Int64, Vector{Int64}, AbstractRange}` - which channel should be highlighted
- `labels::Bool` - plot electrode labels
- `head::Bool` - plot head
- `head_labels::Bool` - plot head labels
- `small::Bool` - draws small plot

# Returns

- `p::Plot`
"""
function eeg_plot_electrodes(eeg::EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head::Bool=true, head_labels::Bool=false, small::Bool=false, kwargs...)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() first."))

    # select channels, default is all channels
    channel == 0 && (channel = 1:eeg_channel_n(eeg))
    length(channel) > 1 && sort!(channel)
    for idx in 1:length(channel)
        (channel[idx] < 1 || channel[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
    typeof(channel) <: AbstractRange && (channel = collect(channel))

    # select channels, default is all channels
    selected == 0 && (selected = 1:eeg_channel_n(eeg))
    length(selected) > 1 && sort!(selected)
    for idx in 1:length(selected)
        (selected[idx] < 1 || selected[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
    typeof(selected) <: AbstractRange && (selected = collect(selected))
    length(selected) > 1 && (intersect(selected, channel) == selected || throw(ArgumentError("channel must include selected.")))
    length(selected) == 1 && (intersect(selected, channel) == [selected] || throw(ArgumentError("channel must include selected.")))

    eeg_tmp = eeg_keep_channel(eeg, channel=channel)

    loc_x = eeg_tmp.eeg_header[:xlocs]
    loc_y = eeg_tmp.eeg_header[:ylocs]
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

    p = plot(grid=false, framestyle=:none, palette=:darktest, size=plot_size, markerstrokewidth=0, border=:none, aspect_ratio=1, margins=-20Plots.px, titlefontsize=10; kwargs...)
    if length(selected) == eeg_tmp.eeg_header[:channel_n]
        for idx in 1:eeg_tmp.eeg_header[:channel_n]
            p = plot!((loc_x[idx], loc_y[idx]), color=idx, seriestype=:scatter, xlims=x_lim, ylims=x_lim, grid=true, label="", markersize=marker_size, markerstrokewidth=0, markerstrokealpha=0)
        end
    else
        p = plot!(loc_x, loc_y, seriestype=:scatter, color=:black, alpha=0.2, xlims=x_lim, ylims=y_lim, grid=true, label="", markersize=marker_size, markerstrokewidth=0, markerstrokealpha=0)
        eeg_tmp = eeg_keep_channel(eeg, channel=selected)
        loc_x = eeg_tmp.eeg_header[:xlocs]
        loc_y = eeg_tmp.eeg_header[:ylocs]
        for idx in 1:eeg_tmp.eeg_header[:channel_n]
            p = plot!((loc_x[idx], loc_y[idx]), color=idx, seriestype=:scatter, xlims=x_lim, ylims=x_lim, grid=true, label="", markersize=marker_size, markerstrokewidth=0, markerstrokealpha=0)
        end
    end
    if labels == true
        for idx in 1:length(eeg_labels(eeg_temp))
        plot!(annotation=(loc_x[idx], loc_y[idx] + 0.05, text(eeg_labels(eeg_temp)[idx], pointsize=font_size)))
        end
        p = plot!()
    end

    if head == true
        # for some reason head is enlarged for channel > 1
        eeg_tmp = eeg_keep_channel(eeg, channel=1)
        loc_x = eeg_tmp.eeg_header[:xlocs]
        loc_y = eeg_tmp.eeg_header[:ylocs]
        hd = eeg_draw_head(p, loc_x, loc_x, head_labels=head_labels)
        plot!(hd)
    end

    plot(p)

    return p
end

"""
    eeg_plot_matrix(eeg, m; epoch, figure="", kwargs...)

Plots matrix `m` of `eeg` signals.

# Arguments

- `eeg:EEG`
- `m::Union{Matrix{Float64}, Array{Float64, 3}}` - channels by channels matrix
- `epoch::Int64` - epoch number to display
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_matrix(eeg::EEG, m::Union{Matrix{Float64}, Array{Float64, 3}}; epoch::Int64=1, figure::String="", kwargs...)
    (epoch < 1 || epoch > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))

    labels = eeg_labels(eeg)
    channel_n = size(m, 1)
    ndims(m) == 3 && (m = m[:, :, epoch])

    p = heatmap(m, xticks=(1:channel_n, labels), yticks=(1:channel_n, eeg_labels(eeg)), xlabelfontsize=8,
             ylabelfontsize=8, xtickfontsize=4, ytickfontsize=4; kwargs...)
    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    eeg_plot_matrix(eeg, cov_m, lags; epoch, figure="", kwargs...)

Plots matrix `m` of `eeg` signals.

# Arguments

- `eeg:EEG`
- `cov_m::Union{Matrix{Float64}, Array{Float64, 3}}` - covariance matrix
- `lags::Union{Vector{Int64}, Vector{Float64}}` - covariance matrix
- `channel::Union{Int64, Vector{Int64}, UnitRange{Int64}, Nothing}` - channel to display
- `epoch::Int64` - epoch number to display
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_covmatrix(eeg::EEG, cov_m::Union{Matrix{Float64}, Array{Float64, 3}}, lags::Union{Vector{Int64}, Vector{Float64}}; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch::Int64=1, figure::String="", kwargs...)

    (epoch < 1 || epoch > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))

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

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_spectrogram(signal; fs, offset=0, norm=true, demean=true, ylim=(0, 0), xlabel="Time  ylabel="Frequency [Hz]", title="Spectrogram", kwargs...)

Plots spectrogram of `signal`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64` - sampling frequency
- `offset::Int64` - displayed segment offset in samples
- `norm::Bool` - normalize powers to dB
- `demean::Bool` - demean signal prior to analysis
- `ylim::Tuple` - y-axis limits
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_spectrogram(signal::Vector{Float64}; fs::Int64, offset::Int64=0, norm::Bool=true, demean::Bool=true, frq_lim::Tuple=(0, 0), xlabel="Time [s]", ylabel="Frequency [Hz]", title="Spectrogram", kwargs...)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1 Hz."))
    frq_lim[1] > fs / 2 || frq_lim[2] > fs / 2 && throw(ArgumentError("frq_lim must be smaller than Nyquist frequency ($(fs/2) Hz)."))
    (frq_lim[1] < 0 || frq_lim[2]) < 0 && throw(ArgumentError("frq_lim must be ≥ 0."))
    frq_lim[1] > frq_lim[2] && (frq_lim = (frq_lim[2], frq_lim[1]))
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
                    xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
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
                    xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
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
    eeg_plot_spectrogram(eeg; epoch=1, channel, offset=0, len=0, norm=true, frq_lim=(0, 0), xlabel="Time [s]", ylabel="Frequency [Hz]", title="Spectrogram", kwargs...)

Plots spectrogram of `eeg` channel.

# Arguments

- `eeg:EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epoch to plot
- `channel::Int64` - channel to plot
- `offset::Int64` - displayed segment offset in samples
- `len::Int64` - displayed segment length in samples, default 1 epoch or 20 seconds
- `norm::Bool` - normalize powers to dB
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `frq_lim::Tuple` - y-axis limits
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_spectrogram(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Int64, offset::Int64=0, len::Int64=0, norm::Bool=true, frq_lim::Tuple=(0, 0), xlabel::String="Time [s]", ylabel::String="Frequency [Hz]", title::String="Spectrogram", figure="", kwargs...)

    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))

    (frq_lim[1] < 0 || frq_lim[1] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    (frq_lim[2] < 0 || frq_lim[2] > eeg_sr(eeg) / 2) && throw(ArgumentError("frq_lim must be > 0 Hz and ≤ $(eeg_sr(eeg))."))
    frq_lim == (0, 0) && (frq_lim = (0, eeg_sr(eeg) / 2))
    frq_lim[1] > frq_lim[2] && (frq_lim = (frq_lim[2], frq_lim[1]))
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
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_signal_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> Int(offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= Int((offset + len) / eeg_sr(eeg))]
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
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                   linestyle=:dash,
                   linewidth=0.2,
                   linecolor=:black,
                   label="")
        for idx in 1:length(epoch_markers)
            p = plot!(annotation=((epoch_markers[idx] - 1), frq_lim[2] * 1.05, text("E$idx", pointsize=6, halign=:center, valign=:center)))
        end
    end

    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_histogram(signal; type=:hist, label="", xlabel="", ylabel="", title="", kwargs...)

Plots histogram of `signal`.

# Arguments

- `signal::Vector{Float64}`
- `type::Symbol[:hist, :kd]` - type of histogram: regular `:hist` or kernel density `:kd`
- `label::String` - channel label
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
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
    signal_plot_histogram(signal, type=:hist, labels=[""], xlabel="", ylabel="", title="", kwargs...)

Plots histogram of `signal`.

# Arguments

- `signal::Matrix{Float64}`
- `type::Symbol[:hist, :kd]` - type of histogram: regular `:hist` or kernel density `:kd`
- `labels::Vector{String}`
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_histogram(signal::Union{Vector{Float64}, Matrix{Float64}}; type::Symbol=:hist, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", kwargs...)

    channel_n = eeg_channel_n(eeg)

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
    eeg_plot_histogram(eeg; epoch=1, channel, offset=0, len=0, labels=[""], xlabel="", ylabel="", title="", figure="", kwargs...)

Plots `eeg` channels histograms.

# Arguments

- `eeg::EEG` - EEG object
- `type::Symbol[:hist, :kd]` - type of histogram: regular `:hist` or kernel density `:kd`
- `epoch::Int64` - epoch number to display
- `channel::Int64` - channel to display
- `offset::Int64` - displayed segment offset in samples
- `len::Int64` - displayed segment length in samples, default 1 epoch or 20 seconds
- `labels::String` - channel label
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_histogram(eeg::EEG; type::Symbol=:hist, epoch::Int64=1, channel::Int64, offset::Int64=0, len::Int64=0, label::String="", xlabel::String="", ylabel::String="", title::String="", figure::String="", kwargs...)

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

    label == "" && (label = eeg_labels(eeg_temp)[channel])

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

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_ica(t, ica; label="", norm=true, xlabel="Time [s]", ylabel="Amplitude [μV]", title="ICA", ylim=nothing, kwargs...)

Plots `ica` against time vector `t`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}` - the time vector
- `ica::Vector{Float64}`
- `label::String` - channel label
- `norm::Bool` - normalize the `ica` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `ylim::Union{Int64, Float64, Nothing}` - y-axis limits (-ylim:ylim)
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_ica(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, ica::Vector{Float64}; label::String="", norm::Bool=true, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="ICA", ylim::Union{Int64, Float64, Nothing}=nothing, kwargs...)

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
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
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
    signal_plot_ica(t, ica, labels="", norm=true, xlabel="Time [s]", ylabel="", title="ICA", kwargs...)

Plots `ica`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
- `ica::Matrix{Float64}`
- `labels::Vector{String}` - labels vector
- `norm::Bool` - normalize the `ica` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_ica(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange,}, ica::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="ICA", kwargs...)

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
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
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
    eeg_plot_ica(eeg; epoch=1, offset=0, len=0, ic=nothing, norm=true, xlabel="Time  ylabel="", title="ICA", figure="", kwargs...)

Plots ICs.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Int64` - epoch number to display
- `offset::Int64` - displayed segment offset in samples
- `len::Int64` - displayed segment length in samples, default 1 epoch or 20 seconds
- `ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}` - which IC to plot
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_ica(eeg::EEG; epoch::Int64=1, offset::Int64=0, len::Int64=0, ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="ICA", figure::String="", kwargs...)

    :ica_activations in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain ICA. Perform eeg_ica!(EEG) first."))
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
    ica_idx = findfirst(isequal(:ica_activations), eeg.eeg_header[:components])
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
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_signal_len(eeg)
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
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        p = vline!(epoch_markers,
                  linestyle=:dash,
                  linewidth=0.5,
                  linecolor=:black,
                  label="")
        for idx in 1:(floor(Int64, (offset + len) / eeg_epoch_len(eeg)))
            p = plot!(annotation=((epoch_markers[idx] - 1), ((length(ic) - 1) * 1.02), text("E$idx", pointsize=6, halign=:center, valign=:center)))
        end
    end

    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    eeg_plot_topo(eeg; offset, len=0, m=:shepard, c=:amp, c_idx=nothing, norm=true, frq_lim=(0,0) head_labels=false, cb_label="", title="", figure="", kwargs...)

Plots topographical view of `eeg` component.

# Arguments

- `eeg::EEG`
- `offset::Int64` - time (in samples) at which to plot
- `len::Int64` - interpolation window
- `m::Symbol[:shepard, :mq, :tp]` - interpolation method: Sherad, Multiquadratic, ThinPlate
- `c::Symbol` - component name (:ica, :pca, :amp, :power)
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange, Tuple, Nothing}` - component index, e.g. ICA number or frequency range
- `norm::Bool` - convert power as dB
- `frq_lim::Tuple` - frequency limit for PSD and spectrogram
- `head_labels::Bool` - plot head labels
- `cb_label::String` - color bar label
- `title::String` - plot title
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `plot`
"""
function eeg_plot_topo(eeg::EEG; offset::Int64, len::Int64=0, m::Symbol=:shepard, c::Symbol=:amp, c_idx::Union{Int64, Vector{Int64}, AbstractRange, Tuple, Nothing}=nothing, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0,0), head_labels::Bool=false, cb_label::String="", title::String="", figure::String="", kwargs...)

    m in [:shepard, :mq, :tp] || throw(ArgumentError("m must be :shepard, :mq or :tp."))
    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() first."))
    offset < 0 || offset > eeg_signal_len(eeg)  && throw(ArgumentError("offset must be ≥ 0 and ≤ $(eeg_signal_len(eeg))."))
    (c === :amp  || c === :power || c in eeg.eeg_header[:components]) || throw(ArgumentError("Component $(c) not found."))
    frq_lim = tuple_order(frq_lim)

    # default length is 100 ms
    len == 0 && (len = round(Int64, eeg_sr(eeg) / 10))
    len < 0 && throw(ArgumentError("len must be > 0."))
    offset + len > eeg_signal_len(eeg) && throw(ArgumentError("offset + len must be ≤ $(eeg_signal_len(eeg))."))

    # get epochs markers for len > epoch_len
    t_tmp = offset
    if offset > (eeg_epoch_len(eeg) - len) && eeg_epoch_n(eeg) > 1
        epoch = 0
        while offset > eeg_epoch_len(eeg)
            epoch += 1
            offset -= eeg_epoch_len(eeg)
        end
    else
        epoch = 1
    end

    typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
    # ignore c_idx for components other than ICA
    if typeof(c_idx) == Vector{Int64} && c === :ica
        component_idx = findfirst(isequal(c), eeg.eeg_header[:components])
        for idx in length(c_idx):-1:1
            c_idx[idx] > size(eeg.eeg_components[component_idx], 1) && throw(ArgumentError("For component $(c) range must be 1:$(size(eeg.eeg_components[component_idx], 1))"))
        end
        p = []
        length(c_idx) <= 10 && (l_row = 2)
        length(c_idx) > 10 && (l_row = 4)
        l = (l_row, ceil(Int64, length(c_idx) / l_row))
        for idx in 1:length(c_idx)
            push!(p, eeg_plot_topo(eeg; offset=offset, len=len, m=m, c=c, c_idx=c_idx[idx], head_labels=head_labels, title=""))
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

    t_1 = eeg.eeg_time[t_tmp]
    t_2 = eeg.eeg_time[t_tmp + len]
    t_1 < 1.0 && (t_s1 = string(round(t_1 * 1000, digits=2)) * " ms")
    t_1 >= 1.0 && (t_s1 = string(round(t_1, digits=2)) * " s")
    t_2 < 1.0 && (t_s2 = string(round(t_2 * 1000, digits=2)) * " ms")
    t_2 >= 1.0 && (t_s2 = string(round(t_2, digits=2)) * " s")

    if c === :amp
        s_non_interpolated = mean(eeg.eeg_signals[:, offset:(offset + len), epoch], dims=2)
        title = "Unweighted amplitude [A.U.]"
    elseif c === :ica
        s = eeg.eeg_signals[:, :, epoch]
        s = reshape(s, size(s, 1), size(s, 2), 1)
        component_idx = findfirst(isequal(c), eeg.eeg_header[:components])
        a = eeg.eeg_components[component_idx][c_idx, :, epoch]
        component_idx = findfirst(isequal(:ica_mw), eeg.eeg_header[:components])
        m_v = eeg.eeg_components[component_idx][:, c_idx, epoch]
        s_w = m_v * a'
        s_non_interpolated = mean(s_w[:, offset:(offset + len), epoch], dims=2)
        title = "$(uppercase(string(c))) #$(c_idx) [A.U.]"
    elseif c === :pca
        s = eeg.eeg_signals[:, :, epoch]
        s = reshape(s, size(s, 1), size(s, 2), 1)
        pca_idx = findfirst(isequal(:pca), eeg.eeg_header[:components])
        pca = eeg.eeg_components[pca_idx][:, :, epoch]
        pca_m_idx = findfirst(isequal(:pca_m), eeg.eeg_header[:components])
        pca_m = eeg.eeg_components[pca_m_idx]
        s_reconstructed = reconstruct(pca_m, pca)
        s_non_interpolated = mean(s_reconstructed[:, offset:(offset + len), epoch], dims=2)
        title = "$(uppercase(string(c))) #1:$(size(pca, 1)) reconstruction [A.U.]"
    elseif c === :power
        if typeof(c_idx) <: Tuple
            s = eeg.eeg_signals[:, offset:(offset + len), epoch]
            s_non_interpolated = zeros(size(s, 1))
            for idx in 1:size(s_non_interpolated, 1)
                s_non_interpolated[idx] = signal_band_power(s[idx, :], fs=eeg_sr(eeg), f=c_idx)
            end
            norm == true && s_non_interpolated[s_non_interpolated .<= 0] .= eps()
            norm == true && (s_non_interpolated = pow2db.(s_non_interpolated))
            title = "Power [A.U.]\n[frequency range: $(c_idx[1]):$(c_idx[end]) Hz]"
        elseif typeof(c_idx) == Int64 || typeof(c_idx) == Float64
            s = eeg.eeg_signals[:, offset:(offset + len), epoch]
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
        throw(ArgumentError("Component $(c) not found."))
    end

    # plot signal at electrodes at time
    loc_x = eeg.eeg_header[:xlocs]
    loc_y = eeg.eeg_header[:ylocs]
    x_lim = (findmin(loc_x)[1] * 1.8, findmax(loc_x)[1] * 1.8)
    y_lim = (findmin(loc_y)[1] * 1.8, findmax(loc_y)[1] * 1.8)

    # interpolate
    x_lim_int = (findmin(loc_x)[1] * 1.4, findmax(loc_x)[1] * 1.4)
    y_lim_int = (findmin(loc_y)[1] * 1.4, findmax(loc_y)[1] * 1.4)
    interpolation_factor = 200
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
    p = plot!(interpolated_x, interpolated_y, s_interpolated, fill=:darktest, seriestype=:contourf,
             colorbar_title=cb_label, clims=(-1, 1))
    p = plot!((loc_x, loc_y), color=:black, seriestype=:scatter, xlims=x_lim, ylims=x_lim, grid=true, label="", markersize=4, markerstrokewidth=0, markerstrokealpha=0)
    # draw head
    pts = Plots.partialcircle(0, 2π, 100, maximum(loc_x))
    x, y = Plots.unzip(pts)
    x = x .* 1.2
    y = y .* 1.2
    head = Shape(x, y)
    nose = Shape([(-0.1, maximum(y)), (0, maximum(y) + 0.1 * maximum(y)), (0.1, maximum(y))])
    ear_l = Shape([(minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), 0.1), (minimum(x), 0.1)])
    ear_r = Shape([(maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), 0.1), (maximum(x), 0.1)])
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
        frq_lim == (0, 0) && (frq_lim = (0, div(eeg_sr(eeg), 2)))
        s = eeg_plot_butterfly(eeg, epoch=epoch, channel=0, len=len, offset=offset, title="Butterfly plot\n[time window: $t_s1:$t_s2]", legend=false)
        ps = eeg_plot_psd(eeg, epoch=epoch, channel=0, len=len, offset=offset, frq_lim=frq_lim, title="PSD\n[frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]", norm=norm, legend=false)
        h = eeg_plot_electrodes(eeg, channel=0, selected=0, labels=true, head_labels=false, title="Channels\n[1:$(length(eeg_labels(eeg)))]")
        l = @layout [a{0.5h} b{0.3w}; c{0.5h}; d{0.3w}]
        p = plot(ps, p, s, h, layout=(2, 2))
    end

    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_bands(signal; fs, band=:all, type)

Plots absolute/relative band powers of `signal`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64` - sampling rate
- `band:Vector{Symbols}` - band name, e.g. :delta (see `eeg_band()`)
- `type::Symbol[:abs, :rel]` - plots absolute or relative power
- `norm::Bool` - convert power to dB if true
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `plot`
"""
function signal_plot_bands(signal::Vector{Float64}; fs::Int64, band::Union{Symbol, Vector{Symbol}}=:all, type::Symbol, norm::Bool=true, xlabel::String="", ylabel::String="", title::String="", kwargs...)

    fs < 0 && throw(ArgumentError("fs must be ≥ 0."))
    type in [:abs, :rel] || throw(ArgumentError("type must be :abs or :rel."))
    band === :all && (band = [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher])
    band_frq = Array{Tuple{Float64, Float64}}(undef, length(band))
    for idx in 1:length(band)
        band[idx] in [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher] || throw(ArgumentError("Available bands: :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher."))
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
    eeg_plot_bands(eeg; epoch=1, channel, offset=0, len=0, labels=[""], xlabel="Time  ylabel="Channels", title="", figure="", kwargs...)

Plots `eeg` channels. If signal is multichannel, only channel amplitudes are plotted. For single-channel signal, the histogram, amplitude, power density and spectrogram are plotted.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epochs to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channels to display
- `offset::Int64` - displayed segment offset in samples
- `len::Int64` - displayed segment length in samples, default 1 epoch or 20 seconds
- `band:Vector{Symbols}` - band name, e.g. :delta (see `eeg_band()`)
- `type::Symbol[:abs, :rel]` - plots absolute or relative power
- `norm::Bool` - convert power to dB if true
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_bands(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Int64, offset::Int64=0, len::Int64=0, band::Union{Symbol, Vector{Symbol}}=:all, type::Symbol, norm::Bool=true, xlabel::String="", ylabel::String="", title::String="", figure::String="", kwargs...)

    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    len < 0 && throw(ArgumentError("len must be > 0."))
    (channel < 1 || channel > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))

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
        epoch_len = eeg_signal_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> Int(offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= Int(len / eeg_sr(eeg))]
    else
        eeg_tmp = eeg
    end

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]

    (offset < 0 || offset > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg_tmp))."))
    (offset + len > eeg_epoch_len(eeg_tmp)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg_tmp))."))

    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + length(t)), epoch]
    signal = vec(signal)

    title == "" && (title = "Band powers\n[epoch: $epoch, channel: $channel ($(eeg_labels(eeg)[channel])), offset: $offset samples, length: $len samples]")

    p = signal_plot_bands(signal,
                          band=band,
                          fs=eeg_sr(eeg),
                          type=type,
                          norm=norm,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title;
                          kwargs...)

    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure will be overwritten."
        savefig(p, figure)
    end

    return p
end