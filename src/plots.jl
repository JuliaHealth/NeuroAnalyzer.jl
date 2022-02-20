"""
    signal_plot(t, signal; label="", norm=true, xlabel="Time [s]", ylabel="Amplitude [μV]", title="Signal plot", ylim=nothing, kwargs...)

Plots `signal` against time vector `t`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}` - the time vector
- `signal::Vector{Float64}`
- `label::String` - channel label
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `ylim::Union{Int64, Float64, Nothing}` - y-axis limits (-ylim:ylim)
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Vector{Float64}; label::String="", norm::Bool=true, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Signal plot", ylim::Union{Int64, Float64, Nothing}=nothing, kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    if ylim === nothing
        ylim = maximum(signal) * 1.5
        ylim = ceil(Int64, ylim)
    end

    p = plot(t,
             signal[1:length(t)],
             xlabel=xlabel,
             ylabel=ylabel,
             label=label,
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
             ylims=(-ylim, ylim),
             title=title,
             palette=:darktest;
             kwargs...)

    return p
end

"""
    signal_plot(t, signal, labels="", norm=true, xlabel="Time [s]", ylabel="Channels", title="Signal plot", kwargs...)

Plots `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
- `signal::Matrix{Float64}`
- `labels::Vector{String}` - labels vector
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=true, xlabel::String="Time [s]", ylabel::String="Channels", title::String="Signal plot", kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    channel_n = size(signal, 1)

    # reverse so 1st channel is on top
    channel_color = channel_n:-1:1
    signal = reverse(signal[:, :], dims = 1)
    s_normalized = zeros(size(signal))

    if norm == true
        # normalize and shift so all channels are visible
        variances = var(signal, dims=2)
        mean_variance = mean(variances)
        for idx in 1:channel_n
            s = @view signal[idx, :]
            s_normalized[idx, :] = (s .- mean(s)) ./ mean_variance .+ (idx - 1)
        end
    else
        s_normalized = signal
    end

    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
             ylims=(-0.5, channel_n-0.5),
             title=title,
             palette=:darktest;
             kwargs...)
    for idx in 1:channel_n
        p = plot!(t,
                  s_normalized[idx, 1:length(t)],
                  label="", color=channel_color[idx])
    end
    p = plot!(p, yticks=((channel_n - 1):-1:0, labels))

    return p
end

"""
    eeg_plot(eeg; epoch=1, channel=nothing, offset=0, len=nothing, labels=[""], norm=true, xlabel="Time  ylabel="Channels", title="Signal plot", head=false, figure="", kwargs...)

Plots `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Int64` - epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}` - channel to display
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Nothing}` - displayed segment length in samples, default 10 seconds
- `labels::Vector{String}` - channel labels vector
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `head::Bool` - add head with electrodes
- `hist::Bool` - add channel(s) histogram(s)
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot(eeg::EEG; epoch::Int64=1, channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, offset::Int64=0, len::Union{Int64, Nothing}=nothing, labels::Vector{String}=[""], norm::Bool=true, xlabel::String="Time [s]", ylabel::String="Channels", title::String="Signal plot", head::Bool=true, hist::Bool=true, figure::String="", kwargs...)
    offset < 0 && throw(ArgumentError("Offset must be ≥ 0."))
    (len !== nothing && len <= 0) && throw(ArgumentError("len must be > 0."))
    len === nothing && (len = 10 * eeg_sr(eeg))
    (epoch < 1 || epoch > eeg.eeg_header[:epoch_n]) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg.eeg_header[:epoch_n])."))
    (typeof(channel) == Int64 && (channel < 1 || channel > eeg.eeg_header[:channel_n])) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg.eeg_header[:chansnel_n])."))
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    if typeof(channel) == Vector{Int64}
        sort!(channel)
        for idx in 1:length(channel)
            (channel[idx] < 1 || channel[idx] > eeg.eeg_header[:channel_n]) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
        end
    end

    # select channels, default is 1:20 or all channels
    if channel === nothing
        if eeg.eeg_header[:channel_n] >= 20
            channel = 1:20
        else
            channel = 1:eeg.eeg_header[:channel_n]
        end
    end

    # get epochs markers for len > epoch_len
    if (((len / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds] && eeg.eeg_header[:epoch_n] > 1)
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
        eeg_tmp = eeg_epochs(eeg_tmp, epoch_n=1)
        epoch_len = size(eeg.eeg_signals, 2)
        epoch_n = size(eeg.eeg_signals, 3)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
    else
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
    end

    labels = eeg_tmp.eeg_header[:labels]

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]

    (offset < 0 || offset > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))
    (offset + len > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset + len must be ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))

    signal = eeg_tmp.eeg_signals[:, (1 + offset):(offset + length(t)), epoch]

    p = signal_plot(t,
                    signal,
                    labels=labels,
                    norm=norm,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    title=title;
                    kwargs...)

    # add epochs markers
    if (((len / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds] && eeg.eeg_header[:epoch_n] > 1)
        p = vline!(p,
                   epoch_markers,
                   timeseries=:vline,
                   linestyle=:dash,
                   linewidth=0.5,
                   linecolor=:grey,
                   label="")
    end

    # cannot plot electrodes without locations
    eeg.eeg_header[:channel_locations] == false && (head = false)

    e = plot(grid=false, framestyle=:none)
    head == true && (hd = eeg_plot_electrodes(eeg, labels=false, selected=channel, small=true))
    head == false && (hd = e)
    hist == true && (ht = signal_plot_histogram(signal, type=:hist, labels=[""]))
    hist == false && (ht = e)
    l = @layout [a b{0.1w}; c{0.9h} d]
    (head == true || hist == true) && (p = plot!(e, hd, p, ht, layout=l, left_margin=20px, bottom_margin=10px))
    head == false && hist == false && (p = plot!(p))

    if figure !== ""
        isfile(figure) && @warn "File $figure cannot be saved."
        savefig(p, figure)
    end

    return p
end

"""
    eeg_draw_head(p, loc_x, loc_y; head_labels=true)

Draws head over a topographical plot `p`.

# Arguments

- `p::Plot` - topographical plot
- `loc_x::Vector{Float64}` - vector of x electrode position
- `loc_y::Vector{Float64}` - vector of y electrode position
- `head_labels::Bool` - add text labels to the plot

# Returns

- `p::Plot`
"""
function eeg_draw_head(p, loc_x::Vector{Float64}, loc_y::Vector{Float64}; head_labels::Bool=true, kwargs...)
    pts = Plots.partialcircle(0, 2π, 100, maximum(loc_x))
    x, y = Plots.unzip(pts)
    x = x .* 1.2
    y = y .* 1.2
    head = Shape(x, y)
    nose = Shape([(-0.1, maximum(y)), (0, maximum(y) + 0.1 * maximum(y)), (0.1, maximum(y))])
    ear_l = Shape([(minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), 0.1), (minimum(x), 0.1)])
    ear_r = Shape([(maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), 0.1), (maximum(x), 0.1)])
    plot!(p, head, fill=nothing, label="")
    plot!(p, nose, fill=nothing, label="")
    plot!(p, ear_l, fill=nothing, label="")
    plot!(p, ear_r, fill=nothing, label="")
    if head_labels == true
        plot!(p, annotation=(0, 1 - maximum(y) / 5, text("Inion", pointsize=12, halign=:center, valign=:center)))
        plot!(p, annotation=(0, -1 - minimum(y) / 5, text("Nasion", pointsize=12, halign=:center, valign=:center)))
        plot!(p, annotation=(-1 - minimum(x) / 5, 0, text("Left", pointsize=12, halign=:center, valign=:center, rotation=90)))
        plot!(p, annotation=(1 - maximum(x) / 5, 0, text("Right", pointsize=12, halign=:center, valign=:center, rotation=-90)))
    end
    plot!(p; kwargs...)
end

"""
    plot_filter_response(fprototype, ftype, cutoff, fs, order, rp, rs, window, figure, kwargs...)

Returns zero phase distortion filter response.

# Arguments

- `fprototype::Symbol[:butterworth, :chebyshev1, :chebyshev2, :elliptic]
- `ftype::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `fs::Int64` - sampling rate
- `order::Int64` - filter order
- `rp::Union{Int64, Float64, Nothing}` - dB ripple in the passband
- `rs::Union{Int64, Float64, Nothing}` - dB attenuation in the stopband
- `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function plot_filter_response(;fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}, fs::Int64, order::Int64, rp::Union{Int64, Float64, Nothing}=nothing, rs::Union{Int64, Float64, Nothing}=nothing, window::Union{Vector{Float64}, Nothing}=nothing, figure::String="", kwargs...)
    ftype in [:lp, :hp, :bp, :bs] || throw(ArgumentError("Filter type must be :bp, :hp, :bp or :bs."))
    fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic] || throw(ArgumentError("Filter prototype must be :butterworth, :chebyshev1:, :chebyshev2 or :elliptic."))

    if ftype === :lp
        length(cutoff) > 1 && throw(ArgumentError("For low-pass filter one frequency must be given."))
        responsetype = Lowpass(cutoff; fs=fs)
    elseif ftype === :hp
        length(cutoff) > 1 && throw(ArgumentError("For high-pass filter one frequency must be given."))
        responsetype = Highpass(cutoff; fs=fs)
    elseif ftype === :bp
        responsetype = Bandpass(cutoff[1], cutoff[2]; fs=fs)
    elseif ftype === :bs
        length(cutoff) < 2 && throw(ArgumentError("For band-stop filter two frequencies must be given."))
        responsetype = Bandstop(cutoff[1], cutoff[2]; fs=fs)
    end

    fprototype === :butterworth && (prototype = Butterworth(order))
    if fprototype === :chebyshev1
        rs == nothing && throw(ArgumentError("For Chebyshev1 filter rs must be given."))
        prototype = Chebyshev1(order, rs)
    end
    if fprototype === :chebyshev2
        rp == nothing && throw(ArgumentError("For Chebyshev2 filter rp must be given."))
        prototype = Chebyshev2(order, rp)
    end
    if fprototype === :elliptic
        rs == nothing && throw(ArgumentError("For Elliptic filter rs must be given."))
        rp == nothing && throw(ArgumentError("For Elliptic filter rp must be given."))
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
              title="Frequency response\nfilter: $(titlecase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order",
              xlims=(0, x_max),
              ylabel="Magnitude [dB]",
              xlabel="Frequency [Hz]",
              label="")
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
              title="Phase response\nfilter: $(titlecase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order",
              ylims=(-180, 180),
              xlims=(0, x_max),
              ylabel="Phase [°]",
              xlabel="Frequency [Hz]",
              label="")
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
              title="Group delay\nfilter: $(titlecase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order", xlims=(0, x_max),
              ylabel="Group delay [samples]",
              xlabel="Frequency [Hz]",
              label="")
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
        isfile(figure) && @warn "File $figure cannot be saved."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_avg(t, signal; labels=[""], norm=true, xlabel"Time [s]", ylabel="Chann title="Averaged signal and 95% CI plot", ylim=nothing)

Plots averaged `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange`
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `ylim::Union{Int64, Float64, Nothing}` - y-axis limits (-ylim:ylim)
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_avg(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Matrix{Float64}; norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Averaged signal and 95% CI plot", ylim::Union{Int64, Float64, Nothing}=nothing, kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))


    if norm == true
        s_normalized = signal_normalize_zscore(signal)
    else
        s_normalized = signal
    end
    
    s_normalized_m, s_normalized_s, s_normalized_u, s_normalized_l = signal_ci95(s_normalized)

    if ylim === nothing
        ylim = maximum(s_normalized_u)
        ylim = ceil(Int64, ylim)
    end

    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
             ylims=(-ylim, ylim),
             title=title,
             palette=:darktest;
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
    eeg_plot_avg(eeg; epoch=1, channel=nothing, offset=0, len=nothing, labels=[""], norm=true, xlabel="Time  ylabel="Channels", title="Averaged signal and 95% CI plot", figure="", kwargs...)

Plots averaged `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Int64` - epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}` - channel to display
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Nothing}` - displayed segment length in samples, default 10 seconds
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `ylim::Union{Int64, Float64, Nothing}` - y-axis limits (-ylim:ylim)
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_avg(eeg::EEG; epoch::Int64=1, channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, offset::Int64=0, len::Union{Int64, Nothing}=nothing, norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Averaged signal and 95% CI plot", ylim::Union{Int64, Float64, Nothing}=nothing, figure::String="", kwargs...)
    offset < 0 && throw(ArgumentError("Offset must be ≥ 0."))
    (len !== nothing && len <= 0) && throw(ArgumentError("len must be > 0."))
    len === nothing && (len = 10 * eeg_sr(eeg))
    (epoch < 1 || epoch > eeg.eeg_header[:epoch_n]) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg.eeg_header[:epoch_n])."))
    (typeof(channel) == Int64 && (channel < 1 || channel > eeg.eeg_header[:channel_n])) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg.eeg_header[:chansnel_n])."))
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    if typeof(channel) == Vector{Int64}
        sort!(channel)
        for idx in 1:length(channel)
            (channel[idx] < 1 || channel[idx] > eeg.eeg_header[:channel_n]) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
        end
    end

    # select channels, default is 1:20 or all channels
    if channel === nothing
        if eeg.eeg_header[:channel_n] >= 20
            channel = 1:20
        else
            channel = 1:eeg.eeg_header[:channel_n]
        end
    end

    # get epochs markers for len > epoch_len
    if (((len / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds] && eeg.eeg_header[:epoch_n] > 1)
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
        eeg_tmp = eeg_epochs(eeg_tmp, epoch_n=1)
        epoch_len = size(eeg.eeg_signals, 2)
        epoch_n = size(eeg.eeg_signals, 3)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
    else
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
    end

    signal = eeg_tmp.eeg_signals[:, :, epoch]
    labels = eeg_tmp.eeg_header[:labels]

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]

    (offset < 0 || offset > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))
    (offset + len > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset + len must be ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))

    signal = eeg_tmp.eeg_signals[:, (1 + offset):(offset + length(t)), epoch]

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
    if (((len / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds] && eeg.eeg_header[:epoch_n] > 1)
        p = vline!(p,
                   epoch_markers,
                   timeseries=:vline,
                   linestyle=:dash,
                   linewidth=0.5,
                   linecolor=:grey,
                   label="")
    end

    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure cannot be saved."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_butterfly(t, signal; labels=[""], norm=true, xlabel"Time [s]", ylabel="Chann title="Butterfly plot", ylim=nothing, kwargs...)

Butterfly plot of `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange`
- `signal::Matrix{Float64}`
- `labels::Vector{String}` - channel labels vector
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `ylim::Union{Int64, Float64, Nothing}` - y-axis limits (-ylim:ylim)
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_butterfly(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Butterfly plot", ylim::Union{Int64, Float64, Nothing}=nothing, kwargs...)

    typeof(t) <: AbstractRange && (t = float(collect(t)))

    channel_n = size(signal, 1)

    if norm == true
        s_normalized = signal_normalize_zscore(reshape(signal, size(signal, 1), size(signal, 2), 1))
    else
        s_normalized = signal
    end

    if ylim === nothing
        ylim = maximum(s_normalized) * 1.5
        ylim = ceil(Int64, ylim)
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
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
             ylims=(-ylim, ylim),
             title=title,
             palette=:darktest;
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
    eeg_plot_butterfly(eeg; epoch=1, channel=nothing, offset=0, len=nothing, labels=[""], norm=true, xlabel="Time  ylabel="Channels", title="Butterfly plot", ylim=nothing, head=true, figure="", kwargs...)

Butterfly plot of `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Int64` - epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}` - channel to display
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Nothing}` - displayed segment length in samples, default 10 seconds
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
function eeg_plot_butterfly(eeg::EEG; epoch::Int64=1, channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, offset::Int64=0, len::Union{Int64, Nothing}=nothing, labels::Vector{String}=[""], norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Butterfly plot", ylim::Union{Int64, Float64, Nothing}=nothing, head::Bool=true, hist::Bool=true, figure::String="", kwargs...)
    offset < 0 && throw(ArgumentError("Offset must be ≥ 0."))
    (len !== nothing && len <= 0) && throw(ArgumentError("len must be > 0."))
    len === nothing && (len = 10 * eeg_sr(eeg))
    (epoch < 1 || epoch > eeg.eeg_header[:epoch_n]) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg.eeg_header[:epoch_n])."))
    (typeof(channel) == Int64 && (channel < 1 || channel > eeg.eeg_header[:channel_n])) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg.eeg_header[:chansnel_n])."))
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    if typeof(channel) == Vector{Int64}
        sort!(channel)
        for idx in 1:length(channel)
            (channel[idx] < 1 || channel[idx] > eeg.eeg_header[:channel_n]) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
        end
    end

    # select channels, default is 1:20 or all channels
    if channel === nothing
        if eeg.eeg_header[:channel_n] >= 20
            channel = 1:20
        else
            channel = 1:eeg.eeg_header[:channel_n]
        end
    end

    # get epochs markers for len > epoch_len
    if (((len / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds] && eeg.eeg_header[:epoch_n] > 1)
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
        eeg_tmp = eeg_epochs(eeg_tmp, epoch_n=1)
        epoch_len = size(eeg.eeg_signals, 2)
        epoch_n = size(eeg.eeg_signals, 3)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
    else
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
    end

    signal = eeg_tmp.eeg_signals[:, :, epoch]
    labels = eeg_tmp.eeg_header[:labels]

    t = collect(0:(1 / eeg_sr(eeg_tmp)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg_tmp))
    t = t[1:(end - 1)]

    (offset < 0 || offset > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))
   (offset + len > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset + len must be ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))

    signal = eeg_tmp.eeg_signals[:, (1 + offset):(offset + length(t)), epoch]

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
    if (((len / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds] && eeg.eeg_header[:epoch_n] > 1)
        p = vline!(p,
                   epoch_markers,
                   timeseries=:vline,
                   linestyle=:dash,
                   linewidth=0.5,
                   linecolor=:grey,
                   label="")
    end

    # cannot plot electrodes without locations
    eeg.eeg_header[:channel_locations] == false && (head = false)

    e = plot(grid=false, framestyle=:none)
    head == true && (hd = eeg_plot_electrodes(eeg, labels=false, selected=channel, small=true))
    head == false && (hd = e)
    hist == true && (ht = signal_plot_histogram(signal, type=:hist, labels=[""]))
    hist == false && (ht = e)
    l = @layout [a b{0.1w}; c{0.9h} d]
    (head == true || hist == true) && (p = plot!(e, hd, p, ht, layout=l, left_margin=20px, bottom_margin=10px))
    head == false && hist == false && (p = plot!(p))

    if figure !== ""
        isfile(figure) && @warn "File $figure cannot be saved."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_psd(s_powers, s_freqs; frq_lim=nothing, xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="PSD", kwargs...)

Plots power spectrum density.

# Arguments

- `s_powers::Vector{Float64}`
- `s_freqs::Vector{Float64}`
- `frq_lim::Union{Int64, Float64, Nothing}` - x-axis limit
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_psd(s_powers::Vector{Float64}, s_freqs::Vector{Float64}; frq_lim::Union{Int64, Float64, Nothing}=nothing, xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="PSD", kwargs...)

    frq_lim === nothing && (frq_lim = s_freqs[end])

    p = plot(s_freqs,
             s_powers,
             xlabel=xlabel,
             ylabel=ylabel,
             xlims=(0, frq_lim),
             legend=false,
             t=:line,
             c=:black,
             title=title,
             palette=:darktest;
             kwargs...)

    plot(p)

    return p
end

"""
    signal_plot_psd(signal; fs, norm=false, frq_lim=nothing, xlabel="Frequency [Hz]", ylabel="Power [μV^2/ title="PSD", kwargs...)

Plots power spectrum density.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64` - sampling frequency
- `norm::Bool` - normalize the `signal` prior to calculations
- `frq_lim::Union{Int64, Float64, Nothing}` - x-axis limit
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_psd(signal::Vector{Float64}; fs::Int64, norm::Bool=false, frq_lim::Union{Int64, Float64, Nothing}=nothing, xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="PSD", kwargs...)
    fs < 0 && throw(ArgumentError("Sampling rate must be ≥ 0."))
    frq_lim === nothing || (frq_lim < 0 && throw(ArgumentError("Frequency limit rate must be ≥ 0.")))

    s_freqs, s_powers = signal_psd(signal, fs=fs, norm=norm)
    norm == true && (ylabel::String="Power [dB]")

    frq_lim === nothing && (frq_lim = s_freqs[end])

    p = plot(s_freqs,
             s_powers,
             xlabel=xlabel,
             ylabel=ylabel,
             xlims=(0, frq_lim),
             legend=false,
             t=:line,
             c=:black,
             title=title,
             palette=:darktest)

    plot(p)

    return p
end

"""
    signal_plot_psd(signal; fs, norm=false, average=false, frq_lim=nothing, labels=[""], xlabel="Frequency [ ylabel="Power [μV^2/Hz]", title="PSD", kwargs...)

Plots power spectrum density.

# Arguments

- `signal::Matrix{Float64}`
- `fs::Int64` - sampling rate
- `norm::Bool` - normalize the `signal` prior to calculations
- `average::Bool` - plots average power and 95%CI for all channels
- `frq_lim::Union{Int64, Float64, Nothing}` - x-axis limit
- `labels::Vector{String}` - channel labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false, average::Bool=false, frq_lim::Union{Int64, Float64, Nothing}=nothing, labels::Vector{String}=[""], xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="PSD", kwargs...)
    fs < 0 && throw(ArgumentError("Sampling rate must be ≥ 0."))
    frq_lim === nothing || (frq_lim < 0 && throw(ArgumentError("Frequency limit rate must be ≥ 0.")))

    norm == true && (ylabel="Power [dB]")

    channel_n = size(signal, 1)
    signal = reshape(signal, size(signal, 1), size(signal, 2), 1)
    s_powers, s_freqs = signal_psd(signal, fs=fs, norm=norm)
    s_powers = s_powers[:, :, 1]
    s_freqs = s_freqs[:, :, 1]

    frq_lim === nothing && (frq_lim = s_freqs[1, end])
    if average == true
        s_powers_m, s_powers_s, s_powers_u, s_powers_l = signal_ci95(s_powers)
        s_freqs = s_freqs[1, :]
        channel_n = 1
        labels == [""]
        title = "Average PSD with 95%CI"
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
             xlims=(0, frq_lim),
             title=title,
             palette=:darktest;
             kwargs...)
    if average == true
    p = plot!(s_freqs,
              s_powers_u,
              fillrange = s_powers_l,
              fillalpha = 0.35,
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
    eeg_plot_psd(eeg; epoch=1, channel=nothing, offset=0, len=nothing, labels=[""], norm=false, average=false, frq_lim=not xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="PSD", head=false, figure="", kwargs...)

Plots power spectrum density.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Int64` - epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}` - channel to display
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Nothing}` - displayed segment length in samples, default 10 seconds
- `labels::Vector{String}` - channel labels vector
- `norm::Bool` - normalize the `signal` prior to calculations
- `average::Bool` - plots average power and 95%CI for all channels
- `frq_lim::Union{Int64, Float64, Nothing}` - x-axis limit
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `head::Bool` - add head with electrodes
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_psd(eeg::EEG; epoch::Int64=1, channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, offset::Int64=0, len::Union{Int64, Nothing}=nothing, labels::Vector{String}=[""], norm::Bool=false, average::Bool=false, frq_lim::Union{Int64, Float64, Nothing}=nothing, xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="PSD", head::Bool=false, figure::String="", kwargs...)
    frq_lim === nothing || (frq_lim < 0 && throw(ArgumentError("frq_lim must be ≥ 0 Hz.")))
    (epoch < 1 || epoch > eeg.eeg_header[:epoch_n]) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg.eeg_header[:epoch_n])."))
    (typeof(channel) == Int64 && (channel < 1 || channel > eeg.eeg_header[:channel_n])) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg.eeg_header[:chansnel_n])."))
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    if typeof(channel) == Vector{Int64}
        sort!(channel)
        for idx in 1:length(channel)
            (channel[idx] < 1 || channel[idx] > eeg.eeg_header[:channel_n]) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
        end
    end
    (len !== nothing && len <= 0) && throw(ArgumentError("len must be > 0."))
    len === nothing && (len = 10 * eeg_sr(eeg))

    # select channels, default is all channels
    channel === nothing && (channel = 1:eeg.eeg_header[:channel_n])
    if (((len / eeg_sr(eeg)) / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds]
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
        eeg_tmp = eeg_epochs(eeg_tmp, epoch_n=1)
    else
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
    end

    (offset < 0 || offset > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))
    (offset + len > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset + len must be ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))

    fs = eeg_sr(eeg_tmp)
    signal = eeg_tmp.eeg_signals[:, (1 + offset):(offset + len), epoch]
    labels = eeg_tmp.eeg_header[:labels]

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

    # cannot plot electrodes without locations
    eeg_tmp.eeg_header[:channel_locations] == false && (head = false)

    if head == true
        p = plot(p)
        h = eeg_plot_electrodes(eeg, labels=false, selected=channel, small=true)
        e = plot(grid=false, framestyle=:none)
        l=@layout @layout [a b{0.1w}; c{0.9h}]
        p = plot(e, h, p, layout=l, left_margin=20px, bottom_margin=10px)
        plot(p)
    else
        plot(p)
    end

    if figure !== ""
        isfile(figure) && @warn "File $figure cannot be saved."
        savefig(p, figure)
    end

    return p
end

"""
    eeg_plot_electrodes(eeg; channel=nothing, selected-nothing, labels=true, head=true, head_labels=false, small=false)

Plots electrodes.

# Arguments

- `eeg:EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}` - channel to display
- `selected::Union{Int64, Vector{Int64}, AbstractRange, Nothing}` - which channel should be highlighted
- `labels::Bool` - plot electrode labels
- `head::Bool` - plot head
- `head_labels::Bool` - plot head labels
- `small::Bool` - draws small plot

# Returns

- `p::Plot`
"""
function eeg_plot_electrodes(eeg::EEG; channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, selected::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, labels::Bool=true, head::Bool=true, head_labels::Bool=false, small::Bool=false)
    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available."))
    (typeof(channel) == Int64 && (channel < 1 || channel > eeg.eeg_header[:channel_n])) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    if typeof(channel) == Vector{Int64}
        sort!(channel)
        for idx in 1:length(channel)
            (channel[idx] < 1 || channel[idx] > eeg.eeg_header[:channel_n]) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
        end
    end
    (typeof(selected) == Int64 && (selected < 1 || selected > eeg.eeg_header[:channel_n])) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
    if typeof(selected) == Vector{Int64}
        sort!(selected)
        for idx in 1:length(selected)
            (selected[idx] < 1 || selected[idx] > eeg.eeg_header[:channel_n]) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
        end
    end

    # select channels, default is all channels
    selected === nothing && (selected = 1:eeg.eeg_header[:channel_n])
    channel === nothing && (channel = 1:eeg.eeg_header[:channel_n])
    eeg_tmp = eeg_keep_channel(eeg, channel=channel)

    loc_x = eeg_tmp.eeg_header[:xlocs]
    loc_y = eeg_tmp.eeg_header[:ylocs]
    x_lim = (findmin(loc_x)[1] * 1.8, findmax(loc_x)[1] * 1.8)
    y_lim = (findmin(loc_y)[1] * 1.8, findmax(loc_y)[1] * 1.8)
    if small == true
        plot_size = (400, 400)
        marker_size = 4
        labels = false
    else
        plot_size = (800, 800)
        marker_size = 8
        font_size = 8
    end

    p = plot(grid=false, framestyle=:none, palette=:darktest, size=plot_size, markerstrokewidth=0, border=:none, margins=0px, aspect_ratio=1)
    head == true && eeg_draw_head(p, loc_x, loc_x, head_labels=head_labels)
    if length(selected) >= eeg_tmp.eeg_header[:channel_n]
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
        for idx in 1:length(eeg_tmp.eeg_header[:labels])
        plot!(p, annotation=(loc_x[idx] + 0.00, loc_y[idx] + 0.04, text(eeg_tmp.eeg_header[:labels][idx], pointsize=font_size)))
        end
        p = plot!()
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
    (epoch < 1 || epoch > eeg.eeg_header[:epoch_n]) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg.eeg_header[:epoch_n])."))

    labels = eeg_labels(eeg)
    channel_n = size(m, 1)
    ndims(m) == 3 && (m = m[:, :, epoch])

    p = heatmap(m, xticks=(1:channel_n, labels), yticks=(1:channel_n, eeg_labels(eeg)); kwargs...)
    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure cannot be saved."
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
function eeg_plot_covmatrix(eeg::EEG, cov_m::Union{Matrix{Float64}, Array{Float64, 3}}, lags::Union{Vector{Int64}, Vector{Float64}}; channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, epoch::Int64=1, figure::String="", kwargs...)
    (epoch < 1 || epoch > eeg.eeg_header[:epoch_n]) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg.eeg_header[:epoch_n])."))
    (typeof(channel) == Int64 && (channel < 1 || channel > eeg.eeg_header[:channel_n])) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg.eeg_header[:chansnel_n])."))
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    if typeof(channel) == Vector{Int64}
        sort!(channel)
        for idx in 1:length(channel)
            (channel[idx] < 1 || channel[idx] > eeg.eeg_header[:channel_n]) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
        end
    end

    # select channels
    channel === nothing && (channel = 1:eeg.eeg_header[:channel_n])

    labels = eeg_labels(eeg)
    ndims(cov_m) == 3 && (cov_m = cov_m[:, :, epoch])
    p = []
    for idx in channel
        push!(p, plot(lags, cov_m[idx, :], title="ch: $(labels[idx])", label="", titlefontsize=6, xtickfontsize=4, ytickfontsize=4, lw=0.5))
    end
    p = plot(p...; kwargs...)

    if figure !== ""
        isfile(figure) && @warn "File $figure cannot be saved."
        savefig(p, figure)
    end

    return p
end

"""
    signal_plot_spectrogram(signal; fs, offset=0, norm=true, demean=true, ylim=nothing, xlabel="Time  ylabel="Frequency [Hz]", title="Spectrogram", kwargs...)

Plots spectrogram of `signal`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64` - sampling frequency
- `offset::Int64` - displayed segment offset in samples
- `norm::Bool` - normalize powers to dB
- `demean::Bool` - demean signal prior to analysis
- `ylim::Union{Int64, Float64, Nothing}` - y-axis limits (0:ylim) Hz
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function signal_plot_spectrogram(signal::Vector{Float64}; fs::Int64, offset::Int64=0, norm::Bool=true, demean::Bool=true, ylim::Union{Int64, Float64, Nothing}=nothing, xlabel="Time [s]", ylabel="Frequency [Hz]", title="Spectrogram", kwargs...)
    fs < 1 && throw(ArgumentError("Sampling rate must be ≥ 1 Hz."))
    ylim > fs / 2 && throw(ArgumentError("ylim must be smaller than Nyquist frequency ($(fs/2) Hz)."))

    demean == true && (signal = signal_demean(signal))

    nfft = length(signal)
    interval = fs
    overlap = round(Int64, fs * 0.85)

    ylim === nothing && (ylim = fs/2)
    spec = spectrogram(signal, interval, overlap, nfft=nfft, fs=fs, window=hanning)
    t = collect(spec.time) .+ (offset / fs)

    if norm == false
        p = heatmap(t,
                    spec.freq,
                    spec.power,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    ylim=(0, ylim),
                    xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
                    title=title,
                    colorbar_title="Power/frequency [μV^2/Hz]";
                    kwargs...)
    else
        # in dB
        p = heatmap(t,
                    spec.freq,
                    pow2db.(spec.power),
                    xlabel=xlabel,
                    ylabel=ylabel,
                    ylim=(0, ylim),
                    xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
                    title=title,
                    colorbar_title="Power/frequency [dB/Hz]";
                    kwargs...)
    end

    return p
end

"""
    eeg_plot_spectrogram(eeg; channel, epoch=1, offset=0, len=nothing, norm=true, title="Signal plot", ylim=nothing, kwargs...)

Plots spectrogram of `eeg` channel.

# Arguments

- `eeg:EEG`
- `channel::Int64` - channel to plot
- `epoch::Int64` - epoch to plot
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Nothing}` - displayed segment length in samples, default 10 seconds
- `norm::Bool` - normalize powers to dB
- `ylim::Union{Int64, Float64, Nothing}` - y-axis limits (0:ylim) Hz
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `ylim::Union{Int64, Float64, Nothing}` - y-axis limits (0:ylim)
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_spectrogram(eeg::EEG; channel::Int64, epoch::Int64=1, offset::Int64=0, len::Union{Int64, Nothing}=nothing, norm::Bool=true, ylim::Union{Int64, Float64, Nothing}=nothing, xlabel="Time [s]", ylabel="Frequency [Hz]", title="Spectrogram", figure="", kwargs...)
    offset < 0 && throw(ArgumentError("Offset must be ≥ 0."))
    (len !== nothing && len <= 0) && throw(ArgumentError("len must be > 0."))
    len === nothing && (len = 10 * eeg_sr(eeg))
    (epoch < 1 || epoch > eeg.eeg_header[:epoch_n]) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg.eeg_header[:epoch_n])."))
    (channel < 1 || channel > eeg.eeg_header[:channel_n]) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))

    fs = eeg_sr(eeg)
    ylim === nothing && (ylim = fs/2)

    # get epochs markers for len > epoch_len
    if len + offset > eeg.eeg_header[:epoch_duration_samples]
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = size(eeg.eeg_signals, 2)
        epoch_n = size(eeg.eeg_signals, 3)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
    else
        eeg_tmp = eeg
    end

    (offset < 0 || offset > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))
    (offset + len > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset + len must be ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))

    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + len), epoch]

    p = signal_plot_spectrogram(signal,
                                fs=fs,
                                offset=offset,
                                norm=norm,
                                xlabel=xlabel,
                                ylabel=ylabel,
                                ylim=ylim,
                                title=title;
                                kwargs...)

    # add epochs markers
    if (len + offset) > eeg.eeg_header[:epoch_duration_samples]
        p = vline!(p,
                   epoch_markers[floor(Int64, 1 + (offset / epoch_len)):floor(Int64, (len + offset) / epoch_len)],
                   timeseries=:vline,
                   linestyle=:dash,
                   linewidth=0.5,
                   linecolor=:black,
                   label="")
    end

    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure cannot be saved."
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
             label=label,
             xlims=(floor(t[1]), ceil(t[end])),
             xticks=floor(t[1]):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end]),
             ylims=(-ylim, ylim),
             title=title,
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
function signal_plot_ica(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, ica::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="ICA", kwargs...)

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
             palette=:darktest;
             kwargs...)
    for idx in 1:ica_n
        p = plot!(t,
                  ica_normalized[idx, 1:length(t)],
                  label="", color=ica_color[idx])
    end
    p = plot!(p, yticks=((ica_n - 1):-1:0, labels))

    return p
end

"""
    eeg_plot_ica(eeg; epoch=1, ic=nothing, offset=0, len=nothing, norm=true, xlabel="Time  ylabel="", title="ICA", figure="", kwargs...)

Plots ICs.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Int64` - epoch number to display
- `ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}` - which IC to plot
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Nothing}` - displayed segment length in samples, default 10 seconds
- `norm::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_ica(eeg::EEG; epoch::Int64=1, ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, offset::Int64=0, len::Union{Int64, Nothing}=nothing, norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="ICA", figure::String="", kwargs...)
    typeof(ic) <: AbstractRange && (ic = collect(ic))
    offset < 0 && throw(ArgumentError("offset must be ≥ 0."))
    (len !== nothing && len <= 0) && throw(ArgumentError("len must be > 0."))
    (typeof(ic) == Int64 && ic <= 0) && throw(ArgumentError("ic must be > 0."))
    len === nothing && (len = 10 * eeg_sr(eeg))
    len > eeg.eeg_header[:eeg_duration_samples] && throw(ArgumentError("len must be < $(eeg.eeg_header[:eeg_duration_samples])."))
    :ica in eeg.eeg_header[:components] || throw(ArgumentError("EEG does not contain ICA. Perform eeg_ica!(EEG) first."))
    (epoch < 1 || epoch > eeg.eeg_header[:epoch_n]) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg.eeg_header[:epoch_n])."))

    # select ICs, default is all
    if typeof(ic) == Vector{Int64}
        sort(ic)
        for idx in 1:length(ic)
            ic[idx] > size(eeg.eeg_components[ica_idx], 1) && throw(ArgumentError("ic must be ≥ 1 and ≤ $(size(eeg.eeg_components[ica_idx], 1))."))
        end
    end
    ica_idx = findfirst(isequal(:ica), eeg.eeg_header[:components])
    (typeof(ic) == Int64 && ic > size(eeg.eeg_components[ica_idx], 1)) && throw(ArgumentError("ic must be ≥ 1 and ≤ $(size(eeg.eeg_components[ica_idx], 1))."))
    ic === nothing && (ic = 1:size(eeg.eeg_components[ica_idx], 1))

    # get epochs markers for len > epoch_len
    if (((len / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds] && eeg.eeg_header[:epoch_n] > 1)
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = size(eeg.eeg_signals, 2)
        epoch_n = size(eeg.eeg_signals, 3)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
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

    (offset < 0 || offset > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))
    (offset + len > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset + len must be ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))

    ica_m = eeg_tmp.eeg_components[ica_idx][ic, (1 + offset):(offset + Int(len)), epoch]

    p = signal_plot_ica(t,
                        ica_m,
                        norm=norm,
                        labels=labels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        title=title;
                        kwargs...)

    # add epochs markers
    if (((len / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds] && eeg.eeg_header[:epoch_n] > 1)
        p = vline!(p,
                   epoch_markers,
                   timeseries=:vline,
                   linestyle=:dash,
                   linewidth=0.5,
                   linecolor=:grey,
                   label="")
    end

    plot(p)

    if figure !== ""
        isfile(figure) && @warn "File $figure cannot be saved."
        savefig(p, figure)
    end

    return p
end

"""
    eeg_topoplot(eeg)

Plots topographical view of `eeg` signals.

# Arguments

- `eeg::EEG`

# Returns

- `plot`
"""
function eeg_topoplot(eeg::EEG; c::Symbol=:amplitude, c_idx::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, t::Int64)
    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available."))
    c in [:amplitude] || throw(ArgumentError("c must be: :amplitude."))
    # c in [:amplitude, :psd, :ica] || throw(ArgumentError("c must be: :amplitude, :psd or :ica."))
    (c !== :amplitude && c in eeg.eeg_header[:components]) || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view available components."))
    
    c === :amplitude && (s_non_interpolated = @view eeg.eeg_signals[:, t, epoch])

    # plot signal at electrodes at time
    loc_x = eeg.eeg_header[:xlocs]
    loc_y = eeg.eeg_header[:ylocs]
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
    time_point_idx = eeg.eeg_time[t]
    s_interpolated = zeros(interpolation_factor, interpolation_factor)
    electrode_locations = [loc_x loc_y]'
    itp = ScatteredInterpolation.interpolate(Shepard(), electrode_locations, s_non_interpolated)
    for idx1 in 1:interpolation_factor
        for idx2 in 1:interpolation_factor
            s_interpolated[idx1, idx2] = ScatteredInterpolation.evaluate(itp, [interpolation_m[idx1, idx2][1]; interpolation_m[idx1, idx2][2]])[1]
        end
    end

    p = plot(grid=false, framestyle=:none, border=:none, margins=0px, aspect_ratio=1)
    p = plot!(interpolated_x, interpolated_y, s_interpolated, fill=:lighttest, seriestype=:contourf)
    p = plot!((loc_x, loc_y), color=:black, seriestype=:scatter, xlims=x_lim, ylims=x_lim, grid=true, label="", markersize=4, markerstrokewidth=0, markerstrokealpha=0)
    h = eeg_draw_head(p, loc_x, loc_y, head_labels=false; lw=2)
    p = plot!(h)

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

    if type === :hist
        p = plot(signal,
                 seriestype=:hist,
                 xlabel=xlabel,
                 ylabel=ylabel,
                 label=label,
                 title=title,
                 palette=:darktest,
                 xticks=nothing,
                 yticks=nothing,
                 grid=false;
                 kwargs...)
    else
        p = plot(signal,
                 seriestype=:density,
                 xlabel=xlabel,
                 ylabel=ylabel,
                 label=label,
                 title=title,
                 palette=:darktest,
                 xticks=nothing,
                 yticks=nothing,
                 grid=false;
                 kwargs...)
    end

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
function signal_plot_histogram(signal::Matrix{Float64}; type::Symbol=:hist, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", kwargs...)

    channel_n = size(signal, 1)

    # reverse so 1st channel is on top
    signal = reverse(signal[:, :], dims = 1)

    if labels == [""]
        labels = Vector{String}(undef, size(signal, 1))
        for idx in 1:size(signal, 1)
            labels[idx] = ""
        end
    end

    # plot channels
    if channel_n == 1
        p = signal_plot_histogram(signal[1, :],
                                  type=type,
                                  label=labels[1],
                                  xlabel=xlabel,
                                  ylabel=ylabel,
                                  linecolor=1,
                                  fillcolor=1,
                                  linewidth=0.5,
                                  top_margin=0px,
                                  bottom_margin=0px,
                                  orientation=:h,
                                  palette=:darktest,
                                  yticks=true,
                                  xaxis=false,
                                  kwargs...)
    else
        p = []
        push!(p, plot(framestyle=:none, margins=0px))
        for idx in 1:channel_n
            if type === :hist
                push!(p, signal_plot_histogram(signal[idx, :],
                                               type=type,
                                               label=labels[idx],
                                               xlabel=xlabel,
                                               ylabel=ylabel,
                                               linecolor=idx,
                                               fillcolor=idx,
                                               linewidth=0.5,
                                               xticks=(range(floor(minimum(signal[idx, :]), digits=1), ceil(maximum(signal[idx, :]), digits=1), length=2)),
                                               yaxis=false))
            else
                push!(p, signal_plot_histogram(signal[idx, :],
                                               type=type,
                                               label=labels[idx],
                                               xlabel=xlabel,
                                               ylabel=ylabel,
                                               linecolor=idx,
                                               linewidth=0.5,
                                               xticks=(range(floor(minimum(signal[idx, :]), digits=1), ceil(maximum(signal[idx, :]), digits=1), length=2)),
                                               yaxis=false))
            end
        end
        push!(p, plot(framestyle=:none, margins=0px))
        p = plot(p...,
                 top_margin=0px,
                 bottom_margin=0px,
                 palette=:darktest,
                 layout=(channel_n+2, 1);
                 kwargs...)
    end
    plot(p)

    return p
end

"""
    eeg_plot_histogram(eeg; epoch=1, channel=nothing, offset=0, len=nothing, labels=[""], xlabel="", ylabel="", title="", figure="", kwargs...)

Plots `eeg` channels histograms.

# Arguments

- `eeg::EEG` - EEG object
- `type::Symbol[:hist, :kd]` - type of histogram: regular `:hist` or kernel density `:kd`
- `epoch::Int64` - epoch number to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}` - channel to display
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Nothing}` - displayed segment length in samples, default 10 seconds
- `labels::Vector{String}` - channel labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `figure::String` - name of the output figure file
- `kwargs` - other arguments for plot() function

# Returns

- `p::Plot`
"""
function eeg_plot_histogram(eeg::EEG; type::Symbol=:hist, epoch::Int64=1, channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, offset::Int64=0, len::Union{Int64, Nothing}=nothing, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", figure::String="", kwargs...)
    offset < 0 && throw(ArgumentError("Offset must be ≥ 0."))
    (len !== nothing && len <= 0) && throw(ArgumentError("len must be > 0."))
    len === nothing && (len = 10 * eeg_sr(eeg))
    (epoch < 1 || epoch > eeg.eeg_header[:epoch_n]) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg.eeg_header[:epoch_n])."))
    (typeof(channel) == Int64 && (channel < 1 || channel > eeg.eeg_header[:channel_n])) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg.eeg_header[:chansnel_n])."))
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    if typeof(channel) == Vector{Int64}
        sort!(channel)
        for idx in 1:length(channel)
            (channel[idx] < 1 || channel[idx] > eeg.eeg_header[:channel_n]) && throw(ArgumentError("selected must be ≥ 1 and ≤ $(eeg.eeg_header[:channel_n])."))
        end
    end

    # select channels, default is 1:20 or all channels
    if channel === nothing
        if eeg.eeg_header[:channel_n] >= 20
            channel = 1:20
        else
            channel = 1:eeg.eeg_header[:channel_n]
        end
    end

    # get epochs markers for len > epoch_len
    if (((len / eeg_sr(eeg)) + (offset / eeg_sr(eeg))) > eeg.eeg_header[:epoch_duration_seconds] && eeg.eeg_header[:epoch_n] > 1)
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
        eeg_tmp = eeg_epochs(eeg_tmp, epoch_n=1)
    else
        eeg_tmp = eeg_keep_channel(eeg, channel=channel)
    end

    labels = eeg_tmp.eeg_header[:labels]

    (offset < 0 || offset > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))
    (offset + len > eeg_tmp.eeg_header[:epoch_duration_samples]) && throw(ArgumentError("offset + len must be ≤ $(eeg_tmp.eeg_header[:epoch_duration_samples])."))

    signal = eeg_tmp.eeg_signals[:, (1 + offset):(offset + len), epoch]

    p = signal_plot_histogram(signal,
                              type=type,
                              labels=labels,
                              xlabel=xlabel,
                              ylabel=ylabel,
                              title=title;
                              kwargs...)

    if figure !== ""
        isfile(figure) && @warn "File $figure cannot be saved."
        savefig(p, figure)
    end

    return p
end