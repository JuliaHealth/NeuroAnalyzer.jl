"""
    signal_plot(t, signal; offset=0, labels=[], normalize=true, xlabel="Time [s]", ylabel="Amplitude [μV]", title="Signal plot", yamp=nothing)

Plots `signal` against time vector `t`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}` - the time vector
- `signal::Vector{Float64}` - the signal vector
- `offset::Int64` - displayed segment offset in samples
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `yamp::Union{Nothing, Int64, Float64}` - y-axis limits (-yamp:yamp)
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}, signal::Vector{Float64}; offset::Int64=0, labels::Vector{String}=[], normalize::Bool=true, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Signal plot", yamp::Union{Nothing, Int64, Float64}=nothing)

    if typeof(t) == UnitRange{Int64} || typeof(t) == StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
        t = float(collect(t))
    end

    if yamp === nothing
        yamp = maximum(signal) * 1.5
        yamp = ceil(Int64, yamp)
    end

    p = plot(t, signal[1+offset:(offset + length(t))], xlabel=xlabel, ylabel=ylabel, label="", xlims=(ceil(t[1]), floor(t[end])), ylims=(-yamp, yamp), title=title)

    plot(p)

    return p
end

"""
    signal_plot(t, signal; offset=0, labels=[""], normalize=true, xlabel"Time [s]", ylabel="Channels", title="Signal plot")

Plots `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}`
- `signal::Matrix{Float64}`
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Float64}` - length in seconds
- `labels::Vector{String}` - labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}, signal::Matrix{Float64}; offset::Int64=0, len::Union{Int64, Float64}=10.0, labels::Vector{String}=[""], normalize::Bool=true, xlabel::String="Time [s]", ylabel::String="Channels", title::String="Signal plot")
    
    if typeof(t) == UnitRange{Int64} || typeof(t) == StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
        t = float(collect(t))
    end
    
    channels_no = size(signal, 1)

    # reverse so 1st channel is on top
    signal = reverse(signal[:, :], dims = 1)
    signal_normalized = zeros(size(signal))

    if normalize == true
        # normalize and shift so all channels are visible
        variances = var(signal, dims=2)
        mean_variance = mean(variances)
        for idx in 1:channels_no
            signal_normalized[idx, :] = (signal[idx, :] .- mean(signal[idx, :])) ./ mean_variance .+ (idx - 1)
        end
    else
        signal_normalized = signal
    end

    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(ceil(t[1]), floor(t[end])),
             ylims=(-0.5, channels_no-0.5),
             title=title)
    for idx in 1:channels_no
        p = plot!(t,
                  signal_normalized[idx, (1 + offset):(offset + length(t))],
                  label="")
    end
    p = plot!(p, yticks = ((channels_no - 1):-1:0, labels))

    return p
end

"""
    eeg_plot(eeg; t, epoch=1, channels=nothing, offset=0, len=10.0, labels=[""], normalize=true, xlabel="Time [s]", ylabel="Channels", title="Signal plot", figure="")

Plots `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `t::Union{Nothing, Vector{Float64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}` - the time vector
- `epoch::Int64` - epoch number to display
- `channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}` - channels to display
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Float64}` - length in seconds
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `figure::String` - name of the output figure file
"""
function eeg_plot(eeg::EEG; t::Union{Vector{Float64}, UnitRange{Int64}, Nothing}=nothing, epoch::Int64=1, channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}=nothing, offset::Int64=0, len::Union{Int64, Float64}=10.0, labels::Vector{String}=[""], normalize::Bool=true, xlabel::String="Time [s]", ylabel::String="Channels", title::String="Signal plot", figure::String="")

    if epoch < 1 || epoch > eeg.eeg_header[:epochs_no]
        throw(ArgumentError("Epoch index out of range."))
    end

    if typeof(t) == UnitRange{Int64}
        t = collect(t)
    end

    # select channels, default is 1:20 or all channels
    if channels === nothing
        if eeg.eeg_header[:channels_no] >= 20
            channels = 1:20
        else
            channels = 1:eeg.eeg_header[:channels_no]
        end
    end

    eeg_temp = eeg_keep_channel(eeg, channels)
    signal = eeg_temp.eeg_signals[:, :, epoch]
    labels = eeg_temp.eeg_header[:labels]

    # default time is 10 seconds or epoch_duration_seconds
    fs = eeg_samplingrate(eeg_temp)
    len > eeg_temp.eeg_header[:epoch_duration_seconds] && (len = eeg_temp.eeg_header[:epoch_duration_seconds])
    t === nothing && (t = (offset / fs) .+ collect(0:1/fs:len))
    t = t[1:(end - 1)]
    if offset < 0 || offset > eeg_temp.eeg_header[:epoch_duration_samples]
        throw(ArgumentError("Offset value out of range."))
    end
    if offset + (len * eeg_samplingrate(eeg_temp)) > eeg_temp.eeg_header[:epoch_duration_samples]
        throw(ArgumentError("Offset value or length value out of range."))
    end

    p = signal_plot(t,
                    signal,
                    offset=offset,
                    labels=labels,
                    normalize=normalize,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    title=title)

    plot(p)

    if figure !== ""
        try
            savefig(p, figure)
        catch error
            throw(ArgumentError("File $figure cannot be saved."))
            return false
        end
    end

    return p
end

"""
    eeg_draw_head(p, loc_x, loc_y, add_labels=true)

Draws head over a topographical plot `p`.

# Arguments

- `p::Plot` - topographical plot
- `loc_x::Vector{Float64}` - vector of x electrode position
- `loc_y::Vector{Float64}` - vector of y electrode position
- `add_labels::Bool` - add text labels to the plot
"""
function eeg_draw_head(p, loc_x::Vector{Float64}, loc_y::Vector{Float64}, add_labels::Bool=true)
    pts = Plots.partialcircle(0, 2π, 100, maximum(loc_x))
    x, y = Plots.unzip(pts)
    x = x .* 1.1
    y = y .* 1.1
    head = Shape(x, y)
    nose = Shape([(-0.1, maximum(y)), (0, maximum(y) + 0.1 * maximum(y)), (0.1, maximum(y))])
    ear_l = Shape([(minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), 0.1), (minimum(x), 0.1)])
    ear_r = Shape([(maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), 0.1), (maximum(x), 0.1)])
    plot!(p, head, fill=nothing, label="")
    plot!(p, nose, fill=nothing, label="")
    plot!(p, ear_l, fill=nothing, label="")
    plot!(p, ear_r, fill=nothing, label="")
    if add_labels == true
        plot!(p, annotation=(0, 1 - maximum(y) / 5, text("Inion", pointsize=12, halign=:center, valign=:center)))
        plot!(p, annotation=(0, -1 - minimum(y) / 5, text("Nasion", pointsize=12, halign=:center, valign=:center)))
        plot!(p, annotation=(-1 - minimum(x) / 5, 0, text("Left", pointsize=12, halign=:center, valign=:center, rotation=90)))
        plot!(p, annotation=(1 - maximum(x) / 5, 0, text("Right", pointsize=12, halign=:center, valign=:center, rotation=-90)))
    end
end

"""
    filter_response(fprototype, ftype, cutoff, fs, order, rp, rs, window, figure)

Returns zero phase distortion filter response.

# Arguments

- `fprototype::Symbol[:butterworth, :chebyshev1, :chebyshev2, :elliptic]
- `ftype::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `fs::Int64` - sampling rate
- `order::Int64` - filter order
- `rp::Union{Nothing, Int64, Float64}` - dB ripple in the passband
- `rs::Union{Nothing, Int64, Float64}` - dB attenuation in the stopband
- `window::Union{Nothing, Vector{Float64}} - window, required for FIR filter
- `figure::String` - name of the output figure file
"""
function filter_response(;fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}, fs::Int64, order::Int64, rp::Union{Nothing, Int64, Float64}=nothing, rs::Union{Nothing, Int64, Float64}=nothing, window::Union{Nothing, Vector{Float64}}=nothing, figure::String="")
    ftype in [:lp, :hp, :bp, :bs] || throw(ArgumentError("""Filter type must be ":bp", ":hp", ":bp" or ":bs"."""))
    fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic] || throw(ArgumentError("""Filter prototype must be ":butterworth", ":chebyshev1:, ":chebyshev2" or ":elliptic"."""))

    if ftype == :lp
        length(cutoff) > 1 && throw(ArgumentError("For low-pass filter one frequency must be given."))
        responsetype = Lowpass(cutoff; fs=fs)
    elseif ftype == :hp
        length(cutoff) > 1 && throw(ArgumentError("For high-pass filter one frequency must be given."))
        responsetype = Highpass(cutoff; fs=fs)
    elseif ftype == :bp
        responsetype = Bandpass(cutoff[1], cutoff[2]; fs=fs)
    elseif ftype == :bs
        length(cutoff) < 2 && throw(ArgumentError("For band-stop filter two frequencies must be given."))
        responsetype = Bandstop(cutoff[1], cutoff[2]; fs=fs)
    end

    fprototype == :butterworth && (prototype = Butterworth(order))
    if fprototype == :chebyshev1
        rs == nothing && throw(ArgumentError("For Chebyshev1 filter rs must be given."))
        prototype = Chebyshev1(order, rs)
    end
    if fprototype == :chebyshev2
        rp == nothing && throw(ArgumentError("For Chebyshev2 filter rp must be given."))
        prototype = Chebyshev2(order, rp)
    end
    if fprototype == :elliptic
        rs == nothing && throw(ArgumentError("For Elliptic filter rs must be given."))
        rp == nothing && throw(ArgumentError("For Elliptic filter rp must be given."))
        prototype = Elliptic(order, rp, rs)
    end

    ffilter = digitalfilter(responsetype, prototype)

    H, w = freqresp(ffilter)
    H = 20 * log10.(abs.(H))
    # convert to dB
    # convert rad/sample to Hz
    w = w .* fs / 2 / pi
    x_max = w[end]
    ftype == :hp && (x_max = cutoff * 10)
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
    ftype == :hp && (x_max = cutoff * 10)
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
    ftype == :hp && (x_max = cutoff * 10)
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

    p = plot(p1, p2, p3, layout=(3, 1))

    if figure !== ""
        try
            savefig(p, figure)
        catch error
            throw(ArgumentError("File $figure cannot be saved."))
            return false
        end
    end

    return p
end

"""
    signal_plot_avg(t, signal; offset=0, len=10.0, labels=[""], normalize=true, xlabel"Time [s]", ylabel="Channels", title="Averaged signal and 95% CI plot", yamp=nothing)

Plots averaged `signal` channels.
# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}`
- `signal::Matrix{Float64}`
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Float64}` - length in seconds
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `yamp::Union{Nothing, Int64, Float64}` - y-axis limits (-yamp:yamp)
"""
function signal_plot_avg(t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}, signal::Matrix{Float64}; offset::Int64=0, len::Union{Int64, Float64}=10.0, normalize::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Averaged signal and 95% CI plot", yamp::Union{Nothing, Int64, Float64}=nothing)
    
    if typeof(t) == UnitRange{Int64} || typeof(t) == StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
        t = float(collect(t))
    end

    if normalize == true
        signal_normalized = signal_normalize_mean(signal)
    else
        signal_normalized = signal
    end
    
    m, s, u, l = signal_ci95(signal_normalized)

    if yamp === nothing
        yamp = maximum(u)
        yamp = ceil(Int64, yamp)
    end

    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(ceil(t[1]),
                    floor(t[end])),
             ylims=(-yamp, yamp),
             title=title)
    p = plot!(t,
              m[(1 + offset):(offset + length(t))],
              label=false,
              t=:line,
              c=:black)
    p = plot!(t,
              u[(1 + offset):(offset + length(t))],
              label=false,
              t=:line,
              c=:grey,
              lw=0.5)
    p = plot!(t,
              l[(1 + offset):(offset + length(t))],
              label=false,
              t=:line,
              c=:grey,
              lw=0.5)

    return p
end

"""
    eeg_plot_avg(eeg; t, epoch=1, channels=nothing, offset=0, len=10.0, labels=[""], normalize=true, xlabel="Time [s]", ylabel="Channels", title="Averaged signal and 95% CI plot", figure="")

Plots averaged `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `t::Union{Nothing, Vector{Float64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}` - the time vector
- `epoch::Int64` - epoch number to display
- `channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}` - channels to display
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Float64}` - length in seconds
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `yamp::Union{Nothing, Int64, Float64}` - y-axis limits (-yamp:yamp)
- `figure::String` - name of the output figure file
"""
function eeg_plot_avg(eeg::EEG; t::Union{Vector{Float64}, UnitRange{Int64}, Nothing}=nothing, epoch::Int64=1, channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}=nothing, offset::Int64=0, len::Union{Int64, Float64}=10.0, normalize::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Averaged signal and 95% CI plot", yamp::Union{Nothing, Int64, Float64}=nothing, figure::String="")

    if epoch < 1 || epoch > eeg.eeg_header[:epochs_no]
        throw(ArgumentError("Epoch index out of range."))
    end

    if typeof(t) == UnitRange{Int64}
        t = collect(t)
    end

    # select channels, default is all channels
    channels === nothing && (channels = 1:eeg.eeg_header[:channels_no])
    eeg_temp = eeg_keep_channel(eeg, channels)

    # default time is 10 seconds or epoch_duration_seconds
    fs = eeg_samplingrate(eeg_temp)
    len > eeg_temp.eeg_header[:epoch_duration_seconds] && (len = eeg_temp.eeg_header[:epoch_duration_seconds])
    t === nothing && (t = (offset / fs) .+ collect(0:1/fs:len))
    t = t[1:(end - 1)]
    if offset < 0 || offset > eeg_temp.eeg_header[:epoch_duration_samples]
        throw(ArgumentError("Offset value out of range."))
    end
    if offset + (len * eeg_samplingrate(eeg_temp)) > eeg_temp.eeg_header[:epoch_duration_samples]
        throw(ArgumentError("Offset value or length value out of range."))
    end

    signal = Matrix(eeg_temp.eeg_signals[:, :, epoch])

    p = signal_plot_avg(t,
                        signal,
                        offset=offset,
                        normalize=normalize,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        title=title,
                        yamp=yamp)

    plot(p)

    if figure !== ""
        try
            savefig(p, figure)
        catch error
            throw(ArgumentError("File $figure cannot be saved."))
            return false
        end
    end

    return p
end

"""
    signal_plot_butterfly(t, signal; offset=0, labels=[""], normalize=true, xlabel"Time [s]", ylabel="Channels", title="Butterfly plot", yamp=nothing)

Butterfly plot of `signal` channels.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}`
- `signal::Matrix{Float64}`
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Float64}` - length in seconds
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `yamp::Union{Nothing, Int64, Float64}` - y-axis limits (-yamp:yamp)
"""
function signal_plot_butterfly(t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}, signal::Matrix{Float64}; offset::Int64=0, len::Union{Int64, Float64}=10.0, labels::Vector{String}=[""], normalize::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Butterfly plot", yamp::Union{Nothing, Int64, Float64}=nothing)
    
    if typeof(t) == UnitRange{Int64} || typeof(t) == StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
        t = float(collect(t))
    end
    
    channels_no = size(signal, 1)

    if normalize == true
        signal_normalized = signal_normalize_mean(reshape(signal, size(signal, 1), size(signal, 2), 1))
    else
        signal_normalized = signal
    end

    if yamp === nothing
        yamp = maximum(signal_normalized) * 1.5
        yamp = ceil(Int64, yamp)
    end

    if labels == [""]
        labels = Vector{String}(undef, channels_no)
        for idx in 1:channels_no
            labels[idx] = ""
        end
    end
    
    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(ceil(t[1]), floor(t[end])),
             ylims=(-yamp, yamp),
             title=title)
    for idx in 1:channels_no
        p = plot!(t,
                  signal_normalized[idx, (1 + offset):(offset + length(t))],
                  t=:line,
                  label=labels[idx])
    end

    return p
end

"""
    eeg_plot_butterfly(eeg; t, epoch=1, channels=nothing, offset=0, len=10.0, labels=[""], normalize=true, xlabel="Time [s]", ylabel="Channels", title="Butterfly plot", yamp=nothing, figure="")

Butterfly plot of `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `t::Union{Nothing, Vector{Float64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}` - the time vector
- `epoch::Int64` - epoch number to display
- `channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}` - channels to display
- `offset::Int64` - displayed segment offset in samples
- `len::Union{Int64, Float64}` - length in seconds
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `yamp::Union{Nothing, Int64, Float64}` - y-axis limits (-yamp:yamp)
- `figure::String` - name of the output figure file
"""
function eeg_plot_butterfly(eeg::EEG; t::Union{Vector{Float64}, UnitRange{Int64}, Nothing}=nothing, epoch::Int64=1, channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}=nothing, offset::Int64=0, len::Union{Int64, Float64}=10.0, labels::Vector{String}=[""], normalize::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Butterfly plot", yamp::Union{Nothing, Int64, Float64}=nothing, figure::String="")

    if epoch < 1 || epoch > eeg.eeg_header[:epochs_no]
        throw(ArgumentError("Epoch index out of range."))
    end

    if typeof(t) == UnitRange{Int64}
        t = collect(t)
    end

    # select channels, default is all channels
    channels === nothing && (channels = 1:eeg.eeg_header[:channels_no])
    eeg_temp = eeg_keep_channel(eeg, channels)

    eeg_temp = eeg_keep_channel(eeg, channels)
    signal = eeg_temp.eeg_signals[:, :, epoch]
    labels = eeg_temp.eeg_header[:labels]

    # default time is 10 seconds or epoch_duration_seconds
    fs = eeg_samplingrate(eeg_temp)
    len > eeg_temp.eeg_header[:epoch_duration_seconds] && (len = eeg_temp.eeg_header[:epoch_duration_seconds])
    t === nothing && (t = (offset / fs) .+ collect(0:1/fs:len))
    t = t[1:(end - 1)]
    if offset < 0 || offset > eeg_temp.eeg_header[:epoch_duration_samples]
        throw(ArgumentError("Offset value out of range."))
    end
    if offset + (len * eeg_samplingrate(eeg_temp)) > eeg_temp.eeg_header[:epoch_duration_samples]
        throw(ArgumentError("Offset value or length value out of range."))
    end

    p = signal_plot_butterfly(t,
                              signal,
                              offset=offset,
                              labels=labels,
                              normalize=normalize,
                              xlabel=xlabel,
                              ylabel=ylabel,
                              title=title,
                              yamp=yamp)

    plot(p)

    if figure !== ""
        try
            savefig(p, figure)
        catch error
            throw(ArgumentError("File $figure cannot be saved."))
            return false
        end
    end

    return p
end

"""
    signal_plot_psd(signal_powers::Vector{Float64}, signal_freqs::Vector{Float64}, normalize=false, xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="PSD")

Plots power spectrum density.

# Arguments

- `signal_powers::Vector{Float64}`
- `signal_freqs::Vector{Float64}`
- `frq_lim::Union{Nothing, Int64, Float64}` - x-axis limit
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
"""
function signal_plot_psd(signal_powers::Vector{Float64}, signal_freqs::Vector{Float64}; frq_lim::Union{Nothing, Int64, Float64}=nothing, xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="PSD")

    frq_lim === nothing && (frq_lim = signal_freqs[end])

    p = plot(signal_freqs,
             signal_powers,
             xlabel=xlabel,
             ylabel=ylabel,
             xlims=(0, frq_lim),
             legend=false,
             t=:line,
             c=:black,
             title=title)

    plot(p)

    return p
end

"""
    signal_plot_psd(signal::Vector{Float64}; fs::Int64, normalize=false, frq_lim=nothing, xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="PSD")

Plots power spectrum density.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64` - sampling frequency
- `normalize::Bool` - normalize the `signal` prior to calculations
- `frq_lim::Union{Nothing, Int64, Float64}` - x-axis limit
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
"""
function signal_plot_psd(signal::Vector{Float64}; fs::Int64, normalize::Bool=false, frq_lim::Union{Nothing, Int64, Float64}=nothing, xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="PSD")

    signal_freqs, signal_powers = signal_psd(signal, fs=fs, normalize=normalize)
    normalize == true && (ylabel::String="Power [dB]")

    frq_lim === nothing && (frq_lim = signal_freqs[end])

    p = plot(signal_freqs,
             signal_powers,
             xlabel=xlabel,
             ylabel=ylabel,
             xlims=(0, frq_lim),
             legend=false,
             t=:line,
             c=:black,
             title=title)

    plot(p)

    return p
end

"""
    signal_plot_psd(signal; fs, normalize=false, average=false, frq_lim=nothing, labels=[""], xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="PSD")

Plots power spectrum density.

# Arguments

- `signal::Matrix{Float64}`
- `fs::Int64` - sampling rate
- `normalize::Bool` - normalize the `signal` prior to calculations
- `average::Bool` - plots average power and 95%CI for all channels
- `frq_lim::Union{Nothing, Int64, Float64}` - x-axis limit
- `labels::Vector{String}` - channel labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
"""
function signal_plot_psd(signal::Matrix{Float64}; fs::Int64, normalize::Bool=false, average::Bool=false, frq_lim::Union{Nothing, Int64, Float64}=nothing, labels::Vector{String}=[""], xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="PSD")

    normalize == true && (ylabel="Power [dB]")

    channels_no = size(signal, 1)
    signal = reshape(signal, size(signal, 1), size(signal, 2), 1)
    signal_powers, signal_freqs = signal_psd(signal, fs=fs, normalize=normalize)
    signal_powers = signal_powers[:, :, 1]
    signal_freqs = signal_freqs[:, :, 1]

    frq_lim === nothing && (frq_lim = signal_freqs[1, end])
    if average == true
        signal_powers_m, signal_powers_s, signal_powers_u, signal_powers_l = signal_ci95(signal_powers)
        signal_freqs = signal_freqs[1, :]
        channels_no = 1
        labels == [""]
        title = "Average PSD with 95%CI"
    end

    if labels == [""]
        labels = Vector{String}(undef, channels_no)
        for idx in 1:channels_no
            labels[idx] = ""
        end
    end

    # plot channels
    p = plot(xlabel=xlabel,
             ylabel=ylabel,
             xlims=(0, frq_lim),
             title=title)
    if average == true
    p = plot!(signal_freqs,
              signal_powers_m,
              label=false,
              t=:line,
              c=:black)
    p = plot!(signal_freqs,
              signal_powers_u,
              label=false,
              t=:line,
              c=:grey,
              lw=0.5)
    p = plot!(signal_freqs,
              signal_powers_l,
              label=false,
              t=:line,
              c=:grey,
              lw=0.5)
    else
        for idx in 1:channels_no
            p = plot!(signal_freqs[idx, :],
                      signal_powers[idx, :],
                      label=labels[idx],
                      t=:line)
        end
    end

    return p
end

"""
    eeg_plot_psd(eeg; t, epoch=1, channels=nothing, labels=[""], normalize=false, average=false, frq_lim=nothing, xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="PSD", figure="")

Plots power spectrum density.

# Arguments

- `eeg::EEG` - EEG object
- `epoch::Int64` - epoch number to display
- `channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}` - channels to display
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `average::Bool` - plots average power and 95%CI for all channels
- `frq_lim::Union{Nothing, Int64, Float64}` - x-axis limit
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `title::String` - plot title
- `figure::String` - name of the output figure file
"""
function eeg_plot_psd(eeg::EEG; epoch::Int64=1, channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}=nothing, labels::Vector{String}=[""], normalize::Bool=false, average::Bool=false, frq_lim::Union{Nothing, Int64, Float64}=nothing, xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="PSD", figure::String="")

    if epoch < 1 || epoch > eeg.eeg_header[:epochs_no]
        throw(ArgumentError("Epoch index out of range."))
    end

    # select channels, default is all channels
    channels === nothing && (channels = 1:eeg.eeg_header[:channels_no])
    eeg_temp = eeg_keep_channel(eeg, channels)

    eeg_temp = eeg_keep_channel(eeg, channels)
    signal = eeg_temp.eeg_signals[:, :, epoch]
    fs = eeg_samplingrate(eeg_temp)
    labels = eeg_temp.eeg_header[:labels]

    p = signal_plot_psd(signal,
                        fs=fs,
                        labels=labels,
                        normalize=normalize,
                        average=average,
                        frq_lim=frq_lim,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        title=title)
    plot(p)

    if figure !== ""
        try
            savefig(p, figure)
        catch error
            throw(ArgumentError("File $figure cannot be saved."))
            return false
        end
    end

    return p
end

"""
    eeg_plot_electrodes(eeg; channels, labels=true, head=true, head_labels=false)

Plots electrodes.

# Arguments

- `eeg:EEG`
- `channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}` - channels to display
- `selected::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}` - which channels should be highlighted
- `labels::Bool` - plot electrode labels
- `head::Bool` - plot head
- `head_labels::Bool` - plot head labels
"""
function eeg_plot_electrodes(eeg::EEG; channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}=nothing, selected::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}=nothing, labels::Bool=true, head::Bool=true, head_labels::Bool=false)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available."))

    # select channels, default is all channels
    selected === nothing && (selected = 1:eeg.eeg_header[:channels_no])
    channels === nothing && (channels = 1:eeg.eeg_header[:channels_no])
    eeg_temp = eeg_keep_channel(eeg, channels)

    loc_x = eeg_temp.eeg_header[:xlocs]
    loc_y = eeg_temp.eeg_header[:ylocs]

    p = plot()
    head == true && eeg_draw_head(p, loc_x, loc_x, head_labels)
    if length(selected) >= eeg_temp.eeg_header[:channels_no]
        p = plot!(loc_x, loc_y, seriestype=:scatter, xlims=(-1, 1), ylims=(-1, 1), grid=true, label="")
    else
        p = plot!(loc_x, loc_y, seriestype=:scatter, color=:black, alpha=0.2, xlims=(-1, 1), ylims=(-1, 1), grid=true, label="")
        eeg_temp = eeg_keep_channel(eeg, selected)
        loc_x = eeg_temp.eeg_header[:xlocs]
        loc_y = eeg_temp.eeg_header[:ylocs]
        p = plot!(loc_x, loc_y, seriestype=:scatter, color=:red, xlims=(-1, 1), ylims=(-1, 1), grid=true, label="")
    end
    if labels == true
        for idx in 1:length(eeg_temp.eeg_header[:labels])
        plot!(p, annotation=(loc_x[idx] + 0.05, loc_y[idx] + 0.05, text(eeg_temp.eeg_header[:labels][idx], pointsize=8)))
        end
        p = plot!()
    end
    plot(p)

    return p
end