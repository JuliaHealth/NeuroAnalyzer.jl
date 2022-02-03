"""
    signal_plot(t, signal; offset=1, labels=[], normalize=true, xlabel="Time [s]", ylabel="Amplitude [μV]", yamp=nothing)

Plots `signal` against time vector `t`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}` - the time vector
- `signal::Vector{Float64}` - the signal vector
- `offset::Int64` - displayed segment offset in samples
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `average::Bool` - plot all channels averaged with 95%CI
- `butterfly::Bool` - plot all channels in butterfly mode
- `yamp::Union{Int64, Float64, Nothing}` - y-axis limits (-yamp:yamp)
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}, signal::Union{Vector{Float64}, Matrix{Float64}}; offset::Int64=1, labels::Vector{String}=[], normalize::Bool=true, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", average::Bool=false, butterfly::Bool=false, yamp::Union{Int64, Float64, Nothing}=nothing)

    if typeof(t) == UnitRange{Int64} || typeof(t) == StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
        t = float(collect(t))
    end

    if yamp === nothing
        yamp = maximum(signal)
        yamp = ceil(Int64, yamp)
    end

    if average == false
        p = plot(t, signal[offset:(offset + length(t))], xlabel=xlabel, ylabel=ylabel, legend=false, t=:line, c=:black, ylims=(-yamp, yamp))
    else
        m, s, u, l = signal_ci95(signal)
        p = plot(t, m[offset:(offset + length(t) - 1)], xlabel=xlabel, ylabel=ylabel, legend=false, t=:line, c=:black, ylims=(-yamp, yamp))
        p = plot!(t, u[offset:(offset + length(t) - 1)], c=:grey, lw=0.5)
        p = plot!(t, l[offset:(offset + length(t) - 1)], c=:grey, lw=0.5)
    end

    plot(p)

    return p
end

"""
    signal_plot(t, signal; offset=1, labels=[], normalize=true, xlabel="Time [s]", ylabel="Channels")

Plots `signal` matrix against time vector `t`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}` - the time vector
- `signal::Matrix{Float64}` - the signal matrix
- `offset::Int64` - displayed segment offset in samples
- `len::Float64` - length in seconds
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `average::Bool` - plot all channels averaged with 95%CI
- `butterfly::Bool` - plot all channels in butterfly mode
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}, signal::Matrix{Float64}; offset::Int64=1, labels::Vector{String}=[""], normalize::Bool=true, xlabel::String="Time [s]", ylabel::String="Channels", average::Bool=false, butterfly::Bool=false)
    
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
    if butterfly == false
        p = plot(xlabel=xlabel, ylabel=ylabel, ylim=(-0.5, channels_no-0.5))
        for idx in 1:channels_no
            p = plot!(t, signal_normalized[idx, offset:(offset + length(t))], legend=false, t=:line, c=:black)
        end
        p = plot!(p, yticks = (channels_no-1:-1:0, labels))
    else
        p = plot(t, signal[:, offset:(offset + length(t))]', xlabel=xlabel, ylabel="Amplitude [μV]")    
    end

    return p
end

"""
    eeg_plot(eeg; t=nothing, epoch=1, channels=nothing, offset=0, labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels", average=false, butterfly=false, figure=nothing)

Plots `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `t::Union{Nothing, Vector{Float64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}` - the time vector
- `epoch::Int64` - epoch number to display
- `channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}` - channels to display
- `offset::Int64` - displayed segment offset in samples
- `len::Float64` - length in seconds
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `average::Bool` - plot all channels averaged with 95%CI
- `butterfly::Bool` - plot all channels in butterfly mode
- `figure::String` - name of the output figure file
"""
function eeg_plot(eeg::EEG; t::Union{Vector{Float64}, UnitRange{Int64}, Nothing}=nothing, epoch::Int64=1, channels::Union{Nothing, Int64, Vector{Float64}, UnitRange{Int64}}=nothing, offset::Int64=1, len::Float64=10.0, labels::Vector{String}=[""], normalize::Bool=true, xlabel::String="Time [s]", ylabel::String="Channels", average::Bool=false, butterfly::Bool=false, figure::String="", overwrite::Bool=false)

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

    fs = eeg_temp.eeg_header[:sampling_rate][1]

    if average == false
        signal = eeg_temp.eeg_signals[:, :, epoch]
        labels = eeg_temp.eeg_header[:labels]
    else
        signal = Matrix(eeg_temp.eeg_signals[:, :, epoch])
        labels = [""]
    end

    # default time is 10 seconds or epoch_duration_seconds
    len > eeg_temp.eeg_header[:epoch_duration_seconds] && (len = eeg_temp.eeg_header[:epoch_duration_seconds])
    t === nothing && (t = collect(0:1/fs:len))
    t = t[1:(end - 2)]

    if offset < 1 || offset + (len * eeg_samplingrate(eeg)) > eeg.eeg_header[:epoch_duration_samples]
        throw(ArgumentError("Offset value out of range."))
    end

    p = signal_plot(t, signal, offset=offset, labels=labels, normalize=normalize, xlabel=xlabel, ylabel=ylabel, average=average, butterfly=butterfly)

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

Returns zero phase distortion filter response.  While saving, it does not check for overwrite.

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
    p1 = plot(w, H, title="Frequency response\nfilter: $(titlecase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order", xlims=(0, x_max), ylabel="Magnitude [dB]", xlabel="Frequency [Hz]", label="")
    if length(cutoff) == 1
        p1 = plot!((0, cutoff), seriestype=:vline, linestyle=:dash, label="")
    else
        p1 = plot!((0, cutoff[1]), seriestype=:vline, linestyle=:dash, label="")
        p1 = plot!((0, cutoff[2]), seriestype=:vline, linestyle=:dash, label="")
    end

    phi, w = phaseresp(ffilter)
    phi = rad2deg.(angle.(phi))
    # convert rad/sample to Hz
    w = w .* fs / 2 / pi
    x_max = w[end]
    ftype == :hp && (x_max = cutoff * 10)
    p2 = plot(w, phi, title="Phase response\nfilter: $(titlecase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order", ylims=(-180, 180), xlims=(0, x_max), ylabel="Phase [°]", xlabel="Frequency [Hz]", label="")
    if length(cutoff) == 1
        p2 = plot!((0, cutoff), seriestype=:vline, linestyle=:dash, label="")
    else
        p2 = plot!((0, cutoff[1]), seriestype=:vline, linestyle=:dash, label="")
        p2 = plot!((0, cutoff[2]), seriestype=:vline, linestyle=:dash, label="")
    end

    tau, w = grpdelay(ffilter)
    tau = abs.(tau)
    # convert rad/sample to Hz
    w = w .* fs / 2 / pi
    x_max = w[end]
    ftype == :hp && (x_max = cutoff * 10)
    p3 = plot(w, tau, title="Group delay\nfilter: $(titlecase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order", xlims=(0, x_max), ylabel="Group delay [samples]", xlabel="Frequency [Hz]", label="")
    if length(cutoff) == 1
        p3 = plot!((0, cutoff), seriestype=:vline, linestyle=:dash, label="")
    else
        p3 = plot!((0, cutoff[1]), seriestype=:vline, linestyle=:dash, label="")
        p3 = plot!((0, cutoff[2]), seriestype=:vline, linestyle=:dash, label="")
    end

    p = plot(p1, p2, p3, layout=(3, 1))

    if figure !== ""
        try
            savefig(p, figure)
        catch error
            throw(SystemError("File $figure cannot be saved."))
            return false
        end
    end
    
    return p
end