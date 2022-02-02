"""
    signal_plot(t, signal; offset=1, labels=[], normalize=false, xlabel="Time [s]", ylabel="Amplitude [μV]", yamp=nothing)

Plots `signal` against time vector `t`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}` - the time vector
- `signal::Vector{Float64}` - the signal vector
- `offset::Int64` - displayed segment offset in samples
- `labels::Vector{String}` - channel labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis lable
- `yamp::Float64` - y-axis limits (-yamp:yamp)
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}, signal::Vector{Float64}; offset::Int64=1, labels::Vector{String}=[], xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", yamp::Union{Float64, Nothing}=nothing)

    if typeof(t) == UnitRange{Int64} || typeof(t) == StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
        t = float(collect(t))
    end

    if yamp === nothing
        yamp = maximum(signal)
        yamp = ceil(Int64, yamp)
    end

    p = plot(t, signal[ofset:(offset + length(t))], xlabel=xlabel, ylabel=ylabel, legend=false, t=:line, c=:black, ylims=(-yamp, yamp))

    plot(p)

    # TO DO: catching error while saving
    figure !== "" && (savefig(p, figure))

    return p
end

"""
    signal_plot(t, signal; offset=1, labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels")

Plots `signal` matrix against time vector `t`.

# Arguments

- `t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}` - the time vector
- `signal::Matrix{Float64}` - the signal matrix
- `offset::Int64` - displayed segment offset in samples
- `len::Float64` - length in seconds
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis lable
"""
function signal_plot(t::Union{Vector{Float64}, Vector{Int64}, UnitRange{Int64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}, signal::Matrix{Float64}; offset::Int64=1, labels::Vector{String}=[""], normalize::Bool=true, xlabel::String="Time [s]", ylabel::String="Channels")
    
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
    p = plot(xlabel=xlabel, ylabel=ylabel, ylim=(-0.5, channels_no-0.5))
    for idx in 1:channels_no
        p = plot!(t, signal_normalized[idx, offset:(offset + length(t))], legend=false, t=:line, c=:black)
    end
    p = plot!(p, yticks = (channels_no-1:-1:0, labels))
    return p
end

"""
    eeg_plot(eeg; t=nothing, epoch=1, offset=0, labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels", figure=nothing)

Plots `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `t::Union{Vector{Float64}, UnitRange{Int64}, Nothing}` - the time vector
- `epoch::Int64` - epoch number to display
- `offset::Int64` - displayed segment offset in samples
- `len::Float64` - length in seconds
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis lable
- `figure::String` - name of the output figure file
"""
function eeg_plot(eeg::EEG; t::Union{Vector{Float64}, UnitRange{Int64}, Nothing}=nothing, epoch::Int64=1, offset::Int64=1, len::Float64=10.0, labels::Vector{String}=[""], normalize::Bool=true, xlabel::String="Time [s]", ylabel::String="Channels", figure::String="")

    if epoch < 1 || epoch > eeg.eeg_header[:epochs_no]
        throw(ArgumentError("Epoch index out of range."))
    end

    if typeof(t) == UnitRange{Int64}
        t = collect(t)
    end

    signal = eeg.eeg_signals[:, :, epoch]
    labels = eeg.eeg_header[:labels]
    fs = eeg.eeg_header[:sampling_rate][1]

    # default time is 10 seconds or epoch_duration_seconds
    len > eeg.eeg_header[:epoch_duration_seconds] && (len = eeg.eeg_header[:epoch_duration_seconds])
    t === nothing && (t = collect(0:1/fs:len))
    t = t[1:(end - 2)]

    p = signal_plot(t, signal, offset=offset, labels=labels, normalize=normalize, xlabel=xlabel, ylabel=ylabel)

    plot(p)

    # TO DO: catching error while saving
    figure !== "" && (savefig(p, figure))
    
    return p
end

"""
    eeg_draw_head(p, loc_x, loc_y, add_labels=true)

Draws head over a topographical plot `p`.

# Arguments

- `p::Plot` - toppgraphical plot
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