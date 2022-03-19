"""
    eeg_plot_signal(eeg; <keyword arguments>)

Plot `eeg` channels. If signal is multi-channel, only channel amplitudes are plotted.

# Arguments

- `eeg::NeuroJ.EEG`: EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=0`: epochs to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
- `offset::Int64=0`: displayed segment offset in samples
- `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function; <keyword arguments>

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
epoch -> offset & len 

function eeg_plot_signal(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=0, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

    (epoch != 0 && len != 0) && throw(ArgumentError("Both epoch and len must not be specified."))

    if epoch != 0
        typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
        _check_epochs(eeg, epoch)
        if length(epoch) > 1
            sort!(epoch)
            len = eeg_epoch_len(eeg) * length(epoch)
            offset = eeg_epoch_len(eeg) * (epoch[1] - 1)
            epoch = epoch[1]
        end
        # default length is one epoch or 20 seconds
    else
        len = _len(eeg, len, 20)
    end

    # select channels, default is 20 channels
    channel != 0 && (channel = _select_channels(eeg, channel, 20))

    # set epoch markers if len > epoch_len
    eeg_tmp, epoch_markers = _get_epoch_markers(eeg, offset, len)

    labels = eeg_labels(eeg_tmp)

    # get time vector
    if length(epoch) == 1 && len <= eeg_epoch_len(eeg)
        t = eeh.eeg_epochs_time[, epoch]
    else
        t = _get_t(eeg_tmp, offset, len)
    end

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

    plot(p)

    return p
end