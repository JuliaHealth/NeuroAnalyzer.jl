"""
    eeg_plot_signal(eeg; <keyword arguments>)

Plot `eeg` channels. If signal is multi-channel, only channel amplitudes are plotted.

# Arguments

- `eeg::NeuroJ.EEG`: EEG object
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epochs to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: other arguments for plot() function; <keyword arguments>

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_signal_epoch(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

    # select channels, default is 20 channels
    channel != 0 && (channel = _select_channels(eeg, channel, 20))
    # get labels
    labels = eeg_labels(eeg)

    # get epochs
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    _check_epochs(eeg, epoch)
    length(epoch) > 1 && (epoch = sort!(epoch))
    offset = (epoch[1] - 1) * eeg_epoch_len(eeg)
    len = (epoch[end] * eeg_epoch_len(eeg)) - 1
    eeg_tmp, epoch_markers = _get_epoch_markers(eeg, offset, len)

    # get time vector
    if length(epoch) == 1 && len <= eeg_epoch_len(eeg)
        t = eeg.eeg_epochs_time[(1 + offset):length, epoch]
    else
        t = _get_t(eeg_tmp, offset, len)
    end

    _check_offset_len(eeg_tmp, offset, len)

    signal = eeg_tmp.eeg_signals[channel, (1 + offset):(offset + length(t)), epoch]

    if length(channel) == 1
        title = ""
        ylabel = "Amplitude [μV]"
        channel_name = labels[1]
        labels = [""]
        signal = vec(signal)
    end

    t_1, t_s1, t_2, t_s2 = _convert_t(t)

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