"""
    eeg_plot_signal(eeg; <keyword arguments>)

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
function eeg_plot_signal(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], xlabel::String="Time [s]", ylabel::String="", title::String="", head::Bool=true, hist::Symbol=:hist, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

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
    
    # select channels, default is 20 channels
    channel = _select_channels(eeg, channel, 20)

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
        _, _, _, s_phase = spectrum(signal)
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