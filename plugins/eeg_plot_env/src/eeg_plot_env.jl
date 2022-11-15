export eeg_plot_env

"""
    eeg_plot_env(eeg; <keyword arguments>)

Plot envelope.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol`: envelope type: :amp (amplitude over time), :pow (power over frequencies), :spec (frequencies over time), :hamp (Hilbert spectrum amplitude)
- `average::Symbol`: averaging method: :no, :mean or :median
- `dims::Union{Int64, Nothing}=nothing`: average over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `epoch::Int64`: epoch number to display
- `channel::Int64`: channel to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String="default"`: plot title
- `y_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for PSD and spectrogram
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_env(eeg::NeuroAnalyzer.EEG; type::Symbol, average::Symbol=:no, dims::Union{Int64, Nothing}=nothing, d::Int64=32, epoch::Int64, channel::Int64, xlabel::String="", ylabel::String="", title::String="default", y_lim::Tuple{Real, Real}=(0, 0), frq_lim::Tuple{Real, Real}=(0, 0), mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    type in [:amp, :pow, :spec, :hamp] || throw(ArgumentError("type must be :amp, :pow, :spec or :hamp."))

    type === :amp && (d = 32)
    type === :hamp && (d = 8)
    type === :pow && (d = 8)
    type === :spec && (d = 8)

    type === :amp && (type = :amplitude)
    type === :hamp && (type = :hamplitude)
    type === :pow && (type = :power)
    type === :spec && (type = :spectrogram)

    average in [:no, :mean, :median] || throw(ArgumentError("average must be :no, :mean or :median."))
    (epoch < 1 || epoch > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    (channel < 1 || epoch > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    average === :no && (dims = nothing)
    (average !== :no && dims === nothing) && throw(ArgumentError("dims must be ≥ 1 and ≤ 3."))
    (average !== :no && (dims < 1 || dims > 3)) && throw(ArgumentError("dims must be ≥ 1 and ≤ 3."))

    fs = eeg_sr(eeg)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_lim = tuple_order(frq_lim)
    (frq_lim[1] < 0 || frq_lim[2] > fs / 2) && throw(ArgumentError("frq_lim must be ≥ 0 and ≤ $(fs / 2)."))

    t = eeg.eeg_epochs_time
    t[1] = floor(t[1], digits=2)
    t[end] = ceil(t[end], digits=2)

    if average === :no
        type === :amplitude && ((e, t) = eeg_tenv(eeg, d=d))
        type === :hamplitude && ((e, t) = eeg_henv(eeg, d=d))
        type === :power && ((e, t) = eeg_penv(eeg, d=d))
        type === :spectrogram && ((e, t) = eeg_senv(eeg, d=d))
    elseif average === :mean
        type === :amplitude && ((e, e_u, e_l, t) = eeg_tenv_mean(eeg, dims=dims, d=d))
        type === :hamplitude && ((e, e_u, e_l, t) = eeg_henv_mean(eeg, dims=dims, d=d))
        type === :power && ((e, e_u, e_l, t) = eeg_penv_mean(eeg, dims=dims, d=d))
        type === :spectrogram && ((e, e_u, e_l, t) = eeg_senv_mean(eeg, dims=dims, d=d))
    elseif average === :median
        type === :amplitude && ((e, e_u, e_l, t) = eeg_tenv_median(eeg, dims=dims, d=d))
        type === :hamplitude && ((e, e_u, e_l, t) = eeg_henv_median(eeg, dims=dims, d=d))
        type === :power && ((e, e_u, e_l, t) = eeg_penv_median(eeg, dims=dims, d=d))
        type === :spectrogram && ((e, e_u, e_l, t) = eeg_senv_median(eeg, dims=dims, d=d))
    end

    type === :amplitude && (xlabel == "" && (xlabel = "Time [s]"))
    type === :amplitude && (ylabel == "" && (ylabel = "Amplitude [μV]"))
    (type === :amplitude && y_lim == (0,0)) && (y_lim = (-200, 200))
    type === :amplitude && (x_lim = _xlims(t))
    type === :amplitude && (x_ticks = _ticks(t))

    type === :hamplitude && (xlabel == "" && (xlabel = "Time [s]"))
    type === :hamplitude && (ylabel == "" && (ylabel = "Amplitude"))
    if type === :hamplitude && y_lim == (0,0)
        if average === :no
            y_lim = (minimum(e[channel, :, epoch]) - 0.1 * minimum(e[channel, :, epoch]), maximum(e[channel, :, epoch]) + 0.1 * maximum(e[channel, :, epoch]))
        else
            y_lim = (minimum(e_l) - 0.1 * minimum(e_l), maximum(e_u) + 0.1 * maximum(e_u))
        end
    end
    type === :hamplitude && (x_lim = _xlims(t))
    type === :hamplitude && (x_ticks = _ticks(t))

    type === :power && (xlabel == "" && (xlabel = "Frequency [Hz]"))
    type === :power && (t = linspace(t[1], t[end], length(t)))
    type === :power && (ylabel == "" && (ylabel = "Power [dB/Hz]"))
    type === :power && (x_lim = (frq_lim[1], frq_lim[end]))
    type === :power && (x_ticks = round.(linspace(frq_lim[1], frq_lim[end], 10), digits=1))
    (type === :power && y_lim == (0,0)) && (y_lim = (-50, 50))

    type === :spectrogram && (xlabel == "" && (xlabel = "Time [s]"))
    type === :spectrogram && (ylabel == "" && (ylabel = "Frequency [Hz]"))
    type === :spectrogram && (x_lim = _xlims(t))
    (type === :spectrogram && y_lim == (0,0)) && (y_lim = frq_lim)
    type === :spectrogram && (x_ticks = _ticks(t))

    channel_name = eeg_labels(eeg)[channel]
    t_1, t_s1, t_2, t_s2 = _convert_t(t[1], t[end])

    if dims == 1
        e = e[:, epoch]
        average !== :no && (e_u = e_u[:, epoch]; e_l = e_l[:, epoch])
        if type === :hamplitude
            title == "default" && (title = "Envelope: Hilbert spectrum amplitude\n[$average of averaged channels, epoch: $epoch, time window: $t_s1:$t_s2]")
        else
            title == "default" && (title = "Envelope: $type\n[$average of averaged channels, epoch: $epoch, time window: $t_s1:$t_s2]")
        end
    elseif dims == 2
        e = e[:, channel]
        e_u = e_u[:, channel]
        e_l = e_l[:, channel]
        if type === :hamplitude
            title == "default" && (title = "Envelope: Hilbert spectrum amplitude\n[$average of averaged epochs, channel: $channel_name, time window: $t_s1:$t_s2]")
        else
            title == "default" && (title = "Envelope: $type\n[$average of averaged epochs, channel: $channel_name, time window: $t_s1:$t_s2]")
        end
    elseif dims == 3
        if type === :hamplitude
            title == "default" && (title = "Envelope: Hilbert spectrum amplitude\n[$average of averaged channels and epochs, time window: $t_s1:$t_s2]")
        else
            title == "default" && (title = "Envelope: $type\n[$average of averaged channels and epochs, time window: $t_s1:$t_s2]")
        end
    else
        e = e[channel, :, epoch]
        if type === :hamplitude
            title == "default" && (title = "Envelope: Hilbert spectrum amplitude\n[channel: $channel_name, epoch: $epoch, time window: $t_s1:$t_s2]")
        else
            title == "default" && (title = "Envelope: $type\n[channel: $channel_name, epoch: $epoch, time window: $t_s1:$t_s2]")
        end
    end

    p = Plots.plot(t,
                   e,
                   label="",
                   legend=false,
                   title=title,
                   xlabel=xlabel,
                   xlims=x_lim,
                   xticks=x_ticks,
                   ylabel=ylabel,
                   ylims=y_lim,
                   yguidefontrotation=0,
                   linewidth=0.5,
                   color=:black,
                   grid=true,
                   titlefontsize=8,
                   xlabelfontsize=6,
                   ylabelfontsize=6,
                   xtickfontsize=4,
                   ytickfontsize=4,
                   palette=pal;
                   kwargs...)
    if average !== :no
        p = Plots.plot!(t,
                        e_u,
                        fillrange=e_l,
                        fillalpha=0.35, 
                        label=false,
                        t=:line,
                        c=:grey,
                        linewidth=0.5)
        p = Plots.plot!(t,
                        e_l,
                        label=false,
                        t=:line,
                        c=:grey,
                        lw=0.5)
    end

    Plots.plot(p)

    return p
end