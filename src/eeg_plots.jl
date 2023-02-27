"""
    plot_save(p; file_name::String)

Saves plot as file (PDF/PNG/TIFF). File format is determined using `file_name` extension.

# Arguments

- `p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}`
- `file_name::String`
"""
function plot_save(p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}; file_name::String)

    ext = splitext(file_name)[2]
    _check_var(ext, [".png", ".pdf", ".jpg", ".tiff"], "File format")
    (isfile(file_name) && verbose == true) && _info("File $file_name will be overwritten.")
    if typeof(p) == Plots.Plot{Plots.GRBackend}
        savefig(p, file_name)
    else
        save(file_name, p)
    end

    nothing
end

"""
    plot_signal(t, signal; <keyword arguments>)

Plot amplitude of single- or multi-channel `signal`.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::Union{AbstractVector, AbstractArray}`: data to plot
- `labels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, signal::Union{AbstractVector, AbstractArray}; labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="μV", kwargs...)

    # convert single-channel signal to single-row matrix
    ndims(signal) == 1 && (signal = reshape(signal, 1, length(signal)))
    channel_n = size(signal, 1)

    # reverse so 1st channel is on top
    signal = @views reverse(signal[:, 1:length(t)], dims = 1)
    # also, reverse colors if palette is not mono
    if mono == true
        pal = :grays
        channel_color = Vector{Symbol}()
        for idx in 1:channel_n
            push!(channel_color, :black)
        end
    else
        pal = :darktest
        channel_color = channel_n:-1:1
    end

    # get range of the original signal for the scale
    range = _get_range(signal)

    # normalize and shift so all channels are visible
    # each channel is between -1.0 and +1.0
    for idx in 1:channel_n
        # scale by 0.5 so maxima do not overlap
        signal[idx, :] = @views s_normalize(signal[idx, :], method=:minmax) .* 0.5 .+ (idx - 1)
    end

    # prepare plot
    channel_n == 1 && (plot_size = (1200, 500))
    channel_n > 1 && (plot_size = (1200, 800))
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=(-1, channel_n),
                   title=title,
                   palette=pal,
                   size=plot_size,
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot zero line
    p = Plots.hline!(collect((channel_n - 1):-1:0),
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot channels
    for idx in 1:channel_n
        p = @views Plots.plot!(t,
                               signal[idx, :],
                               linewidth=1,
                               label="",
                               color=channel_color[idx])
    end

    # plot labels
    p = Plots.plot!(yticks=((channel_n - 1):-1:0, labels))

    # draw scale
    if scale == true
        p = Plots.plot!([t[1], t[1]], [(channel_n - 1.5), (channel_n - 0.5)], color=:red, linewidth=5, label="")
        p = Plots.plot!(annotation=(t[1], (channel_n - 1), Plots.text("$range$units", pointsize=6, halign=:center, valign=:bottom, rotation=90)), label=false)
    end

    return p
end

"""
    plot_signal(t, signal, bad; <keyword arguments>)

Plot amplitude of single- or multi-channel `signal`.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::Union{AbstractVector, AbstractArray}`: data to plot
- `norm::Bool=false`: normalize signal for butterfly and averaged plots
- `bad::Vector{Bool}}`: list of bad channels
- `labels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, signal::Union{AbstractVector, AbstractArray}, bad::Vector{Bool}; labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", scale::Bool=true, units::String="μV", kwargs...)

    length(bad) == size(signal, 1) || throw(ArgumentError("Length of bad channels vector and number of channels must be equal."))

    # convert single-channel signal to single-row matrix
    ndims(signal) == 1 && (signal = reshape(signal, 1, length(signal)))
    channel_n = size(signal, 1)

    # reverse so 1st channel is on top
    signal = @views reverse(signal[:, 1:length(t)], dims = 1)
    bad = reverse(bad)

    pal = mono == true ? :grays : :darktest

    # get range of the original signal for the scale
    range = _get_range(signal)

    # normalize and shift so all channels are visible
    # each channel is between -1.0 and +1.0
    for idx in 1:channel_n
        # scale by 0.5 so maxima do not overlap
        signal[idx, :] = @views s_normalize(signal[idx, :], method=:minmax) .* 0.5 .+ (idx - 1)
    end

    # prepare plot
    channel_n == 1 && (plot_size = (1200, 500))
    channel_n > 1 && (plot_size = (1200, 800))
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=(-1, channel_n),
                   title=title,
                   palette=pal,
                   size=plot_size,
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)
    
    # plot zero line
    p = Plots.hline!(collect((channel_n - 1):-1:0),
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot channels
    for idx in 1:channel_n
        if bad[idx] == true
            p = @views Plots.plot!(t,
                                   signal[idx, :],
                                   linewidth=1,
                                   label="",
                                   color=:red)
        else
            p = @views Plots.plot!(t,
                                   signal[idx, :],
                                   linewidth=1,
                                   label="",
                                   color=:black)
        end
    end

    # plot labels
    p = Plots.plot!(yticks=((channel_n - 1):-1:0, labels))

    # draw scale
    if scale == true
        p = Plots.plot!([t[1], t[1]], [(channel_n - 1.5), (channel_n - 0.5)], color=:red, linewidth=5, label="")
        p = Plots.plot!(annotation=(t[1], (channel_n - 1), Plots.text("$range$units", pointsize=6, halign=:center, valign=:bottom, rotation=90)), label=false)
    end

    return p
end

"""
    plot_signal_avg(t, signal; <keyword arguments>)

Plot amplitude mean and ±95% CI of averaged `signal` channels.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::AbstractArray`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `norm::Bool=false`: normalize to -1 .. +1
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal_avg(t::Union{AbstractVector, AbstractRange}, signal::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="μV", norm::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get range of the original signal for the scale
    range = _get_range(signal)

    # get mean and 95%CI
    s_m, _, s_u, s_l = s_msci95(signal)

    # get limits
    if norm != true
        ylim = (floor(minimum(s_l), digits=0), ceil(maximum(s_u), digits=0))
        ylim = _tuple_max(ylim)
        yticks = [ylim[1], 0, ylim[2]]
    else
        s_m = s_normalize(s_m, method=:minmax)
        s_u = s_normalize(s_u, method=:minmax)
        s_l = s_normalize(s_l, method=:minmax)
        ylim = (-1.0, 1.0)
        yticks = [0]
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)
    # plot upper 95% CI
    p = Plots.plot!(t,
                    s_u,
                    fillrange=s_l,
                    fillalpha=0.35, 
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot lower 95% CI
    p = Plots.plot!(t,
                    s_l,
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot mean
    p = Plots.plot!(t,
                    s_m,
                    label=false,
                    t=:line,
                    c=:black,
                    lw=0.5)

    # draw scale
    if norm == true && scale == true
        p = Plots.plot!([t[1], t[1]], [-1, 1], color=:red, linewidth=5, label=false)
        p = Plots.plot!(annotation=(t[1], 0, Plots.text("$range$units", pointsize=6, halign=:center, valign=:bottom, rotation=90)), label=false)
    end

    return p
end

"""
    plot_signal_butterfly(t, signal; <keyword arguments>)

Butterfly plot of `signal` channels.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::AbstractArray`: data to plot
- `labels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `mono::Bool=false`: use color or grey palette
- `norm::Bool=false`: normalize to -1 .. +1
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal_butterfly(t::Union{AbstractVector, AbstractRange}, signal::AbstractArray; labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="μV", norm::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get range of the original signal for the scale
    range = _get_range(signal)

    channel_n = size(signal, 1)

    # get limits
    if norm != true
        ylim = (floor(minimum(signal), digits=0), ceil(maximum(signal), digits=0))
        ylim = _tuple_max(ylim)
        yticks = [ylim[1], 0, ylim[2]]
    else
        signal = s_normalize(signal, method=:minmax)
        ylim = (-1.0, 1.0)
        yticks = [0]
    end

    # channel labels
    labels == [""] && (labels = repeat([""], channel_n))
    
    # plot channels
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)
    for idx in 1:channel_n
        p = Plots.plot!(t,
                        signal[idx, :],
                        t=:line,
                        linecolor=idx,
                        linewidth=0.5,
                        label=labels[idx],
                        legend=true)
    end
    
    # draw scale
    if norm == true && scale == true
        p = Plots.plot!([t[1], t[1]], [-1, 1], color=:red, linewidth=5, label=false)
        p = Plots.plot!(annotation=(t[1], 0, Plots.text("$range$units", pointsize=6, halign=:center, valign=:bottom, rotation=90)), label=false)
    end

    return p
end

"""
    eeg_plot(eeg; <keyword arguments>)

Plot signal.

# Arguments

- `eeg::NeuroAnalyzer.EEG`: EEG object
- `epoch::Union{Int64, AbstractRange}=0`: epoch to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: channel(s) to plot, default is all channels
- `segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title, default is Amplitude [channels: 1:2, epochs: 1:2, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or grey palette
- `emarkers::Bool`: draw epoch markers if available
- `markers::Bool`: draw markers if available
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`)
- `norm::Bool=false`: normalize signal for butterfly and averaged plots
- `bad::Union{Bool, Matrix{Bool}}=false`: list of bad channels; if not empty - plot bad channels using this list
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, AbstractRange}=0, channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg)), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, markers::Bool=true, scale::Bool=true, units::String="μV", type::Symbol=:normal, norm::Bool=false, bad::Union{Bool, Matrix{Bool}}=false, kwargs...)

    eeg_signal_len(eeg) < 10 * eeg_sr(eeg) && segment == (1, 10*eeg_sr(eeg)) && (segment=(1, eeg_signal_len(eeg)))

    _check_var(type, [:normal, :butterfly, :mean], "type")
    _check_segment(eeg, segment[1], segment[2])

    if epoch != 0
        _check_epochs(eeg, epoch)
        segment = (((epoch[1] - 1) * eeg_epoch_len(eeg) + 1), segment[2])
        if typeof(epoch) == Int64
            segment = (segment[1], (segment[1] + eeg_epoch_len(eeg) - 1))
        else
            segment = (segment[1], (epoch[end] * eeg_epoch_len(eeg)))
        end
    end

    # do not show epoch markers if there are no epochs
    eeg_epoch_n(eeg) == 1 && (emarkers = false)
    if emarkers == true
        epoch_markers = _get_epoch_markers(eeg)
    end

    # check channels
    _check_channels(eeg, channel)
    labels = eeg_labels(eeg)[channel]
    length(channel) == 1 && (labels = [labels])

    # get time vector
    if segment[2] <= eeg_epoch_len(eeg)
        signal = eeg.eeg_signals[channel, segment[1]:segment[2], 1]
    else
        signal = eeg_epoch(eeg, epoch_n=1).eeg_signals[channel, segment[1]:segment[2], 1]
    end
    t = _get_t(segment[1], segment[2], eeg_sr(eeg))

    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    epoch = _s2epoch(eeg, segment[1], segment[2])

    if type === :normal
        if bad == false
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Channel$(_pl(length(channel))) $(_channel2channel_name(channel)) amplitude\n[epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
            p = plot_signal(t,
                            signal,
                            labels=labels,
                            xlabel=xlabel,
                            ylabel=ylabel,
                            title=title,
                            scale=scale,
                            units=units,
                            mono=mono;
                            kwargs...)
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Bad channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
            length(channel) > size(bad, 1) && throw(ArgumentError("Number of channels cannot be larger than number of bad channels rows."))
            epoch > size(bad, 2) && throw(ArgumentError("Epoch number cannot be larger than number of bad channels columns."))
            p = plot_signal(t,
                            signal,
                            bad[channel, epoch],
                            labels=labels,
                            xlabel=xlabel,
                            ylabel=ylabel,
                            title=title,
                            scale=scale,
                            units=units;
                            kwargs...)
        end
    elseif type === :butterfly
        size(signal, 1) == 1 && throw(ArgumentError("For type=:butterfly plot the signal must contain ≥ 2 channels."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [μV]", "Channels $(_channel2channel_name(channel)) amplitude\n[epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
        p = plot_signal_butterfly(t,
                                  signal,
                                  labels=labels,
                                  xlabel=xlabel,
                                  ylabel=ylabel,
                                  title=title,
                                  scale=scale,
                                  units=units,
                                  norm=norm,
                                  mono=mono;
                                  kwargs...)
    elseif type === :mean
        size(signal, 1) == 1 && throw(ArgumentError("For type=:mean plot the signal must contain ≥ 2 channels."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [μV]", "Averaged channels $(_channel2channel_name(channel)) amplitude [mean ± 95%CI]\n [epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
        p = plot_signal_avg(t,
                            signal,
                            xlabel=xlabel,
                            ylabel=ylabel,
                            title=title,
                            scale=scale,
                            units=units,
                            norm=norm,
                            mono=mono;
                            kwargs...)
    end

    # add epochs markers
    # TODO: draw epoch numbers
    if emarkers == true
        p = Plots.vline!(epoch_markers,
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:blue,
                         label="")
    end

    # plot markers if available
    # TODO: draw markers length
    if markers == true && eeg.eeg_header[:markers] == true
        markers_pos = eeg.eeg_markers[!, :start] ./ eeg_sr(eeg)
        markers_desc = eeg.eeg_markers[!, :description]
        p = Plots.vline!(markers_pos,
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:black,
                         label=false)
        for idx in eachindex(markers_desc)
            p = Plots.plot!(annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
        end
    end

    Plots.plot(p)

    return p
end

"""
    eeg_plot(eeg, c; <keyword arguments>)

Plot embedded or external component.

# Arguments

- `eeg::NeuroAnalyzer.EEG`: EEG object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `epoch::Union{Int64, AbstractRange}=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title, default is Amplitude [channels: 1:2, epochs: 1:2, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or grey palette
- `emarkers::Bool`: draw epoch markers if available
- `markers::Bool`: draw markers if available
- `scale::Bool=true`: draw scale
- `units::String=""`: units of the scale
- `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`)
- `norm::Bool=false`: normalize signal for butterfly and averaged plots
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot(eeg::NeuroAnalyzer.EEG, c::Union{Symbol, AbstractArray}; epoch::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg)), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, markers::Bool=true, scale::Bool=true, units::String="", type::Symbol=:normal, norm::Bool=false, kwargs...)

    eeg_signal_len(eeg) < 10 * eeg_sr(eeg) && segment == (1, 10*eeg_sr(eeg)) && (segment=(1, eeg_signal_len(eeg)))

    _check_var(type, [:normal, :butterfly, :mean], "type")
    _check_segment(eeg, segment[1], segment[2])

    if epoch != 0
        _check_epochs(eeg, epoch)
        segment = (((epoch[1] - 1) * eeg_epoch_len(eeg) + 1), segment[2])
        if typeof(epoch) == Int64
            segment = (segment[1], (segment[1] + eeg_epoch_len(eeg) - 1))
        else
            segment = (segment[1], (epoch[end] * eeg_epoch_len(eeg)))
        end
    end

    # do not show epoch markers if there are no epochs
    eeg_epoch_n(eeg) == 1 && (emarkers = false)
    if emarkers == true
        epoch_markers = _get_epoch_markers(eeg)
    end

    # select component channels, default is all channels
    typeof(c) == Symbol && (c = _get_component(eeg, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    labels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (labels = [labels])

    # get time vector
    if segment[2] <= eeg_epoch_len(eeg)
        signal = c[c_idx, segment[1]:segment[2], 1]
    else
        signal = _make_epochs(c, epoch_n=1)[c_idx, segment[1]:segment[2], 1]
    end
    t = _get_t(segment[1], segment[2], eeg_sr(eeg))

    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    epoch = _s2epoch(eeg, segment[1], segment[2])

    if type === :normal
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx)) amplitude\n[epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
        p = plot_signal(t,
                        signal,
                        labels=labels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        title=title,
                        scale=scale,
                        units=units,
                        mono=mono;
                        kwargs...)
    elseif type === :butterfly
        size(signal, 1) == 1 && throw(ArgumentError("For type=:butterfly plot the signal must contain ≥ 2 channels."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [μV]", "Components $(_channel2channel_name(c_idx)) amplitude\n[epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
        p = plot_signal_butterfly(t,
                                  signal,
                                  labels=labels,
                                  xlabel=xlabel,
                                  ylabel=ylabel,
                                  title=title,
                                  scale=scale,
                                  units=units,
                                  norm=norm,
                                  mono=mono;
                                  kwargs...)
    elseif type === :mean
        size(signal, 1) == 1 && throw(ArgumentError("For type=:mean plot the signal must contain ≥ 2 channels."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [μV]", "Averaged components $(_channel2channel_name(c_idx)) amplitude [mean ± 95%CI]\n[epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
        p = plot_signal_avg(t,
                            signal,
                            labels=labels,
                            xlabel=xlabel,
                            ylabel=ylabel,
                            title=title,
                            scale=scale,
                            units=units,
                            norm=norm,
                            mono=mono;
                            kwargs...)
    end

    # add epochs markers
    # TODO: draw epoch numbers
    if emarkers == true
        p = Plots.vline!(epoch_markers,
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:blue,
                         label="")
    end

    # plot markers if available
    # TODO: draw markers length
    if markers == true && eeg.eeg_header[:markers] == true
        markers_pos = eeg.eeg_markers[!, :start] ./ eeg_sr(eeg)
        markers_desc = eeg.eeg_markers[!, :description]
        p = Plots.vline!(markers_pos,
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:black,
                         label=false)
        for idx in eachindex(markers_desc)
            p = Plots.plot!(annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
        end
    end

    Plots.plot(p)

    return p
end

"""
    plot_psd(s_frq, s_pow; <keyword arguments>)

Plot PSD (power spectrum density).

# Arguments

- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Vector{Float64}`: powers
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd(s_frq::Vector{Float64}, s_pow::Vector{Float64}; norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    length(s_pow) == length(s_frq) || throw(ArgumentError("Length of powers vector must equal length of frequencies vector."))
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    frq_lim == (0, 0) && (frq_lim = (s_frq[1], s_frq[end]))
    frq_lim = tuple_order(frq_lim)

    pal = mono == true ? :grays : :darktest

    if ax === :linlin
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = norm == false ? :log10 : :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6)

    # plot powers
    p = Plots.plot!(s_frq,
                    s_pow,
                    xticks=xticks,
                    xscale=xscale,
                    yscale=yscale;
                    kwargs...)

    return p
end


"""
    plot_psd(s_frq, s_pow; <keyword arguments>)

Plot multi-channel PSD (power spectrum density).

# Arguments

- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Matrix{Float64}`: powers
- `labels::Vector{String}=[""]`: signal channel labels vector
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd(s_frq::Vector{Float64}, s_pow::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    channel_n = size(s_pow, 1)
    size(s_pow, 2) == length(s_frq) || throw(ArgumentError("Length of powers vector must equal length of frequencies vector."))
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    # reverse so 1st channel is on top
    s_pow = @views reverse(s_pow[:, 1:length(s_frq)], dims = 1)
    # also, reverse colors if palette is not mono
    if mono == true
        pal = :grays
        channel_color = Vector{Symbol}()
        for idx in 1:channel_n
            push!(channel_color, :black)
        end
    else
        pal = :darktest
        channel_color = channel_n:-1:1
    end

    # channel labels
    labels == [""] && (labels = repeat([""], size(s_pow, 1)))

    # get range of the original s_pow for the scale
    range = _get_range(s_pow)

    # normalize and shift so all channels are visible
    # each channel is between -1.0 and +1.0
    for idx in 1:channel_n
        # scale by 0.5 so maxima do not overlap
        s_pow[idx, :] = @views s_normalize(s_pow[idx, :], method=:minmax) .* 0.5 .+ (idx - 1)
    end

    frq_lim == (0, 0) && (frq_lim = (s_frq[1], s_frq[end]))
    frq_lim = tuple_order(frq_lim)

    if ax === :linlin
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        _info("For multi-channel PSD plots, y-axis log-scale is ignored.")
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = :identity
    elseif ax === :loglog
        _info("For multi-channel PSD plots, y-axis log-scale is ignored.")
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 800),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6)

    # plot zero line
    p = Plots.hline!(collect((channel_n - 1):-1:0),
                     color=:grey,
                     lw=0.5,
                     label="")

    # plot channels
    for idx in 1:channel_n
        p = @views Plots.plot!(s_frq,
                               s_pow[idx, :],
                               linewidth=1,
                               label="",
                               xticks=xticks,
                               xscale=xscale,
                               color=channel_color[idx])
    end

    # plot labels
    p = Plots.plot!(yticks=((channel_n - 1):-1:0, labels))

    return p
end

"""
    plot_psd_avg(s_frq, s_pow; <keyword arguments>)

Plot PSD mean and ±95% CI of averaged channels.

# Arguments

- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Array{Float64, 3}`: powers
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd_avg(s_frq::Vector{Float64}, s_pow::Array{Float64, 2}; norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    size(s_pow, 2) == length(s_frq) || throw(ArgumentError("Length of powers vector must equal length of frequencies vector."))
    _check_var(ax,[:linlin, :loglin, :linlog, :loglog], "ax")

    frq_lim == (0, 0) && (frq_lim = (s_frq[1], s_frq[end]))
    frq_lim = tuple_order(frq_lim)

    pal = mono == true ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = s_msci95(s_pow)

    if ax === :linlin
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = :identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = norm == false ? :log10 : :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   xticks=xticks,
                   xscale=xscale,
                   yscale=yscale,
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot upper 95% CI
    p = Plots.plot!(s_frq,
                    s_u,
                    fillrange=s_l,
                    fillalpha=0.35, 
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot lower 95% CI
    p = Plots.plot!(s_frq,
                    s_l,
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot mean
    p = Plots.plot!(s_frq,
                    s_m,
                    label=false,
                    t=:line,
                    c=:black,
                    lw=0.5)

    return p
end

"""
    plot_psd_butterfly(s_frq, s_pow; <keyword arguments>)

Butterfly PSD plot.

# Arguments

- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Array{Float64, 3}`: powers
- `labels::Vector{String}=[""]`: signal channel labels vector
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_psd_butterfly(s_frq::Vector{Float64}, s_pow::Array{Float64, 2}; labels::Vector{String}=[""], norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    size(s_pow, 2) == length(s_frq) || throw(ArgumentError("Length of powers vector must equal length of frequencies vector."))
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    frq_lim == (0, 0) && (frq_lim = (s_frq[1], s_frq[end]))
    frq_lim = tuple_order(frq_lim)

    pal = mono == true ? :grays : :darktest
    
    # channel labels
    labels == [""] && (labels = repeat([""], size(s_pow, 1)))

    if ax === :linlin
        xticks=_ticks(frq_lim)
        xscale=:identity
        yscale=:identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = norm == false ? :log10 : :identity
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   legend=false,
                   xlims=frq_lim,
                   xticks=xticks,
                   xscale=xscale,
                   yscale=yscale;
                   title=title,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6)

    # plot powers
    for idx in 1:size(s_pow, 1)
        p = Plots.plot!(s_frq,
                        s_pow[idx, :],
                        t=:line,
                        linecolor=idx,
                        linewidth=0.5,
                        label=labels[idx],
                        legend=true;
                        kwargs...)
    end

    return p

end

"""
    plot_psd_w3d(s_frq, s_pow; <keyword arguments>)

Plot 3-d waterfall PSD plot.

# Arguments

- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Array{Float64, 3}`: powers
- `labels::Vector{String}=[""]`: signal channel labels vector
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `zlabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `variant::Symbol`: waterfall (`:w`) or surface (`:s`)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}`
"""
function plot_psd_3d(s_frq::Vector{Float64}, s_pow::Array{Float64, 2}; labels::Vector{String}=[""], norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", zlabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, variant::Symbol, kwargs...)

    _check_var(variant, [:w, :s], "variant")
    size(s_pow, 2) == length(s_frq) || throw(ArgumentError("Length of powers vector must equal length of frequencies vector."))
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    frq_lim == (0, 0) && (frq_lim = (s_frq[1], s_frq[end]))
    frq_lim = tuple_order(frq_lim)

    channel_n = size(s_pow, 1)

    pal = mono == true ? :grays : :darktest
    
    # channel labels
    labels == [""] && (labels = repeat([""], channel_n))

    if ax === :linlin
        xticks=_ticks(frq_lim)
        xscale=:identity
        zscale=:identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        zscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        zscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        zscale = norm == false ? :log10 : :identity
    end

    # prepare plot
    if variant === :w
        p = Plots.plot(s_frq,
                       ones(length(s_frq)),
                       s_pow[1, :],
                       xlabel=xlabel,
                       ylabel=ylabel,
                       zlabel=zlabel,
                       legend=false,
                       xlims=frq_lim,
                       xticks=xticks,
                       xscale=xscale,
                       zscale=zscale;
                       title=title,
                       palette=pal,
                       st=:line,
                       lc=:black,
                       size=(1200, 800),
                       margins=20Plots.px,
                       titlefontsize=8,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6)

        # plot powers
        for idx in 2:channel_n
            p = Plots.plot!(s_frq,
                            ones(length(s_frq)) .* idx,
                            s_pow[idx, :],
                            st=:line,
                            linecolor=idx,
                            linewidth=0.5,
                            kwargs...)
        end
    else
        f1 = vsearch(frq_lim[1], s_frq)
        f2 = vsearch(frq_lim[2], s_frq)
        p = Plots.plot(s_frq[f1:f2],
                       1:length(labels),
                       s_pow[:, f1:f2],
                       xlabel=xlabel,
                       ylabel=ylabel,
                       zlabel=zlabel,
                       legend=false,
                       xlims=frq_lim,
                       xticks=xticks,
                       xscale=xscale,
                       zscale=zscale;
                       title=title,
                       palette=pal,
                       st=:surface,
                       lc=:black,
                       size=(1200, 800),
                       margins=20Plots.px,
                       titlefontsize=8,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6)
    end

    p = Plots.plot!(yticks=(1:channel_n, labels))

    return p
end

"""
    plot_psd_topo(locs, s_frq, s_pow; <keyword arguments>)

Plot topographical map PSDs. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Array{Float64, 3}`: powers
- `channel::Union{Vector{Int64}, AbstractRange}`: which channels to plot
- `labels::Vector{String}=[""]`: signal channel labels vector
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `kwargs`: optional arguments for plot() function

# Returns

- `fig::GLMakie.Figure`
"""
function plot_psd_topo(locs::DataFrame, s_frq::Vector{Float64}, s_pow::Array{Float64, 2}; channel=Union{Vector{Int64}, AbstractRange}, labels::Vector{String}=[""], norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

    size(s_pow, 2) == length(s_frq) || throw(ArgumentError("Length of powers vector must equal length of frequencies vector."))
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")

    length(channel) > nrow(locs) && throw(ArgumentError("Some channels do not have locations."))

    frq_lim == (0, 0) && (frq_lim = (s_frq[1], s_frq[end]))
    frq_lim = tuple_order(frq_lim)

    pal = mono == true ? :grays : :darktest
    
    # channel labels
    labels == [""] && (labels = repeat([""], size(s_pow, 1)))

    if ax === :linlin
        xticks=_ticks(frq_lim)
        xscale=:identity
        yscale=:identity
    elseif ax === :loglin
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = :identity
    elseif ax === :linlog
        xticks = _ticks(frq_lim)
        xscale = :identity
        yscale = norm == false ? :log10 : :identity
    elseif ax === :loglog
        if frq_lim[1] == 0
            frq_lim = (0.1, frq_lim[2])
            _info("Lower frequency bound truncated to 0.1 Hz")
        end
        s_frq[1] == 0 && (s_frq[1] = 0.1)
        xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
        xscale = :log10
        yscale = norm == false ? :log10 : :identity
    end

    # plot parameters
    plot_size = 1200
    marker_size = (150, 75)
    
    # get locations
    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = loc_x[channel]
    loc_y = loc_y[channel]
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)
    # get marker centers
    loc_x .*= ((plot_size / 2) - marker_size[1] / 2)
    loc_y .*= ((plot_size / 2) - marker_size[2] / 2)

    fig = Figure(; resolution=(plot_size, plot_size))
    fig_axis = Axis(fig[1, 1])
    fig_axis.aspect = AxisAspect(1)
    fig_axis.title = title
    GLMakie.xlims!(fig_axis, [-plot_size / 1.75, plot_size / 1.75])
    GLMakie.ylims!(fig_axis, [-plot_size / 1.75, plot_size / 1.75])
    hidedecorations!(fig_axis, grid=true, ticks=true)

    for idx in 1:size(s_pow, 1)
        p = Plots.plot(s_frq,
                       s_pow[idx, :],
                       t=:line,
                       c=:black,
                       linewidth=0.5,
                       xlabel=xlabel,
                       ylabel=ylabel,
                       legend=false,
                       xlims=frq_lim,
                       xticks=false,
                       yticks=false,
                       xscale=xscale,
                       yscale=yscale,
                       title=labels[idx],
                       palette=pal,
                       size=marker_size,
                       #left_margin=20Plots.px,
                       titlefontsize=8,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6;
                       kwargs...)
        marker_img = tempname() * ".png"
        savefig(p, marker_img)
        marker = load(marker_img)
        GLMakie.scatter!(fig_axis, (loc_x[idx], loc_y[idx]), marker=marker, markersize=marker_size)
        rm(marker_img)
    end

    return fig
end

"""
    eeg_plot_psd(eeg; <keyword arguments>)

Plot power spectrum density.

# Arguments

- `eeg::NeuroAnalyzer.EEG`: EEG object
- `epoch::Int64`: epoch to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:welch`: method of calculating PSD:
    - `:welch`: Welch's periodogram
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to EEG band: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher` 
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is Power [dB] or Power [μV^2/Hz]
- `zlabel::String="default"`: z-axis label for 3-d plots, default is Power [dB] or Power [μV^2/Hz]
- `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or grey palette
- `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}`
"""
function eeg_plot_psd(eeg::NeuroAnalyzer.EEG; epoch::Int64, channel::Union{Int64, Vector{Int64}, AbstractRange}, norm::Bool=true, method::Symbol=:welch, nt::Int64=8, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, kwargs...)

    _check_var(type, [:normal, :butterfly, :mean, :w3d, :s3d, :topo], "type")
    _check_var(method, [:welch, :mt, :mw], "method")
    _check_var(ref, [:abs, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "ref")
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")
    ref !== :abs && method === :mw && throw(ArgumentError("For relative PSD, method must be :welch or :mt."))

    _check_epochs(eeg, epoch)
    _check_channels(eeg, channel)

    labels = eeg_labels(eeg)[channel]
    length(channel) == 1 && (labels = [labels])

    ref !== :abs && (f = eeg_band(eeg, band=ref))

    # get frequency range
    fs = eeg_sr(eeg)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_lim = tuple_order(frq_lim)
    (frq_lim[1] < 0 || frq_lim[2] > fs / 2) && throw(ArgumentError("frq_lim must be ≥ 0 and ≤ $(fs / 2)."))

    # calculate PSD
    signal = eeg.eeg_signals[channel, :, epoch]

    # get time vector
    _, t_s1, _, t_s2 = _convert_t(eeg.eeg_epoch_time[1], eeg.eeg_epoch_time[end])

    if ref === :abs
        if method === :welch
            s_pow, s_frq = s_psd(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mt
            s_pow, s_frq = s_psd(signal, fs=fs, norm=norm, mt=true, nt=nt)
            title == "default" && (title = "Absolute PSD (multi-tapered) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mw
            s_pow, s_frq = s_mwpsd(signal, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc)
            title == "default" && (title = "Absolute PSD (Morlet wavelet convolution) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        end
    elseif ref === :total
        if method === :welch
            s_pow, s_frq = s_rel_psd(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "PSD (Welch's periodogram) relative to total power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mt
            s_pow, s_frq = s_rel_psd(signal, fs=fs, norm=norm, mt=true)
            title == "default" && (title = "PSD (multi-tapered) relative to total power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        end
    else
        if method === :welch
            s_pow, s_frq = s_rel_psd(signal, fs=fs, norm=norm, mt=false, f=f)
            title == "default" && (title = "PSD (Welch's periodogram) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mt
            s_pow, s_frq = s_rel_psd(signal, fs=fs, norm=norm, mt=true, f=f)
            title == "default" && (title = "PSD (multi-tapered) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        end
    end

    # set labels
    if type !== :w3d && type !== :s3d && type !== :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if ref !== :abs
            ylabel == "default" && (ylabel = "Power ratio")
        end
        if norm == true
            ylabel == "default" && (ylabel = "Power [dB]")
        else
            ylabel == "default" && (ylabel = "Power [μV^2/Hz]")
        end
    end

    if type === :normal
        # ndims(s_pow) > 1 && throw(ArgumentError("For type=:normal the signal must contain 1 channel."))
        if ndims(s_pow) == 1
            p = plot_psd(s_frq,
                     s_pow,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     title=title,
                     norm=norm,
                     frq_lim=frq_lim,
                     ax=ax,
                     mono=mono;
                     kwargs...)
        else
            p = plot_psd(s_frq,
                     s_pow,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     labels=labels,
                     title=title,
                     norm=norm,
                     frq_lim=frq_lim,
                     ax=ax,
                     mono=mono;
                     kwargs...)
        end
    elseif type === :butterfly
        ndims(s_pow) < 2 && throw(ArgumentError("For type=:butterfly plot the signal must contain ≥ 2 channels."))
        title = replace(title, "channel" => "channels")
        p = plot_psd_butterfly(s_frq,
                               s_pow,
                               labels=labels,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               norm=norm,
                               frq_lim=frq_lim,
                               ax=ax,
                               mono=mono;
                               kwargs...)
    elseif type === :mean
        ndims(s_pow) < 2 && throw(ArgumentError("For type=:mean plot the signal must contain ≥ 2 channels."))
        title = replace(title, "PSD" => "PSD [mean ± 95%CI]")
        title = replace(title, "channel" => "averaged channels")
        p = plot_psd_avg(s_frq,
                         s_pow,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         norm=norm,
                         frq_lim=frq_lim,
                         ax=ax,
                         mono=mono;
                         kwargs...)
    elseif type === :w3d
        ndims(s_pow) < 2 && throw(ArgumentError("For type=:w3d plot the signal must contain ≥ 2 channels."))
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        zlabel == "default" && (zlabel = norm == true ? "Power [dB]" : "Power [μV^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(s_frq,
                        s_pow,
                        labels=labels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        norm=norm,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:w;
                        kwargs...)
    elseif type === :s3d
        ndims(s_pow) < 2 && throw(ArgumentError("For type=:w3d plot the signal must contain ≥ 2 channels."))
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        zlabel == "default" && (zlabel = norm == true ? "Power [dB]" : "Power [μV^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(s_frq,
                        s_pow,
                        labels=labels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        norm=norm,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:s;
                        kwargs...)
    elseif type === :topo
        eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available."))
        ndims(s_pow) == 1 && (s_pow = reshape(s_pow, 1, length(s_pow)))
        xlabel == "default" && (xlabel = "")
        ylabel == "default" && (ylabel = "")
        title = replace(title, "channel" => "channels")
        p = plot_psd_topo(eeg.eeg_locs,
                          s_frq,
                          s_pow,
                          channel=channel,
                          labels=labels,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          norm=norm,
                          frq_lim=frq_lim,
                          ax=ax,
                          mono=mono;
                          kwargs...)
    end

    if typeof(p) == Plots.Plot{Plots.GRBackend}
        Plots.plot(p)
    else
        p
    end

    return p
end

"""
    eeg_plot_psd(eeg::NeuroAnalyzer.EEG, c::Union{Symbol, AbstractArray}; epoch::Int64, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, norm::Bool=true, method::Symbol=:welch, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, kwargs...)

Plot power spectrum density of embedded or external component.

# Arguments

- `eeg::NeuroAnalyzer.EEG`: EEG object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `epoch::Int64`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:welch`: method of calculating PSD:
    - `:welch`: Welch's periodogram
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
- `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to EEG band: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher` 
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
- `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
- `ylabel::String="default"`: y-axis label, default is Power [dB] or Power [μV^2/Hz]
- `zlabel::String="default"`: z-axis label for 3-d plots, default is Power [dB] or Power [μV^2/Hz]
- `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or grey palette
- `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_psd(eeg::NeuroAnalyzer.EEG, c::Union{Symbol, AbstractArray}; epoch::Int64, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, norm::Bool=true, method::Symbol=:welch, nt::Int64=8, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, kwargs...)

    _check_var(type, [:normal, :butterfly, :mean, :w3d, :s3d, :topo], "type")
    _check_var(method, [:welch, :mt, :mw], "method")
    _check_var(ref, [:abs, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "ref")
    _check_var(ax, [:linlin, :loglin, :linlog, :loglog], "ax")
    ref !== :abs && method === :mw && throw(ArgumentError("For relative PSD, method must be :welch or :mt."))
    
    _check_epochs(eeg, epoch)

    # select component c_idxs, default is all c_idxs
    typeof(c) == Symbol && (c = _get_component(eeg, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    labels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (labels = [labels])

    ref !== :abs && (f = eeg_band(eeg, band=ref))

    # get frequency range
    fs = eeg_sr(eeg)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_lim = tuple_order(frq_lim)
    (frq_lim[1] < 0 || frq_lim[2] > fs / 2) && throw(ArgumentError("frq_lim must be ≥ 0 and ≤ $(fs / 2)."))

    # calculate PSD
    signal = c[c_idx, :, epoch]

    # get time vector
    _, t_s1, _, t_s2 = _convert_t(eeg.eeg_epoch_time[1], eeg.eeg_epoch_time[end])

    if ref === :abs
        if method === :welch
            s_pow, s_frq = s_psd(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mt
            s_pow, s_frq = s_psd(signal, fs=fs, norm=norm, mt=true, nt=nt)
            title == "default" && (title = "Absolute PSD (multi-tapered) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mw
            s_pow, s_frq = s_mwpsd(signal, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc)
            title == "default" && (title = "Absolute PSD (Morlet wavelet convolution) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")
        end
    elseif ref === :total
        if method === :welch
            s_pow, s_frq = s_rel_psd(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "PSD (Welch's periodogram) relative to total power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mt
            s_pow, s_frq = s_rel_psd(signal, fs=fs, norm=norm, mt=true)
            title == "default" && (title = "PSD (multi-tapered) relative to total power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")
        end
    else
        if method === :welch
            s_pow, s_frq = s_rel_psd(signal, fs=fs, norm=norm, mt=false, f=f)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mt
            s_pow, s_frq = s_rel_psd(signal, fs=fs, norm=norm, mt=true, f=f)
            title == "default" && (title = "Absolute PSD (multi-tapered) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")
        end
    end

    # set labels
    if type !== :w3d && type !== :s3d
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if norm == true
            ylabel == "default" && (ylabel = "Power [dB]")
        else
            ylabel == "default" && (ylabel = "Power [μV^2/Hz]")
        end
    end

    if type === :normal
        ndims(s_pow) > 1 && throw(ArgumentError("For type=:normal the signal must contain 1 c_idx."))
        p = plot_psd(s_frq,
                     s_pow,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     title=title,
                     norm=norm,
                     frq_lim=frq_lim,
                     ax=ax,
                     mono=mono;
                     kwargs...)
    elseif type === :butterfly
        ndims(s_pow) < 2 && throw(ArgumentError("For type=:butterfly plot the signal must contain ≥ 2 c_idxs."))
        title = replace(title, "component" => "components")
        p = plot_psd_butterfly(s_frq,
                               s_pow,
                               labels=labels,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               norm=norm,
                               frq_lim=frq_lim,
                               ax=ax,
                               mono=mono;
                               kwargs...)
    elseif type === :mean
        ndims(s_pow) < 2 && throw(ArgumentError("For type=:mean plot the signal must contain ≥ 2 c_idxs."))
        title = replace(title, "PSD" => "PSD [mean ± 95%CI]")
        title = replace(title, "component" => "averaged components")
        p = plot_psd_avg(s_frq,
                         s_pow,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         norm=norm,
                         frq_lim=frq_lim,
                         ax=ax,
                         mono=mono;
                         kwargs...)
    elseif type === :w3d
        ndims(s_pow) < 2 && throw(ArgumentError("For type=:w3d plot the signal must contain ≥ 2 channels."))
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        zlabel == "default" && (zlabel = norm == true ? "Power [dB]" : "Power [μV^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(s_frq,
                        s_pow,
                        labels=labels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        norm=norm,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:w;
                        kwargs...)
    elseif type === :s3d
        ndims(s_pow) < 2 && throw(ArgumentError("For type=:w3d plot the signal must contain ≥ 2 channels."))
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        ylabel == "default" && (ylabel = "Channels")
        zlabel == "default" && (zlabel = norm == true ? "Power [dB]" : "Power [μV^2/Hz]")
        title = replace(title, "channel" => "channels")
        p = plot_psd_3d(s_frq,
                        s_pow,
                        labels=labels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        zlabel=zlabel,
                        title=title,
                        norm=norm,
                        frq_lim=frq_lim,
                        ax=ax,
                        mono=mono,
                        variant=:s;
                        kwargs...)
    end

    Plots.plot(p)

    return p
end

"""
    plot_spectrogram(s_t, s_frq, s_pow; <keyword arguments>)

Plot single-channel spectrogram.

# Arguments

- `s_t::Vector{Float64}`: time
- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Array{Float64, 2}`: powers
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(s_t::Vector{Float64}, s_frq::Vector{Float64}, s_pow::Array{Float64, 2}; norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    size(s_pow, 2) == length(s_t) || throw(ArgumentError("Size of powers $(size(s_pow, 2)) and time vector $(length(s_t)) do not match."))
    size(s_pow, 1) == length(s_frq) || throw(ArgumentError("Size of powers $(size(s_pow, 1)) and frequencies vector $(length(s_frq)) do not match."))

    pal = mono == true ? :grays : :darktest
    cb_title = norm == true ? "[dB/Hz]" : "[μV^2/Hz]"

    p = Plots.heatmap(s_t,
                      s_frq,
                      s_pow,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      ylims=frq_lim,
                      xticks=_ticks(s_t),
                      yticks=_ticks(frq_lim),
                      title=title,
                      size=(1200, 800),
                      margins=20Plots.px,
                      seriescolor=pal,
                      colorbar_title=cb_title,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      xtickfontsize=6,
                      ytickfontsize=6;
                      kwargs...)

    return p
end

"""
    plot_spectrogram(s_ch, s_frq, s_pow; <keyword arguments>)

Plot multiple-channel spectrogram.

# Arguments

- `s_ch::Vector{String}`: channel labels
- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Array{Float64, 2}`: powers
- `norm::Bool=true`: whether powers are normalized to dB
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_spectrogram(s_ch::Vector{String}, s_frq::Vector{Float64}, s_pow::Array{Float64, 2}; norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    size(s_pow, 1) == length(s_ch) || throw(ArgumentError("Size of powers $(size(s_pow, 1)) and channels vector $(length(s_ch)) do not match."))
    size(s_pow, 2) == length(s_frq) || throw(ArgumentError("Size of powers $(size(s_pow, 2)) and frequencies vector $(length(s_frq)) do not match."))

    pal = mono == true ? :grays : :darktest
    cb_title = norm == true ? "[dB/Hz]" : "[μV^2/Hz]"
    
    ch = collect(1:length(s_ch)) .- 0.5
    p = Plots.heatmap(s_frq,
                      ch,
                      s_pow,
                      xlabel=xlabel,
                      xticks=_ticks(s_frq),
                      ylabel=ylabel,
                      yticks=(ch, s_ch),
                      title=title,
                      size=(1200, 800),
                      margins=20Plots.px,
                      seriescolor=pal,
                      colorbar_title=cb_title,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      xtickfontsize=6,
                      ytickfontsize=6;
                      kwargs...)

    return p
end

"""
    eeg_plot_spectrogram(eeg; <keyword arguments>)

Plots spectrogram.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `epoch::Int64`: epoch to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or grey palette
- `markers::Bool`: draw markers if available
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_spectrogram(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, AbstractRange}=0, channel::Union{Int64, Vector{Int64}, AbstractRange}, norm::Bool=true, method::Symbol=:standard, nt::Int64=8, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, kwargs...)

    _check_var(method, [:standard, :stft, :mt, :mw], "method")

    _check_epochs(eeg, epoch)
    _check_channels(eeg, channel)

    labels = eeg_labels(eeg)[channel]
    length(channel) == 1 && (labels = [labels])

    # get frequency range
    fs = eeg_sr(eeg)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_lim = tuple_order(frq_lim)
    (frq_lim[1] < 0 || frq_lim[2] > fs / 2) && throw(ArgumentError("frq_lim must be ≥ 0 and ≤ $(fs / 2)."))

    # calculate spectrogram
    signal = eeg.eeg_signals[channel, :, epoch]
    length(channel) > 1 && length(signal) / length(channel) < 4 * eeg_sr(eeg) && throw(ArgumentError("For multi-channel plot, signal length must be ≥ 4 × EEG sampling rate (4 × $(eeg_sr(eeg)) samples)."))

    # get time vector
    _, t_s1, _, t_s2 = _convert_t(eeg.eeg_epoch_time[1], eeg.eeg_epoch_time[end])

    if length(channel) == 1
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        xlabel == "default" && (xlabel = "Time [s]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $(_channel2channel_name(channel)), epoch: $epoch, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f, s_t = s_spectrogram(signal, fs=fs, norm=false, mt=false, st=false, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f, s_t = s_spectrogram(signal, fs=fs, norm=false, mt=true, st=false, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            s_p, s_f, s_t = s_spectrogram(signal, fs=fs, norm=false, mt=false, st=true, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(short-time Fourier transform)")
        elseif method === :mw
            _, s_p, _, s_f = s_wspectrogram(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(s_t, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, kwargs=kwargs)

        # plot markers if available
        # TODO: draw markers length
        if markers == true && eeg.eeg_header[:markers] == true
            markers_pos = eeg.eeg_markers[!, :start] ./ eeg_sr(eeg)
            markers_desc = eeg.eeg_markers[!, :description]
            p = Plots.vline!(markers_pos,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:black,
                             label=false)
            for idx in eachindex(markers_desc)
                p = Plots.plot!(annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
            end
        end

    else
        ylabel == "default" && (ylabel = "Channel")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channels: $(_channel2channel_name(channel)), epoch: $epoch, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=true, nt=nt)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            _info("Method :stft is not available for multi-channel spectrogram, using standard periodogram.")
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mw
            s_p, s_f = s_mwpsd(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_f = linspace(0, frq_lim[2], size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(labels, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, kwargs=kwargs)
    end

    Plots.plot(p)

    return p
end

"""
    eeg_plot_spectrogram(eeg, c; <keyword arguments>)

Plots spectrogram of embedded or external component.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c::Union{Symbol, AbstractArray}`: component to plot
- `epoch::Int64`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `norm::Bool=true`: normalize powers to dB
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[component: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or grey palette
- `markers::Bool`: draw markers if available
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_spectrogram(eeg::NeuroAnalyzer.EEG, c::Union{Symbol, AbstractArray}; epoch::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, AbstractRange}, norm::Bool=true, method::Symbol=:standard, nt::Int64=8, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, kwargs...)

    _check_var(method, [:standard, :stft, :mt, :mw], "method")

    _check_epochs(eeg, epoch)

    # select component c_idxs, default is all c_idxs
    typeof(c) == Symbol && (c = _get_component(eeg, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    labels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (labels = [labels])

    # get frequency range
    fs = eeg_sr(eeg)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_lim = tuple_order(frq_lim)
    (frq_lim[1] < 0 || frq_lim[2] > fs / 2) && throw(ArgumentError("frq_lim must be ≥ 0 and ≤ $(fs / 2)."))

    # calculate spectrogram
    signal = c[c_idx, :, epoch]
    length(c_idx) > 1 && length(signal) / length(c_idx) < 4 * eeg_sr(eeg) && throw(ArgumentError("For multi-channel plot, signal length must be ≥ 4 × EEG sampling rate (4 × $(eeg_sr(eeg)) samples)."))

    # get time vector
    _, t_s1, _, t_s2 = _convert_t(eeg.eeg_epoch_time[1], eeg.eeg_epoch_time[end])

    if length(c_idx) == 1
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        xlabel == "default" && (xlabel = "Time [s]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f, s_t = s_spectrogram(signal, fs=fs, norm=false, mt=false, st=false, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f, s_t = s_spectrogram(signal, fs=fs, norm=false, mt=true, st=false, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            s_p, s_f, s_t = s_spectrogram(signal, fs=fs, norm=false, mt=false, st=true, demean=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[f1:f2, :]
            # s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(short-time Fourier transform)")
        elseif method === :mw
            _, s_p, _, s_f = s_wspectrogram(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_t = linspace(0, (length(signal) / fs), size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(s_t, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, kwargs=kwargs)

        # plot markers if available
        # TODO: draw markers length
        if markers == true && eeg.eeg_header[:markers] == true
            markers_pos = eeg.eeg_markers[!, :start] ./ eeg_sr(eeg)
            markers_desc = eeg.eeg_markers[!, :description]
            p = Plots.vline!(markers_pos,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:black,
                             label=false)
            for idx in eachindex(markers_desc)
                p = Plots.plot!(annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
            end
        end

    else
        ylabel == "default" && (ylabel = "Component")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        title == "default" && (title = "Spectrogram method [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[components: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")

        if method === :standard
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mt
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=true, nt=nt)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            _info("Method :stft is not available for multi-channel spectrogram, using standard periodogram.")
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mw
            s_p, s_f = s_mwpsd(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_f = linspace(0, frq_lim[2], size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
        s_p[s_p .== -Inf] .= minimum(s_p[s_p .!== -Inf])
        p = plot_spectrogram(labels, s_f, s_p, norm=norm, frq_lim=frq_lim, xlabel=xlabel, ylabel=ylabel, title=title, mono=mono, kwargs=kwargs)
    end

    Plots.plot(p)

    return p
end

"""
    plot_electrodes(locs; <keyword arguments>)

Preview of electrode locations. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
- `labels::Bool=true`: plot electrode labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `grid::Bool=false`: draw grid, useful for locating electrode positions
- `plot_size::Int64=400`: plot dimensions in pixels (size × size)

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_electrodes(locs::DataFrame; channel::Union{Int64, Vector{Int64}, AbstractRange}, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head_labels::Bool=true, mono::Bool=false, head_details::Bool=true, grid::Bool=false, plot_size::Int64=400)

    pal = mono == true ? :grays : :darktest

    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    if plot_size > 300
        marker_size = plot_size ÷ 75
        font_size = plot_size ÷ 75
    else
        marker_size = plot_size ÷ 50
        font_size = plot_size ÷ 50
        labels = false
    end

    length(channel) > 64 && (font_size = plot_size ÷ 100)

    if grid == false
        p = Plots.plot(grid=false,
                       framestyle=:none,
                       palette=pal,
                       size=(plot_size, plot_size),
                       border=:none,
                       aspect_ratio=1,
                       right_margin=-30 * Plots.px,
                       bottom_margin=-20 * Plots.px,
                       top_margin=-30 * Plots.px,
                       left_margin=-50 * Plots.px,
                       xlim=(-1.22, 1.23),
                       ylim=(-1.1, 1.2))
    else
        p = Plots.plot(grid=true,
                       palette=pal,
                       size=(plot_size, plot_size),
                       aspect_ratio=1,
                       right_margin=-30 * Plots.px,
                       bottom_margin=-50 * Plots.px,
                       top_margin=-50 * Plots.px,
                       left_margin=-5 * Plots.px,
                       xticks=-1:0.1:1,
                       yticks=-1:0.1:1,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       xlim=(-1.22, 1.23),
                       ylim=(-1.1, 1.2))
    end
    hd = _draw_head(p, head_labels=head_labels, head_details=head_details)
    p = Plots.plot!(hd)

    for idx in eachindex(locs[!, :labels])
        if idx in channel
            if selected != 0
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                #seriestype=:scatter,
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            else
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                #seriestype=:scatter,
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            end
        end
        if idx in selected
            if mono != true
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=idx,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                #seriestype=:scatter,
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            else
                #p = Plots.plot!((loc_x[idx], loc_y[idx]),
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                #seriestype=:scatter,
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            end
        end
    end

    if labels == true
        for idx in eachindex(locs[!, :labels])
            if idx in channel
                Plots.plot!(annotation=(loc_x[idx], loc_y[idx] + 0.075, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
            if idx in selected
                Plots.plot!(annotation=(loc_x[idx], loc_y[idx] + 0.075, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
        end
    end

    Plots.plot(p)

    return p
end

"""
    plot_electrodes3d(locs; <keyword arguments>)

3D interactive preview of electrode locations. It uses spherical :loc_radius_sph, :loc_theta_sph and :loc_phi_sph locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
- `labels::Bool=true`: plot electrode labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `plot_size::Int64=800`: plot dimensions in pixels (plot_size×plot_size)
c

# Returns

- `fig::GLMakie.Figure`
"""
function plot_electrodes3d(locs::DataFrame; channel::Union{Int64, Vector{Int64}, AbstractRange}, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head_labels::Bool=true, mono::Bool=false, plot_size::Int64=800)

    # selected != 0 && length(intersect(channel, selected)) < length(selected) && throw(ArgumentError("channel must include selected."))
    # channel = setdiff(channel, selected)

    pal = mono == true ? :grays : :darktest

    loc_x = zeros(nrow(locs))
    loc_y = zeros(nrow(locs))
    loc_z = zeros(nrow(locs))

    for idx in 1:nrow(locs)
        loc_x[idx], loc_y[idx], loc_z[idx] = sph2cart(locs[idx, :loc_radius_sph], locs[idx, :loc_theta_sph], locs[idx, :loc_phi_sph])
    end

    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    x_lim = (-1.1, 1.1)
    y_lim = (-1.1, 1.1)
    z_lim = extrema(loc_z)

    marker_size = plot_size ÷ 40
    font_size = plot_size ÷ 40

    fig = Figure(; resolution=(plot_size, plot_size))
    ax = Axis3(fig[1, 1]; aspect=(1, 1, 0.5), perspectiveness=0.5, limits = (x_lim, y_lim, z_lim))
    # hidedecorations!(ax, grid=true, ticks=true)

    GLMakie.scatter!(ax, loc_x[channel], loc_y[channel], loc_z[channel], markersize=marker_size, color=:gray)
    if selected != 0
        if mono == true
            GLMakie.scatter!(ax, loc_x[selected], loc_y[selected], loc_z[selected], markersize=marker_size, color=:gray)
        else
            GLMakie.scatter!(ax, loc_x[selected], loc_y[selected], loc_z[selected], markersize=marker_size, color=:red)
        end
    end

    if labels == true
        for idx in eachindex(locs[!, :labels])
            if idx in channel
                GLMakie.text!(ax, locs[!, :labels][idx], position=(loc_x[idx], loc_y[idx], loc_z[idx]), fontsize=font_size)
            end
            if idx in selected
                GLMakie.text!(ax, locs[!, :labels][idx], position=(loc_x[idx], loc_y[idx], loc_z[idx]), fontsize=font_size)
            end
        end
    end

    if head_labels == true
        GLMakie.text!(ax, "Nz", position=(0, 1.025, 0), fontsize = font_size)
        GLMakie.text!(ax, "Iz", position=(0, -1.025, 0), fontsize = font_size)
        GLMakie.text!(ax, "LPA", position=(-1.025, 0, 0), fontsize = font_size)
        GLMakie.text!(ax, "RPA", position=(1.025, 0, 0), fontsize = font_size)
        GLMakie.text!(ax, "top", position=(0, 0, 1.025), fontsize = font_size)
    end
    fig

    return fig
end

"""
    eeg_plot_electrodes(eeg; <keyword arguments>)

Preview of electrode locations.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG/MEG channels
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: which channel should be highlighted
- `labels::Bool=true`: plot electrode labels
- `head::Bool`=true: plot head
- `head_labels::Bool=false`: plot head labels
- `plot_size::Int64=400`: plot dimensions in pixels (plot_size×plot_size)
- `head_details::Bool=true`: draw nose and ears
- `mono::Bool=false`: use color or grey palette
- `threed::Bool=false`: 3-dimensional plot
- `grid::Bool=false`: draw grid, useful for locating electrode positions
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_electrodes(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head::Bool=true, head_labels::Bool=false, plot_size::Int64=400, head_details::Bool=true, mono::Bool=false, threed::Bool=false, grid::Bool=false, kwargs...)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    # select channels, default is all channels
    _check_channels(eeg, channel, Symbol(eeg.eeg_header[:signal_type]))
    selected != 0 && _check_channels(eeg, selected)

    if threed == false
        p = plot_electrodes(eeg.eeg_locs, channel=channel, selected=selected, labels=labels, head_labels=head_labels, mono=mono, head_details=head_details, plot_size=plot_size, grid=grid)
    else
        p = plot_electrodes3d(eeg.eeg_locs, channel=channel, selected=selected, labels=labels, head_labels=head_labels, mono=mono, plot_size=plot_size)
    end

    return p
end

"""
    plot_matrix(m; <keyword arguments>)

Plot matrix.

# Arguments

- `m::Array{<:Real, 2}`
- `xlabels::Vector{String}`
- `ylabels::Vector{String}`
- `xlabel::String=""`
- `ylabel::String=""`
- `title::String=""`
- `cb_title::String=""`: color bar title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_matrix(m::Array{<:Real, 2}; xlabels::Vector{String}, ylabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", cb_title::String="", mono::Bool=false, kwargs...)

    size(m, 1) == size(m, 2) || throw(ArgumentError("Matrix is not square."))
    length(xlabels) == length(ylabels) || throw(ArgumentError("Lengths of xlabels and ylabels must be equal."))
    length(xlabels) == size(m, 1) || throw(ArgumentError("Length of xlabels and matrix size must be equal."))
    length(ylabels) == size(m, 2) || throw(ArgumentError("Length of ylabels and matrix size must be equal."))

    n = size(m, 1)
    r = maximum(length.(xlabels)) > 10 ? 45 : 0
    mar = maximum(length.(xlabels)) > 10 ? 40 : 0
    pal = mono == true ? :grays : :darktest

    p = Plots.heatmap(m,
                      title=title,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      xaxis=(tickfontrotation=r),
                      xticks=(1:n, xlabels),
                      yticks=(1:n, xlabels),
                      seriescolor=pal,
                      colorbar_title=cb_title,
                      size=(1200, 800),
                      left_margin=mar * Plots.px,
                      bottom_margin=mar * Plots.px,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      xtickfontsize=6,
                      ytickfontsize=6;
                      kwargs...)

    return p
end

"""
    plot_covmatrix(m, lags; <keyword arguments>)

Plot cross/auto-covariance matrix.

# Arguments

- `m::Abstractvector`: covariance matrix
- `lags::AbstractVector`: covariance lags, lags will be displayed in ms
- `xlabel::String="lag"`
- `ylabel::String=""`
- `title::String=""`
- `cb_title::String=""`: color bar title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_covmatrix(m::AbstractVector, lags::AbstractVector; xlabel::String="lag [ms]", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest
    r = length(lags) > 10 ? 90 : 0
    p = Plots.plot(lags,
                   m,
                   title=title,
                   xlabel=xlabel,
                   ylabel=ylabel,
                   xticks=round.(lags, digits=1),
                   xaxis=(tickfontrotation=r),
                   yticks=false,
                   palette=pal,
                   size=(600, 200),
                   lw=0.5,
                   grid=false,
                   legend=false,
                   bottom_margin=10 * Plots.px,
                   titlefontsize=5,
                   xlabelfontsize=5,
                   ylabelfontsize=5,
                   xtickfontsize=4,
                   ytickfontsize=4;
                   kwargs...)

    return p
end

"""
    plot_histogram(signal; <keyword arguments>)

Plot histogram.

# Arguments

- `signal::AbstractVector`
- `type::Symbol`: type of histogram: regular (`:hist`) or kernel density (`:kd`)
- `bins::Union{Int64, Symbol, AbstractVector}=(length(signal) ÷ 10)`: histogram bins: number of bins, range or `:sturges`, `:sqrt`, `:rice`, `:scott` or `:fd`)
- `label::String=""`: channel label
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_histogram(signal::AbstractVector; type::Symbol=:hist, bins::Union{Int64, Symbol, AbstractVector}=(length(signal) ÷ 10), label::String="", xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    _check_var(type, [:hist, :kd], "type")

    type === :kd && (type = :density)

    pal = mono == true ? :grays : :darktest

    if mean(signal) < median(signal)
        xticks = [floor(minimum(signal), digits=1), round(mean(signal), digits=1), round(median(signal), digits=1), ceil(maximum(signal), digits=1)]
    else
        xticks = [floor(minimum(signal), digits=1), round(median(signal), digits=1), round(mean(signal), digits=1), ceil(maximum(signal), digits=1)]
    end        

    p = Plots.plot(signal,
                   seriestype=type,
                   xlabel=xlabel,
                   ylabel=ylabel,
                   label=label,
                   title=title,
                   palette=pal,
                   bins=bins,
                   grid=false,
                   linecolor=:black,
                   fillcolor=:grey,
                   fillalpha=0.5,
                   linewidth=1,
                   size=(800, 800),
                   left_margin=30*Plots.px,
                   xaxis=(tickfontrotation=90),
                   margins=10Plots.px,
                   xticks=xticks,
                   yticks=false,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=5,
                   ytickfontsize=5;
                   kwargs...)
    p = Plots.vline!([round(mean(signal), digits=1)], lw=1, ls=:dot, lc=:black, label="mean")
    p = Plots.vline!([round(median(signal), digits=1)], lw=0.5, ls=:dash, lc=:grey, alpha=0.5, label="median")

    Plots.plot(p)

    return p
end

"""
    plot_bar(signal; <keyword arguments>)

Bar plot.

# Arguments

- `signal::AbstractVector`
- `labels::Vector{String}`: x-ticks labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_bar(signal::AbstractVector; labels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    length(signal) == length(labels) || throw(ArgumentError("signal length ($(length(signal))) must be equal to labels length ($(length(labels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :lightblue

    p = Plots.plot(signal,
                   seriestype=:bar,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(labels), labels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)

    Plots.plot(p)

    return p
end

"""
    plot_line(signal; <keyword arguments>)

Line plot.

# Arguments

- `signal::AbstractVector`
- `labels::Vector{String}`: x-ticks labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_line(signal::AbstractVector; labels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    length(signal) == length(labels) || throw(ArgumentError("signal length ($(length(signal))) must be equal to x-ticks length ($(length(labels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(signal,
                   seriestype=:line,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(labels), labels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)

    Plots.plot(p)

    return p
end

"""
    plot_line(signal; <keyword arguments>)

Line plot.

# Arguments

- `signal::AbstractArray`
- `labels::Vector{String}`: signal rows labels
- `xlabels::Vector{String}`: x-ticks labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_line(signal::AbstractArray; labels::Vector{String}, xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    ndims(signal) != 2 && throw(ArgumentError("signal must have 2-dimensions."))
    size(signal, 1) == length(labels) || throw(ArgumentError("Number of signal columns ($(size(signal, 1))) must be equal to labels length ($(length(labels)))."))
    size(signal, 2) == length(xlabels) || throw(ArgumentError("Number of signal columns ($(size(signal, 2))) must be equal to x-ticks length ($(length(xlabels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(signal[1, :],
                   seriestype=:line,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=:topright,
                   label=labels[1],
                   xticks=(1:length(xlabels), xlabels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)
    for idx in 2:size(signal, 1)
        p = Plots.plot!(signal[idx, :],
                        seriestype=:line,
                        color=idx,
                        label=labels[idx])
    end
    Plots.plot(p)

    return p
end

"""
    plot_box(signal; <keyword arguments>)

Box plot.

# Arguments

- `signal::AbstractArray`
- `labels::Vector{String}`: group labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_box(signal::AbstractArray; labels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    ndims(signal) != 2 && throw(ArgumentError("signal must have 2-dimensions."))
    size(signal, 1) == length(labels) || throw(ArgumentError("Number of signal columns ($(size(signal, 1))) must be equal to x-ticks length ($(length(xlabels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(signal',
                   seriestype=:box,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(labels), labels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)
    Plots.plot(p)

    return p
end

"""
    plot_violin(signal; <keyword arguments>)

Violin plot.

# Arguments

- `signal::AbstractArray`
- `labels::Vector{String}`: group labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_violin(signal::AbstractArray; labels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    ndims(signal) != 2 && throw(ArgumentError("signal must have 2-dimensions."))
    size(signal, 1) == length(labels) || throw(ArgumentError("Number of signal columns ($(size(signal, 1))) must be equal to x-ticks length ($(length(xlabels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(signal',
                   seriestype=:violin,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(labels), labels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)
    Plots.plot(p)

    return p
end

"""
    plot_dots(signal; <keyword arguments>)

Dots plot.

# Arguments

- `signal::Vector{Vector{Float64}}`
- `labels::Vector{String}`: group labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_dots(signal::Vector{Vector{Float64}}; labels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    size(signal, 1) == length(labels) || throw(ArgumentError("Number of signal columns ($(size(signal, 1))) must be equal to x-ticks length ($(length(xlabels)))."))

    pal = mono == true ? :grays : :darktest

    p = Plots.plot(size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(labels), labels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs...)
    for idx1 in eachindex(labels)
        for idx2 in eachindex(signal[idx1])
            if mono == false
                p = Plots.scatter!((idx1, signal[idx1][idx2]),
                                   color=idx1)
            else
                p = Plots.scatter!((idx1, signal[idx1][idx2]),
                                   color=:black)
            end
        end
    end

    Plots.plot(p)

    return p
end

"""
    plot_paired(signal; <keyword arguments>)

Plot paired data.

# Arguments

- `signal::Vector{Vector{Float64}}`
- `labels::Vector{String}`: group labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_paired(signal::Vector{Vector{Float64}}; labels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    size(signal, 1) == length(labels) || throw(ArgumentError("Number of signal columns ($(size(signal, 1))) must be equal to x-ticks length ($(length(xlabels)))."))
    ll = Vector{Int64}()
    for idx in eachindex(labels)
        push!(ll, length(signal[idx]))
    end
    length(unique(ll)) == 1 || throw(ArgumentError("Each group must have the same number of values."))

    pal = mono == true ? :grays : :darktest

    p = Plots.plot(size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(labels), labels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs...)
    for idx1 in eachindex(signal[1])
        c_tmp = zeros(length(labels))
        for idx2 in eachindex(labels)
            c_tmp[idx2] = signal[idx2][idx1]
        end
        p = Plots.plot!(c_tmp,
                        color=:black)
    end
    for idx1 in eachindex(labels)
        for idx2 in eachindex(signal[idx1])
            if mono == false
                p = Plots.scatter!((idx1, signal[idx1][idx2]),
                                   color=idx1)
            else
                p = Plots.scatter!((idx1, signal[idx1][idx2]),
                                   color=:black)
            end
        end
    end

    Plots.plot(p)

    return p
end

"""
    plot_polar(signal; <keyword arguments>)

Polar plot.

# Arguments

- `signal::Union{AbstractVector, AbstractArray}`
- `m::Tuple{Real, Real}=(0, 0)`: major value to plot
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_polar(signal::Union{AbstractVector, AbstractArray}; m::Tuple{Real, Real}=(0, 0), title::String="", mono::Bool=false, kwargs...)

    length(m) > 2 && throw(ArgumentError("m must have exactly 2 values: phases and lengths."))
    ndims(signal) > 1 && size(signal, 2) > 2 && throw(ArgumentError("signal must have exactly 2 columns: phases and lengths."))

    pal = mono == true ? :grays : :darktest

    if ndims(signal) == 1
        p = Plots.plot([0, signal[1]], [0, 1],
                       size=(800, 800),
                       projection=:polar,
                       left_margin=50Plots.px,
                       right_margin=50Plots.px,
                       bottom_margin=30Plots.px,
                       legend=false,
                       xticks=false,
                       yticks=false, 
                       title=title,
                       color=:black,
                       palette=pal,
                       linewidth=0.2,
                       titlefontsize=8,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       kwargs...)
        for idx in 2:length(signal)
            p = Plots.plot!([0, signal[idx]], [0, 1],
                            projection=:polar,
                            color=:black)
        end
    else
        p = Plots.plot([0, signal[1, 1]], [0, signal[1, 2]],
                       size=(800, 800),
                       projection=:polar,
                       left_margin=50Plots.px,
                       right_margin=50Plots.px,
                       bottom_margin=30Plots.px,
                       legend=false,
                       xticks=false,
                       yticks=false, 
                       title=title,
                       color=:black,
                       palette=pal,
                       linewidth=0.2,
                       titlefontsize=8,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       kwargs...)
        for idx in 2:size(signal, 1)
            p = Plots.plot!([0, signal[idx, 1]], [0, signal[idx, 2]],
                            projection=:polar,
                            color=:black)
        end
    end
    if m != (0, 0)
        p = Plots.plot!([0, m[1]], [0, m[2]],
                        lw=2,
                        projection=:polar,
                        color=:red)
    end

    Plots.plot(p)

    return p
end

"""
    plot_filter_response(<keyword arguments>)

Plot filter response.

# Arguments

- `fs::Int64`: sampling rate
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`: second-order IIR notch filter
    - `:remez`: Remez FIR filter
- `ftype::Union{Symbol, Nothing}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
- `n::Int64`: signal length in samples
- `fs::Int64`: sampling rate
- `order::Int64=8`: filter order (6 dB/octave), number of taps for `:remez`, attenuation (× 4 dB) for `:fir` filters
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Real=-1`: bandwidth for `:iirnotch` and :remez filters
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter; default is Hamming window, number of taps is calculated using fred harris' rule-of-thumb
- `mono::Bool=false`: use color or grey palette
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_filter_response(; fs::Int64, n::Int64, fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple}, order::Int64=8, rp::Real=8, rs::Real=-1, bw::Real=-1, window::Union{Vector{Float64}, Nothing}=nothing, mono::Bool=false, frq_lim::Tuple{Real, Real}=(0, 0), kwargs...)

    pal = mono == true ? :grays : :darktest

    frq_lim == (0, 0) && (frq_lim = (0, fs / 2))
    frq_lim = tuple_order(frq_lim)

    flt = s_filter_create(fprototype=fprototype, ftype=ftype, cutoff=cutoff, n=n, fs=fs, order=order, rp=rp, rs=rs, bw=bw, window=window)

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch]
        H, w = freqresp(flt)
        # convert to dB
        H = 20 * log10.(abs.(H))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        
        if fprototype !== :iirnotch
            fname = titlecase(String(fprototype))
            title = "Filter: $(fname), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order\n\nFrequency response"
        else
            fname = "IIR notch"
            title = "Filter: $(fname), cutoff: $cutoff Hz, band width: $bw\n\nFrequency response"
        end

        p1 = Plots.plot(w,
                        H,
                        title=title,
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylims=(-100, 10),
                        ylabel="Magnitude\n[dB]",
                        xlabel="Frequency [Hz]",
                        label="",
                        top_margin=10Plots.px,
                        bottom_margin=10Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            p1 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p1 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p1 = Plots.plot!((0, cutoff[2]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:green)
        end

        phi, w = phaseresp(flt)
        phi = rad2deg.(angle.(phi))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        p2 = Plots.plot(w,
                        phi,
                        title="Phase response",
                        ylims=(-180, 180),
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylabel="Phase\n[°]",
                        xlabel="Frequency [Hz]",
                        label="",
                        bottom_margin=10Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            p2 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p2 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p2 = Plots.plot!((0, cutoff[2]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:green)
        end

        tau, w = grpdelay(flt)
        tau = abs.(tau)
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        p3 = Plots.plot(w,
                        tau,
                        title="Group delay",
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylabel="Group delay\n[samples]",
                        xlabel="Frequency [Hz]",
                        label="",
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            p3 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p3 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p3 = Plots.plot!((0, cutoff[2]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:green)
        end

        p = Plots.plot(p1, p2, p3, size=(1200, 800), margins=20Plots.px, layout=(3, 1), palette=pal; kwargs...)
    else
        @show "FIR"
        w = range(0, stop=pi, length=1024)
        H = _fir_response(flt, w)
        # convert to dB
        H = 20 * log10.(abs.(H))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        if fprototype === :fir
            title = "Filter: FIR, type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, taps: $(length(flt)), attenuation: $(order * 15) dB\nFrequency response"
        elseif fprototype === :remez
            title = "Filter: Remez, type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, taps: $(length(flt))\nFrequency response"
        end
        p1 = Plots.plot(w,
                        H,
                        title=title,
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylims=(-100, 10),
                        ylabel="Magnitude\n[dB]",
                        xlabel="Frequency [Hz]",
                        label="",
                        bottom_margin=10Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            p1 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p1 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p1 = Plots.plot!((0, cutoff[2]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:green)
        end
        w = range(0, stop=pi, length=1024)
        phi = _fir_response(flt, w)
        phi = DSP.unwrap(-atan.(imag(phi), real(phi)))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        p2 = Plots.plot(w,
                        phi,
                        title="Phase response",
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylabel="Phase\n[rad]",
                        xlabel="Frequency [Hz]",
                        label="",
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            p2 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p2 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p2 = Plots.plot!((0, cutoff[2]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:green)
        end

        p = Plots.plot(p1, p2, size=(1200, 800), margins=20*Plots.px, layout=(2, 1), palette=pal; kwargs...)
    end

    return p
end

"""
    plot_weights(locs; <keyword arguments>)

Plot weights at electrode positions. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
- `weights::Vector{<:Real}=[]`: weights vector
- `labels::Bool=true`: plot electrode labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `plot_size::Int64=400`: plot dimensions in pixels (size × size)

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_weights(locs::DataFrame; channel::Union{Int64, Vector{Int64}, AbstractRange}, weights::Vector{<:Real}=[], labels::Bool=true, head_labels::Bool=true, mono::Bool=false, head_details::Bool=true, plot_size::Int64=400)

    length(weights) > length(channel) && throw(ArgumentError("Number of weights must be ≤ number of channels to plot ($(length(channel)))."))
    length(weights) < 1 && throw(ArgumentError("weights must contain at least one value."))

    # selected != 0 && length(intersect(channel, selected)) < length(selected) && throw(ArgumentError("channel must include selected."))
    # channel = setdiff(channel, selected)

    pal = mono == true ? :grays : :darktest

    marker_size = plot_size ÷ 100
    font_size = plot_size ÷ 100

    p = Plots.plot(grid=true,
                   framestyle=:none,
                   palette=pal,
                   size=(plot_size, plot_size),
                   markerstrokewidth=0,
                   border=:none,
                   aspect_ratio=1,
                   margins=-plot_size * Plots.px,
                   titlefontsize=plot_size ÷ 50)

    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    for idx in eachindex(locs[!, :labels])
        if idx in channel
            p = Plots.plot!((loc_x[idx], loc_y[idx]),
                            color=:black,
                            seriestype=:scatter,
                            grid=true,
                            label="",
                            markersize=marker_size,
                            markerstrokewidth=0,
                            markerstrokealpha=0)
        end
    end

    for idx in eachindex(locs[!, :labels])
        if idx in channel
            Plots.plot!(annotation=(loc_x[idx], loc_y[idx] + 0.05, Plots.text(string(weights[idx]), pointsize=font_size)))
        end
    end

    hd = _draw_head(p, head_labels=head_labels, head_details=head_details)
    p = Plots.plot!(hd)

    Plots.plot(p)

    return p
end

"""
    eeg_plot_weights(eeg; <keyword arguments>)

Plot weights at electrode positions. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG/MEG channels
- `weights::Matrix{<:Real}`: matrix of weights
- `labels::Bool=false`: plot electrode labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `title::String=""`: plot title
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_weights(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), weights::Vector{<:Real}, labels::Bool=true, head_labels::Bool=false, mono::Bool=false, head_details::Bool=true, plot_size::Int64=800, title::String="", kwargs...)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    # remove non-EEG/MEG channels
    eeg_tmp = deepcopy(eeg)
    eeg_keep_channel_type!(eeg_tmp, type=Symbol(eeg_tmp.eeg_header[:signal_type]))

    _check_channels(eeg, channel, Symbol(eeg.eeg_header[:signal_type]))
    typeof(channel) == Int64 && throw(ArgumentError("≥ 2 channels are required."))

    p = plot_weights(eeg_tmp.eeg_locs, weights=weights, channel=channel, labels=labels, head_labels=head_labels, mono=mono, plot_size=plot_size, head_details=head_details)

    Plots.plot!(p, title=title; kwargs)

    return p
end

"""
    eeg_plot_connections(eeg; <keyword arguments>)

Plot weights at electrode positions. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG/MEG channels
- `connections::Matrix{<:Real}`: matrix of connections weights
- `threshold::Real`: plot all connection above threshold
- `threshold_type::Symbol=:g`: rule for thresholding: = (`:eq`), ≥ (`:geq`), ≤ (`:leq`), > (`:g`), < (`:l`)
- `weights::Bool=true`: weight line widths and alpha based on connection value
- `labels::Bool=false`: plot electrode labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `title::String=""`: plot title
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_connections(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), connections::Matrix{<:Real}, threshold::Real, threshold_type::Symbol=:g, weights::Bool=true, labels::Bool=true, head_labels::Bool=false, mono::Bool=false, head_details::Bool=true, plot_size::Int64=800, title::String="", kwargs...)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    _check_var(threshold_type, [:eq, :geq, :leq, :g, :l], "threshold_type")

    # remove non-EEG/MEG channels
    eeg_tmp = deepcopy(eeg)
    eeg_keep_channel_type!(eeg_tmp, type=Symbol(eeg_tmp.eeg_header[:signal_type]))

    _check_channels(eeg, channel, Symbol(eeg.eeg_header[:signal_type]))
    typeof(channel) == Int64 && throw(ArgumentError("≥ 2 channels are required."))

    p = plot_connections(eeg_tmp.eeg_locs, connections=connections, channel=channel, threshold=threshold, threshold_type=threshold_type, weights=weights, labels=labels, head_labels=head_labels, mono=mono, plot_size=plot_size, head_details=head_details)

    Plots.plot!(p, title=title; kwargs)

    return p
end

"""
    eeg_plot_connections(eeg; <keyword arguments>)

Plot connections between channels. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Vector{Int64}, AbstractRange}`: channel(s) to plot
- `connections::Matrix{<:Real}`: matrix of connections weights
- `threshold::Real`: plot all connection above threshold
- `threshold_type::Symbol=:g`: rule for thresholding: = (`:eq`), ≥ (`:geq`), ≤ (`:leq`), > (`:g`), < (`:l`)
- `weights::Bool=true`: weight line widths and alpha based on connection value
- `labels::Bool=false`: plot electrode labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
- `head_details::Bool=true`: draw nose and ears
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_connections(locs::DataFrame; channel::Union{Vector{Int64}, AbstractRange}, connections::Matrix{<:Real}, threshold::Real, threshold_type::Symbol=:g, weights::Bool=true, labels::Bool=true, head_labels::Bool=false, mono::Bool=false, head_details::Bool=true, plot_size::Int64=800, kwargs...)

    size(connections, 1) == length(channel) || throw(ArgumentError("Length of channel and number of connections rows must be equal."))
    _check_var(threshold_type, [:eq, :geq, :leq, :g, :l], "threshold_type")
    pal = mono == true ? :grays : :darktest
    marker_size = plot_size ÷ 100
    font_size = plot_size ÷ 100

    p = Plots.plot(grid=true,
                   framestyle=:none,
                   palette=pal,
                   size=(plot_size, plot_size),
                   markerstrokewidth=0,
                   border=:none,
                   aspect_ratio=1,
                   margins=-plot_size * Plots.px,
                   titlefontsize=plot_size ÷ 50)

    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    for idx in eachindex(locs[!, :labels])
        if idx in channel
            p = Plots.plot!((loc_x[idx], loc_y[idx]),
                            color=:gray,
                            seriestype=:scatter,
                            grid=true,
                            label="",
                            markersize=marker_size,
                            markerstrokewidth=0,
                            markerstrokealpha=0)
        end
    end

    if labels == true
        for idx in eachindex(locs[!, :labels])
            if idx in channel
                Plots.plot!(annotation=(loc_x[idx], loc_y[idx] + 0.05, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
        end
    end

    hd = _draw_head(p, head_labels=head_labels, head_details=head_details)
    p = Plots.plot!(hd)

    m_tmp = s_normalize_n(connections)

    for idx1 in 1:size(connections, 1)
        for idx2 in 1:size(connections, 1)
            if threshold_type === :g
                if connections[idx1, idx2] > threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            elseif threshold_type === :l
                if connections[idx1, idx2] < threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            elseif threshold_type === :eq
                if connections[idx1, idx2] == threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            elseif threshold_type === :leq
                if connections[idx1, idx2] <= threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            elseif threshold_type === :geq
                if connections[idx1, idx2] >= threshold
                    if weights == true
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=6 * m_tmp[idx1, idx2], alpha=0.25 * m_tmp[idx1, idx2], lc=:black, legend=false)
                    else
                        p = Plots.plot!([loc_x[idx1], loc_x[idx2]], [loc_y[idx1], loc_y[idx2]], lw=0.2, lc=:black, legend=false)
                    end
                end
            end
        end
    end

    return p
end

"""
    plot_topo(c; <keyword arguments>)

Plot topographical view.

# Arguments

- `signal::Vector{<:Real}`: values to plot (one value per channel)
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `cb::Bool=true`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `s_normalize()`
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `head_labels::Bool=false`: plot head labels
- `head_details::Bool=true`: draw nose and ears
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_topo(signal::Vector{<:Real}; channel::Union{Int64, Vector{Int64}, AbstractRange}, locs::DataFrame, cb::Bool=true, cb_label::String="[A.U.]", title::String="default", mono::Bool=false, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, plot_size::Int64=800, head_labels::Bool=false, head_details::Bool=true, kwargs...)
    
    pal = mono == true ? :grays : :darktest
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")

    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = loc_x[channel]
    loc_y = loc_y[channel]
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    s_interpolated, interpolated_x, interpolated_y = _interpolate(signal, loc_x, loc_y, 100, imethod, nmethod)

    p = Plots.plot(grid=true,
                   framestyle=:none,
                   palette=pal,
                   size=(plot_size, plot_size),
                   border=:none,
                   aspect_ratio=1,
                   left_margin=-20 * Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=6,
                   ylabelfontsize=6,
                   xtickfontsize=4,
                   ytickfontsize=4,
                   title=title;
                   kwargs...)

    p = Plots.plot!(interpolated_x,
                    interpolated_y,
                    s_interpolated,
                    fill=:darktest,
                    seriestype=:heatmap,
                    seriescolor=pal,
                    colorbar=cb,
                    colorbar_title=cb_label,
                    levels=10,
                    linewidth=0)
    if plot_contours
        p = Plots.plot!(interpolated_x,
                        interpolated_y,
                        s_interpolated,
                        fill=:darktest,
                        seriestype=:contour,
                        seriescolor=pal,
                        colorbar=cb,
                        colorbar_title=cb_label,
                        levels=5,
                        linecolor=:black,
                        linewidth=0.2)
    end
    if plot_electrodes
        p = Plots.plot!((loc_x, loc_y),
                        color=:black,
                        seriestype=:scatter,
                        grid=true,
                        label="",
                        markersize=2,
                        markeralpha=0.5,
                        markerstrokewidth=0,
                        markerstrokealpha=0)
    end

    # draw head
    hd = _draw_head(p, head_labels=head_labels, head_details=head_details, topo=true)
    p = Plots.plot!(hd)

    return p

end

"""
    eeg_plot_topo(eeg; <keyword arguments>)

Topographical plot.

# Arguments

- `eeg::NeuroAnalyzer.EEG`: EEG object
- `epoch::Union{Int64, AbstractRange}=0`: epoch to display
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG/MEG channels
- `segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
- `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or grey palette
- `cb::Bool=true`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `s_normalize()`
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `head_labels::Bool=false`: plot head labels
- `head_details::Bool=true`: draw nose and ears
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_topo(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, AbstractRange}=0, channel::Union{Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg)), title::String="default", mono::Bool=false, cb::Bool=true, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, plot_size::Int64=800, head_labels::Bool=false, head_details::Bool=true, kwargs...)

    eeg_signal_len(eeg) < 10 * eeg_sr(eeg) && segment == (1, 10*eeg_sr(eeg)) && (segment=(1, eeg_signal_len(eeg)))

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(amethod, [:mean, :median], "amethod")
    _check_segment(eeg, segment[1], segment[2])

    if epoch != 0
        _check_epochs(eeg, epoch)
        if eeg_epoch_n(eeg) == 1
            epoch = 0
        else
            segment = (((epoch[1] - 1) * eeg_epoch_len(eeg) + 1), segment[2])
            if typeof(epoch) == Int64
                segment = (segment[1], (segment[1] + eeg_epoch_len(eeg) - 1))
            else
                segment = (segment[1], (epoch[end] * eeg_epoch_len(eeg)))
            end
            epoch = 0
        end
    end

    # remove non-EEG/MEG channels
    eeg_tmp = deepcopy(eeg)
    eeg_keep_channel_type!(eeg_tmp, type=Symbol(eeg_tmp.eeg_header[:signal_type]))

    length(channel) < 2 && throw(ArgumentError("eeg_plot_topo() requires ≥ 2 channels."))
    _check_channels(eeg_tmp, channel)

    length(channel) > nrow(eeg_tmp.eeg_locs) && throw(ArgumentError("Some channels do not have locations."))

    # get time vector
    if segment[2] <= eeg_epoch_len(eeg_tmp)
        signal = eeg_tmp.eeg_signals[channel, segment[1]:segment[2], 1]
    else
        signal = eeg_epoch(eeg_tmp, epoch_n=1).eeg_signals[channel, segment[1]:segment[2], 1]
    end
    t = _get_t(segment[1], segment[2], eeg_sr(eeg_tmp))
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    epoch = _s2epoch(eeg_tmp, segment[1], segment[2])
    
    # average signal and convert to vector
    if size(signal, 2) > 1
        if amethod === :mean
            signal = vec(mean(signal, dims=2))
        elseif amethod === :median
            signal = vec(median(signal, dims=2))
        end
    else
        signal = vec(signal)
    end

    if segment[2] != segment[1] + 1
        title == "default" && (title = "Amplitude topographical plot\n[channel$(_pl(length(channel))): $(_channel2channel_name(channel)), epoch$(_pl(length(epoch))): $epoch, averaged ($(string(amethod))) over time window: $t_s1:$t_s2]")
    else
        title == "default" && (title = "Amplitude topographical plot\n[channel$(_pl(length(channel))): $(_channel2channel_name(channel)), epoch$(_pl(length(epoch))): $epoch, time point: $t_s1]")
    end
    cb_label == "default" && (cb_label = "[A.U.]")

    p = plot_topo(signal, channel=channel, locs=eeg_tmp.eeg_locs, cb=cb, cb_label=cb_label, title=title, mono=mono, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, plot_size=plot_size, head_labels=head_labels, head_details=head_details, kwargs=kwargs)

    Plots.plot(p)

    return p
end

"""
    eeg_plot_topo(eeg; <keyword arguments>)

Topographical plot of embedded or external component.

# Arguments

- `eeg::NeuroAnalyzer.EEG`: EEG object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `epoch::Union{Int64, AbstractRange}=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
- `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or grey palette
- `cb::Bool=true`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `s_normalize()`
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `head_labels::Bool=false`: plot head labels
- `head_details::Bool=true`: draw nose and ears
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_topo(eeg::NeuroAnalyzer.EEG, c::Union{Symbol, AbstractArray}; epoch::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg)), title::String="default", mono::Bool=false, cb::Bool=true, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, plot_size::Int64=800, head_labels::Bool=false, head_details::Bool=true, kwargs...)

    eeg_signal_len(eeg) < 10 * eeg_sr(eeg) && segment == (1, 10*eeg_sr(eeg)) && (segment=(1, eeg_signal_len(eeg)))

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(amethod, [:mean, :median], "amethod")

    no_timepoint = false
    if typeof(c) == Matrix{Float64}
        c = reshape(c, size(c, 1), size(c, 2), 1)
        segment = (1, size(c, 2))
        no_timepoint = true
    elseif typeof(c) == Vector{Float64}
        c = reshape(c, length(c), 1, 1)
        segment = (1, 1)
        no_timepoint = true
    elseif size(c, 2) == 1
        no_timepoint = true
        segment = (1, 1)
    end

    _check_segment(eeg, segment[1], segment[2])

    if epoch != 0
        _check_epochs(eeg, epoch)
        if eeg_epoch_n(eeg) == 1
            epoch = 0
        else
            segment = (((epoch[1] - 1) * eeg_epoch_len(eeg) + 1), segment[2])
            if typeof(epoch) == Int64
                segment = (segment[1], (segment[1] + eeg_epoch_len(eeg) - 1))
            else
                segment = (segment[1], (epoch[end] * eeg_epoch_len(eeg)))
            end
            epoch = 0
        end
    end

    # remove non-EEG/MEG channels
    eeg_tmp = deepcopy(eeg)
    eeg_keep_channel_type!(eeg_tmp, type=Symbol(eeg_tmp.eeg_header[:signal_type]))

    # select component channels, default is all channels
    typeof(c) == Symbol && (c = _get_component(eeg_tmp, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    labels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (labels = [labels])

    length(c_idx) < 2 && throw(ArgumentError("eeg_plot_topo() requires ≥ 2 channels."))
    length(channel) > nrow(eeg_tmp.eeg_locs) && throw(ArgumentError("Some channels do not have locations."))

    # get time vector
    if segment[2] <= eeg_epoch_len(eeg_tmp)
        signal = c[c_idx, segment[1]:segment[2], 1]
    else
        signal = _make_epochs(c, epoch_n=1)[c_idx, segment[1]:segment[2], 1]
    end
    if segment[1] != segment[2]
        t = _get_t(segment[1], segment[2], eeg_sr(eeg_tmp))
    else
        t = _get_t(segment[1], segment[2] + 1, eeg_sr(eeg_tmp))
    end
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    epoch = _s2epoch(eeg_tmp, segment[1], segment[2])
    
    # average signal and convert to vector
    if size(signal, 2) > 1
        if amethod === :mean
            signal = vec(mean(signal, dims=2))
        elseif amethod === :median
            signal = vec(median(signal, dims=2))
        end
    else
        signal = vec(signal)
    end

    if segment[2] != segment[1]
        if no_timepoint != true
            title == "default" && (title = "Amplitude topographical plot\n[component$(_pl(length(c_idx))): $(_channel2channel_name(c_idx)), epoch$(_pl(length(epoch))): $epoch, averaged ($(string(amethod))) over time window: $t_s1:$t_s2]")
        else
            title == "default" && (title = "Amplitude topographical plot\n[component$(_pl(length(c_idx))): $(_channel2channel_name(c_idx)), averaged ($(string(amethod))) over $(length(segment[1]:segment[2])) time point$(_pl(length(segment)))]")
        end
    else
        if no_timepoint != true
            title == "default" && (title = "Amplitude topographical plot\n[component$(_pl(length(c_idx))): $(_channel2channel_name(c_idx)), epoch$(_pl(length(epoch))): $epoch, time point: $t_s1]")
        else
            title == "default" && (title = "Amplitude topographical plot\n[component$(_pl(length(c_idx))): $(_channel2channel_name(c_idx)), $(size(c, 2)) time point$(_pl(size(c, 2)))]")
        end
    end
    cb_label == "default" && (cb_label = "[A.U.]")

    p = plot_topo(signal, channel=c_idx, locs=eeg_tmp.eeg_locs, cb=cb, cb_label=cb_label, title=title, mono=mono, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, plot_size=plot_size, head_labels=head_labels, head_details=head_details, kwargs=kwargs)

    Plots.plot(p)

    return p
end

"""
    plot_compose(p; <keyword arguments>)

Compose a complex plot of various plots contained in vector `p` using layout `layout`. Layout scheme is:
- `(2, 2)`: 2 × 2 plots, regular layout
- `@layout [a{0.2w} b{0.8w};_ c{0.6}]`: complex layout using Plots.jl `@layout` macro

# Arguments

- `p::Vector{Plots.Plot{Plots.GRBackend}}`: vector of plots
- `layout::Union(Matrix{Any}, Tuple{Int64, Int64}}`: layout
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for `p` vector plots

# Returns

- `pc::Plots.Plot{Plots.GRBackend}`
"""
function plot_compose(p::Vector{Plots.Plot{Plots.GRBackend}}; layout::Union{Matrix{Any}, Tuple{Int64, Int64}}, mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest
    if typeof(layout) == Tuple{Int64, Int64} && length(p) < layout[1] * layout[2]
        for idx in 1:(layout[1] * layout[2]) - length(p)
            push!(p, Plots.plot(border=:none, title=""))
        end
    end
    pc = Plots.plot(grid=false,
                    framestyle=:none,
                    border=:none,
                    margins=0Plots.px)
    pc = Plots.plot!(p..., layout=layout, palette=pal; kwargs...)
    Plots.plot(pc)

    return pc
end

"""
    plot_empty()

Return an empty plot, useful for filling matrices of plots. 

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_empty()
    return Plots.plot(grid=false, border=:none, title="")
end

"""
    plot_erp(t, signal, bad; <keyword arguments>)

Plot ERP.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::AbstractVector`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp(t::Union{AbstractVector, AbstractRange}, signal::AbstractVector; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get limits
    ylim = (floor(minimum(signal) * 1.1, digits=0), ceil(maximum(signal) * 1.1, digits=0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_erpticks(t), string.(_erpticks(t) .* 1000)),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)
    
    # plot 0 h-line
    p = Plots.hline!([0],
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot ERP
    p = Plots.plot!(t,
                    signal,
                    linewidth=1,
                    label="",
                    color=:black)

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    return p
end

"""
    plot_erp_avg(t, signal; <keyword arguments>)

Plot ERP amplitude mean and ±95% CI.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::AbstractArray`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_avg(t::Union{AbstractVector, AbstractRange}, signal::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = s_msci95(signal')

    # get limits
    ylim = (floor(minimum(s_l) * 1.1, digits=0), ceil(maximum(s_u) * 1.1, digits=0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_erpticks(t), string.(_erpticks(t) .* 1000)),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot 0 h-line
    p = Plots.hline!([0],
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot upper 95% CI
    p = Plots.plot!(t,
                    s_u,
                    fillrange=s_l,
                    fillalpha=0.35, 
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot lower 95% CI
    p = Plots.plot!(t,
                    s_l,
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot mean
    p = Plots.plot!(t,
                    s_m,
                    label=false,
                    t=:line,
                    c=:black,
                    lw=0.5)

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    return p
end

"""
    plot_erp_butterfly(t, signal; <keyword arguments>)

Butterfly plot of ERP.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::AbstractArray`: data to plot
- `labels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `avg::Bool=false`: plot average ERP
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_butterfly(t::Union{AbstractVector, AbstractRange}, signal::AbstractArray; labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, avg::Bool=true, kwargs...)

    pal = mono == true ? :grays : :darktest

    channel_n = size(signal, 1)

    # get limits
    ylim = (floor(minimum(signal) * 1.1, digits=0), ceil(maximum(signal) * 1.1, digits=0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # plot channels
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_erpticks(t), string.(_erpticks(t) .* 1000)),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot 0 h-line
    p = Plots.hline!([0],
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot signals
    for idx in 1:channel_n
        if labels == [""]
            p = Plots.plot!(t,
                            signal[idx, :],
                            t=:line,
                            linecolor=idx,
                            linewidth=0.2,
                            alpha=0.2,
                            legend=false)
        else
            if labels == repeat([""], channel_n)
                p = Plots.plot!(t,
                                signal[idx, :],
                                t=:line,
                                legend=false,
                                linecolor=idx,
                                linewidth=0.5,
                                alpha=0.5)
            else
                p = Plots.plot!(t,
                                signal[idx, :],
                                t=:line,
                                label=labels[idx],
                                linecolor=idx,
                                linewidth=0.5,
                                alpha=0.5)
            end
        end
    end

    # plot averaged ERP
    if avg == true
        if channel_n == 1
            signal = mean(signal, dims=2)[:]
        else
            signal = mean(signal, dims=1)[:]
        end
        p = Plots.plot!(t,
                        signal,
                        linewidth=1,
                        linecolor=:black,
                        label=false)
    end

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    return p
end

"""
    plot_erp_topo(locs, t, erp; <keyword arguments>)

Plot topographical map ERPs. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `t::Vector{Float64}`: time vector
- `signal::Array{Float64, 2}`: ERPs
- `channels::Union{Vector{Int64}, AbstractRange}`: which channels to plot
- `labels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `fig::GLMakie.Figure`
"""
function plot_erp_topo(locs::DataFrame, t::Vector{Float64}, signal::Array{Float64, 2}; channel=Union{Vector{Int64}, AbstractRange}, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    size(signal, 2) == length(t) || throw(ArgumentError("Length of powers vector must equal length of frequencies vector."))
    length(channel) > nrow(locs) && throw(ArgumentError("Some channels do not have locations."))

    pal = mono == true ? :grays : :darktest
    
    # channel labels
    labels == [""] && (labels = repeat([""], size(s_pow, 1)))

    # get limits
    ylim = (floor(minimum(signal) * 1.1, digits=0), ceil(maximum(signal) * 1.1, digits=0))
    ylim = _tuple_max(ylim)

    # plot parameters
    plot_size = 1200
    marker_size = (150, 75)
    
    # get locations
    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = loc_x[channel]
    loc_y = loc_y[channel]
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)
    # get marker centers
    loc_x .*= ((plot_size / 2) - marker_size[1] / 2)
    loc_y .*= ((plot_size / 2) - marker_size[2] / 2)

    fig = Figure(; resolution=(plot_size, plot_size))
    fig_axis = Axis(fig[1, 1])
    fig_axis.aspect = AxisAspect(1)
    fig_axis.title = title
    GLMakie.xlims!(fig_axis, [-plot_size / 1.75, plot_size / 1.75])
    GLMakie.ylims!(fig_axis, [-plot_size / 1.75, plot_size / 1.75])
    hidedecorations!(fig_axis, grid=true, ticks=true)

    for idx in 1:size(signal, 1)
        p = Plots.plot(xlabel=xlabel,
                       ylabel=ylabel,
                       legend=false,
                       xticks=false,
                       yticks=false,
                       grid=false,
                       border=:none, 
                       xlims=(t[1], t[end]),
                       ylims=ylim,
                       title=labels[idx],
                       palette=pal,
                       size=marker_size,
                       titlefontsize=8,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6;
                       kwargs...)
        # plot 0 h-line
        p = Plots.hline!([0],
                         color=:grey,
                         linewidth=0.5,
                         labels="")

        # plot ERP
        p = Plots.plot!(t,
                        signal[idx, :],
                        t=:line,
                        color=:black,
                        linewidth=0.5,
                        alpha=0.75)

        # plot 0 v-line
        p = Plots.vline!([0],
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:black,
                         label=false)

        marker_img = tempname() * ".png"
        savefig(p, marker_img)
        marker = load(marker_img)
        GLMakie.scatter!(fig_axis, (loc_x[idx], loc_y[idx]), marker=marker, markersize=marker_size)
        rm(marker_img)
    end

    return fig
end

"""
    eeg_plot_erp(eeg; <keyword arguments>)

Plot ERP.

# Arguments

- `eeg::NeuroAnalyzer.EEG`: EEG object
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points 
- `xlabel::String="default"`: x-axis label, default is Time [ms]
- `ylabel::String="default"`: y-axis label, default is Amplitude [μV] 
- `title::String="default"`: plot title, default is ERP amplitude [channel: 1, epochs: 1:2, time window: -0.5 s:1.5 s]
- `cb::Bool=true`: plot color bar
- `cb_title::String="default"`: color bar title, default is Amplitude [μV] 
- `mono::Bool=false`: use color or grey palette
- `peaks::Bool=true`: draw peaks
- `labels::Bool=true`: draw labels legend (using EEG channel labels) for multi-channel `:butterfly` plot
- `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`), topographical plot of ERPs (`:topo`) or stacked epochs/channels (`:stack`)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_erp(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}, tm::Union{Int64, Vector{Int64}}=0, xlabel::String="default", ylabel::String="default", title::String="default", cb::Bool=true, cb_title::String="default", mono::Bool=false, peaks::Bool=true, labels::Bool=true, type::Symbol=:normal, kwargs...)

    _check_var(type, [:normal, :butterfly, :mean, :topo, :stack], "type")

    type in [:normal, :mean] && length(channel) > 1 && throw(ArgumentError("For :normal and :mean plot types, only one channel must be specified."))

    # check channels
    _check_channels(eeg, channel)

    # average all epochs
    epoch = 1:eeg_epoch_n(eeg)

    signal = eeg.eeg_signals[channel, :, epoch]

    # get time vector
    t = eeg.eeg_epoch_time
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])

    if tm != 0
        for tm_idx in eachindex(tm)
            tm[tm_idx] / 1000 < t[1] && throw(ArgumentError("tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."))
            tm[tm_idx] / 1000 > t[end] && throw(ArgumentError("tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."))
            tm[tm_idx] = vsearch(tm[tm_idx] / 1000, t)
        end
    end

    if type === :normal
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude channel $(_channel2channel_name(channel))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
        signal = mean(signal, dims=2)[:]
        p = plot_erp(t,
                     signal,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     title=title,
                     mono=mono;
                     kwargs...)
    elseif type === :butterfly
        if length(channel) > 1
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = mean(signal, dims=3)[:, :]
            if labels == true
                labels = eeg_labels(eeg)[channel]
            else
                labels = repeat([""], length(channel))
            end
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = signal'
            labels = [""]
        end
        p = plot_erp_butterfly(t,
                               signal,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               labels=labels,
                               mono=mono;
                               kwargs...)
    elseif type === :mean
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude [mean ± 95%CI] channel $(_channel2channel_name(channel))\n[averaged epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
        p = plot_erp_avg(t,
                         signal,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         mono=mono;
                         kwargs...)
    elseif type === :topo
        eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "", "", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
        peaks = false
        signal = mean(signal, dims=3)[:, :]
        ndims(signal) == 1 && (signal = reshape(signal, 1, length(signal)))
        labels = eeg_labels(eeg)[channel]
        typeof(labels) == String && (labels = [labels])
        p = plot_erp_topo(eeg.eeg_locs,
                          t,
                          signal,
                          channel=channel,
                          labels=labels,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          mono=mono;
                          kwargs...)
    elseif type === :stack
        peaks = false
        cb_title == "default" && (cb_title = "Amplitude [μV]")
        if length(channel) > 1
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = mean(signal, dims=3)[:, :]
            if labels == true
                labels = eeg_labels(eeg)[channel]
            else
                labels = repeat([""], length(channel))
            end
            ylabel = "Channel"
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = signal'
            labels = [""]
            ylabel = "Epoch"
        end
        p = plot_erp_stack(t,
                           signal,
                           xlabel=xlabel,
                           ylabel=ylabel,
                           title=title,
                           labels=labels,
                           cb=cb,
                           cb_title=cb_title,
                           mono=mono;
                           kwargs...)
    end

    # draw time markers
    if tm != 0
        for tm_idx in eachindex(tm)
            p = Plots.vline!([t[tm[tm_idx]]],
                             linewidth=1,
                             linecolor=:black,
                             label=false)
        end
    end

    # draw peaks
    if peaks == true
        if length(channel) == 1
            erp = eeg_erp(eeg).eeg_signals
            pp = eeg_erp_peaks(eeg)
            if mono == false
                Plots.scatter!((t[pp[channel, 1]], erp[channel, pp[channel, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[channel, 2]], erp[channel, pp[channel, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[channel, 1]], erp[channel, pp[channel, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[channel, 2]], erp[channel, pp[channel, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[channel, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(erp[channel, pp[channel, 1]], digits=2)) μV")
            _info("Negative peak time: $(round(t[pp[channel, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(erp[channel, pp[channel, 2]], digits=2)) μV")
        else
            erp = mean(eeg_erp(eeg).eeg_signals[channel, :], dims=1)[:]
            eeg_tmp = eeg_keep_channel(eeg, channel=1)
            eeg_tmp.eeg_signals = reshape(erp, 1, length(erp), 1)
            pp = eeg_erp_peaks(eeg_tmp)
            if mono == false
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(erp[pp[1, 1]], digits=2)) μV")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(erp[pp[1, 2]], digits=2)) μV")
        end
    end

    if type !== :topo
        Plots.plot(p)
    else
        GLMakie.show(p)
    end

    return p
end

"""
    plot_erp_stack(signal; <keyword arguments>)

Plot EPRs stacked by channels or by epochs.

# Arguments

- `t::AbstractVector`: x-axis values
- `signal::AbstractArray`
- `labels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `cb::Bool=true`: plot color bar
- `cb_title::String=""`: color bar title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_stack(t::AbstractVector, signal::AbstractArray; labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", cb::Bool=true, cb_title::String="", mono::Bool=false, kwargs...)

    ndims(signal) == 2 || throw(ArgumentError("signal must have 2 dimensions."))
    length(t) == size(signal, 2) || throw(ArgumentError("Number of signal columns ($(size(signal, 2))) must be equal to length of x-axis values ($(length(t)))."))

    pal = mono == true ? :grays : :darktest

    if labels == [""]
        yticks = round.(Int64, range(1, size(signal, 1), length=10))
    else
        yticks = (1:size(signal, 1), labels)
    end
    p = Plots.heatmap(t,
                      1:size(signal, 1),
                      signal,
                      size=(1200, 500),
                      margins=20Plots.px,
                      legend=false,
                      xticks=(_erpticks(t), string.(_erpticks(t) .* 1000)),
                      yticks=yticks,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      cb=cb,
                      cbtitle=cb_title,
                      title=title,
                      seriescolor=pal,
                      linewidth=0.5,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      xtickfontsize=8,
                      ytickfontsize=8;
                      kwargs...)

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    Plots.plot(p)

    return p
end

"""
    eeg_plot_erp(eeg, c; <keyword arguments>)

Plot ERP.

# Arguments

- `eeg::NeuroAnalyzer.EEG`: EEG object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points 
- `xlabel::String="default"`: x-axis label, default is Time [ms]
- `ylabel::String="default"`: y-axis label, default is Amplitude [μV] 
- `title::String="default"`: plot title, default is ERP amplitude [component: 1, epochs: 1:2, time window: -0.5 s:1.5 s]
- `cb::Bool=true`: plot color bar
- `cb_title::String="default"`: color bar title
- `mono::Bool=false`: use color or grey palette
- `peaks::Bool=true`: draw peaks
- `labels::Bool=true`: draw labels legend (using EEG component labels) for multi-channel `:butterfly` plot
- `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`), topographical plot of ERPs (`:topo`) or stacked epochs/channels (`:stack`)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_erp(eeg::NeuroAnalyzer.EEG, c::Union{Symbol, AbstractArray}; c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, tm::Union{Int64, Vector{Int64}}=0, xlabel::String="default", ylabel::String="default", title::String="default", cb::Bool=true, cb_title::String="default", mono::Bool=false, peaks::Bool=true, labels::Bool=true, type::Symbol=:normal, kwargs...)

    _check_var(type, [:normal, :butterfly, :mean, :topo, :stack], "type")

    # select component channels, default is all channels
    typeof(c) == Symbol && (c = _get_component(eeg, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    if labels == true
        labels = _gen_clabels(c)[c_idx]
    else
        labels = repeat([""], length(c_idx))
    end
    length(c_idx) == 1 && (labels = [labels])

    type in [:normal, :mean] && length(c_idx) > 1 && throw(ArgumentError("For :normal and :mean plot types, only one component channel must be specified."))

    # average all epochs
    epoch = 1:eeg_epoch_n(eeg)

    signal = c[c_idx, :, epoch]

    # get time vector
    t = eeg.eeg_epoch_time
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])

    if tm != 0
        for tm_idx in eachindex(tm)
            tm[tm_idx] / 1000 < t[1] && throw(ArgumentError("tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."))
            tm[tm_idx] / 1000 > t[end] && throw(ArgumentError("tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."))
            tm[tm_idx] = vsearch(tm[tm_idx] / 1000, t)
        end
    end

    if type === :normal
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude component $(_channel2channel_name(c_idx))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
        signal = mean(signal, dims=2)[:]
        p = plot_erp(t,
                     signal,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     title=title,
                     mono=mono;
                     kwargs...)
    elseif type === :butterfly
        if length(c_idx) > 1
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = mean(signal, dims=3)[:, :]
            if labels == true
                labels = _gen_clabels(c)[c_idx]
            else
                labels = repeat([""], length(c_idx))
            end
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = signal'
            labels = [""]
        end
        p = plot_erp_butterfly(t,
                               signal,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               labels=labels,
                               mono=mono;
                               kwargs...)
    elseif type === :mean
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude [mean ± 95%CI] component $(_channel2channel_name(c_idx))\n[averaged epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
        p = plot_erp_avg(t,
                         signal,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         mono=mono;
                         kwargs...)
    elseif type === :topo
        eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "", "", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
        peaks = false
        signal = mean(signal, dims=3)[:, :]
        ndims(signal) == 1 && (signal = reshape(signal, 1, length(signal)))
        typeof(labels) == String && (labels = [labels])
        p = plot_erp_topo(eeg.eeg_locs,
                          t,
                          signal,
                          channel=c_idx,
                          labels=labels,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          mono=mono;
                          kwargs...)
    elseif type === :stack
        cb_title == "default" && (cb_title = "Amplitude [μV]")
        peaks = false
        if length(c_idx) > 1
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = mean(signal, dims=3)[:, :]
            ylabel = "Component channel"
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = signal'
            labels = [""]
            ylabel = "Epoch"
        end
        p = plot_erp_stack(t,
                           signal,
                           xlabel=xlabel,
                           ylabel=ylabel,
                           title=title,
                           labels=labels,
                           cb=cb,
                           cb_title=cb_title,
                           mono=mono;
                           kwargs...)
    end

    # draw time markers
    if tm != 0
        for tm_idx in eachindex(tm)
            p = Plots.vline!([t[tm[tm_idx]]],
                             linewidth=1,
                             linecolor=:black,
                             label=false)
        end
    end

    # draw peaks
    if peaks == true
        signal = c[c_idx, :, epoch]
        if length(c_idx) == 1
            erp = mean(signal, dims=2)[:]
            eeg_tmp = eeg_keep_channel(eeg, channel=1)
            eeg_tmp.eeg_signals = reshape(erp, 1, length(erp), 1)
            pp = eeg_erp_peaks(eeg_tmp)
            if mono == false
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(erp[pp[1, 1]], digits=2)) μV")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(erp[pp[1, 2]], digits=2)) μV")
        else         
            erp = mean(mean(signal, dims=3), dims=1)[:]
            eeg_tmp = eeg_keep_channel(eeg, channel=1)
            eeg_tmp.eeg_signals = reshape(erp, 1, length(erp), 1)
            pp = eeg_erp_peaks(eeg_tmp)
            if mono == false
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(erp[pp[1, 1]], digits=2)) μV")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(erp[pp[1, 2]], digits=2)) μV")
        end
    end

    if type !== :topo
        Plots.plot(p)
    else
        GLMakie.show(p)
    end

    return p
end

"""
    plot_dipole3d(eeg, c; <keyword arguments>)

Plot dipole in 3D.

# Arguments

- `d::NeuroAnalyzer.DIPOLE)`

# Returns
- `p::GLMakie.Figure`
"""
function plot_dipole3d(d::NeuroAnalyzer.DIPOLE)

    brain_top = Point3f[[-1.5,-1.5,-0.5],[1.5,-1.5,-0.5], [1.5,1.5,-0.5], [-1.5,1.5,-0.5]]
    brain_top_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_top_fs = GLTriangleFace[(1, 2, 3), (1, 3, 4)]
    brain_top_mesh = GeometryBasics.Mesh(meta(brain_top, uv = brain_top_uvs, normals = normals(brain_top, brain_top_fs)), brain_top_fs)

    brain_side = Point3f[[-1.5,-1.5,-0.5],[-1.5,1.5,-0.5], [-1.5,1.5,1.0], [-1.5,-1.5,1.0]]
    brain_side_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_side_fs = GLTriangleFace[(1, 2, 3), (1, 3, 4)]
    brain_side_mesh = GeometryBasics.Mesh(meta(brain_side, uv = brain_side_uvs, normals = normals(brain_side, brain_side_fs)), brain_side_fs)

    brain_front = Point3f[[-1.5,1.5,-0.5],[-1.5,1.5,1.0], [1.5,1.5,1.0], [1.5,1.5,-0.5]]
    brain_front_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_front_fs = GLTriangleFace[(1, 2, 3), (1, 3, 4)]
    brain_front_mesh = GeometryBasics.Mesh(meta(brain_front, uv = brain_front_uvs, normals = normals(brain_front, brain_front_fs)), brain_front_fs)

    brain_top_texture = FileIO.load("images/brain_top.png")
    brain_side_texture = FileIO.load("images/brain_side.png")
    brain_front_texture = FileIO.load("images/brain_front.png")

    x = d.loc[1]
    y = d.loc[2]
    z = d.loc[3]

    p = Figure()
    ax = Axis3(p[1, 1])
    mesh!(ax, brain_side_mesh, color=brain_side_texture)
    mesh!(ax, brain_top_mesh, color=brain_top_texture)
    mesh!(ax, brain_front_mesh, color=brain_front_texture)

    # draw dipole
    GLMakie.scatter!(ax, x, y, z, markersize=20, color=:red)

    # project at top-plane
    GLMakie.lines!(ax, [x, x], [y, y], [z, -0.5], linestyle=:dash, color=:blue)
    # project at side-axis
    GLMakie.lines!(ax, [x, -1.5], [y, y], [z, z], linestyle=:dash, color=:blue)
    # project at front-axis
    GLMakie.lines!(ax, [x, x], [y, 1.5], [z, z], linestyle=:dash, color=:blue)
    
    GLMakie.show(p)
    return p
end