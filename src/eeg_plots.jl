"""
    eeg_plot_save(p; file_name::String)

Saves plot as file (PDF/PNG/TIFF). File format is determined using `file_name` extension.

# Arguments

- `p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}`
- `file_name::String`
"""
function eeg_plot_save(p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}; file_name::String)

    ext = splitext(file_name)[2]
    _check_var(ext, [".png", ".pdf", ".jpg", ".tiff"], "File format")
    (isfile(file_name) && verbose == true) && @info "File $file_name will be overwritten."
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
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=(-1, channel_n),
                   title=title,
                   palette=pal,
                   size=(1200, 800),
                   left_margin=20Plots.px,
                   bottom_margin=20Plots.px,
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
function plot_signal_avg(t::Union{AbstractVector, AbstractRange}, signal::AbstractArray; xlabel::String="default", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="μV", norm::Bool=false, kwargs...)

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
                   size=(1200, 800),
                   left_margin=20Plots.px,
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
                   size=(1200, 800),
                   left_margin=20Plots.px,
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
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, AbstractRange}=0, channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg)), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, markers::Bool=true, scale::Bool=true, units::String="μV", type::Symbol=:normal, norm::Bool=false, kwargs...)

    _check_var(type, [:normal, :butterfly, :mean], "type")
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
        signal = eeg_epochs(eeg, epoch_n=1).eeg_signals[channel, segment[1]:segment[2], 1]
    end
    t = _get_t(segment[1], segment[2], eeg_sr(eeg))

    t_1, t_s1, t_2, t_s2 = _convert_t(t[1], t[end])
    epoch = _t2epoch(eeg, segment[1], segment[2])

    if type === :normal
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
    elseif type === :butterfly
        size(signal, 1) == 1 && throw(ArgumentError("For type=:butterfly plot the signal must contain ≥ 2 channels."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Channels $(_channel2channel_name(channel)) amplitude\n[epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
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
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Averaged channels $(_channel2channel_name(channel)) amplitude [mean ± 95%CI]\n [epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
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
        for idx in 1:length(markers_desc)
            p = Plots.plot!(annotation=(markers_pos[idx], -0.95, Plots.text("$(markers_desc[idx])", pointsize=4, halign=:left, valign=:top, rotation=90)), label=false)
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

    _check_var(type, [:normal, :butterfly, :mean], "type")
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

    # get time vector
    if segment[2] <= eeg_epoch_len(eeg)
        signal = c[c_idx, segment[1]:segment[2], 1]
    else
        signal = _make_epochs(c, epoch_n=1)[c_idx, segment[1]:segment[2], 1]
    end
    t = _get_t(segment[1], segment[2], eeg_sr(eeg))

    t_1, t_s1, t_2, t_s2 = _convert_t(t[1], t[end])
    epoch = _t2epoch(eeg, segment[1], segment[2])

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
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Components $(_channel2channel_name(c_idx)) amplitude\n[epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
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
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Averaged components $(_channel2channel_name(c_idx)) amplitude [mean ± 95%CI]\n[epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
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
        for idx in 1:length(markers_desc)
            p = Plots.plot!(annotation=(markers_pos[idx], -0.95, Plots.text("$(markers_desc[idx])", pointsize=4, halign=:left, valign=:top, rotation=90)), label=false)
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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
                   size=(1200, 800),
                   left_margin=20Plots.px,
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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
                   size=(1200, 800),
                   left_margin=20Plots.px,
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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
                   size=(1200, 800),
                   left_margin=20Plots.px,
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

- `p::Plots.Plot{Plots.GRBackend}`
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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
                       left_margin=20Plots.px,
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
                       left_margin=20Plots.px,
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
    plot_psd_topo(s_frq, s_pow; <keyword arguments>)

Plot topographical map `eeg` PSD. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `s_frq::Vector{Float64}`: frequencies
- `s_pow::Array{Float64, 3}`: powers
- `Union{Vector{Int64}, AbstractRange}`: which channels to plot
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
function plot_psd_topo(locs::DataFrame, s_frq::Vector{Float64}, s_pow::Array{Float64, 2}; channels=Union{Vector{Int64}, AbstractRange}, labels::Vector{String}=[""], norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, ax::Symbol=:linlin, kwargs...)

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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
            verbose == true && @info "Lower frequency bound truncated to 0.1 Hz"
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
    loc_x, loc_y = _locnorm(loc_x, loc_y)
    # get marker centers
    loc_x .*= ((plot_size / 2) - marker_size[1] / 2)
    loc_y .*= ((plot_size / 2) - marker_size[2] / 2)
    loc_x = loc_x[channels]
    loc_y = loc_y[channels]

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
- `method::Symbol=:welch`: method of calculating PSD: Welch's periodogram, (`:welch`), multi-tapered periodogram (`:mt`), Morlet wavelet convolution (`:mw`)
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

- `p::Plots.Plot{Plots.GRBackend} | GLMakie.Figure`
"""
function eeg_plot_psd(eeg::NeuroAnalyzer.EEG; epoch::Int64, channel::Union{Int64, Vector{Int64}, AbstractRange}, norm::Bool=true, method::Symbol=:welch, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, kwargs...)

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
    t_1, t_s1, t_2, t_s2 = _convert_t(eeg.eeg_epochs_time[1], eeg.eeg_epochs_time[end])

    if ref === :abs
        if method === :welch
            s_pow, s_frq = s_psd(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mt
            s_pow, s_frq = s_psd(signal, fs=fs, norm=norm, mt=true)
            title == "default" && (title = "Absolute PSD (multi-tapered) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mw
            s_pow, s_frq = s_wspectrum(signal, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc)
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
            title == "default" && (title = "Absolute PSD (Welch's periodogram) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mt
            s_pow, s_frq = s_rel_psd(signal, fs=fs, norm=norm, mt=true, f=f)
            title == "default" && (title = "Absolute PSD (multi-tapered) relative to $ref power [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[channel: $channel, epoch: $epoch, time window: $t_s1:$t_s2]")
        end
    end

    # set labels
    if type !== :w3d && type !== :s3d && type !== :topo
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if norm == true
            ylabel == "default" && (ylabel = "Power [dB]")
        else
            ylabel == "default" && (ylabel = "Power [μV^2/Hz]")
        end
    end

    if type === :normal
        ndims(s_pow) > 1 && throw(ArgumentError("For type=:normal the signal must contain 1 channel."))
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
                          channels=channel,
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
- `method::Symbol=:welch`: method of calculating PSD: Welch's periodogram, (`:welch`), multi-tapered periodogram (`:mt`), Morlet wavelet convolution (`:mw`)
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
function eeg_plot_psd(eeg::NeuroAnalyzer.EEG, c::Union{Symbol, AbstractArray}; epoch::Int64, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, norm::Bool=true, method::Symbol=:welch, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, kwargs...)

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
    t_1, t_s1, t_2, t_s2 = _convert_t(eeg.eeg_epochs_time[1], eeg.eeg_epochs_time[end])

    if ref === :abs
        if method === :welch
            s_pow, s_frq = s_psd(signal, fs=fs, norm=norm, mt=false)
            title == "default" && (title = "Absolute PSD (Welch's periodogram) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mt
            s_pow, s_frq = s_psd(signal, fs=fs, norm=norm, mt=true)
            title == "default" && (title = "Absolute PSD (multi-tapered) [frequency limit: $(frq_lim[1])-$(frq_lim[2]) Hz]\n[component: $(_channel2channel_name(c_idx)), epoch: $epoch, time window: $t_s1:$t_s2]")
        elseif method === :mw
            s_pow, s_frq = s_wspectrum(signal, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc)
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

    palette = mono == true ? :grays : :darktest
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
                      left_margin=20Plots.px,
                      seriescolor=palette,
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

    palette = mono == true ? :grays : :darktest
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
                      left_margin=20Plots.px,
                      seriescolor=:darktest,
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
- `method::Symbol=:standard`: method of calculating spectrogram: standard (`:standard`), short-time Fourier transform (`:stft`), multi-tapered periodogram (`:mt`), Morlet wavelet convolution (`:mw`)
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
function eeg_plot_spectrogram(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, AbstractRange}=0, channel::Union{Int64, Vector{Int64}, AbstractRange}, norm::Bool=true, method::Symbol=:standard, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, kwargs...)

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
    t_1, t_s1, t_2, t_s2 = _convert_t(eeg.eeg_epochs_time[1], eeg.eeg_epochs_time[end])

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
            for idx in 1:length(markers_desc)
                p = Plots.plot!(annotation=(markers_pos[idx], -0.95, Plots.text("$(markers_desc[idx])", pointsize=4, halign=:left, valign=:top, rotation=90)), label=false)
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
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=true)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            @info "Method :stft is not available for multi-channel spectrogram, using standard periodogram."
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mw
            s_p, s_f = s_wspectrum(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_f = linspace(0, frq_lim[2], size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
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
- `method::Symbol=:standard`: method of calculating spectrogram: standard (`:standard`), short-time Fourier transform (`:stft`), multi-tapered periodogram (`:mt`), Morlet wavelet convolution (`:mw`)
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
function eeg_plot_spectrogram(eeg::NeuroAnalyzer.EEG, c::Union{Symbol, AbstractArray}; epoch::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, AbstractRange}, norm::Bool=true, method::Symbol=:standard, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, kwargs...)

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
    t_1, t_s1, t_2, t_s2 = _convert_t(eeg.eeg_epochs_time[1], eeg.eeg_epochs_time[end])

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
            for idx in 1:length(markers_desc)
                p = Plots.plot!(annotation=(markers_pos[idx], -0.95, Plots.text("$(markers_desc[idx])", pointsize=4, halign=:left, valign=:top, rotation=90)), label=false)
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
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=true)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(multi-tapered periodogram)")
        elseif method === :stft
            @info "Method :stft is not available for multi-channel spectrogram, using standard periodogram."
            s_p, s_f = s_psd(signal, fs=fs, norm=false, mt=false)
            f1 = vsearch(frq_lim[1], s_f)
            f2 = vsearch(frq_lim[2], s_f)
            s_f = s_f[f1:f2]
            s_p = s_p[:, f1:f2]
            title = replace(title, "method" => "(standard periodogram)")
        elseif method === :mw
            s_p, s_f = s_wspectrum(signal, fs=fs, frq_lim=frq_lim, frq_n=length(frq_lim[1]:frq_lim[2]), ncyc=ncyc, norm=false)
            s_f = linspace(0, frq_lim[2], size(s_p, 2))
            title = replace(title, "method" => "(Morlet-wavelet transform)")
        end

        norm == true && (s_p = pow2db.(s_p))
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
- `plot_size::Int64=400`: plot dimensions in pixels (size × size)

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_electrodes(locs::DataFrame; channel::Union{Int64, Vector{Int64}, AbstractRange}, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head_labels::Bool=true, mono::Bool=false, head_details::Bool=true, plot_size::Int64=400)

    pal = mono == true ? :grays : :darktest

    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    loc_x, loc_y = _locnorm(loc_x, loc_y)

    if plot_size > 300
        marker_size = plot_size ÷ 75
        font_size = plot_size ÷ 75
    else
        marker_size = plot_size ÷ 50
        font_size = plot_size ÷ 50
        labels = false
    end

    loc_x .*= 0.9
    loc_y .*= 0.9
    p = Plots.plot(grid=true,
                   framestyle=:none,
                   palette=pal,
                   size=(plot_size, plot_size),
                   border=:none,
                   aspect_ratio=1,
                   right_margin=-30 * Plots.px,
                   bottom_margin=-30 * Plots.px,
                   top_margin=-30 * Plots.px,
                   left_margin=-50 * Plots.px,
                   xlim=(-1.22, 1.23),
                   ylim=(-1.1, 1.2))

    hd = _draw_head(p, head_labels=head_labels, head_details=head_details)
    p = Plots.plot!(hd)

    for idx in 1:length(locs[!, :labels])
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
        for idx in 1:length(locs[!, :labels])
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

3D interactive preview of electrode locations. It uses spherical :loc_x, :loc_y and :loc_z locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
- `labels::Bool=true`: plot electrode labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or grey palette
 `plot_size::Int64=800`: plot dimensions in pixels (plot_size×plot_size)

# Returns

- `fig::GLMakie.Figure`
"""
function plot_electrodes3d(locs::DataFrame; channel::Union{Int64, Vector{Int64}, AbstractRange}, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head_labels::Bool=true, mono::Bool=false, plot_size::Int64=800)

    # selected != 0 && length(intersect(channel, selected)) < length(selected) && throw(ArgumentError("channel must include selected."))
    # channel = setdiff(channel, selected)

    pal = mono == true ? :grays : :darktest

    loc_x = locs[!, :loc_x]
    loc_y = locs[!, :loc_y]
    loc_z = locs[!, :loc_z]

    loc_x, loc_y = _locnorm(loc_x, loc_y)
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
        for idx in 1:length(locs[!, :labels])
            if idx in channel
                GLMakie.text!(ax, locs[!, :labels][idx], position=(loc_x[idx], loc_y[idx], loc_z[idx]), textsize=font_size)
            end
            if idx in selected
                GLMakie.text!(ax, locs[!, :labels][idx], position=(loc_x[idx], loc_y[idx], loc_z[idx]), textsize=font_size)
            end
        end
    end

    if head_labels == true
        GLMakie.text!(ax, "Nz", position=(0, 1.025, 0), textsize = font_size)
        GLMakie.text!(ax, "Iz", position=(0, -1.025, 0), textsize = font_size)
        GLMakie.text!(ax, "LPA", position=(-1.025, 0, 0), textsize = font_size)
        GLMakie.text!(ax, "RPA", position=(1.025, 0, 0), textsize = font_size)
        GLMakie.text!(ax, "top", position=(0, 0, 1.025), textsize = font_size)
    end
    fig

    return fig
end

"""
    eeg_plot_electrodes(eeg; <keyword arguments>)

Preview of electrode locations.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: which channel should be highlighted
- `labels::Bool=true`: plot electrode labels
- `head::Bool`=true: plot head
- `head_labels::Bool=false`: plot head labels
- `plot_size::Int64=400`: plot dimensions in pixels (plot_size×plot_size)
- `head_details::Bool=true`: draw nose and ears
- `mono::Bool=false`: use color or grey palette
- `threed::Bool=false`: 3-dimensional plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_electrodes(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type])), selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head::Bool=true, head_labels::Bool=false, plot_size::Int64=400, head_details::Bool=true, mono::Bool=false, threed::Bool=false, kwargs...)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    # select channels, default is all channels
    _check_channels(eeg, channel, Symbol(eeg.eeg_header[:signal_type]))
    selected != 0 && _check_channels(eeg, selected)

    if threed == false
        p = plot_electrodes(eeg.eeg_locs, channel=channel, selected=selected, labels=labels, head_labels=head_labels, mono=mono, head_details=head_details, plot_size=plot_size)
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

Plot histogram of `signal`.

# Arguments

- `signal::AbstractVector`
- `type::Symbol`: type of histogram: regular (`:hist`) or kernel density (`:kd`)
- `label::String=""`: channel label
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_histogram(signal::AbstractVector; type::Symbol=:hist, label::String="", xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

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
                   grid=false,
                   linecolor=:black,
                   fillcolor=:grey,
                   fillalpha=0.5,
                   linewidth=1,
                   size=(600, 400),
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
    plot_filter_response(<keyword arguments>)

Plot filter response.

# Arguments

- `fs::Int64`: sampling rate
- `fprototype::Symbol`: filter class: `:fir`, `:butterworth`, `:chebyshev1`, `:chebyshev2`, `:elliptic`
- `ftype::Symbol`: filter type: low-pass (`:lp`), high-pass (`:hp`), band-pass (`:bp`), band-stop (`:bs`)
- `cutoff::Union{Real, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64`: filter order
- `rp::Real`: dB ripple in the passband
- `rs::Real`: dB attenuation in the stopband
- `window::window::Union{Vector{Float64}, Nothing}`: window, required for FIR filter
- `mono::Bool=false`: use color or grey palette
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_filter_response(; fs::Int64, fprototype::Symbol, ftype::Symbol, cutoff::Union{Real, Tuple}, order::Int64=-1, rp::Real=-1, rs::Real=-1, window::Union{Vector{Float64}, Nothing}=nothing, mono::Bool=false, frq_lim::Tuple{Real, Real}=(0, 0), kwargs...)

    _check_var(fprototype, [:fir, :butterworth, :chebyshev1, :chebyshev2, :elliptic], "fprototype")
    _check_var(ftype, [:lp, :hp, :bp, :bs], "ftype")

    fprototype !== :fir && order < 1 && throw(ArgumentError("order must be > 0."))

    pal = mono == true ? :grays : :darktest

    frq_lim == (0, 0) && (frq_lim = (0, fs / 2))
    frq_lim = tuple_order(frq_lim)

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
    if fprototype === :fir
        if window === nothing
            verbose == true && @info "Using default window for :fir filter: hanning($(3 * floor(Int64, fs / cutoff[1])))."
            window = hanning(3 * floor(Int64, fs / cutoff[1]))
        end
        if ftype === :hp || ftype === :bp || ftype === :bs
            mod(length(window), 2) == 0 && (window = vcat(window[1:((length(window) ÷ 2) - 1)], window[((length(window) ÷ 2) + 1):end]))
        end
        prototype = FIRWindow(window)
    end
    if fprototype === :chebyshev1
        (rs < 0 || rs > eeg_sr(eeg) / 2) && throw(ArgumentError("For :chebyshev1 filter rs must be ≥ 0 and ≤ $(eeg_sr(eeg) / 2)."))
        prototype = Chebyshev1(order, rs)
    end
    if fprototype === :chebyshev2
        (rp < 0 || rp > eeg_sr(eeg) / 2) && throw(ArgumentError("For :chebyshev2 filter rp must be ≥ 0 and ≤ $(eeg_sr(eeg) / 2)."))
        prototype = Chebyshev2(order, rp)
    end
    if fprototype === :elliptic
        (rs < 0 || rs > eeg_sr(eeg) / 2) && throw(ArgumentError("For :elliptic filter rs must be ≥ 0 and ≤ $(eeg_sr(eeg) / 2)."))
        (rp < 0 || rp > eeg_sr(eeg) / 2) && throw(ArgumentError("For :elliptic filter rp must be ≥ 0 and ≤ $(eeg_sr(eeg) / 2)."))
        prototype = Elliptic(order, rp, rs)
    end

    ffilter = digitalfilter(responsetype, prototype)

    if fprototype !== :fir
        H, w = freqresp(ffilter)
        # convert to dB
        H = 20 * log10.(abs.(H))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        p1 = Plots.plot(w,
                        H,
                        title="Filter: $(titlecase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order\nFrequency response",
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylims=(-100, 0),
                        ylabel="Magnitude\n[dB]",
                        xlabel="Frequency [Hz]",
                        label="",
                        titlefontsize=5,
                        xlabelfontsize=4,
                        ylabelfontsize=4,
                        xtickfontsize=3,
                        ytickfontsize=3,
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

        phi, w = phaseresp(ffilter)
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
                        titlefontsize=5,
                        xlabelfontsize=4,
                        ylabelfontsize=4,
                        xtickfontsize=3,
                        ytickfontsize=3,
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

        tau, w = grpdelay(ffilter)
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
                        titlefontsize=5,
                        xlabelfontsize=4,
                        ylabelfontsize=4,
                        xtickfontsize=3,
                        ytickfontsize=3,
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

        p = Plots.plot(p1, p2, p3, size=(1200, 800), left_margin=20*Plots.px, layout=(3, 1), palette=pal; kwargs...)
    else
        w = range(0, stop=pi, length=1024)
        H = _fir_response(ffilter, w)
        # convert to dB
        H = 20 * log10.(abs.(H))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        p1 = Plots.plot(w,
                        H,
                        title="Filter: $(uppercase(String(fprototype))), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz\nFrequency response",
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylims=(-100, 0),
                        ylabel="Magnitude\n[dB]",
                        xlabel="Frequency [Hz]",
                        label="",
                        titlefontsize=5,
                        xlabelfontsize=4,
                        ylabelfontsize=4,
                        xtickfontsize=3,
                        ytickfontsize=3,
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
        phi = _fir_response(ffilter, w)
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
                        titlefontsize=5,
                        xlabelfontsize=4,
                        ylabelfontsize=4,
                        xtickfontsize=3,
                        ytickfontsize=3,
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

        p = Plots.plot(p1, p2, size=(1200, 800), left_margin=20*Plots.px, layout=(2, 1), palette=pal; kwargs...)
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
    loc_x, loc_y = _locnorm(loc_x, loc_y)

    for idx in 1:length(locs[!, :labels])
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

    for idx in 1:length(locs[!, :labels])
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
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
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
function eeg_plot_weights(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type])), weights::Vector{<:Real}, labels::Bool=true, head_labels::Bool=false, mono::Bool=false, head_details::Bool=true, plot_size::Int64=800, title::String="", kwargs...)

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
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
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
function eeg_plot_connections(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type])), connections::Matrix{<:Real}, threshold::Real, threshold_type::Symbol=:g, weights::Bool=true, labels::Bool=true, head_labels::Bool=false, mono::Bool=false, head_details::Bool=true, plot_size::Int64=800, title::String="", kwargs...)

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
    loc_x, loc_y = _locnorm(loc_x, loc_y)

    for idx in 1:length(locs[!, :labels])
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
        for idx in 1:length(locs[!, :labels])
            if idx in channel
                Plots.plot!(annotation=(loc_x[idx], loc_y[idx] + 0.05, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
        end
    end

    hd = _draw_head(p, head_labels=head_labels, head_details=head_details)
    p = Plots.plot!(hd)

    m_tmp = s_normalize_max(connections)

#    loc_x = loc_x[channel]
#    loc_y = loc_y[channel]

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
    eeg_plot(eeg, c; <keyword arguments>)

Plot channel/epoch data.

# Arguments

- `eeg::NeuroAnalyzer.EEG`: EEG object
- `c::Union{Vector{<:Real}, Matrix{<:Real}, Symbol, Dict}`: component to plot
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: epoch to display
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `plot_by::Symbol`: c values refer to: :labels, :channels or :epochs
- `type::Symbol`: plot type: histogram (`:hist`), kernel density (`:kd`), bar plot (`:bar`), box plot (`:box`), violin plot (`:violin`), paired (`:paired`) or polar (`:polar`); for `:box` and `:violin` `c` must contain ≥ 2 values per channel
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`

# Notes

Labeled matrix is a dictionary of labels and vectors associated with these labels. This way plotting of vectors non-equal length is possible. Labeled matrix is created using `_labeled_matrix2dict(l::Vector{String}, v::Vector{Vector{<:Real}})` and converted back to keys and values using `_dict2labeled_matrix(d::Dict)`.

For `:polar plot` if `c` is a vector, than it contains phases in radians. If `c` is a two column matrix, than first column contains phases in radians and second column contains lengths.
"""
function eeg_plot_stats(eeg::NeuroAnalyzer.EEG, c::Union{Vector{<:Real}, Matrix{<:Real}, Symbol, Dict}; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", plot_by::Symbol, type::Symbol, mono::Bool=false, kwargs...)

    _check_var(plot_by, [:channels, :epochs, :labels], "plot_by")
    _check_var(type, [:hist, :kd, :line, :bar, :box, :violin, :dots, :paired, :polar], "type")

    typeof(c) == Symbol && (c = _get_component(eeg, c).c)

    if plot_by === :channels
        channel == 0 && throw(ArgumentError("channel must be specified for plot by :channels."))
        epoch == 0 && (epoch = _select_epochs(eeg, epoch))
        _check_channels(eeg, channel)
        _check_epochs(eeg, epoch)
        labels = eeg_labels(eeg)[channel]
        if ndims(c) == 2
            c = c[channel, epoch]
        else
            length(c) < length(channel) && throw(ArgumentError("Length of c ($(length(c))) does not match the number of channels ($(length(channel)))"))
            c = c[channel]
        end
    elseif plot_by === :epochs
        epoch == 0 && throw(ArgumentError("epoch must be specified for plot by :epochs."))
        channel == 0 && (channel = _select_channels(eeg, channel))
        _check_epochs(eeg, epoch)
        _check_channels(eeg, channel)
        labels = "e" .* string.(collect(epoch))
        if ndims(c) == 2
            c = c[channel, epoch]
        else
            length(c) < length(epoch) && throw(ArgumentError("Length of c ($(length(c))) does not match the number of epochs ($(length(epoch)))"))
            c = c[epoch]
        end
    elseif plot_by === :labels
        if typeof(c) <: Dict
            l_tmp, c = _dict2labeled_matrix(c)
            l_tmp = reverse!(l_tmp)
            c = reverse!(c)
            labels == [""] && (labels = l_tmp)
            length(c) == length(labels) || throw(ArgumentError("Number of rows of c ($(length(c))) does not match number of labels ($(length(labels)))"))
        else
            if ndims(c) == 1
                length(c) == length(labels) || throw(ArgumentError("Length c ($(length(c))) does not match the number of labels ($(length(epoch)))"))
            else
                size(c, 1) == length(labels) || throw(ArgumentError("Number of rows of c ($(size(c, 1))) does not match number of labels ($(length(labels)))"))
            end
        end
    end

    pal = mono == true ? :grays : :darktest

    if type in [:bar, :hist, :kd]
        if plot_by === :epochs
            length(channel) > 1 && throw(ArgumentError("For :bar, :hist and :kd plots and plot by :epochs only one channel may be specified."))
            length(epoch) == 1 && throw(ArgumentError("More than 1 epoch must be specified."))
        elseif plot_by === :channels
            plot_by === :channels && length(epoch) > 1 && throw(ArgumentError("For :bar, :hist and :kd plots and plot by :channels only one epoch may be specified."))
            length(channel) == 1 && throw(ArgumentError("More than 1 channel must be specified."))
        end
    end

    if type === :hist
        p = plot_histogram(c,
                           type=:hist,
                            mono=mono,
                            xlabel=xlabel,
                            title=title;
                            kwargs)
    elseif type === :kd
        p = plot_histogram(c,
                           type=:kd,
                           mono=mono,
                           xlabel=xlabel,
                           title=title;
                           kwargs)
    elseif type === :line
        if ndims(c) == 1
            color = mono == true ? :lightgrey : :lightblue
            p = Plots.plot(c,
                           seriestype=:line,
                           size=(1200, 800),
                           left_margin=20Plots.px,
                           legend=false,
                           xticks=(1:length(labels), labels),
                           xlabel=xlabel,
                           ylabel=ylabel,
                           color=color,
                           title=title,
                           palette=pal,
                           titlefontsize=8,
                           xlabelfontsize=8,
                           ylabelfontsize=8,
                           xtickfontsize=8,
                           ytickfontsize=8,
                           kwargs=kwargs)
        else
            if plot_by === :channels
                p = Plots.plot(c[:, 1],
                               seriestype=:line,
                               size=(1200, 800),
                               left_margin=20Plots.px,
                               label="e 1",
                               legend=true,
                               xticks=(1:length(labels), labels),
                               xlabel=xlabel,
                               ylabel=ylabel,
                               color=1,
                               title=title,
                               palette=pal,
                               titlefontsize=8,
                               xlabelfontsize=8,
                               ylabelfontsize=8,
                               xtickfontsize=8,
                               ytickfontsize=8,
                               kwargs=kwargs)
                for idx in 2:size(c, 2)
                    p = Plots.plot!(c[:, idx],
                                    seriestype=:line,
                                    label="e $idx",
                                    color=idx)
                end
            elseif plot_by === :epochs
                p = Plots.plot(c[1, :],
                               seriestype=:line,
                               size=(1200, 800),
                               left_margin=20Plots.px,
                               label="ch 1",
                               legend=true,
                               xticks=(1:length(labels), labels),
                               xlabel=xlabel,
                               ylabel=ylabel,
                               color=1,
                               title=title,
                               palette=pal,
                               titlefontsize=8,
                               xlabelfontsize=8,
                               ylabelfontsize=8,
                               xtickfontsize=8,
                               ytickfontsize=8,
                               kwargs=kwargs)
                for idx in 2:size(c, 1)
                    p = Plots.plot!(c[idx, :],
                                    seriestype=:line,
                                    label="ch $idx",
                                    color=idx)
                end
            else
                p = Plots.plot(c,
                               seriestype=:line,
                               size=(1200, 800),
                               left_margin=20Plots.px,
                               legend=false,
                               xticks=(1:length(labels), labels),
                               xlabel=xlabel,
                               ylabel=ylabel,
                               color=1,
                               title=title,
                               palette=pal,
                               titlefontsize=8,
                               xlabelfontsize=8,
                               ylabelfontsize=8,
                               xtickfontsize=8,
                               ytickfontsize=8,
                               kwargs=kwargs)
            end
        end
    elseif type === :bar
        color = mono == true ? :lightgrey : :lightblue
        p = Plots.plot(c,
                       seriestype=:bar,
                       size=(1200, 800),
                       left_margin=20Plots.px,
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
                       ytickfontsize=8,
                       kwargs=kwargs)
    elseif type in [:box, :violin]
        color = mono == true ? :lightgrey : :auto
        if plot_by !== :labels
            p = Plots.plot(c',
                           seriestype=type,
                           size=(1200, 800),
                           left_margin=20Plots.px,
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
                           kwargs...)
        else
            p = Plots.plot(c,
                           seriestype=type,
                           size=(1200, 800),
                           left_margin=20Plots.px,
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
                           kwargs...)
        end
    elseif type === :dots
        p = Plots.plot(size=(1200, 800),
                       left_margin=20Plots.px,
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
        for idx1 in 1:length(labels)
            for idx2 in 1:length(c[idx1])
                if mono == false
                    p = Plots.scatter!((idx1, c[idx1][idx2]),
                                       color=idx1)
                else
                    p = Plots.scatter!((idx1, c[idx1][idx2]),
                                       color=:black)
                end
            end
        end
    elseif type === :paired
        ll = Vector{Int64}()
        for idx in 1:length(labels)
            push!(ll, length(c[idx]))
        end
        length(unique(ll)) == 1 || throw(ArgumentError("For :paired plot each label must have the same number of values."))
        p = Plots.plot(size=(1200, 800),
                       left_margin=20Plots.px,
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
        for idx1 in 1:length(labels)
            for idx2 in 1:length(c[idx1])
                if mono == false
                    p = Plots.scatter!((idx1, c[idx1][idx2]),
                                       color=idx2)
                else
                    p = Plots.scatter!((idx1, c[idx1][idx2]),
                                       color=:black)
                end
            end
        end
        for idx1 in 1:length(c[1])
            c_tmp = zeros(length(labels))
            for idx2 in 1:length(labels)
                c_tmp[idx2] = c[idx2][idx1]
            end
            p = Plots.plot!(c_tmp,
                            color=:black)
        end
    elseif type === :polar
        if ndims(c) == 1
            p = Plots.plot([0, c[1]], [0, 1],
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
                           linewidth=0.5,
                           titlefontsize=8,
                           xlabelfontsize=8,
                           ylabelfontsize=8,
                           xtickfontsize=8,
                           ytickfontsize=8;
                           kwargs...)
            for idx in 2:length(c)
                p = Plots.plot!([0, c[idx]], [0, 1],
                                projection=:polar,
                                color=:black)
            end
        else
            size(c, 2) > 2 && throw(ArgumentError("c must have exactly 2 columns: phases and lengths."))
            p = Plots.plot([0, c[1, 1]], [0, c[1, 2]],
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
                           linewidth=0.5,
                           titlefontsize=8,
                           xlabelfontsize=8,
                           ylabelfontsize=8,
                           xtickfontsize=8,
                           ytickfontsize=8;
                           kwargs...)
            for idx in 2:size(c, 1)
                p = Plots.plot!([0, c[idx, 1]], [0, c[idx, 2]],
                                projection=:polar,
                                color=:black)
            end
        end
    end

    Plots.plot(p)

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
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
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
    loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = loc_x[channel]
    loc_y = loc_y[channel]

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
                        markersize=5,
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
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
- `segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
- `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or grey palette
- `cb::Bool=true`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method: `:mean`, `:median`
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
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
function eeg_plot_topo(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, AbstractRange}=0, channel::Union{Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type])), segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg)), title::String="default", mono::Bool=false, cb::Bool=true, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, plot_size::Int64=800, head_labels::Bool=false, head_details::Bool=true, kwargs...)

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

    # get time vector
    if segment[2] <= eeg_epoch_len(eeg_tmp)
        signal = eeg_tmp.eeg_signals[channel, segment[1]:segment[2], 1]
    else
        signal = eeg_epochs(eeg_tmp, epoch_n=1).eeg_signals[channel, segment[1]:segment[2], 1]
    end
    t = _get_t(segment[1], segment[2], eeg_sr(eeg_tmp))
    t_1, t_s1, t_2, t_s2 = _convert_t(t[1], t[end])
    epoch = _t2epoch(eeg_tmp, segment[1], segment[2])
    
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
- `amethod::Symbol=:mean`: averaging method: `:mean`, `:median`
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
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
    t_1, t_s1, t_2, t_s2 = _convert_t(t[1], t[end])
    epoch = _t2epoch(eeg_tmp, segment[1], segment[2])
    
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
    eeg_plot_compose(p; <keyword arguments>)

Compose a complex plot of various plots contained in vector `p` using layout `layout`. Layout scheme is:
- `(2, 2)`: 2 × 2 plots, regular layout
- `@layout [a{0.2w} b{0.8w};_ c{0.6}]`: complex layout using Plots.jl `@layout` macro

# Arguments

- `p::Vector{Plots.Plot{Plots.GRBackend}}`: vector of plots
- `layout::Union(Matrix{Any}, Tuple{Int64, Int64}}`: layout
- `mono::Bool=false`: use color or grey palette
- `title::String=""`: plot title
- `kwargs`: optional arguments for `p` vector plots

# Returns

- `pc::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_compose(p::Vector{Plots.Plot{Plots.GRBackend}}; title::String="", layout::Union{Matrix{Any}, Tuple{Int64, Int64}}, mono::Bool=false, kwargs...)

    palette = mono == true ? :grays : :darktest

    pc = Plots.plot(grid=false,
                    framestyle=:none,
                    border=:none,
                    margins=0Plots.px)
    pc = Plots.plot!(p..., title=title, layout=layout, palette=palette; kwargs...)
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