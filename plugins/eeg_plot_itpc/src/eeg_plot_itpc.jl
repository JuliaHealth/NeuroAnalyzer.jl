export eeg_plot_itpc
export eeg_plot_itpc_s
export eeg_plot_itpc_f

"""
    eeg_plot_itpc(eeg; <keyword arguments>)

Plot ITPC (Inter-Trial-Phase Clustering) at time `t` over epochs/trials of `channel` of `eeg`.

# Arguments

- `eeg:NeuroAnalyzer.EEG`
- `channel::Int64`: channel to plot
- `t::Int64`: time point to plot
- `z::Bool=false`: plot ITPCz instead of ITPC
- `w::Union{AbstractVector, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_itpc(eeg::NeuroAnalyzer.EEG; channel::Int64, t::Int64, z::Bool=false, w::Union{AbstractVector, Nothing}=nothing, mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    itpc, itpcz, itpc_angle, itpc_phases = eeg_itpc(eeg, channel=channel, t=t, w=w)
    itpc = round(itpc, digits=2)
    itpcz = round(itpcz, digits=2)
    t = eeg_s2t(eeg, t=t)

    p1 = Plots.plot(itpc_phases,
                    seriestype=:histogram,
                    bins=(length(itpc_phases) ÷ 10),
                    xticks=[-3.14, 0, 3.14],
                    xlims=(-pi, pi),
                    fill=:lightgrey,
                    title="Phase angles across trials\nchannel: $channel",
                    xlabel="Phase angle [rad]",
                    ylabel="Count/bin")

    if z == false
        if w === nothing
            p2 = @views Plots.plot([0, itpc_phases[1]], [0, 1], projection=:polar, yticks=false, color=:black, lw=0.2, legend=nothing, title="Phase differences\nITPC at $t s = $itpc")
        else
            p2 = @views Plots.plot([0, itpc_phases[1]], [0, 1], projection=:polar, yticks=false, color=:black, lw=0.2, legend=nothing, title="Phase differences\nwITPC at $t s = $itpc")
        end
    else
        if w === nothing
            p2 = @views Plots.plot([0, itpc_phases[1]], [0, 1], projection=:polar, yticks=false, color=:black, lw=0.2, legend=nothing, title="Phase differences\nITPCz at $t s = $itpcz")
        else
            nothing && (p2 = @views Plots.plot([0, itpc_phases[1]], [0, 1], projection=:polar, yticks=false, color=:black, lw=0.2, legend=nothing, title="Phase differences\nwITPCz at $t s = $itpcz"))
        end
    end
    for idx in 2:length(itpc_phases)
        p2 = @views Plots.plot!([0, itpc_phases[idx]], [0, 1], projection=:polar, color=:black, lw=0.2)
    end
    itpcz > 1 && (itpcz = 1)
    z == false && (p2 = Plots.plot!([0, itpc_angle], [0, itpc], lw=1, color=:red))
    z == true && (p2 = Plots.plot!([0, itpc_angle], [0, itpcz], lw=1, color=:red))

    p = Plots.plot(p1, p2,
                   legend=false,
                   titlefontsize=8,
                   xlabelfontsize=6,
                   ylabelfontsize=6,
                   xtickfontsize=4,
                   ytickfontsize=4,
                   margins=10Plots.px,
                   palette=pal;
                   kwargs...)

    return p
end

"""
    eeg_plot_itpc_s(eeg; <keyword arguments>)

Plot spectrogram of ITPC (Inter-Trial-Phase Clustering).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`: channel to plot
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:lin`: linear (:lin) or logarithmic (:log) frequencies
- `z::Bool=false`: plot ITPCz instead of ITPC
- `w::Union{AbstractVector, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Frequency [Hz]"`: y-axis label
- `title::String="ITPC spectrogram"`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_itpc_s(eeg::NeuroAnalyzer.EEG; channel::Int64, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:lin, z::Bool=false, w::Union{AbstractVector, Nothing}=nothing, xlabel::String="Time [s]", ylabel::String="Frequency [Hz]", title::String="default", mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest
        title == "default" && (title = "ITPC spectrogram\nchannel: $channel")
    
    itpc_s, itpc_z_s, frq_list = eeg_itpc_s(eeg, channel=channel, frq_lim=frq_lim, frq_n=frq_n, frq=frq)

    s = z == true ? itpc_z_s : itpc_s

    p = eeg_plot_spectrogram(eeg,
                             frq_list,
                             s,
                             title=title,
                             xlabel=xlabel,
                             ylabel=ylabel;
                             kwargs...)

    return p
end

"""
    eeg_plot_itpc_f(eeg; <keyword arguments>)

Plot time-frequency plot of ITPC (Inter-Trial-Phase Clustering) for `channel` of `eeg` for frequency `f`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`: channel to plot
- `f::Int64`: frequency to plot
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:lin`: linear (:lin) or logarithmic (:log) frequencies
- `z::Bool=false`: plot ITPCz instead of ITPC
- `w::Union{AbstractVector, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Frequency [Hz]"`: y-axis label
- `title::String="default"`: plot title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_itpc_f(eeg::NeuroAnalyzer.EEG; channel::Int64, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:lin, f::Int64, z::Bool=false, w::Union{AbstractVector, Nothing}=nothing, xlabel::String="Time [s]", ylabel::String="ITPC", title::String="default", mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest
    f < 0 && throw(ArgumentError("f must be ≥ 0."))
    f > eeg_sr(eeg) ÷ 2 && throw(ArgumentError("f must be ≤ $(eeg_sr(eeg) ÷ 2)."))

    itpc_s, itpc_z_s, frq_list = eeg_itpc_s(eeg, channel=channel, frq_lim=frq_lim, frq_n=frq_n, frq=frq, w=w)
    title == "default" && (title = "ITPC at frequency $(vsearch(f, frq_list)) Hz\nchannel: $channel")

    s = z == true ? itpc_z_s[vsearch(f, frq_list), :] : itpc_s[vsearch(f, frq_list), :]

    p = eeg_plot(eeg,
                 s,
                 title=title,
                 xlabel=xlabel,
                 ylabel=ylabel;
                 kwargs...)

    return p
end