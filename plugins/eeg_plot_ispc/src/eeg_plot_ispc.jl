export eeg_plot_ispc

"""
    eeg_plot_ispc(eeg1, eeg2; <keyword arguments>)

Plot ISPC `eeg1` and `eeg2` channels/epochs.

# Arguments

- `eeg1:NeuroAnalyzer.EEG`
- `eeg2:NeuroAnalyzer.EEG`
- `channel1::Int64`: channel to plot
- `channel2::Int64`: channel to plot
- `epoch1::Int64`: epoch to plot
- `epoch2::Int64`: epoch to plot
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function eeg_plot_ispc(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Int64, channel2::Int64, epoch1::Int64, epoch2::Int64, mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    ispc, ispc_angle, signal_diff, phase_diff, s1_phase, s2_phase = eeg_ispc(eeg1, eeg2, channel1=channel1, channel2=channel2, epoch1=epoch1, epoch2=epoch2)
    ispc = round.(ispc, digits=2)

    t = eeg1.eeg_epochs_time

    p1 = @views Plots.plot(t, eeg1.eeg_signals[channel1, :, epoch1], color=:black, lw=0.2)
    p1 = @views Plots.plot!(t, eeg2.eeg_signals[channel2, :, epoch2], color=:grey, lw=0.2, title="Signals", legend=false, xlabel="Time [s]", ylabel="Amplitude [μv]")

    p2 = @views Plots.plot(t, signal_diff[1, :, 1], color=:black, lw=0.2, title="Signals difference", legend=false, xlabel="Time [s]", ylabel="Amplitude [μv]")

    p3 = @views Plots.plot(t, s1_phase[1, :, 1], color=:black, lw=0.2)
    p3 = @views Plots.plot!(t, s2_phase[1, :, 1], color=:grey, lw=0.2, title="Phases", legend=false, xlabel="Time [s]", ylabel="Angle [rad]")

    p4 = @views Plots.plot(t, phase_diff[1, :, 1], color=:black, lw=0.2, title="Phases difference", legend=false, xlabel="Time [s]", ylabel="Angle [rad]")

    p5 = @views Plots.plot([0, s1_phase[1, :, 1][1]], [0, 1], projection=:polar, yticks=false, color=:black, lw=0.2, legend=nothing, title="Phases")
    @views for idx in 2:length(phase_diff[1, :, 1])
        p5 = @views Plots.plot!([0, s1_phase[1, :, 1][idx]], [0, 1], projection=:polar, color=:black, lw=0.2)
    end

    p5 = @views Plots.plot!([0, s2_phase[1, :, 1][1]], [0, 1], projection=:polar, yticks=false, color=:grey, lw=0.2, legend=nothing)
    @views for idx in 2:length(phase_diff[1, :, 1])
        p5 = @views Plots.plot!([0, s2_phase[1, :, 1][idx]], [0, 1], projection=:polar, color=:grey, lw=0.2)
    end

    p6 = @views Plots.plot([0, phase_diff[1, :, 1][1]], [0, 1], projection=:polar, yticks=false, color=:black, lw=0.2, legend=nothing, title="Phases difference and ISPC = $(ispc[1])")
    @views for idx in 2:length(phase_diff[1, :, 1])
        p6 = @views Plots.plot!([0, phase_diff[1, :, 1][idx]], [0, 1], projection=:polar, color=:black, lw=0.2)
    end
    p6 = @views Plots.plot!([0, ispc_angle[1]], [0, ispc[1]], lw=1, color=:red)
    
    p = Plots.plot(p1, p2, p3, p4, p5, p6,
                   layout=(3, 2),
                   titlefontsize=8,
                   xlabelfontsize=6,
                   ylabelfontsize=6,
                   xtickfontsize=4,
                   ytickfontsize=4,
                   palette=pal;
                   kwargs...)

    return p
end