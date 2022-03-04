# NeuroJ.jl Documentation

Welcome fellow researcher! NeuroJ.jl is a [Julia](https://julialang.org) package for analyzing of EEG data. Future versions will also process NIRS data and use MRI data for source localization techniques. Also, various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will be included.

This is a non-commercial projected, aimed for researchers in psychiatry, neurology and neuroscience.

Initially NeuroJ.jl will be focused on resting-state EEG analysis, but ERP and other type of analyses will be developed in future versions.

Every contribution (bug reports, fixes, new ideas, feature requests or additions, documentation improvements, etc.) to the project is highly welcomed.

## EEG plots

```@docs
signal_plot(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Vector{Float64}; ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", kwargs...)

signal_plot(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::AbstractArray; labels::Vector{String}=[""], xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

eeg_plot(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], xlabel::String="Time [s]", ylabel::String="", title::String="", head::Bool=true, hist::Symbol=:hist, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

signal_plot_avg(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Matrix{Float64}; norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

eeg_plot_avg(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), hist::Symbol=:hist, head::Bool=true, kwargs...)

signal_plot_butterfly(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, signal::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

eeg_plot_butterfly(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], norm::Bool=false, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="Butterfly plot", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), head::Bool=true, hist::Bool=true, average::Bool=false, kwargs...)

signal_plot_psd(s_powers::Vector{Float64}, s_freqs::Vector{Float64}; frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="", kwargs...)

signal_plot_psd(signal::Vector{Float64}; fs::Int64, norm::Bool=false, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel="Frequency [Hz]", ylabel="Power [μV^2/Hz]", title="", kwargs...)

signal_plot_psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false, average::Bool=false, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), labels::Vector{String}=[""], xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="", kwargs...)

eeg_plot_psd(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}=0, offset::Int64=0, len::Int64=0, labels::Vector{String}=[""], norm::Bool=false, average::Bool=false, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel::String="Frequency [Hz]", ylabel::String="Power [μV^2/Hz]", title::String="", head::Bool=false, kwargs...)

signal_plot_spectrogram(signal::Vector{Float64}; fs::Int64, offset::Int64=0, norm::Bool=true, demean::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel="Time [s]", ylabel="Frequency [Hz]", title="Spectrogram", kwargs...)

eeg_plot_spectrogram(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Int64, offset::Int64=0, len::Int64=0, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel::String="Time [s]", ylabel::String="Frequency [Hz]", title::String="", kwargs...)



eeg_plot_matrix(eeg::NeuroJ.EEG, m::Union{Matrix{Float64}, Array{Float64, 3}}; epoch::Int64=1, kwargs...)

eeg_plot_covmatrix(eeg::NeuroJ.EEG, cov_m::Union{Matrix{Float64}, Array{Float64, 3}}, lags::Union{Vector{Int64}, Vector{Float64}}; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch::Int64=1, kwargs...)

eeg_plot_electrodes(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head::Bool=true, head_labels::Bool=false, small::Bool=false, kwargs...)

eeg_draw_head(p::Plots.Plot{Plots.GRBackend}, loc_x::Vector{Float64}, loc_y::Vector{Float64}; head_labels::Bool=true, kwargs...)

eeg_plot_filter_response(eeg::NeuroJ.EEG; fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Tuple}, order::Int64, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, window::Union{Vector{Float64}, Nothing}=nothing, kwargs...)
```
