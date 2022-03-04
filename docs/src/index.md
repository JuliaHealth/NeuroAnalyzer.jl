# NeuroJ.jl Documentation

This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).

## NeuroJ

```@docs
neuroj_version()

neuroj_reload_plugins()
```

## EEG io

```@docs
eeg_import_edf(file_name::String; read_annotations::Bool=true, clean_labels::Bool=true)

eeg_import_ced(file_name::String)

eeg_import_locs(file_name::String)

eeg_import_elc(file_name::String)

eeg_load_electrodes(eeg::NeuroJ.EEG; file_name::String)

eeg_load_electrodes!(eeg::NeuroJ.EEG; file_name::String)

eeg_load(file_name::String)

eeg_save(eeg::NeuroJ.EEG; file_name::String, overwrite::Bool=false)

eeg_export_csv(eeg::NeuroJ.EEG; file_name::String, header::Bool=false, components::Bool=false, overwrite::Bool=false)
```

## EEG edit

## EEG process

```@docs
eeg_reference_channel(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

eeg_reference_channel!(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

eeg_reference_car(eeg::NeuroJ.EEG)

eeg_reference_car!(eeg::NeuroJ.EEG)

eeg_derivative(eeg::NeuroJ.EEG)

eeg_derivative!(eeg::NeuroJ.EEG)

eeg_detrend(eeg::NeuroJ.EEG; type::Symbol=:linear)

eeg_detrend!(eeg::NeuroJ.EEG; type::Symbol=:linear)

eeg_taper(eeg::NeuroJ.EEG; taper::Vector)

eeg_taper!(eeg::NeuroJ.EEG; taper::Vector)

eeg_demean(eeg::NeuroJ.EEG)

eeg_demean!(eeg::NeuroJ.EEG)

eeg_normalize_zscore(eeg::NeuroJ.EEG)

eeg_normalize_zscore!(eeg::NeuroJ.EEG)

eeg_normalize_minmax(eeg::NeuroJ.EEG)

eeg_normalize_minmax!(eeg::NeuroJ.EEG)

eeg_average(eeg::NeuroJ.EEG)

eeg_average!(eeg::NeuroJ.EEG)

eeg_resample(eeg::NeuroJ.EEG; new_sr::Int64)

eeg_resample!(eeg::NeuroJ.EEG; new_sr::Int64)

eeg_upsample(eeg::NeuroJ.EEG; new_sr::Int64)

eeg_upsample!(eeg::NeuroJ.EEG; new_sr::Int64)

eeg_downsample(eeg::NeuroJ.EEG; new_sr::Int64)

eeg_downsample!(eeg::NeuroJ.EEG; new_sr::Int64)

eeg_filter(eeg::NeuroJ.EEG; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

eeg_filter!(eeg::NeuroJ.EEG; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

eeg_tconv(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

eeg_tconv!(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

eeg_fconv(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

eeg_fconv!(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

eeg_pca(eeg::NeuroJ.EEG; n::Int64)

eeg_pca!(eeg::NeuroJ.EEG; n::Int64)

eeg_ica(eeg::NeuroJ.EEG; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

eeg_ica!(eeg::NeuroJ.EEG; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

eeg_ica_reconstruct(eeg::NeuroJ.EEG; ica::Union{Int64, Vector{Int64}, AbstractRange})

eeg_ica_reconstruct!(eeg::NeuroJ.EEG; ica::Union{Int64, Vector{Int64}, AbstractRange})
```

## EEG analyze

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

signal_plot_histogram(signal::Vector{Float64}; type::Symbol=:hist, label::String="", xlabel::String="", ylabel::String="", title::String="", kwargs...)

signal_plot_histogram(signal::Union{Vector{Float64}, Matrix{Float64}}; type::Symbol=:hist, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", kwargs...)

eeg_plot_histogram(eeg::NeuroJ.EEG; type::Symbol=:hist, epoch::Int64=1, channel::Int64, offset::Int64=0, len::Int64=0, label::String="", xlabel::String="", ylabel::String="", title::String="", kwargs...)

eeg_plot_matrix(eeg::NeuroJ.EEG, m::Union{Matrix{Float64}, Array{Float64, 3}}; epoch::Int64=1, kwargs...)

eeg_plot_covmatrix(eeg::NeuroJ.EEG, cov_m::Union{Matrix{Float64}, Array{Float64, 3}}, lags::Union{Vector{Int64}, Vector{Float64}}; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch::Int64=1, kwargs...)

signal_plot_ica(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, ica::Vector{Float64}; label::String="", norm::Bool=true, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

signal_plot_ica(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange,}, ica::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

eeg_plot_ica(eeg::NeuroJ.EEG; epoch::Int64=1, offset::Int64=0, len::Int64=0, ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

eeg_plot_topo(eeg::NeuroJ.EEG; offset::Int64, len::Int64=0, m::Symbol=:shepard, c::Symbol=:amp, c_idx::Union{Int64, Vector{Int64}, AbstractRange, Tuple, Nothing}=nothing, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0,0), head_labels::Bool=false, cb::Bool=false, cb_label::String="", average::Bool=true, title::String="", kwargs...)

signal_plot_bands(signal::Vector{Float64}; fs::Int64, band::Union{Symbol, Vector{Symbol}}=:all, type::Symbol, norm::Bool=true, xlabel::String="", ylabel::String="", title::String="", kwargs...)

eeg_plot_bands(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Int64, offset::Int64=0, len::Int64=0, band::Union{Symbol, Vector{Symbol}}=:all, type::Symbol, norm::Bool=true, xlabel::String="", ylabel::String="", title::String="", kwargs...)

eeg_plot_electrodes(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head::Bool=true, head_labels::Bool=false, small::Bool=false, kwargs...)

eeg_draw_head(p::Plots.Plot{Plots.GRBackend}, loc_x::Vector{Float64}, loc_y::Vector{Float64}; head_labels::Bool=true, kwargs...)

eeg_plot_filter_response(eeg::NeuroJ.EEG; fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Tuple}, order::Int64, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, window::Union{Vector{Float64}, Nothing}=nothing, kwargs...)

eeg_plot_save(p::Plots.Plot{Plots.GRBackend}; file_name::String)
```

## Signal

## Misc

## NSTIM

```@docs
tes_dose(current::Union{Int64, Float64}, pad_area::Union{Int64, Float64}, duration::Int64)
```