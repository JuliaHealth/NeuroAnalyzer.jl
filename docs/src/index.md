# NeuroJ.jl Documentation

This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).

## NeuroJ

```@docs
neuroj_version()

neuroj_reload_plugins()
```

---

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

---

## EEG edit

```@docs
eeg_delete_channel(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

eeg_delete_channel!(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

eeg_keep_channel(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

eeg_keep_channel!(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

eeg_get_channel(eeg::NeuroJ.EEG; channel::Union{Int64, String})

eeg_rename_channel(eeg::NeuroJ.EEG; channel::Union{Int64, String}, new_name::String)

eeg_rename_channel!(eeg::NeuroJ.EEG; channel::Union{Int64, String}, new_name::String)

eeg_extract_channel(eeg::NeuroJ.EEG; channel::Union{Int64, String})

eeg_history(eeg::NeuroJ.EEG)

eeg_labels(eeg::NeuroJ.EEG)

eeg_sr(eeg::NeuroJ.EEG)

eeg_channel_n(eeg::NeuroJ.EEG; type::Symbol=:all)

eeg_epoch_n(eeg::NeuroJ.EEG)

eeg_signal_len(eeg::NeuroJ.EEG)

eeg_epoch_len(eeg::NeuroJ.EEG)

eeg_info(eeg::NeuroJ.EEG)

eeg_epochs(eeg::NeuroJ.EEG; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

eeg_epochs!(eeg::NeuroJ.EEG; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

eeg_extract_epoch(eeg::NeuroJ.EEG; epoch::Int64)

eeg_trim(eeg::NeuroJ.EEG; len::Int64, offset::Int64=1, from::Symbol=:start, keep_epochs::Bool=true)

eeg_trim!(eeg::NeuroJ.EEG; len::Int64, offset::Int64=1, from::Symbol=:start, keep_epochs::Bool=true)

eeg_edit_header(eeg::NeuroJ.EEG; field::Symbol, value::Any)

eeg_edit_header!(eeg::NeuroJ.EEG; field::Symbol, value::Any)

eeg_show_header(eeg::NeuroJ.EEG)

eeg_delete_epoch(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

eeg_delete_epoch!(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

eeg_keep_epoch(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

eeg_keep_epoch!(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

eeg_detect_bad_epochs(eeg::NeuroJ.EEG; method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p], ch_t::Float64=0.1)

eeg_delete_bad_epochs(eeg::NeuroJ.EEG; bad_epochs::Vector{Int64}, confirm::Bool=true)

eeg_delete_bad_epochs!(eeg::NeuroJ.EEG; bad_epochs::Vector{Int64}, confirm::Bool=true)

eeg_add_labels(eeg::NeuroJ.EEG, labels::Vector{String})

eeg_add_labels!(eeg::NeuroJ.EEG, labels::Vector{String})

eeg_edit_channel(eeg::NeuroJ.EEG; channel::Int64, field::Any, value::Any)

eeg_edit_channel!(eeg::NeuroJ.EEG; channel::Int64, field::Any, value::Any)

eeg_keep_eeg_channels(eeg::NeuroJ.EEG)

eeg_keep_eeg_channels!(eeg::NeuroJ.EEG)

eeg_list_components(eeg::NeuroJ.EEG)

eeg_extract_component(eeg::NeuroJ.EEG; c::Symbol)

eeg_delete_component(eeg::NeuroJ.EEG; c::Symbol)

eeg_delete_component!(eeg::NeuroJ.EEG; c::Symbol)

eeg_reset_components(eeg::NeuroJ.EEG)

eeg_reset_components!(eeg::NeuroJ.EEG)
```

---

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

---

## EEG analyze

```@docs
eeg_total_power(eeg::NeuroJ.EEG)

eeg_total_power!(eeg::NeuroJ.EEG)

eeg_band_power(eeg::NeuroJ.EEG; f::Tuple)

eeg_cov(eeg::NeuroJ.EEG; norm=true)

eeg_cov!(eeg::NeuroJ.EEG; norm=true)

eeg_cor(eeg::NeuroJ.EEG)

eeg_cor!(eeg::NeuroJ.EEG)

eeg_autocov(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

eeg_autocov!(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

eeg_crosscov(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

eeg_crosscov!(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

eeg_crosscov(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

eeg_psd(eeg::NeuroJ.EEG; norm::Bool=false)

eeg_psd!(eeg::NeuroJ.EEG; norm::Bool=false)

eeg_stationarity(eeg::NeuroJ.EEG; window::Int64=10, method::Symbol=:hilbert)

eeg_stationarity!(eeg::NeuroJ.EEG; window::Int64=10, method::Symbol=:hilbert)

eeg_mi(eeg::NeuroJ.EEG)

eeg_mi!(eeg::NeuroJ.EEG)

eeg_mi(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG)

eeg_entropy(eeg::NeuroJ.EEG)

eeg_entropy!(eeg::NeuroJ.EEG)

eeg_band(eeg; band::Symbol)

eeg_coherence(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG)

eeg_coherence(eeg::NeuroJ.EEG; channel1::Int64, channel2::Int64, epoch1::Int64, epoch2::Int64)

eeg_freqs(eeg::NeuroJ.EEG)

eeg_freqs!(eeg::NeuroJ.EEG)

eeg_difference(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG; n::Int64=3, method::Symbol=:absdiff)

eeg_pick(eeg::NeuroJ.EEG; pick::Union{Symbol, Vector{Symbol}})

eeg_epochs_stats(eeg::NeuroJ.EEG)

eeg_epochs_stats!(eeg::NeuroJ.EEG)

eeg_spectrogram(eeg::NeuroJ.EEG; norm::Bool=true, demean::Bool=true)

eeg_spectrogram!(eeg::NeuroJ.EEG; norm::Bool=true, demean::Bool=true)

eeg_spectrum(eeg::NeuroJ.EEG; pad::Int64=0)

eeg_spectrum!(eeg::NeuroJ.EEG; pad::Int64=0)

eeg_s2t(eeg::NeuroJ.EEG; t::Int64)

eeg_t2s(eeg::NeuroJ.EEG; t::Union{Int64, Float64})
```

---

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

---

## Signal

---

## Misc

---

## NSTIM

```@docs
tes_dose(current::Union{Int64, Float64}, pad_area::Union{Int64, Float64}, duration::Int64)
```

---