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

eeg_add_component(eeg::NeuroJ.EEG; c::Symbol, v::Any)

eeg_add_component!(eeg::NeuroJ.EEG; c::Symbol, v::Any)

eeg_reset_components(eeg::NeuroJ.EEG)

eeg_reset_components!(eeg::NeuroJ.EEG)
```

## EEG process

```@docs
eeg_invert_polarity(eeg::NeuroJ.EEG; channel::Int64)

eeg_invert_polarity!(eeg::NeuroJ.EEG; channel::Int64)

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

eeg_channels_stats(eeg::NeuroJ.EEG)

eeg_channels_stats!(eeg::NeuroJ.EEG)

eeg_epochs_stats(eeg::NeuroJ.EEG)

eeg_epochs_stats!(eeg::NeuroJ.EEG)

eeg_spectrogram(eeg::NeuroJ.EEG; norm::Bool=true, demean::Bool=true)

eeg_spectrogram!(eeg::NeuroJ.EEG; norm::Bool=true, demean::Bool=true)

eeg_spectrum(eeg::NeuroJ.EEG; pad::Int64=0)

eeg_spectrum!(eeg::NeuroJ.EEG; pad::Int64=0)

eeg_s2t(eeg::NeuroJ.EEG; t::Int64)

eeg_t2s(eeg::NeuroJ.EEG; t::Union{Int64, Float64})

eeg_snr(eeg::NeuroJ.EEG)

eeg_snr!(eeg::NeuroJ.EEG)

eeg_standardize(eeg::NeuroJ.EEG)

eeg_standardize!(eeg::NeuroJ.EEG)
```

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

eeg_plot_spectrogram(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}, offset::Int64=0, len::Int64=0, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), xlabel::String="Time [s]", ylabel::String="Frequency [Hz]", title::String="", kwargs...)

signal_plot_histogram(signal::Vector{Float64}; type::Symbol=:hist, label::String="", xlabel::String="", ylabel::String="", title::String="", kwargs...)

signal_plot_histogram(signal::Union{Vector{Float64}, Matrix{Float64}}; type::Symbol=:hist, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", kwargs...)

eeg_plot_histogram(eeg::NeuroJ.EEG; type::Symbol=:hist, epoch::Int64=1, channel::Int64, offset::Int64=0, len::Int64=0, label::String="", xlabel::String="", ylabel::String="", title::String="", kwargs...)

eeg_plot_matrix(eeg::NeuroJ.EEG, m::Union{Matrix{Float64}, Array{Float64, 3}}; epoch::Int64=1, kwargs...)

eeg_plot_covmatrix(eeg::NeuroJ.EEG, cov_m::Union{Matrix{Float64}, Array{Float64, 3}}, lags::Union{Vector{Int64}, Vector{Float64}}; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch::Int64=1, kwargs...)

signal_plot_ica(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}, ica::Vector{Float64}; label::String="", norm::Bool=true, xlabel::String="Time [s]", ylabel::String="Amplitude [μV]", title::String="", ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0), kwargs...)

signal_plot_ica(t::Union{Vector{Float64}, Vector{Int64}, AbstractRange,}, ica::Matrix{Float64}; labels::Vector{String}=[""], norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

eeg_plot_ica(eeg::NeuroJ.EEG; epoch::Int64=1, offset::Int64=0, len::Int64=0, ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing, norm::Bool=true, xlabel::String="Time [s]", ylabel::String="", title::String="", kwargs...)

eeg_plot_topo(eeg::NeuroJ.EEG; offset::Int64, len::Int64=0, m::Symbol=:shepard, c::Symbol=:amp, c_idx::Union{Int64, Vector{Int64}, AbstractRange, Tuple, Nothing}=nothing, norm::Bool=true, frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0,0), head_labels::Bool=false, cb::Bool=false, cb_label::String="", average::Bool=true, title::String="", kwargs...)

signal_plot_band(signal::Vector{Float64}; fs::Int64, band::Vector{Symbol}=[:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], type::Symbol, norm::Bool=true, xlabel::String="", ylabel::String="", title::String="", kwargs...)

signal_plot_band(signal::Matrix{Float64}; fs::Int64, band::Symbol, type::Symbol, norm::Bool=true, labels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", kwargs...)

eeg_plot_band(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange}=1, channel::Union{Int64, Vector{Int64}, AbstractRange}, offset::Int64=0, len::Int64=0, band::Union{Symbol, Vector{Symbol}}=:all, type::Symbol, norm::Bool=true, xlabel::String="", ylabel::String="", title::String="", kwargs...)

eeg_plot_electrodes(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, selected::Union{Int64, Vector{Int64}, AbstractRange}=0, labels::Bool=true, head::Bool=true, head_labels::Bool=false, small::Bool=false, kwargs...)

eeg_draw_head(p::Plots.Plot{Plots.GRBackend}, loc_x::Vector{Float64}, loc_y::Vector{Float64}; head_labels::Bool=true, kwargs...)

eeg_plot_filter_response(eeg::NeuroJ.EEG; fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Tuple}, order::Int64, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, window::Union{Vector{Float64}, Nothing}=nothing, kwargs...)

eeg_plot_save(p::Plots.Plot{Plots.GRBackend}; file_name::String)
```

## Signal

```@docs
signal_derivative(signal::AbstractArray)

signal_derivative(signal::Array{Float64, 3})

signal_total_power(signal::AbstractArray; fs::Int64)

signal_total_power(signal::Array{Float64, 3}; fs::Int64)

signal_band_power(signal::AbstractArray; fs::Int64, f::Tuple)

signal_band_power(signal::Array{Float64, 3}; fs::Int64, f::Tuple)

signal_make_spectrum(signal::AbstractArray; fs::Int64)

signal_make_spectrum(signal::Array{Float64, 3}; fs::Int64)

signal_detrend(signal::AbstractArray; type::Symbol=:linear)

signal_detrend(signal::Array{Float64, 3}; type::Symbol=:linear)

signal_ci95(signal::Vector{Float64}; n::Int64=3, method::Symbol=:normal)

signal_ci95(signal::AbstractArray; n::Int64=3, method::Symbol=:normal)

signal_ci95(signal::Array{Float64, 3}; n::Int64=3, method::Symbol=:normal)

signal_mean(signal1::Vector{Float64}, signal2::Vector{Float64})

signal_mean(signal1::Array{Float64, 3}, signal2::Array{Float64, 3})

signal_difference(signal1::AbstractArray, signal2::AbstractArray; n::Int64=3, method::Symbol=:absdiff)

signal_difference(signal1::Array{Float64, 3}, signal2::Array{Float64, 3}; n::Int64=3, method::Symbol=:absdiff)

signal_autocov(signal::AbstractArray; lag::Int64=1, demean::Bool=false, norm::Bool=false)

signal_autocov(signal::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, norm::Bool=false)

signal_crosscov(signal1::AbstractArray, signal2::AbstractArray; lag::Int64=1, demean::Bool=false, norm::Bool=false)

signal_crosscov(signal::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, norm::Bool=false)

signal_crosscov(signal1::Array{Float64, 3}, signal2::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, norm::Bool=false)

signal_spectrum(signal::AbstractArray; pad::Int64=0)

signal_spectrum(signal::Array{Float64, 3}; pad::Int64=0)

signal_epochs(signal::Vector{Float64}; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

signal_epochs(signal::Matrix{Float64}; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

signal_delete_channel(signal::Matrix{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange})

signal_delete_channel(signal::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange})

signal_reference_channel(signal::Matrix{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange})

signal_reference_channel(signal::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange})

signal_reference_car(signal::Matrix{Float64})

signal_reference_car(signal::Array{Float64, 3})

signal_taper(signal::AbstractArray; taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

signal_taper(signal::Array{Float64, 3}; taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

signal_demean(signal::AbstractArray)

signal_demean(signal::Array{Float64, 3})

signal_normalize_zscore(signal::AbstractArray)

signal_normalize_zscore(signal::Array{Float64, 3})

signal_normalize_minmax(signal::AbstractArray)

signal_normalize_minmax(signal::Array{Float64, 3})

signal_cov(signal1::AbstractArray, signal2::AbstractArray; norm::Bool=false)

signal_cov(signal::AbstractArray; norm::Bool=false)

signal_cov(signal::Array{Float64, 3}; norm::Bool=false)

signal_cor(signal::Array{Float64, 3})

signal_add_noise(signal::AbstractArray)

signal_add_noise(signal::Array{Float64, 3})

signal_upsample(signal::AbstractArray; t::AbstractRange, new_sr::Int64)

signal_upsample(signal::Array{Float64, 3}; t::AbstractRange, new_sr::Int64)

signal_tconv(signal::AbstractArray; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

signal_tconv(signal::Array{Float64, 3}; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

signal_filter(signal::AbstractArray; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, fs::Int64=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

signal_filter(signal::Array{Float64, 3}; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, fs::Int64=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

signal_downsample(signal::AbstractArray; t::AbstractRange, new_sr::Int64)

signal_downsample(signal::Array{Float64, 3}; t::AbstractRange, new_sr::Int64)

signal_psd(signal::AbstractArray; fs::Int64, norm::Bool=false)

signal_psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false)

signal_psd(signal::Array{Float64, 3}; fs::Int64, norm::Bool=false)

signal_stationarity_hilbert(signal::AbstractArray)

signal_stationarity_mean(signal::AbstractArray; window::Int64)

signal_stationarity_var(signal::AbstractArray; window::Int64)

signal_stationarity(signal::Array{Float64, 3}; window::Int64=10, method::Symbol=:hilbert)

signal_trim(signal::AbstractArray; len::Int64, offset::Int64=1, from::Symbol=:start)

signal_trim(signal::Array{Float64, 3}; len::Int64, offset::Int64=1, from::Symbol=:start)

signal_mi(signal1::AbstractArray, signal2::AbstractArray)

signal_mi(signal::Array{Float64, 3})

signal_mi(signal1::Array{Float64, 3}, signal2::Array{Float64, 3})

signal_entropy(signal::AbstractArray)

signal_entropy(signal::Array{Float64, 3})

signal_average(signal1::AbstractArray, signal2::AbstractArray)

signal_average(signal::Array{Float64, 3})

signal_average(signal1::Array{Float64, 3}, signal2::Array{Float64, 3})

signal_coherence(signal1::AbstractArray, signal2::AbstractArray)

signal_coherence(signal1::Matrix{Float64}, signal2::Matrix{Float64})

signal_coherence(signal1::Array{Float64, 3}, signal2::Array{Float64, 3})

signal_pca(signal::Array{Float64, 3}; n::Int64)

signal_fconv(signal::AbstractArray; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

signal_fconv(signal::Array{Float64, 3}; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

signal_ica(signal::Array{Float64, 3}; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

signal_ica_reconstruct(signal::Array{Float64, 3}; ic_activations::Array{Float64, 3}, ic_mw::Array{Float64, 3}, ic_v::Union{Int64, Vector{Int64}, AbstractRange})

signal_epochs_stats(signal::Array{Float64, 3})

signal_channels_stats(signal::AbstractArray)

signal_channels_stats(signal::Array{Float64, 3})

signal_spectrogram(signal::AbstractArray; fs::Int64, norm::Bool=true, demean::Bool=true)

signal_spectrogram(signal::Array{Float64, 3}; fs::Int64, norm::Bool=true, demean::Bool=true)

signal_band(fs::Union{Int64, Float64}, band::Symbol)

signal_detect_epoch_flat(signal::Array{Float64, 3})

signal_detect_epoch_rmse(signal::Array{Float64, 3})

signal_detect_epoch_rmsd(signal::Array{Float64, 3})

signal_detect_epoch_euclid(signal::Array{Float64, 3})

signal_detect_epoch_p2p(signal::Array{Float64, 3})

signal_invert_polarity(signal::AbstractArray)

signal_invert_polarity(signal::Array{Float64, 3})

signal_snr(signal::AbstractArray)

signal_snr(signal::Array{Float64, 3})

signal_standardize(signal::Array{Float64, 3})
```

## Misc

```@docs
pad0m(m::Union{Matrix{Int64}, Matrix{Float64}, Matrix{ComplexF64}})

vsearch(y::Union{Int64, Float64}, x::Union{Vector{Int64}, Vector{Float64}}; return_distance::Bool=false)

vsearch(y::Union{Vector{Int64}, Vector{Float64}}, x::Union{Vector{Int64}, Vector{Float64}}; return_distance=false)

jaccard_similarity(x::Union{Vector{Int64}, Vector{Float64}}, y::Union{Vector{Int64}, Vector{Float64}})

fft0(x::AbstractArray, n::Int64)

ifft0(x::AbstractArray, n::Int64)

nextpow2(x::Int64)

vsplit(x::Union{Vector{Int64}, Vector{Float64}}, n::Int64=1)

freqs(t::Union{Vector{Int64}, Vector{Float64}, AbstractRange})

freqs(signal::Vector{Float64}, fs::Union{Int64, Float64})

matrix_sortperm(m::Matrix; rev::Bool=false, dims::Int64=1)

matrix_sort(m::Matrix, m_idx::Vector{Int64}; rev::Bool=false, dims::Int64=1)

pad0(x::Union{Vector{Int64}, Vector{Float64}}, n::Int64, sym::Bool=false)

generate_window(type::Symbol, n::Int64; even::Bool=false)

generate_sinc(t::AbstractRange=-2:0.01:2; f::Union{Int64, Float64}=1, peak::Union{Int64, Float64}=0, norm::Bool=true)

generate_morlet(fs::Int64, wt::Union{Int64, Float64}, wf::Union{Int64, Float64}; ncyc::Int64=5, complex::Bool=false)

generate_gaussian(fs::Int64, gt::Union{Int64, Float64}, gw::Union{Int64, Float64})

tuple_order(t::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}, rev::Bool=false)

rmse(signal1::Vector{Float64}, signal2::Vector{Float64})
```

## NSTIM

```@docs
tes_dose(current::Union{Int64, Float64}, pad_area::Union{Int64, Float64}, duration::Int64)
```