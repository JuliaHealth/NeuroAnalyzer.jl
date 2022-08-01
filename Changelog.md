# NeuroJ.jl changes

Changes for v0.22.8 (01/08/2022)

- initial CUDA support
- new reference mode: simple (planar) Laplacian
- plugins use git
- many minor additions and fixes

# Git changelog

```
2a8a915 update: k_categories()
7bb6f58 update: seg_tcmp() -> seg_cmp() update: seg_cmp() use parametric or non-parametric test
c22dcff add: s_specseg() add: seg_tcmp()
abcf9bd add: s_cbp(), eeg_cbp(), eeg_cbp!() add: eeg_plot_electrode() update: new method for s_detrend(): HP filter fix: eeg_plot_electrode() margins update: eeg_plot_electrodes() code cleanup
3e7803b update: use fft0() and ifft0() insted of fft() and ifft()
c70a773 add: s_normalize_gauss() update: eeg_normalize() now use all available methods of normalization
5abdc38 update: ADF test added to eeg_stationarity()
3186f4e update: CUDA version in neuroj_version()
00f689f update: fft0() and ifft0() use CUDA if available
26381cb add: s_wbp() add: eeg_wbp() add: eeg_wbp!() update: threshold for eeg_senv()
4a223dd update: reduce edge artifacts in s_wspectrogram()
fee6e07 add: eeg_chdiff()
84dc71b add: eeg_zero() add: eeg_zero!()
6614910 add: eeg_fbsplit()
7287553 update: README.md
244a506 fix: for short signal (<4 * fs) always use multitaper PSD
d2ac664 fix: ticks for log x-axis in PSD plots fix: for short signal (<4 * fs) always use multitaper PSD
0ec3843 update: Tutorial.md
931a6df fix: lower frequency bound 0.1 hz for log xaxis PSD plots fix: ticks for log xaxis PSD plots
20bf018 fix: bug in plot_spectrogram() for mw=true update: log-lin, lin-log and log-log axes in PSD plots
f3adb41 update: eeg_plot_signal_psd() plots relative PSD
71ecbdc update: Tutorial.md
a936fa6 update: power at 0 Hz frequency = power at next frequency in PSD
6b6621c update: remove power at frequency = 0 Hz in PSD add: s_rel_psd() add: eeg_rel_psd() update: idx names in eeg_analyze.jl update: docstrings fix: s_total_power() docstring
4f9dabe add: eeg_loc_swap_axes() add: eeg_loc_swap_axes!()
415fe8d update: plugins use git add: neuroj_plugins_remove() add: neuroj_plugins_install() add: neuroj_plugins_update()
6414f37 update: s_signal_filter() uses signal-reflecting to reduce edge artifacts
965c849 add: eeg_band_mpower() add: s_band_mpower()
f13c468 update: eeg_interpolate_channel()
323d236 add: eeg_interpolate_channel() add: f_nearest()
20ee476 add: eeg_replace_channel() add: eeg_replace_channel!() add: eeg_add_note() add: eeg_add_note!() update: eeg_view_comment() -> eeg_view_note() add: eeg_delete_note() add: eeg_delete_note!()
ed8672e add: eeg_vartest()
545bdb7 add: eeg_normalize_log() add: eeg_normalize_log!() add: s_normalize_log()
099ac2e add: eeg_normalize_max() add: eeg_normalize_max!() add: s_normalize_max() update: grubbs()
47ee9f7 add: outlier_detect()
889445f update: statistic functions moved to statistic.jl submodule
1209dff add: grubbs()
492b760 add: infcrit()
073e8bf update: eeg_reference_slap() -> eeg_reference_plap()
5124b27 add: eeg_reference_slap() add: eeg_reference_slap!()
bebb495 add: brain.stl update: neuroj_version()
47a672e add: generate_morlet_fwhm()
1527c6a add: neuroj_plugins_list() update: neuroj_reload_plugins() -> neuroj_plugins_reload()
a40d910 update: README.md
8ac146b Merge branch 'master' into dev
eaf1591 update: versioning scheme changed to SymVer
03a7f2c update: coherence functions also return MSC
fdff918 add: s_gfp(), s_gfp_norm(), s2_diss()
a416bb9 Merge remote-tracking branch 'origin/dev' into dev
a85a439 Merge remote-tracking branch 'origin/dev'
fbc6b4b update: README.md
dc1b51c first official stable monthly release: 2022.06 add: s_cums()
481c59f update: install required packages
2fb582a add: a2_l2()
888eb3b add: plot_psd_3ds() update: eeg_plot_signal_psd_3d() option to plot 3d surface
d7dcdea add: s_fcoherence() add: eeg_fcoherence()
b4a7326 add: s_cmp()
34f4fb3 update: generate_window() :exp and :triangle windows added update: s_filter() :mavg and :mmed filters may use weighting window
850241e update: :mavg and :mmed filters use _reflect() and _chop()
9e9bba6 add: plot_psd_3d() add: eeg_plot_signal_psd_3d() update: version: 0.0.26
a599429 Merge remote-tracking branch 'origin/dev' into dev
117a3ff add: s_tkeo() add: eeg_tkeo() add: s_wspectrum() add: eeg_wspectrum()
45e49e7 update: coif2, coif4, coif8 wavelets added to s_wdenoise()
6087b26 update: eeg_reference_car() and eeg_reference_car!() options to exclude Fp1, Fp2, O1, O2 and current electrode from common average calculation
c3bd9ce Merge remote-tracking branch 'origin/dev' into dev
6f31036 update: new filters added to s_filter(): iirnotch and remez update: version: 0.0.25
3197511 update: s_gfilter() update: s_ghspectrogram()
5b6f5ce add: s_gfilter() add: s_ghspectrogram()
50d36e6 update: s_wdenoise(), eeg_wdenoise() add: s_fftdenoise(), eeg_fftdenoise()
1e78fe2 Merge remote-tracking branch 'origin/dev' into dev
e73b292 update: plot_spectrogram() option to use Morlet wavelet convolution
36a83e3 update: ITPCz, wITPC
195640d add: t2f() add: f2t()
d535ef8 add: eeg_itpc_s() add: eeg_plot_itpc_f() update: s_itpc() speed improved by 25% update: eeg_plot_spectrogram_itpc() -> eeg_plot_itpc_s()
29e4519 Merge remote-tracking branch 'origin/dev' into dev
6868b64 add: s_hspectrum() fix: s_ispc(), s_itpc(), s_pli() and plots now use Hilbert transform for phase angles fix: pad0() for n=0 update: eeg_spectrum() option to use Hilbert transform or FFT
eba5262 add: s_hspectrum(signal::AbstractArray; pad::Int64=0) fix: s_ispc(), s_itpc(), s_pli() and plots now use Hilbert transform for phase angles fix: pad0() for n=0 update: eeg_spectrum() option to use Hilbert transform or FFT
97b33c4 add: mono parameter for plots
aa901a6 add: s_frqinst() add: eeg_frqinst() add: eeg_ged()
d7d8b3e add: eeg_plot_connections() update: ear shapes
8e89807 add: eeg_aec()
a9fcbcb add: eeg_plot_spectrogram_itpc() fix: eeg_keep_epoch(), eeg_keep_epoch!() error while checking epoch value
b56040f add: s_pli(), eeg_pli(), eeg_plot_pli() add: eeg_pli_m() add: eeg_stpc_m()
d463072 update: s_itpc(), eeg_itpc(), eeg_plot_itpc()
d60df31 update: s_ispc(), eeg_ispc(), eeg_plot_ispc() over trials/epochs
9f61e37 add: s_ispc() add: eeg_ispc() add: eeg_plot_ispc() update: version 0.0.24
f01759a update: tutorial.md
b3ea2f6 update: eeg_plot_signal() default number of channels is 10 update: eeg_plot_signal() option to plot channels in mono
ca60830 add: head and brain images
f923343 add: effsize()
5fed213 update: version 0.0.23 fix: eeg_reference_m() and eeg_reference_a() use eeg_pick() for :c and :i referencing
52953fb add: eeg_reference_a(), eeg_reference_a!() add: eeg_reference_m(), eeg_reference_m!() fix: eeg_delete_channel!() fix: eeg_epochs!() update: eeg channel types update: eeg_reference_channel() -> eeg_reference_ch() update: eeg_reference_channel!() -> eeg_reference_ch!()
05d09c7 update: mutators return nothing
1f4fe23 fix: generate_sinc() NaN detection update: generate_morlet() and generate_gaussian()
054db76 update: README.md
20007b1 Merge remote-tracking branch 'origin/dev' into dev
2b5b03c update: README.md
cea5ad8 add: s_wt_denoise() add: eeg_wt_denoise()
d50380c update: eeg_plot_filter_response() plots amplitude and phase response for FIR filters
d32ca23 add: plot_signal_scaled()
0359384 add: multi-taper periodograms and spectrograms
7f283ab update: tutorial.md
d9af957 update: README.md
2e3950a add: eeg_negentropy() add: s_negentropy() add: plot_signal_scaled() update: eeg_epochs_stats() update: eeg_channels_stats()
34097b9 update: plot_signal() now plots unscaled signal
ff1fab4 add: eeg_plot_env()
897ed70 update: version 0.0.21 add: eeg_tenv_mean() add: eeg_tenv_median() add: eeg_penv_mean() add: eeg_penv_median() add: eeg_senv_mean() add: eeg_senv_median()
7fe11f0 add: eeg_tenv() add: eeg_penv() add: eeg_senv()
2a9476c add: s_findpeaks() add: eeg_tenv()
65ebd45 update: documentation
2d5ba7f update: plotting functions rewritten update: version 0.0.20
c713eb1 update: tutorial.md
b9eb498 add: eeg_plot_tile()
787253f update: cb parameter for _topo() plots
b04f468 update: s_pca() check for maximum number of PC
91e18a2 Merge remote-tracking branch 'origin/dev' into dev
495daf9 add: eeg_plot_ica_topo()
326c74b add: eeg_plot_ica_topo
1d9dc89 add: eeg_plot_weights_topo() add: eeg_plot_mcomponent_topo() update: eeg_plot_component_topo() -> eeg_plot_acomponent_topo()
2e14668 add: s_pca_reconstruct() add: eeg_pca_reconstruct() add: eeg_pca_reconstruct!()
6b215f8 add eeg_plot_component_spectrogram() add eeg_plot_component_spectrogram_avg()
b36fd50 add eeg_plot_component_idx_psd() add eeg_plot_component_idx_psd_avg() add eeg_plot_component_idx_psd_butterfly()
a7e68bb add: eeg_plot_component_idx() add: eeg_plot_component_idx_avg() add: eeg_plot_component_idx_butterfly()
f902e8f add: eeg_plot_signal_spectrogram_avg() add: eeg_plot_component_spectrogram_avg() add: eeg_plot_signal_topo()
b5ef694 Merge remote-tracking branch 'origin/dev' into dev
cf95722 add: eeg_plot_component_spectrogram() update: function inputs types simplified
25dda4e update: eeg_plot_sigal_spectrogram()
8087b93 add: eeg_plot_component_psd() add: eeg_plot_component_psd_avg() add: eeg_plot_component_psd_psd()
54c64ca fix: signal_plot_psd() frequencies and powers were misplaced add: eeg_plot_signal_psd_avg() add: eeg_plot_signal_psd_butterfly()
0ce5f42 add: eeg_plot_component_butterfly()
5e36f8c add: eeg_plot_signal_butterfly_details()
16b95d5 Merge remote-tracking branch 'origin/dev' into dev
b231f49 add: eeg_plot_signal_avg_details()
6edf85c add: eeg_plot_signal_avg_details(edf)
93aaf33  add: eeg_plot_signal_details()
79ed2e4 add: eeg_plot_compose()
213b382 add: eeg_plot_component()
d14efc7 update: eeg_plot()
06ffc1e update: eeg now contains time vector for each epoch (eeg.eeg_epochs_time) update: epochs time can start from negative value add: eeg_epochs_time() add: eeg_epochs_time!()
e33aa3d rewriting of eeg functions, removal of intermediate signal functions
dfb4fd8  update: version: 0.0.19  rewriting of eeg functions, removal of intermediate signal functions  update: documentation, README.md  add: internal low-level functions
b20c578 rewriting of eeg functions, removal of intermediate signal functions
fc6b0d2 update: many mutators removed add: eeg_rename_component() add: eeg_component_type() add: eeg_component_idx()
485e8f0 update: documentation
af58677 update: version 0.0.18 add: eeg_plot_epochs() update: eeg_plot_channels() plots external values and embedded components
933e0c4 add: eeg_plot_channels()
38690fb add: eeg_standardize(), eeg_standardize!() add: signal_standardize()
02326ee add: signal_snr() add: eeg_snr() add: eeg_snr!()
7b915cf fix: generate_window() add: generate_window() optional generating even-long window add: _pl()
dea4f9d update: performance improvement for signal_crosscov() and signal_autocov()
42dfd02 update: eeg_import_edf() detects channel type from channel name update: pad0() pads symmetrically or at the end update: generate_sinc() generates normalized or unnormalized sinc update: signal_upsample() uses DSP.resample()
ef22e95 update: generate_window() add: window types: blackmanharris, bohman, flattop, nutall update: tutorial
3b3d413 add: signal_channels_stats(), eeg_channels_stats() update: signal_epochs_stats(), eeg_epochs_stats() update: :fir filter default window length update: rp and rs defaults for signal_filter()
bbcd35e add: eeg_comment() add: eeg_delete() update: README.md
f16e652 update: tutorial
06960a5 add: eeg_plot_spectrogram() now plots single- (time vs. frq) and multi-channel (channels vs frq) signals
d5acecc add: eeg_add_component() add: eeg_add_component!() add: signal_invert_polarity() add: eeg_invert_polarity()
b4c32de fix: signal_detrending() add: :poly and :loess detrending
a4bbb32 fix: time axis in plots
e7b026c fix: epoch markers positioning
85f9ca5 update: eeg_plot_band() update: Tutorial.md
a5a7816 update: Tutorial.md
58a9254 update: channel numbers/labels are better formatted/positioned
2cf62a7 update: version 0.0.17
7a1793c update: repository moved to codeberg.org
98b6bcd Merge remote-tracking branch 'origin/master'
478c02e update: docstrings to follow Julia recommendations
73aa7ec update: Tutorial.md
fd3a5b3 update: documentation
df70726 add: Markdown documentation started
2757ad5 add: HTML documentation started
98bec7e update: eeg_ functions splitted into eeg_edit.jl, eeg_processing.jl and eeg_analysis.jl submodules
5cc8787 update: EEG is check for non-eeg channels update: code cleanup
22b49cf Merge remote-tracking branch 'origin/dev' into dev
ec0d472 update: version 0.0.16 add: eeg_edit_channel() add: eeg_edit_channel!() add: eeg_keep_eeg_channels() add: eeg_keep_eeg_channels!() update: eeg locations cleaned up
2e21049 update: massive performance increase
578b3a7 add: eeg_import_ced() add: eeg_import_elc() add: eeg_import_locs() add: pol2cart_sph() add: spherical locs
14ca64f add: plugins framework add: neuroj_reload_plugins()
6650811 add: rmse() add: signal_detect_epoch_flat(),signal_detect_epoch_rmse(), signal_detect_epoch_rmsd(), signal_detect_epoch_euclid(), signal_detect_epoch_p2p() add: eeg_detect_bad_epochs() add: eeg_delete_bad_epochs(), eeg_delete_bad_epochs!()
bc902ca update: version 0.0.15 add: eeg_ica_reconstruct(), eeg_ica_reconstruct!(), signal_ica_reconstruct()
d33ecc2 update: eeg_epochs_stats(), signal_epoch_stats() calculate median
4d3a2be update: eeg_plot_topo() now plots averaged or butterfly signal and PSD
4870e1a update: eeg_plot_ica() now plots weights distribution update: eeg_plot_electrodes() accepts no selected electrodes
dc0046a update: error messages more detailed update: code cleanup
07722e5 update: code cleanup
208094b fix: eeg_trim(), eeg_trim!(), signal_trim() fix: epoch_markers collecting add: tuple_order()
66c3732 update: version 0.0.14
d4adc35 update: version 0.0.13 update: eeg_plot_topo() accepts unweighted amplitude, ICA, PCA, power update: eeg_band() add: signal_band() update: vsearch() update: signal_ica() update: signal_pca() update: signal_normalize_minmax()
4923eb5 add: signal_plotbands(), eeg_plot_bands()
9d1f1aa add: eeg_spectrum(), eeg_spectrum!() update: eeg_plot() and eeg_plot_avg() now plot phase histogram
b904ce1 update: for eeg_plot_*() epoch as range or offset+len update: eeg_plot_avg() rewritten for single/multichannel plots add: eeg_average() update: plots font size update: version 0.0.12
c996cf5 update: code cleanup
35b763b update: eeg_plot() rewritten for single / multichannel plots add: epoch markers in plots
5e9e7ca add: eeg_plot_histogram() add: signal_plot_histogram() update: eeg_plot() also plots histograms
0f1dcba update: eeg_plot*() segment length set in samples
da8d2ff add: signal_spectrogram() add: eeg_spectrogram() update: signal_ica() parameters iter, f
d1ec111 fix: xticks spacing
0eacce6 add: eeg_ mutator (!) functions add: eeg_list_components() add: eeg_delete_component() add: eeg_extract_component() add: eeg_reset_components()
f46995e add: signal_plot_ica() update: mutators replace component if already present
d99511a add: components (e.g. ICA, PCA, PSD) are stored within EEG object update: version 0.0.12
5d557c6 add: eeg_ica() add: signal_epochs_stats(), eeg_epoch_stats()
d307f50 add: eeg_export_csv()
b55de0a add: eeg_pick()
4e6b19e add: eeg_delete_epoch() add: eeg_keep_epoch() update: eeg_save() calculates file size and keeps new name fix: eeg_delete_channel() deleting all channels
d36e1b9 add: eeg_edit()
90dd78c update: memory allocation
aa3b4de update: version 0.0.11 add: kwargs for plotting functions update: code cleanup
f9d3aa3 fix: eeg_crosscov()
4dd2105 update: multi-threading
0013c70 update: multi-threading
f694efb add: signal_pca(), eeg_pca() update: code cleanup
a7d3f8e add: signal_fconv(), eeg_fconv()
a30221e update: nfft calculation for signal_plot_spectrogram()
8121f18 add: signal_plot_spectrogram(), eeg_plot_spectrogram()
fa71c5b add: eeg_difference()
66a0d40 update: version 0.0.10
308aa9a add: signal_coherence(), eeg_coherence()
c638238 add: polynomial filter fix: eeg duration values in header for eeg_upsample() and eeg_downsample() update: mavg and mmed filters calculated for the whole signal
d759ddc add: threshold for :mavg and :mmed filters
f3701a2 add: signal_average(), eeg_average() add: generate_gausian() add: mavg filter with window
ec00a45 add: mavg and mmed filters
d7c7ef6 add: signal_mi(), eeg_mi()
d7130ff add: plot epoch markers
dd3543c add: offset parameter for signal_trim(), eeg_trim()
47c883b add: signal_trim(), eeg_trim()
8d10f3d add: signal_stationarity_hilbert() add: signal_stationarity_mean() add: signal_stationarity_var() add: signal_stationarity() add: eeg_stationarity()
cef7179 add: filter direction: oneway, oneway_reverse, twoway
3a57a36 fix: eeg_rename_channel() labels of the parent object also renamed
45c2d46 fix: lot of bugfixes update: all missing tests added update: __doc__ returns
eb67e35 add: tests for misc.jl
db7f425 update: tests added
5ea6c6d add: eeg_plot_matrix()
39c303a update: version 0.0.7
200cff2 add: eeg_load_electrode_positions() eeg_plot_electrodes()
0e25bdb update: plot line colors
e592c69 add: signal_plot_psd, eeg_plot_psd()
a2a34f9 add: signal_psd() eeg_psd()
b7a02ef add: signal_plot_butterfly(), eeg_plot_butterfly()
2da614f fix: time ticks
011b25a update: x axis time display on plot
f98c260 add: signal_plot_avg() and eeg_plot_avg()
f492932 add: average and butterfly plots
028c01c add: eeg_keep_channel() update: eeg_drop_channel() -> eeg_delete_channel()
c727b24 update: error handling while saving plots
ad7e0e3 add: multi-threading for epoch processing
b2bd74a add: eeg_crosscov() eeg_autocov()
f84405e add: signal_downsample() eeg_downsample()
1725df5 add: filter_response() phase response
bb7a1a2 fix: eeg_get_epoch()
8e4014d add: Chebyshev1, Chebyshev2, Elliptic filters
8672754 update: filter_response()
bd8a7d4 fix: filter_response()
612a3f4 Merge remote-tracking branch 'origin/master'
7ffdc91 fix: filter_response()
0da9b36 add: signal_filer() plots filter response fix: eeg_get_channel() should return vector
471ee16 update: eeg_header[:reference]
fc3a45a fix: signal_filter() window length for :fir
2ec7529 update
820df0f add: FIR filter update: eeg_filter_butter() -> eeg_filter(), signal_filter_butter() -> signal_filter() update: frequencies()
bb182c7 Merge remote-tracking branch 'origin/master'
fb96869 fix: fft0() for complex numbers update: ifft0(), fft0() now pads to n or up to n
3e6ac3e fix: fft0 for complex numbers update: ifft0, fft0 for pads to n or up to n
34ed62b add: morlet() add: signal_tconv(), eeg_tconv() update: plotting functions moved to plots.jl update: __doc__ cleanup update: function types cleanup update: eeg_show_processing_history() -> eeg_history() update: tutorial.md
a106adc generated package NeuroJ
ca4489b add: eeg.jl add: tes.jl
78ec22f first commit
```