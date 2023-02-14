using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

eeg = eeg_import_bdf("eeg-test-bdfplus.bdf")
eeg_delete_marker!(eeg, n=1)
@test size(eeg.eeg_markers) == (1, 5)
eeg_add_marker!(eeg, id="event", start=1, len=1, desc="test", channel=0)
@test size(eeg.eeg_markers) == (2, 5)
eeg_edit_marker!(eeg, n=2, id="event2", start=1, len=1, desc="test2", channel=0)

eeg = eeg_import_edf("eeg-test-edf.edf")
@test size(eeg.eeg_signals) == (24, 309760, 1)

ecg = eeg_extract_channel(eeg, channel=24)
eog2 = eeg_extract_channel(eeg, channel=23)
eog1 = eeg_extract_channel(eeg, channel=22)
eeg_delete_channel!(eeg, channel=22:24)

eeg1 = eeg_reference_a(eeg)
@test size(eeg1.eeg_signals) == (21, 309760, 1)
a1 = eeg_extract_channel(eeg, channel=20)
a2 = eeg_extract_channel(eeg, channel=21)
eeg_delete_channel!(eeg, channel=[20, 21])

eeg1 = eeg_delete_channel(eeg, channel=1)
@test eeg1.eeg_header[:channel_n] == 18

eeg1 = eeg_keep_channel(eeg, channel=1)
@test eeg1.eeg_header[:channel_n] == 1

eeg1 = eeg_derivative(eeg)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

tbp = eeg_total_power(eeg)
@test size(tbp) == (19, 1)

abp = eeg_band_power(eeg, f=(2, 4))
@test size(abp) == (19, 1)

eeg1 = eeg_detrend(eeg)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

eeg1 = eeg_reference_ch(eeg, channel=1)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

eeg1 = eeg_reference_car(eeg)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

eeg1 = eeg_extract_channel(eeg, channel="Cz")
@test size(eeg1) == (1, 309760, 1)

eeg1 = eeg_extract_channel(eeg, channel=18)
@test size(eeg1) == (1, 309760, 1)

@test eeg_get_channel(eeg, channel=1) == "Fp1"
@test eeg_get_channel(eeg, channel="Fp1") == 1

eeg1 = eeg_rename_channel(eeg, channel="Cz", name="CZ")
@test eeg1.eeg_header[:labels][18] == "CZ"
eeg1 = eeg_rename_channel(eeg, channel=1, name="FP1")
@test eeg1.eeg_header[:labels][1] == "FP1"

eeg1 = eeg_taper(eeg, taper=eeg.eeg_signals[1, :, 1])
@test size(eeg1.eeg_signals) == (19, 309760, 1)

eeg1 = eeg_demean(eeg)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

eeg1 = eeg_normalize(eeg, method=:zscore)
@test size(eeg1.eeg_signals) == (19, 309760, 1)
eeg1 = eeg_normalize(eeg, method=:minmax)
@test size(eeg1.eeg_signals) == (19, 309760, 1)
eeg1 = eeg_normalize(eeg, method=:log)
@test size(eeg1.eeg_signals) == (19, 309760, 1)
eeg1 = eeg_normalize(eeg, method=:gauss)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

cov_m = eeg_cov(eeg)
@test size(cov_m) == (19, 19, 1)

cor_m = eeg_cor(eeg)
@test size(cor_m) == (19, 19, 1)

eeg1 = eeg_upsample(eeg, new_sr=512)
@test size(eeg1.eeg_signals) == (19, 619519, 1)

@test typeof(eeg_history(eeg)) == Vector{String}

@test eeg_labels(eeg)[1] == "Fp1"

@test eeg_sr(eeg) == 256

eeg1 = eeg_epoch(eeg, epoch_len=1000)
eeg_erp!(eeg1)
@test size(eeg1.eeg_signals) == (19, 1000, 1)

eeg10 = eeg_epoch(eeg, epoch_n=10)
eeg1 = eeg_extract_epoch(eeg, epoch=1)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

f, s = eeg_dft(eeg)
@test size(f) == (19, 309760, 1)

m, _, _, _ = eeg_msci95(eeg)
@test size(m) == (1, 309760)

m, _, _, _ = eeg_mean(eeg, eeg)
@test m == zeros(1, 309760)

s, ss, p = eeg_difference(eeg, eeg)
@test p == [1.0]

eeg1 = eeg_filter(eeg, fprototype=:butterworth, ftype=:lp, cutoff=2, order=8)
@test size(eeg1.eeg_signals) == (19, 309760, 1)
eeg1 = eeg_filter(eeg, fprototype=:mavg, order=10)
@test size(eeg1.eeg_signals) == (19, 309760, 1)
eeg1 = eeg_filter(eeg, fprototype=:mmed, order=10)
@test size(eeg1.eeg_signals) == (19, 309760, 1)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

eeg1 = eeg_downsample(eeg, new_sr=128)
@test size(eeg1.eeg_signals) == (19, 154880, 1)

acov_m, _ = eeg_acov(eeg)
@test size(acov_m) == (19, 3, 1)
xcov_m, _ = eeg_xcov(eeg)
@test size(xcov_m) == (361, 3, 1)

p, f = eeg_psd(eeg)
@test size(p, 1) == 19

p = eeg_stationarity(eeg, method=:mean)
@test size(p) == (19, 10, 1)
p = eeg_stationarity(eeg, method=:var)
@test size(p) == (19, 10, 1)
p = eeg_stationarity(eeg, method=:hilbert)
@test size(p) == (19, 309759, 1)
p = eeg_stationarity(eeg, window=10000, method=:cov)
@test size(p) == (32, 1)

e = eeg_trim(eeg, segment=(10 * eeg_sr(eeg), 20 * eeg_sr(eeg)), remove_epochs=false)
@test eeg_signal_len(e) == 307199

m = eeg_mi(eeg)
@test size(m) == (19, 19, 1)
m = eeg_mi(eeg, eeg)
@test size(m) == (19, 19, 1)

e = eeg_entropy(eeg)
@test length(e) == 3
e = eeg_negentropy(eeg)
@test size(e) == (19, 1)

a = eeg_band(eeg, band=:alpha)
@test a == (8, 13)

c, msc, ic = eeg_tcoherence(eeg, eeg)
@test size(c) == (19, 309760, 1)

hz, nyq = eeg_freqs(eeg)
@test nyq == 128.0
@test length(hz) == 154880

e10 = eeg_epoch(eeg, epoch_len=2560)
s_conv = eeg_fconv(e10, kernel=generate_window(:hann, 256))
@test size(s_conv) == (19, 2560, 121)
s_conv = eeg_tconv(e10, kernel=generate_window(:hann, 256))
@test size(s_conv) == (19, 2560, 121)

p, v, m, pca = eeg_pca(eeg, n=2)
@test size(p) == (2, 309760, 1)
@test size(v) == (2, 1)
e1 = eeg_add_component(eeg, c=:pc, v=p)
eeg_add_component!(e1, c=:pca, v=pca)
e2 = eeg_pca_reconstruct(e1)
e2 = eeg_pca_reconstruct(eeg, p, pca)
@test size(e2.eeg_signals) == (19, 309760, 1)

e = eeg_edit_header(eeg, field=:patient, value="unknown")
@test e.eeg_header[:patient] == "unknown"

e = eeg_epoch(eeg, epoch_n=10)
e9 = eeg_delete_epoch(e, epoch=10)
@test size(e9.eeg_signals) == (19, 30976, 9)
e1 = eeg_keep_epoch(e, epoch=1)
@test size(e1.eeg_signals) == (19, 30976, 1)

@test length(eeg_channel_pick(eeg, pick=:left)) == 8

e = eeg_epoch(eeg, epoch_len=20*256)
v = eeg_epoch_stats(e)
@test length(v) == 10

e = eeg_epoch(eeg, epoch_len=20)
eeg_erp!(e)
i, _ = eeg_ica(e, n=5, tol=1.0)
@test size(i) == (5, 20, 1)

e = eeg_copy(eeg)
e_stats = eeg_epoch_stats(e)
@test length(e_stats) == (10)
eeg_add_component!(e, c=:epochs_mean, v=e_stats[1])
v = eeg_extract_component(e, c=:epochs_mean)
@test size(v) == (1, )
eeg_rename_component!(e, c_old=:epochs_mean, c_new=:epochs_m)
c = eeg_list_components(e)
@test size(c) == (1, )
c = eeg_component_type(e, c=:epochs_m)
@test c == Vector{Float64}
eeg_delete_component!(e, c=:epochs_m)
c = eeg_list_components(e)
@test size(c) == (0, )
eeg_reset_components!(e)
c = eeg_list_components(e)
@test size(c) == (0, )

e = eeg_epoch(eeg, epoch_len=2560)
eeg_erp!(e)
p, f, t = eeg_spectrogram(e)
@test size(p) == (1281, 37, 19, 1)
p, f, t = eeg_spectrogram(e, method=:mt)
@test size(p) == (257, 15, 19, 1)
p, f, t = eeg_spectrogram(e, method=:mw)
@test size(p) == (129, 2560, 19, 1)
p, f, t = eeg_spectrogram(e, method=:stft)
@test size(p) == (1281, 37, 19, 1)
p, f, t = eeg_spectrogram(e, method=:gh)
@test size(p) == (129, 2560, 19, 1)
p, f, t = eeg_spectrogram(e, method=:cwt)
@test size(p) == (18, 2560, 19, 1)

f, a, p, ph = eeg_spectrum(e)
@test size(p) == (19, 1280, 1)

e = eeg_copy(eeg)
i, iw = eeg_ica(e, tol=1.0, n=10)
eeg_add_component!(e, c=:ica, v=i)
eeg_add_component!(e, c=:ica_mw, v=iw)
@test size(e.eeg_components[1]) == (10, 309760, 1)
e2 = eeg_ica_reconstruct(e, ic=1)
@test size(e2.eeg_signals) == (19, 309760, 1)

b = eeg_detect_bad(eeg)
@test length(b) == 2

@test eeg_t2s(eeg, t=10) == 2561
@test eeg_s2t(eeg, t=10) == 0.04

e = eeg_keep_channel_type(eeg)
@test size(e.eeg_signals) == (19, 309760, 1)
eeg_edit_channel!(e, channel=19, field=:channel_type, value="ecg")
eeg_keep_channel_type!(e, type=:eeg)
@test size(e.eeg_signals) == (18, 309760, 1)

e = eeg_invert_polarity(eeg, channel=1)
@test e.eeg_signals[1, 1, 1] == -eeg.eeg_signals[1, 1, 1]

v = eeg_channel_stats(eeg)
@test length(v) == 10

eeg = eeg_import_edf("eeg-test-edf.edf")
eeg_delete_channel!(eeg, channel=20:24)
eeg_load_electrodes!(eeg, file_name="../locs/standard-10-20-cap19-elmiko.ced")

s, h = eeg_snr(e10)
@test size(s) == (19, 1280)

s, _ = eeg_standardize(eeg)
@test size(s.eeg_signals) == (19, 309760, 1)

eeg1 = eeg_epoch_time(eeg, ts=-10.0)
eeg1.eeg_epoch_time[1, 1] == -10.0

e10 = eeg_epoch(eeg, epoch_len=10*256)
@test size(eeg_tenv(e10)[1]) == (19, 2560, 121)
@test size(eeg_tenv_mean(e10, dims=1)[1]) == (2560, 121)
@test size(eeg_tenv_median(e10, dims=1)[1]) == (2560, 121)
@test size(eeg_penv(e10)[1]) == (19, 513, 121)
@test size(eeg_penv_mean(e10, dims=1)[1]) == (513, 121)
@test size(eeg_penv_median(e10, dims=1)[1]) == (513, 121)
@test size(eeg_senv(e10)[1]) == (19, 37, 121)
@test size(eeg_senv_mean(e10, dims=1)[1]) == (37, 121)
@test size(eeg_senv_median(e10, dims=1)[1]) == (37, 121)
@test size(eeg_wdenoise(eeg, wt=wavelet(WT.haar)).eeg_signals) == (19, 309760, 1)
@test length(eeg_ispc(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)) == 6
@test length(eeg_itpc(e10, channel=1, t=12)) == 4
@test length(eeg_pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)) == 5
@test size(eeg_pli(e10)) == (19, 19, 121)
@test size(eeg_ispc(e10)) == (19, 19, 121)
@test length(eeg_ec(eeg, eeg, channel1=1, channel2=2, epoch1=1, epoch2=1)) == 2
@test length(eeg_ged(eeg, eeg)) == 3
@test size(eeg_frqinst(eeg)) == size(eeg.eeg_signals)
@test size(eeg_fftdenoise(eeg).eeg_signals) == (19, 309760, 1)
@test size(eeg_tkeo(eeg)) == (19, 309760, 1)
@test length(eeg_mwpsd(eeg, frq_lim=(0, 20), frq_n=21)) == 2

c, msc, f = eeg_fcoherence(eeg, eeg, channel1=1, channel2=2, epoch1=1, epoch2=1)
@test length(c) == 262145

eeg1 = eeg_reference_plap(eeg)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

f, p = eeg_vartest(eeg)
@test size(f) == (19, 19, 1)

eeg1 = eeg_add_note(eeg, note="test")
@test eeg_view_note(eeg1) == "test"
eeg_delete_note!(eeg1)
@test eeg_view_note(eeg1) == ""

eeg1 = eeg_epoch(eeg, epoch_len=2560)
new_channel = zeros(1, eeg_epoch_len(eeg1), eeg_epoch_n(eeg1))
eeg1 = eeg_replace_channel(eeg1, channel=1, signal=new_channel);
@test eeg1.eeg_signals[1, :, :] == zeros(eeg_epoch_len(eeg1), eeg_epoch_n(eeg1))
eeg2 = eeg_plinterpolate_channel(eeg1, channel=1, epoch=1)
@test eeg2.eeg_signals[1, :, 1] != zeros(eeg_epoch_len(eeg1))
eeg2 = eeg_lrinterpolate_channel(eeg1, channel=1, epoch=1)
@test eeg2.eeg_signals[1, :, 1] != zeros(eeg_epoch_len(eeg1))

@test length(eeg_band_mpower(eeg, f=(1,4))) == 3

p, f = eeg_rel_psd(eeg, f=(8,12))
@test size(p) == (19, 513, 1)

_, _, ss = eeg_fbsplit(eeg)
@test size(ss) == (10, 19, 309760, 1)

eeg1 = eeg_zero(eeg)
@test eeg1.eeg_signals[1, 1, 1] == 0

c = eeg_chdiff(eeg, eeg, channel1=1, channel2=2)
@test size(c) == (1, 309760, 1)

eeg1 = eeg_wbp(eeg, frq=10)
@test size(eeg1.eeg_signals) == (19, 309760, 1)
eeg1 = eeg_cbp(eeg, frq=10)
@test size(eeg1.eeg_signals) == (19, 309760, 1)
eeg1 = eeg_denoise_wien(eeg)
@test size(eeg1.eeg_signals) == (19, 309760, 1)

p, _, _ = eeg_cps(eeg, eeg, channel1=1, channel2=2, epoch1=1, epoch2=1)
@test length(p) == 262145

eeg2 = eeg_channel_type(eeg, channel=1, type="eog")
@test eeg2.eeg_header[:channel_type][1] == "eog"
eeg2 = eeg_edit_electrode(eeg, channel=1, x=2)
@test eeg2.eeg_locs[!, :loc_x][1] == 2.0
_, _, x, _, _, _, _, _ = eeg_electrode_loc(eeg2, channel=1, output=false)
@test x == 2.0

ch1 = eeg_electrode_loc(eeg, channel=1, output=false)
@test ch1[1] == 108.0

locs = locs_import_ced("../locs/standard-10-20-cap19-elmiko.ced")
locs2 = locs_flipx(locs)
@test locs2[1, 3] == 72.0
locs2 = locs_flipy(locs)
@test locs2[1, 3] == 252.0
locs2 = locs_flipz(locs)
@test locs2[1, 3] == 108.0
locs2 = locs_swapxy(locs)
@test locs2[1, 3] == 198.0
locs2 = locs_sph2cart(locs)
@test locs2[1, 5] == -0.309
locs2 = locs_cart2sph(locs)
@test locs2[1, 3] == 108.0
locs2 = locs_maximize(locs)
@test locs2[1, :loc_radius] == 1.0

@test size(eeg_phdiff(eeg)) == (19, 309760, 1)
@test size(eeg_scale(eeg, channel=1, factor=0.1).eeg_signals) == (19, 309760, 1)

_, _, f = eeg_psdslope(eeg)
@test length(f) == 513

@test size(eeg_vch(e10, f="fp1 + fp2")) == (1, 2560, 121)
@test size(eeg_dwt(e10, wt=wavelet(WT.haar), type=:sdwt)) == (19, 10, 2560, 121)
@test size(eeg_cwt(e10, wt=wavelet(Morlet(π), β=2))) == (19, 33, 2560, 121)

@test size(eeg_henv(e10)[1]) == (19, 2560, 121)
@test size(eeg_henv_mean(e10, dims=1)[1]) == (2560, 121)
@test size(eeg_henv_median(e10, dims=1)[1]) == (2560, 121)
@test size(eeg_apply(e10, f="mean(eeg, dims=1)")) == (19, 1, 121)

@test eeg_channel_cluster(e10, cluster=:f1) == [1, 3, 11]

e1 = eeg_copy(eeg)
eeg_add_marker!(e1, id="1", start=100, len=1, desc="test")
eeg_add_marker!(e1, id="1", start=1000, len=1, desc="test")
eeg_add_marker!(e1, id="1", start=2000, len=1, desc="test")
eeg_add_marker!(e1, id="1", start=3000, len=1, desc="test")
eeg_add_marker!(e1, id="1", start=4000, len=1, desc="test")
eeg_add_marker!(e1, id="1", start=5000, len=1, desc="test")
e2 = eeg_trim(e1, segment=(1, 400), remove_epochs=false)
@test e2.eeg_markers[1, :start] == 600
e2 = eeg_epoch(e1, epoch_len=200)
eeg_delete_epoch!(e2, epoch=1)
@test e2.eeg_markers[1, :start] == 800

eeg1, g, h = eeg_slaplacian(eeg)
@test size(eeg1.eeg_signals) == (24, 309760, 1)
@test size(g) == (19, 19)
@test size(h) == (19, 19)

b = eeg_bands_dwt(eeg, channel=1, wt=wavelet(WT.db2), type=:sdwt, n=5)
@test size(b) == (5, 309760)

true