using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

eeg = import_bdf("eeg-test-bdfplus.bdf")
delete_marker!(eeg, n=1)
@test size(eeg.markers) == (1, 5)
add_marker!(eeg, id="event", start=1, len=1, desc="test", channel=0)
@test size(eeg.markers) == (2, 5)
edit_marker!(eeg, n=2, id="event2", start=1, len=1, desc="test2", channel=0)

eeg = import_edf("eeg-test-edf.edf")
@test size(eeg.data) == (24, 309760, 1)

ecg = extract_channel(eeg, channel=24)
eog2 = extract_channel(eeg, channel=23)
eog1 = extract_channel(eeg, channel=22)
delete_channel!(eeg, channel=22:24)

eeg1 = reference_a(eeg)
@test size(eeg1.data) == (21, 309760, 1)
a1 = extract_channel(eeg, channel=20)
a2 = extract_channel(eeg, channel=21)
delete_channel!(eeg, channel=[20, 21])

eeg1 = delete_channel(eeg, channel=1)
@test eeg1.header.recording[:channel_n] == 18

eeg1 = keep_channel(eeg, channel=1)
@test eeg1.header.recording[:channel_n] == 1

eeg1 = derivative(eeg)
@test size(eeg1.data) == (19, 309760, 1)

tbp = total_power(eeg)
@test size(tbp) == (19, 1)

abp = band_power(eeg, f=(2, 4))
@test size(abp) == (19, 1)

eeg1 = detrend(eeg)
@test size(eeg1.data) == (19, 309760, 1)

eeg1 = reference_ch(eeg, channel=1)
@test size(eeg1.data) == (19, 309760, 1)

eeg1 = reference_car(eeg)
@test size(eeg1.data) == (19, 309760, 1)

eeg1 = extract_channel(eeg, channel="Cz")
@test size(eeg1) == (1, 309760, 1)

eeg1 = extract_channel(eeg, channel=18)
@test size(eeg1) == (1, 309760, 1)

@test get_channel(eeg, channel=1) == "Fp1"
@test get_channel(eeg, channel="Fp1") == 1

eeg1 = rename_channel(eeg, channel="Cz", name="CZ")
@test eeg1.header.recording[:labels][18] == "CZ"
eeg1 = rename_channel(eeg, channel=1, name="FP1")
@test eeg1.header.recording[:labels][1] == "FP1"

eeg1 = taper(eeg, t=eeg.data[1, :, 1])
@test size(eeg1.data) == (19, 309760, 1)

eeg1 = remove_dc(eeg)
@test size(eeg1.data) == (19, 309760, 1)

eeg1 = normalize(eeg, method=:zscore)
@test size(eeg1.data) == (19, 309760, 1)
eeg1 = normalize(eeg, method=:minmax)
@test size(eeg1.data) == (19, 309760, 1)
eeg1 = normalize(eeg, method=:log)
@test size(eeg1.data) == (19, 309760, 1)
eeg1 = normalize(eeg, method=:gauss)
@test size(eeg1.data) == (19, 309760, 1)

cov_m = covm(eeg)
@test size(cov_m) == (19, 19, 309760, 1)

cor_m = corm(eeg)
@test size(cor_m) == (19, 19, 309760, 1)

eeg1 = NeuroAnalyzer.upsample(eeg, new_sr=512)
@test size(eeg1.data) == (19, 619519, 1)

@test typeof(history(eeg)) == Vector{String}

@test labels(eeg)[1] == "Fp1"

@test sr(eeg) == 256

eeg1 = epoch(eeg, ep_len=1000)
erp!(eeg1)
@test size(eeg1.data) == (19, 1000, 1)

eeg10 = epoch(eeg, ep_n=10)
eeg1 = extract_epoch(eeg, epoch=1)
@test size(eeg1.data) == (19, 309760, 1)

f, s = dft(eeg)
@test size(f) == (19, 309760, 1)

m, _, _, _ = msci95(eeg)
@test size(m) == (1, 309760)

m, _, _, _ = msci95(eeg, eeg)
@test m == zeros(1, 309760)

s, ss, p = difference(eeg, eeg)
@test p == [1.0]

eeg1 = NeuroAnalyzer.filter(eeg, fprototype=:butterworth, ftype=:lp, cutoff=2, order=8)
@test size(eeg1.data) == (19, 309760, 1)
eeg1 = NeuroAnalyzer.filter(eeg, fprototype=:mavg, order=10)
@test size(eeg1.data) == (19, 309760, 1)
eeg1 = NeuroAnalyzer.filter(eeg, fprototype=:mmed, order=10)
@test size(eeg1.data) == (19, 309760, 1)

eeg1 = NeuroAnalyzer.downsample(eeg, new_sr=128)
@test size(eeg1.data) == (19, 154880, 1)

acov_m, _ = acov(eeg)
@test size(acov_m) == (19, 3, 1)
xcov_m, _ = xcov(eeg)
@test size(xcov_m) == (361, 3, 1)

p, f = psd(eeg)
@test size(p, 1) == 19

p = stationarity(eeg, method=:mean)
@test size(p) == (19, 10, 1)
p = stationarity(eeg, method=:var)
@test size(p) == (19, 10, 1)
p = stationarity(eeg, method=:hilbert)
@test size(p) == (19, 309759, 1)
p = stationarity(eeg, window=10000, method=:cov)
@test size(p) == (32, 1)

e = trim(eeg, segment=(10 * sr(eeg), 20 * sr(eeg)), remove_epochs=false)
@test signal_len(e) == 307199

m = mutual_information(eeg)
@test size(m) == (19, 19, 1)
m = mutual_information(eeg, eeg)
@test size(m) == (19, 19, 1)

e = entropy(eeg)
@test length(e) == 3
e = negentropy(eeg)
@test size(e) == (19, 1)

a = band_frq(eeg, band=:alpha)
@test a == (8, 13)

c, msc, ic = tcoherence(eeg, eeg)
@test size(c) == (19, 309760, 1)

hz, nyq = freqs(eeg)
@test nyq == 128.0
@test length(hz) == 154880

e10 = epoch(eeg, ep_len=2560)
s_conv = fconv(e10, kernel=generate_window(:hann, 256))
@test size(s_conv) == (19, 2560, 121)
s_conv = tconv(e10, kernel=generate_window(:hann, 256))
@test size(s_conv) == (19, 2560, 121)

p, v, m, pc_model = pca(eeg, n=2)
@test size(p) == (2, 309760, 1)
@test size(v) == (2, 1)
e1 = add_component(eeg, c=:pc, v=p)
add_component!(e1, c=:pc_model, v=pc_model)
e2 = pca_reconstruct(e1)
e2 = pca_reconstruct(eeg, p, pc_model)
@test size(e2.data) == (19, 309760, 1)

eeg.header.subject[:subject_last_name] = "unknown"
@test eeg.header.subject[:subject_last_name] == "unknown"

e = epoch(eeg, ep_n=10)
e9 = delete_epoch(e, epoch=10)
@test size(e9.data) == (19, 30976, 9)
e1 = keep_epoch(e, epoch=1)
@test size(e1.data) == (19, 30976, 1)

@test length(pick(eeg, p=:left)) == 8

e = epoch(eeg, ep_len=20*256)
v = epoch_stats(e)
@test length(v) == 10

e = epoch(eeg, ep_len=20)
erp!(e)
i, _ = ica(e, n=5, tol=1.0)
@test size(i) == (5, 20, 1)

e = deepcopy(eeg)
e_stats = epoch_stats(e)
@test length(e_stats) == (10)
add_component!(e, c=:epochs_mean, v=e_stats[1])
v = extract_component(e, c=:epochs_mean)
@test size(v) == (1, )
rename_component!(e, c_old=:epochs_mean, c_new=:epochs_m)
c = list_component(e)
@test size(c) == (1, )
c = component_type(e, c=:epochs_m)
@test c == Vector{Float64}
delete_component!(e, c=:epochs_m)
c = list_component(e)
@test size(c) == (0, )
reset_components!(e)
c = list_component(e)
@test size(c) == (0, )

e = epoch(eeg, ep_len=2560)
erp!(e)
p, f, t = spectrogram(e)
@test size(p) == (1281, 37, 19, 1)
p, f, t = spectrogram(e, method=:mt)
@test size(p) == (257, 15, 19, 1)
p, f, t = spectrogram(e, method=:mw)
@test size(p) == (129, 2560, 19, 1)
p, f, t = spectrogram(e, method=:stft)
@test size(p) == (1281, 37, 19, 1)
p, f, t = spectrogram(e, method=:gh)
@test size(p) == (129, 2560, 19, 1)
p, f, t = spectrogram(e, method=:cwt)
@test size(p) == (18, 2560, 19, 1)

f, a, p, ph = spectrum(e)
@test size(p) == (19, 1280, 1)

e = deepcopy(eeg)
i, iw = ica(e, tol=1.0, n=10)
add_component!(e, c=:ic, v=i)
add_component!(e, c=:ic_mw, v=iw)
@test size(e.components[1]) == (10, 309760, 1)
e2 = ica_reconstruct(e, ic_idx=1)
@test size(e2.data) == (19, 309760, 1)

b = detect_bad(eeg)
@test length(b) == 2

@test t2s(eeg, t=10) == 2561
@test s2t(eeg, t=10) == 0.04

e = keep_channel_type(eeg)
@test size(e.data) == (19, 309760, 1)
edit_channel!(e, channel=19, field=:channel_type, value="ecg")
keep_channel_type!(e, type=:eeg)
@test size(e.data) == (18, 309760, 1)

e = invert_polarity(eeg, channel=1)
@test e.data[1, 1, 1] == -eeg.data[1, 1, 1]

v = channel_stats(eeg)
@test length(v) == 10

eeg = import_edf("eeg-test-edf.edf")
delete_channel!(eeg, channel=20:24)
load_locs!(eeg, file_name="../locs/standard-10-20-cap19-elmiko.ced")
e10 = epoch(eeg, ep_len=10*256)

s, h = snr(e10)
@test size(s) == (19, 1280)

s, _ = standardize(eeg)
@test size(s.data) == (19, 309760, 1)

eeg1 = epoch_time(eeg, ts=-10.0)
eeg1.epoch_time[1, 1] == -10.0

@test size(tenv(e10)[1]) == (19, 2560, 121)
@test size(tenv_mean(e10, dims=1)[1]) == (2560, 121)
@test size(tenv_median(e10, dims=1)[1]) == (2560, 121)
@test size(penv(e10)[1]) == (19, 513, 121)
@test size(penv_mean(e10, dims=1)[1]) == (513, 121)
@test size(penv_median(e10, dims=1)[1]) == (513, 121)
@test size(senv(e10)[1]) == (19, 37, 121)
@test size(senv_mean(e10, dims=1)[1]) == (37, 121)
@test size(senv_median(e10, dims=1)[1]) == (37, 121)
@test size(denoise_wavelet(eeg, wt=wavelet(WT.haar)).data) == (19, 309760, 1)
@test length(ispc(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)) == 6
@test length(itpc(e10, channel=1, t=12)) == 4
@test length(pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)) == 5
@test size(pli(e10)) == (19, 19, 121)
@test size(ispc(e10)) == (19, 19, 121)
@test length(env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)) == 2
@test length(ged(e10, e10)) == 3
@test size(frqinst(eeg)) == (19, 2560, 121)
@test size(fftdenoise(eeg).data) == (19, 309760, 1)
@test size(tkeo(eeg)) == (19, 309760, 1)
@test length(psd_mw(eeg, frq_lim=(0, 20), frq_n=21)) == 2

c, msc, f = fcoherence(eeg, eeg, channel1=1, channel2=2, epoch1=1, epoch2=1)
@test length(c) == 262145

eeg1 = reference_plap(eeg)
@test size(eeg1.data) == (19, 309760, 1)

f, p = vartest(eeg)
@test size(f) == (19, 19, 1)

eeg1 = add_note(eeg, note="test")
@test view_note(eeg1) == "test"
delete_note!(eeg1)
@test view_note(eeg1) == ""

eeg1 = epoch(eeg, ep_len=2560)
new_channel = zeros(1, epoch_len(eeg1), epoch_n(eeg1))
eeg1 = replace_channel(eeg1, channel=1, signal=new_channel)
@test eeg1.data[1, :, :] == zeros(epoch_len(eeg1), epoch_n(eeg1))
eeg2 = plinterpolate_channel(eeg1, channel=1, epoch=1);
@test eeg2.data[1, :, 1] != zeros(ep_len(eeg1))

eeg1 = epoch(eeg, ep_len=2560);
new_channel = zeros(1, epoch_len(eeg1), 1)
eeg1.data[1, :, 1] = zeros(epoch_len(eeg1))
eeg2 = lrinterpolate_channel(eeg1, channel=1, epoch=1);
@test eeg2.data[1, :, 1] != zeros(epoch_len(eeg1))

@test length(band_mpower(eeg, f=(1,4))) == 3

p, f = rel_psd(eeg, f=(8,12))
@test size(p) == (19, 513, 1)

_, _, ss = fbsplit(eeg)
@test size(ss) == (10, 19, 309760, 1)

eeg1 = zero(eeg)
@test eeg1.data[1, 1, 1] == 0

c = chdiff(eeg, eeg, channel1=1, channel2=2)
@test size(c) == (1, 309760, 1)

eeg1 = wbp(eeg, frq=10)
@test size(eeg1.data) == (19, 309760, 1)
eeg1 = cbp(eeg, frq=10)
@test size(eeg1.data) == (19, 309760, 1)
eeg1 = denoise_wien(eeg)
@test size(eeg1.data) == (19, 309760, 1)

p, _, _ = cps(eeg, eeg, channel1=1, channel2=2, epoch1=1, epoch2=1)
@test length(p) == 262145

eeg2 = channel_type(eeg, channel=1, type="eog")
@test eeg2.header[:channel_type][1] == "eog"

@test size(phdiff(eeg)) == (19, 309760, 1)
@test size(scale(eeg, channel=1, factor=0.1).data) == (19, 309760, 1)

_, _, f = psdslope(eeg)
@test length(f) == 513

@test size(vch(e10, f="fp1 + fp2")) == (1, 2560, 121)
@test size(dwt(e10, wt=wavelet(WT.haar), type=:sdwt)) == (19, 10, 2560, 121)
@test size(cwt(e10, wt=wavelet(Morlet(π), β=2))) == (19, 33, 2560, 121)

@test size(henv(e10)[1]) == (19, 2560, 121)
@test size(henv_mean(e10, dims=1)[1]) == (2560, 121)
@test size(henv_median(e10, dims=1)[1]) == (2560, 121)
@test size(apply(e10, f="mean(eeg, dims=1)")) == (19, 1, 121)

@test channel_cluster(e10, cluster=:f1) == [1, 3, 11]

e1 = copy(eeg)
add_marker!(e1, id="1", start=100, len=1, desc="test")
add_marker!(e1, id="1", start=1000, len=1, desc="test")
add_marker!(e1, id="1", start=2000, len=1, desc="test")
add_marker!(e1, id="1", start=3000, len=1, desc="test")
add_marker!(e1, id="1", start=4000, len=1, desc="test")
add_marker!(e1, id="1", start=5000, len=1, desc="test")
e2 = trim(e1, segment=(1, 400), remove_epochs=false)
@test e2.markers[1, :start] == 600
e2 = epoch(e1, ep_len=200)
delete_epoch!(e2, epoch=1)
@test e2.markers[1, :start] == 800

eeg1, g, h = slaplacian(eeg)
@test size(eeg1.data) == (19, 309760, 1)
@test size(g) == (19, 19)
@test size(h) == (19, 19)

b = bands_dwt(eeg, channel=1, wt=wavelet(WT.db2), type=:sdwt, n=5)
@test size(b) == (5, 309760, 1)

r = reflect(eeg)
c = chop(r)
@test size(eeg.data) == size(c.data)

@test size(extract_data(eeg, channel=1:channel_n(eeg))) == size(eeg.data)
@test length(extract_time(eeg)) == length(eeg.time_pts)
@test length(extract_etime(eeg)) == length(eeg.epoch_time)

true