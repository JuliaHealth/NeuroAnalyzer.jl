using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@info "Initializing"
eeg = import_edf("files/eeg-test-edf.edf")
e10 = epoch(eeg, ep_len=10*sr(eeg))
keep_epoch!(e10, ep=1:10)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "test 1/19: add_signal()"
@test add_signal(v1, v2) == v1 + v2
x = rand(epoch_len(e10))
e10_tmp = add_signal(e10, s=x)
e10_tmp.data[1, :, 1] == e10.data[1, :, 1] + x

@info "test 2/19: average()"
@test average(a1) == ones(1, 3, 2)
@test average(a1, a2) == 0.5 .* ones(2, 1, 2)
e10_tmp = average(e10)
@test size(e10_tmp.data) == (1, 2560, 10)
e10_tmp = average(e10, e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 3/19: cbp()"
@test length(cbp(rand(100), fs=10, frq=4)) == 100
e10_tmp = cbp(e10, frq=4)
@test size(e10_tmp.data) == size(e10.data)

@info "test 4/19: ch_zero()"
e10_tmp = ch_zero(e10)
@test e10_tmp.data[1, 1, 1] == 0
@test e10_tmp.data[1, end, 1] == 0

@info "test 5/19: cw_trans()"
s = rand(100)
ct = cw_trans(s, wt=wavelet(Morlet(π), β=2))
@test size(ct) == (14, 100)
ct = cw_trans(e10, wt=wavelet(Morlet(π), β=2))
@test size(ct) == (19, 33, 2560, 10)

@info "test 6/19: icw_trans()"
ct = cw_trans(s, wt=wavelet(Morlet(π), β=2))
s_new = icw_trans(ct, wt=wavelet(Morlet(π), β=2), type=:nd)
@test length(s_new) == 100
s_new = icw_trans(ct, wt=wavelet(Morlet(π), β=2), type=:pd)
@test length(s_new) == 100
s_new = icw_trans(ct, wt=wavelet(Morlet(π), β=2), type=:df)
@test length(s_new) == 100

@info "test 7/19: denoise_fft()"
s = rand(100)
s2, f = denoise_fft(s)
@test length(s2) == 100
@test length(f) == 100
e10_tmp = denoise_fft(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 8/19: denoise_wavelet())"
s = denoise_wavelet(rand(100), wt=wavelet(WT.haar))
@test length(s) == 100
e10_tmp = denoise_wavelet(e10, wt=wavelet(WT.haar))
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 9/19: denoise_wien()"
s = denoise_wien(a1)
@test size(s) == (2, 3, 2)
e10_tmp = denoise_wien(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 10/19: derivative()"
@test derivative(v1) == ones(5)
@test derivative(a1) == zeros(2, 3, 2)
e10_tmp = derivative(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 11/19: detrend()"
@test detrend(v1) == zeros(5)
@test round.(detrend(a1)) == zeros(2, 3, 2)
e10_tmp = detrend(e10, type=:ls)
@test size(e10_tmp.data) == (24, 2560, 10)
e10_tmp = detrend(e10, type=:linear)
@test size(e10_tmp.data) == (24, 2560, 10)
e10_tmp = detrend(e10, type=:constant)
@test size(e10_tmp.data) == (24, 2560, 10)
e10_tmp = detrend(e10, type=:poly)
@test size(e10_tmp.data) == (24, 2560, 10)
e10_tmp = detrend(e10, type=:loess)
@test size(e10_tmp.data) == (24, 2560, 10)
e10_tmp = detrend(e10, type=:hp)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 12/19: dw_trans()"
s = rand(100)
dt = dw_trans(s, wt=wavelet(WT.haar), type=:sdwt)
@test size(dt) == (3, 100)
dt = dw_trans(s, wt=wavelet(WT.haar), type=:acdwt)
@test size(dt) == (3, 100)
dt = dw_trans(e10, wt=wavelet(WT.haar), type=:sdwt)
@test size(dt) == (19, 10, 2560, 10)
dt = dw_trans(e10, wt=wavelet(WT.haar), type=:acdwt)
@test size(dt) == (19, 10, 2560, 10)

@info "test 13/19: idw_trans()"
dt = dw_trans(s, wt=wavelet(WT.haar), type=:sdwt)
s_new = idw_trans(dt, wt=wavelet(WT.haar), type=:sdwt)
@test length(s_new) == 100
s_new = idw_trans(dt, wt=wavelet(WT.haar), type=:acdwt)
@test length(s_new) == 100

@info "test 14/19: dwtsplit()"
s = dwtsplit(e10, ch=1, wt = wavelet(WT.haar), type=:sdwt)
@test size(s) == (10, 2560, 10)

@info "test 15/19: erp()"
e = erp(e10)
@test size(e.data) == (24, 2560, 1)
@test e.time_pts == e.epoch_time

@info "test 16/19: bpsplit()"
s, bn, bf = bpsplit(e10)
@test length(bn) == 10
@test length(bf) == 10
@test size(s) == (10, 19, 2560, 10)
Plots.plot(s[1, 1, :, 1])

@info "test 17/19: fconv()"
@test fconv(v1, kernel=v2) == [2.0, 2.0, 3.0, 3.0, 2.0]
@test fconv(a1, kernel=[0.5, 1.0, 0.5]) == ones(2, 3, 2)
e10_tmp = fconv(e10, kernel=[0.0, 0.5, 1.0, 0.5, 0.0])
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 18/19: filter_mavg()"
@test filter_mavg(v1, k=2) == [0.0, 0.0, 3.0, 0.0, 0.0]
e10_tmp = filter_mavg(e10, k=2)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 19/19: filter_mmed()"
@test filter_mmed(v1, k=2) == [0.0, 0.0, 3.0, 0.0, 0.0]
e10_tmp = filter_mmed(e10, k=2)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 20/19: filter_poly()"
s = filter_poly(rand(20))
@test length(s) == 20
e10_tmp = filter_poly(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 21/19: filter_sg()"
s = filter_sg(rand(20))
@test length(s) == 20
e10_tmp = filter_poly(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 22/19: filter()"
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:lp, cutoff=40, order=8)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:hp, cutoff=1, order=12)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:bs, cutoff=(49, 51), order=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:bp, cutoff=(49, 51), order=4)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev1, ftype=:lp, cutoff=40, order=8, rs=35)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev1, ftype=:hp, cutoff=1, order=12, rs=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev1, ftype=:bs, cutoff=(49, 51), order=2, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev1, ftype=:bp, cutoff=(49, 51), order=4, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev2, ftype=:lp, cutoff=40, order=8, rs=35)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev2, ftype=:hp, cutoff=1, order=12, rs=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev2, ftype=:bs, cutoff=(49, 51), order=2, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev2, ftype=:bp, cutoff=(49, 51), order=4, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:elliptic, ftype=:lp, cutoff=40, order=8, rs=35)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:elliptic, ftype=:hp, cutoff=1, order=12, rs=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:elliptic, ftype=:bs, cutoff=(49, 51), order=2, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:elliptic, ftype=:bp, cutoff=(49, 51), order=4, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:lp, cutoff=40, order=8)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:hp, cutoff=1, order=12)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:bs, cutoff=(49, 51), order=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:bp, cutoff=(49, 51), order=4)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:lp, cutoff=40, order=8, bw=10)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:hp, cutoff=1, order=12, bw=0.5)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:bs, cutoff=(49, 51), order=2, bw=0.5)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:bp, cutoff=(49, 51), order=4, bw=0.5)
@test size(eeg_tmp.data) == (24, 2560, 10)

@info "test 23/19: filter_g()"
s = filter_g(rand(20), fs=2, f=4)
@test length(s) == 20
e10_tmp = filter_g(e10, f=20)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 24/19: invert_polarity()"
e10_tmp = invert_polarity(e10)
@test e10_tmp.data == .-(e10.data)

@info "test 25/19: lrinterpolate_channel()"
e10_tmp = deepcopy(e10)
e10_tmp.data[1, :, 1] = zeros(epoch_len(e10))
e10_int = lrinterpolate_channel(e10_tmp, ch=1, ep=1)
@test e10_int.data[1, :, 1] != e10_tmp.data[1, :, 1]

@info "test 26/19: normalize()"
@show normalize(v1, method=:none) == v1
@show normalize_zscore(v1) == [-1.2649110640673518, -0.6324555320336759, 0.0, 0.6324555320336759, 1.2649110640673518]
@show normalize_minmax(v1) == [-1.0, -0.5, 0.0, 0.5, 1.0]
@show normalize_n(v1) == [0.0, 0.25, 0.5, 0.75, 1.0]
@show normalize_log(v1) == [1.0986122886681098, 1.3862943611198906, 1.6094379124341003, 1.791759469228055, 1.9459101490553132]
@show normalize_gauss(v1) == [-0.8047189562170504, -0.3465735902799727, 0.0, 0.3465735902799726, 0.8047189562170504]
@show normalize_log10(v1) == [0.47712125471966244, 0.6020599913279624, 0.6989700043360189, 0.7781512503836436, 0.8450980400142568]
@show normalize_neglog(v1) == [-0.0, -0.6931471805599453, -1.0986122886681098, -1.3862943611198906, -1.6094379124341003]
@show normalize_neglog10(v1) == [-0.0, -0.3010299956639812, -0.47712125471966244, -0.6020599913279624, -0.6989700043360189]
@show normalize_neg(v1) == [-4, -3, -2, -1, 0]
@show normalize_pos(v1) == [2, 3, 4, 5, 6]
@show normalize_perc(v1) == [0.0, 0.25, 0.5, 0.75, 1.0]
@show normalize_invroot(v1) == [0.7071067811865475, 0.5773502691896258, 0.5, 0.4472135954999579, 0.4082482904638631]
e10_tmp = normalize(e10, method=:zscore)
@test size(e10_tmp.data) == (24, 2560, 10)




eeg1 = reference_a(eeg)
@test size(eeg1.data) == (21, 309760, 1)
a1 = extract_channel(eeg, channel=20)
a2 = extract_channel(eeg, channel=21)
delete_channel!(eeg, channel=[20, 21])

eeg1 = derivative(eeg)
@test size(eeg1.data) == (19, 309760, 1)

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

@test length(NeuroAnalyzer.pick(eeg, p=:left)) == 8

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

@test size(denoise_wavelet(eeg, wt=wavelet(WT.haar)).data) == (19, 309760, 1)
@test size(fftdenoise(eeg).data) == (19, 309760, 1)

eeg1 = reference_plap(eeg)
@test size(eeg1.data) == (19, 309760, 1)

eeg1 = epoch(eeg, ep_len=2560)
new_channel = zeros(1, epoch_len(eeg1), epoch_n(eeg1))
eeg1 = replace_channel(eeg1, channel=1, signal=new_channel)
@test eeg1.data[1, :, :] == zeros(epoch_len(eeg1), epoch_n(eeg1))
eeg2 = plinterpolate_channel(eeg1, channel=1, epoch=1);
@test eeg2.data[1, :, 1] != zeros(epoch_len(eeg1))

eeg1 = epoch(eeg, ep_len=2560);
new_channel = zeros(1, epoch_len(eeg1), 1)
eeg1.data[1, :, 1] = zeros(epoch_len(eeg1))
eeg2 = lrinterpolate_channel(eeg1, channel=1, epoch=1);
@test eeg2.data[1, :, 1] != zeros(epoch_len(eeg1))

@test length(band_mpower(eeg, f=(1,4))) == 3

p, f = psd_rel(eeg, f=(8,12))
@test size(p) == (19, 513, 1)

_, _, ss = fbsplit(eeg)
@test size(ss) == (10, 19, 309760, 1)

eeg1 = ch_zero(eeg)
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
@test eeg2.header.recording[:channel_type][1] == "eog"

@test size(phdiff(eeg)) == (19, 309760, 1)
@test size(scale(eeg, channel=1, factor=0.1).data) == (19, 309760, 1)

_, _, f = psdslope(eeg)
@test length(f) == 513

@test size(vch(e10, f="fp1 + fp2")) == (1, 2560, 121)
@test size(dw_trans(e10, wt=wavelet(WT.haar), type=:sdwt)) == (19, 10, 2560, 121)
@test size(cw_trans(e10, wt=wavelet(Morlet(π), β=2))) == (19, 33, 2560, 121)

@test size(henv(e10)[1]) == (19, 2560, 121)
@test size(henv_mean(e10, dims=1)[1]) == (2560, 121)
@test size(henv_median(e10, dims=1)[1]) == (2560, 121)
@test size(apply(e10, f="mean(obj, dims=1)")) == (19, 1, 121)

eeg1, g, h = slaplacian(eeg)
@test size(eeg1.data) == (19, 309760, 1)
@test size(g) == (19, 19)
@test size(h) == (19, 19)

b = dwtsplit(eeg, channel=1, wt=wavelet(WT.db2), type=:sdwt, n=5)
@test size(b) == (5, 309760, 1)

r = reflect(eeg)
c = NeuroAnalyzer.chop(r)
@test size(eeg.data) == size(c.data)

@test size(extract_data(eeg, channel=1:channel_n(eeg))) == size(eeg.data)
@test length(extract_time(eeg)) == length(eeg.time_pts)
@test length(extract_etime(eeg)) == length(eeg.epoch_time)

true