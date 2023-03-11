using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

eeg = import_edf("files/eeg-test-edf.edf")
e10 = epoch(eeg, ep_len=10*sr(eeg))
keep_epoch!(e10, epoch=1:10)
v = [1, 2, 3, 4, 5]
m = [1 2 3; 4 5 6]
a1 = ones(2, 3, 2)
a0 = zeros(2, 3, 2)

@test NeuroAnalyzer.acov(v) == (acov = [40.0, 55.0, 40.0], lags = [-1, 0, 1])
a, l = NeuroAnalyzer.acov(e10)
@test size(a) == 19, 3, 10
@test l == [-3.90625, 0.0, 3.90625]

NeuroAnalyzer.band_power(e10, f=(10, 20))
NeuroAnalyzer.band_power(e10, f=(10, 20), mt=true)
NeuroAnalyzer.band_mpower(e10, f=(10, 20))

NeuroAnalyzer.corm(e10)
NeuroAnalyzer.covm(e10)
NeuroAnalyzer.cps(e10)
NeuroAnalyzer.cps(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)
NeuroAnalyzer.cw_trans(e10, wt=wavelet(Morlet(π), β=2))
NeuroAnalyzer.ampdiff(e10)
NeuroAnalyzer.chdiff(e10, e10) 
NeuroAnalyzer.chdiff(e10, e10, channel1=1, channel2=2)
NeuroAnalyzer.phdiff(e10) 
NeuroAnalyzer.denoise_fft(e10)
NeuroAnalyzer.denoise_wien(e10)
NeuroAnalyzer.derivative(e10)
NeuroAnalyzer.detrend(e10, type=:constant)
NeuroAnalyzer.dft(e10)
NeuroAnalyzer.difference(e10, e10 method=:absdiff)
NeuroAnalyzer.difference(e10, e10 method=:diff2int)
NeuroAnalyzer.dw_trans(e10, wt=wavelet(WT.haar), type=:sdwt)
NeuroAnalyzer.dwtsplit(e10, channel=1, wt=wavelet(WT.db2), type=:sdwt, n=5)
NeuroAnalyzer.e10_ica = add_component(e10, c=:ic, v=ic)
NeuroAnalyzer.entropy(e10)
NeuroAnalyzer.env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:amp)
NeuroAnalyzer.env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:hamp)
NeuroAnalyzer.env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:pow)
NeuroAnalyzer.env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:spec)
NeuroAnalyzer.channel_stats(e10)
NeuroAnalyzer.epoch_stats(e10)
NeuroAnalyzer.fbsplit(e10)
NeuroAnalyzer.fcoherence(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)
NeuroAnalyzer.fconv(e10, kernel=generate_morlet(256, 1, 32, complex=true))
NeuroAnalyzer.frqinst(e10)
NeuroAnalyzer.ged(e10, e10)
NeuroAnalyzer.henv(e10)
NeuroAnalyzer.henv_mean(e10, dims=1)
NeuroAnalyzer.henv_median(e10, dims=1)
NeuroAnalyzer.ica(e10, n=15, tol=1.0)
NeuroAnalyzer.ica_reconstruct(e10_ica, ic_idx=1)
NeuroAnalyzer.ispc(e10)
NeuroAnalyzer.ispc(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)
NeuroAnalyzer.itpc(e10, channel=1, t=256)
NeuroAnalyzer.itpc_s(e10, channel=1, frq_lim=(10, 20), frq_n=11)
NeuroAnalyzer.msci95(e10)
NeuroAnalyzer.msci95(e10, e10)
NeuroAnalyzer.msci95(e10, method=:boot)
NeuroAnalyzer.mutual_information(e10)
NeuroAnalyzer.mutual_information(e10, e10)
NeuroAnalyzer.negentropy(e10)
NeuroAnalyzer.normalize(e10, method=:zscore)
NeuroAnalyzer.pca(e10, n=4)
NeuroAnalyzer.penv(e10)
NeuroAnalyzer.penv_mean(e10, dims=1)
NeuroAnalyzer.penv_median(e10, dims=1)
NeuroAnalyzer.pli(e10)
NeuroAnalyzer.pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)
NeuroAnalyzer.psd(e10)
NeuroAnalyzer.psd(e10, mt=true)
NeuroAnalyzer.psd_mw(e10)
NeuroAnalyzer.psd_rel(e10 f=(10, 20))
NeuroAnalyzer.psdslope(e10)
NeuroAnalyzer.remove_dc(e10)
NeuroAnalyzer.senv(e10)
NeuroAnalyzer.senv_mean(e10, dims=1)
NeuroAnalyzer.senv_median(e10, dims=1)
NeuroAnalyzer.snr(e10)
NeuroAnalyzer.spectrogram(e10)
NeuroAnalyzer.spectrogram(e10, method=:cwt)
NeuroAnalyzer.spectrogram(e10, method=:mt)
NeuroAnalyzer.spectrogram(e10, method=:mw)
NeuroAnalyzer.spectrogram(e10, method=:stft)
NeuroAnalyzer.spectrum(e10)
NeuroAnalyzer.spectrum(e10, h=true)
NeuroAnalyzer.standardize(e10)
NeuroAnalyzer.stationarity(e10, method=:adf)
NeuroAnalyzer.stationarity(e10, method=:cov)
NeuroAnalyzer.stationarity(e10, method=:hilbert)
NeuroAnalyzer.stationarity(e10, method=:mean)
NeuroAnalyzer.stationarity(e10, method=:var)
NeuroAnalyzer.taper(e10, t=generate_window(:hann, epoch_len(e10)))
NeuroAnalyzer.tcoherence(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=2)
NeuroAnalyzer.tconv(e10, kernel=generate_morlet(256, 1, 32, complex=true))
NeuroAnalyzer.tenv(e10)
NeuroAnalyzer.tenv_mean(e10, dims=1)
NeuroAnalyzer.tenv_median(e10, dims=1)
NeuroAnalyzer.tkeo(e10)
NeuroAnalyzer.total_power(e10)
NeuroAnalyzer.total_power(e10, mt=true)
NeuroAnalyzer.vartest(e10)
NeuroAnalyzer.vartest(e10, e10)
NeuroAnalyzer.xcov(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=2, lag=10, demean=true)
NeuroAnalyzer.xcov(e10, lag=10, demean=true)


acov_m, _ = acov(eeg)
@test size(acov_m) == (19, 3, 1)

tbp = total_power(eeg)
@test size(tbp) == (19, 1)

abp = band_power(eeg, f=(2, 4))
@test size(abp) == (19, 1)

cov_m = covm(eeg)
@test size(cov_m) == (19, 19, 309760, 1)

cor_m = corm(eeg)
@test size(cor_m) == (19, 19, 309760, 1)

m, _, _, _ = msci95(eeg)
@test size(m) == (1, 309760)

m, _, _, _ = msci95(eeg, eeg)
@test m == zeros(1, 309760)

s, ss, p = difference(eeg, eeg)
@test p == [1.0]

f, s = dft(eeg)
@test size(f) == (19, 309760, 1)

xcov_m, _ = xcov(eeg)
@test size(xcov_m) == (361, 3, 1)

p = stationarity(eeg, method=:mean)
@test size(p) == (19, 10, 1)
p = stationarity(eeg, method=:var)
@test size(p) == (19, 10, 1)
p = stationarity(eeg, method=:hilbert)
@test size(p) == (19, 309759, 1)
p = stationarity(eeg, window=10000, method=:cov)
@test size(p) == (32, 1)

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

e = epoch(eeg, ep_len=20*256)
v = epoch_stats(e)
@test length(v) == 10

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

v = channel_stats(eeg)
@test length(v) == 10

s, h = snr(e10)
@test size(s) == (19, 1280)

@test size(tenv(e10)[1]) == (19, 2560, 121)
@test size(tenv_mean(e10, dims=1)[1]) == (2560, 121)
@test size(tenv_median(e10, dims=1)[1]) == (2560, 121)
@test size(penv(e10)[1]) == (19, 513, 121)
@test size(penv_mean(e10, dims=1)[1]) == (513, 121)
@test size(penv_median(e10, dims=1)[1]) == (513, 121)
@test size(senv(e10)[1]) == (19, 37, 121)
@test size(senv_mean(e10, dims=1)[1]) == (37, 121)
@test size(senv_median(e10, dims=1)[1]) == (37, 121)
@test length(env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)) == 2

@test length(ispc(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)) == 6
@test length(itpc(e10, channel=1, t=12)) == 4
@test length(pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)) == 5
@test size(pli(e10)) == (19, 19, 121)
@test size(ispc(e10)) == (19, 19, 121)
@test length(ged(e10, e10)) == 3
@test size(frqinst(eeg)) == (19, 309760, 1)
@test size(denoise_fft(eeg).data) == (19, 309760, 1)
@test size(tkeo(eeg)) == (19, 309760, 1)
@test length(psd_mw(eeg, frq_lim=(0, 20), frq_n=21)) == 2

c, msc, f = fcoherence(eeg, eeg, channel1=1, channel2=2, epoch1=1, epoch2=1)
@test length(c) == 262145

f, p = vartest(eeg)
@test size(f) == (19, 19, 1)

@test length(band_mpower(eeg, f=(1,4))) == 3

p, f = psd_rel(eeg, f=(8,12))
@test size(p) == (19, 513, 1)

c = chdiff(eeg, eeg, channel1=1, channel2=2)
@test size(c) == (1, 309760, 1)

p, _, _ = cps(eeg, eeg, channel1=1, channel2=2, epoch1=1, epoch2=1)
@test length(p) == 262145

_, _, f = psdslope(eeg)
@test length(f) == 513

@test size(henv(e10)[1]) == (19, 2560, 121)
@test size(henv_mean(e10, dims=1)[1]) == (2560, 121)
@test size(henv_median(e10, dims=1)[1]) == (2560, 121)

####

eeg1 = wbp(eeg, frq=10)
@test size(eeg1.data) == (19, 309760, 1)
eeg1 = cbp(eeg, frq=10)
@test size(eeg1.data) == (19, 309760, 1)

@test size(phdiff(eeg)) == (19, 309760, 1)

true