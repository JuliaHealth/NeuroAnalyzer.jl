using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@info "- Initializing"
eeg = import_edf("files/eeg-test-edf.edf")
e10 = epoch(eeg, ep_len=10*sr(eeg))
keep_epoch!(e10, epoch=1:10)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "test 1/ : acov()"
@test NeuroAnalyzer.acov(v) == (ac = [40.0, 55.0, 40.0], l = [-1, 0, 1])
@test NeuroAnalyzer.acov(a1) == (ac = [2.0 3.0 2.0; 2.0 3.0 2.0;;; 2.0 3.0 2.0; 2.0 3.0 2.0], l = [-1, 0, 1])
ac, l = NeuroAnalyzer.acov(e10)
@test size(ac) == (19, 3, 10)
@test l == [-3.90625, 0.0, 3.90625]

@info "test 2/ : ampdiff()"
@test size(NeuroAnalyzer.ampdiff(a1)) == (2, 3, 2)
ad = NeuroAnalyzer.ampdiff(e10)
@test size(ad) == (19, 2560, 10)

@info "test 3/ : band_power()"
@test NeuroAnalyzer.band_power(v, fs=10, f=(1, 2)) == 0
@test NeuroAnalyzer.band_power(a1, fs=10, f=(1, 2)) == zeros(2, 2)
@test size(NeuroAnalyzer.band_power(e10, f=(10, 20), mt=false)) == (19, 10)
@test size(NeuroAnalyzer.band_power(e10, f=(10, 20), mt=true)) == (19, 10)

@info "test 4/ : band_mpower()"
@test NeuroAnalyzer.band_mpower(v, fs=10, f=(1, 2)) == (mbp = 0.0, maxfrq = 1.0, maxbp = 0.0)
@test NeuroAnalyzer.band_mpower(a1, fs=10, f=(1, 2)) == (mbp = [0.0 0.0; 0.0 0.0], maxfrq = [1.0 1.0; 1.0 1.0], maxbp = [0.0 0.0; 0.0 0.0])
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, f=(10, 20), mt=false)
@test size(mbp) == (19, 10)
@test size(maxf) == (19, 10)
@test size(maxbp) == (19, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, f=(10, 20), mt=true)
@test size(mbp) == (19, 10)
@test size(maxf) == (19, 10)
@test size(maxbp) == (19, 10)

@info "test 5/ : corm()"
@test NeuroAnalyzer.corm(v) ≈ ones(5, 5)
@test size(NeuroAnalyzer.corm(a1)) == (2, 2, 3, 2)
@test size(NeuroAnalyzer.corm(e10)) == (19, 19, 2560, 10)

@info "test 6/ : covm()"
@test NeuroAnalyzer.covm(v) == [ 2.5  5.0  7.5 10.0 12.5;
                                 5.0 10.0 15.0 20.0 25.0;
                                 7.5 15.0 22.5 30.0 37.5;
                                10.0 20.0 30.0 40.0 50.0;
                                12.5 25.0 37.5 50.0 62.5]
@test size(NeuroAnalyzer.covm(a1)) == (2, 2, 3, 2)
@test size(NeuroAnalyzer.covm(e10)) == (19, 19, 2560, 10)

@info "test 7/ : cps()"
cp, cph, cf = NeuroAnalyzer.cps(rand(10), rand(10), fs=1)
@test length(cp) == 9
@test length(cph) == 9
@test length(cf) == 9
cp, cph, cf = NeuroAnalyzer.cps(rand(10, 10, 2), fs=1)
@test size(cp) == (10, 10, 9, 2) 
@test size(cph) == (10, 10, 9, 2) 
@test length(cf) == 9
cp, cph, cf = NeuroAnalyzer.cps(e10)
@test size(cp) == (19, 19, 2049, 10) 
@test size(cph) == (19, 19, 2049, 10) 
@test length(cf) == 2049
cp, cph, cf = NeuroAnalyzer.cps(rand(10, 10, 2), rand(10, 10, 2), fs=1)
@test size(cp) == (10, 9, 2) 
@test size(cph) == (10, 9, 2) 
@test length(cf) == 180
cp, cph, cf = NeuroAnalyzer.cps(e10, e10, ch1=1:2, ch2=2:3, ep1=1, ep2=1)
@test size(cp) == (2, 2049, 1) 
@test size(cph) == (2, 2049, 1) 
@test length(cf) == 4098

@info "test 8/ : diss()"
@test NeuroAnalyzer.diss(v1, v2) == (gd = 0.21320071635561044, sc = 0.9772727272727273)
@test NeuroAnalyzer.diss(a1) == (gd = [0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0], sc = [1.0 1.0; 1.0 1.0;;; 1.0 1.0; 1.0 1.0])
gd, sc = NeuroAnalyzer.diss(a1, a2)
@test size(gd) == (2, 2)
@test size(sc) == (2, 2)
gd, sc = NeuroAnalyzer.diss(e10)
@test size(gd) == (19, 19, 10)
@test size(sc) == (19, 19, 10)

@info "test 9/ : entropy()"
e, s, l = entropy(rand(10))
@test e < s < l
e, s, l = entropy(rand(10, 10))
@test size(e) == (10, 1)
@test size(s) == (10, 1)
@test size(l) == (10, 1)
e, s, l = NeuroAnalyzer.entropy(e10)
@test size(e) == (19, 10)
@test size(s) == (19, 10)
@test size(l) == (19, 10)

@info "test 10/ : negentropy()"
n = NeuroAnalyzer.negentropy(rand(10))
@test n < 0
n = NeuroAnalyzer.negentropy(rand(10, 10))
@test size(n) == (10, 1)
n = NeuroAnalyzer.negentropy(eeg)
@test size(n) == (19, 1)

@info "test 11/ : tenv()"
e, t = NeuroAnalyzer.tenv(e10)
@test size(e) == (19, 2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_mean(e10, dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_median(e10, dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560

@info "test 12/ : senv()"
e, t = NeuroAnalyzer.senv(e10)
@test size(e) == (19, 37, 10)
@test length(t) == 37
em, eu, el, t = NeuroAnalyzer.senv_mean(e10, dims=1)
@test size(em) == (37, 10)
@test size(eu) == (37, 10)
@test size(el) == (37, 10)
@test length(t) == 37
em, eu, el, t = NeuroAnalyzer.senv_median(e10, dims=1)
@test size(em) == (37, 10)
@test size(eu) == (37, 10)
@test size(el) == (37, 10)
@test length(t) == 37

@info "test 13/ : penv()"
e, t = NeuroAnalyzer.penv(e10)
@test size(e) == (19, 513, 10)
@test length(t) == 513
em, eu, el, t = NeuroAnalyzer.penv_mean(e10, dims=1)
@test size(em) == (513, 10)
@test size(eu) == (513, 10)
@test size(el) == (513, 10)
@test length(t) == 513
em, eu, el, t = NeuroAnalyzer.penv_median(e10, dims=1)
@test size(em) == (513, 10)
@test size(eu) == (513, 10)
@test size(el) == (513, 10)
@test length(t) == 513

@info "test 14/ : henv()"
e, t = NeuroAnalyzer.henv(e10)
@test size(e) == (19, 2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_mean(e10, dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_median(e10, dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560

@info "test 15/ : env_cor()"
ec, p = NeuroAnalyzer.env_cor(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1, type=:amp)
@test ec[1] <= 1.0
@test p[1] <= 1.0
ec, p = NeuroAnalyzer.env_cor(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1, type=:pow)
@test ec[1] <= 1.0
@test p[1] <= 1.0
ec, p = NeuroAnalyzer.env_cor(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1, type=:spec)
@test ec[1] <= 1.0
@test p[1] <= 1.0
ec, p = NeuroAnalyzer.env_cor(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1, type=:hamp)
@test ec[1] <= 1.0
@test p[1] <= 1.0

@info "test 15/ : erp_peaks()"
p = NeuroAnalyzer.erp_peaks(e10)
@test size(p) == (19, 2)

#=
NeuroAnalyzer.mdiff(e10, e10, method=:absdiff)
NeuroAnalyzer.mdiff(e10, e10, method=:diff2int)
NeuroAnalyzer.phdiff(e10) 
NeuroAnalyzer.e10_ica = add_component(e10, c=:ic, v=ic)
NeuroAnalyzer.channel_stats(e10)
NeuroAnalyzer.epoch_stats(e10)
NeuroAnalyzer.fbsplit(e10)
NeuroAnalyzer.fcoherence(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)
NeuroAnalyzer.fconv(e10, kernel=generate_morlet(256, 1, 32, complex=true))
NeuroAnalyzer.frqinst(e10)
NeuroAnalyzer.ged(e10, e10)
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
NeuroAnalyzer.pli(e10)
NeuroAnalyzer.pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)
NeuroAnalyzer.psd(e10)
NeuroAnalyzer.psd(e10, mt=true)
NeuroAnalyzer.psd_mw(e10)
NeuroAnalyzer.psd_rel(e10 f=(10, 20))
NeuroAnalyzer.psdslope(e10)
NeuroAnalyzer.remove_dc(e10)
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

####

eeg1 = wbp(eeg, frq=10)
@test size(eeg1.data) == (19, 309760, 1)
eeg1 = cbp(eeg, frq=10)
@test size(eeg1.data) == (19, 309760, 1)

@test size(phdiff(eeg)) == (19, 309760, 1)

=#

true