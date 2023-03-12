using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@info "- Initializing"
eeg = import_edf("test/files/eeg-test-edf.edf")
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

@info "test 16/ : fcoherence()"
c, msc, f = fcoherence(rand(10, 100), fs=10)
@test size(c) == (10, 10, 65)
@test size(msc) == (10, 10, 65)
@test length(f) == 65
c, msc, f = fcoherence(rand(10, 100), rand(10, 100), fs=10)
@test length(c) == 9
@test length(msc) == 9
@test length(f) == 9
c, msc, f = fcoherence(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1)
@test length(c) == 2049
@test length(msc) == 2049
@test length(f) == 2049

@info "test 17/ : frqinst()"
f = frqinst(rand(100), fs=10)
@test length(f) == 100
f = frqinst(rand(10, 100, 10), fs=10)
@test size(f) == (10, 100, 10)

@info "test 18/ : ged()"
s, r, rn = ged(rand(10, 10), rand(10, 10))
@test length(s) == 100
@test length(r) == 10
@test length(rn) == 10
s, r, rn = ged(e10, e10)
@test length(s) == 486400
@test length(r) == 190
@test length(rn) == 190

@info "test 19/ : ica()"
ic, ic_mw = NeuroAnalyzer.ica(rand(10, 1000), n=5, tol=1.0)
@test size(ic) == (5, 1000)
@test size(ic_mw) == (10, 5)
ic, ic_mw = NeuroAnalyzer.ica(rand(10, 1000, 10), n=5, tol=1.0)
@test size(ic) == (5, 1000, 10)
@test size(ic_mw) == (10, 5, 10)
ic, ic_mw = NeuroAnalyzer.ica(e10, n=5, tol=1.0)
@test size(ic) == (5, 1000, 10)
@test size(ic_mw) == (10, 5, 10)

@info "test 20/ : ica_reconstruct()"
ic, ic_mw = NeuroAnalyzer.ica(rand(10, 1000), n=5, tol=1.0)
s = NeuroAnalyzer.ica_reconstruct(rand(10, 1000), ic=ic, ic_mw=ic_mw, ic_idx=5)
@test size(s) == (10, 1000)
ic, ic_mw = NeuroAnalyzer.ica(rand(10, 1000, 10), n=5, tol=1.0)
s = NeuroAnalyzer.ica_reconstruct(rand(10, 1000, 10), ic=ic, ic_mw=ic_mw, ic_idx=5)
@test size(s) == (10, 1000, 10)
ic, ic_mw = NeuroAnalyzer.ica(e10, n=5, tol=1.0)
e10_tmp = NeuroAnalyzer.ica_reconstruct(e10, ic, ic_mw; ic_idx=1)
@test size(e10_tmp.data) == (24, 2560, 10)
add_component!(e10, c=:ic, v=ic)
add_component!(e10, c=:ic_mw, v=ic_mw)
e10_tmp = NeuroAnalyzer.ica_reconstruct(e10, ic_idx=1)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 21/ : ispc()"
iv, ia, sd, pd, s1p, s2p = NeuroAnalyzer.ispc(v1, v2)
@test iv == 0.6125992852305387
@test ia == -0.0017801930770334254
@test sd == [5, 3, 1, -1, -3]
@test pd == [-1.3157044982273682, 0.8713795327960081, 0.3743702916488456, 0.7615478999167377, -1.0328436072470222]
@test s1p == [1.039406675134543, -0.6027563879589182, -0.21331750626000984, -0.33140501474755424, 0.3279718365439915]
@test s2p == [-0.27629782309282525, 0.26862314483709, 0.16105278538883577, 0.4301428851691834, -0.7048717707030308]
iv, ia = NeuroAnalyzer.ispc(e10)
@test size(iv) == (19, 19, 10)
@test size(ia) == (19, 19, 10)
iv, ia, sd, pd, s1p, s2p = NeuroAnalyzer.ispc(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1)
@test iv == [0.934947546068702;;]
@test ia == [0.0013110745656857103;;]
@test size(sd) == (1, 2560, 1)
@test size(pd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "test 22/ : itpc()"
iv, izv, ia, ip = NeuroAnalyzer.itpc(ones(1, 10, 10), t=1)
@test iv == 1.0
@test izv == 10.0
@test ia == 0.0
@test ip == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
iv, izv, ia, ip = NeuroAnalyzer.itpc(e10, ch=1, t=256)
@test iv == [0.9742650309891113]
@test izv == [9.491923506082141]
@test ia == [-0.09471374965342466]
@test ip[1] == -0.34043930845454473

@info "test 23/ : itpc_spec()"
iv, izv, f = NeuroAnalyzer.itpc_spec(e10, ch=1, frq_lim=(0, 4), frq_n=5)
@test size(iv) == (5, 2560)
@test size(izv) == (5, 2560)
@test f == [0.01, 0.044721359549995794, 0.20000000000000004, 0.8944271909999159, 4.0]

@info "test 24/ : mdiff()"
st, sts, p = NeuroAnalyzer.mdiff(m1, m2, method=:absdiff)
@test length(st) == 6
@test sts == 3.0
@test p in [0.0, 1.0]
st, sts, p = NeuroAnalyzer.mdiff(a1, a2, method=:absdiff)
@test size(st) == (2, 6)
@test sts == [1.0, 1.0]
@test p == [0.0, 0.0]
@test NeuroAnalyzer.mdiff(m1, m2, method=:diff2int) == (st = [2.6666666666666665, 8.666666666666666, 4.666666666666666, 1.1666666666666665, 15.166666666666666, 4.5], sts = 4.666666666666666, p = 1.0)
st, sts, p = NeuroAnalyzer.mdiff(a1, a2, method=:diff2int)
@test size(st) == (2, 6)
@test sts == [2.0, 2.0]
@test p == [0.0, 0.0]
st, sts, p = NeuroAnalyzer.mdiff(e10, e10, method=:absdiff)
@test size(st) == (10, 57)
@test sts == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test p == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
st, sts, p = NeuroAnalyzer.mdiff(e10, e10, method=:diff2int)
@test size(st) == (10, 57)
@test sts == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test p == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

@info "test 25/ : mutual_information()"
@test NeuroAnalyzer.mutual_information(v1, v2) == 0.4199730940219748
@test NeuroAnalyzer.mutual_information(a1) == [0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0]
@test NeuroAnalyzer.mutual_information(a1, a2) == [0.0 0.0; 0.0 0.0] 
m = NeuroAnalyzer.mutual_information(e10)
@test size(m) == (19, 19, 10) 
m = NeuroAnalyzer.mutual_information(e10, e10, ch1=1, ch2=2)
@test size(m) == (1, 1, 10)

#=
NeuroAnalyzer.phdiff(e10) 
NeuroAnalyzer.e10_ica = add_component(e10, c=:ic, v=ic)
NeuroAnalyzer.channel_stats(e10)
NeuroAnalyzer.epoch_stats(e10)
NeuroAnalyzer.fbsplit(e10)
NeuroAnalyzer.fconv(e10, kernel=generate_morlet(256, 1, 32, complex=true))
NeuroAnalyzer.ged(e10, e10)
NeuroAnalyzer.itpc_s(e10, channel=1, frq_lim=(10, 20), frq_n=11)
NeuroAnalyzer.msci95(e10)
NeuroAnalyzer.msci95(e10, e10)
NeuroAnalyzer.msci95(e10, method=:boot)
NeuroAnalyzer.negentropy(e10)
NeuroAnalyzer.normalize(e10, method=:zscore)
NeuroAnalyzer.pca(e10, n=4)
NeuroAnalyzer.pli(e10)
NeuroAnalyzer.pli(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1)
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
NeuroAnalyzer.tcoherence(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2)
NeuroAnalyzer.tconv(e10, kernel=generate_morlet(256, 1, 32, complex=true))
NeuroAnalyzer.tkeo(e10)
NeuroAnalyzer.total_power(e10)
NeuroAnalyzer.total_power(e10, mt=true)
NeuroAnalyzer.vartest(e10)
NeuroAnalyzer.vartest(e10, e10)
NeuroAnalyzer.xcov(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2, lag=10, demean=true)
NeuroAnalyzer.xcov(e10, lag=10, demean=true)


tbp = total_power(eeg)
@test size(tbp) == (19, 1)

m, _, _, _ = msci95(eeg)
@test size(m) == (1, 309760)

m, _, _, _ = msci95(eeg, eeg)
@test m == zeros(1, 309760)

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

v = channel_stats(eeg)
@test length(v) == 10

s, h = snr(e10)
@test size(s) == (19, 1280)

@test length(env_cor(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1)) == 2

@test length(pli(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1)) == 5
@test size(pli(e10)) == (19, 19, 121)
@test size(ispc(e10)) == (19, 19, 121)
@test length(ged(e10, e10)) == 3
@test size(denoise_fft(eeg).data) == (19, 309760, 1)
@test size(tkeo(eeg)) == (19, 309760, 1)
@test length(psd_mw(eeg, frq_lim=(0, 20), frq_n=21)) == 2

f, p = vartest(eeg)
@test size(f) == (19, 19, 1)

p, f = psd_rel(eeg, f=(8,12))
@test size(p) == (19, 513, 1)

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