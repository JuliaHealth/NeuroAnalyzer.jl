using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

ntests = 61

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
NeuroAnalyzer.filter!(e10, fprototype=:fir, ftype=:lp, cutoff=40, order=8)
NeuroAnalyzer.filter!(e10, fprototype=:fir, ftype=:hp, cutoff=1, order=8)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "Test 1/$ntests: acov()"
@test NeuroAnalyzer.acov(v) == [-0.8 -0.8 -0.2 0.8 2.0 0.8 -0.2 -0.8 -0.8;;;]
ac, l = NeuroAnalyzer.acov(e10)
@test size(ac) == (23, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acov(e10, biased=false)
@test size(ac) == (23, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acov(e10, method=:cov)
@test size(ac) == (23, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acov(e10, method=:cov, biased=false)
@test size(ac) == (23, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acov(e10, method=:stat)
@test size(ac) == (23, 3, 10)
@test length(l) == 3

@info "Test 2/$ntests: ampdiff()"
@test size(NeuroAnalyzer.ampdiff(a1)) == (2, 3, 2)
ad = NeuroAnalyzer.ampdiff(e10)
@test size(ad) == (23, 2560, 10)

@info "Test 3/$ntests: band_power()"
@test size(NeuroAnalyzer.band_power(e10, frq_lim=(10, 20))) == (23, 10)
@test size(NeuroAnalyzer.band_power(e10, frq_lim=(10, 20), method=:welch)) == (23, 10)
@test size(NeuroAnalyzer.band_power(e10, frq_lim=(10, 20), method=:fft)) == (23, 10)
@test size(NeuroAnalyzer.band_power(e10, frq_lim=(10, 20), method=:stft)) == (23, 10)
@test size(NeuroAnalyzer.band_power(e10, frq_lim=(10, 20), method=:mt)) == (23, 10)
@test size(NeuroAnalyzer.band_power(e10, frq_lim=(10, 20), method=:mw)) == (23, 10)
@test size(NeuroAnalyzer.band_power(e10, frq_lim=(10, 20), method=:gh)) == (23, 10)
@test size(NeuroAnalyzer.band_power(e10, frq_lim=(10, 20), method=:cwt)) == (23, 10)

@info "Test 4/$ntests: band_mpower()"
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, frq_lim=(10, 20))
@test size(mbp) == (23, 10)
@test size(maxf) == (23, 10)
@test size(maxbp) == (23, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, frq_lim=(10, 20), method=:mt)
@test size(mbp) == (23, 10)
@test size(maxf) == (23, 10)
@test size(maxbp) == (23, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, frq_lim=(10, 20), method=:stft)
@test size(mbp) == (23, 10)
@test size(maxf) == (23, 10)
@test size(maxbp) == (23, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, frq_lim=(10, 20), method=:fft)
@test size(mbp) == (23, 10)
@test size(maxf) == (23, 10)
@test size(maxbp) == (23, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, frq_lim=(10, 20), method=:mw)
@test size(mbp) == (23, 10)
@test size(maxf) == (23, 10)
@test size(maxbp) == (23, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, frq_lim=(10, 20), method=:gh)
@test size(mbp) == (23, 10)
@test size(maxf) == (23, 10)
@test size(maxbp) == (23, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, frq_lim=(10, 20), method=:cwt)
@test size(mbp) == (23, 10)
@test size(maxf) == (23, 10)
@test size(maxbp) == (23, 10)

@info "Test 5/$ntests: corm()"
@test NeuroAnalyzer.corm(v) ≈ ones(5, 5)
@test size(NeuroAnalyzer.corm(a1)) == (2, 2, 3, 2)
@test size(NeuroAnalyzer.corm(e10)) == (23, 23, 2560, 10)

@info "Test 6/$ntests: covm()"
@test NeuroAnalyzer.covm(v) == [ 2.5  5.0  7.5 10.0 12.5;
                                 5.0 10.0 15.0 20.0 25.0;
                                 7.5 15.0 22.5 30.0 37.5;
                                10.0 20.0 30.0 40.0 50.0;
                                12.5 25.0 37.5 50.0 62.5]
@test size(NeuroAnalyzer.covm(a1)) == (2, 2, 3, 2)
@test size(NeuroAnalyzer.covm(e10)) == (23, 23, 2560, 10)

@info "Test 7/$ntests: cph()"
ph, f = NeuroAnalyzer.cph(rand(10), rand(10), fs=1)
@test length(ph) == 9
@test length(f) == 9
ph, f = NeuroAnalyzer.cph(rand(10, 10, 2), fs=1)
@test size(ph) == (10, 10, 9, 2)
@test length(f) == 9
ph, f = NeuroAnalyzer.cph(e10)
@test size(ph) == (23, 23, 2049, 10)
@test length(f) == 2049
ph, f = NeuroAnalyzer.cph(rand(10, 10, 2), rand(10, 10, 2), fs=1)
@test size(ph) == (10, 9, 2)
@test length(f) == 9
ph, f = NeuroAnalyzer.cph(e10, e10, ch1=1:2, ch2=2:3, ep1=1, ep2=1)
@test size(ph) == (2, 2049, 1)
@test length(f) == 2049

@info "Test 8/$ntests: diss()"
@test NeuroAnalyzer.diss(v1, v2) == (gd = 0.21320071635561044, sc = 0.9772727272727273)
@test NeuroAnalyzer.diss(a1) == (gd = [0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0], sc = [1.0 1.0; 1.0 1.0;;; 1.0 1.0; 1.0 1.0])
gd, sc = NeuroAnalyzer.diss(a1, a2)
@test size(gd) == (2, 2)
@test size(sc) == (2, 2)
gd, sc = NeuroAnalyzer.diss(e10)
@test size(gd) == (23, 23, 10)
@test size(sc) == (23, 23, 10)

@info "Test 9/$ntests: entropy()"
e, s, l = NeuroAnalyzer.entropy(rand(10))
@test e < l
@test s < l
e, s, l = NeuroAnalyzer.entropy(rand(10, 10))
@test size(e) == (10, 1)
@test size(s) == (10, 1)
@test size(l) == (10, 1)
e, s, l = NeuroAnalyzer.entropy(e10)
@test size(e) == (23, 10)
@test size(s) == (23, 10)
@test size(l) == (23, 10)

@info "Test 10/$ntests: negentropy()"
n = NeuroAnalyzer.negentropy(rand(10))
@test n < 0
n = NeuroAnalyzer.negentropy(rand(10, 10))
@test size(n) == (10, 1)
n = NeuroAnalyzer.negentropy(eeg)
@test size(n) == (23, 1)

@info "Test 11/$ntests: tenv()"
e, t = NeuroAnalyzer.tenv(e10)
@test size(e) == (23, 2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_mean(e10, dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_mean(e10, dims=2)
@test size(em) == (2560, 23)
@test size(eu) == (2560, 23)
@test size(el) == (2560, 23)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_mean(e10, dims=3)
@test size(em) == (2560,)
@test size(eu) == (2560,)
@test size(el) == (2560,)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_median(e10, dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_median(e10, dims=2)
@test size(em) == (2560, 23)
@test size(eu) == (2560, 23)
@test size(el) == (2560, 23)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_median(e10, dims=3)
@test size(em) == (2560,)
@test size(eu) == (2560,)
@test size(el) == (2560,)
@test length(t) == 2560

@info "Test 12/$ntests: senv()"
e, t = NeuroAnalyzer.senv(e10)
@test size(e) == (23, 289, 10)
@test length(t) == 289
e, t = NeuroAnalyzer.senv(e10, method=:mt)
@test size(e) == (23, 15, 10)
@test length(t) == 15
e, t = NeuroAnalyzer.senv(e10, method=:mw)
@test size(e) == (23, 2560, 10)
@test length(t) == 2560
e, t = NeuroAnalyzer.senv(e10, method=:gh)
@test size(e) == (23, 2560, 10)
@test length(t) == 2560
e, t = NeuroAnalyzer.senv(e10, method=:cwt)
@test size(e) == (23, 2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.senv_mean(e10, dims=1)
@test size(em) == (289, 10)
@test size(eu) == (289, 10)
@test size(el) == (289, 10)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_mean(e10, dims=2)
@test size(em) == (289, 23)
@test size(eu) == (289, 23)
@test size(el) == (289, 23)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_mean(e10, dims=3)
@test size(em) == (289,)
@test size(eu) == (289,)
@test size(el) == (289,)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_median(e10, dims=1)
@test size(em) == (289, 10)
@test size(eu) == (289, 10)
@test size(el) == (289, 10)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_median(e10, dims=2)
@test size(em) == (289, 23)
@test size(eu) == (289, 23)
@test size(el) == (289, 23)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_median(e10, dims=3)
@test size(em) == (289,)
@test size(eu) == (289,)
@test size(el) == (289,)
@test length(t) == 289

@info "Test 13/$ntests: penv()"
e, t = NeuroAnalyzer.penv(e10, method=:welch)
@test size(e) == (23, 129, 10)
@test length(t) == 129
e, t = NeuroAnalyzer.penv(e10, method=:fft)
@test size(e) == (23, 1281, 10)
@test length(t) == 1281
e, t = NeuroAnalyzer.penv(e10, method=:stft)
@test size(e) == (23, 129, 10)
@test length(t) == 129
e, t = NeuroAnalyzer.penv(e10, method=:mt)
@test size(e) == (23, 1281, 10)
@test length(t) == 1281
e, t = NeuroAnalyzer.penv(e10, method=:mw)
@test size(e) == (23, 129, 10)
@test length(t) == 129
e, t = NeuroAnalyzer.penv(e10, method=:gh)
@test size(e) == (23, 129, 10)
@test length(t) == 129
e, t = NeuroAnalyzer.penv(e10, method=:cwt)
@test size(e) == (23, 131, 10)
@test length(t) == 131
em, eu, el, t = NeuroAnalyzer.penv_mean(e10, dims=1)
@test size(em) == (129, 10)
@test size(eu) == (129, 10)
@test size(el) == (129, 10)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_mean(e10, dims=2)
@test size(em) == (129, 23)
@test size(eu) == (129, 23)
@test size(el) == (129, 23)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_mean(e10, dims=3)
@test size(em) == (129,)
@test size(eu) == (129,)
@test size(el) == (129,)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_median(e10, dims=1)
@test size(em) == (129, 10)
@test size(eu) == (129, 10)
@test size(el) == (129, 10)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_median(e10, dims=2)
@test size(em) == (129, 23)
@test size(eu) == (129, 23)
@test size(el) == (129, 23)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_median(e10, dims=3)
@test size(em) == (129,)
@test size(eu) == (129,)
@test size(el) == (129,)
@test length(t) == 129

@info "Test 14/$ntests: henv()"
e, t = NeuroAnalyzer.henv(e10)
@test size(e) == (23, 2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_mean(e10, dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_mean(e10, dims=2)
@test size(em) == (2560, 23)
@test size(eu) == (2560, 23)
@test size(el) == (2560, 23)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_mean(e10, dims=3)
@test size(em) == (2560,)
@test size(eu) == (2560,)
@test size(el) == (2560,)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_median(e10, dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_median(e10, dims=2)
@test size(em) == (2560, 23)
@test size(eu) == (2560, 23)
@test size(el) == (2560, 23)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_median(e10, dims=3)
@test size(em) == (2560,)
@test size(eu) == (2560,)
@test size(el) == (2560,)
@test length(t) == 2560

@info "Test 15/$ntests: erp_peaks()"
e = NeuroAnalyzer.erp(e10)
p = NeuroAnalyzer.erp_peaks(e)
@test size(p) == (23, 2)

@info "Test 16/$ntests: coherence()"
c, msc, f = NeuroAnalyzer.coherence(rand(100), rand(100), fs=10, method=:mt)
@test length(c) == 65
@test length(msc) == 65
@test length(f) == 65
c, msc, f = NeuroAnalyzer.coherence(rand(100), rand(100), fs=10, method=:fft)
@test length(c) == 65
@test length(msc) == 65
@test length(f) == 65
c, msc, f = NeuroAnalyzer.coherence(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1, method=:mt)
@test size(c) == (1, 2049, 1)
@test size(msc) == (1, 2049, 1)
@test length(f) == 2049
c, msc, f = NeuroAnalyzer.coherence(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1, method=:fft)
@test size(c) == (1, 2049, 1)
@test size(msc) == (1, 2049, 1)
@test length(f) == 2049

@info "Test 17/$ntests: frqinst()"
f = NeuroAnalyzer.frqinst(rand(100))
@test length(f) == 100
f = NeuroAnalyzer.frqinst(rand(10, 100, 10))
@test size(f) == (10, 100, 10)
f = NeuroAnalyzer.frqinst(e10)
@test size(f) == (23, 2560, 10)

@info "Test 18/$ntests: ged()"
s, r, rn = NeuroAnalyzer.ged(rand(10, 10), rand(10, 10))
@test length(s) == 100
@test length(r) == 10
@test length(rn) == 10
s, r, rn = NeuroAnalyzer.ged(e10, e10)
@test length(s) == 588800
@test length(r) == 230
@test length(rn) == 230

@info "Test 19/$ntests: erop()"
p, f = erop(e10, ch=1, method=:welch)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10, ch=1, method=:fft)
@test size(p) == (1281, 1)
@test length(f) == 1281
p, f = erop(e10, ch=1, method=:stft)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10, ch=1, method=:mt)
@test size(p) == (1281, 1)
@test length(f) == 1281
p, f = erop(e10, ch=1, method=:mw)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10, ch=1, method=:gh)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10, ch=1, method=:cwt)
@test size(p) == (131, 1)
@test length(f) == 131

@info "Test 20/$ntests: acor()"
@test NeuroAnalyzer.acor(v) == [-0.32 -0.32 -0.08 0.32 0.8 0.32 -0.08 -0.32 -0.32;;;]
ac, l = NeuroAnalyzer.acor(e10)
@test size(ac) == (23, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acor(e10, biased=false)
@test size(ac) == (23, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acor(e10, method=:cor)
@test size(ac) == (23, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acor(e10, method=:stat)
@test size(ac) == (23, 3, 10)
@test length(l) == 3

@info "Test 21/$ntests: ispc()"
iv, ia, sd, pd, s1p, s2p = NeuroAnalyzer.ispc(v1, v2)
@test iv ≈ 0.6125992852305387
@test ia ≈ -0.0017801930770334254
@test sd == [5, 3, 1, -1, -3]
@test pd ≈ [-1.3157044982273682, 0.8713795327960081, 0.3743702916488456, 0.7615478999167377, -1.0328436072470222]
@test s1p ≈ [1.039406675134543, -0.6027563879589182, -0.21331750626000984, -0.33140501474755424, 0.3279718365439915]
@test s2p ≈ [-0.27629782309282525, 0.26862314483709, 0.16105278538883577, 0.4301428851691834, -0.7048717707030308]
iv, ia = NeuroAnalyzer.ispc(e10)
@test size(iv) == (23, 23, 10)
@test size(ia) == (23, 23, 10)
iv, ia, sd, pd, s1p, s2p = NeuroAnalyzer.ispc(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1)
@test iv ≈ [0.1649363573346739;;]
@test ia ≈ [-2.757457461687163;;]
@test size(sd) == (1, 2560, 1)
@test size(pd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test 22/$ntests: itpc()"
iv, izv, ia, ip = NeuroAnalyzer.itpc(ones(1, 10, 10), t=1)
@test iv == 1.0
@test izv == 10.0
@test ia == 0.0
@test ip == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
iv, izv, ia, ip = NeuroAnalyzer.itpc(e10, ch=1, t=256)
@test iv ≈ [0.3987026452273473]
@test izv ≈ [1.5896379931128393]
@test ia ≈ [2.636121475719976]
@test ip[1] ≈ 1.6921361543218993

@info "Test 23/$ntests: itpc_spec()"
iv, izv, f = NeuroAnalyzer.itpc_spec(e10, ch=1, frq_lim=(0, 4), frq_n=5)
@test size(iv) == (5, 2560)
@test size(izv) == (5, 2560)
@test f == [0.01, 0.045, 0.2, 0.894, 4.0]

@info "Test 24/$ntests: mdiff()"
st, sts, p = NeuroAnalyzer.mdiff(m1, m2, method=:absdiff)
@test length(st) == 6
@test sts == 3.0
@test p in [0.0, 1.0]
st, sts, p = NeuroAnalyzer.mdiff(a1, a2, method=:absdiff)
@test size(st) == (2, 6)
@test sts == [1.0, 1.0]
@test p == [0.0, 0.0]
st, sts, p = NeuroAnalyzer.mdiff(m1, m2, method=:diff2int)
@test length(st) == 6
@test sts == 4.666666666666666
@test p == 1.0 || p == 0.0
st, sts, p = NeuroAnalyzer.mdiff(a1, a2, method=:diff2int)
@test size(st) == (2, 6)
@test sts == [2.0, 2.0]
@test p == [0.0, 0.0]
st, sts, p = NeuroAnalyzer.mdiff(e10, e10, method=:absdiff)
@test size(st) == (10, 69)
@test sts == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test p == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
st, sts, p = NeuroAnalyzer.mdiff(e10, e10, method=:diff2int)
@test size(st) == (10, 69)
@test sts == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test p == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

@info "Test 25/$ntests: mutual_information()"
@test NeuroAnalyzer.mutual_information(v1, v2) ≈ 0.4199730940219748
@test NeuroAnalyzer.mutual_information(a1) == [0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0]
@test NeuroAnalyzer.mutual_information(a1, a2) == [0.0 0.0; 0.0 0.0]
m = NeuroAnalyzer.mutual_information(e10)
@test size(m) == (23, 23, 10)
m = NeuroAnalyzer.mutual_information(e10, e10, ch1=1, ch2=2)
@test size(m) == (1, 10)

@info "Test 26/$ntests: msci95()"
@test NeuroAnalyzer.msci95(v1) == (sm = 3.0, ss = 0.7071067811865476, su = 4.385929291125633, sl = 1.6140707088743669)
@test NeuroAnalyzer.msci95(v2) == (sm = 4.0, ss = 0.7071067811865476, su = 5.385929291125633, sl = 2.614070708874367)
@test NeuroAnalyzer.msci95(m1) == (sm = [2.5, 3.5, 4.5], ss = [1.4999999999999998, 1.4999999999999998, 1.4999999999999998], su = [5.4399999999999995, 6.4399999999999995, 7.4399999999999995], sl = [-0.4399999999999995, 0.5600000000000005, 1.5600000000000005])
@test NeuroAnalyzer.msci95(a1) == (sm = [1.0 1.0 1.0; 1.0 1.0 1.0], ss = [0.0 0.0 0.0; 0.0 0.0 0.0], su = [1.0 1.0 1.0; 1.0 1.0 1.0], sl = [1.0 1.0 1.0; 1.0 1.0 1.0])
sm, ss, su, sl = NeuroAnalyzer.msci95(e10)
@test size(sm) == (10, 2560)
@test size(ss) == (10, 2560)
@test size(su) == (10, 2560)
@test size(sl) == (10, 2560)
sm, ss, su, sl = NeuroAnalyzer.msci95(e10, method=:boot)
@test size(sm) == (10, 2560)
@test size(ss) == (10, 2560)
@test size(su) == (10, 2560)
@test size(sl) == (10, 2560)
@test NeuroAnalyzer.msci95(v1, v2) == (sm = -1.0, ss = 1.0, su = 0.96, sl = -2.96)
@test NeuroAnalyzer.msci95(m1, m2) == (sm = [-4.0; 2.0;;], ss = [0.8164965809277261; 0.8164965809277261;;], su = [-2.3996667013816566; 3.6003332986183434;;], sl = [-5.600333298618343; 0.39966670138165683;;])
@test NeuroAnalyzer.msci95(a1, a2) == (sm = [1.0 1.0; 1.0 1.0], ss = [0.0 0.0; 0.0 0.0], su = [1.0 1.0; 1.0 1.0], sl = [1.0 1.0; 1.0 1.0])
sm, ss, su, sl = NeuroAnalyzer.msci95(e10, e10)
@test size(sm) == (23, 10)
@test size(ss) == (23, 10)
@test size(su) == (23, 10)
@test size(sl) == (23, 10)

@info "Test 27/$ntests: eros()"
s, f, t = eros(e10, ch=1, method=:stft)
@test size(s) == (129, 289, 1)
@test length(f) == 129
@test length(t) == 289
s, f, t = eros(e10, ch=1, method=:mt)
@test size(s) == (257, 15, 1)
@test length(f) == 257
@test length(t) == 15
s, f, t = eros(e10, ch=1, method=:mw)
@test size(s) == (129, 2560, 1)
@test length(f) == 129
@test length(t) == 2560
s, f, t = eros(e10, ch=1, method=:gh)
@test size(s) == (129, 2560, 1)
@test length(f) == 129
@test length(t) == 2560
s, f, t = eros(e10, ch=1, method=:cwt)
@test size(s) == (131, 2560, 1)
@test length(f) == 131
@test length(t) == 2560

@info "Test 28/$ntests: phdiff()"
@test NeuroAnalyzer.phdiff(a1, avg=:phase, h=true) == zeros(2, 3, 2)
p = NeuroAnalyzer.phdiff(e10, avg=:phase)
@test size(p) == (23, 1281, 10)
p = NeuroAnalyzer.phdiff(e10, avg=:phase)
@test size(p) == (23, 1281, 10)
p = NeuroAnalyzer.phdiff(e10, avg=:phase, h=true)
@test size(p) == (23, 2560, 10)
p = NeuroAnalyzer.phdiff(e10, avg=:phase, h=true)
@test size(p) == (23, 2560, 10)

@info "Test 29/$ntests: pli()"
pv, phd, s1ph, s2ph = pli(v1, v2)
@test pv == 0.2
@test phd == [5, 3, 1, -1, -3]
@test length(s1ph) == 5
@test length(s2ph) == 5
pv = NeuroAnalyzer.pli(e10);
@test size(pv) == (23, 23, 10)
pv, sd, phd, s1p, s2p = NeuroAnalyzer.pli(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1)
@test pv == [0.1390625;;]
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test 30/$ntests: psd()"
p, f = psd(rand(100), fs=10, wlen=10, woverlap=0)
@test length(p) == 6
@test f == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
p, f = psd(rand(10, 100), fs=10, wlen=10, woverlap=0)
@test size(p) == (10, 6)
p, f = psd(rand(10, 100, 10), fs=10, wlen=10, woverlap=0)
@test size(p) == (10, 6, 10)
p, f = psd(rand(100), fs=10, method=:mt)
@test length(p) == 51
@test round.(f, digits=3) == 0.0:0.1:5.0
p, f = psd(rand(10, 100), fs=10, method=:mt)
@test size(p) == (10, 51)
p, f = psd(rand(10, 100, 10), fs=10, method=:mt)
@test size(p) == (10, 51, 10)
p, f = NeuroAnalyzer.psd(e10)
@test size(p) == (23, 129, 10)
@test f == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10, method=:fft)
@test size(p) == (23, 1281, 10)
@test round.(f, digits=3) == 0.0:0.1:128.0
p, f = NeuroAnalyzer.psd(e10, method=:stft)
@test size(p) == (23, 129, 10)
@test round.(f, digits=3) == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10, method=:mt)
@test size(p) == (23, 1281, 10)
@test round.(f, digits=3) == 0.0:0.1:128.0
p, f = NeuroAnalyzer.psd(e10, method=:mw)
@test size(p) == (23, 129, 10)
@test round.(f, digits=3) == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10, method=:gh)
@test size(p) == (23, 129, 10)
@test round.(f, digits=3) == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10, method=:cwt)
@test size(p) == (23, 131, 10)
@test length(f) == 131

@info "Test 31/$ntests: env_cor()"
e1, t = NeuroAnalyzer.tenv(e10)
e2, t = NeuroAnalyzer.tenv(e10)
ec, p = NeuroAnalyzer.env_cor(e1, e2)
@test ec[1] <= 1.0
@test p[1] <= 1.0

@info "Test 32/$ntests: psd_rel()"
p, f = psd_rel(rand(100), fs=10, frq_lim=(0, 1), wlen=10, woverlap=0)
@test length(p) == 6
@test f == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
p, f = psd_rel(rand(10, 100), fs=10, frq_lim=(0, 1), wlen=10, woverlap=0)
@test size(p) == (10, 6)
p, f = psd_rel(rand(10, 100, 10), fs=10, frq_lim=(0, 1), wlen=10, woverlap=0)
@test size(p) == (10, 6, 10)
p, f = psd_rel(rand(100), fs=10, method=:mt, frq_lim=(0, 1), wlen=10, woverlap=0)
@test length(p) == 51
@test round.(f, digits=3) == 0.0:0.1:5.0
p, f = psd_rel(rand(10, 100), fs=10, method=:mt, frq_lim=(0, 1), wlen=10, woverlap=0)
@test size(p) == (10, 51)
p, f = psd_rel(rand(10, 100, 10), fs=10, method=:mt, frq_lim=(0, 1), wlen=10, woverlap=0)
@test size(p) == (10, 51, 10)
p, f = NeuroAnalyzer.psd_rel(e10, frq_lim=(0, 1))
@test size(p) == (23, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, method=:mt, frq_lim=(0, 1))
@test size(p) == (23, 1281, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, method=:fft, frq_lim=(0, 1))
@test size(p) == (23, 1281, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, method=:stft, frq_lim=(0, 1))
@test size(p) == (23, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, method=:mw, frq_lim=(0, 1))
@test size(p) == (23, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, method=:gh, frq_lim=(0, 1))
@test size(p) == (23, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, method=:cwt, frq_lim=(0, 1))
@test size(p) == (23, 131, 10)
@test f[1] == 0.0
@test f[end] == 86.63

@info "Test 33/$ntests: psd_slope()"
lf, ls, pf = psd_slope(rand(100), fs=10, wlen=10, woverlap=0)
@test length(lf) == 6
@test pf == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
lf, ls, pf = psd_slope(rand(10, 100), fs=10, wlen=10, woverlap=0)
@test size(lf) == (10, 6, 1)
@test length(ls) == 10
lf, ls, pf = psd_slope(rand(10, 100, 10), fs=10, wlen=10, woverlap=0)
@test size(lf) == (10, 6, 10)
@test size(ls) == (10, 10)
lf, ls, pf = psd_slope(e10)
@test size(lf) == (23, 129, 10)
@test size(ls) == (23, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, method=:stft)
@test size(lf) == (23, 129, 10)
@test size(ls) == (23, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, method=:fft)
@test size(lf) == (23, 1281, 10)
@test size(ls) == (23, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, method=:mt)
@test size(lf) == (23, 1281, 10)
@test size(ls) == (23, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, method=:mw)
@test size(lf) == (23, 129, 10)
@test size(ls) == (23, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, method=:gh)
@test size(lf) == (23, 129, 10)
@test size(ls) == (23, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, method=:cwt)
@test size(lf) == (23, 131, 10)
@test size(ls) == (23, 10)
@test pf[1] == 0.0
@test pf[end] == 86.63

@info "Test 34/$ntests: amp()"
p, r, p2p, semi_p2p, msa, rmsa, nrg, rms = NeuroAnalyzer.amp(e10)
@test size(p) == (23, 10)
@test size(r) == (23, 10)
@test size(p2p) == (23, 10)
@test size(semi_p2p) == (23, 10)
@test size(msa) == (23, 10)
@test size(rmsa) == (23, 10)
@test size(nrg) == (23, 10)
@test size(rms) == (23, 10)

@info "Test 35/$ntests: rmse()"
@test NeuroAnalyzer.rmse(v1, v2) == 1.0
@test NeuroAnalyzer.rmse(m1, m2) == [0.0; 0.0;;]
@test NeuroAnalyzer.rmse(a1, a2) == [0.0 0.0; 0.0 0.0]
@test NeuroAnalyzer.rmse(e10, e10) == zeros(23, 10)

@info "Test 36/$ntests: snr()"
@test NeuroAnalyzer.snr(v1) == 1.8973665961010275
@test NeuroAnalyzer.snr2(v1) == 1.2060453783110545
sn, f = NeuroAnalyzer.snr(e10, type=:rms)
sn, f = NeuroAnalyzer.snr(e10, type=:mean)
@test size(sn) == (23, 1281)
@test length(f) == 1281

@info "Test 37/$ntests: spectrogram()"
sp, sf, st = NeuroAnalyzer.spectrogram(e10, method=:stft)
@test size(sp) == (129, 289, 23, 10)
@test length(sf) == 129
@test length(st) == 289
sp, sf, st = NeuroAnalyzer.spectrogram(e10, method=:mt)
@test size(sp) == (257, 15, 23, 10)
@test length(sf) == 257
@test length(st) == 15
sp, sf, st = NeuroAnalyzer.spectrogram(e10, method=:mw)
@test size(sp) == (129, 2560, 23, 10)
@test length(sf) == 129
@test length(st) == 2560
sp, sf, st = NeuroAnalyzer.spectrogram(e10, method=:gh)
@test size(sp) == (129, 2560, 23, 10)
@test length(sf) == 129
@test length(st) == 2560
sp, sf, st = NeuroAnalyzer.spectrogram(e10, method=:cwt)
@test size(sp) == (131, 2560, 23, 10)
@test length(sf) == 131
@test length(st) == 2560

@info "Test 38/$ntests: spec_seg()"
sp, sf, st = NeuroAnalyzer.spectrogram(e10)
sp, sst, t, f = spec_seg(sp, sf, st, ch=1, t=(0, 1), f=(0, 10))
@test size(sp) == (11, 30, 10)
@test t == (1, 30)
@test f == (1, 11)

@info "Test 39/$ntests: spectrum()"
c, sa, sp, sph = NeuroAnalyzer.spectrum(rand(100))
@test length(c) == 51
@test length(sa) == 51
@test length(sp) == 51
@test length(sph) == 51
c, sa, sp, sph = NeuroAnalyzer.hspectrum(rand(100))
@test length(c) == 100
@test length(sa) == 100
@test length(sp) == 100
@test length(sph) == 100
c, sa, sp, sph = NeuroAnalyzer.spectrum(rand(10, 100, 10))
@test size(c) == (10, 51, 10)
@test size(sa) == (10, 51, 10)
@test size(sp) == (10, 51, 10)
@test size(sph) == (10, 51, 10)
c, sa, sp, sph = NeuroAnalyzer.hspectrum(rand(10, 100, 10))
@test size(c) == (10, 100, 10)
@test size(sa) == (10, 100, 10)
@test size(sp) == (10, 100, 10)
@test size(sph) == (10, 100, 10)
c, sa, sp, sph = NeuroAnalyzer.spectrum(rand(10, 100, 10), h=true)
@test size(c) == (10, 100, 10)
@test size(sa) == (10, 100, 10)
@test size(sp) == (10, 100, 10)
@test size(sph) == (10, 100, 10)
c, sa, sp, sph = NeuroAnalyzer.spectrum(e10)
@test size(c) == (23, 1281, 10)
@test size(sa) == (23, 1281, 10)
@test size(sp) == (23, 1281, 10)
@test size(sph) == (23, 1281, 10)
c, sa, sp, sph = NeuroAnalyzer.spectrum(e10, h=true)
@test size(c) == (23, 2560, 10)
@test size(sa) == (23, 2560, 10)
@test size(sp) == (23, 2560, 10)
@test size(sph) == (23, 2560, 10)

@info "Test 40/$ntests: stationarity()"
s = NeuroAnalyzer.stationarity(e10, method=:adf)
@test size(s) == (23, 2, 10)
s = NeuroAnalyzer.stationarity(e10, method=:cov)
@test size(s) == (257, 10)
s = NeuroAnalyzer.stationarity(e10, method=:hilbert)
@test size(s) == (23, 2559, 10)
s = NeuroAnalyzer.stationarity(e10, method=:mean)
@test size(s) == (23, 10, 10)
s = NeuroAnalyzer.stationarity(e10, method=:var)
@test size(s) == (23, 10, 10)

@info "Test 41/$ntests: channel_stats()"
c = NeuroAnalyzer.channel_stats(e10)
for idx in 1:length(c)
    @test size(c[idx]) == (24, 10)
end

@info "Test 42/$ntests: epoch_stats()"
e = NeuroAnalyzer.epoch_stats(e10)
for idx in 1:length(e)
    @test length(e[idx]) == 10
end

@info "Test 43/$ntests: cpsd()"
pxy, f = cpsd(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1, method=:mt)
@test size(pxy) == (1, 2049, 1)
@test length(f) == 2049
pxy, f = cpsd(e10, e10, ch1=1, ch2=2, ep1=1, ep2=1, method=:fft)
@test size(pxy) == (1, 1290, 1)
@test length(f) == 1290

@info "Test 44/$ntests: tkeo()"
@test NeuroAnalyzer.tkeo(v1) == [1.0, 1.0, 1.0, 1.0, 5.0]
@test NeuroAnalyzer.tkeo(a1) == [1.0 0.0 1.0; 1.0 0.0 1.0;;; 1.0 0.0 1.0; 1.0 0.0 1.0]
t = NeuroAnalyzer.tkeo(e10, method=:pow)
@test size(t) == (23, 2560, 10)
t = NeuroAnalyzer.tkeo(e10, method=:der)
@test size(t) == (23, 2560, 10)
t = NeuroAnalyzer.tkeo(e10, method=:amp)
@test size(t) == (23, 2560, 10)

@info "Test 45/$ntests: total_power()"
tp = NeuroAnalyzer.total_power(e10)
@test size(tp) == (23, 10)
tp = NeuroAnalyzer.total_power(e10, method=:welch)
@test size(tp) == (23, 10)
tp = NeuroAnalyzer.total_power(e10, method=:fft)
@test size(tp) == (23, 10)
tp = NeuroAnalyzer.total_power(e10, method=:stft)
@test size(tp) == (23, 10)
tp = NeuroAnalyzer.total_power(e10, method=:mt)
@test size(tp) == (23, 10)
tp = NeuroAnalyzer.total_power(e10, method=:mw)
@test size(tp) == (23, 10)
tp = NeuroAnalyzer.total_power(e10, method=:gh)
@test size(tp) == (23, 10)
tp = NeuroAnalyzer.total_power(e10, method=:cwt)
@test size(tp) == (23, 10)

@info "Test 46/$ntests: pacor()"
pac, l = pacor(e10, l=2)
@test size(pac) == (23, 5, 10)
@test length(l) == 5

@info "Test 47/$ntests: xcov()"
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2, biased=false)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2, method=:cov, biased=false)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2, method=:cov, biased=true)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2, method=:stat)
@test size(xc) == (1, 3, 1)
@test length(l) == 3

@info "Test 48/$ntests: xcor()"
xc, l = NeuroAnalyzer.xcor(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcor(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2, biased=false)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcor(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2, method=:cor)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcor(e10, e10, ch1=1, ch2=2, ep1=1, ep2=2, method=:stat)
@test size(xc) == (1, 3, 1)
@test length(l) == 3

@info "Test 49/$ntests: amp_at()"
e = NeuroAnalyzer.erp(e10)
@test size(amp_at(e, t=2)) == (23, 11)

@info "Test 50/$ntests: avgamp_at()"
@test size(avgamp_at(e, t=(2, 2.5))) == (23, 11)

@info "Test 51/$ntests: maxamp_at()"
@test size(maxamp_at(e, t=(2, 2.5))) == (23, 11)

@info "Test 52/$ntests: minamp_at()"
@test size(minamp_at(e, t=(2, 2.5))) == (23, 11)

@info "Test 53/$ntests: env_up()"
x = rand(-10:0.1:10, 1000)
t = linspace(0, 10, 1000)
@test length(env_up(x, t)) == 1000

@info "Test 54/$ntests: env_lo()"
x = rand(-10:0.1:10, 1000)
t = linspace(0, 10, 1000)
@test length(env_lo(x, t)) == 1000

@info "Test 55/$ntests: henv_up()"
x = rand(-10:0.1:10, 1000)
@test length(henv_up(x)) == 1000

@info "Test 56/$ntests: henv_lo()"
x = rand(-10:0.1:10, 1000)
@test length(henv_lo(x)) == 1000

@info "Test 57/$ntests: axc2frq()"
x = [1, -2, 3, -4, 5]
y = [-1, 2, -3, 4, -5]
xc = xcor(x, y, l=4, demean=false)
l = collect(-4:4)
f = axc2frq(xc[1, :, 1], l)
@test length(f) == 1

@info "Test 58/$ntests: hjorth()"
h_act, h_mob, h_comp = hjorth(v1)
@test h_act == 2.5
@test h_mob == 0.0
@test isnan(h_comp)
h_act, h_mob, h_comp = hjorth(a1)
@test h_act == zeros(2, 2)
h_act, h_mob, h_comp = hjorth(e10)
@test size(h_act) == (23, 10)
@test size(h_mob) == (23, 10)
@test size(h_comp) == (23, 10)

@info "Test 58/$ntests: hjorth()"
pf = peak_frq(e10, f=(8, 13))
@test size(pf) == (23, 10)

@info "Test 59/$ntests: phsd()"
ph, f = phsd(e10)
@test size(ph) == (23, 1281, 10)
@test length(f) == 1281

@info "Test 60/$ntests: band_asymmetry()"
@test band_asymmetry(e10, ch1=1, ch2=1, frq_lim=(0, 10)) == (ba = 0.0, ba_norm = 0.0)

@info "Test 61/$ntests: symmetry()"
@test symmetry(v) == 5
@test symmetry(m1) == [3.0; 3.0;;]
@test symmetry(e10, ch=1) == [1.001563721657545 0.9527078565980168 1.0285261489698891 0.9452887537993921 1.0 0.970746728252502 0.9219219219219219 1.0496397117694156 1.0173364854215918 0.9393939393939394]

true
