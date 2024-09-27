using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
NeuroAnalyzer.filter!(e10, ch="all", fprototype=:fir, ftype=:lp, cutoff=40, order=8)
NeuroAnalyzer.filter!(e10, ch="all", fprototype=:fir, ftype=:hp, cutoff=1, order=8)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "Test: acov()"
@test NeuroAnalyzer.acov(v) == [-0.8 -0.8 -0.2 0.8 2.0 0.8 -0.2 -0.8 -0.8;;;]
ac, l = NeuroAnalyzer.acov(e10, ch="all")
@test size(ac) == (24, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acov(e10, ch="all", biased=false)
@test size(ac) == (24, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acov(e10, ch="all", method=:cov)
@test size(ac) == (24, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acov(e10, ch="all", method=:cov, biased=false)
@test size(ac) == (24, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acov(e10, ch="all", method=:stat)
@test size(ac) == (24, 3, 10)
@test length(l) == 3

@info "Test: ampdiff()"
@test size(NeuroAnalyzer.ampdiff(a1)) == (2, 3, 2)
ad = NeuroAnalyzer.ampdiff(e10, ch="all")
@test size(ad) == (24, 2560, 10)

@info "Test: band_power()"
@test size(NeuroAnalyzer.band_power(e10, ch="Fp1", frq_lim=(10, 20))) == (1, 10)
@test size(NeuroAnalyzer.band_power(e10, ch="Fp1", frq_lim=(10, 20), method=:welch)) == (1, 10)
@test size(NeuroAnalyzer.band_power(e10, ch="Fp1", frq_lim=(10, 20), method=:fft)) == (1, 10)
@test size(NeuroAnalyzer.band_power(e10, ch="Fp1", frq_lim=(10, 20), method=:stft)) == (1, 10)
@test size(NeuroAnalyzer.band_power(e10, ch="Fp1", frq_lim=(10, 20), method=:mt)) == (1, 10)
@test size(NeuroAnalyzer.band_power(e10, ch="Fp1", frq_lim=(10, 20), method=:mw)) == (1, 10)
@test size(NeuroAnalyzer.band_power(e10, ch="Fp1", frq_lim=(10, 20), method=:gh)) == (1, 10)
@test size(NeuroAnalyzer.band_power(e10, ch="Fp1", frq_lim=(10, 20), method=:cwt)) == (1, 10)

@info "Test: band_mpower()"
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, ch="Fp1", frq_lim=(10, 20))
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, ch="Fp1", frq_lim=(10, 20), method=:mt)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, ch="Fp1", frq_lim=(10, 20), method=:stft)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, ch="Fp1", frq_lim=(10, 20), method=:fft)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, ch="Fp1", frq_lim=(10, 20), method=:mw)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, ch="Fp1", frq_lim=(10, 20), method=:gh)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
mbp, maxf, maxbp = NeuroAnalyzer.band_mpower(e10, ch="Fp1", frq_lim=(10, 20), method=:cwt)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)

@info "Test: corm()"
@test NeuroAnalyzer.corm(v) ≈ ones(5, 5)
@test size(NeuroAnalyzer.corm(a1)) == (2, 2, 3, 2)
@test size(NeuroAnalyzer.corm(e10, ch="all")) == (24, 24, 2560, 10)

@info "Test: covm()"
@test NeuroAnalyzer.covm(v) == [ 2.5  5.0  7.5 10.0 12.5;
                                 5.0 10.0 15.0 20.0 25.0;
                                 7.5 15.0 22.5 30.0 37.5;
                                10.0 20.0 30.0 40.0 50.0;
                                12.5 25.0 37.5 50.0 62.5]
@test size(NeuroAnalyzer.covm(a1)) == (2, 2, 3, 2)
@test size(NeuroAnalyzer.covm(e10, ch="all")) == (24, 24, 2560, 10)

@info "Test: cph()"
ph, f = NeuroAnalyzer.cph(rand(10), rand(10), fs=1)
@test length(ph) == 9
@test length(f) == 9
ph, f = NeuroAnalyzer.cph(rand(10, 10, 2), fs=1)
@test size(ph) == (10, 10, 9, 2)
@test length(f) == 9
ph, f = NeuroAnalyzer.cph(e10, ch="all")
@test size(ph) == (24, 24, 2049, 10)
@test length(f) == 2049
ph, f = NeuroAnalyzer.cph(rand(10, 10, 2), rand(10, 10, 2), fs=1)
@test size(ph) == (10, 9, 2)
@test length(f) == 9
ph, f = NeuroAnalyzer.cph(e10, e10, ch1=["Fp1", "Fp2"], ch2=["Fp1", "Fp2"], ep1=1, ep2=1)
@test size(ph) == (2, 2049, 1)
@test length(f) == 2049

@info "Test: diss()"
@test NeuroAnalyzer.diss(v1, v2) == (gd = 0.21320071635561044, sc = 0.9772727272727273)
@test NeuroAnalyzer.diss(a1) == (gd = [0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0], sc = [1.0 1.0; 1.0 1.0;;; 1.0 1.0; 1.0 1.0])
gd, sc = NeuroAnalyzer.diss(a1, a2)
@test size(gd) == (2, 2)
@test size(sc) == (2, 2)
gd, sc = NeuroAnalyzer.diss(e10, ch="all")
@test size(gd) == (24, 24, 10)
@test size(sc) == (24, 24, 10)

@info "Test: entropy()"
e, s, l = NeuroAnalyzer.entropy(rand(10))
@test e < l
@test s < l
e, s, l = NeuroAnalyzer.entropy(e10, ch="all")
@test size(e) == (24, 10)
@test size(s) == (24, 10)
@test size(l) == (24, 10)

@info "Test: negentropy()"
n = NeuroAnalyzer.negentropy(rand(10))
@test n < 0
n = NeuroAnalyzer.negentropy(eeg, ch="all")
@test size(n) == (24, 1)

@info "Test: tenv()"
e, t = NeuroAnalyzer.tenv(e10, ch="all")
@test size(e) == (24, 2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_mean(e10, ch="all", dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_mean(e10, ch="all", dims=2)
@test size(em) == (2560, 24)
@test size(eu) == (2560, 24)
@test size(el) == (2560, 24)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_mean(e10, ch="all", dims=3)
@test size(em) == (2560,)
@test size(eu) == (2560,)
@test size(el) == (2560,)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_median(e10, ch="all", dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_median(e10, ch="all", dims=2)
@test size(em) == (2560, 24)
@test size(eu) == (2560, 24)
@test size(el) == (2560, 24)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.tenv_median(e10, ch="all", dims=3)
@test size(em) == (2560,)
@test size(eu) == (2560,)
@test size(el) == (2560,)
@test length(t) == 2560

@info "Test: senv()"
e, t = NeuroAnalyzer.senv(e10, ch="Fp1")
@test size(e) == (1, 289, 10)
@test length(t) == 289
e, t = NeuroAnalyzer.senv(e10, ch="Fp1", method=:mt)
@test size(e) == (1, 15, 10)
@test length(t) == 15
e, t = NeuroAnalyzer.senv(e10, ch="Fp1", method=:mw)
@test size(e) == (1, 2560, 10)
@test length(t) == 2560
e, t = NeuroAnalyzer.senv(e10, ch="Fp1", method=:gh)
@test size(e) == (1, 2560, 10)
@test length(t) == 2560
e, t = NeuroAnalyzer.senv(e10, ch="Fp1", method=:cwt)
@test size(e) == (1, 2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.senv_mean(e10, ch="all", dims=1)
@test size(em) == (289, 10)
@test size(eu) == (289, 10)
@test size(el) == (289, 10)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_mean(e10, ch="all", dims=2)
@test size(em) == (289, 24)
@test size(eu) == (289, 24)
@test size(el) == (289, 24)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_mean(e10, ch="all", dims=3)
@test size(em) == (289,)
@test size(eu) == (289,)
@test size(el) == (289,)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_median(e10, ch="all", dims=1)
@test size(em) == (289, 10)
@test size(eu) == (289, 10)
@test size(el) == (289, 10)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_median(e10, ch="all", dims=2)
@test size(em) == (289, 24)
@test size(eu) == (289, 24)
@test size(el) == (289, 24)
@test length(t) == 289
em, eu, el, t = NeuroAnalyzer.senv_median(e10, ch="all", dims=3)
@test size(em) == (289,)
@test size(eu) == (289,)
@test size(el) == (289,)
@test length(t) == 289

@info "Test: penv()"
e, t = NeuroAnalyzer.penv(e10, ch="Fp1", method=:welch)
@test size(e) == (1, 129, 10)
@test length(t) == 129
e, t = NeuroAnalyzer.penv(e10, ch="Fp1", method=:fft)
@test size(e) == (1, 1281, 10)
@test length(t) == 1281
e, t = NeuroAnalyzer.penv(e10, ch="Fp1", method=:stft)
@test size(e) == (1, 129, 10)
@test length(t) == 129
e, t = NeuroAnalyzer.penv(e10, ch="Fp1", method=:mt)
@test size(e) == (1, 1281, 10)
@test length(t) == 1281
e, t = NeuroAnalyzer.penv(e10, ch="Fp1", method=:mw)
@test size(e) == (1, 129, 10)
@test length(t) == 129
e, t = NeuroAnalyzer.penv(e10, ch="Fp1", method=:gh)
@test size(e) == (1, 129, 10)
@test length(t) == 129
e, t = NeuroAnalyzer.penv(e10, ch="Fp1", method=:cwt)
@test size(e) == (1, 131, 10)
@test length(t) == 131
em, eu, el, t = NeuroAnalyzer.penv_mean(e10, ch="all", dims=1)
@test size(em) == (129, 10)
@test size(eu) == (129, 10)
@test size(el) == (129, 10)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_mean(e10, ch="all", dims=2)
@test size(em) == (129, 24)
@test size(eu) == (129, 24)
@test size(el) == (129, 24)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_mean(e10, ch="all", dims=3)
@test size(em) == (129,)
@test size(eu) == (129,)
@test size(el) == (129,)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_median(e10, ch="all", dims=1)
@test size(em) == (129, 10)
@test size(eu) == (129, 10)
@test size(el) == (129, 10)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_median(e10, ch="all", dims=2)
@test size(em) == (129, 24)
@test size(eu) == (129, 24)
@test size(el) == (129, 24)
@test length(t) == 129
em, eu, el, t = NeuroAnalyzer.penv_median(e10, ch="all", dims=3)
@test size(em) == (129,)
@test size(eu) == (129,)
@test size(el) == (129,)
@test length(t) == 129

@info "Test: henv()"
e, t = NeuroAnalyzer.henv(e10, ch="all")
@test size(e) == (24, 2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_mean(e10, ch="all", dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_mean(e10, ch="all", dims=2)
@test size(em) == (2560, 24)
@test size(eu) == (2560, 24)
@test size(el) == (2560, 24)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_mean(e10, ch="all", dims=3)
@test size(em) == (2560,)
@test size(eu) == (2560,)
@test size(el) == (2560,)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_median(e10, ch="all", dims=1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_median(e10, ch="all", dims=2)
@test size(em) == (2560, 24)
@test size(eu) == (2560, 24)
@test size(el) == (2560, 24)
@test length(t) == 2560
em, eu, el, t = NeuroAnalyzer.henv_median(e10, ch="all", dims=3)
@test size(em) == (2560,)
@test size(eu) == (2560,)
@test size(el) == (2560,)
@test length(t) == 2560

@info "Test: erp_peaks()"
e = NeuroAnalyzer.average_epochs(e10)
p = NeuroAnalyzer.erp_peaks(e)
@test size(p) == (19, 2)

@info "Test: coherence()"
c, msc, f = NeuroAnalyzer.coherence(rand(100), rand(100), fs=10, method=:mt)
@test length(c) == 65
@test length(msc) == 65
@test length(f) == 65
c, msc, f = NeuroAnalyzer.coherence(rand(100), rand(100), fs=10, method=:fft)
@test length(c) == 65
@test length(msc) == 65
@test length(f) == 65
c, msc, f = NeuroAnalyzer.coherence(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=1, method=:mt)
@test size(c) == (1, 2049, 1)
@test size(msc) == (1, 2049, 1)
@test length(f) == 2049
c, msc, f = NeuroAnalyzer.coherence(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=1, method=:fft)
@test size(c) == (1, 2049, 1)
@test size(msc) == (1, 2049, 1)
@test length(f) == 2049

@info "Test: frqinst()"
f = NeuroAnalyzer.frqinst(rand(100))
@test length(f) == 100
f = NeuroAnalyzer.frqinst(rand(10, 100, 10))
@test size(f) == (10, 100, 10)
f = NeuroAnalyzer.frqinst(e10, ch="all")
@test size(f) == (24, 2560, 10)

@info "Test: ged()"
s, r, rn = NeuroAnalyzer.ged(rand(10, 10), rand(10, 10))
@test length(s) == 100
@test length(r) == 10
@test length(rn) == 10
s, r, rn = NeuroAnalyzer.ged(e10, e10, ch1="Fp1", ch2="Fp2")
@test length(s) == 25600
@test length(r) == 10
@test length(rn) == 10

@info "Test: erop()"
p, f = erop(e10, ch="Fp1", method=:welch)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10, ch="Fp1", method=:fft)
@test size(p) == (1281, 1)
@test length(f) == 1281
p, f = erop(e10, ch="Fp1", method=:stft)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10, ch="Fp1", method=:mt)
@test size(p) == (1281, 1)
@test length(f) == 1281
p, f = erop(e10, ch="Fp1", method=:mw)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10, ch="Fp1", method=:gh)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10, ch="Fp1", method=:cwt)
@test size(p) == (131, 1)
@test length(f) == 131

@info "Test: acor()"
@test NeuroAnalyzer.acor(v) == [-0.32 -0.32 -0.08 0.32 0.8 0.32 -0.08 -0.32 -0.32;;;]
ac, l = NeuroAnalyzer.acor(e10, ch="all")
@test size(ac) == (24, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acor(e10, ch="all", biased=false)
@test size(ac) == (24, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acor(e10, ch="all", method=:cor)
@test size(ac) == (24, 3, 10)
@test length(l) == 3
ac, l = NeuroAnalyzer.acor(e10, ch="all", method=:stat)
@test size(ac) == (24, 3, 10)
@test length(l) == 3

@info "Test: ispc()"
iv, ia, sd, pd, s1p, s2p = NeuroAnalyzer.ispc(v1, v2)
@test iv ≈ 0.6125992852305387
@test ia ≈ -0.0017801930770334254
@test sd == [5, 3, 1, -1, -3]
@test pd ≈ [-1.3157044982273682, 0.8713795327960081, 0.3743702916488456, 0.7615478999167377, -1.0328436072470222]
@test s1p ≈ [1.039406675134543, -0.6027563879589182, -0.21331750626000984, -0.33140501474755424, 0.3279718365439915]
@test s2p ≈ [-0.27629782309282525, 0.26862314483709, 0.16105278538883577, 0.4301428851691834, -0.7048717707030308]
iv, ia = NeuroAnalyzer.ispc(e10, ch="all")
@test size(iv) == (24, 24, 10)
@test size(ia) == (24, 24, 10)
iv, ia, sd, pd, s1p, s2p = NeuroAnalyzer.ispc(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=1)
@test iv ≈ [0.1649363573346739;;]
@test ia ≈ [-2.757457461687163;;]
@test size(sd) == (1, 2560, 1)
@test size(pd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: itpc()"
iv, izv, ia, ip = NeuroAnalyzer.itpc(ones(1, 10, 10), t=1)
@test iv == 1.0
@test izv == 10.0
@test ia == 0.0
@test ip == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
iv, izv, ia, ip = NeuroAnalyzer.itpc(e10, ch="Fp1", t=256)
@test iv ≈ [0.3987026452273473]
@test izv ≈ [1.5896379931128393]
@test ia ≈ [2.636121475719976]
@test ip[1] ≈ 1.6921361543218993

@info "Test: itpc_spec()"
iv, izv, f = NeuroAnalyzer.itpc_spec(e10, ch="Fp1", frq_lim=(0, 4), frq_n=5)
@test size(iv) == (5, 2560)
@test size(izv) == (5, 2560)
@test f == [0.01, 0.045, 0.2, 0.894, 4.0]

@info "Test: mdiff()"
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
st, sts, p = NeuroAnalyzer.mdiff(e10, e10, ch1="Fp1", ch2="Fp1", method=:absdiff)
@test size(st) == (10, 3)
@test sts == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test p == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
st, sts, p = NeuroAnalyzer.mdiff(e10, e10, ch1="Fp1", ch2="Fp1", method=:diff2int)
@test size(st) == (10, 3)
@test sts == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test p == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

@info "Test: mutual_information()"
@test NeuroAnalyzer.mutual_information(v1, v2) ≈ 0.4199730940219748
@test NeuroAnalyzer.mutual_information(a1) == [0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0]
@test NeuroAnalyzer.mutual_information(a1, a2) == [0.0 0.0; 0.0 0.0]
m = NeuroAnalyzer.mutual_information(e10, ch="all")
@test size(m) == (24, 24, 10)
m = NeuroAnalyzer.mutual_information(e10, e10, ch1="Fp1", ch2="Fp2")
@test size(m) == (1, 10)

@info "Test: msci95()"
@test NeuroAnalyzer.msci95(v1) == (sm = 3.0, ss = 0.7071067811865476, su = 4.385929291125633, sl = 1.6140707088743669)
@test NeuroAnalyzer.msci95(v2) == (sm = 4.0, ss = 0.7071067811865476, su = 5.385929291125633, sl = 2.614070708874367)
@test NeuroAnalyzer.msci95(m1) == (sm = [2.5, 3.5, 4.5], ss = [1.4999999999999998, 1.4999999999999998, 1.4999999999999998], su = [5.4399999999999995, 6.4399999999999995, 7.4399999999999995], sl = [-0.4399999999999995, 0.5600000000000005, 1.5600000000000005])
@test NeuroAnalyzer.msci95(a1) == (sm = [1.0 1.0 1.0; 1.0 1.0 1.0], ss = [0.0 0.0 0.0; 0.0 0.0 0.0], su = [1.0 1.0 1.0; 1.0 1.0 1.0], sl = [1.0 1.0 1.0; 1.0 1.0 1.0])
sm, ss, su, sl = NeuroAnalyzer.msci95(e10, ch="all")
@test size(sm) == (10, 2560)
@test size(ss) == (10, 2560)
@test size(su) == (10, 2560)
@test size(sl) == (10, 2560)
sm, ss, su, sl = NeuroAnalyzer.msci95(e10, ch="all", method=:boot)
@test size(sm) == (10, 2560)
@test size(ss) == (10, 2560)
@test size(su) == (10, 2560)
@test size(sl) == (10, 2560)
@test NeuroAnalyzer.msci95(v1, v2) == (sm = -1.0, ss = 1.0, su = 0.96, sl = -2.96)
@test NeuroAnalyzer.msci95(m1, m2) == (sm = [-4.0; 2.0;;], ss = [0.8164965809277261; 0.8164965809277261;;], su = [-2.3996667013816566; 3.6003332986183434;;], sl = [-5.600333298618343; 0.39966670138165683;;])
@test NeuroAnalyzer.msci95(a1, a2) == (sm = [1.0 1.0; 1.0 1.0], ss = [0.0 0.0; 0.0 0.0], su = [1.0 1.0; 1.0 1.0], sl = [1.0 1.0; 1.0 1.0])
sm, ss, su, sl = NeuroAnalyzer.msci95(e10, e10, ch1="Fp1", ch2="Fp2")
@test size(sm) == (1, 10)
@test size(ss) == (1, 10)
@test size(su) == (1, 10)
@test size(sl) == (1, 10)

@info "Test: eros()"
s, f, t = eros(e10, ch="Fp1", method=:stft)
@test size(s) == (129, 289, 1)
@test length(f) == 129
@test length(t) == 289
s, f, t = eros(e10, ch="Fp1", method=:mt)
@test size(s) == (257, 15, 1)
@test length(f) == 257
@test length(t) == 15
s, f, t = eros(e10, ch="Fp1", method=:mw)
@test size(s) == (129, 2560, 1)
@test length(f) == 129
@test length(t) == 2560
s, f, t = eros(e10, ch="Fp1", method=:gh)
@test size(s) == (129, 2560, 1)
@test length(f) == 129
@test length(t) == 2560
s, f, t = eros(e10, ch="Fp1", method=:cwt)
@test size(s) == (131, 2560, 1)
@test length(f) == 131
@test length(t) == 2560

@info "Test: phdiff()"
@test NeuroAnalyzer.phdiff(a1, avg=:phase, h=true) == zeros(2, 3, 2)
p = NeuroAnalyzer.phdiff(e10, ch="all", avg=:phase)
@test size(p) == (24, 1281, 10)
p = NeuroAnalyzer.phdiff(e10, ch="all", avg=:phase)
@test size(p) == (24, 1281, 10)
p = NeuroAnalyzer.phdiff(e10, ch="all", avg=:phase, h=true)
@test size(p) == (24, 2560, 10)
p = NeuroAnalyzer.phdiff(e10, ch="all", avg=:phase, h=true)
@test size(p) == (24, 2560, 10)

@info "Test: pli()"
pv, phd, s1ph, s2ph = pli(v1, v2)
@test pv == 0.2
@test phd == [5, 3, 1, -1, -3]
@test length(s1ph) == 5
@test length(s2ph) == 5
pv = NeuroAnalyzer.pli(e10, ch="all");
@test size(pv) == (24, 24, 10)
pv, sd, phd, s1p, s2p = NeuroAnalyzer.pli(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=1)
@test pv == [0.1390625;;]
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: psd()"
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
p, f = NeuroAnalyzer.psd(e10, ch="Fp1")
@test size(p) == (1, 129, 10)
@test f == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10, ch="Fp1", method=:fft)
@test size(p) == (1, 1281, 10)
@test round.(f, digits=3) == 0.0:0.1:128.0
p, f = NeuroAnalyzer.psd(e10, ch="Fp1", method=:stft)
@test size(p) == (1, 129, 10)
@test round.(f, digits=3) == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10, ch="Fp1", method=:mt)
@test size(p) == (1, 1281, 10)
@test round.(f, digits=3) == 0.0:0.1:128.0
p, f = NeuroAnalyzer.psd(e10, ch="Fp1", method=:mw)
@test size(p) == (1, 129, 10)
@test round.(f, digits=3) == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10, ch="Fp1", method=:gh)
@test size(p) == (1, 129, 10)
@test round.(f, digits=3) == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10, ch="Fp1", method=:cwt)
@test size(p) == (1, 131, 10)
@test length(f) == 131

@info "Test: env_cor()"
e1, t = NeuroAnalyzer.tenv(e10, ch="all")
e2, t = NeuroAnalyzer.tenv(e10, ch="all")
ec, p = NeuroAnalyzer.env_cor(e1, e2)
@test ec[1] <= 1.0
@test p[1] <= 1.0

@info "Test: psd_rel()"
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
p, f = NeuroAnalyzer.psd_rel(e10, ch="Fp1", frq_lim=(0, 1))
@test size(p) == (1, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, ch="Fp1", method=:mt, frq_lim=(0, 1))
@test size(p) == (1, 1281, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, ch="Fp1", method=:fft, frq_lim=(0, 1))
@test size(p) == (1, 1281, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, ch="Fp1", method=:stft, frq_lim=(0, 1))
@test size(p) == (1, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, ch="Fp1", method=:mw, frq_lim=(0, 1))
@test size(p) == (1, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, ch="Fp1", method=:gh, frq_lim=(0, 1))
@test size(p) == (1, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = NeuroAnalyzer.psd_rel(e10, ch="Fp1", method=:cwt, frq_lim=(0, 1))
@test size(p) == (1, 131, 10)
@test f[1] == 0.0
@test f[end] == 86.63

@info "Test: psd_slope()"
lf, ls, pf = psd_slope(rand(100), fs=10, wlen=10, woverlap=0)
@test length(lf) == 6
@test pf == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
lf, ls, pf = psd_slope(rand(10, 100, 10), fs=10, wlen=10, woverlap=0)
@test size(lf) == (10, 6, 10)
@test size(ls) == (10, 10)
lf, ls, pf = psd_slope(e10, ch="Fp1")
@test size(lf) == (1, 129, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, ch="Fp1", method=:stft)
@test size(lf) == (1, 129, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, ch="Fp1", method=:fft)
@test size(lf) == (1, 1281, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, ch="Fp1", method=:mt)
@test size(lf) == (1, 1281, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, ch="Fp1", method=:mw)
@test size(lf) == (1, 129, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, ch="Fp1", method=:gh)
@test size(lf) == (1, 129, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10, ch="Fp1", method=:cwt)
@test size(lf) == (1, 131, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 86.63

@info "Test: amp()"
p, r, p2p, semi_p2p, msa, rmsa, nrg, rms = NeuroAnalyzer.amp(e10, ch="all")
@test size(p) == (24, 10)
@test size(r) == (24, 10)
@test size(p2p) == (24, 10)
@test size(semi_p2p) == (24, 10)
@test size(msa) == (24, 10)
@test size(rmsa) == (24, 10)
@test size(nrg) == (24, 10)
@test size(rms) == (24, 10)

@info "Test: rmse()"
@test NeuroAnalyzer.rmse(v1, v2) == 1.0
@test NeuroAnalyzer.rmse(a1, a2) == [0.0 0.0; 0.0 0.0]
@test NeuroAnalyzer.rmse(e10, e10, ch1="Fp1", ch2="Fp2") == zeros(1, 10)

@info "Test: snr()"
@test NeuroAnalyzer.snr(v1) == 1.8973665961010275
@test NeuroAnalyzer.snr2(v1) == 1.2060453783110545
sn, f = NeuroAnalyzer.snr(e10, ch="all", type=:rms)
sn, f = NeuroAnalyzer.snr(e10, ch="all", type=:mean)
@test size(sn) == (24, 1281)
@test length(f) == 1281

@info "Test: spectrogram()"
sp, sf, st = NeuroAnalyzer.spectrogram(e10, ch="Fp1", method=:stft)
@test size(sp) == (129, 289, 1, 10)
@test length(sf) == 129
@test length(st) == 289
sp, sf, st = NeuroAnalyzer.spectrogram(e10, ch="Fp1", method=:mt)
@test size(sp) == (257, 15, 1, 10)
@test length(sf) == 257
@test length(st) == 15
sp, sf, st = NeuroAnalyzer.spectrogram(e10, ch="Fp1", method=:mw)
@test size(sp) == (129, 2560, 1, 10)
@test length(sf) == 129
@test length(st) == 2560
sp, sf, st = NeuroAnalyzer.spectrogram(e10, ch="Fp1", method=:gh)
@test size(sp) == (129, 2560, 1, 10)
@test length(sf) == 129
@test length(st) == 2560
sp, sf, st = NeuroAnalyzer.spectrogram(e10, ch="Fp1", method=:cwt)
@test size(sp) == (131, 2560, 1, 10)
@test length(sf) == 131
@test length(st) == 2560

@info "Test: spec_seg()"
sp, sf, st = NeuroAnalyzer.spectrogram(e10, ch="all")
sp, sst, t, f = spec_seg(sp, sf, st, ch=1, t=(0, 1), f=(0, 10))
@test size(sp) == (11, 30, 10)
@test t == (1, 30)
@test f == (1, 11)

@info "Test: spectrum()"
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
c, sa, sp, sph = NeuroAnalyzer.spectrum(e10, ch="all")
@test size(c) == (24, 1281, 10)
@test size(sa) == (24, 1281, 10)
@test size(sp) == (24, 1281, 10)
@test size(sph) == (24, 1281, 10)
c, sa, sp, sph = NeuroAnalyzer.spectrum(e10, ch="all", h=true)
@test size(c) == (24, 2560, 10)
@test size(sa) == (24, 2560, 10)
@test size(sp) == (24, 2560, 10)
@test size(sph) == (24, 2560, 10)

@info "Test: stationarity()"
s = NeuroAnalyzer.stationarity(e10, ch="all", method=:adf)
@test size(s) == (24, 2, 10)
s = NeuroAnalyzer.stationarity(e10, ch="all", method=:cov)
@test size(s) == (257, 10)
s = NeuroAnalyzer.stationarity(e10, ch="all", method=:hilbert)
@test size(s) == (24, 2559, 10)
s = NeuroAnalyzer.stationarity(e10, ch="all", method=:mean)
@test size(s) == (24, 10, 10)
s = NeuroAnalyzer.stationarity(e10, ch="all", method=:var)
@test size(s) == (24, 10, 10)

@info "Test: channel_stats()"
c = NeuroAnalyzer.channel_stats(e10)
for idx in eachindex(c)
    @test size(c[idx]) == (24, 10)
end

@info "Test: epoch_stats()"
e = NeuroAnalyzer.epoch_stats(e10)
for idx in eachindex(e)
    @test length(e[idx]) == 10
end

@info "Test: cpsd()"
pxy, f = cpsd(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=1, method=:mt)
@test size(pxy) == (1, 2049, 1)
@test length(f) == 2049
pxy, f = cpsd(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=1, method=:fft)
@test size(pxy) == (1, 1290, 1)
@test length(f) == 1290

@info "Test: tkeo()"
@test NeuroAnalyzer.tkeo(v1) == [1.0, 1.0, 1.0, 1.0, 5.0]
@test NeuroAnalyzer.tkeo(a1) == [1.0 0.0 1.0; 1.0 0.0 1.0;;; 1.0 0.0 1.0; 1.0 0.0 1.0]
t = NeuroAnalyzer.tkeo(e10, ch="all", method=:pow)
@test size(t) == (24, 2560, 10)
t = NeuroAnalyzer.tkeo(e10, ch="all", method=:der)
@test size(t) == (24, 2560, 10)
t = NeuroAnalyzer.tkeo(e10, ch="all", method=:amp)
@test size(t) == (24, 2560, 10)

@info "Test: total_power()"
tp = NeuroAnalyzer.total_power(e10, ch="Fp1")
@test size(tp) == (1, 10)
tp = NeuroAnalyzer.total_power(e10, ch="Fp1", method=:welch)
@test size(tp) == (1, 10)
tp = NeuroAnalyzer.total_power(e10, ch="Fp1", method=:fft)
@test size(tp) == (1, 10)
tp = NeuroAnalyzer.total_power(e10, ch="Fp1", method=:stft)
@test size(tp) == (1, 10)
tp = NeuroAnalyzer.total_power(e10, ch="Fp1", method=:mt)
@test size(tp) == (1, 10)
tp = NeuroAnalyzer.total_power(e10, ch="Fp1", method=:mw)
@test size(tp) == (1, 10)
tp = NeuroAnalyzer.total_power(e10, ch="Fp1", method=:gh)
@test size(tp) == (1, 10)
tp = NeuroAnalyzer.total_power(e10, ch="Fp1", method=:cwt)
@test size(tp) == (1, 10)

@info "Test: pacor()"
pac, l = pacor(e10, ch="all", l=2)
@test size(pac) == (24, 5, 10)
@test length(l) == 5

@info "Test: xcov()"
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=2)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=2, biased=false)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=2, method=:cov, biased=false)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=2, method=:cov, biased=true)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcov(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=2, method=:stat)
@test size(xc) == (1, 3, 1)
@test length(l) == 3

@info "Test: xcor()"
xc, l = NeuroAnalyzer.xcor(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=2)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcor(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=2, biased=false)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcor(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=2, method=:cor)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = NeuroAnalyzer.xcor(e10, e10, ch1="Fp1", ch2="Fp2", ep1=1, ep2=2, method=:stat)
@test size(xc) == (1, 3, 1)
@test length(l) == 3

@info "Test: amp_at()"
e = NeuroAnalyzer.average_epochs(e10)
@test size(amp_at(e, t=2)) == (19, 11)

@info "Test: avgamp_at()"
@test size(avgamp_at(e, t=(2, 2.5))) == (19, 11)

@info "Test: maxamp_at()"
@test size(maxamp_at(e, t=(2, 2.5))) == (19, 11)

@info "Test: minamp_at()"
@test size(minamp_at(e, t=(2, 2.5))) == (19, 11)

@info "Test: env_up()"
x = rand(-10:0.1:10, 1000)
t = linspace(0, 10, 1000)
@test length(env_up(x, t)) == 1000

@info "Test: env_lo()"
x = rand(-10:0.1:10, 1000)
t = linspace(0, 10, 1000)
@test length(env_lo(x, t)) == 1000

@info "Test: henv_up()"
x = rand(-10:0.1:10, 1000)
@test length(henv_up(x)) == 1000

@info "Test: henv_lo()"
x = rand(-10:0.1:10, 1000)
@test length(henv_lo(x)) == 1000

@info "Test: axc2frq()"
x = [1, -2, 3, -4, 5]
y = [-1, 2, -3, 4, -5]
xc = xcor(x, y, l=4, demean=false)
l = collect(-4:4)
f = axc2frq(xc[1, :, 1], l)
@test length(f) == 1

@info "Test: hjorth()"
h_act, h_mob, h_comp = hjorth(v1)
@test h_act == 2.5
@test h_mob == 0.0
@test isnan(h_comp)
h_act, h_mob, h_comp = hjorth(a1)
@test h_act == zeros(2, 2)
h_act, h_mob, h_comp = hjorth(e10, ch="all")
@test size(h_act) == (24, 10)
@test size(h_mob) == (24, 10)
@test size(h_comp) == (24, 10)

@info "Test: hjorth()"
pf = peak_frq(e10, ch="all", f=(8, 13))
@test size(pf) == (24, 10)

@info "Test: phsd()"
ph, f = phsd(e10, ch="all")
@test size(ph) == (24, 1281, 10)
@test length(f) == 1281

@info "Test: band_asymmetry()"
@test band_asymmetry(e10, ch1="Fp1", ch2="Fp1", frq_lim=(0, 10)) == (ba = 0.0, ba_norm = 0.0)

@info "Test: symmetry()"
@test symmetry(v) == 5
@test symmetry(e10, ch="Fp1") == [1.001563721657545 0.9527078565980168 1.0285261489698891 0.9452887537993921 1.0 0.970746728252502 0.9219219219219219 1.0496397117694156 1.0173364854215918 0.9393939393939394]

true
