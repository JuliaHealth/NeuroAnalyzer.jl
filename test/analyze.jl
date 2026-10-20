using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets
using DSP

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)
NeuroAnalyzer.filter!(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :lp,
    cutoff = 40,
    order = 91,
)
NeuroAnalyzer.filter!(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :hp,
    cutoff = 1,
    order = 91,
)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "Test: acov()"
@test NeuroAnalyzer._acov(v) ≈ [-0.8, -0.8, -0.2, 0.8, 2.0, 0.8, -0.2, -0.8, -0.8]
ac, l = acov(e10; ch = "all")
@test size(ac) == (24, 69, 10)
@test length(l) == 69
ac, l = acov(e10; ch = "all", biased = false)
@test size(ac) == (24, 69, 10)
@test length(l) == 69
ac, l = acov(e10; ch = "all", method = :cov)
@test size(ac) == (24, 69, 10)
@test length(l) == 69
ac, l = acov(e10; ch = "all", method = :cov, biased = false)
@test size(ac) == (24, 69, 10)
@test length(l) == 69
ac, l = acov(e10; ch = "all", method = :stat)
@test size(ac) == (24, 69, 10)
@test length(l) == 69

@info "Test: ampdiff()"
@test size(NeuroAnalyzer._ampdiff(a1)) == (2, 3, 2)
ad = ampdiff(e10; ch = "all")
@test size(ad) == (24, 2560, 10)

@info "Test: band_power()"
@test size(band_power(e10, ch = "Fp1", flim = (10, 20))) == (1, 10)
@test size(band_power(e10, ch = "Fp1", flim = (10, 20), method = :welch)) == (1, 10)
@test size(band_power(e10, ch = "Fp1", flim = (10, 20), method = :fft)) == (1, 10)
@test size(band_power(e10, ch = "Fp1", flim = (10, 20), method = :stft)) == (1, 10)
@test size(band_power(e10, ch = "Fp1", flim = (10, 20), method = :mt)) == (1, 10)
@test size(band_power(e10, ch = "Fp1", flim = (10, 20), method = :mw)) == (1, 10)
@test size(band_power(e10, ch = "Fp1", flim = (10, 20), method = :gh)) == (1, 10)

@info "Test: band_mpower()"
mbp, maxf, maxbp, maxba = band_mpower(e10; ch = "Fp1", flim = (10, 20))
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
@test size(maxba) == (1, 10)
mbp, maxf, maxbp, maxba = band_mpower(e10; ch = "Fp1", flim = (10, 20), method = :mt)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
@test size(maxba) == (1, 10)
mbp, maxf, maxbp, maxba = band_mpower(e10; ch = "Fp1", flim = (10, 20), method = :stft)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
@test size(maxba) == (1, 10)
mbp, maxf, maxbp, maxba = band_mpower(e10; ch = "Fp1", flim = (10, 20), method = :fft)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
@test size(maxba) == (1, 10)
mbp, maxf, maxbp, maxba = band_mpower(e10; ch = "Fp1", flim = (10, 20), method = :mw)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
@test size(maxba) == (1, 10)
mbp, maxf, maxbp, maxba = band_mpower(e10; ch = "Fp1", flim = (10, 20), method = :gh)
@test size(mbp) == (1, 10)
@test size(maxf) == (1, 10)
@test size(maxbp) == (1, 10)
@test size(maxba) == (1, 10)

@info "Test: corm()"
@test size(NeuroAnalyzer._corm(a1)) == (2, 2, 2)
@test size(corm(e10, ch = "all")) == (24, 24, 10)

@info "Test: covm()"
@test NeuroAnalyzer._covm(v) == [2.5;;]
@test size(NeuroAnalyzer._covm(a1)) == (2, 2, 2)
@test size(covm(e10, ch = "all")) == (24, 24, 10)

@info "Test: cph()"
ph, f = NeuroAnalyzer._cph(rand(10), rand(10); fs = 1)
@test length(ph) == 9
@test length(f) == 9
ph, f = NeuroAnalyzer._cph(rand(10, 10, 2); fs = 1)
@test size(ph) == (10, 10, 9, 2)
@test length(f) == 9
ph, f = cph(e10; ch = "all")
@test size(ph) == (24, 24, 2049, 10)
@test length(f) == 2049
ph, f = NeuroAnalyzer._cph(rand(10, 10, 2), rand(10, 10, 2); fs = 1)
@test size(ph) == (10, 9, 2)
@test length(f) == 9
ph, f = cph(e10, e10; ch1 = ["Fp1", "Fp2"], ch2 = ["Fp1", "Fp2"], ep1 = 1, ep2 = 1)
@test size(ph) == (2, 2049, 1)
@test length(f) == 2049

@info "Test: entropy()"
e, sh, l, s, ns = NeuroAnalyzer._entropy(rand(100))
@test e < l
@test sh < l
@test s > ns

e, sh, l, s, ns = NeuroAnalyzer.entropy(e10; ch = "all")
@test size(e) == (24, 10)
@test size(sh) == (24, 10)
@test size(l) == (24, 10)
@test size(s) == (24, 10)
@test size(ns) == (24, 10)

@info "Test: negentropy()"
n = NeuroAnalyzer._negentropy(rand(10000))
@test n > 0
n = NeuroAnalyzer.negentropy(eeg; ch = "all")
@test size(n) == (24, 1)

@info "Test: tenv()"
e, t = tenv(e10; ch = "all")
@test size(e) == (24, 2560, 10)
@test length(t) == 2560
em, eu, el, t = tenv_mean(e10; ch = "all", dims = 1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = tenv_mean(e10; ch = "all", dims = 2)
@test size(em) == (2560, 24)
@test size(eu) == (2560, 24)
@test size(el) == (2560, 24)
@test length(t) == 2560
em, eu, el, t = tenv_mean(e10; ch = "all", dims = 3)
@test size(em) == (2560, 1)
@test size(eu) == (2560, 1)
@test size(el) == (2560, 1)
@test length(t) == 2560
em, eu, el, t = tenv_median(e10; ch = "all", dims = 1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = tenv_median(e10; ch = "all", dims = 2)
@test size(em) == (2560, 24)
@test size(eu) == (2560, 24)
@test size(el) == (2560, 24)
@test length(t) == 2560
em, eu, el, t = tenv_median(e10; ch = "all", dims = 3)
@test size(em) == (2560, 1)
@test size(eu) == (2560, 1)
@test size(el) == (2560, 1)
@test length(t) == 2560

@info "Test: senv()"
e, t = senv(e10; ch = "Fp1")
@test size(e) == (1, 89, 10)
@test length(t) == 89
e, t = senv(e10; ch = "Fp1", method = :mt)
@test size(e) == (1, 89, 10)
@test length(t) == 89
e, t = senv(e10; ch = "Fp1", method = :mw)
@test size(e) == (1, 2560, 10)
@test length(t) == 2560
e, t = senv(e10; ch = "Fp1", method = :gh)
@test size(e) == (1, 2560, 10)
@test length(t) == 2560
em, eu, el, t = senv_mean(e10; ch = "all", dims = 1)
@test size(em) == (89, 10)
@test size(eu) == (89, 10)
@test size(el) == (89, 10)
@test length(t) == 89
em, eu, el, t = senv_mean(e10; ch = "all", dims = 2)
@test size(em) == (89, 24)
@test size(eu) == (89, 24)
@test size(el) == (89, 24)
@test length(t) == 89
em, eu, el, t = senv_mean(e10; ch = "all", dims = 3)
@test size(em) == (89, 1)
@test size(eu) == (89, 1)
@test size(el) == (89, 1)
@test length(t) == 89
em, eu, el, t = senv_median(e10; ch = "all", dims = 1)
@test size(em) == (89, 10)
@test size(eu) == (89, 10)
@test size(el) == (89, 10)
@test length(t) == 89
em, eu, el, t = senv_median(e10; ch = "all", dims = 2)
@test size(em) == (89, 24)
@test size(eu) == (89, 24)
@test size(el) == (89, 24)
@test length(t) == 89
em, eu, el, t = senv_median(e10; ch = "all", dims = 3)
@test size(em) == (89, 1)
@test size(eu) == (89, 1)
@test size(el) == (89, 1)
@test length(t) == 89

@info "Test: penv()"
e, t = penv(e10; ch = "Fp1", method = :welch)
@test size(e) == (1, 129, 10)
@test length(t) == 129
e, t = penv(e10; ch = "Fp1", method = :fft)
@test size(e) == (1, 1281, 10)
@test length(t) == 1281
e, t = penv(e10; ch = "Fp1", method = :stft)
@test size(e) == (1, 129, 10)
@test length(t) == 129
e, t = penv(e10; ch = "Fp1", method = :mt)
@test size(e) == (1, 1281, 10)
@test length(t) == 1281
e, t = penv(e10; ch = "Fp1", method = :mw)
@test size(e) == (1, 128, 10)
@test length(t) == 128
e, t = penv(e10; ch = "Fp1", method = :gh)
@test size(e) == (1, 129, 10)
@test length(t) == 129
em, eu, el, t = penv_mean(e10; ch = "all", dims = 1)
@test size(em) == (129, 10)
@test size(eu) == (129, 10)
@test size(el) == (129, 10)
@test length(t) == 129
em, eu, el, t = penv_mean(e10; ch = "all", dims = 2)
@test size(em) == (129, 24)
@test size(eu) == (129, 24)
@test size(el) == (129, 24)
@test length(t) == 129
em, eu, el, t = penv_mean(e10; ch = "all", dims = 3)
@test size(em) == (129, 1)
@test size(eu) == (129, 1)
@test size(el) == (129, 1)
@test length(t) == 129
em, eu, el, t = penv_median(e10; ch = "all", dims = 1)
@test size(em) == (129, 10)
@test size(eu) == (129, 10)
@test size(el) == (129, 10)
@test length(t) == 129
em, eu, el, t = penv_median(e10; ch = "all", dims = 2)
@test size(em) == (129, 24)
@test size(eu) == (129, 24)
@test size(el) == (129, 24)
@test length(t) == 129
em, eu, el, t = penv_median(e10; ch = "all", dims = 3)
@test size(em) == (129, 1)
@test size(eu) == (129, 1)
@test size(el) == (129, 1)
@test length(t) == 129

@info "Test: henv()"
e, t = henv(e10; ch = "all")
@test size(e) == (24, 2560, 10)
@test length(t) == 2560
em, eu, el, t = henv_mean(e10; ch = "all", dims = 1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = henv_mean(e10; ch = "all", dims = 2)
@test size(em) == (2560, 24)
@test size(eu) == (2560, 24)
@test size(el) == (2560, 24)
@test length(t) == 2560
em, eu, el, t = henv_mean(e10; ch = "all", dims = 3)
@test size(em) == (2560, 1)
@test size(eu) == (2560, 1)
@test size(el) == (2560, 1)
@test length(t) == 2560
em, eu, el, t = henv_median(e10; ch = "all", dims = 1)
@test size(em) == (2560, 10)
@test size(eu) == (2560, 10)
@test size(el) == (2560, 10)
@test length(t) == 2560
em, eu, el, t = henv_median(e10; ch = "all", dims = 2)
@test size(em) == (2560, 24)
@test size(eu) == (2560, 24)
@test size(el) == (2560, 24)
@test length(t) == 2560
em, eu, el, t = henv_median(e10; ch = "all", dims = 3)
@test size(em) == (2560, 1)
@test size(eu) == (2560, 1)
@test size(el) == (2560, 1)
@test length(t) == 2560

@info "Test: erp_peaks()"
e = average_epochs(e10)
p = erp_peaks(e)
@test size(p) == (19, 2)

@info "Test: erp_auc()"
auc_data = erp_auc(e; ch = "all")
@test length(auc_data) == 19

@info "Test: coherence()"
coh_data = NeuroAnalyzer._coherence(rand(100), rand(100); fs = 10, wstep = 5, method = :mt)
@test length(coh_data.coh) == 51
@test length(coh_data.imcoh) == 51
@test length(coh_data.msc) == 51
@test length(coh_data.f) == 51
coh_data = NeuroAnalyzer._coherence(rand(100), rand(100); fs = 10, wstep = 5, method = :fft)
@test length(coh_data.coh) == 51
@test length(coh_data.imcoh) == 51
@test length(coh_data.msc) == 51
@test length(coh_data.f) == 51
coh_data = NeuroAnalyzer._coherence(rand(100), rand(100); fs = 10, wstep = 5, method = :stft)
@test length(coh_data.coh) == 9
@test length(coh_data.imcoh) == 9
@test length(coh_data.msc) == 9
@test length(coh_data.f) == 9
coh_data = NeuroAnalyzer.coherence(
    e10,
    e10;
    ch1 = "Fp1",
    ch2 = "Fp2",
    method = :mt,
)
@test size(coh_data.coh) == (1, 1281)
@test size(coh_data.imcoh) == (1, 1281)
@test size(coh_data.msc) == (1, 1281)
@test length(coh_data.f) == 1281
coh_data = NeuroAnalyzer.coherence(
    e10,
    e10;
    ch1 = "Fp1",
    ch2 = "Fp2",
    method = :fft,
)
@test size(coh_data.coh) == (1, 1281)
@test size(coh_data.imcoh) == (1, 1281)
@test size(coh_data.msc) == (1, 1281)
@test length(coh_data.f) == 1281
coh_data = NeuroAnalyzer.coherence(
    e10,
    e10;
    ch1 = "Fp1",
    ch2 = "Fp2",
    method = :stft,
)
@test size(coh_data.coh) == (1, 129)
@test size(coh_data.imcoh) == (1, 129)
@test size(coh_data.msc) == (1, 129)
@test length(coh_data.f) == 129

@info "Test: frqinst()"
f = NeuroAnalyzer._frqinst(rand(100))
@test length(f) == 100
f = NeuroAnalyzer._frqinst(rand(10, 100, 10))
@test size(f) == (10, 100, 10)
f = NeuroAnalyzer.frqinst(e10; ch = "all")
@test size(f) == (24, 2560, 10)

@info "Test: ged()"
s, r, rn = NeuroAnalyzer._ged(rand(10, 10), rand(10, 10))
@test length(s) == 100
@test length(r) == 10
@test length(rn) == 10
s, r, rn = ged(e10, e10; ch1 = "Fp1", ch2 = "Fp2")
@test length(s) == 25600
@test length(r) == 10
@test length(rn) == 10

@info "Test: erop()"
p, f = erop(e10; ch = "Fp1", method = :welch)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10; ch = "Fp1", method = :fft)
@test size(p) == (1281, 1)
@test length(f) == 1281
p, f = erop(e10; ch = "Fp1", method = :stft)
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10; ch = "Fp1", method = :mt)
@test size(p) == (1281, 1)
@test length(f) == 1281
p, f = erop(e10; ch = "Fp1", method = :mw)
@test size(p) == (128, 1)
@test length(f) == 128
p, f = erop(e10; ch = "Fp1", method = :gh)
@test size(p) == (129, 1)
@test length(f) == 129

@info "Test: acor()"
# normalized to 1 at lag 0
@test NeuroAnalyzer._acor(v) ≈ [-0.4, -0.4, -0.1, 0.4, 1.0, 0.4, -0.1, -0.4, -0.4]
ac, l = acor(e10; ch = "all")
@test size(ac) == (24, 69, 10)
@test length(l) == 69
ac, l = acor(e10; ch = "all", biased = false)
@test size(ac) == (24, 69, 10)
@test length(l) == 69
ac, l = acor(e10; ch = "all", method = :cor)
@test size(ac) == (24, 69, 10)
@test length(l) == 69
ac, l = acor(e10; ch = "all", method = :stat)
@test size(ac) == (24, 69, 10)
@test length(l) == 69

@info "Test: ispc()"
iv, ia, sd, pd, s1p, s2p = NeuroAnalyzer._ispc(v1, v2)
@test iv ≈ 0.61259928523
@test ia ≈ 0.00178019307703
@test sd == [-5.0, -3.0, -1.0, 1.0, 3.0]
@test pd ≈ [
    1.31570449822,
    -0.87137953279,
    -0.37437029164,
    -0.76154789991,
    1.03284360724,
]
@test s1p ≈ [
    1.0394066751,
    -0.60275638795,
    -0.213317506260,
    -0.331405014747,
    0.32797183654,
]
@test s2p ≈ [
    -0.276297823092,
    0.268623144,
    0.161052785388,
    0.43014288516,
    -0.70487177070,
]
iv, ia = ispc(e10; ch = "all")
@test size(iv) == (24, 24, 10)
@test size(ia) == (24, 24, 10)
iv, ia, sd, pd, s1p, s2p = ispc(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test iv ≈ [0.9920495704;;]
@test ia ≈ [0.000134242680366;;]
@test size(sd) == (1, 2560, 1)
@test size(pd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: itpc()"
iv, izv, ia, ip = NeuroAnalyzer._itpc(ones(1, 10, 10); t = 1)
@test iv == 1.0
@test izv == 10.0
@test ia == 0.0
@test ip == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
iv, izv, ia, ip = itpc(e10; ch = "Fp1", t = 1)
@test iv ≈ [0.99979974993]
@test izv ≈ [9.9959953997]
@test ia ≈ [-0.00360014005561]
@test ip[1] ≈ 0.0306443702042

@info "Test: itpc_spec()"
iv, izv, f = itpc_spec(e10; ch = "Fp1", flim = (0, 4), nfrq = 5)
@test size(iv) == (5, 2560)
@test size(izv) == (5, 2560)
@test f == [0.0, 1.0, 2.0, 3.0, 4.0]

@info "Test: mdiff()"
st, sts, p = NeuroAnalyzer._mdiff(m1, m2; method = :absdiff)
@test length(st) == 6
@test sts == 3.0
@test p >= 0.0 && p <= 1.0
st, sts, p = NeuroAnalyzer._mdiff(a1, a2; method = :absdiff)
@test size(st) == (2, 6)
@test sts == [1.0, 1.0]
@test p == [0.0, 0.0]
st, sts, p = NeuroAnalyzer._mdiff(m1, m2; method = :diff2int)
@test length(st) == 6
@test sts ≈ 4.6666666666
@test p >= 0.0 && p <= 1.0
st, sts, p = NeuroAnalyzer._mdiff(a1, a2; method = :diff2int)
@test size(st) == (2, 6)
@test sts == [2.0, 2.0]
@test p == [0.0, 0.0]
st, sts, p = mdiff(e10, e10; ch1 = "Fp1", ch2 = "Fp1", method = :absdiff)
@test size(st) == (10, 3)
@test sts == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test p == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
st, sts, p = mdiff(e10, e10; ch1 = "Fp1", ch2 = "Fp1", method = :diff2int)
@test size(st) == (10, 3)
@test sts == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test p == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

@info "Test: mutual_information()"
@test NeuroAnalyzer._mutual_information(v1, v2) ≈ 0.41997309402
@test NeuroAnalyzer._mutual_information(a1) == [0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0]
@test NeuroAnalyzer._mutual_information(a1, a2) == [0.0 0.0; 0.0 0.0]
m = NeuroAnalyzer.mutual_information(e10; ch = "all")
@test size(m) == (24, 24, 10)
m = NeuroAnalyzer.mutual_information(e10, e10; ch1 = "Fp1", ch2 = "Fp2")
@test size(m) == (1, 10)

@info "Test: msci95()"
@test NeuroAnalyzer._msci95(v1).sm == 3.0
@test NeuroAnalyzer._msci95(v1).se ≈ 0.70710678118
@test NeuroAnalyzer._msci95(v1).ll ≈ 1.61407070887
@test NeuroAnalyzer._msci95(v1).ul ≈ 4.3859292911
@test NeuroAnalyzer._msci95(v2).sm == 4.0
@test NeuroAnalyzer._msci95(v2).se ≈ 0.70710678118
@test NeuroAnalyzer._msci95(v2).ll ≈ 2.6140707088
@test NeuroAnalyzer._msci95(v2).ul ≈ 5.3859292911
@test NeuroAnalyzer._msci95(m1).sm == [2.5, 3.5, 4.5]
@test NeuroAnalyzer._msci95(m1).se ≈ [1.49999999999, 1.49999999999, 1.49999999999]
@test NeuroAnalyzer._msci95(m1).ll ≈ [-0.43999999999, 0.56000000000, 1.56000000000]
@test NeuroAnalyzer._msci95(m1).ul ≈ [5.43999999999, 6.43999999999, 7.43999999999]
@test NeuroAnalyzer._msci95(a1).sm == [1.0 1.0 1.0; 1.0 1.0 1.0]
@test NeuroAnalyzer._msci95(a1).se == [0.0 0.0 0.0; 0.0 0.0 0.0]
@test NeuroAnalyzer._msci95(a1).ll == [1.0 1.0 1.0; 1.0 1.0 1.0]
@test NeuroAnalyzer._msci95(a1).ul == [1.0 1.0 1.0; 1.0 1.0 1.0]

@test NeuroAnalyzer._msci95(v1, v2).sm == -1.0
@test NeuroAnalyzer._msci95(v1, v2).se == 1.0
@test NeuroAnalyzer._msci95(v1, v2).ll == -2.96
@test NeuroAnalyzer._msci95(v1, v2).ul == 0.96
@test NeuroAnalyzer._msci95(m1, m2).sm == [-4.0; 2.0;;]
@test NeuroAnalyzer._msci95(m1, m2).se ≈ [0.81649658092; 0.81649658092;;]
@test NeuroAnalyzer._msci95(m1, m2).ll ≈ [-5.6003332986; 0.399666701381;;]
@test NeuroAnalyzer._msci95(m1, m2).ul ≈ [-2.39966670138; 3.60033329861;;]

@test NeuroAnalyzer._msci95(a1, a2).sm == [1.0 1.0; 1.0 1.0]
@test NeuroAnalyzer._msci95(a1, a2).se == [0.0 0.0; 0.0 0.0]
@test NeuroAnalyzer._msci95(a1, a2).ll == [1.0 1.0; 1.0 1.0]
@test NeuroAnalyzer._msci95(a1, a2).ul == [1.0 1.0; 1.0 1.0]

sm, ss, su, sl = msci95(e10; ch = "all")
@test size(sm) == (10, 2560)
@test size(ss) == (10, 2560)
@test size(su) == (10, 2560)
@test size(sl) == (10, 2560)
sm, ss, su, sl = msci95(e10; ch = "all", method = :boot)
@test size(sm) == (10, 2560)
@test size(ss) == (10, 2560)
@test size(su) == (10, 2560)
@test size(sl) == (10, 2560)
sm, ss, su, sl = msci95(e10, e10; ch1 = "Fp1", ch2 = "Fp2")
@test size(sm) == (1, 10)
@test size(ss) == (1, 10)
@test size(su) == (1, 10)
@test size(sl) == (1, 10)

@info "Test: eros()"
s, f, t = eros(e10; ch = "Fp1", method = :stft)
@test size(s) == (129, 89, 1)
@test length(f) == 129
@test length(t) == 89
s, f, t = eros(e10; ch = "Fp1", method = :mt)
@test size(s) == (129, 89, 1)
@test length(f) == 129
@test length(t) == 89
s, f, t = eros(e10; ch = "Fp1", method = :mw)
@test size(s) == (129, 2560, 1)
@test length(f) == 129
@test length(t) == 2560
s, f, t = eros(e10; ch = "Fp1", method = :gh)
@test size(s) == (129, 2560, 1)
@test length(f) == 129
@test length(t) == 2560
s, f, t = eros(e10; ch = "Fp1", method = :cwt)
@test size(s) == (30, 2560, 1)
@test length(f) == 30
@test length(t) == 2560

@info "Test: phdiff()"
@test NeuroAnalyzer._phdiff(a1, avg = :phase, h = true) == zeros(2, 3, 2)
p = phdiff(e10; ch = "all", avg = :phase)
@test size(p) == (24, 1281, 10)
p = phdiff(e10; ch = "all", avg = :phase)
@test size(p) == (24, 1281, 10)
p = phdiff(e10; ch = "all", avg = :phase, h = true)
@test size(p) == (24, 2560, 10)
p = phdiff(e10; ch = "all", avg = :phase, h = true)
@test size(p) == (24, 2560, 10)

@info "Test: pli()"
pv, sd, phd, s1ph, s2ph = NeuroAnalyzer._pli(v1, v2)
@test 0 <= pv <= 1
@test sd == [-5.0, -3.0, -1.0, 1.0, 3.0]
@test all(-π .< phd .<= π)
@test length(s1ph) == 5
@test length(s2ph) == 5
# constant phase lag: PLI = 1 and wPLI = 1
x1 = sin.(2π * 8 .* (0:511) ./ 256)
x2 = sin.(2π * 8 .* (0:511) ./ 256 .- π / 2)
@test NeuroAnalyzer._pli(x1, x2).pv ≈ 1 atol = 0.05
@test NeuroAnalyzer._wpli(x1, x2).pv ≈ 1 atol = 0.05
pv = pli(e10; ch = "all");
@test size(pv) == (24, 24, 10)
pv, sd, phd, s1p, s2p = pli(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test pv == [0.265625;;]
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: plv()"
pv, phd, s1ph, s2ph = NeuroAnalyzer._plv(v1, v2)
@test pv ≈ 0.612599285
@test phd == [-5.0, -3.0, -1.0, 1.0, 3.0]
@test length(s1ph) == 5
@test length(s2ph) == 5
pv = plv(e10; ch = "all");
@test size(pv) == (24, 24, 10)
pv, sd, phd, s1p, s2p = plv(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test pv ≈ [0.99204957;;]
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: psd()"
p, f = NeuroAnalyzer._psd(rand(100); fs = 10, wlen = 10, woverlap = 5)
@test length(p) == 6
@test f == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
p, f = NeuroAnalyzer._psd(rand(10, 100); fs = 10, wlen = 10, woverlap = 5)
@test size(p) == (10, 6)
p, f = NeuroAnalyzer._psd(rand(10, 100, 10); fs = 10, wlen = 10, woverlap = 5)
@test size(p) == (10, 6, 10)
p, f = NeuroAnalyzer._psd(rand(100); fs = 10, method = :mt, woverlap = 5)
@test length(p) == 51
@test round.(f, digits = 3) == 0.0:0.1:5.0
p, f = NeuroAnalyzer._psd(rand(10, 100); fs = 10, method = :mt, woverlap = 5)
@test size(p) == (10, 51)
p, f = NeuroAnalyzer._psd(rand(10, 100, 10); fs = 10, method = :mt, woverlap = 5)
@test size(p) == (10, 51, 10)
p, f = NeuroAnalyzer.psd(e10; ch = "Fp1")
@test size(p) == (1, 129, 10)
@test f == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10; ch = "Fp1", method = :fft)
@test size(p) == (1, 1281, 10)
@test round.(f, digits = 3) == 0.0:0.1:128.0
p, f = NeuroAnalyzer.psd(e10; ch = "Fp1", method = :stft)
@test size(p) == (1, 129, 10)
@test round.(f, digits = 3) == 0.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10; ch = "Fp1", method = :mt)
@test size(p) == (1, 1281, 10)
@test round.(f, digits = 3) == 0.0:0.1:128.0
p, f = NeuroAnalyzer.psd(e10; ch = "Fp1", method = :mw)
@test size(p) == (1, 128, 10)
@test round.(f, digits = 3) == 1.0:1.0:128.0
p, f = NeuroAnalyzer.psd(e10; ch = "Fp1", method = :gh)
@test size(p) == (1, 129, 10)
@test round.(f, digits = 3) == 0.0:1.0:128.0

@info "Test: env_cor()"
e1, t = tenv(e10; ch = "all")
e2, t = tenv(e10; ch = "all")
ec, p = env_cor(e1, e2)
@test ec[1] <= 1.0
@test p[1] <= 1.0

@info "Test: psd_rel()"
p, f = NeuroAnalyzer._psd_rel(rand(100); fs = 10, flim = (0, 1), wlen = 10, woverlap = 5)
@test length(p) == 6
@test f == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
p, f = NeuroAnalyzer._psd_rel(rand(10, 100); fs = 10, flim = (0, 1), wlen = 10, woverlap = 5)
@test size(p) == (10, 6)
p, f = NeuroAnalyzer._psd_rel(rand(10, 100, 10); fs = 10, flim = (0, 1), wlen = 10, woverlap = 5)
@test size(p) == (10, 6, 10)
p, f = NeuroAnalyzer._psd_rel(rand(100); fs = 10, method = :mt, flim = (0, 1), wlen = 10, woverlap = 5)
@test length(p) == 51
@test round.(f, digits = 3) == 0.0:0.1:5.0
p, f = NeuroAnalyzer._psd_rel(rand(10, 100); fs = 10, method = :mt, flim = (0, 1), wlen = 10, woverlap = 5)
@test size(p) == (10, 51)
p, f = NeuroAnalyzer._psd_rel(
    rand(10, 100, 10);
    fs = 10,
    method = :mt,
    flim = (0, 1),
    wlen = 10,
    woverlap = 5,
)
@test size(p) == (10, 51, 10)
p, f = psd_rel(e10; ch = "Fp1", flim = (0, 1))
@test size(p) == (1, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = psd_rel(e10; ch = "Fp1", method = :mt, flim = (0, 1))
@test size(p) == (1, 1281, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = psd_rel(e10; ch = "Fp1", method = :fft, flim = (0, 1))
@test size(p) == (1, 1281, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = psd_rel(e10; ch = "Fp1", method = :stft, flim = (0, 1))
@test size(p) == (1, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = psd_rel(e10; ch = "Fp1", method = :mw, flim = (0, 1))
@test size(p) == (1, 128, 10)
@test f[1] == 1.0
@test f[end] == 128.0
p, f = psd_rel(e10; ch = "Fp1", method = :gh, flim = (0, 1))
@test size(p) == (1, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0

@info "Test: psd_slope()"
lf, ls, pf = NeuroAnalyzer._psd_slope(rand(100); fs = 10, wlen = 10, woverlap = 5)
@test length(lf) == 6
@test pf == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
lf, ls, pf = NeuroAnalyzer._psd_slope(rand(10, 100, 10); fs = 10, wlen = 10, woverlap = 5)
@test size(lf) == (10, 6, 10)
@test size(ls) == (10, 10)
lf, ls, pf = psd_slope(e10; ch = "Fp1")
@test size(lf) == (1, 129, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10; ch = "Fp1", method = :stft)
@test size(lf) == (1, 129, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10; ch = "Fp1", method = :fft)
@test size(lf) == (1, 1281, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10; ch = "Fp1", method = :mt)
@test size(lf) == (1, 1281, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10; ch = "Fp1", method = :mw)
@test size(lf) == (1, 128, 10)
@test size(ls) == (1, 10)
@test pf[1] == 1.0
@test pf[end] == 128.0
lf, ls, pf = psd_slope(e10; ch = "Fp1", method = :gh)
@test size(lf) == (1, 129, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
@test pf[end] == 128.0

@info "Test: amp()"
amp_data = NeuroAnalyzer.amp(e10; ch = "all")
@test size(amp_data.peak_amp) == (24, 10)
@test size(amp_data.rms_amp) == (24, 10)
@test size(amp_data.p2p_amp) == (24, 10)
@test size(amp_data.semi_p2p_amp) == (24, 10)
@test size(amp_data.ms_amp) == (24, 10)
@test size(amp_data.te_signal) == (24, 10)

@info "Test: rms()"
@test NeuroAnalyzer._rms(v1) == 3.3166247903554
@test NeuroAnalyzer._rms(a1) == [1.0 1.0; 1.0 1.0]
@test NeuroAnalyzer.rms(e10, ch = "Fp1") ≈
      [
    363.813292624 349.851835483 336.525088954 329.744566991 317.292431492 311.970593552 285.781759746 267.596794688 259.64042809 243.790799956
]

@info "Test: rmse()"
@test NeuroAnalyzer._rmse(v1, v2) == 3.0
@test NeuroAnalyzer._rmse(a1, a2) == [1.0 1.0; 1.0 1.0]
@test length(NeuroAnalyzer.rmse(e10, e10, ch1 = "Fp1", ch2 = "Fp2")) == 10

@info "Test: snr()"
@test NeuroAnalyzer._snr(v1, v2) ≈ 3.01029995663
@test NeuroAnalyzer._snr(v1) ≈ 1.89736659610
@test NeuroAnalyzer.snr2(v1) ≈ 1.20604537831
sn, f = NeuroAnalyzer.snr(e10; ch = "all", type = :rms)
@test size(sn) == (24, 1281)
@test length(f) == 1281
sn, f = NeuroAnalyzer.snr(e10; ch = "all", type = :mean)
@test size(sn) == (24, 1281)
@test length(f) == 1281

@info "Test: spectrogram()"
sp, sf, st = NeuroAnalyzer.spectrogram(e10; ch = "Fp1", method = :stft)
@test size(sp) == (129, 89, 1, 10)
@test length(sf) == 129
@test length(st) == 89
sp, sf, st = NeuroAnalyzer.spectrogram(e10; ch = "Fp1", method = :mt)
@test size(sp) == (129, 89, 1, 10)
@test length(sf) == 129
@test length(st) == 89
sp, sf, st = NeuroAnalyzer.spectrogram(e10; ch = "Fp1", method = :mw)
@test size(sp) == (129, 2560, 1, 10)
@test length(sf) == 129
@test length(st) == 2560
sp, sf, st = NeuroAnalyzer.spectrogram(e10; ch = "Fp1", method = :gh)
@test size(sp) == (129, 2560, 1, 10)
@test length(sf) == 129
@test length(st) == 2560
sp, sf, st = NeuroAnalyzer.spectrogram(e10; ch = "Fp1", method = :cwt)
@test size(sp) == (30, 2560, 1, 10)
@test length(sf) == 30
@test length(st) == 2560

@info "Test: spec_seg()"
spec_data = NeuroAnalyzer.spectrogram(e10; ch = "Fp1")
sp, sf, st = spec_data.p, spec_data.f, spec_data.t
# time points are the window centers (the first one at wlen / 2 = 0.5 s)
sp, sst, t, f = spec_seg(sp, sf, st; ch = 1, t = (1, 2), f = (0, 10))
@test size(sp) == (11, 11, 10)
@test t == (6, 16)
@test f == (1, 11)

@info "Test: stationarity()"
s = NeuroAnalyzer.stationarity(e10; ch = "all", method = :adf)
@test size(s) == (24, 2, 10)
s = NeuroAnalyzer.stationarity(e10; ch = "all", method = :cov)
# 2560 samples / 10 = 256 windows, 255 distances
@test size(s) == (255, 10)
s = NeuroAnalyzer.stationarity(e10; ch = "all", method = :hilbert)
@test size(s) == (24, 2559, 10)
s = NeuroAnalyzer.stationarity(e10; ch = "all", method = :mean)
@test size(s) == (24, 10, 10)
s = NeuroAnalyzer.stationarity(e10; ch = "all", method = :var)
@test size(s) == (24, 10, 10)

@info "Test: channel_stats()"
c = channel_stats(e10)
for idx in eachindex(c)
    @test size(c[idx]) == (24, 10)
end

@info "Test: epoch_stats()"
e = epoch_stats(e10)
for idx in eachindex(e)
    @test length(e[idx]) == 10
end

@info "Test: cpsd()"
pxy, f = cpsd(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1, method = :mt)
@test size(pxy) == (1, 1281, 1)
@test length(f) == 1281
pxy, f = cpsd(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1, method = :fft)
@test size(pxy) == (1, 1281, 1)
@test length(f) == 1281
pxy, f = cpsd(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1, method = :stft)
@test size(pxy) == (1, 129, 1)
@test length(f) == 129

@info "Test: tkeo()"
@test NeuroAnalyzer._tkeo(v1) == [1.0, 1.0, 1.0, 1.0, 1.0]
@test NeuroAnalyzer._tkeo(a1) == zeros(2, 3, 2)
# :der: ẋ² - x·ẍ of sin(ωn) is ω² (central differences)
x = sin.(0.1 .* (1:200))
@test NeuroAnalyzer._tkeo(x; method = :der)[10:190] ≈ fill(0.01, 181) rtol = 0.01
t = tkeo(e10; ch = "all", method = :pow)
@test size(t) == (24, 2560, 10)
t = tkeo(e10; ch = "all", method = :der)
@test size(t) == (24, 2560, 10)
t = tkeo(e10; ch = "all", method = :amp)
@test size(t) == (24, 2560, 10)

@info "Test: total_power()"
tp = total_power(e10; ch = "Fp1")
@test size(tp) == (1, 10)
tp = total_power(e10; ch = "Fp1", method = :welch)
@test size(tp) == (1, 10)
tp = total_power(e10; ch = "Fp1", method = :fft)
@test size(tp) == (1, 10)
tp = total_power(e10; ch = "Fp1", method = :stft)
@test size(tp) == (1, 10)
tp = total_power(e10; ch = "Fp1", method = :mt)
@test size(tp) == (1, 10)
tp = total_power(e10; ch = "Fp1", method = :mw)
@test size(tp) == (1, 10)
tp = total_power(e10; ch = "Fp1", method = :gh)
@test size(tp) == (1, 10)

@info "Test: pacor()"
pac, l = pacor(e10; ch = "all", l = 2)
@test size(pac) == (24, 1025, 10)
@test length(l) == 1025

@info "Test: xcov()"
xc, l = xcov(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2)
@test size(xc) == (1, 513, 1)
@test length(l) == 513
xc, l = xcov(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, biased = false)
@test size(xc) == (1, 513, 1)
@test length(l) == 513
xc, l = xcov(
    e10,
    e10;
    ch1 = "Fp1",
    ch2 = "Fp2",
    ep1 = 1,
    ep2 = 2,
    method = :cov,
    biased = false,
)
@test size(xc) == (1, 513, 1)
@test length(l) == 513
xc, l =
    xcov(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, method = :cov, biased = true)
@test size(xc) == (1, 513, 1)
@test length(l) == 513
xc, l = xcov(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, method = :stat)
@test size(xc) == (1, 513, 1)
@test length(l) == 513

@info "Test: xcor()"
xc, l = xcor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2)
@test size(xc) == (1, 513, 1)
@test length(l) == 513
xc, l = xcor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, biased = false)
@test size(xc) == (1, 513, 1)
@test length(l) == 513
xc, l = xcor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, method = :cov)
@test size(xc) == (1, 513, 1)
@test length(l) == 513
xc, l = xcor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, method = :stat)
@test size(xc) == (1, 513, 1)
@test length(l) == 513

@info "Test: amp_at()"
e = average_epochs(e10)
@test size(amp_at(e, t = 2)) == (19, 11)

@info "Test: avgamp_at()"
@test size(avgamp_at(e, t = (2, 2.5))) == (19, 11)

@info "Test: maxamp_at()"
@test size(maxamp_at(e, t = (2, 2.5))) == (19, 11)

@info "Test: minamp_at()"
@test size(minamp_at(e, t = (2, 2.5))) == (19, 11)

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
xc = NeuroAnalyzer._xcor(x, y; l = 4, demean = false)
l = collect(-4:4)
f = axc2frq(xc[1, :, 1], l)
@test length(f) == 1

@info "Test: hjorth()"
h_act, h_mob, h_comp = NeuroAnalyzer._hjorth(v1)
@test h_act == 2.5
@test h_mob == 0.17320508075688773
@test h_comp == 5.270462766947299
h_act, h_mob, h_comp = NeuroAnalyzer._hjorth(a1)
@test h_act == zeros(2, 2)
h_act, h_mob, h_comp = hjorth(e10; ch = "all")
@test size(h_act) == (24, 10)
@test size(h_mob) == (24, 10)
@test size(h_comp) == (24, 10)

@info "Test: peak_frq()"
pf = peak_frq(e10; ch = "all", flim = (8, 13))
@test size(pf) == (24, 10)

@info "Test: peak_amp()"
pa = peak_amp(e10; ch = "all", flim = (8, 13))
@test size(pa) == (24, 10)

@info "Test: peak_pow()"
pp = peak_pow(e10; ch = "all", flim = (8, 13))
@test size(pp) == (24, 10)

@info "Test: phs()"
ph, f = phs(e10; ch = "all")
@test size(ph) == (24, 1281, 10)
@test length(f) == 1281

@info "Test: asy_idx()"
@test asy_idx(e10, ch1 = "Fp1", ch2 = "Fp1", flim = (0, 10), method = :welch) ==
      (asi = 0.0, nasi = 0.0)
@test asy_idx(e10, ch1 = "Fp1", ch2 = "Fp1", flim = (0, 10), method = :fft) ==
      (asi = 0.0, nasi = 0.0)
@test asy_idx(e10, ch1 = "Fp1", ch2 = "Fp1", flim = (0, 10), method = :mt) ==
      (asi = 0.0, nasi = 0.0)
@test asy_idx(e10, ch1 = "Fp1", ch2 = "Fp1", flim = (0, 10), method = :stft) ==
      (asi = 0.0, nasi = 0.0)
@test asy_idx(e10, ch1 = "Fp1", ch2 = "Fp1", flim = (0, 10), method = :mw) ==
      (asi = 0.0, nasi = 0.0)
@test asy_idx(e10, ch1 = "Fp1", ch2 = "Fp1", flim = (0, 10), method = :gh) ==
      (asi = 0.0, nasi = 0.0)

@info "Test: sym_idx()"
@test NeuroAnalyzer._sym_idx(v) == 5.0
@test sym_idx(e10, ch = "Fp1") ==
      [
    852.3333333333334 1279.0 852.3333333333334 852.3333333333334 1279.0 852.3333333333334 639.0 852.3333333333334 1279.0 1279.0
]

@info "Test: lat_idx()"
@test lat_idx(e10, frq = (1, 3.5)) isa Float64
# right hemisphere with more power: positive index
e_tmp = deepcopy(e10)
e_tmp.data[get_channel(e_tmp; ch = channel_pick(e_tmp; pick = :r)), :, :] .*= 2
@test lat_idx(e_tmp, frq = (1, 3.5)) > 0

@info "Test: vartest()"
f, p = NeuroAnalyzer.vartest(e10; ch = "all")
@test size(f) == (24, 24, 10)
@test size(p) == (24, 24, 10)
f, p = NeuroAnalyzer.vartest(e10, e10; ch1 = "all", ch2 = "all")
@test size(f) == (24, 24, 10)
@test size(p) == (24, 24, 10)

@info "Test: seg_mean()"
@test seg_mean(ones(5, 5, 5)) == ones(5)
@test seg_mean(ones(5, 5, 5), ones(5, 5, 5)) == (seg1 = ones(5), seg2 = ones(5))

@info "Test: spec_flim()"
p = ones(10, 100, 5)
f = collect(1:100)
p2, f2 = spec_flim(p, f; flim = (5, 10))
@test size(p2) == (10, 6, 5)
@test length(f2) == 6
p = ones(100, 200, 10, 5)
f = collect(1:100)
p2, f2 = spec_flim(p, f; flim = (5, 10))
@test size(p2) == (6, 200, 10, 5)
@test length(f2) == 6

@info "Test: spec_tlim()"
p = ones(100, 200, 10, 5)
t = collect(1:200)
p2, t2 = spec_tlim(p, t; seg = (5, 10))
@test size(p2) == (100, 6, 10, 5)
@test length(t2) == 6

@info "Test: seg_extract()"
x = ones(100, 100)
@test seg_extract(x, (10, 10, 20, 20)) == ones(11, 11)
@test seg_extract(x, (10, 10, 20, 20), v = true) == ones(11 * 11)
@test seg_extract(x, (10, 10, 20, 20), c = true) == ones(496)
x = ones(10, 10)
x[4:6, 4:6] .+= 1
@test length(seg_extract(x, threshold = 1, threshold_type = :eq)[1]) == 91
@test length(seg_extract(x, threshold = 1, threshold_type = :neq)[1]) == 9
@test length(seg_extract(x, threshold = 1.1, threshold_type = :leq)[1]) == 91
@test length(seg_extract(x, threshold = 1.1, threshold_type = :geq)[1]) == 9
@test length(seg_extract(x, threshold = 1.1, threshold_type = :g)[1]) == 9
@test length(seg_extract(x, threshold = 1.1, threshold_type = :l)[1]) == 91
@test x[seg_extract(x, threshold = 1, threshold_type = :neq)[1]] == ones(9) .* 2

@info "Test: std()"
s = std(e10)
@test size(s) == (24, 2560)

@info "Test: topo_var()"
erp = average_epochs(e10)
s = topo_var(erp; ch = "eeg")
@test length(s) == 2560

@info "Test: erp_gfp()"
erp = average_epochs(e10)
s = erp_gfp(erp; ch = "all")
@test length(s) == 2560
s = erp_gfp(erp; ch = "all", norm = true)
@test size(s) == (19, 2560)

@info "Test: diss()"
erp = average_epochs(e10)
gd, sc = diss(erp, erp; ch1 = "all", ch2 = "all")
@test gd == zeros(2560)
@test sc == ones(2560)

@info "Test: sumsim()"
@test NeuroAnalyzer._sumsim(v1, v2, theta = 1) ≈ 0.00122085489442
@test NeuroAnalyzer._sumsim(a1, a2, theta = 1) ≈
      [0.176921206317 0.176921206317; 0.176921206317 0.176921206317]
@test sumsim(e10, e10, ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1, theta = 0.0001) ≈
      [0.33211228;;]

@info "Test: hfd()"
@test NeuroAnalyzer._hfd([1.0, 2.0, 3.5, 2.0, 5.0, 11.0, 2.0, 11.0]) ≈ 0.86044864760
@test size(hfd(e10, ch = "all")) == (24, 10)

@info "Test: dirinrg()"
@test NeuroAnalyzer._dirinrg([1.0, 2.0, 3.5, 2.0, 5.0, 11.0, 2.0, 11.0]) == 212.5
@test size(dirinrg(e10, ch = "all")) == (24, 10)

@info "Test: transform()"
c, a, p, ph = NeuroAnalyzer.ftransform(rand(100))
@test length(c) == 51
@test length(a) == 51
@test length(p) == 51
@test length(ph) == 51
c, a, p, ph = NeuroAnalyzer.ftransform(rand(100); nf = true)
@test length(c) == 100
@test length(a) == 100
@test length(p) == 100
@test length(ph) == 100
c, a, p, ph = NeuroAnalyzer.htransform(rand(100))
@test length(c) == 100
@test length(a) == 100
@test length(p) == 100
@test length(ph) == 100
c, a, p, ph = NeuroAnalyzer._transform(rand(10, 100, 10); h = false, nf = false)
@test size(c) == (10, 51, 10)
@test size(a) == (10, 51, 10)
@test size(p) == (10, 51, 10)
@test size(ph) == (10, 51, 10)
c, a, p, ph = NeuroAnalyzer._transform(rand(10, 100, 10); h = false, nf = true)
@test size(c) == (10, 100, 10)
@test size(a) == (10, 100, 10)
@test size(p) == (10, 100, 10)
@test size(ph) == (10, 100, 10)
c, a, p, ph = NeuroAnalyzer._transform(rand(10, 100, 10); h = true)
@test size(c) == (10, 100, 10)
@test size(a) == (10, 100, 10)
@test size(p) == (10, 100, 10)
@test size(ph) == (10, 100, 10)
c, a, p, ph = NeuroAnalyzer.transform(e10; ch = "all")
@test size(c) == (24, 1281, 10)
@test size(a) == (24, 1281, 10)
@test size(p) == (24, 1281, 10)
@test size(ph) == (24, 1281, 10)
c, a, p, ph = NeuroAnalyzer.transform(e10; ch = "all", nf = true)
@test size(c) == (24, 2560, 10)
@test size(a) == (24, 2560, 10)
@test size(p) == (24, 2560, 10)
@test size(ph) == (24, 2560, 10)
c, a, p, ph = NeuroAnalyzer.transform(e10; ch = "all", h = true)
@test size(c) == (24, 2560, 10)
@test size(a) == (24, 2560, 10)
@test size(p) == (24, 2560, 10)
@test size(ph) == (24, 2560, 10)

@info "Test: zipratio()"
e10_tmp = deepcopy(e10)
z1 = zipratio(e10_tmp)
e10_tmp.data = rand(size(e10.data, 1), size(e10.data, 2), size(e10.data, 3))
z2 = zipratio(e10_tmp)
e10_tmp.data = zeros(size(e10.data))
z3 = zipratio(e10_tmp)
@test z1 > z3
@test z2 > z3

@info "Test: hhtspectrogram()"
imf = emd(e10; ch = "Fp1", ep = 1)[1:(end - 1), :]
hht_data = hhtspectrogram(imf, e10.epoch_time; fs = sr(e10))
@test size(hht_data.p) == (257, 2560)
@test length(hht_data.f) == 257
@test length(hht_data.t) == 2560

@info "Test: hmspectrum()"
hms_data = hmspectrum(e10; ch = "Fp1")
@test size(hms_data.p) == (257, 10)
@test length(hms_data.f) == 257

@info "Test: ghexp()"
g = ghexp(e10; ch = "Fp1", tau_range = 1:10)
@test size(g) == (1, 1, 2, 10)
g = ghexp(e10; ch = "Fp1", tau_range = 1:10, q_range = 0.1:0.1:1.0)
@test size(g) == (1, 10, 2, 10)

@info "Test: wpli()"
pv, sd, phd, s1p, s2p = wpli(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test size(pv) == (1, 1)
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)
pv, sd, phd, s1p, s2p =
    wpli(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1, debiased = true)
@test size(pv) == (1, 1)
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: dpli()"
pv, sd, phd, s1p, s2p = dpli(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test size(pv) == (1, 1)
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: iplv()"
pv, sd, phd, s1p, s2p = iplv(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test size(pv) == (1, 1)
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: aecor()"
aec = aecor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test size(aec) == (1, 1)

@info "Test: escor()"
esc = escor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test size(esc) == (1, 1)

@info "Test: cosim()"
cs = cosim(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test size(cs) == (1, 1)

@info "Test: psa()"
ps = psa(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test size(ps) == (1, 1)

@info "Test: corr()"
cr = corr(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test size(cr) == (1, 1)

@info "Test: psi()"
pv = psi(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test pv[1][1] == -pv[1][2]

@info "Test: fooof()"
# reference values from fooof 1.1.0 (Python) on the same noise-free spectra:
#   fixed: ap = [1.00846, 1.50140], peaks = [10.0016 0.78872 2.91857; 22.0024 0.38969 3.80960], R² = 0.99986
#   knee:  ap = [1.0, 25.0, 2.0], no peaks
f = collect(1:0.5:40)
p = 10 .^ (1 .- 1.5 .* log10.(f) .+ 0.8 .* exp.(-(f .- 10) .^ 2 ./ (2 * 1.5^2)) .+ 0.4 .* exp.(-(f .- 22) .^ 2 ./ (2 * 2.0^2)))
r = NeuroAnalyzer._fooof(f, p)
@test isapprox(r.ap, [1.00846, 0.0, 1.50140]; atol = 1e-3)
@test size(r.pe) == (2, 3)
@test isapprox(r.pe, [10.0016 0.78872 2.91857; 22.0024 0.38969 3.80960]; atol = 1e-3)
@test isapprox(r.r2, 0.99986; atol = 1e-4)
@test length(r.model) == length(r.f) == length(f)
pk = 10 .^ (1 .- log10.(25 .+ f .^ 2))
r = NeuroAnalyzer._fooof(f, pk; aperiodic_mode = :knee)
@test isapprox(r.ap, [1.0, 25.0, 2.0]; rtol = 1e-3)
@test size(r.pe) == (0, 3)
r = NeuroAnalyzer._fooof(f, p; flim = (2, 30), max_n_peaks = 1)
@test r.f[1] == 2.0
@test r.f[end] == 30.0
@test size(r.pe, 1) == 1
@test_throws ArgumentError NeuroAnalyzer._fooof(f, DSP.pow2db.(p))
@test_throws ArgumentError NeuroAnalyzer._fooof(vcat(0.0, f), vcat(1.0, p); flim = (0, 40))
@test_throws ArgumentError NeuroAnalyzer._fooof(f, p; aperiodic_mode = :foo)
r = NeuroAnalyzer._fooof(rand(10, 1000, 10); fs = 100)
@test size(r.ap) == (10, 3, 10)
@test size(r.pe) == (10, 10)
@test size(r.r2) == (10, 10)
@test r.f[1] == 1.0
@test r.f[end] == 50.0
r = NeuroAnalyzer._fooof(rand(10, 1000, 10); fs = 100, method = :stft)
@test size(r.ap) == (10, 3, 10)
@test r.f[1] == 1.0
# non-integer resolution: psd() rounds frequencies to 3 decimals
r = NeuroAnalyzer._fooof(rand(2, 2500, 2); fs = 250, wlen = 256)
@test size(r.ap) == (2, 3, 2)
# flat channel: NaN result instead of an error
s = rand(3, 1000, 2)
s[2, :, :] .= 0
r = NeuroAnalyzer._fooof(s; fs = 100)
@test all(isnan, r.ap[2, :, :])
@test !any(isnan, r.ap[1, :, :])
r = fooof(e10; ch = "Fp1")
@test size(r.ap) == (1, 3, 10)
@test size(r.model) == (1, length(r.f), 10)

@info "Test: regressions"
# GFP: population standard deviation; columns with GFP = 0 are left unchanged
m_tmp = [1.0 2.0 5.0; 3.0 2.0 1.0]
@test NeuroAnalyzer._erp_gfp(m_tmp) ≈ [1.0, 0.0, 2.0]
@test erp_gfp_norm(m_tmp)[:, 2] == m_tmp[:, 2]
# ftt_analyze: intervals within trials only
t_tmp = (taps = [3, 1], tap_t = [[0.0, 100.0, 300.0], [50.0]], tap_d = [[1.0, 1.0, 1.0], [1.0]],
    taps_int = [0, 0], tap_t_int = [Float64[], Float64[]], tap_d_int = [Float64[], Float64[]])
r_tmp = ftt_analyze(t_tmp)
@test r_tmp.t_mean == 150.0 && r_tmp.t_rmssd == 100.0 && isnan(r_tmp.t_sdsd)
# cosim within one object (used an undefined variable)
@test size(cosim(e10; ch = ["Fp1", "Fp2"])) == (2, 2, 10)
# dPLI within one object: directed matrix
d_tmp = dpli(e10; ch = ["Fp1", "Fp2", "F3"])
@test d_tmp[1, 2, 1] ≈ 1 - d_tmp[2, 1, 1] && d_tmp[1, 1, 1] == 0.5
# dPLI: s1 leads s2 by a quarter cycle
x1 = sin.(2π * 8 .* (0:511) ./ 256)
x2 = sin.(2π * 8 .* (0:511) ./ 256 .- π / 2)
@test NeuroAnalyzer._dpli(x1, x2).pv > 0.9
# pacor() with :reg (always threw)
@test size(pacor(e10; ch = "Fp1", l = 0.05, method = :reg).pac, 2) == 2 * round(Int, 0.05 * sr(e10)) + 1
# phase differences are wrapped to (-π, π]
x1 = sin.(2π * 5 .* (0:255) ./ 256)
x2 = sin.(2π * 5 .* (0:255) ./ 256 .+ 3.0)
@test all(abs.(NeuroAnalyzer._phdiff(x1, x2; h = true)) .<= π)
# DISS: identical maps (up to scaling and offset) → 0, inverted maps → 2
m_tmp = rand(10, 20)
d_tmp = NeuroAnalyzer._diss(m_tmp, 3 .* m_tmp .+ 1)
@test all(d_tmp.gd .< 1e-10) && all(d_tmp.sc .≈ 1)
d_tmp = NeuroAnalyzer._diss(m_tmp, -m_tmp)
@test all(d_tmp.gd .≈ 2) && all(d_tmp.sc .≈ -1)
# SEF of a sine is its frequency, independent of the frequency resolution
x = sin.(2π * 10 .* (0:2559) ./ 256)
@test abs(NeuroAnalyzer._sef(x; fs = 256, wlen = 512) - 10) <= 1
# xcor: a signal with itself is 1 at lag 0 (:sum and :cov)
x = rand(200)
@test NeuroAnalyzer._xcor(x, x; l = 5)[1, 6, 1] ≈ 1
@test NeuroAnalyzer._xcor(x, x; l = 5, method = :cov)[1, 6, 1] ≈ 1
# aecor()/escor() within one object (used an undefined variable); ESC is not symmetric
@test size(aecor(e10; ch = ["Fp1", "Fp2"])) == (2, 2, 10)
@test size(escor(e10; ch = ["Fp1", "Fp2"])) == (2, 2, 10)
# Granger causality: x drives y (lag 1)
rng_tmp = NeuroAnalyzer.Random.MersenneTwister(1)
x_gc = zeros(2, 2000, 1)
for t = 2:2000
    x_gc[1, t, 1] = 0.5 * x_gc[1, t - 1, 1] + randn(rng_tmp)
    x_gc[2, t, 1] = 0.3 * x_gc[2, t - 1, 1] + 0.6 * x_gc[1, t - 1, 1] + randn(rng_tmp)
end
g_tmp = NeuroAnalyzer._granger(x_gc; order = 1)
@test g_tmp.gc[1, 2] > 0.1 && g_tmp.gc[2, 1] < 0.01 && g_tmp.p[1, 2] < 0.001
# spectral GC: 2 channels, conditional = pairwise; frequency average = time-domain GC
gs_p = NeuroAnalyzer._granger_spectral(x_gc; fs = 100, order = 1, nfrq = 128)
gs_c = NeuroAnalyzer._granger_spectral(x_gc; fs = 100, order = 1, nfrq = 128, conditional = true)
@test maximum(abs.(gs_p.gc .- gs_c.gc)) < 1e-6
gm = sum((gs_p.gc[1, 2, 1:(end - 1)] .+ gs_p.gc[1, 2, 2:end]) ./ 2) / (length(gs_p.frq) - 1)
@test gm ≈ g_tmp.gc[1, 2] rtol = 0.05
# NEURO methods
@test size(granger(e10; ch = ["Fp1", "Fp2"], order = 2).gc) == (2, 2)
@test size(granger_spectral(e10; ch = ["Fp1", "Fp2"], order = 2, nfrq = 16).gc) == (2, 2, 16)
# psd_slope(): slope per Hz (0.5 Hz resolution)
x = randn(NeuroAnalyzer.Random.MersenneTwister(2), 5120)
r_tmp = NeuroAnalyzer._psd_slope(x; fs = 256, wlen = 512)
@test r_tmp.ls ≈ (r_tmp.lf[end] - r_tmp.lf[1]) / (r_tmp.f[end] - r_tmp.f[1])
# cph(): antisymmetric cross-phase matrix
c_tmp = cph(e10; ch = ["Fp1", "Fp2"])
@test c_tmp.ph[1, 2, :, 1] ≈ -c_tmp.ph[2, 1, :, 1]
# erp_auc(): selected channels (not the first ones); positive + negative = total
e_erp = average_epochs(e10)
a_all = erp_auc(e_erp; ch = "Fp2")
@test a_all ≈ erp_auc(e_erp; ch = "Fp2", type = :pos) .+ erp_auc(e_erp; ch = "Fp2", type = :neg)
@test a_all != erp_auc(e_erp; ch = "Fp1")
# negentropy of a Gaussian signal ≈ 0 (bits vs bits); the object is not modified
x = randn(NeuroAnalyzer.Random.MersenneTwister(3), 20000)
@test abs(NeuroAnalyzer._negentropy(x; norm = false)) < 0.05
# negentropy does not depend on the signal scale (differential entropy over the whole KDE support)
xu = rand(NeuroAnalyzer.Random.MersenneTwister(4), 10000)
@test NeuroAnalyzer._negentropy(xu; norm = false) > 0
@test NeuroAnalyzer._negentropy(xu; norm = false) ≈ NeuroAnalyzer._negentropy(1e-3 .* xu; norm = false) rtol = 1e-4
@test NeuroAnalyzer._negentropy(xu) ≈ NeuroAnalyzer._negentropy(xu; norm = false) rtol = 1e-4
e_tmp = deepcopy(e10)
negentropy(e_tmp; ch = "Fp1", demean = false)
@test e_tmp.data == e10.data
# msci95() bootstrap: SE ≈ analytical SE; few resamples do not fail
x = randn(NeuroAnalyzer.Random.MersenneTwister(4), 400)
r_n = NeuroAnalyzer._msci95(x)
r_b = NeuroAnalyzer._msci95(x; method = :boot, n = 20)
@test r_b.se ≈ r_n.se rtol = 0.2
@test NeuroAnalyzer._msci95(rand(5); method = :boot, n = 1).ll isa Float64
# ftransform(): one-sided amplitude in the last bin of an odd-length FFT (not Nyquist)
x = 3 .* cos.(2π * 50 .* (0:100) ./ 101)
@test maximum(ftransform(x).a) ≈ 3
# spectrograms: Gaussian-Hilbert power of a stationary sine = squared amplitude, flat in time
x = 2 .* cos.(2π * 10 .* (0:1023) ./ 128)
g_tmp = ghtspectrogram(x; fs = 128, db = false)
fi = argmin(abs.(g_tmp.f .- 10))
@test g_tmp.p[fi, 300] ≈ 4 rtol = 0.05
@test g_tmp.p[fi, 100] ≈ g_tmp.p[fi, 500] rtol = 0.05
# STFT time axis: centers of the segments (not 0 … end)
s_tmp = NeuroAnalyzer._spectrogram(x; fs = 128, wlen = 128, woverlap = 64, db = false)
@test s_tmp.t[1] ≈ 0.5 rtol = 0.02
@test length(s_tmp.t) == size(s_tmp.p, 2)
# Hilbert spectrum of IMFs (the matrix method does not decompose again)
h_tmp = hhtspectrogram(reshape(x, 1, :), collect(0:1023) ./ 128; fs = 128, db = false)
@test argmax(vec(sum(h_tmp.p; dims = 2))) == argmin(abs.(h_tmp.f .- 10))
# envelope CI: 95% CI of the mean across channels at each point; median CI: lower ≤ upper
a_tmp = randn(NeuroAnalyzer.Random.MersenneTwister(5), 40, 3, 2)
em, el, eu = NeuroAnalyzer._env_mean_ci(a_tmp, 1)
@test em[2, 1] ≈ NeuroAnalyzer.Statistics.mean(a_tmp[:, 2, 1])
@test eu[2, 1] - em[2, 1] ≈ 1.96 * NeuroAnalyzer.Statistics.std(a_tmp[:, 2, 1]) / sqrt(40)
em, el, eu = NeuroAnalyzer._env_median_ci(a_tmp, 3)
@test size(em) == (3, 1) && all(el .<= em .<= eu)
# lat_idx() / channel_pick() for an ERP object
@test lat_idx(average_epochs(e10); frq = (8, 12)) isa Float64
# band_power() 3-D method default demean
@test size(NeuroAnalyzer._band_power(rand(2, 512, 1); fs = 256, flim = (1, 10))) == (2, 1)
# PSI: x leads y → PSI(x, y) > 0; a signal shorter than 1 s uses its own frequency bins
x = randn(NeuroAnalyzer.Random.MersenneTwister(1), 256 * 20)
y = circshift(x, 5) .+ 0.5 .* randn(NeuroAnalyzer.Random.MersenneTwister(2), length(x))
p = NeuroAnalyzer._psi(x, y; fs = 256, flim = (1, 40))
@test p[1] > 0 && p[2] ≈ -p[1]
@test NeuroAnalyzer._psi(x[1:200], y[1:200]; fs = 256, flim = (2, 40))[1] > 0
@test_throws ArgumentError NeuroAnalyzer._psi(x[1:8], y[1:8]; fs = 256, flim = (1, 2))
# vendored Simpson / FIRLS: values from SciPy 1.17 (scipy.integrate.simpson, scipy.signal.firls)
@test NeuroAnalyzer.Simpson.simpson(collect(0:0.1:2) .^ 3; dx = 0.1) ≈ 4.0
h = NeuroAnalyzer.FIRLSFilterDesign.firls_design(100, [0, 9, 11, 128.0], [1, 1, 0, 0.0], ones(4), false; fs = 256)
@test length(h) == 101 && h ≈ reverse(h)
@test sum(h) ≈ 0.9627675365974948
@test h[51] ≈ 0.078125
@test h[1] ≈ -0.0014178653128872954


##
## values: SciPy 1.17 references, closed-form results on synthetic signals
##

fs = 256
tv = (0:1023) ./ fs
xa = sin.(2π * 10 .* tv) .+ 0.5 .* sin.(2π * 23 .* tv .+ 0.3) .+ 0.2 .* cos.(2π * 47 .* tv) .+ 0.1 .* sin.(2π * 3.7 .* tv) .^ 3
ya = 0.8 .* sin.(2π * 10 .* tv .- 0.7) .+ 0.3 .* sin.(2π * 31 .* tv) .+ 0.1 .* cos.(2π * 47 .* tv .+ 1) .+ 0.05 .* cos.(2π * 71 .* tv) .* sin.(2π * 0.9 .* tv)
s2 = 2 .* sin.(2π * 10 .* tv)

@info "Test: psd() / spectrogram() / coherence() / htransform() vs SciPy"
# scipy.signal.welch(x - mean, window = hann(256, sym = True), nperseg = 256, noverlap = 230, detrend = False)
p = NeuroAnalyzer._psd(xa; fs = fs, method = :welch, wlen = 256, woverlap = 230)
@test length(p.p) == 129
@test p.p[[11, 24, 48]] ≈ [0.3323734482782393, 0.08300780760898713, 0.013281255151670541] rtol = 1e-10
# scipy.signal.periodogram(x - mean, window = hann(1024, sym = True), detrend = False)
p = NeuroAnalyzer._psd(xa; fs = fs, method = :fft)
@test p.f[[41, 93]] == [10.0, 23.0]
@test p.p[[41, 93]] ≈ [1.3321094937218387, 0.33300812959494636] rtol = 1e-10
# scipy.signal.spectrogram(..., mode = "psd"): 30 segments, centers 0.5, 0.6015625, … s
sp = NeuroAnalyzer._spectrogram(xa; fs = fs, method = :stft, wlen = 256, woverlap = 230, db = false)
@test size(sp.p) == (129, 30)
@test sp.t[1:2] == [0.5, 0.6015625]
@test sp.p[11, 1] ≈ 0.3388363680892147 rtol = 1e-6
@test sp.p[24, 6] ≈ 0.08300121883523325 rtol = 1e-6
# scipy.signal.coherence(x, y, window = hann(256, sym = True), nperseg = 256, noverlap = 192, nfft = 256)
c = NeuroAnalyzer._coherence(xa, ya; fs = fs, method = :stft, wlen = 256, wstep = 64, demean = true)
@test c.f[[11, 32, 48, 61]] == [10.0, 31.0, 47.0, 60.0]
@test c.msc[[11, 32, 48, 61]] ≈ [0.9998861347266094, 0.9611892034248443, 0.9999999991801273, 0.2826798537576459] rtol = 1e-9
@test c.imcoh ≈ imag.(c.coh)
# a single segment (:fft): |coherence| = 1 at every frequency
@test all(NeuroAnalyzer._coherence(xa, ya; fs = fs, method = :fft).msc .≈ 1)
# |scipy.signal.hilbert(x)|
@test htransform(xa).a[[501, 101]] ≈ [0.3079868348251104, 1.249903120202579] rtol = 1e-10
@test NeuroAnalyzer._hanalytic(xa) ≈ NeuroAnalyzer.DSP.hilbert(xa)
@test NeuroAnalyzer._hanalytic(reshape(xa, 1, :, 1))[1, :, 1] ≈ NeuroAnalyzer.DSP.hilbert(xa)
@test hanalytic(e10; ch = "Fp1")[1, :, 2] ≈ NeuroAnalyzer.DSP.hilbert(e10.data[1, :, 2])

@info "Test: spectral estimates of a sine (closed form)"
# Parseval: the PSD integrates to the variance (Welch, periodogram, multitaper)
for m in (:welch, :fft, :mt)
    p = NeuroAnalyzer._psd(s2; fs = fs, method = m)
    @test sum(p.p) * (p.f[2] - p.f[1]) ≈ 2.0 rtol = 1e-3
    @test p.f[argmax(p.p)] == 10.0
    @test NeuroAnalyzer._total_power(s2; fs = fs, method = m) ≈ 2.0 rtol = 1e-3
end
@test NeuroAnalyzer._total_power(xa; fs = fs) ≈ sum(abs2, xa .- sum(xa) / length(xa)) / length(xa) rtol = 1e-3
# Morlet wavelets have unit gain at their center frequency: power of A sin = (A/2)²
mw = mwspectrogram(s2; fs = fs, db = false)
@test mw.f[argmax(vec(sum(mw.p; dims = 2)))] == 10.0
@test mw.p[findfirst(==(10.0), mw.f), 512] ≈ 1.0 rtol = 1e-2
r = mwpsd(s2; fs = fs)
@test r.f[argmax(r.p)] == 10.0
@test r.f[1] == 1.0
# Gaussian filter + Hilbert: envelope² of A sin = A²
g = ghtspectrogram(s2; fs = fs, db = false)
@test g.p[findfirst(==(10.0), g.f), 512] ≈ 4.0 rtol = 1e-6
r = ghpsd(s2; fs = fs)
@test r.f[argmax(r.p)] == 10.0
cw = cwtspectrogram(s2; fs = fs)
@test abs(cw.f[argmax(vec(sum(cw.m; dims = 2)))] - 10) < 1
@test size(cw.m) == (length(cw.f), 1024) && issorted(cw.f)
@test NeuroAnalyzer._pow2db_finite([0.0, 1.0, 10.0]) == [0.0, 0.0, 10.0]
@test NeuroAnalyzer._pow2db_finite([0.0, 0.0]) == [-Inf, -Inf]

@info "Test: amp() / band_mpower() / peak_*() / phs() / h_mob / sef values"
a = NeuroAnalyzer._amp(s2)
@test a.peak_amp ≈ 2.0 && a.p2p_amp ≈ 4.0 && a.semi_p2p_amp ≈ 2.0
@test a.rms_amp ≈ sqrt(2) && a.ms_amp ≈ 2.0 && a.te_signal ≈ 2048.0
bm = NeuroAnalyzer._band_mpower(s2; fs = fs, flim = (8, 12))
pw = NeuroAnalyzer._psd(s2; fs = fs)
@test bm.maxfrq == 10.0
@test bm.maxbp == maximum(pw.p)
@test bm.maxba == sqrt(bm.maxbp)
@test bm.mbp ≈ sum(pw.p[9:13]) / 5
@test NeuroAnalyzer._peak_frq(s2; fs = fs, flim = (5, 15)) == 10.0
@test NeuroAnalyzer._peak_frq(xa; fs = fs, flim = (15, 30)) == 23.0
@test NeuroAnalyzer._peak_pow(s2; fs = fs, flim = (5, 15)) == maximum(pw.p)
@test NeuroAnalyzer._peak_amp(s2; fs = fs, flim = (5, 15)) == sqrt(maximum(pw.p))
# phase spectrum: cos → 0, sin → -π/2 at its frequency
@test NeuroAnalyzer._phs(cos.(2π * 10 .* tv); fs = fs).ph[41] ≈ 0 atol = 1e-10
@test NeuroAnalyzer._phs(sin.(2π * 10 .* tv); fs = fs).ph[41] ≈ -π / 2
@test NeuroAnalyzer._phs(cos.(2π * 10 .* tv); fs = fs).f[41] == 10.0
# Hjorth mobility of a sine (central differences): ≈ sin(2πf / fs)
@test NeuroAnalyzer._h_mob(s2) ≈ sin(2π * 10 / fs) rtol = 0.02
# spectral edge: 40% of the power of two equal sines (10, 40 Hz) lies at or below 10 Hz; 95% below ~40 Hz
x2 = sin.(2π * 10 .* tv) .+ sin.(2π * 40 .* tv)
@test NeuroAnalyzer._sef(x2; fs = fs, x = 0.4) in (10.0, 11.0)
@test NeuroAnalyzer._sef(x2; fs = fs, x = 0.95) in (40.0, 41.0)
@test size(sef(e10; ch = "Fp1")) == (1, 10)

@info "Test: connectivity values (known phase lag)"
x1 = sin.(2π * 10 .* tv)
x2 = sin.(2π * 10 .* tv .- π / 4)
@test NeuroAnalyzer._iplv(x1, x2).ipl ≈ sin(π / 4)
@test NeuroAnalyzer._psa(x1, x2) ≈ cos(π / 4)
@test NeuroAnalyzer._plv(x1, x2).pv ≈ 1.0
@test NeuroAnalyzer._pli(x1, x2).pv ≈ 1.0
@test NeuroAnalyzer._cosim(x1, x1) ≈ 1.0
@test NeuroAnalyzer._cosim(x1, -x1) ≈ -1.0
@test NeuroAnalyzer._cosim(x1, cos.(2π * 10 .* tv)) ≈ 0 atol = 1e-10
@test_throws ArgumentError NeuroAnalyzer._cosim(x1, x1[1:10])
# amplitude envelope correlation: the same 1 Hz modulation of different carriers
mod1 = 1 .+ 0.5 .* sin.(2π * 1 .* tv)
@test NeuroAnalyzer._aecor(mod1 .* sin.(2π * 20 .* tv), mod1 .* cos.(2π * 30 .* tv)) > 0.99
# envelope-signal correlation: raw signal vs the envelope of an AM signal
@test NeuroAnalyzer._escor(mod1, mod1 .* sin.(2π * 30 .* tv)) > 0.99
@test_throws ArgumentError NeuroAnalyzer._aecor(x1, x1[1:10])
# ITPC: identical epochs → 1
ep3 = repeat(reshape(xa, 1, :, 1), 1, 1, 5)
r = NeuroAnalyzer._itpc_spec(ep3)
@test all(r.itpcv .≈ 1)
@test r.itpcz ≈ 5 .* r.itpcv .^ 2
@test_throws ArgumentError NeuroAnalyzer._itpc_spec(repeat(ep3, 2, 1, 1))

@info "Test: xcov() / pacor() values"
xc = NeuroAnalyzer._xcov(xa, xa; l = 5)
# biased auto-covariance: lag 0 = variance (1/n), symmetric
@test xc[1, 6, 1] ≈ sum(abs2, xa .- sum(xa) / 1024) / 1024
@test xc[1, 1:5, 1] ≈ reverse(xc[1, 7:11, 1])
# a delayed copy: the cross-covariance peaks at the delay
xd = vcat(zeros(3), xa[1:(end - 3)])
xc = NeuroAnalyzer._xcov(xd, xa; l = 6)
@test argmax(xc[1, :, 1]) - 7 == 3
@test NeuroAnalyzer._xcov(xa, ya; l = 4, method = :cov)[1, 5, 1] ≈ NeuroAnalyzer.Statistics.cov(xa, ya; corrected = false)
@test_throws ArgumentError NeuroAnalyzer._xcov(xa, ya[1:10])
# AR(1) with φ = 0.7: partial autocorrelation 0.7 at lag 1, ≈ 0 beyond
rng_ = NeuroAnalyzer.Random.MersenneTwister(11)
ar = zeros(20000)
for i = 2:20000
    ar[i] = 0.7 * ar[i - 1] + randn(rng_)
end
pc = NeuroAnalyzer._pacor(ar; l = 3)
@test pc[1, 5, 1] ≈ 0.7 atol = 0.02
@test abs(pc[1, 6, 1]) < 0.03
@test NeuroAnalyzer._pacor(ar; l = 3, method = :reg)[1, 5, 1] ≈ 0.7 atol = 0.02

@info "Test: ghexp() values"
# Brownian motion: Hurst exponent ≈ 0.5
bm_ = cumsum(randn(NeuroAnalyzer.Random.MersenneTwister(2), 20000))
h = NeuroAnalyzer._ghexp(bm_; tau_range = 1:20)
@test h[1] ≈ 0.5 atol = 0.05

@info "Test: stationarity_*() values"
s = vcat(fill(1.0, 10), fill(3.0, 10), [5.0])
@test stationarity_mean(s; window = 2) == [1.0, 3.0]
@test stationarity_var(vcat(collect(1.0:4), collect(1.0:4) .* 2); window = 2) ≈ [NeuroAnalyzer.Statistics.var(1:4), NeuroAnalyzer.Statistics.var(2:2:8)]
@test_throws ArgumentError stationarity_mean(s; window = 0)
@test_throws ArgumentError stationarity_var(s; window = 100)
# a sine: constant phase increments of 2πf / fs
@test all(stationarity_hilbert(s2)[100:900] .≈ 2π * 10 / fs)

@info "Test: hrv_detect() / hrv_analyze() values"
e_ecg = deepcopy(keep_epoch(e10; ep = 1))
ecg = findfirst(==("ecg"), e_ecg.header.recording[:channel_type])
e_ecg.data[ecg, :, 1] .= 0.0
e_ecg.data[ecg, 1:256:end, 1] .= 1000.0
h = hrv_detect(e_ecg)
@test h.r_idx == collect(257:256:2305)
@test all(h.nn_seg .== 1000.0)
ha = hrv_analyze([800.0, 850.0, 900.0, 1000.0])
@test ha.menn == 887.5 && ha.mdnn == 875.0
@test ha.vnn ≈ NeuroAnalyzer.Statistics.var([800.0, 850.0, 900.0, 1000.0]) && ha.sdnn ≈ NeuroAnalyzer.Statistics.std([800.0, 850.0, 900.0, 1000.0])
@test ha.rmssd ≈ sqrt((50^2 + 50^2 + 100^2) / 3)
@test ha.sdsd ≈ NeuroAnalyzer.Statistics.std([50.0, 50.0, 100.0])
@test ha.nn50 == 1
@test ha.pnn50 ≈ 1 / 4
e_tmp = keep_channel(e10; ch = "Fp1")
@test_throws ArgumentError hrv_detect(e_tmp)

@info "Test: tpt_detect() / tpt_analyze()"
tp = create_object(; data_type = "tpt")
tp.header.recording[:sampling_rate] = 50
add_channel!(tp; data = zeros(6, 1000, 1), label = ["x", "y", "z", "p1", "p2", "p3"], type = repeat(["other"], 6), unit = repeat([""], 6))
tp.time_pts, tp.epoch_time = NeuroAnalyzer._get_t(tp)
for k in 60:100:960
    tp.data[4, k:(k + 2), 1] .= 100.0
end
pi_ = tpt_detect(tp)
@test length(pi_) >= 1
ta = tpt_analyze(tp)
@test isnothing(ta) || ta.n == length(pi_)
@test_throws ArgumentError tpt_detect(e10)

@info "Test: mep_peaks() values"
mep = import_duomag(joinpath(testfiles_path, "mep", "mep-duomag.ascii"))
pk = mep_peaks(mep)
ss = mep.header.recording[:stimulation_sample] .+ 20
for ch in 1:nchannels(mep)
    @test mep.data[ch, pk[ch, 1], 1] == maximum(mep.data[ch, ss[ch]:end, 1])
    @test mep.data[ch, pk[ch, 2], 1] == minimum(mep.data[ch, ss[ch]:end, 1])
    @test all(pk[ch, :] .>= ss[ch])
end
@test_throws ArgumentError mep_peaks(e10)

@info "Test: envelope CI helper"
a3 = reshape(collect(1.0:24), 2, 3, 4)
em, el, eu = NeuroAnalyzer._env_ci(a3, 1, false)
@test em[:, 1] ≈ vec(sum(a3[:, :, 1]; dims = 1)) ./ 2
@test eu .- em ≈ em .- el
@test (eu .- em)[1, 1] ≈ 1.96 * NeuroAnalyzer.Statistics.std(a3[:, 1, 1]) / sqrt(2)
em, el, eu = NeuroAnalyzer._env_ci(a3, 1, true)
@test em[1, 1] == NeuroAnalyzer.Statistics.median(a3[:, 1, 1])
@test_throws ArgumentError NeuroAnalyzer._env_ci(a3, 4, false)

@info "Test: cpsd() values"
r = NeuroAnalyzer._cpsd(xa, xa; fs = fs, method = :stft, wlen = 256, wstep = 64)
@test all(abs.(imag.(r.pxy)) .< 1e-12)
r2 = NeuroAnalyzer._cpsd(ya, xa; fs = fs, method = :stft, wlen = 256, wstep = 64)
r1 = NeuroAnalyzer._cpsd(xa, ya; fs = fs, method = :stft, wlen = 256, wstep = 64)
@test r2.pxy ≈ conj.(r1.pxy)
@test_throws ArgumentError NeuroAnalyzer._cpsd(xa, ya; fs = fs, method = :x)

@info "Test: FOOOF internals"
fq = collect(1.0:0.5:40.0)
@test NeuroAnalyzer._log10nan(100.0) == 2.0
@test isnan(NeuroAnalyzer._log10nan(0.0))
@test NeuroAnalyzer._fooof_fres(fq) == 0.5
@test NeuroAnalyzer._fooof_ap(fq, [1.0, 2.0]) ≈ 1 .- 2 .* log10.(fq)
@test NeuroAnalyzer._fooof_ap(fq, [1.0, 5.0, 2.0]) ≈ 1 .- log10.(5 .+ fq .^ 2)
@test NeuroAnalyzer._fooof_gauss(fq, [10.0, 1.0, 2.0])[findfirst(==(10.0), fq)] == 1.0
# analytic Jacobians = finite differences
for (fun, jac, p0) in ((NeuroAnalyzer._fooof_ap, NeuroAnalyzer._fooof_ap_jac, [1.0, 2.0]),
                       (NeuroAnalyzer._fooof_ap, NeuroAnalyzer._fooof_ap_jac, [1.0, 5.0, 2.0]),
                       (NeuroAnalyzer._fooof_gauss, NeuroAnalyzer._fooof_gauss_jac, [10.0, 1.0, 2.0, 25.0, 0.5, 3.0]))
    J = jac(fq, p0)
    for k in eachindex(p0)
        d = zeros(length(p0)); d[k] = 1e-6
        @test J[:, k] ≈ (fun(fq, p0 .+ d) .- fun(fq, p0 .- d)) ./ 2e-6 rtol = 1e-5 atol = 1e-8
    end
end
# noise-free aperiodic spectrum: exact recovery of offset and exponent
pap = 10 .^ NeuroAnalyzer._fooof_ap(fq, [2.0, 1.5])
@test NeuroAnalyzer._fooof_ap_fit(fq, log10.(pap), NeuroAnalyzer._fooof_ap_guess(fq, log10.(pap), log10.(pap), false)) ≈ [2.0, 1.5] atol = 1e-6
@test NeuroAnalyzer._fooof_robust_ap_fit(fq, log10.(pap), false) ≈ [2.0, 1.5] atol = 1e-6
@test NeuroAnalyzer._fooof_lm(NeuroAnalyzer._fooof_ap, NeuroAnalyzer._fooof_ap_jac, fq, log10.(pap), [0.0, 0.0], fill(-Inf, 2), fill(Inf, 2)) ≈ [2.0, 1.5] atol = 1e-6
# a Gaussian peak on a flat spectrum: center, height and width
pks = NeuroAnalyzer._fooof_fit_peaks(fq, NeuroAnalyzer._fooof_gauss(fq, [10.0, 0.8, 1.5]), (0.5, 12.0), 5, 0.0, 2.0)
@test size(pks, 1) == 1
@test pks[1, :] ≈ [10.0, 0.8, 1.5] rtol = 1e-3
r = NeuroAnalyzer._fooof_fit(fq, pap .* 10 .^ NeuroAnalyzer._fooof_gauss(fq, [10.0, 0.8, 1.5]);
                             peak_width_limits = (0.5, 12.0), max_n_peaks = 5, min_peak_height = 0.0,
                             peak_threshold = 2.0, aperiodic_mode = :fixed)
# fixed model: ap = [offset, 0 (no knee), exponent]
@test r.ap[[1, 3]] ≈ [2.0, 1.5] atol = 0.01
@test r.ap[2] == 0.0
@test r.r2 ≈ 1 atol = 1e-4
rn = NeuroAnalyzer._fooof_nan(fq, pap)
@test all(isnan, rn.ap) && isnan(rn.r2) && size(rn.pe) == (0, 3)
@test NeuroAnalyzer._fooof_check_f(fq) === nothing
@test_throws ArgumentError NeuroAnalyzer._fooof_check_f([1.0, 2.0, 4.0, 5.0])
@test_throws ArgumentError NeuroAnalyzer._fooof_check_f([1.0, 2.0])
@test NeuroAnalyzer._fooof_check_p(pap) === nothing
@test_throws ArgumentError NeuroAnalyzer._fooof_check_p(-pap)
@test_throws ArgumentError NeuroAnalyzer._fooof_check_p([NaN])
@test NeuroAnalyzer._fooof_check_settings((0.5, 12.0), 5, 2.0, :fixed) === nothing
@test_throws ArgumentError NeuroAnalyzer._fooof_check_settings((0.5, 12.0), 5, 2.0, :x)
# 0.5 Hz < 2 × resolution (1 Hz): a warning is printed
wout = mktemp() do path, io
    redirect_stdout(io) do
        NeuroAnalyzer._fooof_width_warn(fq, (0.5, 12.0))
    end
    flush(io)
    read(path, String)
end
@test occursin("overfitting", wout) || !NeuroAnalyzer.verbose

@info "Test: Granger internals"
# VAR(2): channel 2 is driven by channel 1
rng_ = NeuroAnalyzer.Random.MersenneTwister(4)
xv = zeros(2, 4000, 1)
for i = 3:4000
    xv[1, i, 1] = 0.5 * xv[1, i - 1, 1] - 0.3 * xv[1, i - 2, 1] + randn(rng_)
    xv[2, i, 1] = 0.4 * xv[2, i - 1, 1] + 0.6 * xv[1, i - 1, 1] + randn(rng_)
end
Y, Z = NeuroAnalyzer._gc_design(xv, 2)
@test size(Y) == (3998, 2) && size(Z) == (3998, 4)
@test Y[1, :] == xv[:, 3, 1]
@test NeuroAnalyzer._gc_rss(Y[:, 1], Z) ≈ sum(abs2, Y[:, 1] .- Z * (Z \ Y[:, 1]))
@test NeuroAnalyzer._gc_order(xv, 6) == 2
A, Σ = NeuroAnalyzer._gc_var_fit(xv, 2)
@test A[1, 1, 1] ≈ 0.5 atol = 0.05
@test A[1, 1, 2] ≈ -0.3 atol = 0.05
@test A[2, 1, 1] ≈ 0.6 atol = 0.05
@test A[1, 2, 1] ≈ 0.0 atol = 0.05
@test Σ ≈ [1.0 0.0; 0.0 1.0] atol = 0.1
@test_throws ArgumentError NeuroAnalyzer._gc_design(xv[:, 1:2, :], 2)
# Benjamini-Hochberg adjusted p-values: monotone in the raw p-values, ≥ raw, ≤ 1, diagonal untouched
pv = [1.0 0.01 0.04; 0.03 1.0 0.5; 0.02 0.9 1.0]
adj = NeuroAnalyzer._gc_fdr(pv)
off = [CartesianIndex(i, j) for i in 1:3 for j in 1:3 if i != j]
@test all(adj[k] >= pv[k] for k in off)
@test all(adj[k] <= 1 for k in off)
@test sortperm([adj[k] for k in off]) == sortperm([pv[k] for k in off]) || issorted([adj[k] for k in off][sortperm([pv[k] for k in off])])
@test adj[1, 2] ≈ 0.01 * 6 / 1
# Wilson factorization: S = H Σ Hᴴ at every frequency
N2 = 64
S = zeros(ComplexF64, 2, 2, N2)
for k in 1:N2
    ω = 2π * (k - 1) / N2
    Hk = inv([1.0 0.0; 0.0 1.0] .- [0.5 0.0; 0.6 0.4] .* exp(-im * ω))
    S[:, :, k] = Hk * Hk'
end
H, Σw = NeuroAnalyzer._gc_wilson(S)
@test all(H[:, :, k] * Σw * H[:, :, k]' ≈ S[:, :, k] for k in 1:N2)
@test Σw ≈ [1.0 0.0; 0.0 1.0] atol = 1e-6
gp = NeuroAnalyzer._gc_plus(S)
@test size(gp) == size(S)

true
