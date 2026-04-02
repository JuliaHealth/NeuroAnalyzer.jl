using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)
NeuroAnalyzer.filter!(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :lp,
    cutoff = 40,
    order = 91,
    bw = 0.5,
)
NeuroAnalyzer.filter!(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :hp,
    cutoff = 1,
    order = 91,
    bw = 0.5,
)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)
#=
@info "Test: acov()"
@test acov(v) == [-0.8, -0.8, -0.2, 0.8, 2.0, 0.8, -0.2, -0.8, -0.8]
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
@test size(ampdiff(a1)) == (2, 3, 2)
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
@test size(corm(a1)) == (2, 2, 2)
@test size(corm(e10, ch = "all")) == (24, 24, 10)

@info "Test: covm()"
@test covm(v) == [2.5;;]
@test size(covm(a1)) == (2, 2, 2)
@test size(covm(e10, ch = "all")) == (24, 24, 10)

@info "Test: cph()"
ph, f = cph(rand(10), rand(10); fs = 1)
@test length(ph) == 9
@test length(f) == 9
ph, f = cph(rand(10, 10, 2); fs = 1)
@test size(ph) == (10, 10, 9, 2)
@test length(f) == 9
ph, f = cph(e10; ch = "all")
@test size(ph) == (24, 24, 2049, 10)
@test length(f) == 2049
ph, f = cph(rand(10, 10, 2), rand(10, 10, 2); fs = 1)
@test size(ph) == (10, 9, 2)
@test length(f) == 9
ph, f = cph(e10, e10; ch1 = ["Fp1", "Fp2"], ch2 = ["Fp1", "Fp2"], ep1 = 1, ep2 = 1)
@test size(ph) == (2, 2049, 1)
@test length(f) == 2049

@info "Test: entropy()"
e, sh, l, s, ns = NeuroAnalyzer.entropy(rand(100))
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
n = NeuroAnalyzer.negentropy(rand(10000))
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
@test size(e) == (1, 15, 10)
@test length(t) == 15
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
@test size(e) == (1, 129, 10)
@test length(t) == 129
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
v = erp_auc(e; ch = "all")
@test length(v) == 19
=#
@info "Test: coherence()"
coh_data = NeuroAnalyzer.coherence(rand(100), rand(100); fs = 10, wstep = 5, method = :mt)
@test length(coh_data.coh) == 51
@test length(coh_data.imcoh) == 51
@test length(coh_data.msc) == 51
@test length(coh_data.f) == 51
coh_data = NeuroAnalyzer.coherence(rand(100), rand(100); fs = 10, wstep = 5, method = :fft)
@test length(coh_data.coh) == 51
@test length(coh_data.imcoh) == 51
@test length(coh_data.msc) == 51
@test length(coh_data.f) == 51
coh_data = NeuroAnalyzer.coherence(rand(100), rand(100); fs = 10, wstep = 5, method = :stft)
@test length(coh_data.coh) == 6
@test length(coh_data.imcoh) == 6
@test length(coh_data.msc) == 6
@test length(coh_data.f) == 6
coh_data = NeuroAnalyzer.coherence(
    e10,
    e10;
    ch1 = "Fp1",
    ch2 = "Fp2",
    ep1 = 1,
    ep2 = 1,
    method = :mt,
)
@test size(coh_data.coh) == (1, 1281, 1)
@test size(coh_data.imcoh) == (1, 1281, 1)
@test size(coh_data.msc) == (1, 1281, 1)
@test length(coh_data.f) == 1281
coh_data = NeuroAnalyzer.coherence(
    e10,
    e10;
    ch1 = "Fp1",
    ch2 = "Fp2",
    ep1 = 1,
    ep2 = 1,
    method = :fft,
)
@test size(coh_data.coh) == (1, 1281, 1)
@test size(coh_data.imcoh) == (1, 1281, 1)
@test size(coh_data.msc) == (1, 1281, 1)
@test length(coh_data.f) == 1281
coh_data = NeuroAnalyzer.coherence(
    e10,
    e10;
    ch1 = "Fp1",
    ch2 = "Fp2",
    ep1 = 1,
    ep2 = 1,
    method = :stft,
)
@test size(coh_data.coh) == (1, 129, 1)
@test size(coh_data.imcoh) == (1, 129, 1)
@test size(coh_data.msc) == (1, 129, 1)
@test length(coh_data.f) == 129

@info "Test: frqinst()"
f = NeuroAnalyzer.frqinst(rand(100))
@test length(f) == 100
f = NeuroAnalyzer.frqinst(rand(10, 100, 10))
@test size(f) == (10, 100, 10)
f = NeuroAnalyzer.frqinst(e10; ch = "all")
@test size(f) == (24, 2560, 10)

@info "Test: ged()"
s, r, rn = ged(rand(10, 10), rand(10, 10))
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
@test size(p) == (129, 1)
@test length(f) == 129
p, f = erop(e10; ch = "Fp1", method = :gh)
@test size(p) == (129, 1)
@test length(f) == 129

@info "Test: acor()"
@test acor(v) == [-0.32, -0.32, -0.08, 0.32, 0.8, 0.32, -0.08, -0.32, -0.32]
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
iv, ia, sd, pd, s1p, s2p = ispc(v1, v2)
@test iv == 0.6125992852305386
@test ia == 0.0017801930770334259
@test sd == [-5.0, -3.0, -1.0, 1.0, 3.0]
@test pd == [1.3157044982273685, -0.8713795327960081, -0.3743702916488455, -0.7615478999167378, 1.0328436072470222]
@test s1p ≈ [1.039406675134543, -0.6027563879589182, -0.21331750626000984, -0.33140501474755424, 0.3279718365439915]
@test s2p ≈ [-0.27629782309282525, 0.26862314483709, 0.16105278538883577, 0.4301428851691834, -0.7048717707030308]
iv, ia = ispc(e10; ch = "all")
@test size(iv) == (24, 24, 10)
@test size(ia) == (24, 24, 10)
iv, ia, sd, pd, s1p, s2p = ispc(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test iv ≈ [0.992049536181781;;]
@test ia ≈ [0.00013424580155568734;;]
@test size(sd) == (1, 2560, 1)
@test size(pd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: itpc()"
iv, izv, ia, ip = itpc(ones(1, 10, 10); t = 1)
@test iv == 1.0
@test izv == 10.0
@test ia == 0.0
@test ip == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
iv, izv, ia, ip = itpc(e10; ch = "Fp1", t = 1)
@test iv ≈ [0.9997997474350548]
@test izv ≈ [9.995995349711995]
@test ia ≈ [-0.0036001408903651292]
@test ip[1] ≈ 0.030644553057378048

@info "Test: itpc_spec()"
iv, izv, f = itpc_spec(e10; ch = "Fp1", flim = (0, 4), nfrq = 5)
@test size(iv) == (5, 2560)
@test size(izv) == (5, 2560)
@test f == [0.0, 1.0, 2.0, 3.0, 4.0]

@info "Test: mdiff()"
st, sts, p = mdiff(m1, m2; method = :absdiff)
@test length(st) == 6
@test sts == 3.0
@test p >= 0.0 && p <= 1.0
st, sts, p = mdiff(a1, a2; method = :absdiff)
@test size(st) == (2, 6)
@test sts == [1.0, 1.0]
@test p == [0.0, 0.0]
st, sts, p = mdiff(m1, m2; method = :diff2int)
@test length(st) == 6
@test sts == 4.666666666666666
@test p >= 0.0 && p <= 1.0
st, sts, p = mdiff(a1, a2; method = :diff2int)
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
@test NeuroAnalyzer.mutual_information(v1, v2) ≈ 0.4199730940219748
@test NeuroAnalyzer.mutual_information(a1) == [0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0]
@test NeuroAnalyzer.mutual_information(a1, a2) == [0.0 0.0; 0.0 0.0]
m = NeuroAnalyzer.mutual_information(e10; ch = "all")
@test size(m) == (24, 24, 10)
m = NeuroAnalyzer.mutual_information(e10, e10; ch1 = "Fp1", ch2 = "Fp2")
@test size(m) == (1, 10)

@info "Test: msci95()"
@test msci95(v1) == (sm = 3.0, se = 0.7071067811865476, ll = 1.6140707088743669, ul = 4.385929291125633)
@test msci95(v2) == (sm = 4.0, se = 0.7071067811865476, ll = 2.614070708874367, ul = 5.385929291125633)
@test msci95(m1) == (sm = [2.5, 3.5, 4.5], se = [1.4999999999999998, 1.4999999999999998, 1.4999999999999998], ll = [-0.4399999999999995, 0.5600000000000005, 1.5600000000000005], ul = [5.4399999999999995, 6.4399999999999995, 7.4399999999999995])
@test msci95(a1) == (sm = [1.0 1.0 1.0; 1.0 1.0 1.0], se = [0.0 0.0 0.0; 0.0 0.0 0.0], ll = [1.0 1.0 1.0; 1.0 1.0 1.0], ul = [1.0 1.0 1.0; 1.0 1.0 1.0])
@test msci95(v1, v2) == (sm = -1.0, se = 1.0, ll = -2.96, ul = 0.96)
@test msci95(m1, m2) == (sm = [-4.0; 2.0;;], se = [0.8164965809277261; 0.8164965809277261;;], ll = [-5.600333298618343; 0.39966670138165683;;], ul = [-2.3996667013816566; 3.6003332986183434;;])
@test msci95(a1, a2) == (sm = [1.0 1.0; 1.0 1.0], se = [0.0 0.0; 0.0 0.0], ll = [1.0 1.0; 1.0 1.0], ul = [1.0 1.0; 1.0 1.0])
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
@test size(s) == (257, 15, 1)
@test length(f) == 257
@test length(t) == 15
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
@test phdiff(a1, avg = :phase, h = true) == zeros(2, 3, 2)
p = phdiff(e10; ch = "all", avg = :phase)
@test size(p) == (24, 1281, 10)
p = phdiff(e10; ch = "all", avg = :phase)
@test size(p) == (24, 1281, 10)
p = phdiff(e10; ch = "all", avg = :phase, h = true)
@test size(p) == (24, 2560, 10)
p = phdiff(e10; ch = "all", avg = :phase, h = true)
@test size(p) == (24, 2560, 10)

@info "Test: pli()"
pv, phd, s1ph, s2ph = pli(v1, v2)
@test pv == 0.2
@test phd == [-5.0, -3.0, -1.0, 1.0, 3.0]
@test length(s1ph) == 5
@test length(s2ph) == 5
pv = pli(e10; ch = "all");
@test size(pv) == (24, 24, 10)
pv, sd, phd, s1p, s2p = pli(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test pv == [0.265625;;]
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: plv()"
pv, phd, s1ph, s2ph = plv(v1, v2)
@test pv == 0.6125992852305386
@test phd == [-5.0, -3.0, -1.0, 1.0, 3.0]
@test length(s1ph) == 5
@test length(s2ph) == 5
pv = plv(e10; ch = "all");
@test size(pv) == (24, 24, 10)
pv, sd, phd, s1p, s2p = plv(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
@test pv == [0.992049536181781;;]
@test size(sd) == (1, 2560, 1)
@test size(phd) == (1, 2560, 1)
@test size(s1p) == (1, 2560, 1)
@test size(s2p) == (1, 2560, 1)

@info "Test: psd()"
p, f = NeuroAnalyzer.psd(rand(100); fs = 10, wlen = 10, woverlap = 5)
@test length(p) == 6
@test f == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
p, f = NeuroAnalyzer.psd(rand(10, 100); fs = 10, wlen = 10, woverlap = 5)
@test size(p) == (10, 6)
p, f = NeuroAnalyzer.psd(rand(10, 100, 10); fs = 10, wlen = 10, woverlap = 5)
@test size(p) == (10, 6, 10)
p, f = NeuroAnalyzer.psd(rand(100); fs = 10, method = :mt, woverlap = 5)
@test length(p) == 51
@test round.(f, digits = 3) == 0.0:0.1:5.0
p, f = NeuroAnalyzer.psd(rand(10, 100); fs = 10, method = :mt, woverlap = 5)
@test size(p) == (10, 51)
p, f = NeuroAnalyzer.psd(rand(10, 100, 10); fs = 10, method = :mt, woverlap = 5)
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
@test size(p) == (1, 129, 10)
@test round.(f, digits = 3) == 0.0:1.0:128.0
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
p, f = psd_rel(rand(100); fs = 10, flim = (0, 1), wlen = 10, woverlap = 5)
@test length(p) == 6
@test f == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
p, f = psd_rel(rand(10, 100); fs = 10, flim = (0, 1), wlen = 10, woverlap = 5)
@test size(p) == (10, 6)
p, f = psd_rel(rand(10, 100, 10); fs = 10, flim = (0, 1), wlen = 10, woverlap = 5)
@test size(p) == (10, 6, 10)
p, f = psd_rel(rand(100); fs = 10, method = :mt, flim = (0, 1), wlen = 10, woverlap = 5)
@test length(p) == 51
@test round.(f, digits = 3) == 0.0:0.1:5.0
p, f = psd_rel(rand(10, 100); fs = 10, method = :mt, flim = (0, 1), wlen = 10, woverlap = 5)
@test size(p) == (10, 51)
p, f = psd_rel(
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
@test size(p) == (1, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0
p, f = psd_rel(e10; ch = "Fp1", method = :gh, flim = (0, 1))
@test size(p) == (1, 129, 10)
@test f[1] == 0.0
@test f[end] == 128.0

@info "Test: psd_slope()"
lf, ls, pf = psd_slope(rand(100); fs = 10, wlen = 10, woverlap = 5)
@test length(lf) == 6
@test pf == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
lf, ls, pf = psd_slope(rand(10, 100, 10); fs = 10, wlen = 10, woverlap = 5)
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
@test size(lf) == (1, 129, 10)
@test size(ls) == (1, 10)
@test pf[1] == 0.0
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
@test NeuroAnalyzer.rms(v1) == 3.3166247903554
@test NeuroAnalyzer.rms(a1) == [1.0 1.0; 1.0 1.0]
@test NeuroAnalyzer.rms(e10, ch = "Fp1") ≈ [363.81145362705735 349.8499943131811 336.52324682075516 329.7427258961967 317.29058937434763 311.9687522577657 285.7799209936517 267.59495573180266 259.6385867764095 243.7889621794278]

@info "Test: rmse()"
@test NeuroAnalyzer.rmse(v1, v2) == 3.0
@test NeuroAnalyzer.rmse(a1, a2) == [1.0 1.0; 1.0 1.0]
@test length(NeuroAnalyzer.rmse(e10, e10, ch1 = "Fp1", ch2 = "Fp2")) == 10

@info "Test: snr()"
@test NeuroAnalyzer.snr(v1, v2) == 3.0102999566398116
@test NeuroAnalyzer.snr(v1) == 1.8973665961010275
@test NeuroAnalyzer.snr2(v1) == 1.2060453783110545
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
@test size(sp) == (257, 15, 1, 10)
@test length(sf) == 257
@test length(st) == 15
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
sp, sst, t, f = spec_seg(sp, sf, st; ch = 1, t = (0, 1), f = (0, 10))
@test size(sp) == (11, 10, 10)
@test t == (1, 10)
@test f == (1, 11)

@info "Test: stationarity()"
s = NeuroAnalyzer.stationarity(e10; ch = "all", method = :adf)
@test size(s) == (24, 2, 10)
s = NeuroAnalyzer.stationarity(e10; ch = "all", method = :cov)
@test size(s) == (257, 10)
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
@test size(pxy) == (1, 2049, 1)
@test length(f) == 2049
pxy, f = cpsd(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1, method = :fft)
@test size(pxy) == (1, 2049, 1)
@test length(f) == 2049
pxy, f = cpsd(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1, method = :stft)
@test size(pxy) == (1, 257, 1)
@test length(f) == 257

@info "Test: tkeo()"
@test tkeo(v1) == [1.0, 1.0, 1.0, 1.0, 5.0]
@test tkeo(a1) == [1.0 0.0 1.0; 1.0 0.0 1.0;;; 1.0 0.0 1.0; 1.0 0.0 1.0]
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
@test size(pac) == (24, 5, 10)
@test length(l) == 5

@info "Test: xcov()"
xc, l = xcov(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = xcov(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, biased = false)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
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
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l =
    xcov(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, method = :cov, biased = true)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = xcov(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, method = :stat)
@test size(xc) == (1, 3, 1)
@test length(l) == 3

@info "Test: xcor()"
xc, l = xcor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = xcor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, biased = false)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = xcor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, method = :cor)
@test size(xc) == (1, 3, 1)
@test length(l) == 3
xc, l = xcor(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 2, method = :stat)
@test size(xc) == (1, 3, 1)
@test length(l) == 3

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
xc = xcor(x, y; l = 4, demean = false)
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

@info "Test: phsd()"
ph, f = phsd(e10; ch = "all")
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
@test sym_idx(v) == 19.0
@test sym_idx(e10, ch = "Fp1") ==
      [852.3333333333334 1279.0 852.3333333333334 852.3333333333334 1279.0 852.3333333333334 639.0 852.3333333333334 1279.0 1279.0]

@info "Test: lat_idx()"
@test lat_idx(e10) isa Float64
@test lat_idx(e10, frq = (1, 3.5)) isa Float64

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
@test sumsim(v1, v2, theta = 1) == 0.0012208548944264495
@test sumsim(a1, a2, theta = 1) ==
      [0.17692120631776423 0.17692120631776423; 0.17692120631776423 0.17692120631776423]
@test sumsim(e10, e10, ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1, theta = 0.0001) ≈
      [0.33211228;;]

@info "Test: hfd()"
@test hfd([1.0, 2.0, 3.5, 2.0, 5.0, 11.0, 2.0, 11.0]) == 0.8604486476012065
@test size(hfd(e10, ch = "all")) == (24, 10)

@info "Test: dirinrg()"
@test dirinrg([1.0, 2.0, 3.5, 2.0, 5.0, 11.0, 2.0, 11.0]) == 212.5
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
c, a, p, ph = NeuroAnalyzer.transform(rand(10, 100, 10); h = false, nf = false)
@test size(c) == (10, 51, 10)
@test size(a) == (10, 51, 10)
@test size(p) == (10, 51, 10)
@test size(ph) == (10, 51, 10)
c, a, p, ph = NeuroAnalyzer.transform(rand(10, 100, 10); h = false, nf = true)
@test size(c) == (10, 100, 10)
@test size(a) == (10, 100, 10)
@test size(p) == (10, 100, 10)
@test size(ph) == (10, 100, 10)
c, a, p, ph = NeuroAnalyzer.transform(rand(10, 100, 10); h = true)
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
p, ph, f, t = hhtspectrogram(imf; fs = sr(e10))
@test size(p) == (128, 2560)
@test size(ph)[2] == 2560
@test length(f) == 128
@test length(t) == 2560

@info "Test: hmspectrum()"
p, t = hmspectrum(e10; ch = "Fp1")
@test size(p) == (1, 2560, 10)
@test length(t) == 2560

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

true
