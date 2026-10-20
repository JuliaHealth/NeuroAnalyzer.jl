using NeuroAnalyzer
using Test

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
n = import_snirf(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "Test: apply()"
@test size(apply(e10; ch = "all", f = "mean(obj, dims=1)")) == (24, 1, 10)

@info "Test: l1()"
@test l1(a1, a2) == 12

@info "Test: l2()"
@test l2(a1, a2) == 3.4641016151377544

@info "Test: perm_cmp()"
z, b = perm_cmp(a1, a2)
@test round.(z) == -2.0 .* ones(2, 3)
@test size(b) == (2, 3)

@info "Test: vsearch()"
@test vsearch(2.1, [1, 2, 3, 4]) == 2
@test vsearch([2.1, 2.9], [1, 2, 3, 4]) == [2.0, 3.0]
@test vsearch("a", ["a", "b", "c"]) == 1
@test vsearch("d", ["a", "b", "c"]) === nothing

@info "Test: trim()"
@test vsplit(1:10, 2) == [1:2, 3:4, 5:6, 7:8, 9:10]

@info "Test: fft0()"
x = fft0(v1)
@test length(x) == 5
x = fft0(v1, 10)
@test length(x) == 15
x = rfft0(v1)
@test length(x) == 3
x = rfft0(v1, 10)
@test length(x) == div(length(v1) + 10, 2) + 1

@info "Test: ifft0()"
x = fft0(v1, 2)
@test round.(real.(ifft0(x, 2))) == [1.0, 2.0, 3.0, 4.0, 5.0]

@info "Test: fft2()"
x = fft2(v1)
@test length(x) == 8
x = rfft2(v1)
@test length(x) == 5

@info "Test: nextpow()"
@test nextpow2(5) == 8

@info "Test: gradient()"
g, g_l = NeuroAnalyzer.gradient(rand(10))
@test size(g) == (10,)
@test size(g_l) == (10,)
g, g_l = NeuroAnalyzer.gradient(rand(10, 10))
@test size(g) == (10, 10)
@test size(g_l) == (10, 10)
g, g_l = NeuroAnalyzer.gradient(rand(10, 10, 10))
@test size(g) == (10, 10, 10)
@test size(g_l) == (10, 10, 10)

@info "Test: findpeaks()"
x = rand(100)
x[10] *= 1000
@test 10 in findpeaks(x)

@info "Test: hz2rads()"
@test hz2rads(1) == 2pi

@info "Test: rads2hz()"
@test rads2hz(2pi) == 1

@info "Test: t2f()"
@test t2f(1000) == 1.0

@info "Test: f2t()"
@test f2t(1.0) == 1000.0

@info "Test: e2t()"
@test e2t(e10; ep = 1) == (0.0, 2559 / 256)
@test e2t(e10; ep = 1:10) == (0.0, 25599 / 256)

@info "Test: freqs()"
f, nf = NeuroAnalyzer._freqs(0:(1 / 10):10)
@test length(f) == 51
@test nf == 5
f, nf = NeuroAnalyzer._freqs(0:(1 / 10):10; nf = true)
@test length(f) == 101
@test nf == 5
f, nf = NeuroAnalyzer._freqs(rand(100), 10)
@test length(f) == 51
@test nf == 5
f, nf = NeuroAnalyzer._freqs(rand(100), 10; nf = true)
@test length(f) == 100
@test nf == 5
f, nf = NeuroAnalyzer.freqs(100, 10)
@test length(f) == 51
@test nf == 5
f, nf = NeuroAnalyzer.freqs(100, 10; nf = true)
@test length(f) == 100
@test nf == 5
f, nf = NeuroAnalyzer.freqs(e10)
@test length(f) == 1281
@test nf == 128
f, nf = NeuroAnalyzer.freqs(e10; nf = true)
@test length(f) == 2560
@test nf == 128

@info "Test: generate_window()"
s = generate_window(:hann, 100)
@test length(s) == 100
s = generate_window(:bh, 100)
@test length(s) == 100
s = generate_window(:bohman, 100)
@test length(s) == 100
s = generate_window(:flat, 100)
@test length(s) == 100
s = generate_window(:bn, 100)
@test length(s) == 100
s = generate_window(:nutall, 100)
@test length(s) == 100
s = generate_window(:triangle, 100)
@test length(s) == 100
s = generate_window(:exp, 100)
@test length(s) == 100

@info "Test: generate_sine()"
s = generate_sine(10, 1:100)
@test length(s) == 100

@info "Test: generate_cosine()"
s = generate_cosine(10, 1:100)
@test length(s) == 100

@info "Test: generate_csine()"
s = generate_csine(10, 1:100)
@test length(s) == 100

@info "Test: generate_sinc()"
s = generate_sinc()
@test length(s) == 401

@info "Test: generate_morlet()"
s = generate_morlet(100, 10)
@test length(s) == 201

@info "Test: generate_gaussian()"
s = generate_gaussian(100, 10)
@test length(s) == 201

@info "Test: generate_noise()"
s = generate_noise(100)
@test length(s) == 100

@info "Test: generate_morlet_fwhm()"
s = generate_morlet_fwhm(100, 10)
@test length(s) == 201

@info "Test: sr()"
@test sr(e10) == 256

@info "Test: nchannels()"
@test nchannels(e10) == 24

@info "Test: nepochs()"
@test nepochs(e10) == 10

@info "Test: signal_len()"
@test signal_len(e10) == 25600

@info "Test: epoch_len()"
@test epoch_len(e10) == 2560

@info "Test: get_channel()"
@test length(get_channel(e10; type = ["eeg", "eeg", "ecg", "mrk"])) == 20
@test length(get_channel(e10; ch = r"Fp.*")) == 2

@info "Test: cwtfrq()"
s = rand(100)
@test length(NeuroAnalyzer._cwtfrq(s, fs = 10)) == 12
@test length(cwtfrq(e10)) == 30

@info "Test: history()"
@test NeuroAnalyzer.history(e10) isa Vector{String}

@info "Test: labels()"
@test length(labels(e10)) == 24

@info "Test: channel_cluster()"
@test channel_cluster(e10; cluster = :f1) == ["Fp1", "F3", "F7"]

@info "Test: band_frq()"
@test band_frq(256; band = :alpha) == (8.0, 13.0)
@test band_frq(e10; band = :alpha) == (8.0, 13.0)

@info "Test: m_pad0()"
@test m_pad0(m1) == [1 2 3; 4 5 6; 0 0 0]
@test size(m_pad0(m1, 6, 5)) == (6, 5)

@info "Test: m_sortperm()"
@test m_sortperm(m1) == [1 1 1; 2 2 2]

@info "Test: m_sort()"
@test m_sort(m1, [2, 1]) == [4 5 6; 1 2 3]

@info "Test: m_norm()"
@test m_norm(m1) == [0.5 1.0 1.5; 2.0 2.5 3.0]

@info "Test: linspace()"
@test linspace(1, 10, 10) == 1:10

@info "Test: logspace()"
@test logspace(1, 10, 5) ==
      [1.0, 1.7782794100389228, 3.1622776601683795, 5.623413251903491, 10.0]

@info "Test: cmax()"
@test cmax([1 + 2im, 10 + 10im]) == 10 + 10im

@info "Test: cmin()"
@test cmin([1 + 2im, 10 + 10im]) == 1 + 2im

@info "Test: cums()"
@test cums(a1) == [1.0 2.0 3.0; 1.0 2.0 3.0;;; 1.0 2.0 3.0; 1.0 2.0 3.0]

@info "Test: f_nearest()"
m = [(1.0, 2.0) (3.0, 4.0); (5.0, 6.0) (7.0, 8.0)]
@test f_nearest(m, (3.2, 3.9)) == (1, 2)

@info "Test: view_note()"
@test view_note(e10) == ""

@info "Test: add_note()"
add_note!(e10; note = "test")
@test view_note(e10) == "test"

@info "Test: delete_note()"
delete_note!(e10)
@test view_note(e10) == ""

@info "Test: pad0()"
@test pad0(v1, 5) == [1, 2, 3, 4, 5, 0, 0, 0, 0, 0]
@test pad0(m1, 5) == [1 2 3 0 0 0 0 0; 4 5 6 0 0 0 0 0]

@info "Test: pad2()"
@test pad2(v1) == [1, 2, 3, 4, 5, 0, 0, 0]
@test pad2(m1) == [1 2 3 0; 4 5 6 0]

@info "Test: phases()"
@test length(phases(e10.data[1, :, 1])) == 2560

@info "Test: pick()"
@test channel_pick(e10; pick = [:l, :f]) == ["F3", "F7", "Fp1"]

@info "Test: t2s()"
@test t2s(0, 256) == 1
@test t2s(1.0, 256) == 257
@test t2s(e10; t = 1.0) == 257
@test t2s(0.3, 10) == 4
@test all(t2s(s2t(s, 256), 256) == s for s = 1:5000)

@info "Test: s2t()"
@test s2t(1, 256) == 0.0
@test s2t(2560, 256) == 2559 / 256
@test s2t(e10; s = 257) == 1.0
@test s2t(2, 20_000) != s2t(3, 20_000)

@info "Test: t2n()"
@test t2n(10, 256) == 2560
@test t2n(0, 256) == 0
@test t2n(e10; t = 1.0) == 256
@test length(t2s(0, 256):(t2s(0, 256) + t2n(10, 256) - 1)) == 2560

@info "Test: n2t()"
@test n2t(2560, 256) == 10.0
@test n2t(0, 256) == 0.0
@test n2t(e10; n = 256) == 1.0
@test all(t2n(n2t(n, 256), 256) == n for n = 0:5000)

@info "Test: get_channel(wl)"
@test length(get_channel(n; wl = 760)) == 36

@info "Test: size()"
@test size(e10) == (24, 2560, 10)

@info "Test: generate_signal()"
s = generate_signal(100)
@test length(s) == 100

@info "Test: chtypes()"
@test length(chtypes(e10)) == 24

@info "Test: optode_labels()"
n = import_nirs(joinpath(testfiles_path, "nirs", "fnirs-test-nirs.nirs"))
@test length(optode_labels(n)) == 3
@test length(source_labels(n)) == 1
@test length(detector_labels(n)) == 2

@info "Test: tavg()"
@test tavg(a1) == ones(2, 3, 1)

@info "Test: padm()"
@test padm(ones(2), 2) == ones(4)
@test padm(ones(2, 4), 2; mode = :all) == ones(2, 6)
@test padm(ones(2, 4), 2; mode = :row) == ones(2, 6)
@test padm(ones(2, 4, 3), 2; mode = :all) == ones(2, 6, 3)
@test padm(ones(2, 4, 3), 2; mode = :row) == ones(2, 6, 3)

@info "Test: vec2mat()"
@test size(vec2mat(ones(10); wlen = 2, woverlap = 1)) == (9, 2)

@info "Test: arr2mat()"
@test size(arr2mat(rand(1, 10, 10))) == (10, 10)

@info "Test: minat()"
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
@test minat(v1, v2) == (6, 1)

@info "Test: maxat()"
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
@test maxat(v1, v2) == (2, 5)

@info "Test: paired_labels()"
l = ["ch1", "ch2", "ch3"]
@test length(paired_labels(l; unq = true)) == 6
@test length(paired_labels(l; unq = false)) == 9
@test length(paired_labels(l, l)) == 3

@info "Test: vreduce()"
x = rand(100)
f = linspace(0, 10, 100)
x2, f2 = vreduce(x, f)
@test length(x2) == length(f2) == 21

@info "Test: areduce()"
x = rand(10, 100)
f = linspace(0, 10, 100)
x2, f2 = areduce(x, f)
@test size(x2, 2) == length(f2) == 21
x = rand(10, 100, 5)
f = linspace(0, 10, 100)
x2, f2 = areduce(x, f)
@test size(x2, 2) == length(f2) == 21

@info "Test: info()"
@test isa(info(e10, df = true), DataFrame)

@info "Test: describe()"
@test isa(NeuroAnalyzer.describe(e10; df = true), DataFrame)

@info "Test: ntapers()"
@test ntapers(e10, df = 1) == 9

@info "Test: trtm()"
@test size(trtm(e10, ch = "Fp1")) == (10, 2560)

@info "Test: meshgrid()"
@test length(meshgrid(collect(range(-1, 1, 100)), collect(range(-1, 1, 100)))) == 2

@info "Test: aff_mni2tal()"
@test aff_mni2tal([10, 12, 14]) == [8.0, 8.32, 12.48]

@info "Test: mni2tal()"
@test mni2tal([10, 12, 14]) == [9.9, 12.2696, 12.2826]

@info "Test: aff_tal2mni()"
@test aff_tal2mni([9.9, 12.2692, 12.2826]) ==
      [12.15909090909091, 16.071340206185567, 13.544355670103092]

@info "Test: tal2mni()"
@test tal2mni([9.9, 12.2692, 12.2821]) == [10.0, 11.999613921643126, 13.999435493742183]

@info "Test: fir_order_bw()"
@test fir_order_bw(eeg; bw = 0.2, a = 50) == 2909
@test fir_order_bw(bw = 0.2, a = 50, fs = 256) == 2909

@info "Test: fir_order_f()"
@test fir_order_f(eeg; f = 35) == (32, 40)
@test fir_order_f(f = 35, fs = 256) == (32, 40)

@info "Test: iir_order()"
@test iir_order(eeg; fprototype = :butterworth, ftype = :lp, cutoff = 12, bw = 0.2) == 199
@test iir_order(eeg; fprototype = :butterworth, ftype = :hp, cutoff = 12, bw = 0.2) == 199
@test iir_order(eeg; fprototype = :butterworth, ftype = :bp, cutoff = (12, 15), bw = 0.2) ==
      28
@test iir_order(eeg; fprototype = :butterworth, ftype = :bs, cutoff = (12, 15), bw = 0.2) ==
      28
@test iir_order(;
    fprototype = :butterworth,
    ftype = :lp,
    cutoff = 12,
    bw = 0.2,
    fs = 256,
) == 199
@test iir_order(;
    fprototype = :butterworth,
    ftype = :hp,
    cutoff = 12,
    bw = 0.2,
    fs = 256,
) == 199
@test iir_order(;
    fprototype = :butterworth,
    ftype = :bp,
    cutoff = (12, 15),
    bw = 0.2,
    fs = 256,
) == 28
@test iir_order(;
    fprototype = :butterworth,
    ftype = :bs,
    cutoff = (12, 15),
    bw = 0.2,
    fs = 256,
) == 28

@info "Test: anonymize()"
eeg_tmp = anonymize(eeg)
@test eeg_tmp.header.subject[:first_name] == ""
@test eeg_tmp.header.subject[:middle_name] == ""
@test eeg_tmp.header.subject[:last_name] == ""


##
## values (formulas, synthetic signals, SciPy references)
##

@info "Test: fft0() / rfft0() / ifft0() / fft2() values"
x = [1.0, -2.0, 3.0, 0.5, 4.0]
@test fft0(x) ≈ NeuroAnalyzer.FFTW.fft(x)
@test fft0(x, 3) ≈ NeuroAnalyzer.FFTW.fft(vcat(x, zeros(3)))
@test rfft0(x, 3) ≈ NeuroAnalyzer.FFTW.rfft(vcat(x, zeros(3)))
@test real.(ifft0(fft0(x, 3), 3)) ≈ x
@test fft2(x) ≈ NeuroAnalyzer.FFTW.fft(vcat(x, zeros(3)))
@test rfft2(x) ≈ NeuroAnalyzer.FFTW.rfft(vcat(x, zeros(3)))
# DC bin = sum, Parseval
@test real(fft0(x)[1]) ≈ sum(x)
@test sum(abs2, fft0(x)) / length(x) ≈ sum(abs2, x)
@test nextpow2(8) == 8
@test nextpow2(1) == 1
@test_throws ArgumentError fft0(x, -1)

@info "Test: frequency conversions"
@test rads2hz(hz2rads(13.5)) ≈ 13.5
@test f2t(t2f(250)) ≈ 250
@test t2f(100) == 10.0
# DFT bins: k × fs / n
f, nf = NeuroAnalyzer.freqs(8, 8)
@test f == [0.0, 1.0, 2.0, 3.0, 4.0]
@test nf == 4.0
f, nf = NeuroAnalyzer.freqs(8, 8; nf = true)
# full spectrum, centered (fftshift order)
@test f == [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]

@info "Test: generate_window() values"
# SciPy scipy.signal.windows (sym = True), n = 7
@test generate_window(:hann, 7) ≈ [0.0, 0.25, 0.75, 1.0, 0.75, 0.25, 0.0] atol = 1e-12
@test generate_window(:bh, 7) ≈ [6.0e-5, 0.055645, 0.520575, 1.0, 0.520575, 0.055645, 6.0e-5] atol = 1e-10
@test generate_window(:bohman, 7) ≈ [0.0, 0.10899778104422927, 0.6089977810442293, 1.0, 0.6089977810442293, 0.10899778104422927, 0.0] atol = 1e-12
# flat top: SciPy uses more digits of the coefficients
@test generate_window(:flat, 7) ≈ [-0.000421051, -0.051263156, 0.19821053, 1.0, 0.19821053, -0.051263156, -0.000421051] atol = 1e-4
# Blackman-Nuttall = scipy.signal.windows.nuttall
@test generate_window(:bn, 7) ≈ [0.0003628, 0.0613345, 0.5292298, 1.0, 0.5292298, 0.0613345, 0.0003628] atol = 1e-10
# Nuttall (continuous first derivative): w(0) = a0 - a1 + a2 - a3, w(centre) = a0 + a1 + a2 + a3
w = generate_window(:nutall, 7)
@test w[1] ≈ 0.355768 - 0.487396 + 0.144232 - 0.012604
@test w[4] ≈ 1.0
@test generate_window(:triangle, 5) ≈ [1 / 3, 2 / 3, 1.0, 2 / 3, 1 / 3]
@test generate_window(:exp, 5) ≈ [1 / 3, 1 / 2, 1.0, 1 / 2, 1 / 3]
@test generate_window(:exp, 4) ≈ [1 / 2, 1.0, 1.0, 1 / 2]
@test length(generate_window(:hann, 5; even = true)) == 6
@test generate_window(:hann, 1) == [1.0]
@test_throws ArgumentError generate_window(:x, 5)
@test_throws ArgumentError generate_window(:hann, 0)

@info "Test: generate_*() values"
t = 0:0.001:1
s = generate_sine(5, t, 2, 90)
@test s ≈ 2 .* cos.(2π * 5 .* t)
@test generate_cosine(5, t, 2, -90) ≈ 2 .* sin.(2π * 5 .* t)
cs = generate_csine(5, t, 3)
@test real.(cs) ≈ generate_cosine(5, t, 3)
@test imag.(cs) ≈ generate_sine(5, t, 3)
@test all(abs.(cs) .≈ 3)
# sinc: value 2f at the peak, zeros at multiples of 1 / (2f)
s = generate_sinc(-2:0.25:2; f = 1.0)
@test s[9] == 2.0
@test abs(s[11]) < 1e-12
s = generate_sinc(-2:0.25:2; f = 1.0, peak = 0.5)
@test s[11] == 2.0
@test generate_sinc(-2:0.25:2; f = 1.0, norm = false)[9] == 2.0
# Gaussian: σ = ncyc / (2πf); value at ±σ = exp(-1/2)
g = generate_gaussian(1000, 10, 0.5; ncyc = 5)
σ = 5 / (2π * 10)
@test length(g) == 1001
@test g[501] == 1.0
@test g[501 + round(Int64, σ * 1000)] ≈ exp(-0.5) rtol = 1e-2
@test generate_gaussian(100, 10, 1; a = 3)[101] == 3.0
@test_throws ArgumentError generate_gaussian(0, 10)
@test_throws ArgumentError generate_gaussian(100, -1)
# Morlet: complex sine × Gaussian; the real part = cos × Gaussian
m = generate_morlet(1000, 10, 0.5; complex = true)
@test abs.(m) ≈ g
@test generate_morlet(1000, 10, 0.5) ≈ real.(m)
@test_throws ArgumentError generate_morlet(1000, 10, 0)
# FWHM Morlet: the envelope is 0.5 at ±h/2
mf = generate_morlet_fwhm(1000, 10, 1; h = 0.2)
@test abs(mf[1001]) == 1.0
@test abs(mf[1001 + 100]) ≈ 0.5
@test_throws ArgumentError generate_morlet_fwhm(1000, 10, 1; h = 0.0)
# noise: normalized to [-a, a]
for type in (:whiten, :whiteu, :pink)
    local s = generate_noise(1000, 2; type = type)
    @test extrema(s) == (-2.0, 2.0)
end
@test_throws ArgumentError generate_noise(10; type = :x)
# pink noise: power ∝ 1/f → log-log PSD slope ≈ -1
s = generate_noise(2^16; type = :pink)
p = abs2.(NeuroAnalyzer.FFTW.rfft(s))[2:end]
fr = 1:length(p)
@test -1.3 < ([ones(length(p)) log.(fr)] \ log.(p))[2] < -0.7
s = generate_signal(500, 3)
@test extrema(s) == (-3.0, 3.0)
@test_throws ArgumentError generate_signal(0)
# square wave: period 2, duty cycle (2 - w) / 2, amplitude a, offset
s = generate_square(0:0.5:3.5, 2, 0)
@test s == [0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0]
@test generate_square(0:0.5:3.5, 1, 0; w = 0.5, offset = -0.5) == [-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5]
# triangle: a × |mod(t, 2) - 1|
@test generate_triangle(0:0.5:2, 4) == [4.0, 2.0, 0.0, 2.0, 4.0]

@info "Test: fwhm()"
x = exp.(-0.5 .* ((-50:50) ./ 10) .^ 2)
p1, pk, p2 = fwhm(x)
@test pk == 51
# FWHM of a Gaussian = 2 √(2 ln 2) σ
@test (p2 - p1) ≈ 2 * sqrt(2 * log(2)) * 10 atol = 2
@test_throws ArgumentError fwhm([1.0])

@info "Test: findpeaks() values"
x = sin.(2π .* (0:999) ./ 100)
@test findpeaks(x; d = 10) == 26:100:926
# minimum distance: only one of two close peaks is kept
x = zeros(100)
x[[30, 35, 70]] .= [1.0, 2.0, 1.5]
@test findpeaks(x; d = 10) == [35, 70]
@test_throws ArgumentError findpeaks(x; d = 0)

@info "Test: gradient() values"
g, gl = NeuroAnalyzer.gradient(3.0 .* (1:10))
@test all(x -> x[1] ≈ 3, g)
g, gl = NeuroAnalyzer.gradient(3.0 .* (1:10); rev = true)
@test all(x -> x[1] ≈ -3, g)

@info "Test: l1() / l2() / tavg() / cums() values"
a = reshape(collect(1.0:12), 2, 3, 2)
b = zeros(2, 3, 2)
@test l1(a, b) == sum(1:12)
@test l2(a, b) ≈ sqrt(sum(abs2, 1:12))
@test l1(a, a) == 0.0
@test_throws ArgumentError l1(a, zeros(2, 3, 3))
@test tavg(a)[:, :, 1] == (a[:, :, 1] .+ a[:, :, 2]) ./ 2
@test cums(a)[:, :, 1] == cumsum(a[:, :, 1]; dims = 2)

@info "Test: perm_cmp() values"
# no difference → no significant cells; a large shift → all cells significant
a1 = randn(NeuroAnalyzer.Random.MersenneTwister(1), 2, 3, 20)
a2 = randn(NeuroAnalyzer.Random.MersenneTwister(2), 2, 3, 20) .+ 5
zm, bm = perm_cmp(a1, a2; perm_n = 500)
@test all(bm)
@test all(zm .> 0)
@test_throws ArgumentError perm_cmp(a1, zeros(2, 3, 2))
@test_throws ArgumentError perm_cmp(a1, a1; p = 1.5)

@info "Test: vsearch() / vsplit() / minat() / maxat() / vreduce() / areduce() values"
@test vsearch(3.6, [1.0, 2.0, 3.0, 4.0]) == 4
@test vsearch([0.0, 10.0], [1.0, 2.0, 3.0]) == [1, 3]
@test_throws ArgumentError vsplit(1:10, 3)
@test vsplit(1:6, 3) == [1:3, 4:6]
x = [5.0, 1.0, 9.0]
y = [10.0, 20.0, 30.0]
# (value of y at the extremum of x, index)
@test minat(x, y) == (20.0, 2)
@test maxat(x, y) == (30.0, 3)
f = collect(0.0:0.1:3.0)
x = collect(1.0:31.0)
x2, f2 = vreduce(x, f)
@test f2 == [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
@test x2 == x[[1, 6, 11, 16, 21, 26, 31]]
a = repeat(x', 2)
a2, f2 = areduce(a, f; n = 1.0)
@test f2 == [0.0, 1.0, 2.0, 3.0]
@test a2[1, :] == x[[1, 11, 21, 31]]
@test_throws ArgumentError areduce(rand(2), f)
@test_throws ArgumentError vreduce(x, f[1:3])

@info "Test: cmax() / cmin() / cextrema()"
x = [1 + 1im, -3 + 0im, 0.5im]
@test cextrema(x) == (-3 + 0im, 0 + 0.5im)
@test cextrema(x) == (cmax(x), cmin(x))

@info "Test: linspace() / logspace() values"
@test linspace(0, 1, 5) == [0.0, 0.25, 0.5, 0.75, 1.0]
# logspace: n geometrically spaced points from start to stop
@test logspace(1, 100, 3) ≈ [1.0, 10.0, 100.0]
@test_throws ArgumentError logspace(0, 2, 3)

@info "Test: m_*() / pad*() / vec2mat() / arr2mat() / meshgrid() values"
@test m_sortperm([3 1; 1 2]) == [2 1; 1 2]
@test m_sortperm([3 1; 1 2]; dims = 2) == [2 1; 1 2]
@test m_sort([3 1; 1 2], [2, 1]) == [1 2; 3 1]
@test m_pad0([1 2 3]) == [1 2 3; 0 0 0; 0 0 0]
@test pad0(ones(1, 2, 2), 1) == cat([1.0 1 0], [1.0 1 0]; dims = 3)
@test pad2(ones(1, 3, 1)) == reshape([1.0, 1, 1, 0], 1, 4, 1)
@test padm([1.0, 3.0], 2) == [1.0, 3.0, 2.0, 2.0]
@test vec2mat(collect(1:6); wlen = 3, woverlap = 1) == [1 2 3; 3 4 5]
@test arr2mat(reshape(collect(1.0:6), 1, 3, 2)) == [1.0 2 3; 4 5 6]
mx, my = meshgrid([1.0, 2.0, 3.0], [10.0, 20.0])
@test mx == [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]
@test my == [[10.0, 10.0, 10.0], [20.0, 20.0, 20.0]]
@test f_nearest([(1.0, 2.0) (3.0, 4.0); (5.0, 6.0) (7.0, 8.0)], (7.2, 7.9)) == (2, 2)

@info "Test: phases() values"
# Hilbert phase of a cosine with an integer number of cycles = 2π f t (wrapped)
t = (0:255) ./ 256
p = phases(cos.(2π * 8 .* t))
@test all(abs.(angle.(exp.(im .* (p .- 2π * 8 .* t)))) .< 1e-8)

@info "Test: sr() / signal_duration() / epoch_duration() / datatype() / channel_order()"
@test signal_duration(e10) == 100.0
@test epoch_duration(e10) == 10.0
@test signal_duration(eeg) == signal_len(eeg) / 256
@test datatype(e10) == "eeg"
@test channel_order(e10) == e10.header.recording[:channel_order]
@test size(e10, 3) == 10
@test nchannels(e10; type = "eeg") == length(get_channel(e10; type = "eeg"))

@info "Test: band_frq() values"
@test band_frq(256; band = :total) == (0.1, 128.0)
@test band_frq(256; band = :delta) == (0.5, 4.0)
@test band_frq(256; band = :gamma_higher) == (80.0, 127.9)
# band above the Nyquist frequency: truncated just below it
@test band_frq(100; band = :gamma_higher) == (49.8, 49.9)
@test band_frq(256; band = :list) == (0.0, 128.0)
@test_throws ArgumentError band_frq(256; band = :x)
@test NeuroAnalyzer._band_table(:theta, 128.0) == (4.0, 8.0)
@test isnothing(NeuroAnalyzer._band_table(:list, 128.0))
@test_throws ArgumentError NeuroAnalyzer._band_table(:x, 128.0)
@test NeuroAnalyzer._clamp_band((8.0, 13.0), 128.0, :alpha, "") == (8.0, 13.0)
@test NeuroAnalyzer._clamp_band((8.0, 130.0), 128.0, :alpha, "") == (8.0, 127.9)
@test NeuroAnalyzer._clamp_band((130.0, 150.0), 128.0, :alpha, "") == (127.8, 127.9)

@info "Test: printing (info, header, channel_info, view_header, _band_list)"
function capture(f)
    mktemp() do path, io
        redirect_stdout(io) do
            f()
        end
        flush(io)
        read(path, String)
    end
end
out = capture(() -> NeuroAnalyzer._print_header(e10))
@test occursin("Sampling rate (Hz): 256", out)
@test occursin("Number of epochs: 10", out)
out = capture(() -> NeuroAnalyzer.header(e10))
@test occursin("Data type: EEG", out)
out = capture(() -> info(e10))
@test occursin("Fp1", out)
@test names(info(e10; df = true))[1:2] == ["time", "Fp1"]
@test channel_info(e10; ch = "Fp1", pr = false) == " ch: 1    label: Fp1      type: EEG      unit: μV       bad: false"
out = capture(() -> channel_info(e10; ch = "Fp1"))
@test occursin("label: Fp1", out)
@test_throws ArgumentError channel_info(e10; ch = r"Fp" |> string)
nn = import_snirf(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf"))
@test occursin("wavelength", channel_info(nn; ch = labels(nn)[1], pr = false)) ||
      occursin("type:", channel_info(nn; ch = labels(nn)[1], pr = false))
out = capture(() -> view_header(e10))
@test occursin("Header fields: subject, recording, experiment", out)
@test occursin("header.recording[:sampling_rate]: 256", out)
out = capture(() -> NeuroAnalyzer._band_list([:list, :a, :b]))
@test out == "Available band names: :a, :b\n"
d = NeuroAnalyzer.describe(e10; df = true)
@test d[1, 1] == "Fp1" || d[1, :label] == "Fp1"

@info "Test: to_df()"
df = to_df(e10)
@test size(df) == (25600, 25)
@test names(df)[1:2] == ["time", "Fp1"]
@test df[!, :time] == e10.time_pts
@test df[2561, :Fp1] == e10.data[1, 1, 2]

@info "Test: anonymize!()"
eeg_tmp = deepcopy(eeg)
eeg_tmp.header.subject[:first_name] = "A"
anonymize!(eeg_tmp)
@test eeg_tmp.header.subject[:first_name] == ""
@test eeg_tmp.header.subject[:last_name] == ""
eeg_tmp.header.subject[:last_name] = "B"
e_tmp = anonymize(eeg_tmp)
@test e_tmp.header.subject[:last_name] == ""
@test eeg_tmp.header.subject[:last_name] == "B"

@info "Test: add_note() / delete_note()"
e_tmp = add_note(e10; note = "x")
@test view_note(e_tmp) == "x"
@test view_note(e10) == ""
@test view_note(delete_note(e_tmp)) == ""
@test view_note(e_tmp) == "x"

@info "Test: time conversion values"
@test e2t(e10; ep = 2:3) == (10.0, 30 - 1 / 256)
@test_throws ArgumentError e2t(e10; ep = 11)
@test_throws ArgumentError t2s(-1, 256)
@test_throws ArgumentError s2t(1, 0)
m = DataFrame(:id => ["a"], :start => [257.0], :length => [256.0], :value => ["x"], :channel => [0])
@test markers_s2t(m; fs = 256)[1, :start] == 1.0
@test markers_s2t(m; fs = 256)[1, :length] == 1.0
@test m[1, :start] == 257.0

@info "Test: trtm() / apply() values"
@test trtm(e10; ch = "Fp1") == e10.data[1, :, :]'
r = apply(e10; ch = "Fp1", f = "mean(obj, dims=1)")
@test r[1, 1, 3] ≈ sum(e10.data[1, :, 3]) / 2560

@info "Test: ntapers() values"
# K = T × df - 1 (T: epoch duration)
@test ntapers(e10; df = 2) == 19
@test ntapers(e10; df = 0.1) == 1
@test_throws ArgumentError ntapers(e10; df = 200)

@info "Test: get_channel() values"
@test get_channel(e10; ch = "Fp2") == [2]
# indices in channel order
@test get_channel(e10; ch = ["Fp2", "Fp1"]) == [1, 2]
@test get_channel(e10; type = "ecg") == labels(e10)[e10.header.recording[:channel_type] .== "ecg"]
@test_throws ArgumentError get_channel(e10; ch = "no_such_channel")

@info "Test: cwtfrq() values"
# frequencies within [0, fs / 2] (0: the averaging wavelet), ascending
f = cwtfrq(e10)
@test f[1] == 0
@test all(0 .< f[2:end] .<= 128)
@test issorted(f)

@info "Test: paired_labels() values"
@test paired_labels(["a", "b"]; unq = false) == ["a-a", "a-b", "b-a", "b-b"]
@test paired_labels(["a", "b"], ["c", "d"]) == ["a-c", "b-d"]

@info "Test: play()"
# argument validation only (no audio output in tests)
@test_throws ArgumentError play(e10; ch = "Fp1", seg = (0, 1), ep = 11)
@test_throws ArgumentError play(e10; ch = "Fp1", seg = (1, 0), ep = 1)

@info "Test: mni2tal() / tal2mni() round trips"
@test aff_tal2mni(aff_mni2tal([10.0, -20.0, 30.0])) ≈ [10.0, -20.0, 30.0]
@test tal2mni(mni2tal([10.0, -20.0, 30.0])) ≈ [10.0, -20.0, 30.0] atol = 1e-6
@test tal2mni(mni2tal([10.0, 20.0, -30.0])) ≈ [10.0, 20.0, -30.0] atol = 1e-6

true
