using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
nirs = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
meg = import_fiff(joinpath(testfiles_path, "meg-test-fiff.fif"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)
epoch!(meg; ep_len = 2)
keep_epoch!(meg; ep = 1:10)
load_locs!(e10; file_name = joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "Test: add_signal()"
@test NeuroAnalyzer._add_signal(v1, v2) == v1 + v2
x = rand(epoch_len(e10))
e10_tmp = add_signal(e10; ch = "Fp1", s = x)
e10_tmp.data[1, :, 1] == e10.data[1, :, 1] + x

@info "Test: cbp()"
@test length(NeuroAnalyzer._cbp(rand(100), fs = 10, frq = 4)) == 100
@test length(NeuroAnalyzer._cbp(rand(100), fs = 10, frq = 4, pad = 20)) == 100
# zero-phase, unit gain at frq (away from the edges)
x = cos.(2π * 4 .* (0:999) ./ 100)
@test NeuroAnalyzer._cbp(x, fs = 100, frq = 4)[300:700] ≈ x[300:700] atol = 0.05
e10_tmp = cbp(e10; ch = "all", frq = 4)
@test size(e10_tmp) == size(e10)

@info "Test: cwd()"
s = rand(100)
NeuroAnalyzer._log_off()
ct = NeuroAnalyzer._cwd(s; wt = wavelet(Morlet(π); β = 32, Q = 128))
NeuroAnalyzer._log_on()
@test size(ct) == (130, 100)
ct = NeuroAnalyzer.cwd(e10; ch = "Fp1", wt = wavelet(Morlet(π); β = 32, Q = 128));
@test size(ct) == (1, 131, 2560, 10)

@info "Test: icwd()"
# accurate reconstruction from complex coefficients
x = sin.(2π * 10 .* (0:1023) ./ 256) .+ 0.5 .* sin.(2π * 40 .* (0:1023) ./ 256)
r = icwd(NeuroAnalyzer._cwd(x; complex = true); type = :pd)
@test NeuroAnalyzer.norm(r .- x) / NeuroAnalyzer.norm(x) < 0.05
# artrem_cwd(): removing a band outside the signal leaves it (nearly) unchanged
r = NeuroAnalyzer._artrem_cwd(x, collect(0:1023) ./ 256; fs = 256, tseg = (1, 2), fseg = (70, 80))
@test NeuroAnalyzer.norm(r .- x) / NeuroAnalyzer.norm(x) < 0.05
NeuroAnalyzer._log_off()
ct = NeuroAnalyzer._cwd(s; wt = wavelet(Morlet(π); β = 32, Q = 128))
s_new = icwd(ct; wt = wavelet(Morlet(π); β = 32, Q = 128), type = :nd)
@test length(s_new) == 100
s_new = icwd(ct; wt = wavelet(Morlet(π); β = 32, Q = 128), type = :pd)
@test length(s_new) == 100
s_new = icwd(ct; wt = wavelet(Morlet(π); β = 32, Q = 128), type = :df)
@test length(s_new) == 100
NeuroAnalyzer._log_on()

@info "Test: denoise_fft()"
s = rand(100)
s2, f = NeuroAnalyzer._denoise_fft(s)
@test length(s2) == 100
@test length(f) == 100
# the sign of the signal is kept (no bin exceeds the threshold → the signal is unchanged)
x = sin.(2π .* (1:100) ./ 10)
@test NeuroAnalyzer._denoise_fft(x; t = 1e10).s ≈ x
e10_tmp = denoise_fft(e10; ch = "all")
@test size(e10_tmp) == size(e10)

@info "Test: denoise_dwd()"
s = NeuroAnalyzer._denoise_dwd(rand(100); wt = wavelet(WT.haar))
@test length(s) == 100
e10_tmp = denoise_dwd(e10; ch = "all", wt = wavelet(WT.haar))
@test size(e10_tmp) == size(e10)

@info "Test: denoise_wien()"
s = NeuroAnalyzer._denoise_wien(a1)
@test size(s) == (2, 3, 2)
e10_tmp = denoise_wien(e10; ch = "all")
@test size(e10_tmp) == size(e10)

@info "Test: derivative()"
@test NeuroAnalyzer._derivative(v1) == [0.5, 1.0, 1.0, 1.0, 0.5]
@test NeuroAnalyzer._derivative(a1) == [0.0 0.0 0.0; 0.0 0.0 0.0;;; 0.0 0.0 0.0; 0.0 0.0 0.0]
e10_tmp = NeuroAnalyzer.derivative(e10; ch = "all")
@test size(e10_tmp) == size(e10)

@info "Test: detrend()"
@test round.(NeuroAnalyzer._detrend(v1)) == [0.0, 0.0, 0.0, 0.0, 0.0]
e10_tmp = detrend(e10; ch = "all", type = :ls)
@test size(e10_tmp) == size(e10)
e10_tmp = detrend(e10; ch = "all", type = :linear)
@test size(e10_tmp) == size(e10)
e10_tmp = detrend(e10; ch = "all", type = :constant)
@test size(e10_tmp) == size(e10)
e10_tmp = detrend(e10; ch = "all", type = :poly)
@test size(e10_tmp) == size(e10)
e10_tmp = detrend(e10; ch = "all", type = :loess)
@test size(e10_tmp) == size(e10)

@info "Test: dwd()"
NeuroAnalyzer._log_off()
s = rand(100)
dt = NeuroAnalyzer._dwd(s; type = :sdwt)
@test size(dt) == (7, 100)
dt = NeuroAnalyzer._dwd(s; type = :acdwt)
@test size(dt) == (7, 100)
NeuroAnalyzer._log_on()
dt = dwd(e10; ch = "all", type = :sdwt, l = 2)
@test size(dt) == (24, 7, 2560, 10)
dt = dwd(e10; ch = "all", type = :acdwt, l = 2)
@test size(dt) == (24, 7, 2560, 10)

@info "Test: idwd()"
NeuroAnalyzer._log_off()
s = rand(256)
dt = NeuroAnalyzer._dwd(s; type = :sdwt, l = 2)
s_new = idwd(dt; type = :sdwt)
@test length(s_new) == 256
s = rand(256)
dt = NeuroAnalyzer._dwd(s; type = :acdwt, l = 2)
s_new = idwd(dt; type = :acdwt)
@test length(s_new) == 256
# exact reconstruction from the full packet decomposition
for ty in (:sdwt, :acdwt)
    @test idwd(NeuroAnalyzer._dwd(s; type = ty, l = 3); type = ty) ≈ s
end
NeuroAnalyzer._log_on()

@info "Test: average_epochs()"
e_erp = average_epochs(e10)
@test size(e_erp.data) == (19, 2560, 11)
m_erf = average_epochs(meg)
@test size(m_erf.data) == (306, 2000, 11)

@info "Test: bpsplit()"
s, bn, bf = bpsplit(e10; ch = "all")
@test length(bn) == 13
@test length(bf) == 13
@test size(s) == (13, 24, 2560, 10)

@info "Test: fconv()"
@test NeuroAnalyzer._fconv(v1, kernel = v2) == ComplexF64[
    1.5999999999999996 + 4.13548630602436e-17im,
    2.4999999999999996 + 0.0im,
    3.499999999999999 + 7.236218295908405e-17im,
    2.6999999999999997 - 2.632430811870764e-17im,
    1.8999999999999997 + 0.0im,
]
@test NeuroAnalyzer._fconv(a1, kernel = [0.5, 1.0, 0.5]) == ComplexF64[
    0.75 + 0.0im 1.0 + 0.0im 0.75 + 0.0im; 0.75 + 0.0im 1.0 + 0.0im 0.75 + 0.0im;;;
    0.75 + 0.0im 1.0 + 0.0im 0.75 + 0.0im; 0.75 + 0.0im 1.0 + 0.0im 0.75 + 0.0im
]
s_conv = fconv(e10; ch = "all", kernel = [0.0, 0.5, 1.0, 0.5, 0.0])
@test size(s_conv) == size(e10)

@info "Test: filter_mavg()"
# centered windows, truncated at the edges
@test NeuroAnalyzer._filter_mavg(vcat(v1, v1), k = 2) == [2.0, 2.5, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.5, 4.0]
# weighted mean: a constant signal is unchanged
@test NeuroAnalyzer._filter_mavg(fill(2.0, 20), k = 3, ww = NeuroAnalyzer.DSP.hanning(7)) ≈ fill(2.0, 20)
e10_tmp = filter_mavg(e10; ch = "all", k = 2)
@test size(e10_tmp) == size(e10)

@info "Test: filter_mmed()"
@test NeuroAnalyzer._filter_mmed(vcat(v1, v1), k = 2) == [2.0, 2.5, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.5, 4.0]
# weighted median: values are not rescaled by the weights
@test NeuroAnalyzer._filter_mmed(fill(2.0, 20), k = 3, ww = NeuroAnalyzer.DSP.hanning(7)) ≈ fill(2.0, 20)
e10_tmp = filter_mmed(e10; ch = "all", k = 2)
@test size(e10_tmp) == size(e10)

@info "Test: filter_poly()"
s = NeuroAnalyzer._filter_poly(rand(20))
@test length(s) == 20
e10_tmp = filter_poly(e10; ch = "all")
@test size(e10_tmp) == size(e10)

@info "Test: filter_sg()"
s = NeuroAnalyzer._filter_sg(rand(20))
@test length(s) == 20
e10_tmp = filter_poly(e10; ch = "all")
@test size(e10_tmp) == size(e10)

@info "Test: filter()"
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :butterworth,
    ftype = :lp,
    cutoff = 40,
    order = 8,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :butterworth,
    ftype = :hp,
    cutoff = 1,
    order = 12,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :butterworth,
    ftype = :bs,
    cutoff = (49, 51),
    order = 2,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :butterworth,
    ftype = :bp,
    cutoff = (49, 51),
    order = 4,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :chebyshev1,
    ftype = :lp,
    cutoff = 40,
    order = 8,
    rs = 35,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :chebyshev1,
    ftype = :hp,
    cutoff = 1,
    order = 12,
    rs = 2,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :chebyshev1,
    ftype = :bs,
    cutoff = (49, 51),
    order = 2,
    rs = 50,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :chebyshev1,
    ftype = :bp,
    cutoff = (49, 51),
    order = 4,
    rs = 50,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :chebyshev2,
    ftype = :lp,
    cutoff = 40,
    order = 8,
    rs = 35,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :chebyshev2,
    ftype = :hp,
    cutoff = 1,
    order = 12,
    rs = 2,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :chebyshev2,
    ftype = :bs,
    cutoff = (49, 51),
    order = 2,
    rs = 50,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :chebyshev2,
    ftype = :bp,
    cutoff = (49, 51),
    order = 4,
    rs = 50,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :elliptic,
    ftype = :lp,
    cutoff = 40,
    order = 8,
    rs = 35,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :elliptic,
    ftype = :hp,
    cutoff = 1,
    order = 12,
    rs = 2,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :elliptic,
    ftype = :bs,
    cutoff = (49, 51),
    order = 2,
    rs = 50,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :elliptic,
    ftype = :bp,
    cutoff = (49, 51),
    order = 4,
    rs = 50,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :lp,
    cutoff = 40,
    order = 91,
    bw = 0.5,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :hp,
    cutoff = 1,
    order = 91,
    bw = 0.5,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :bs,
    cutoff = (49, 51),
    order = 91,
    bw = 0.5,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :bp,
    cutoff = (49, 51),
    order = 91,
    bw = 0.5,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :firls,
    ftype = :lp,
    cutoff = 40,
    order = 32,
    bw = 2,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :firls,
    ftype = :hp,
    cutoff = 10,
    order = 31,
    bw = 2,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :firls,
    ftype = :bs,
    cutoff = (49, 51),
    bw = 0.5,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :firls,
    ftype = :bp,
    cutoff = (49, 51),
    bw = 0.5,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :remez,
    ftype = :lp,
    cutoff = 40,
    order = 32,
    bw = 10,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :remez,
    ftype = :hp,
    cutoff = 10,
    order = 31,
    bw = 0.5,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :remez,
    ftype = :bs,
    cutoff = (49, 51),
    order = 31,
    bw = 0.1,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :remez,
    ftype = :bp,
    cutoff = (49, 51),
    order = 31,
    bw = 0.5,
)
@test size(eeg_tmp) == size(e10)

@info "Test: filter_g()"
s = NeuroAnalyzer._filter_g(rand(10); fs = 10, f = 4)
@test length(s) == 10
# envelope of a sine at f equals its amplitude
@test NeuroAnalyzer._filter_g(2 .* cos.(2π * 10 .* (0:999) ./ 100); fs = 100, f = 10)[200:800] ≈ fill(2.0, 601) atol = 1e-3
e10_tmp = filter_g(e10; ch = "all", f = 20)
@test size(e10_tmp) == size(e10)

@info "Test: lrinterpolate_channel()"
e10_tmp = deepcopy(e10)
e10_tmp.data[1, :, 1] = zeros(epoch_len(e10))
e10_int = lrinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1)
@test e10_int.data[1, :, 1] != e10_tmp.data[1, :, 1]

@info "Test: normalize()"
@test NeuroAnalyzer._normalize(v1, method = :none) == v1
@test NeuroAnalyzer._normalize(m1, method = :none) == m1
@test NeuroAnalyzer._normalize(a1, method = :none) == a1
@test normalize_zscore(v1) == [
    -1.2649110640673518,
    -0.6324555320336759,
    0.0,
    0.6324555320336759,
    1.2649110640673518,
]
@test normalize_zscore(m1; bych = true) == [-1 0 1; -1 0 1]
@test normalize_zscore(m1; bych = false) == [
    -1.3363062095621219 -0.8017837257372732 -0.2672612419124244;
    0.2672612419124244 0.8017837257372732 1.3363062095621219
]
@test normalize_minmax(v1) == [-1.0, -0.5, 0.0, 0.5, 1.0]
@test normalize_minmax(m1; bych = true) == [-1 0 1; -1 0 1]
@test normalize_minmax(m1; bych = false) ==
      [-1.0 -0.6 -0.19999999999999996;
      0.19999999999999996 0.6000000000000001 1.0]
@test normalize_n(v1) == [0.0, 0.25, 0.5, 0.75, 1.0]
@test normalize_n(m1; bych = true) == [0.0 0.5 1.0; 0.0 0.5 1.0]
@test normalize_n(m1; bych = false) == [0.0 0.2 0.4; 0.6 0.8 1.0]
@test normalize_log(v1) == [
    1.0986122886681098,
    1.3862943611198906,
    1.6094379124341003,
    1.791759469228055,
    1.9459101490553132,
]
@test normalize_log(m1; bych = true) == [
    1.0986122886681098 1.3862943611198906 1.6094379124341003;
    2.1972245773362196 2.302585092994046 2.3978952727983707
]
@test normalize_log(m1; bych = false) == [
    1.0986122886681098 1.3862943611198906 1.6094379124341003;
    1.791759469228055 1.9459101490553132 2.0794415416798357
]
# rank-based inverse normal transform: Φ⁻¹(rank / (n + 1))
@test normalize_gauss(v1) ≈ norminv.((1:5) ./ 6)
@test normalize_gauss(v1)[1] ≈ -0.9674215661017014
@test normalize_gauss(m1; bych = true) ≈ [
    -0.6744897501960818 0.0 0.6744897501960818;
    -0.6744897501960818 0.0 0.6744897501960818
]
@test normalize_gauss(m1; bych = false) ≈ norminv.(m1 ./ 7)
@test normalize_log10(v1) == [
    0.47712125471966244,
    0.6020599913279624,
    0.6989700043360189,
    0.7781512503836436,
    0.8450980400142568,
]
@test normalize_log10(m1; bych = true) == [
    0.47712125471966244 0.6020599913279624 0.6989700043360189;
    0.9542425094393249 1.0 1.0413926851582251
]
@test normalize_log10(m1; bych = false) == [
    0.47712125471966244 0.6020599913279624 0.6989700043360189;
    0.7781512503836436 0.8450980400142568 0.9030899869919435
]
@test normalize_neglog(v1) == [
    -0.0,
    -0.6931471805599453,
    -1.0986122886681098,
    -1.3862943611198906,
    -1.6094379124341003,
]
@test normalize_neglog(m1; bych = true) == [
    -0.0 -0.6931471805599453 -1.0986122886681098;
    -1.3862943611198906 -1.6094379124341003 -1.791759469228055
]
@test normalize_neglog(m1; bych = false) == [
    -0.0 -0.6931471805599453 -1.0986122886681098;
    -1.3862943611198906 -1.6094379124341003 -1.791759469228055
]
@test normalize_neglog10(v1) == [
    -0.0,
    -0.3010299956639812,
    -0.47712125471966244,
    -0.6020599913279624,
    -0.6989700043360189,
]
@test normalize_neglog10(m1; bych = true) == [
    -0.0 -0.3010299956639812 -0.47712125471966244;
    -0.6020599913279624 -0.6989700043360189 -0.7781512503836436
]
@test normalize_neglog10(m1; bych = false) == [
    -0.0 -0.3010299956639812 -0.47712125471966244;
    -0.6020599913279624 -0.6989700043360189 -0.7781512503836436
]
@test normalize_neg(v1) == [-4, -3, -2, -1, 0]
@test normalize_neg(m1; bych = true) == [-2.0 -1.0 0.0; -2.0 -1.0 0.0]
@test normalize_neg(m1; bych = false) == [-5 -4 -3; -2 -1 0]
@test normalize_pos(v1) == [2, 3, 4, 5, 6]
@test normalize_pos(m1; bych = true) == [2.0 3.0 4.0; 8.0 9.0 10.0]
@test normalize_pos(m1; bych = false) == [2 3 4; 5 6 7]
@test normalize_perc(v1) == [0.0, 0.25, 0.5, 0.75, 1.0]
@test normalize_perc(m1; bych = true) == [0.0 0.5 1.0; 0.0 0.5 1.0]
@test normalize_perc(m1; bych = false) == [0.0 0.2 0.4; 0.6 0.8 1.0]
@test normalize_invroot(v1) ==
      [1.0, 0.7071067811865475, 0.5773502691896258, 0.5, 0.4472135954999579]
@test normalize_invroot(m1; bych = true) ==
      [1.0 0.7071067811865475 0.5773502691896258;
      0.5 0.4472135954999579 0.4082482904638631]
@test normalize_invroot(m1; bych = false) ==
      [1.0 0.7071067811865475 0.5773502691896258;
      0.5 0.4472135954999579 0.4082482904638631]
@test normalize_softmax(v1) ≈ [
    0.011656230956039607,
    0.03168492079612427,
    0.0861285444362687,
    0.23412165725273662,
    0.6364086465588308,
]
@test normalize_softmax(m1; bych = true) ≈ [
    0.0042697785452821095 0.011606461431184656 0.03154963320110001;
    0.08576079462509834 0.23312200962361299 0.6336913225737218
]
@test normalize_softmax(m1; bych = false) ≈ [
    0.0042697785452821095 0.011606461431184656 0.03154963320110001;
    0.08576079462509834 0.23312200962361299 0.6336913225737218
]
@test normalize_sigmoid(v1) == [
    0.7310585786300049,
    0.8807970779778823,
    0.9525741268224334,
    0.9820137900379085,
    0.9933071490757153,
]
@test normalize_sigmoid(m1; bych = true) == [
    0.7310585786300049 0.8807970779778823 0.9525741268224334;
    0.9820137900379085 0.9933071490757153 0.9975273768433653
]
@test normalize_sigmoid(m1; bych = false) == [
    0.7310585786300049 0.8807970779778823 0.9525741268224334;
    0.9820137900379085 0.9933071490757153 0.9975273768433653
]
@test normalize_mad(v1) == [
    -0.9098742077378683,
    -0.45493710386893416,
    0.0,
    0.45493710386893416,
    0.9098742077378683,
]
@test normalize_mad(m1; bych = true) == [
    -0.45493710386893416 0.0 0.45493710386893416;
    -0.45493710386893416 0.0 0.45493710386893416
]
@test normalize_mad(m1; bych = false) == [
    -0.7582285064482237 -0.4549371038689342 -0.15164570128964475;
    0.15164570128964475 0.4549371038689342 0.7582285064482237
]
@test normalize_rank(v1) == 1.0:5.0
@test normalize_rank(m1; bych = true) == [1.0 2.0 3.0; 1.0 2.0 3.0]
@test normalize_rank(m1; bych = false) == [1.0 2.0 3.0; 4.0 5.0 6.0]
@test normalize_fisher(v1) ==
      [-18.36840028483855, -0.5493061443340549, 0.0, 0.5493061443340549, 18.36840028483855]
@test normalize_fisher(m1; bych = true) == [
    -18.36840028483855 -0.6931471805599453 -0.20273255405408216;
    0.20273255405408216 0.6931471805599454 18.36840028483855
]
@test normalize_fisher(m1; bych = false) == [
    -18.36840028483855 -0.6931471805599453 -0.20273255405408216;
    0.20273255405408216 0.6931471805599454 18.36840028483855
]
e10_tmp = NeuroAnalyzer.normalize(e10; ch = "all", method = :zscore)
@test size(e10_tmp) == size(e10)

@info "Test: plinterpolate_channel()"
e10_tmp = deepcopy(e10)
e10_tmp.data[1, :, 1] = zeros(epoch_len(e10))
e10_int = plinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1)
@test e10_int.data[1, :, 1] != e10_tmp.data[1, :, 1]

@info "Test: remove_dc()"
@test NeuroAnalyzer._remove_dc(v1) == [-2.0, -1.0, 0.0, 1.0, 2.0]
e10_tmp = remove_dc(e10; ch = "all")
@test size(e10_tmp) == size(e10)

@info "Test: reference()"
e10_tmp = reference_ce(e10; ch = "Fp1")
@test size(e10_tmp) == size(e10)
e10_tmp = reference_ce(e10; ch = labels(e10)[1:5])
@test size(e10_tmp) == size(e10)
e10_tmp = reference_ce(e10; ch = "Fp1", med = true)
@test size(e10_tmp) == size(e10)
ch_eeg = get_channel(e10; ch = get_channel(e10; type = "eeg"))
fp1 = get_channel(e10; ch = "Fp1")[1]
e10_tmp = reference_ce(e10; ch = "Fp1")
@test e10_tmp.data[ch_eeg, :, :] ≈ e10.data[ch_eeg, :, :] .- e10.data[fp1:fp1, :, :]
@test all(iszero, e10_tmp.data[fp1, :, :])
@test labels(e10_tmp)[ch_eeg[1]] == labels(e10)[ch_eeg[1]] * "-Fp1"
@test e10_tmp.header.recording[:reference] == "common (Fp1)"
r5 = get_channel(e10; ch = labels(e10)[1:5])
e10_tmp = reference_ce(e10; ch = labels(e10)[1:5])
@test e10_tmp.data[ch_eeg, :, :] ≈
      e10.data[ch_eeg, :, :] .- NeuroAnalyzer.mean(e10.data[r5, :, :]; dims = 1)
e10_tmp = reference_ce(e10; ch = labels(e10)[1:5], med = true)
@test e10_tmp.data[ch_eeg, :, :] ≈
      e10.data[ch_eeg, :, :] .- NeuroAnalyzer.median(e10.data[r5, :, :]; dims = 1)
# regex selects all matching labels (Fp1, Fp2)
fp = get_channel(e10; ch = r"^Fp")
@test length(fp) == 2
@test reference_ce(e10; ch = r"^Fp").data ≈ reference_ce(e10; ch = labels(e10)[fp]).data

@info "Test: reference_a()"
e10_tmp = reference_a(e10)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_a(e10; med = true)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_a(e10; type = :c)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_a(e10; type = :i)
@test size(e10_tmp) == size(e10)

@info "Test: reference_m()"
edit_channel!(e10; ch = "A1", field = :label, value = "M1")
edit_channel!(e10; ch = "A2", field = :label, value = "M2")
e10_tmp = reference_m(e10)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_m(e10; med = true)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_m(e10; type = :c)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_m(e10; type = :i)
@test size(e10_tmp) == size(e10)

@info "Test: reference_avg()"
e10_tmp = reference_avg(e10)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_avg(e10; exclude_fpo = true)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_avg(e10; exclude_current = true)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_avg(e10; average = false)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_avg(e10; weighted = true)
@test size(e10_tmp) == size(e10)
# weighted CAR: locs row order does not matter; constant signal → 0
e10_shuf = deepcopy(e10)
e10_shuf.locs = e10_shuf.locs[reverse(1:size(e10_shuf.locs, 1)), :]
for md in (false, true)
    @test reference_avg(e10_shuf; weighted = true, med = md).data ≈
          reference_avg(e10; weighted = true, med = md).data
end
ch_eeg = get_channel(e10; ch = get_channel(e10; type = "eeg"))
e10_sim = deepcopy(e10)
e10_sim.data[ch_eeg, :, :] .= 5.0
@test all(isapprox.(reference_avg(e10_sim; weighted = true).data[ch_eeg, :, :], 0; atol = 1e-10))
# weights follow the polar (topographic) layout: out_i = x_i - Σ_j w_ij x_j / Σ_j w_ij, w_ij = 1/d_ij
e10_sim.data[ch_eeg, :, :] .= 0.0
e10_sim.data[ch_eeg[1], :, :] .= 1.0
out = reference_avg(e10_sim; weighted = true).data[ch_eeg, 1, 1]
locs = e10.locs[NeuroAnalyzer._find_bylabel(e10.locs, labels(e10)[ch_eeg]), :]
px = locs[!, :loc_radius] .* cosd.(locs[!, :loc_theta])
py = locs[!, :loc_radius] .* sind.(locs[!, :loc_theta])
n = length(ch_eeg)
w = [i == j ? 0.0 : 1 / hypot(px[i] - px[j], py[i] - py[j]) for i in 1:n, j in 1:n]
w[.!isfinite.(w)] .= 0.0
expected = [(i == 1 ? 1.0 : 0.0) - w[i, 1] / sum(w[i, :]) for i in 1:n]
@test out ≈ expected
e10_nopol = deepcopy(e10)
e10_nopol.locs[!, :loc_radius] .= 0.0
@test_throws ArgumentError reference_avg(e10_nopol; weighted = true)

@info "Test: reference_plap()"
e10_tmp = reference_plap(e10)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_plap(e10; weighted = true)
@test size(e10_tmp) == size(e10)

@info "Test: reference_slap()"
e10_tmp = reference_slap(e10)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_slap(e10; weighted = true)
@test size(e10_tmp) == size(e10)
# Laplacian: result does not depend on the order of rows in locs
e10_shuf = deepcopy(e10)
e10_shuf.locs = e10_shuf.locs[reverse(1:size(e10_shuf.locs, 1)), :]
for f in (reference_plap, reference_slap), wt in (false, true), md in (false, true)
    @test f(e10_shuf; weighted = wt, med = md).data ≈ f(e10; weighted = wt, med = md).data
end
# constant signal → 0; impulse at one channel → 1 there, 0 or -1/nn elsewhere
ch_eeg = get_channel(e10; ch = get_channel(e10; type = "eeg"))
e10_sim = deepcopy(e10)
e10_sim.data[ch_eeg, :, :] .= 5.0
for f in (reference_plap, reference_slap)
    @test all(isapprox.(f(e10_sim; nn = 3).data[ch_eeg, :, :], 0; atol = 1e-10))
end
e10_sim.data[ch_eeg, :, :] .= 0.0
e10_sim.data[ch_eeg[1], :, :] .= 1.0
for f in (reference_plap, reference_slap)
    out = f(e10_sim; nn = 3).data[ch_eeg, 1, 1]
    @test out[1] ≈ 1.0
    @test all(v -> isapprox(v, 0; atol = 1e-12) || isapprox(v, -1 / 3), out[2:end])
    @test any(v -> isapprox(v, -1 / 3), out[2:end])
end
@test_throws ArgumentError reference_slap(e10; nn = length(ch_eeg))

@info "Test: csd()"
g, h = gh(e10.locs)
@test size(g) == (23, 23)
@test size(h) == (23, 23)
e10_tmp = csd(e10)
@test size(e10_tmp) == size(e10)
# the surface Laplacian of a constant potential is zero
e_tmp = deepcopy(e10)
ch_eeg = get_channel(e_tmp; ch = get_channel(e_tmp; type = "eeg"))
e_tmp.data[ch_eeg, :, :] .= 5.0
@test maximum(abs.(csd(e_tmp).data[ch_eeg, :, :])) < 1e-6

@info "Test: standardize()"
m_s, sc = NeuroAnalyzer._standardize(a1)
@test length(sc) == 2
m_s, sc = NeuroAnalyzer.standardize(e10; ch = "all")
@test length(sc) == 10

@info "Test: taper()"
@test NeuroAnalyzer._taper(v1, t = v1) == [1, 4, 9, 16, 25]
e10_tmp = taper(e10; ch = "all", t = e10.data[1, :, 1])
@test size(e10_tmp) == size(e10)

@info "Test: tconv()"
@test round.(NeuroAnalyzer._tconv(v1, kernel = [0.2, 0.1, 0.2]), digits = 2) == [0.5, 1.0, 1.5, 2.0, 1.3]
@test round.(NeuroAnalyzer._tconv(a1, kernel = [0.2, 0.1, 0.2]), digits = 2) ==
      [0.3 0.5 0.3; 0.3 0.5 0.3;;; 0.3 0.5 0.3; 0.3 0.5 0.3]
e10_tmp = tconv(e10; ch = "all", kernel = [0.2, 0.1, 0.2])
@test size(e10_tmp) == size(e10)

@info "Test: wbp()"
@test length(NeuroAnalyzer._wbp(e10.data[1, :, 1], fs = 10, frq = 4)) == 2560
@test length(NeuroAnalyzer._wbp(e10.data[1, :, 1], fs = 10, frq = 4, pad = 100)) == 2560
# unit gain at frq
x = cos.(2π * 10 .* (0:2559) ./ 256)
@test NeuroAnalyzer._wbp(x, fs = 256, frq = 10)[500:2000] ≈ x[500:2000] atol = 1e-6
e10_tmp = wbp(e10; ch = "all", frq = 4)
@test size(e10_tmp) == size(e10)

@info "Test: intensity2od()"
n_tmp = intensity2od(nirs)
@test size(n_tmp.data) == (16, 9001, 1)  # 25.04 Hz resampled to 25 Hz on import

@info "Test: od2conc()"
n_tmp2 = od2conc(n_tmp)
@test size(n_tmp2.data) == (22, 9001, 1)

@info "Test: npl()"
e = average_epochs(e10)
npl!(e)
@test size(e.data) == (19, 2560, 11)

@info "Test: remove_pops()"
eeg_tmp, pl, ls, rs = remove_pops(eeg; ch = "all")
@test size(eeg_tmp) == size(eeg.data)
# the input signal is not modified
x = eeg.data[1, 1:2560, 1]
x0 = copy(x)
NeuroAnalyzer._remove_pops(x; r = 128)
@test x == x0

@info "Test: remove_powerline()"
e10_tmp = keep_epoch(e10; ep = 1)
remove_powerline!(e10_tmp; pl_frq = 50, ch = "Fp1")
@test size(e10_tmp) == (24, 2560, 1)

@info "Test: ica_decompose()"
ica_data = NeuroAnalyzer._ica_decompose(rand(10, 1000); n = 5)
@test size(ica_data.ic) == (5, 1000)
@test size(ica_data.ic_mw) == (10, 5)
ica_data = ica_decompose(eeg; ch = "all", n = 5, iter = 10)
@test size(ica_data.ic) == (5, 308480)
@test size(ica_data.ic_mw) == (24, 5)
@test length(ica_data.ic_var) == 5
# removing no component keeps the data (including channel means)
e_tmp = ica_remove(eeg; ch = "all", ic_idx = Int64[], ic = ica_data.ic, ic_mw = ica_data.ic_mw)
@test e_tmp.data == eeg.data
# keeping all components keeps the data
e_tmp = ica_remove(eeg; ch = "all", ic_idx = 1:5, keep = true, ic = ica_data.ic, ic_mw = ica_data.ic_mw)
@test e_tmp.data ≈ eeg.data

@info "Test: ica_remove()"
ica_data = ica_decompose(eeg; ch = "all", n = 5, iter = 10)
eeg_tmp =
    ica_remove(eeg; ch = "all", ic_idx = 1:2, ic = ica_data.ic, ic_mw = ica_data.ic_mw)
@test size(eeg_tmp) == size(eeg)
eeg_tmp = ica_remove(
    eeg;
    ch = "all",
    ic_idx = 1:2,
    ic = ica_data.ic,
    ic_mw = ica_data.ic_mw,
    keep = true,
)
@test size(eeg_tmp) == size(eeg)

@info "Test: pca_decompose()"
pc, pcv, pcm, pc_model = NeuroAnalyzer._pca_decompose(rand(4, 4, 2); n = 2)
@test size(pc) == (2, 4, 2)
@test size(pcv) == (2, 2)
@test length(pcm) == 4
pc, pcv, pcm, _ = pca_decompose(e10; ch = "all", n = 4)
@test size(pc) == (4, 2560, 10)
@test size(pcv) == (4, 10)
@test length(pcm) == 24

@info "Test: pca_reconstruct()"
pc, pcv, pcm, pc_model = NeuroAnalyzer._pca_decompose(rand(4, 4, 2); n = 2)
s = NeuroAnalyzer._pca_reconstruct(rand(4, 4, 2); pc = pc, pc_model = pc_model)
@test size(s) == (4, 4, 2)
pc, pcv, pcm, pc_model = pca_decompose(e10; ch = "all", n = 4)
e10_rec = pca_reconstruct(e10; ch = "all", pc = pc, pc_model = pc_model)
@test size(e10_rec.data) == size(e10)
# all components: exact reconstruction of every epoch
x = rand(4, 50, 3)
pc, _, _, pc_model = NeuroAnalyzer._pca_decompose(x; n = 4)
@test NeuroAnalyzer._pca_reconstruct(x; pc = pc, pc_model = pc_model) ≈ x

@info "Test: reference_custom()"
e10_tmp = reference_custom(e10)
@test size(e10_tmp) == (23, 2560, 10)

@info "Test: normpower()"
@test NeuroAnalyzer._rms(NeuroAnalyzer._normpower(1:10)) ≈ 1
@test NeuroAnalyzer._normpower(zeros(5)) == zeros(5)
@test size(normpower(e10, ch = "all")) == size(e10)

@info "Test: sort_epochs()"
e10_erp = average_epochs(e10)
e = sort_epochs(e10_erp; s = collect((nepochs(e10_erp)):-1:2))
@test size(e.data) == (19, 2560, 11)

@info "Test: denoise_cwd())"
s = NeuroAnalyzer._denoise_cwd(rand(100); fs = 10, nf = 2)
@test length(s) == 100
e10_tmp = denoise_cwd(e10; ch = "all", nf = 50)
@test size(e10_tmp) == size(e10)

@info "Test: artrem_cwd())"
e10_tmp = artrem_cwd(e10; ch = "Fp1", ep = 1, tseg = (0.2, 0.4), fseg = (10, 12.5))
@test size(e10_tmp) == size(e10)

@info "Test: emd()"
@test size(NeuroAnalyzer._emd(e10.data[1, :, 1], e10.epoch_time, epsilon = 0.3), 1) in [8, 9]
@test size(emd(e10, ch = "Fp1", ep = 1, epsilon = 0.3), 1) in [8, 9]

@info "Test: process regressions"
# resample: time axes use the new sampling rate; epoch offset is kept; markers unchanged
r = resample(eeg; new_sr = 128)
@test sr(r) == 128 && r.time_pts[2] == 1 / 128 && length(r.time_pts) == size(r.data, 2)
@test r.markers == eeg.markers
e_off = deepcopy(e10); e_off.epoch_time .-= 1.0
@test resample(e_off; new_sr = 512).epoch_time[1] == -1.0
@test resample(e10; new_sr = sr(e10)) !== e10
@test_throws ArgumentError downsample(e10; new_sr = 512)
@test_throws ArgumentError upsample(e10; new_sr = 128)
# average_epochs: markers shifted by one epoch (in seconds); in-place version keeps locs in sync
e_avg = average_epochs(e10)
isempty(e10.markers) || @test e_avg.markers[!, :start] ≈ e10.markers[!, :start] .+ 10
e_tmp = deepcopy(e10); average_epochs!(e_tmp)
@test e_tmp.locs == e_avg.locs
# sort_epochs: s must be a permutation
@test_throws ArgumentError sort_epochs(e_avg; s = fill(2, nepochs(e_avg) - 1))
# plinterpolate_channel: any channel, all channels kept
e_pl = plinterpolate_channel(e10; ch = "Cz", ep = 1)
@test labels(e_pl) == labels(e10)
cz = findfirst(==("Cz"), labels(e10))
@test e_pl.data[cz, :, 1] != e10.data[cz, :, 1] && e_pl.data[cz, :, 2] == e10.data[cz, :, 2]
# lrinterpolate_channel: deterministic
@test lrinterpolate_channel(e10; ch = "Cz", ep = 1).data == lrinterpolate_channel(e10; ch = "Cz", ep = 1).data
# standardize!, tconv! with a complex kernel
e_tmp = deepcopy(e10)
@test standardize!(e_tmp; ch = "Fp1") isa Vector
@test tconv!(deepcopy(e10); ch = "Fp1", kernel = ComplexF64[0.2 + 0.1im, 0.5, 0.2]) isa Array{ComplexF64, 3}
# normalize: :fisher available; rank keeps shape; softmax does not overflow
@test size(NeuroAnalyzer._normalize(rand(3, 4); method = :fisher)) == (3, 4)
@test size(normalize_rank(ones(2, 3))) == (2, 3)
@test sum(normalize_softmax([1000.0, 1001.0])) ≈ 1
# weighted median with equal weights is the median
X = rand(6, 50)
@test NeuroAnalyzer._wmedian_cols(X, ones(6)) ≈ vec(NeuroAnalyzer.Statistics.median(X; dims = 1))

@info "Test: maxwell_filter()"
# reference: MaxFilter output (MNE-Python testing data); per-channel SNR of the MEG channels,
# ‖ref‖ / ‖out - ref‖ (as MNE-Python's assert_meg_snr); maxwell_filter() reproduces
# MNE-Python's output, which deviates from MaxFilter by the same amount
sss_path = joinpath(testfiles_path, "SSS")
function _mf_snr(out, ref_file)
    ref = import_fiff(joinpath(sss_path, ref_file))
    rr = ref.header.recording
    r = Float64[]
    for i in sort(vcat(rr[:magnetometers], rr[:gradiometers]))
        j = findfirst(==(labels(ref)[i]), labels(out))
        n = min(size(out.data, 2), size(ref.data, 2))
        x = out.data[j, 1:n, 1]
        y = ref.data[i, 1:n, 1]
        push!(r, NeuroAnalyzer.norm(y) / NeuroAnalyzer.norm(x .- y))
    end
    return minimum(r), NeuroAnalyzer.Statistics.median(r)
end
raw = import_fiff(joinpath(sss_path, "test_move_anon_raw.fif"); allow_maxshield = true)
# the references are 1 s long
raw = NeuroAnalyzer.trim(raw; seg = (0, 1), keep = true)
erm = import_fiff(joinpath(sss_path, "test_move_anon_erm_raw.fif"); allow_maxshield = true)
o = [0.0, 0.0, 0.04]
e_o = (; coord_frame = :meg, origin = [0.0, 0.013, -0.006])
s = maxwell_filter(raw; origin = o, regularize = :none)
mn, md = _mf_snr(s, "test_move_anon_stdOrigin_raw_sss.fif")
@test mn > 10 && md > 100
@test s.header.recording[:sss] && !s.header.recording[:maxshield]
@test_throws ArgumentError maxwell_filter(s)
mn, md = _mf_snr(maxwell_filter(raw; origin = o), "test_move_anon_regIn_raw_sss.fif")
@test mn > 20 && md > 200
# fine calibration: channel names "MEG 113" (data) and "MEG0113" (file) must match
cal = read_fine_calibration(joinpath(sss_path, "sss_cal_3053.dat"))
@test haskey(cal, NeuroAnalyzer._sss_name("MEG 113")) && NeuroAnalyzer._sss_name("MEG0113") == "MEG113"
mn, md = _mf_snr(
    maxwell_filter(raw; origin = o, regularize = :none, calibration = joinpath(sss_path, "sss_cal_3053.dat")),
    "test_move_anon_fineCal_raw_sss.fif",
)
@test md > 50
mn, md = _mf_snr(
    maxwell_filter(erm; regularize = :none, cross_talk = joinpath(sss_path, "ct_sparse.fif"), e_o...),
    "test_move_anon_erm_crossTalk_raw_sss.fif",
)
@test mn > 40 && md > 200
# tSSS with MaxFilter-like (non-overlapping) buffers
mn, md = _mf_snr(
    maxwell_filter(erm; regularize = :none, st_duration = 1, st_overlap = false, e_o...),
    "test_move_anon_erm_st1_raw_sss.fif",
)
@test mn > 100 && md > 150
# movement compensation with MaxFilter head positions (times from the start of the acquisition)
hp = [parse.(Float64, split(l)) for l in readlines(joinpath(sss_path, "test_move_anon_raw.pos"))[2:end] if !isempty(strip(l))]
hp = reduce(vcat, permutedims.(hp))
hp[:, 1] .-= raw.header.recording[:first_time]
mn, md = _mf_snr(
    maxwell_filter(raw; origin = o, regularize = :none, head_pos = hp, mc_interp = :zero),
    "test_move_anon_movecomp_raw_sss.fif",
)
@test mn > 10 && md > 100
@test maxwell_filter(raw; origin = o, regularize = :none, head_pos = hp).header.recording[:sss_info][:movement_compensation]
# destination: translation without rotation
mn, md = _mf_snr(maxwell_filter(raw; origin = o, regularize = :none, destination = o), "test_move_anon_transDefault_raw_sss.fif")
@test mn > 5 && md > 100
# internal helpers: Hann window is COLA at 50% overlap, quaternion average of one position
@test all(NeuroAnalyzer._sss_hann(8)[1:4] .+ NeuroAnalyzer._sss_hann(8)[5:8] .≈ 1)
R = NeuroAnalyzer._vec_rotation([1.0, 0.0, 0.0], [0.0, 0.0, -1.0])
@test R * [1.0, 0.0, 0.0] ≈ [0.0, 0.0, -1.0] && R' * R ≈ NeuroAnalyzer.I(3)

@info "Test: chpi_head_positions()"
c = import_fiff(joinpath(sss_path, "chpi5_raw.fif"))
c = NeuroAnalyzer.trim(c; seg = (0, 3), keep = true)
hp = chpi_head_positions(c)
@test size(hp, 2) == 10 && all(hp[:, 8] .> 0.98)
ref = [parse.(Float64, split(l)) for l in readlines(joinpath(sss_path, "chpi5_raw_mc.pos"))[2:end] if !isempty(strip(l))]
ref = reduce(vcat, permutedims.(ref))
ref[:, 1] .-= c.header.recording[:first_time]
for r in eachrow(ref)
    0 < r[1] < 3 || continue
    k = argmin(abs.(hp[:, 1] .- r[1]))
    @test NeuroAnalyzer.norm(hp[k, 5:7] .- r[5:7]) < 0.002
end

@info "Test: apply_ssp_projectors()"
m_ssp = import_fiff(joinpath(testfiles_path, "test_raw.fif"))
p_ssp = generate_ssp_projectors(m_ssp)
@test size(p_ssp.U, 2) == 3 && p_ssp.ssp_projectors * p_ssp.ssp_projectors ≈ p_ssp.ssp_projectors
m_p = apply_ssp_projectors(m_ssp)
@test apply_ssp_projectors(m_p).data ≈ m_p.data

true
