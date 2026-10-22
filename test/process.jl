using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
nirs = import_nirs(joinpath(testfiles_path, "nirs", "fnirs-test-nirs.nirs"))
meg = import_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fif"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)
epoch!(meg; ep_len = 2)
keep_epoch!(meg; ep = 1:10)
load_locs!(e10; file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"))
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
r = NeuroAnalyzer._artrem_cwd(x, collect(0:1023) ./ 256; fs = 256, tseg = (1, 2), fseg = (60, 70))
@test NeuroAnalyzer.norm(r .- x) / NeuroAnalyzer.norm(x) < 0.02
# artrem_cwd(): the artifact is reconstructed and subtracted, the signal far from it is not changed
tb = collect(0:1023) ./ 256
xb = x .+ 5 .* exp.(-(tb .- 2) .^ 2 ./ (2 * 0.1^2))
r = NeuroAnalyzer._artrem_cwd(xb, tb; fs = 256, tseg = (1.5, 2.5), fseg = (0, 5))
far = (tb .< 1) .| (tb .> 3)
@test maximum(abs.(r[far] .- xb[far])) < 0.01
@test maximum(abs.(r .- x)) < 0.5
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
# deterministic (the noise model is the expected white-noise spectrum, not a random realization)
xw = randn(NeuroAnalyzer.Random.MersenneTwister(1), 3, 256, 2)
@test NeuroAnalyzer._denoise_wien(xw) == NeuroAnalyzer._denoise_wien(xw)
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
# band clamped to Nyquist (fs = 256 Hz): extracted with a high-pass filter
@test bf[findfirst(==(:gamma), bn)] == (30.0, 127.9)
s, bn, bf = bpsplit(e10; ch = "Fp1", band = [:alpha, :gamma], fprototype = :butterworth, order = 4)
@test bn == [:alpha, :gamma]
@test size(s) == (2, 1, 2560, 10)
# 10 Hz + 20 Hz sines: each band keeps its own component (std of a unit sine = 1/√2)
e_tmp = deepcopy(e10)
e_tmp.data[1, :, :] .= sin.(2π * 10 .* e_tmp.epoch_time) .+ sin.(2π * 20 .* e_tmp.epoch_time)
std_mid(x) = (m = x[600:1960]; sqrt(sum((m .- sum(m) / length(m)) .^ 2) / (length(m) - 1)))
for (p, kw) in ((:fir, (;)), (:firls, (;)), (:remez, (;)), (:butterworth, (order = 6,)))
    local s = bpsplit(e_tmp; ch = "Fp1", band = [:alpha, :beta, :gamma], fprototype = p, kw...).s
    @test std_mid(s[1, 1, :, 1]) ≈ 1 / sqrt(2) atol = 0.01
    @test std_mid(s[2, 1, :, 1]) ≈ 1 / sqrt(2) atol = 0.01
    @test std_mid(s[3, 1, :, 1]) < 0.01
end
@test_throws ArgumentError bpsplit(e10; ch = "Fp1", fprototype = :iirnotch)
@test_throws ArgumentError bpsplit(e10; ch = "Fp1", band = [:foo])
@test_throws ArgumentError bpsplit(e10; ch = "Fp1", fprototype = :butterworth)

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
    order = 4,
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
    order = 4,
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
    bw = 1,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :hp,
    cutoff = 1,
    bw = 0.5,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :bs,
    cutoff = (49, 51),
    bw = 1,
)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(
    e10;
    ch = "all",
    fprototype = :fir,
    ftype = :bp,
    cutoff = (49, 51),
    bw = 1,
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
# segment of epoch 1 only, trained on the rest of the epoch; the samples outside the
# segment are not changed
e10_seg = deepcopy(e10)
r_seg = NeuroAnalyzer._epoch_seg(e10_seg, (2, 4))
e10_seg.data[1, r_seg, 1] .= 0
e10_int = lrinterpolate_channel(e10_seg; ch = "Fp1", ep = 1, ep_ref = Int64[], seg = (2, 4))
@test all(!iszero, e10_int.data[1, r_seg, 1])
@test e10_int.data[1, setdiff(1:epoch_len(e10), r_seg), 1] == e10.data[1, setdiff(1:epoch_len(e10), r_seg), 1]
@test_throws ArgumentError lrinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1, ep_ref = Int64[])
@test_throws ArgumentError lrinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1, seg = (4, 2))
@test_throws ArgumentError lrinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1, ep_ref = 1)
# bad channels are not used as predictors
e10_bad = deepcopy(e10_tmp)
e10_bad.data[2, :, :] .= 1e6 .* randn(size(e10_bad.data, 2), size(e10_bad.data, 3))
e10_bad.header.recording[:bad_channel][2] = true
@test lrinterpolate_channel(e10_bad; ch = "Fp1", ep = 1).data[1, :, 1] ≈ lrinterpolate_channel(delete_channel(e10_tmp; ch = "Fp2"); ch = "Fp1", ep = 1).data[1, :, 1]

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
# (x − median) / (1.4826 × MAD): v1 median 3, MAD 1; rows of m1 median 2/5, MAD 1; m1 median 3.5, MAD 1.5
@test normalize_mad(v1) ≈ (v1 .- 3) ./ 1.4826
@test normalize_mad(m1; bych = true) ≈ [-1 0 1; -1 0 1] ./ 1.4826
@test normalize_mad(m1; bych = false) ≈ (m1 .- 3.5) ./ (1.4826 * 1.5)
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
# all epochs by default; segment only
@test plinterpolate_channel(e10_tmp; ch = "Fp1").data[1, :, :] != e10_tmp.data[1, :, :]
e10_int = plinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1, seg = (2, 4))
r_seg = NeuroAnalyzer._epoch_seg(e10_tmp, (2, 4))
@test all(!iszero, e10_int.data[1, r_seg, 1])
@test all(iszero, e10_int.data[1, setdiff(1:epoch_len(e10), r_seg), 1])
@test_throws ArgumentError plinterpolate_channel(e10_tmp; ch = "Fp1", seg = (0, 100))

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
    local out = f(e10_sim; nn = 3).data[ch_eeg, 1, 1]
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
# synthetic: 10 Hz sine + noise (σ = 2), a positive pop (200, τ = 40 samples) and a negative pop (-150, τ = 20 samples)
fs_p = 256
n_p = 10 * fs_p
t_p = (0:(n_p - 1)) ./ fs_p
clean_p = 10 .* sin.(2π * 10 .* t_p) .+ 2 .* randn(NeuroAnalyzer.Random.Xoshiro(1), n_p)
pops_p = zeros(n_p)
pops_p[1200:end] .+= 200 .* exp.(-(0:(n_p - 1200)) ./ 40)
pops_p[1900:end] .-= 150 .* exp.(-(0:(n_p - 1900)) ./ 20)
x = clean_p .+ pops_p
x0 = copy(x)
p = NeuroAnalyzer._remove_pops(x; r = 128)
@test x == x0
@test p.pop_loc == [1200, 1900]
@test isapprox(p.amp, [200, -150]; rtol = 0.05)
@test isapprox(p.tau, [40, 20]; rtol = 0.1)
# the pops are removed, the underlying signal is preserved (residual below the noise level)
@test maximum(abs.(p.s .- clean_p)) < 10
@test sqrt(sum(abs2, p.s .- clean_p) / n_p) < 1
# detection only
p = NeuroAnalyzer._remove_pops(x; r = 128, repair = false)
@test p.s == x && p.pop_loc == [1200, 1900]
# no pops in a clean signal; a single-sample spike is not a pop
p = NeuroAnalyzer._remove_pops(clean_p; r = 128)
@test isempty(p.pop_loc) && p.s == clean_p
x_sp = copy(clean_p)
x_sp[1500] += 300
@test isempty(NeuroAnalyzer._remove_pops(x_sp; r = 128).pop_loc)
# exponential fit
@test all(isapprox.(values(NeuroAnalyzer._pop_fit(50 .* exp.(-(0:99) ./ 12))), (50, 12); rtol = 0.02))
# a step that does not recover is not a pop: not repaired (a warning is shown)
x_st = copy(clean_p)
x_st[1200:end] .+= 300
p = NeuroAnalyzer._remove_pops(x_st; r = 128)
@test isempty(p.pop_loc) && p.s == x_st
# acceptance: a small pop (relative to the signal SD) or a shape that does not follow the model is rejected
@test isempty(NeuroAnalyzer._remove_pops(x; r = 128, a_min = 100).pop_loc)
x_bl = copy(clean_p)
x_bl[1200:1327] .+= 300 .* sin.(π .* (0:127) ./ 127) .^ 8   # a fast bump that returns (not a step with exponential recovery)
@test isempty(NeuroAnalyzer._remove_pops(x_bl; r = 128).pop_loc)
# a pop at the very beginning of the signal (e.g. a filter transient) is not analyzed
x_b = copy(clean_p)
x_b[2:end] .+= 300 .* exp.(-(0:(n_p - 2)) ./ 40)
@test isempty(NeuroAnalyzer._remove_pops(x_b; r = 128).pop_loc)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x; r = 128, a_min = -1)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x; r = 128, r2_min = 2)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x; r = 4)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x; r = 128, z = 0)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x[1:100]; r = 128)
# NEURO object: channel 1 with the pops, channel 2 clean
o_p = create_object(; data_type = "eeg")
o_p.header.recording[:sampling_rate] = fs_p
add_channel!(o_p; data = reshape(vcat(x', clean_p'), 2, n_p, 1), label = ["P1", "P2"], type = ["eeg", "eeg"], unit = ["μV", "μV"])
create_time!(o_p; fs = fs_p)
o_r, pops = remove_pops(o_p; ch = "all", window = 5)
@test pops[!, :channel] == ["P1", "P1"]
@test pops[!, :sample] == [1200, 1900]
@test pops[!, :time] ≈ [1199, 1899] ./ fs_p
@test isapprox(pops[!, :tau], [40, 20] ./ fs_p; rtol = 0.1)
@test o_r.data[2, :, 1] == clean_p
@test sqrt(sum(abs2, o_r.data[1, :, 1] .- clean_p) / n_p) < 1
@test o_p.data[1, :, 1] == x
# detection only: same pops (amplitudes of later pops differ: measured from the unrepaired level)
@test remove_pops(o_p; ch = "all", repair = false)[!, [:channel, :sample]] == pops[!, [:channel, :sample]]
@test_throws ArgumentError remove_pops(o_p; ch = "all", window = 0.5)
@test_throws ArgumentError remove_pops(o_p; ch = "all", window = 20)
@test_throws ArgumentError remove_pops(epoch(o_p; ep_len = 5); ch = "all")

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
sss_path = joinpath(testfiles_path, "meg", "SSS")
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
# fine calibration: channel names "MEG 0113" (data) and "MEG0113" (file) must match
cal = read_fine_calibration(joinpath(sss_path, "sss_cal_3053.dat"))
@test haskey(cal, NeuroAnalyzer._sss_name("MEG 0113")) && NeuroAnalyzer._sss_name("MEG0113") == "MEG113"
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
    local k = argmin(abs.(hp[:, 1] .- r[1]))
    @test NeuroAnalyzer.norm(hp[k, 5:7] .- r[5:7]) < 0.002
end

@info "Test: apply_ssp_projectors()"
m_ssp = import_fiff(joinpath(testfiles_path, "meg", "test_raw.fif"))
p_ssp = generate_ssp_projectors(m_ssp)
@test size(p_ssp.U, 2) == 3 && p_ssp.ssp_projectors * p_ssp.ssp_projectors ≈ p_ssp.ssp_projectors
m_p = apply_ssp_projectors(m_ssp)
@test apply_ssp_projectors(m_p).data ≈ m_p.data


##
## values: SciPy 1.17 / MNE-Python 1.13.2 references, synthetic signals, formulas
##

eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
eb = epoch(eeg; ep_len = 10)
keep_epoch!(eb; ep = 1:10)
load_locs!(eb; file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"))
fs = 256
tt = (0:(4 * fs - 1)) ./ fs
xm = sin.(2π * 5 .* tt) .+ 0.5 .* sin.(2π * 60 .* tt) .+ 0.3 .* cos.(2π * 90 .* tt)
# amplitude of a sinusoid at f in the middle of a signal (least squares)
sine_amp(s, f) = (mid = 257:768; A = [sin.(2π * f .* tt[mid]) cos.(2π * f .* tt[mid])]; hypot((A \ s[mid])...))

@info "Test: filter_create() / _filter_apply() vs SciPy / MNE"
# MNE filter_data(method = "iir", iir_params = (order = 4, ftype = "butter", output = "sos"), phase = "zero")
f = filter_create(; fprototype = :butterworth, ftype = :lp, cutoff = 30, fs = fs, order = 4)
@test NeuroAnalyzer._filter_apply(xm; flt = f)[[401, 513, 601]] ≈ [-0.9244172619059996, 6.06058617302202e-07, -0.98116451395666] atol = 1e-12
f = filter_create(; fprototype = :butterworth, ftype = :hp, cutoff = 40, fs = fs, order = 4)
@test NeuroAnalyzer._filter_apply(xm; flt = f)[[401, 513, 601]] ≈ [-0.704915844879155, 0.2999917628471555, -0.07129865301596978] atol = 1e-12
# scipy.signal.firwin(101, 30, window = "hamming", fs = 256)
f = filter_create(; fprototype = :fir, ftype = :lp, cutoff = 30, fs = fs, order = 101)
@test length(f) == 101
@test f[1:3] ≈ [-0.00039415085773660165, -0.0005255662713842926, -0.0003925985049683375] atol = 1e-15
@test sum(f) ≈ 1.0
# scipy.signal.remez(101, [0, 28, 32, 128], [1, 0], fs = 256)
f = filter_create(; fprototype = :remez, ftype = :lp, cutoff = 30, fs = fs, order = 101, bw = 4)
@test f[1:3] ≈ [-0.0010783991016224733, -0.011648688737594769, -0.0011394905586340003] atol = 1e-12
# scipy.signal.iirnotch(50, 25, fs = 256)
f = filter_create(; fprototype = :iirnotch, cutoff = 50, fs = fs, bw = 2)
@test [f.b0, f.b1, f.b2] ≈ [0.9760395733504627, -0.6576356575420848, 0.9760395733504627]
@test [f.a1, f.a2] ≈ [-0.6576356575420848, 0.9520791467009253]
# |H(f)| of scipy.signal.cheby1(4, 0.5, 30), cheby2(4, 20, 30), ellip(4, 0.5, 20, 30) at 5, 20, 29, 31, 45, 80 Hz
fr = [5.0, 20, 29, 31, 45, 80]
for (p, ref) in ((:chebyshev1, [0.9629080387446526, 0.9490149558771923, 0.9868240257223284, 0.867296802722799, 0.08728069796786471, 0.0016906547989606268]),
                 (:chebyshev2, [0.9999996629508234, 0.9383994834987681, 0.1648264095457529, 0.04838851547390617, 0.09073270752540165, 0.050573680748433526]),
                 (:elliptic, [0.9541360941244506, 0.9824598868041596, 0.9999988024626586, 0.7061632202106207, 0.07599550817942868, 0.06624646323565964]))
    local f = filter_create(; fprototype = p, ftype = :lp, cutoff = 30, fs = fs, order = 4, rp = 0.5, rs = 20)
    @test abs.(NeuroAnalyzer.DSP.freqresp(f, 2π .* fr ./ fs)) ≈ ref atol = 1e-10
end
@test_throws ArgumentError filter_create(; fprototype = :butterworth, ftype = :lp, cutoff = 30, fs = fs)
@test_throws ArgumentError filter_create(; fprototype = :firls, ftype = :lp, cutoff = 30, fs = fs, order = 101)
@test_throws ArgumentError filter_create(; fprototype = :fir, ftype = :hp, cutoff = 30, fs = fs, order = 100)
@test_throws ArgumentError filter_create(; fprototype = :iirnotch, cutoff = 127, fs = fs, bw = 4)

@info "Test: filter() responses (all prototypes and types)"
# synthetic mixture of 5, 30 and 60 Hz sines: pass bands kept, stop bands removed (two-pass)
xs = sin.(2π * 5 .* tt) .+ sin.(2π * 30 .* tt) .+ sin.(2π * 60 .* tt)
designs = Dict(:lp => (40, [5], [60]), :hp => (20, [30, 60], [5]), :bp => ((20, 40), [30], [5, 60]), :bs => ((20, 40), [5, 60], [30]))
for p in (:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :firls, :remez), (ft, (c, keep, stop)) in designs
    kw = p in (:butterworth, :chebyshev1, :chebyshev2, :elliptic) ? (order = 6, rp = 0.5, rs = 40) :
         p === :fir ? (order = 151,) : (order = 151, bw = 6)
    p in (:butterworth,) && (kw = (order = 6,))
    local f = filter_create(; fprototype = p, ftype = ft, cutoff = c, fs = fs, kw...)
    local y = NeuroAnalyzer._filter_apply(xs; flt = f)
    # Chebyshev I / elliptic: pass-band ripple rp = 0.5 dB (twice for two-pass filtering); Chebyshev II:
    # the cutoff is the stop-band edge, so frequencies near it are in the transition band
    lo = p in (:chebyshev1, :elliptic) ? 10^(-2 * 0.5 / 20) - 0.02 : 0.9
    for k in keep
        p === :chebyshev2 && ft === :hp && k == 30 && continue
        @test lo <= sine_amp(y, k) <= 1.02
    end
    for s in stop
        @test sine_amp(y, s) < 0.05
    end
end

@info "Test: two-pass FIR filtering of long signals"
# DSP.filtfilt(b, x) filters in place and is corrupted when filt! uses FFT overlap-save
# reference: forward and reverse convolution (the middle is unaffected by the edge extension)
xl = randn(NeuroAnalyzer.Random.MersenneTwister(1), 100_000) .+ 100
f = filter_create(; fprototype = :fir, ftype = :hp, cutoff = 0.5, fs = fs, bw = 0.5)
fc(v) = NeuroAnalyzer.DSP.conv(v, f)[((length(f) - 1) ÷ 2 + 1):((length(f) - 1) ÷ 2 + length(v))]
yl = NeuroAnalyzer._filter_apply(xl; flt = f)
@test yl[20_000:80_000] ≈ reverse(fc(reverse(fc(xl))))[20_000:80_000] atol = 1e-9
@test NeuroAnalyzer._fir_filtfilt(f[1:31], xl) ≈ NeuroAnalyzer.DSP.filtfilt(f[1:31], xl)
@test length(NeuroAnalyzer._fir_filtfilt(f, xl[1:500])) == 500

@info "Test: filter() / filter_apply() / filter!() on objects"
f = filter_create(; fprototype = :butterworth, ftype = :lp, cutoff = 30, fs = fs, order = 4)
e_f = filter_apply(eb; ch = "Fp1", flt = f, report = false)
@test e_f.data[1, :, 3] ≈ NeuroAnalyzer._filter_apply(eb.data[1, :, 3]; flt = f)
@test e_f.data[2, :, :] == eb.data[2, :, :]
e_tmp = deepcopy(eb)
filter_apply!(e_tmp; ch = "Fp1", flt = f, report = false)
@test e_tmp.data == e_f.data
@test NeuroAnalyzer._filter_one(eb, 1, 3, f, :twopass) ≈ e_f.data[1, :, 3]
@test NeuroAnalyzer._filter_one(eb, 1, 3, f, :onepass) ≈ NeuroAnalyzer.DSP.filt(f, eb.data[1, :, 3])
e_f2 = NeuroAnalyzer.filter(eb; ch = "Fp1", fprototype = :butterworth, ftype = :lp, cutoff = 30, order = 4, report = false)
@test e_f2.data ≈ e_f.data
e_tmp = deepcopy(eb)
NeuroAnalyzer.filter!(e_tmp; ch = "Fp1", fprototype = :butterworth, ftype = :lp, cutoff = 30, order = 4, report = false)
@test e_tmp.data ≈ e_f.data

@info "Test: filter_order() / filter_report() / filter helpers"
# FIR: N = ⌈k fs / bw⌉ (Hamming k = 3.3), odd; FIRLS/Remez: harris' rule ⌈A/22 · fs / bw⌉
@test filter_order(; fprototype = :fir, fs = 256, bw = 2) == 423
@test filter_order(; fprototype = :fir, fs = 256, bw = 2, window = :blackman) == 705
@test filter_order(; fprototype = :firls, fs = 256, bw = 2) == 309
@test filter_order(; fprototype = :remez, fs = 256, bw = 2, rs = 44) == 257
@test_throws ArgumentError filter_order(; fprototype = :butterworth, fs = 256, bw = 2)
@test_throws ArgumentError filter_order(; fprototype = :fir, fs = 256, bw = 0)
f = filter_create(; fprototype = :butterworth, ftype = :lp, cutoff = 30, fs = fs, order = 4)
rep = filter_report(f; fs = fs, ftype = :lp, cutoff = 30, bw = 10, fcheck = [5.0, 30.0, 60.0], verbose = false)
# two-pass Butterworth: -6 dB at the cutoff
@test rep.gain_at[30.0] ≈ -6.02 atol = 0.05
@test rep.gain_at[5.0] ≈ 0 atol = 0.01
@test rep.gain_at[60.0] < -40
@test rep.summary isa String
@test_throws ArgumentError filter_report(f; fs = fs, fcheck = [200.0], verbose = false)
@test NeuroAnalyzer._fir_window(:hamming, 11) == NeuroAnalyzer.DSP.hamming(11)
@test NeuroAnalyzer._fir_window(:rect, 4) == ones(4)
@test NeuroAnalyzer._transition_edges(:lp, 30, 4) == [(28.0, 32.0)]
@test NeuroAnalyzer._transition_edges(:bp, (10, 20), 2) == [(9.0, 11.0), (19.0, 21.0)]
@test NeuroAnalyzer._check_edges(:lp, 30, 4, 128; strict = true) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_edges(:lp, 126, 8, 128; strict = true)

@info "Test: csd() vs MNE"
# MNE compute_current_source_density(sphere = (0, 0, 0, 1), lambda2 = 1e-5, stiffness = 4,
# n_legendre_terms = 8), EEG channels on the unit sphere; V/m² × 10⁶ (μV)
e2 = keep_epoch(eb; ep = 1:2)
c = csd(e2)
@test c.data[1, [1, 2, 101, 512], 1] ≈ [-1692.0568207927302, -1775.5717618275437, 659.987805703435, 603.4875914557678] rtol = 1e-9
@test c.data[10, [1, 2, 101, 512], 1] ≈ [-1962.2371848290109, -545.3811613543534, -136.083974881286, -1501.146212469763] rtol = 1e-9
@test c.header.recording[:data_type] == "csd"
@test c.header.recording[:unit][1] == "μV/m²"
# electrode positions are projected onto the unit sphere: the radius does not matter
e_tmp = deepcopy(e2)
e_tmp.locs[!, :loc_x] .*= 2; e_tmp.locs[!, :loc_y] .*= 2; e_tmp.locs[!, :loc_z] .*= 2
@test csd(e_tmp).data ≈ c.data
G, H = gh(e2.locs)
@test G ≈ G' && H ≈ H'
@test_throws ArgumentError csd(e2; m = 1)

@info "Test: reference values"
ch_eeg = get_channel(eb; ch = get_channel(eb; type = "eeg"))
e_r = reference_avg(eb)
@test e_r.data[ch_eeg, :, :] ≈ eb.data[ch_eeg, :, :] .- sum(eb.data[ch_eeg, :, :]; dims = 1) ./ length(ch_eeg)
# MNE set_eeg_reference("average"): Fp1, sample 101 of the first 2 s
@test e_r.data[1, 101, 1] ≈ 27.5266006673819 rtol = 1e-10
@test e_r.header.recording[:reference] == "average"
a1 = get_channel(eb; ch = "A1")[1]
a2 = get_channel(eb; ch = "A2")[1]
e_r = reference_a(eb)
@test e_r.data[1, :, :] ≈ eb.data[1, :, :] .- (eb.data[a1, :, :] .+ eb.data[a2, :, :]) ./ 2
e_r = reference_custom(eb; ref_list = ["Fp1-F3", "Cz"], ref_name = "test")
@test e_r.data[1, :, :] ≈ eb.data[1, :, :] .- eb.data[3, :, :]
@test labels(e_r)[1] == "Fp1-F3"
@test e_r.header.recording[:reference] == "test"
@test NeuroAnalyzer._parse_montage_entry("Fp1-F3", ["Fp1", "F3"]) == ("Fp1", "F3")
@test NeuroAnalyzer._parse_montage_entry("Fp1", ["Fp1", "F3"]) == ("Fp1", nothing)
@test NeuroAnalyzer._parse_montage_entry("A-B-C", ["A-B", "C"]) == ("A-B", "C")
@test_throws ArgumentError NeuroAnalyzer._parse_montage_entry("X-Y", ["Fp1"])
@test NeuroAnalyzer._eeg_ch_idx(eb) == ch_eeg
@test NeuroAnalyzer._inv_dist_weights([0.0 2.0; Inf 0.5]) == [0.0 0.5; 0.0 2.0]
@test NeuroAnalyzer._lap_suffix("plap", true, true) == "-wplapmed"
@test NeuroAnalyzer._lap_name("planar", 4, true, true) == "weighted planar Laplacian (median, 4)"
lo = NeuroAnalyzer._locs_ch_order(eb, [3, 1])
@test lo[!, :label] == ["F3", "Fp1"]
d = NeuroAnalyzer._locs_gc_dist(lo)
@test d[1, 1] ≈ 0 atol = 1e-7
@test d[1, 2] ≈ acos(clamp(sum(collect(lo[1, [:loc_x, :loc_y, :loc_z]]) .* collect(lo[2, [:loc_x, :loc_y, :loc_z]])) /
                         (NeuroAnalyzer.norm(collect(lo[1, [:loc_x, :loc_y, :loc_z]])) * NeuroAnalyzer.norm(collect(lo[2, [:loc_x, :loc_y, :loc_z]]))), -1, 1))
px, py = NeuroAnalyzer._locs_polar_xy(lo)
@test px ≈ lo.loc_radius .* cosd.(lo.loc_theta) && py ≈ lo.loc_radius .* sind.(lo.loc_theta)
ch_l, s_l = NeuroAnalyzer._laplacian_reference(eb; nn = 3, weighted = false, med = false, spherical = true)
@test ch_l == ch_eeg && size(s_l) == (length(ch_eeg), 2560, 10)
e_p = NeuroAnalyzer._reference_paired(eb, "A1", "A2", :l, "test")
@test e_p isa NeuroAnalyzer.NEURO
e_tmp = deepcopy(eb)
NeuroAnalyzer._append_label_suffix!(e_tmp, eb, [1, 2], "-x")
@test labels(e_tmp)[1:2] == ["Fp1-x", "Fp2-x"]
@test "Fp1-x" in e_tmp.locs[!, :label]
e_tmp = deepcopy(eb)
NeuroAnalyzer._copy_into!(e_tmp, reference_avg(eb))
@test e_tmp.data == reference_avg(eb).data

@info "Test: resample() values"
# 10 Hz sine: down- and upsampling keep the waveform (away from the edges)
x = sin.(2π * 10 .* (0:2559) ./ 256)
y = NeuroAnalyzer._resample(x; old_sr = 256, new_sr = 128)
@test length(y) == 1280
@test maximum(abs.(y[200:1080] .- sin.(2π * 10 .* (199:1079) ./ 128))) < 1e-3
y = NeuroAnalyzer._resample(x; old_sr = 256, new_sr = 512)
@test length(y) == 5120
@test maximum(abs.(y[400:4700] .- sin.(2π * 10 .* (399:4699) ./ 512))) < 1e-3
@test NeuroAnalyzer._resample(x; old_sr = 256, new_sr = 256) == x
@test_throws ArgumentError NeuroAnalyzer._resample(x; old_sr = 0, new_sr = 256)
r = NeuroAnalyzer._resample_obj(eb, 128, "resample")
@test sr(r) == 128 && size(r.data, 2) == 1280
@test r.history[end] == "resample(obj, new_sr=128)"

@info "Test: detect_powerline()"
x = 0.2 .* randn(NeuroAnalyzer.Random.MersenneTwister(3), 2560) .+ sin.(2π * 50 .* (0:2559) ./ 256)
@test NeuroAnalyzer._detect_powerline(x; fs = 256) == 50.0
e_tmp = keep_epoch(eb; ep = 1)
e_tmp.data[1, :, 1] = x
@test detect_powerline(keep_channel(e_tmp; ch = "Fp1")) == 50.0
# only 50 or 60 Hz: a stronger 10 Hz rhythm is not power-line noise; 0 if there is none
tt_pl = (0:2559) ./ 256
x_pl = 5 .* sin.(2π * 10 .* tt_pl) .+ 0.5 .* randn(NeuroAnalyzer.Random.MersenneTwister(4), 2560)
@test NeuroAnalyzer._detect_powerline(x_pl; fs = 256) == 0.0
@test NeuroAnalyzer._detect_powerline(x_pl .+ 0.5 .* sin.(2π * 60 .* tt_pl); fs = 256) == 60.0
@test NeuroAnalyzer._detect_powerline(x_pl .+ 0.5 .* sin.(2π * 50 .* tt_pl); fs = 256) == 50.0
@test NeuroAnalyzer._detect_powerline(x_pl .+ 0.5 .* sin.(2π * 50 .* tt_pl); fs = 256, pr = 100) == 0.0
@test NeuroAnalyzer._line_prominence(fill(1.0, 21), collect(40.0:60.0), 50) == 0
@test_throws ArgumentError NeuroAnalyzer._detect_powerline(x_pl; fs = 1)
e_tmp.data[1, :, 1] = x_pl
@test detect_powerline(keep_channel(e_tmp; ch = "Fp1")) == 0.0

@info "Test: mlinterpolate_channel()"
# constant regressor: the interpolated epoch = the mean of the channel in the reference epochs
cr = NeuroAnalyzer.MLJ.DeterministicConstantRegressor()
e_ml = mlinterpolate_channel(eb; ch = "Cz", ep = 1, ep_ref = 2:4, model = cr)
cz = get_channel(eb; ch = "Cz")[1]
@test all(e_ml.data[cz, :, 1] .≈ sum(eb.data[cz, :, 2:4]) / (3 * 2560))
@test e_ml.data[cz, :, 2] == eb.data[cz, :, 2]
e_tmp = deepcopy(eb)
mlinterpolate_channel!(e_tmp; ch = "Cz", ep = 1, ep_ref = 2:4, model = cr)
@test e_tmp.data == e_ml.data
@test_throws ArgumentError mlinterpolate_channel(eb; ch = "Cz", ep = 1, ep_ref = 1:3, model = cr)
@test_throws ArgumentError mlinterpolate_channel(eb; ch = "ECG", ep = 1, model = cr)

@info "Test: plinterpolate()"
ploc = eb.locs[NeuroAnalyzer._find_bylabel(eb.locs, labels(eb)[ch_eeg]), :]
pi_ = plinterpolate(eb.data[ch_eeg, 1, 1]; locs = ploc, ch = 1, ifactor = 50)
@test size(pi_.int_s) == (50, 50)
@test length(pi_.int_x) == 50 && all(isfinite, pi_.int_s)

@info "Test: read_cross_talk()"
ct = read_cross_talk(joinpath(testfiles_path, "meg", "SSS", "ct_sparse.fif"))
@test size(ct.decoupler, 1) == size(ct.decoupler, 2) == length(ct.ch_names)
@test all(NeuroAnalyzer.diag(ct.decoupler) .≈ 1)
@test_throws ArgumentError read_cross_talk(joinpath(testfiles_path, "no_such_file.fif"))

@info "Test: edit_montage()"
@test_throws ArgumentError edit_montage(joinpath(testfiles_path, "no_such_file.mnt"))
if Sys.islinux()
    withenv("EDITOR" => "true") do
        @test edit_montage(joinpath(NeuroAnalyzer.PATH, "montages", "bip_long.mnt")) === nothing
    end
end

@info "Test: small processing helpers"
@test NeuroAnalyzer._ssp_unit_scale("fT") == 1e15
@test NeuroAnalyzer._ssp_unit_scale("fT/cm") == 1e13
@test NeuroAnalyzer._ssp_unit_scale("T") == 1.0
s_i = reshape([1.0 2.0; 4.0 8.0], 2, 2, 1)
@test NeuroAnalyzer._intensity2od(s_i)[:, :, 1] ≈ -log.([1 / 1.5 2 / 1.5; 4 / 6 8 / 6])
z = NeuroAnalyzer._filter_g_analytic(2 .* cos.(2π * 10 .* (0:999) ./ 100); fs = 100, f = 10)
@test abs.(z[200:800]) ≈ fill(2.0, 601) atol = 1e-3
@test maximum(abs.(angle.(z[200:800] .* exp.(-im .* 2π * 10 .* (199:799) ./ 100)))) < 1e-3

@info "Test: Maxwell / cHPI helpers"
@test NeuroAnalyzer._sss_ncomp(8) == 80
@test NeuroAnalyzer._sss_ncomp(3) == 15
@test NeuroAnalyzer._sss_idx(1, -1) == 1
@test NeuroAnalyzer._sss_idx(2, 2) == 8
deg, ord = NeuroAnalyzer._sss_degrees_orders(2)
@test deg == [1, 1, 1, 2, 2, 2, 2, 2]
@test ord == [-1, 0, 1, -2, -1, 0, 1, 2]
@test NeuroAnalyzer._lfact(5) ≈ log(120)
@test NeuroAnalyzer._lfact(0) == 0.0
# dual numbers: value and gradient
d = NeuroAnalyzer._d3sqrt(NeuroAnalyzer._D3(4.0, 1.0, 0.0, 2.0))
@test (d.v, d.gx, d.gy, d.gz) == (2.0, 0.25, 0.0, 0.5)
@test NeuroAnalyzer._D3(3).v == 3.0 && NeuroAnalyzer._D3(3).gx == 0.0
# solid harmonics: gradients (dual numbers) = finite differences; regular harmonics of degree l scale as r^l
p0 = [0.03, -0.02, 0.07]
R, Ir = NeuroAnalyzer._sss_harmonics(p0, 3)
h = 1e-7
for k = 1:3
    e_k = zeros(3); e_k[k] = h
    Rp, Ip = NeuroAnalyzer._sss_harmonics(p0 .+ e_k, 3)
    Rm, Im = NeuroAnalyzer._sss_harmonics(p0 .- e_k, 3)
    local g = getfield.(R, (:gx, :gy, :gz)[k])
    @test g ≈ (getfield.(Rp, :v) .- getfield.(Rm, :v)) ./ (2h) rtol = 1e-5
    gi = getfield.(Ir, (:gx, :gy, :gz)[k])
    @test gi ≈ (getfield.(Ip, :v) .- getfield.(Im, :v)) ./ (2h) rtol = 1e-5
end
R2, _ = NeuroAnalyzer._sss_harmonics(2 .* p0, 3)
deg, _ = NeuroAnalyzer._sss_degrees_orders(3)
@test getfield.(R2, :v) ≈ getfield.(R, :v) .* 2.0 .^ deg
# rotations ↔ quaternions
q = [0.1, -0.2, 0.3]
Rq = NeuroAnalyzer._quat2rot(q)
@test Rq' * Rq ≈ NeuroAnalyzer.I(3) && NeuroAnalyzer.det(Rq) ≈ 1
@test NeuroAnalyzer._rot2quat(Rq) ≈ q
@test NeuroAnalyzer._hom(Rq, [1.0, 2.0, 3.0]) == [Rq [1.0, 2.0, 3.0]; 0 0 0 1]
# Kabsch: exact recovery of a rigid transform
A = randn(NeuroAnalyzer.Random.MersenneTwister(5), 6, 3)
B = (Rq * A')' .+ [0.1 -0.2 0.05]
Rk, tk = NeuroAnalyzer._kabsch(A, B)
@test Rk ≈ Rq && tk ≈ [0.1, -0.2, 0.05]
# Nelder-Mead: minimum of a quadratic
xm_, fm_ = NeuroAnalyzer._nelder_mead(x -> sum(abs2, x .- [1.0, -2.0]), [0.0, 0.0])
@test xm_ ≈ [1.0, -2.0] atol = 1e-4
@test fm_ < 1e-8
# orthonormal basis / pseudoinverse
Ao = [1.0 2.0; 2.0 4.0; 0.0 1.0]
Q = NeuroAnalyzer._orth(Ao)
@test Q' * Q ≈ NeuroAnalyzer.I(size(Q, 2))
@test Q * Q' * Ao ≈ Ao
Pi, nrm = NeuroAnalyzer._col_norm_pinv(Ao)
@test Pi ≈ NeuroAnalyzer.pinv(Ao)
# tSSS projector: a time course shared by both inputs is found
tk_ = sin.(2π .* (1:200) ./ 37)
U = NeuroAnalyzer._tsss_projector(0.01 .* randn(NeuroAnalyzer.Random.MersenneTwister(1), 5, 200) .+ [1.0, 2, 3, 4, 5] * tk_',
                                   0.01 .* randn(NeuroAnalyzer.Random.MersenneTwister(2), 5, 200) .+ [5.0, 4, 3, 2, 1] * tk_', 0.98)
@test size(U, 2) >= 1
@test abs(NeuroAnalyzer.norm(U' * tk_) / NeuroAnalyzer.norm(tk_)) > 0.9
# coil integration points: weights of a magnetometer sum to its area normal
meg_raw = import_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fif"))
mr = meg_raw.header.recording
pts, wn, cidx = NeuroAnalyzer._sss_coil_points(mr[:coil_loc][mr[:magnetometers][1:2], :], mr[:coil_id][mr[:magnetometers][1:2]], NeuroAnalyzer.Matrix{Float64}(NeuroAnalyzer.I(4)))
@test size(pts, 2) == 3 && size(wn) == size(pts) && sort(unique(cidx)) == [1, 2]
@test size(NeuroAnalyzer._vv_mag_acc(0.0258), 2) == 4
@test size(NeuroAnalyzer._vv_mag_normal(0.0258), 2) == 4
cp_ = NeuroAnalyzer._chpi_coil_points(mr[:coil_loc][mr[:magnetometers][1:30], :], mr[:coil_id][mr[:magnetometers][1:30]])
@test size(cp_[1], 2) == 3
# dipole forward model and fit quality: the exact field of a dipole has goodness of fit 1
fwd(rq) = NeuroAnalyzer._chpi_mag_dipole_fwd(rq, cp_..., 30)
rq = [0.0, 0.0, 0.04]
b = fwd(rq) * [1e-8, 0.0, 2e-8]
@test NeuroAnalyzer._chpi_gof(rq, b, fwd) ≈ 1.0
@test NeuroAnalyzer._chpi_gof([0.0, 0.03, 0.0], b, fwd) < 1.0

@info "Test: SSS basis and helpers"
# a magnetic dipole at the expansion origin is a pure l = 1 internal multipole: its field on the
# coils lies in the span of the first 3 internal basis vectors
mag_i = mr[:magnetometers][1:20]
grad_i = mr[:gradiometers][1:20]
sel = vcat(mag_i, grad_i)
pts, wn, cidx = NeuroAnalyzer._sss_coil_points(mr[:coil_loc][sel, :], mr[:coil_id][sel], NeuroAnalyzer.Matrix{Float64}(NeuroAnalyzer.I(4)))
S_in, S_out = NeuroAnalyzer._sss_basis(pts, wn, cidx, length(sel), [0.0, 0.0, 0.0], 3, 2)
@test size(S_in) == (40, NeuroAnalyzer._sss_ncomp(3))
@test size(S_out) == (40, NeuroAnalyzer._sss_ncomp(2))
b = NeuroAnalyzer._chpi_mag_dipole_fwd([0.0, 0.0, 0.0], pts, wn, cidx, length(sel)) * [1e-8, -2e-8, 0.5e-8]
B1 = S_in[:, 1:3]
@test NeuroAnalyzer.norm(b .- B1 * (B1 \ b)) / NeuroAnalyzer.norm(b) < 1e-6
# an off-origin dipole needs higher orders
b2 = NeuroAnalyzer._chpi_mag_dipole_fwd([0.0, 0.02, 0.04], pts, wn, cidx, length(sel)) * [1e-8, 0.0, 0.0]
@test NeuroAnalyzer.norm(b2 .- B1 * (B1 \ b2)) / NeuroAnalyzer.norm(b2) > 1e-3
# information-based regularization returns a subset of the internal components
keep = NeuroAnalyzer._sss_regularize_in(hcat(S_in, S_out), size(S_in, 2), 3, collect(1:size(S_out, 2)))
@test !isempty(keep) && all(in(1:size(S_in, 2)), keep)
# fine calibration orientation: unit columns; an orthogonal calibration frame is used as is
o0 = vec([1.0 0 0; 0 1 0; 0 0 1])
@test NeuroAnalyzer._fine_cal_orientation(o0, 2 .* o0) ≈ o0
oc = vec(NeuroAnalyzer._quat2rot([0.05, 0.0, 0.0]))
@test NeuroAnalyzer._fine_cal_orientation(o0, oc) ≈ oc
# calibration file that only rotates ez: the original frame is rotated onto the new ez
ez = [0.0, 0.1, 1.0] ./ NeuroAnalyzer.norm([0.0, 0.1, 1.0])
oz = vcat([1.0, 0, 0], [0.0, 1, 0], ez)
R_ = reshape(NeuroAnalyzer._fine_cal_orientation(o0, oz), 3, 3)
@test R_[:, 3] ≈ ez && R_' * R_ ≈ NeuroAnalyzer.I(3)
# average head transform of a buffer: one position → that position; two equal-length positions with
# opposite translations → zero translation
T1 = NeuroAnalyzer._hom(NeuroAnalyzer._quat2rot([0.1, 0.0, 0.0]), [0.01, 0.0, 0.0])
@test NeuroAnalyzer._sss_avg_trans([0], [T1], [[0.1, 0.0, 0.0, 0.01, 0.0, 0.0]], 10, 20) == T1
Ta = NeuroAnalyzer._sss_avg_trans([0, 15], [T1, T1], [[0.1, 0.0, 0.0, 0.01, 0.0, 0.0], [0.1, 0.0, 0.0, -0.01, 0.0, 0.0]], 10, 20)
@test Ta[1:3, 4] ≈ [0.0, 0.0, 0.0] atol = 1e-12
@test Ta[1:3, 1:3] ≈ NeuroAnalyzer._quat2rot([0.1, 0.0, 0.0])

@info "Test: in-place variants = non-mutating variants"
cmpdata(a, b) = a.data ≈ b.data && labels(a) == labels(b)
erp = average_epochs(eb)
for (f, f!, o, kw) in (
    (add_signal, add_signal!, eb, (ch = "Fp1", s = collect(1.0:2560))),
    (cbp, cbp!, eb, (ch = "Fp1", frq = 4)),
    (csd, csd!, eb, ()),
    (denoise_fft, denoise_fft!, eb, (ch = "Fp1",)),
    (denoise_wien, denoise_wien!, eb, (ch = "Fp1",)),
    (denoise_dwd, denoise_dwd!, eb, (ch = "Fp1", wt = wavelet(WT.haar))),
    (NeuroAnalyzer.derivative, derivative!, eb, (ch = "Fp1",)),
    (detrend, detrend!, eb, (ch = "Fp1", type = :linear)),
    (filter_g, filter_g!, eb, (ch = "Fp1", f = 20)),
    (filter_mavg, filter_mavg!, eb, (ch = "Fp1", k = 2)),
    (filter_mmed, filter_mmed!, eb, (ch = "Fp1", k = 2)),
    (filter_poly, filter_poly!, eb, (ch = "Fp1",)),
    (filter_sg, filter_sg!, eb, (ch = "Fp1",)),
    (NeuroAnalyzer.normalize, NeuroAnalyzer.normalize!, eb, (ch = "Fp1", method = :zscore)),
    (normpower, normpower!, eb, (ch = "Fp1",)),
    (remove_dc, remove_dc!, eb, (ch = "Fp1",)),
    (taper, taper!, eb, (ch = "Fp1", t = NeuroAnalyzer.DSP.hanning(2560))),
    (wbp, wbp!, eb, (ch = "Fp1", frq = 4)),
    (reference_avg, reference_avg!, eb, ()),
    (reference_a, reference_a!, eb, ()),
    (reference_ce, reference_ce!, eb, (ch = "Fp1",)),
    (reference_plap, reference_plap!, eb, ()),
    (reference_slap, reference_slap!, eb, ()),
    (reference_custom, reference_custom!, eb, ()),
    (resample, resample!, eb, (new_sr = 128,)),
    (downsample, downsample!, eb, (new_sr = 128,)),
    (upsample, upsample!, eb, (new_sr = 512,)),
    (lrinterpolate_channel, lrinterpolate_channel!, eb, (ch = "Cz", ep = 1)),
    (plinterpolate_channel, plinterpolate_channel!, eb, (ch = "Cz", ep = 1)),
    (artrem_cwd, artrem_cwd!, eb, (ch = "Fp1", ep = 1, tseg = (0.2, 0.4), fseg = (10, 12.5))),
    (denoise_cwd, denoise_cwd!, eb, (ch = "Fp1", nf = 50)),
    (sort_epochs, sort_epochs!, erp, (s = collect(nepochs(erp):-1:2),)),
    (intensity2od, intensity2od!, nirs, ()),
)
    o2 = deepcopy(o)
    f!(o2; kw...)
    @test cmpdata(o2, f(o; kw...))
end
n_od = intensity2od(nirs)
n_tmp = deepcopy(n_od)
od2conc!(n_tmp)
@test n_tmp.data ≈ od2conc(n_od).data
pc, _, _, pc_model = pca_decompose(eb; ch = "all", n = 4)
e_tmp = deepcopy(eb)
pca_reconstruct!(e_tmp; ch = "all", pc = pc, pc_model = pc_model)
@test e_tmp.data ≈ pca_reconstruct(eb; ch = "all", pc = pc, pc_model = pc_model).data
ica_d = ica_decompose(eeg; ch = "all", n = 5, iter = 10)
e_tmp = deepcopy(eeg)
ica_remove!(e_tmp; ch = "all", ic_idx = 1:2, ic = ica_d.ic, ic_mw = ica_d.ic_mw)
@test e_tmp.data ≈ ica_remove(eeg; ch = "all", ic_idx = 1:2, ic = ica_d.ic, ic_mw = ica_d.ic_mw).data
e_tmp = keep_epoch(eb; ep = 1)
rp_out = remove_pops(o_p; ch = "all")
o_tmp = deepcopy(o_p)
@test remove_pops!(o_tmp; ch = "all") == rp_out[2]
@test o_tmp.data == rp_out[1].data
m_tmp = deepcopy(m_ssp)
apply_ssp_projectors!(m_tmp)
@test m_tmp.data ≈ apply_ssp_projectors(m_ssp).data
r_tmp = deepcopy(raw)
maxwell_filter!(r_tmp; origin = o, regularize = :none)
@test r_tmp.data ≈ maxwell_filter(raw; origin = o, regularize = :none).data
e_tmp = keep_epoch(eb; ep = 1)
detect_powerline!(e_tmp)
@test e_tmp.header.recording[:line_frequency] == detect_powerline(keep_epoch(eb; ep = 1)) ||
      e_tmp.header.recording[:line_frequency] isa Real

@info "Test: process values (synthetic)"
# derivative of a linear ramp; detrend removes a linear trend; remove_dc removes the mean
@test NeuroAnalyzer._derivative(2.0 .* collect(1:10))[2:9] ≈ fill(2.0, 8)
x = 3.0 .* (1:100) .+ sin.(2π .* (1:100) ./ 20)
@test NeuroAnalyzer._detrend(x; type = :linear) ≈ NeuroAnalyzer._detrend(sin.(2π .* (1:100) ./ 20); type = :linear) atol = 1e-10
@test sum(NeuroAnalyzer._remove_dc(x)) ≈ 0 atol = 1e-9
# taper: elementwise product
@test NeuroAnalyzer._taper([1.0, 2.0, 3.0]; t = [0.0, 0.5, 1.0]) == [0.0, 1.0, 3.0]
# moving average of a constant, Savitzky-Golay of a polynomial of its order
@test NeuroAnalyzer._filter_sg(collect(1.0:50) .^ 2; order = 2, window = 11)[10:40] ≈ (collect(10.0:40)) .^ 2 rtol = 1e-8
# normalize: z-score of a NEURO channel
# bych = false: the selected data are normalized jointly (all epochs); bych = true: per channel and epoch
e_n = NeuroAnalyzer.normalize(eb; ch = "Fp1", method = :zscore)
@test sum(e_n.data[1, :, :]) ≈ 0 atol = 1e-6
@test NeuroAnalyzer.Statistics.std(vec(e_n.data[1, :, :])) ≈ 1
e_n = NeuroAnalyzer.normalize(eb; ch = "Fp1", method = :zscore, bych = true)
@test sum(e_n.data[1, :, 3]) ≈ 0 atol = 1e-8
@test NeuroAnalyzer.Statistics.std(e_n.data[1, :, 3]) ≈ 1
# add_signal adds the signal to every epoch
e_a = add_signal(eb; ch = "Fp1", s = collect(1.0:2560))
@test e_a.data[1, :, 4] == eb.data[1, :, 4] .+ (1.0:2560)
@test e_a.data[2, :, 4] == eb.data[2, :, 4]


@info "Test: remove_powerline() values"
# 10 Hz signal + 50 Hz line noise (+ a little noise): the 50 Hz component is removed, 10 Hz kept
e_pl = keep_epoch(epoch(eeg; ep_len = 20); ep = 1)
tp_ = (0:(epoch_len(e_pl) - 1)) ./ sr(e_pl)
e_pl.data[1, :, 1] = sin.(2π * 10 .* tp_) .+ sin.(2π * 50 .* tp_) .+ 0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(1), length(tp_))
e_out, df_pl = remove_powerline(e_pl; ch = "Fp1", pl_frq = 50)
amp_f(s, f) = (mid = 1001:(length(s) - 1000); A = [sin.(2π * f .* tp_[mid]) cos.(2π * f .* tp_[mid])]; hypot((A \ s[mid])...))
@test amp_f(e_out.data[1, :, 1], 50) < 0.05
@test amp_f(e_out.data[1, :, 1], 10) > 0.95
@test df_pl[1, :channel] == "Fp1" && df_pl[1, :frequency] == 50 && df_pl[1, :bandwidth] > 0
# only channels with a peak are filtered; no harmonic in this signal
@test nrow(df_pl) == 1
# other channels are not modified
@test e_out.data[2, :, 1] == e_pl.data[2, :, 1]
# harmonics: 50 Hz + 100 Hz removed; a channel without line noise is not modified
e_pl.data[1, :, 1] = sin.(2π * 10 .* tp_) .+ sin.(2π * 50 .* tp_) .+ 0.3 .* sin.(2π * 100 .* tp_) .+ 0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(1), length(tp_))
e_pl.data[2, :, 1] = sin.(2π * 10 .* tp_) .+ 0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(2), length(tp_))
e_out, df_pl = remove_powerline(e_pl; ch = ["Fp1", "Fp2"], pl_frq = 50)
@test df_pl[!, :channel] == ["Fp1", "Fp1"] && df_pl[!, :frequency] == [50, 100]
@test amp_f(e_out.data[1, :, 1], 50) < 0.05 && amp_f(e_out.data[1, :, 1], 100) < 0.03
@test amp_f(e_out.data[1, :, 1], 10) > 0.95
@test e_out.data[2, :, 1] == e_pl.data[2, :, 1]
# 60 Hz
e_pl.data[1, :, 1] = sin.(2π * 10 .* tp_) .+ sin.(2π * 60 .* tp_) .+ 0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(1), length(tp_))
e_out, df_pl = remove_powerline(e_pl; ch = "Fp1", pl_frq = 60)
@test df_pl[!, :frequency] == [60] && amp_f(e_out.data[1, :, 1], 60) < 0.05
@test_throws ArgumentError remove_powerline(e_pl; ch = "Fp1", pl_frq = 0)
@test_throws ArgumentError remove_powerline(e_pl; ch = "Fp1", pl_frq = 50, d = 0)
@test_throws ArgumentError remove_powerline(epoch(eeg; ep_len = 10); ch = "Fp1")

true
