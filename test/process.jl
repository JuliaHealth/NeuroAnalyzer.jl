using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets
using LinearAlgebra: LinearAlgebra

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
nirs = import_nirs(joinpath(testfiles_path, "nirs", "fnirs-test-nirs.nirs"))
meg = import_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fif"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)
epoch!(meg; ep_len = 2)
keep_epoch!(meg; ep = 1:10)
load_locs!(
    e10;
    file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"),
)
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
r = NeuroAnalyzer._artrem_cwd(
    x,
    collect(0:1023) ./ 256;
    fs = 256,
    tseg = (1, 2),
    fseg = (60, 70),
)
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
s = NeuroAnalyzer._denoise_wien(a1; fs = 10)
@test size(s) == (2, 3, 2)
xw = randn(NeuroAnalyzer.Random.MersenneTwister(1), 3, 256, 2)
@test NeuroAnalyzer._denoise_wien(xw; fs = 256) == NeuroAnalyzer._denoise_wien(xw; fs = 256)
# a signal common to all channels is kept, independent noise is removed, and activity
# present in one channel only is attenuated
rng_w = NeuroAnalyzer.Random.MersenneTwister(2)
t_w = (0:2559) ./ 256
common_w = 10 .* sin.(2π * 10 .* t_w)
xw = zeros(8, 2560, 1)
for idx = 1:8
    xw[idx, :, 1] = common_w .+ 10 .* randn(rng_w, 2560)
end
xw[1, :, 1] .+= 10 .* sin.(2π * 6 .* t_w)
yw = NeuroAnalyzer._denoise_wien(xw; fs = 256)
@test sqrt(sum(abs2, yw[2, :, 1] .- common_w) / 2560) < 3
fw_in = abs.(NeuroAnalyzer.rfft(xw[1, :, 1]))
fw_out = abs.(NeuroAnalyzer.rfft(yw[1, :, 1]))
@test fw_out[61] / fw_in[61] < 0.1      # 6 Hz, channel 1 only
@test fw_out[101] / fw_in[101] > 0.9    # 10 Hz, all channels
# ref = :epochs: noisy trials of the same response
xe = zeros(1, 256, 30)
for idx = 1:30
    xe[1, :, idx] = common_w[1:256] .+ 10 .* randn(rng_w, 256)
end
ye = NeuroAnalyzer._denoise_wien(xe; fs = 256, ref = :epochs)
@test sqrt(sum(abs2, ye[1, :, 1] .- common_w[1:256]) / 256) <
      sqrt(sum(abs2, xe[1, :, 1] .- common_w[1:256]) / 256) / 2
@test_throws ArgumentError NeuroAnalyzer._denoise_wien(xw[1:1, :, :]; fs = 256)
@test_throws ArgumentError NeuroAnalyzer._denoise_wien(xw; fs = 256, ref = :epochs)
@test_throws ArgumentError NeuroAnalyzer._denoise_wien(xw; fs = 256, ref = :x)
@test_throws ArgumentError NeuroAnalyzer._denoise_wien(xw; fs = 256, fsmooth = -1)
e10_tmp = denoise_wien(e10; ch = "all")
@test size(e10_tmp) == size(e10)
@test size(denoise_wien(e10; ch = "all", ref = :epochs)) == size(e10)

@info "Test: derivative()"
@test NeuroAnalyzer._derivative(v1) == [0.5, 1.0, 1.0, 1.0, 0.5]
@test NeuroAnalyzer._derivative(a1) ==
      [0.0 0.0 0.0; 0.0 0.0 0.0;;; 0.0 0.0 0.0; 0.0 0.0 0.0]
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
s, bn, bf =
    bpsplit(e10; ch = "Fp1", band = [:alpha, :gamma], fprototype = :butterworth, order = 4)
@test bn == [:alpha, :gamma]
@test size(s) == (2, 1, 2560, 10)
# 10 Hz + 20 Hz sines: each band keeps its own component (std of a unit sine = 1/√2)
e_tmp = deepcopy(e10)
e_tmp.data[1, :, :] .=
    sin.(2π * 10 .* e_tmp.epoch_time) .+ sin.(2π * 20 .* e_tmp.epoch_time)
std_mid(x) = (m = x[600:1960]; sqrt(sum((m .- sum(m) / length(m)) .^ 2) / (length(m) - 1)))
for (p, kw) in ((:fir, (;)), (:firls, (;)), (:remez, (;)), (:butterworth, (order = 6,)))
    local s =
        bpsplit(e_tmp; ch = "Fp1", band = [:alpha, :beta, :gamma], fprototype = p, kw...).s
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
@test NeuroAnalyzer._filter_mavg(vcat(v1, v1), k = 2) ==
      [2.0, 2.5, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.5, 4.0]
# weighted mean: a constant signal is unchanged
@test NeuroAnalyzer._filter_mavg(fill(2.0, 20), k = 3, ww = NeuroAnalyzer.DSP.hanning(7)) ≈
      fill(2.0, 20)
e10_tmp = filter_mavg(e10; ch = "all", k = 2)
@test size(e10_tmp) == size(e10)

@info "Test: filter_mmed()"
@test NeuroAnalyzer._filter_mmed(vcat(v1, v1), k = 2) ==
      [2.0, 2.5, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.5, 4.0]
# weighted median: values are not rescaled by the weights
@test NeuroAnalyzer._filter_mmed(fill(2.0, 20), k = 3, ww = NeuroAnalyzer.DSP.hanning(7)) ≈
      fill(2.0, 20)
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
# low-pass and high-pass frequencies of the channels follow the filters (EDF: HP 0.18 Hz, LP 104 Hz)
rec0 = eeg.header.recording
@test rec0[:lowpass][1] == 104.0 && rec0[:highpass][1] == 0.18 && !haskey(rec0, :prefiltering)
f_tmp = NeuroAnalyzer.filter(eeg; ch = "Fp1", fprototype = :butterworth, ftype = :lp,
                             cutoff = 40, order = 4)
i1, i2 = get_channel(f_tmp; ch = ["Fp1", "Fp2"])
@test f_tmp.header.recording[:lowpass][i1] == 40.0
@test f_tmp.header.recording[:lowpass][i2] == 104.0 && rec0[:lowpass][i1] == 104.0
f_tmp = NeuroAnalyzer.filter(f_tmp; ch = "Fp1", fprototype = :butterworth, ftype = :bp,
                             cutoff = (1, 80), order = 4)
@test f_tmp.header.recording[:lowpass][i1] == 40.0 && f_tmp.header.recording[:highpass][i1] == 1.0
@test all(isempty, rec0[:notch])
f_tmp2 = NeuroAnalyzer.filter(f_tmp; ch = "Fp1", fprototype = :iirnotch, cutoff = 50, bw = 2)
@test f_tmp2.header.recording[:lowpass] == f_tmp.header.recording[:lowpass]
@test f_tmp2.header.recording[:highpass] == f_tmp.header.recording[:highpass]
# notch filters are recorded in :notch (once per frequency)
@test f_tmp2.header.recording[:notch][i1] == [50.0] && isempty(f_tmp2.header.recording[:notch][i2])
NeuroAnalyzer.filter!(f_tmp2; ch = "Fp1", fprototype = :iirnotch, cutoff = 50, bw = 2)
@test f_tmp2.header.recording[:notch][i1] == [50.0]
# resampling: notches above the new Nyquist frequency are dropped
@test isempty(resample(f_tmp2; new_sr = 64).header.recording[:notch][i1])
f_tmp2 = NeuroAnalyzer.filter(f_tmp; ch = "Fp1", fprototype = :butterworth, ftype = :hp,
                             cutoff = 0.5, order = 4)
@test f_tmp2.header.recording[:highpass][i1] == 1.0
# resampling: no frequencies above the new Nyquist frequency
@test all(==(64.0), resample(eeg; new_sr = 128).header.recording[:lowpass])
@test all(==(104.0), resample(eeg; new_sr = 512).header.recording[:lowpass])
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
@test NeuroAnalyzer._filter_g(2 .* cos.(2π * 10 .* (0:999) ./ 100); fs = 100, f = 10)[200:800] ≈
      fill(2.0, 601) atol = 1e-3
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
@test e10_int.data[1, setdiff(1:epoch_len(e10), r_seg), 1] ==
      e10.data[1, setdiff(1:epoch_len(e10), r_seg), 1]
@test_throws ArgumentError lrinterpolate_channel(
    e10_tmp;
    ch = "Fp1",
    ep = 1,
    ep_ref = Int64[],
)
@test_throws ArgumentError lrinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1, seg = (4, 2))
@test_throws ArgumentError lrinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1, ep_ref = 1)
# bad channels are not used as predictors
e10_bad = deepcopy(e10_tmp)
e10_bad.data[2, :, :] .= 1e6 .* randn(size(e10_bad.data, 2), size(e10_bad.data, 3))
e10_bad.header.recording[:bad_channel][2] = true
@test lrinterpolate_channel(e10_bad; ch = "Fp1", ep = 1).data[1, :, 1] ≈
      lrinterpolate_channel(delete_channel(e10_tmp; ch = "Fp2"); ch = "Fp1", ep = 1).data[
    1,
    :,
    1,
]

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
      [                         -1.0 -0.6 -0.19999999999999996;
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
      [                         1.0 0.7071067811865475 0.5773502691896258;
    0.5 0.4472135954999579 0.4082482904638631]
@test normalize_invroot(m1; bych = false) ==
      [                         1.0 0.7071067811865475 0.5773502691896258;
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
# :sph delegates to sphinterpolate_channel()
@test plinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1, seg = (2, 4), imethod = :sph).data ==
      sphinterpolate_channel(e10_tmp; ch = "Fp1", ep = 1, seg = (2, 4)).data

@info "Test: sphinterpolate_channel()"
# G matrix: the same as in gh() (CSD) for the same parameters
ch_eeg_sph = get_channel(e10; ch = "eeg")
locs_sph = e10.locs[[findfirst(==(l), e10.locs.label) for l in labels(e10)[ch_eeg_sph]], :]
p_sph = Matrix{Float64}(locs_sph[:, [:loc_x, :loc_y, :loc_z]])
p_sph ./= sqrt.(sum(abs2, p_sph; dims = 2))
@test NeuroAnalyzer._sph_g(p_sph, p_sph; m = 4, n = 8) ≈ gh(locs_sph; m = 4, n = 8).G
# a constant field is reproduced exactly (the spline has a constant term)
e10_c = deepcopy(e10)
e10_c.data[ch_eeg_sph, :, :] .= 7.0
@test sphinterpolate_channel(e10_c; ch = "C3").data[get_channel(e10; ch = "C3")[1], :, :] ≈
      fill(7.0, epoch_len(e10), nepochs(e10)) atol = 1e-6
# only the interpolated channels (and epochs) change; several channels at once
e10_tmp = deepcopy(e10)
c34 = get_channel(e10; ch = ["C3", "C4"])
e10_tmp.data[c34, :, :] .= 0
e10_int = sphinterpolate_channel(e10_tmp; ch = ["C3", "C4"], ep = 1:2)
@test all(!iszero, e10_int.data[c34, :, 1:2])
@test all(iszero, e10_int.data[c34, :, 3:end])
@test e10_int.data[setdiff(1:nchannels(e10), c34), :, :] ==
      e10_tmp.data[setdiff(1:nchannels(e10), c34), :, :]
# segment only
e10_int = sphinterpolate_channel(e10_tmp; ch = "C3", ep = 1, seg = (2, 4))
r_seg = NeuroAnalyzer._epoch_seg(e10_tmp, (2, 4))
@test all(!iszero, e10_int.data[c34[1], r_seg, 1])
@test all(iszero, e10_int.data[c34[1], setdiff(1:epoch_len(e10), r_seg), 1])
# bad channels are not used
e10_bad = deepcopy(e10)
fp2 = get_channel(e10; ch = "Fp2")[1]
e10_bad.data[fp2, :, :] .= 1e6
e10_bad.header.recording[:bad_channel][fp2] = true
@test sphinterpolate_channel(e10_bad; ch = "C3").data[c34[1], :, :] ≈
      sphinterpolate_channel(delete_channel(e10; ch = "Fp2"); ch = "C3").data[
    get_channel(delete_channel(e10; ch = "Fp2"); ch = "C3")[1],
    :,
    :,
]
# in-place
e10_tmp2 = deepcopy(e10_tmp)
sphinterpolate_channel!(e10_tmp2; ch = "C3")
@test e10_tmp2.data == sphinterpolate_channel(e10_tmp; ch = "C3").data
@test_throws ArgumentError sphinterpolate_channel(e10; ch = "ECG")
@test_throws ArgumentError sphinterpolate_channel(e10; ch = "C3", m = 1)
@test_throws ArgumentError sphinterpolate_channel(e10; ch = "C3", seg = (0, 100))

@info "Test: interpolation, reset_bads"
e10_rb = deepcopy(e10)
c3 = get_channel(e10_rb; ch = "C3")[1]
e10_rb.header.recording[:bad_channel][c3] = true
for f in (sphinterpolate_channel, plinterpolate_channel)
    @test f(e10_rb; ch = "C3", ep = 1).header.recording[:bad_channel][c3]
    @test !f(e10_rb; ch = "C3", ep = 1, reset_bads = true).header.recording[:bad_channel][c3]
end
@test lrinterpolate_channel(e10_rb; ch = "C3", ep = 1).header.recording[:bad_channel][c3]
@test !lrinterpolate_channel(e10_rb; ch = "C3", ep = 1, reset_bads = true).header.recording[:bad_channel][c3]
for f! in (sphinterpolate_channel!, plinterpolate_channel!)
    e_rb = deepcopy(e10_rb)
    f!(e_rb; ch = "C3", ep = 1, reset_bads = true)
    @test !e_rb.header.recording[:bad_channel][c3]
end
e_tmp = deepcopy(e10_rb)
lrinterpolate_channel!(e_tmp; ch = "C3", ep = 1, reset_bads = true)
@test !e_tmp.header.recording[:bad_channel][c3]

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
@test all(
    isapprox.(reference_avg(e10_sim; weighted = true).data[ch_eeg, :, :], 0; atol = 1e-10),
)
# weights follow the polar (topographic) layout: out_i = x_i - Σ_j w_ij x_j / Σ_j w_ij, w_ij = 1/d_ij
e10_sim.data[ch_eeg, :, :] .= 0.0
e10_sim.data[ch_eeg[1], :, :] .= 1.0
out = reference_avg(e10_sim; weighted = true).data[ch_eeg, 1, 1]
locs = e10.locs[NeuroAnalyzer._find_bylabel(e10.locs, labels(e10)[ch_eeg]), :]
px = locs[!, :loc_radius] .* cosd.(locs[!, :loc_theta])
py = locs[!, :loc_radius] .* sind.(locs[!, :loc_theta])
n = length(ch_eeg)
w = [i == j ? 0.0 : 1 / hypot(px[i] - px[j], py[i] - py[j]) for i = 1:n, j = 1:n]
w[.!isfinite.(w)] .= 0.0
expected = [(i == 1 ? 1.0 : 0.0) - w[i, 1] / sum(w[i, :]) for i = 1:n]
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

@info "Test: reference_rest()"
# shell factors: equal conductivities → homogeneous sphere
@test NeuroAnalyzer._shell_factors([0.87, 0.92, 1.0], [0.33, 0.33, 0.33], 50) ≈ ones(50)
@test NeuroAnalyzer._shell_factors([1.0], [0.33], 10) ≈ ones(10)
e10_tmp = reference_rest(e10; ndip = 500)
@test size(e10_tmp) == size(e10)
@test all(isfinite, e10_tmp.data)
@test startswith(e10_tmp.header.recording[:reference], "REST")
# REST does not depend on the original common reference
@test reference_rest(reference_ce(e10; ch = "Cz"); ndip = 500).data[ch_eeg, :, :] ≈
      e10_tmp.data[ch_eeg, :, :]
@test reference_rest(reference_avg(e10); ndip = 500).data[ch_eeg, :, :] ≈
      e10_tmp.data[ch_eeg, :, :]
# non-EEG channels are not modified
ch_other = setdiff(1:nchannels(e10), ch_eeg)
@test e10_tmp.data[ch_other, :, :] == e10.data[ch_other, :, :]
e10_tmp = reference_rest(e10; ndip = 500, model = :sphere, rdip = 0.9)
@test all(isfinite, e10_tmp.data)
@test_throws ArgumentError reference_rest(e10; rdip = 0.9)
@test_throws ArgumentError reference_rest(e10; model = :foo)
@test_throws ArgumentError reference_rest(e10; ndip = 0)

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
@test round.(NeuroAnalyzer._tconv(v1, kernel = [0.2, 0.1, 0.2]), digits = 2) ==
      [0.5, 1.0, 1.5, 2.0, 1.3]
@test round.(NeuroAnalyzer._tconv(a1, kernel = [0.2, 0.1, 0.2]), digits = 2) ==
      [0.3 0.5 0.3; 0.3 0.5 0.3;;; 0.3 0.5 0.3; 0.3 0.5 0.3]
e10_tmp = tconv(e10; ch = "all", kernel = [0.2, 0.1, 0.2])
@test size(e10_tmp) == size(e10)
# "same" convolution for kernels of any length; a symmetric kernel does not shift the signal
x_imp = zeros(50)
x_imp[25] = 1
for k = 1:6
    @test length(NeuroAnalyzer._tconv(x_imp; kernel = ones(k))) == 50
    @test length(NeuroAnalyzer._fconv(x_imp; kernel = ones(k))) == 50
end
@test argmax(NeuroAnalyzer._tconv(x_imp; kernel = [1.0, 2, 3, 2, 1])) == 25
@test size(tconv(e10; ch = "all", kernel = ones(4))) == size(e10)
@test size(fconv(e10; ch = "all", kernel = ones(4))) == size(e10.data)

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
# intensities: zero → error, negative → absolute value
n_bad = deepcopy(nirs)
n_bad.data[1, 10, 1] = 0
@test_throws ArgumentError intensity2od(n_bad)
n_neg = deepcopy(nirs)
n_neg.data[1, :, :] .*= -1
@test intensity2od(n_neg).data[5:end, :, :] ≈ n_tmp.data[5:end, :, :]
# the optode position unit is read from the file
@test nirs.header.recording[:length_unit] == "mm"
# concentration channels: after the optical channels, before the auxiliary channels, with
# their optode pairs (wavelength index 0)
r2 = n_tmp2.header.recording
n_opt = 8
@test r2[:channel_type][(n_opt + 1):(n_opt + 6)] == repeat(["nirs_hbo", "nirs_hbr", "nirs_hbt"], 2)
@test all(==("nirs_aux"), r2[:channel_type][(n_opt + 7):end])
@test length(r2[:wavelength_index]) == n_opt + 6
@test all(iszero, r2[:wavelength_index][(n_opt + 1):end])
@test r2[:optode_pairs][(n_opt + 1):(n_opt + 3), :] == repeat(r2[:optode_pairs][5:5, :], 3)
@test all(==("μM·mm"), r2[:unit][(n_opt + 1):(n_opt + 6)])
# the OD → concentration conversion can be repeated after the concentration channels
@test od2conc(n_tmp2) isa NeuroAnalyzer.NEURO
# DPF (Scholkmann & Wolf 2013)
@test nirs_dpf(; age = 30, wl = 760) ≈ 6.30 atol = 0.01
@test nirs_dpf(; age = 30, wl = [760, 850]) isa Vector{Float64}
@test_throws ArgumentError nirs_dpf(; age = -1, wl = 760)
# path length: ppf (or DPF) × source-detector distance (3-D, in mm) → μM
src = n_tmp.locs[findfirst(==(r2[:src_labels][r2[:optode_pairs][5, 1]]), n_tmp.locs.label), :]
det = n_tmp.locs[findfirst(==(r2[:det_labels][r2[:optode_pairs][5, 2]]), n_tmp.locs.label), :]
rho = sqrt((src.loc_x - det.loc_x)^2 + (src.loc_y - det.loc_y)^2 + (src.loc_z - det.loc_z)^2)
n_ppf = od2conc(n_tmp; ppf = [2.0, 2.0])
@test n_ppf.data[n_opt + 1, :, 1] ≈ n_tmp2.data[n_opt + 1, :, 1] ./ (2 * rho)
@test all(==("μM"), n_ppf.header.recording[:unit][(n_opt + 1):(n_opt + 6)])
dpf = nirs_dpf(; age = 30, wl = nirs.header.recording[:wavelengths])
@test od2conc(n_tmp; age = 30).data ≈ od2conc(n_tmp; ppf = dpf).data
@test_throws ArgumentError od2conc(n_tmp; ppf = [1.0, 1.0], age = 30)

@info "Test: nirs_distance() / nirs_sci() / nirs_cv() / nirs_quality()"
@test nirs_distance(nirs).d ≈ fill(30.0, 4)
sc = nirs_sci(nirs)
@test sc.pair == ["S1_D1", "S1_D2"] && all(0.9 .< sc.sci .<= 1)
cvn = nirs_cv(nirs)
@test size(cvn.cv) == (4, 1) && all(0 .< cvn.cv .< 15)
q = nirs_quality(nirs)
@test q.pair == ["S1_D1", "S1_D2"] && !any(q.bad)
@test nirs_quality(nirs; sci_t = 0.99).bad == [true, true]
# a pair whose wavelengths do not share the cardiac pulsation (one replaced by noise)
n_q = deepcopy(nirs)
n_q.data[3, :, 1] .= NeuroAnalyzer.Statistics.mean(n_q.data[3, :, 1]) .+ randn(NeuroAnalyzer.Random.MersenneTwister(1), size(n_q.data, 2))
@test nirs_sci(n_q).sci[1] < 0.3

@info "Test: tddr()"
n_m = deepcopy(n_tmp)
# baseline shift (motion) in an OD channel: removed by TDDR
n_m.data[5, 4000:end, 1] .+= 1.0
n_t = tddr(n_m; ch = "S1_D1 690 OD")
@test abs(NeuroAnalyzer.Statistics.mean(n_t.data[5, 4100:4500, 1]) - NeuroAnalyzer.Statistics.mean(n_t.data[5, 3500:3900, 1])) < 0.2
@test n_t.data[6, :, :] == n_m.data[6, :, :]
n_t2 = deepcopy(n_m)
tddr!(n_t2; ch = "S1_D1 690 OD")
@test n_t2.data ≈ n_t.data

# simulated task and short-separation channel: S1_D2 is moved 8 mm from S1
function _nirs_sim(n0)
    n = deepcopy(n0)
    s = findfirst(==("S1"), n.locs.label)
    d = findfirst(==("D2"), n.locs.label)
    n.locs[d, :loc_x] = n.locs[s, :loc_x]
    n.locs[d, :loc_y] = n.locs[s, :loc_y] + 8
    fs = sr(n)
    onsets = 30:45:330
    for o in onsets
        add_marker!(n; type = "stim", start = Float64(o), len = 15.0, value = "task")
    end
    h = nirs_hrf(; fs = fs).h
    u = zeros(epoch_len(n))
    for o in onsets
        u[(round(Int, o * fs) + 1):round(Int, (o + 15) * fs)] .= 1
    end
    reg = NeuroAnalyzer.DSP.conv(u, h)[1:epoch_len(n)]
    ho, hr, ht = get_channel(n; ch = ["S1_D1 HbO", "S1_D1 HbR", "S1_D1 HbT"])
    a = NeuroAnalyzer.Statistics.std(n.data[ho, :, 1])
    n.data[ho, :, 1] .+= a .* reg
    n.data[hr, :, 1] .-= a / 3 .* reg
    n.data[ht, :, 1] .+= 2a / 3 .* reg
    return n, a
end
n_sim, sim_amp = _nirs_sim(n_tmp2)

@info "Test: nirs_hrf()"
hrf = nirs_hrf(; fs = 10)
@test sum(hrf.h) ≈ 1
@test 4 < hrf.t[argmax(hrf.h)] < 6

@info "Test: nirs_glm()"
g = nirs_glm(n_sim)
@test g.channel == ["S1_D1 HbO", "S1_D1 HbR", "S1_D2 HbO", "S1_D2 HbR"]
@test g.p[1] < 1e-6 && g.beta[1] > 0 && g.p[2] < 0.05 && g.beta[2] < 0
@test g.beta[1] ≈ sim_amp rtol = 0.3
@test g.p[3] > 0.01
@test_throws ArgumentError nirs_glm(n_tmp2)
# several marker values modeled as one condition
n_c = deepcopy(n_sim)
n_c.markers[1:2:end, :value] .= "task_b"
g_c = nirs_glm(n_c; condition = Dict("task_b" => "task"))
@test g_c.beta ≈ g.beta && all(==("task"), g_c.condition)
@test nrow(nirs_glm(n_c)) == 2 * nrow(g)
@test_throws ArgumentError nirs_glm(n_c; cond = "x")
@test_throws ArgumentError nirs_glm(n_sim; marker_type = "response")

@info "Test: nirs_ssr()"
n_s = nirs_ssr(n_sim)
# the short channels are not changed, the long channels are
@test n_s.data[get_channel(n_s; ch = "S1_D2 HbO"), :, :] == n_sim.data[get_channel(n_sim; ch = "S1_D2 HbO"), :, :]
@test n_s.data[get_channel(n_s; ch = "S1_D1 HbO"), :, :] != n_sim.data[get_channel(n_sim; ch = "S1_D1 HbO"), :, :]
# the systemic signal is removed: smaller standard error of the activation
g_s = nirs_glm(n_s)
@test g_s.se[1] < g.se[1]
@test_throws ArgumentError nirs_ssr(n_tmp2)
n_s2 = deepcopy(n_sim)
nirs_ssr!(n_s2)
@test n_s2.data == n_s.data

@info "Test: average_epochs() (NIRS)"
e_n = epoch(n_s; marker_value = "task", offset = 5, ep_len = 30)
a_n = average_epochs(e_n; bl = (-5, 0))
@test datatype(a_n) == "nirs"
@test nepochs(a_n) == nepochs(e_n) + 1
@test !("nirs_aux" in a_n.header.recording[:channel_type])
ho = get_channel(a_n; ch = "S1_D1 HbO")[1]
@test maximum(a_n.data[ho, :, 1]) > 0.5 * sim_amp

@info "Test: export_snirf()"
f_snirf = joinpath(mktempdir(), "test.snirf")
export_snirf(n_s; file_name = f_snirf)
@test_throws ArgumentError export_snirf(n_s; file_name = f_snirf)
r_s = import_snirf(f_snirf)
# written with HDF5.jl: strings are variable-length (SNIRF specification)
@test NeuroAnalyzer.HDF5.h5open(f_snirf, "r") do h
    all(
        k -> NeuroAnalyzer.HDF5.API.h5t_is_variable_str(NeuroAnalyzer.HDF5.datatype(h[k])),
        ("formatVersion", "nirs/metaDataTags/SubjectID", "nirs/probe/sourceLabels"),
    )
end
@test labels(r_s) == labels(n_s)
@test r_s.header.recording[:channel_type] == n_s.header.recording[:channel_type]
@test r_s.data ≈ n_s.data
@test r_s.header.recording[:optode_pairs] == n_s.header.recording[:optode_pairs]
@test r_s.markers[!, :start] == n_s.markers[!, :start] && all(==("task"), r_s.markers[!, :value])
@test_throws ArgumentError export_snirf(e_n; file_name = joinpath(mktempdir(), "e.snirf"))

@info "Test: npl()"
e = average_epochs(e10)
npl!(e)
@test size(e.data) == (19, 2560, 11)

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
e_tmp =
    ica_remove(eeg; ch = "all", ic_idx = Int64[], ic = ica_data.ic, ic_mw = ica_data.ic_mw)
@test e_tmp.data == eeg.data
# keeping all components keeps the data
e_tmp = ica_remove(
    eeg;
    ch = "all",
    ic_idx = 1:5,
    keep = true,
    ic = ica_data.ic,
    ic_mw = ica_data.ic_mw,
)
@test e_tmp.data ≈ eeg.data

@info "Test: ica_decompose(method = :infomax)"
# known sources: 2 super-Gaussian (Laplace), 1 sub-Gaussian (uniform), 1 sinusoid (sub-Gaussian),
# mixed into 6 channels
rng_ica = NeuroAnalyzer.Random.MersenneTwister(7)
n_ica = 20_000
lap_ica() = (u = rand(rng_ica, n_ica) .- 0.5; -sign.(u) .* log.(1 .- 2 .* abs.(u)))
src_ica = permutedims(hcat(lap_ica(), lap_ica(), rand(rng_ica, n_ica) .- 0.5,
                           sin.(2π * 50 .* (0:(n_ica - 1)) ./ 500)))
x_ica = randn(rng_ica, 6, 4) * src_ica
best_r(ic) = [maximum(abs(NeuroAnalyzer.Statistics.cor(ic[i, :], src_ica[j, :])) for i in axes(ic, 1)) for j in 1:4]
# extended Infomax separates sub- and super-Gaussian sources (PCA reduction 6 → 4)
r_ica = NeuroAnalyzer._ica_infomax(x_ica; n = 4)
@test size(r_ica.ic) == (4, n_ica) && size(r_ica.ic_mw) == (6, 4)
@test all(best_r(r_ica.ic) .> 0.99)
# the components reconstruct the (demeaned) data
@test r_ica.ic_mw * r_ica.ic ≈ x_ica .- NeuroAnalyzer.Statistics.mean(x_ica; dims = 2)
# logistic Infomax separates only the super-Gaussian sources
r_ica = NeuroAnalyzer._ica_infomax(x_ica; n = 4, extended = false)
@test all(best_r(r_ica.ic)[1:2] .> 0.99) && minimum(best_r(r_ica.ic)[3:4]) < 0.9
# reproducible
@test NeuroAnalyzer._ica_infomax(x_ica; n = 4).ic == NeuroAnalyzer._ica_infomax(x_ica; n = 4).ic
# rank: 6 channels mixing 4 sources
@test NeuroAnalyzer._data_rank(x_ica) == 4
@test_throws ArgumentError NeuroAnalyzer._ica_infomax(x_ica; n = 5)
ica_data = ica_decompose(eeg; ch = "all", n = 5, method = :infomax, iter = 10)
@test size(ica_data.ic) == (5, 308480) && size(ica_data.ic_mw) == (24, 5)
@test issorted(ica_data.ic_var; rev = true)
# rank-deficient data: the number of components is reduced to the rank
e_rk = keep_channel(eeg; ch = ["Fp1", "Fp2", "F3"])
e_rk.data[3, :, 1] = e_rk.data[1, :, 1] .+ e_rk.data[2, :, 1]
ica_data = ica_decompose(e_rk; ch = "all", method = :infomax, iter = 10)
@test size(ica_data.ic, 1) == 2
@test_throws ArgumentError ica_decompose(eeg; ch = "all", n = 5, method = :x)

@info "Test: ica_decompose(method = :amica)"
@test_throws ArgumentError ica_decompose(
    e_rk;
    ch = "all",
    method = :amica,
    amica = "no_such_amica_program",
)
# the AMICA program must be installed by the user: tested if amica15ub is in the PATH
if isnothing(Sys.which("amica15ub"))
    @info "amica15ub not found: AMICA run not tested"
else
    e_am = create_object(; data_type = "eeg")
    e_am.header.recording[:sampling_rate] = 500
    for c in 1:6
        append_channel!(
            e_am;
            data = reshape(x_ica[c, :], 1, n_ica, 1),
            label = "E$c",
            type = "eeg",
            unit = "μV",
        )
    end
    create_time!(e_am; fs = 500)
    am_dir = joinpath(mktempdir(), "amicaout")
    ica_data = ica_decompose(
        e_am;
        ch = "all",
        n = 4,
        method = :amica,
        iter = 300,
        amica_threads = 2,
        amica_outdir = am_dir,
    )
    @test size(ica_data.ic) == (4, n_ica) && size(ica_data.ic_mw) == (6, 4)
    @test all(best_r(ica_data.ic) .> 0.98)
    # the AMICA output is kept and gives the same components with import_amica()
    @test isfile(joinpath(am_dir, "W"))
    @test import_amica(e_am; ch = "all", path = am_dir).ic ≈ ica_data.ic
end

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
@test length(e10_tmp.header.recording[:lowpass]) == length(e10_tmp.header.recording[:highpass]) == 23

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
# 20 Hz + 4 Hz + linear trend: two IMFs and the trend as the residual
t_emd = collect(0:(1 / 256):(4 - 1 / 256))
x_emd = 5 .* sin.(2π * 20 .* t_emd) .+ 10 .* sin.(2π * 4 .* t_emd) .+ 2 .* t_emd
imf_emd = NeuroAnalyzer._emd(x_emd, t_emd)
@test size(imf_emd) == (3, length(t_emd))
@test vec(sum(imf_emd; dims = 1)) ≈ x_emd
@test NeuroAnalyzer.Statistics.cor(imf_emd[1, :], 5 .* sin.(2π * 20 .* t_emd)) > 0.99
@test NeuroAnalyzer.Statistics.cor(imf_emd[2, :], 10 .* sin.(2π * 4 .* t_emd)) > 0.99
@test maximum(abs.(imf_emd[3, :] .- 2 .* t_emd)) < 1
# the rows (IMFs and residual) always sum to the signal
imf_emd = emd(e10; ch = "Fp1", ep = 1)
@test vec(sum(imf_emd; dims = 1)) ≈ e10.data[1, :, 1]
@test size(emd(e10; ch = "Fp1", ep = 1, max_imf = 2), 1) == 3
# a trend (no oscillation) is returned as the residual only
@test size(NeuroAnalyzer._emd(t_emd, t_emd), 1) == 1
@test_throws ArgumentError NeuroAnalyzer._emd(x_emd, t_emd; epsilon = 0)
@test_throws ArgumentError NeuroAnalyzer._emd(x_emd, t_emd[1:10])

@info "Test: process regressions"
# resample: time axes use the new sampling rate; epoch offset is kept; markers unchanged
r = resample(eeg; new_sr = 128)
@test sr(r) == 128 && r.time_pts[2] == 1 / 128 && length(r.time_pts) == size(r.data, 2)
@test r.markers == eeg.markers
e_off = deepcopy(e10);
e_off.epoch_time .-= 1.0
@test resample(e_off; new_sr = 512).epoch_time[1] == -1.0
@test resample(e10; new_sr = sr(e10)) !== e10
@test_throws ArgumentError downsample(e10; new_sr = 512)
@test_throws ArgumentError upsample(e10; new_sr = 128)
# average_epochs: markers shifted by one epoch (in seconds); in-place version keeps locs in sync
e_avg = average_epochs(e10)
isempty(e10.markers) || @test e_avg.markers[!, :start] ≈ e10.markers[!, :start] .+ 10
e_tmp = deepcopy(e10);
average_epochs!(e_tmp)
@test e_tmp.locs == e_avg.locs
# sort_epochs: s must be a permutation
@test_throws ArgumentError sort_epochs(e_avg; s = fill(2, nepochs(e_avg) - 1))
# plinterpolate_channel: any channel, all channels kept
e_pl = plinterpolate_channel(e10; ch = "Cz", ep = 1)
@test labels(e_pl) == labels(e10)
cz = findfirst(==("Cz"), labels(e10))
@test e_pl.data[cz, :, 1] != e10.data[cz, :, 1] && e_pl.data[cz, :, 2] == e10.data[cz, :, 2]
# lrinterpolate_channel: deterministic
@test lrinterpolate_channel(e10; ch = "Cz", ep = 1).data ==
      lrinterpolate_channel(e10; ch = "Cz", ep = 1).data
# standardize!, tconv! with a complex kernel
e_tmp = deepcopy(e10)
@test standardize!(e_tmp; ch = "Fp1") isa Vector
@test tconv!(deepcopy(e10); ch = "Fp1", kernel = ComplexF64[0.2 + 0.1im, 0.5, 0.2]) isa
      Array{ComplexF64, 3}
# normalize: :fisher available; rank keeps shape; softmax does not overflow
@test size(NeuroAnalyzer._normalize(rand(3, 4); method = :fisher)) == (3, 4)
@test size(normalize_rank(ones(2, 3))) == (2, 3)
@test sum(normalize_softmax([1000.0, 1001.0])) ≈ 1
# weighted median with equal weights is the median
X = rand(6, 50)
@test NeuroAnalyzer._wmedian_cols(X, ones(6)) ≈
      vec(NeuroAnalyzer.Statistics.median(X; dims = 1))

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
@test haskey(cal, NeuroAnalyzer._sss_name("MEG 0113")) &&
      NeuroAnalyzer._sss_name("MEG0113") == "MEG113"
mn, md = _mf_snr(
    maxwell_filter(
        raw;
        origin = o,
        regularize = :none,
        calibration = joinpath(sss_path, "sss_cal_3053.dat"),
    ),
    "test_move_anon_fineCal_raw_sss.fif",
)
@test md > 50
mn, md = _mf_snr(
    maxwell_filter(
        erm;
        regularize = :none,
        cross_talk = joinpath(sss_path, "ct_sparse.fif"),
        e_o...,
    ),
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
hp = [
    parse.(Float64, split(l)) for
    l in readlines(joinpath(sss_path, "test_move_anon_raw.pos"))[2:end] if
    !isempty(strip(l))
]
hp = reduce(vcat, permutedims.(hp))
hp[:, 1] .-= raw.header.recording[:first_time]
mn, md = _mf_snr(
    maxwell_filter(raw; origin = o, regularize = :none, head_pos = hp, mc_interp = :zero),
    "test_move_anon_movecomp_raw_sss.fif",
)
@test mn > 10 && md > 100
@test maxwell_filter(raw; origin = o, regularize = :none, head_pos = hp).header.recording[:sss_info][:movement_compensation]
# destination: translation without rotation
mn, md = _mf_snr(
    maxwell_filter(raw; origin = o, regularize = :none, destination = o),
    "test_move_anon_transDefault_raw_sss.fif",
)
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
ref = [
    parse.(Float64, split(l)) for
    l in readlines(joinpath(sss_path, "chpi5_raw_mc.pos"))[2:end] if !isempty(strip(l))
]
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
@test size(p_ssp.U, 2) == 3 &&
      p_ssp.ssp_projectors * p_ssp.ssp_projectors ≈ p_ssp.ssp_projectors
m_p = apply_ssp_projectors(m_ssp)
@test apply_ssp_projectors(m_p).data ≈ m_p.data

##
## values: SciPy 1.17 / MNE-Python 1.13.2 references, synthetic signals, formulas
##

eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
eb = epoch(eeg; ep_len = 10)
keep_epoch!(eb; ep = 1:10)
load_locs!(
    eb;
    file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"),
)
fs = 256
tt = (0:(4 * fs - 1)) ./ fs
xm = sin.(2π * 5 .* tt) .+ 0.5 .* sin.(2π * 60 .* tt) .+ 0.3 .* cos.(2π * 90 .* tt)
# amplitude of a sinusoid at f in the middle of a signal (least squares)
sine_amp(s, f) = (
    mid = 257:768;
    A = [sin.(2π * f .* tt[mid]) cos.(2π * f .* tt[mid])];
    hypot((A \ s[mid])...)
)

@info "Test: filter_create() / _filter_apply() vs SciPy / MNE"
# MNE filter_data(method = "iir", iir_params = (order = 4, ftype = "butter", output = "sos"), phase = "zero")
f = filter_create(; fprototype = :butterworth, ftype = :lp, cutoff = 30, fs = fs, order = 4)
@test NeuroAnalyzer._filter_apply(xm; flt = f)[[401, 513, 601]] ≈
      [-0.9244172619059996, 6.06058617302202e-07, -0.98116451395666] atol = 1e-12
f = filter_create(; fprototype = :butterworth, ftype = :hp, cutoff = 40, fs = fs, order = 4)
@test NeuroAnalyzer._filter_apply(xm; flt = f)[[401, 513, 601]] ≈
      [-0.704915844879155, 0.2999917628471555, -0.07129865301596978] atol = 1e-12
# scipy.signal.firwin(101, 30, window = "hamming", fs = 256)
f = filter_create(; fprototype = :fir, ftype = :lp, cutoff = 30, fs = fs, order = 101)
@test length(f) == 101
@test f[1:3] ≈ [-0.00039415085773660165, -0.0005255662713842926, -0.0003925985049683375] atol =
    1e-15
@test sum(f) ≈ 1.0
# scipy.signal.remez(101, [0, 28, 32, 128], [1, 0], fs = 256)
f = filter_create(;
    fprototype = :remez,
    ftype = :lp,
    cutoff = 30,
    fs = fs,
    order = 101,
    bw = 4,
)
@test f[1:3] ≈ [-0.0010783991016224733, -0.011648688737594769, -0.0011394905586340003] atol =
    1e-12
# scipy.signal.iirnotch(50, 25, fs = 256)
f = filter_create(; fprototype = :iirnotch, cutoff = 50, fs = fs, bw = 2)
@test [f.b0, f.b1, f.b2] ≈ [0.9760395733504627, -0.6576356575420848, 0.9760395733504627]
@test [f.a1, f.a2] ≈ [-0.6576356575420848, 0.9520791467009253]
# |H(f)| of scipy.signal.cheby1(4, 0.5, 30), cheby2(4, 20, 30), ellip(4, 0.5, 20, 30) at 5, 20, 29, 31, 45, 80 Hz
fr = [5.0, 20, 29, 31, 45, 80]
for (p, ref) in (
    (
        :chebyshev1,
        [
            0.9629080387446526,
            0.9490149558771923,
            0.9868240257223284,
            0.867296802722799,
            0.08728069796786471,
            0.0016906547989606268,
        ],
    ),
    (
        :chebyshev2,
        [
            0.9999996629508234,
            0.9383994834987681,
            0.1648264095457529,
            0.04838851547390617,
            0.09073270752540165,
            0.050573680748433526,
        ],
    ),
    (
        :elliptic,
        [
            0.9541360941244506,
            0.9824598868041596,
            0.9999988024626586,
            0.7061632202106207,
            0.07599550817942868,
            0.06624646323565964,
        ],
    ))
    local f = filter_create(;
        fprototype = p,
        ftype = :lp,
        cutoff = 30,
        fs = fs,
        order = 4,
        rp = 0.5,
        rs = 20,
    )
    @test abs.(NeuroAnalyzer.DSP.freqresp(f, 2π .* fr ./ fs)) ≈ ref atol = 1e-10
end
@test_throws ArgumentError filter_create(;
    fprototype = :butterworth,
    ftype = :lp,
    cutoff = 30,
    fs = fs,
)
@test_throws ArgumentError filter_create(;
    fprototype = :firls,
    ftype = :lp,
    cutoff = 30,
    fs = fs,
    order = 101,
)
@test_throws ArgumentError filter_create(;
    fprototype = :fir,
    ftype = :hp,
    cutoff = 30,
    fs = fs,
    order = 100,
)
@test_throws ArgumentError filter_create(;
    fprototype = :iirnotch,
    cutoff = 127,
    fs = fs,
    bw = 4,
)

@info "Test: filter() responses (all prototypes and types)"
# synthetic mixture of 5, 30 and 60 Hz sines: pass bands kept, stop bands removed (two-pass)
xs = sin.(2π * 5 .* tt) .+ sin.(2π * 30 .* tt) .+ sin.(2π * 60 .* tt)
designs = Dict(
    :lp => (40, [5], [60]),
    :hp => (20, [30, 60], [5]),
    :bp => ((20, 40), [30], [5, 60]),
    :bs => ((20, 40), [5, 60], [30]),
)
for p in (:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :firls, :remez),
    (ft, (c, keep, stop)) in designs

    kw =
        p in (:butterworth, :chebyshev1, :chebyshev2, :elliptic) ?
        (order = 6, rp = 0.5, rs = 40) :
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
# DSP.filtfilt(b, x) filtered in place and was corrupted when filt! used FFT overlap-save (DSP.jl ≤ 0.8.6)
# reference: forward and reverse convolution (the middle is unaffected by the edge extension)
xl = randn(NeuroAnalyzer.Random.MersenneTwister(1), 100_000) .+ 100
f = filter_create(; fprototype = :fir, ftype = :hp, cutoff = 0.5, fs = fs, bw = 0.5)
fc(v) =
    NeuroAnalyzer.DSP.conv(v, f)[((length(f) - 1) ÷ 2 + 1):((length(f) - 1) ÷ 2 + length(
        v,
    ))]
yl = NeuroAnalyzer._filter_apply(xl; flt = f)
@test yl[20_000:80_000] ≈ reverse(fc(reverse(fc(xl))))[20_000:80_000] atol = 1e-9
@test_throws ArgumentError NeuroAnalyzer._filter_apply(xl[1:(length(f) - 1)]; flt = f)

@info "Test: filter() / filter_apply() / filter!() on objects"
f = filter_create(; fprototype = :butterworth, ftype = :lp, cutoff = 30, fs = fs, order = 4)
e_f = filter_apply(eb; ch = "Fp1", flt = f, report = false)
@test e_f.data[1, :, 3] ≈ NeuroAnalyzer._filter_apply(eb.data[1, :, 3]; flt = f)
@test e_f.data[2, :, :] == eb.data[2, :, :]
e_tmp = deepcopy(eb)
filter_apply!(e_tmp; ch = "Fp1", flt = f, report = false)
@test e_tmp.data == e_f.data
@test NeuroAnalyzer._filter_one(eb, 1, 3, f, :twopass) ≈ e_f.data[1, :, 3]
@test NeuroAnalyzer._filter_one(eb, 1, 3, f, :onepass) ≈
      NeuroAnalyzer.DSP.filt(f, eb.data[1, :, 3])
e_f2 = NeuroAnalyzer.filter(
    eb;
    ch = "Fp1",
    fprototype = :butterworth,
    ftype = :lp,
    cutoff = 30,
    order = 4,
    report = false,
)
@test e_f2.data ≈ e_f.data
e_tmp = deepcopy(eb)
NeuroAnalyzer.filter!(
    e_tmp;
    ch = "Fp1",
    fprototype = :butterworth,
    ftype = :lp,
    cutoff = 30,
    order = 4,
    report = false,
)
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
rep = filter_report(
    f;
    fs = fs,
    ftype = :lp,
    cutoff = 30,
    bw = 10,
    fcheck = [5.0, 30.0, 60.0],
    verbose = false,
)
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
@test c.data[1, [1, 2, 101, 512], 1] ≈
      [-1692.0568207927302, -1775.5717618275437, 659.987805703435, 603.4875914557678] rtol =
    1e-9
@test c.data[10, [1, 2, 101, 512], 1] ≈
      [-1962.2371848290109, -545.3811613543534, -136.083974881286, -1501.146212469763] rtol =
    1e-9
@test c.header.recording[:data_type] == "csd"
@test c.header.recording[:unit][1] == "μV/m²"
# electrode positions are projected onto the unit sphere: the radius does not matter
e_tmp = deepcopy(e2)
e_tmp.locs[!, :loc_x] .*= 2;
e_tmp.locs[!, :loc_y] .*= 2;
e_tmp.locs[!, :loc_z] .*= 2
@test csd(e_tmp).data ≈ c.data
G, H = gh(e2.locs)
@test G ≈ G' && H ≈ H'
@test_throws ArgumentError csd(e2; m = 1)

@info "Test: reference values"
ch_eeg = get_channel(eb; ch = get_channel(eb; type = "eeg"))
e_r = reference_avg(eb)
@test e_r.data[ch_eeg, :, :] ≈
      eb.data[ch_eeg, :, :] .- sum(eb.data[ch_eeg, :, :]; dims = 1) ./ length(ch_eeg)
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
e_r = reference_custom(e10; ref_list = ["Fp1-F3", "Cz"])
lfp1 = e10.locs[NeuroAnalyzer._find_bylabel(e10.locs, "Fp1"), :]
lf3 = e10.locs[NeuroAnalyzer._find_bylabel(e10.locs, "F3"), :]
lm = e_r.locs[NeuroAnalyzer._find_bylabel(e_r.locs, "Fp1-F3"), :]
pm = [lm.loc_x, lm.loc_y, lm.loc_z]
p1 = [lfp1.loc_x, lfp1.loc_y, lfp1.loc_z]
p2 = [lf3.loc_x, lf3.loc_y, lf3.loc_z]
@test LinearAlgebra.norm(pm) ≈ (LinearAlgebra.norm(p1) + LinearAlgebra.norm(p2)) / 2
@test LinearAlgebra.normalize(pm) ≈ LinearAlgebra.normalize(p1 .+ p2)
@test lm.loc_radius_sph ≈ LinearAlgebra.norm(pm)
@test e_r.locs[NeuroAnalyzer._find_bylabel(e_r.locs, "Cz"), :loc_x] ==
      e10.locs[NeuroAnalyzer._find_bylabel(e10.locs, "Cz"), :loc_x]
@test isempty(NeuroAnalyzer._find_bylabel(e_r.locs, "Fp1"))
e_r = reference_custom(e10; ref_list = ["Fp1-F3", "Cz"], loc = :first)
lm = e_r.locs[NeuroAnalyzer._find_bylabel(e_r.locs, "Fp1-F3"), :]
@test lm.loc_x == lfp1.loc_x && lm.loc_theta == lfp1.loc_theta
@test_throws ArgumentError reference_custom(e10; loc = :last)
@test NeuroAnalyzer._parse_montage_entry("Fp1-F3", ["Fp1", "F3"]) == ("Fp1", "F3")
@test NeuroAnalyzer._parse_montage_entry("Fp1", ["Fp1", "F3"]) == ("Fp1", nothing)
@test NeuroAnalyzer._parse_montage_entry("A-B-C", ["A-B", "C"]) == ("A-B", "C")
@test_throws ArgumentError NeuroAnalyzer._parse_montage_entry("X-Y", ["Fp1"])
@test NeuroAnalyzer._eeg_ch_idx(eb) == ch_eeg
@test NeuroAnalyzer._inv_dist_weights([0.0 2.0; Inf 0.5]) == [0.0 0.5; 0.0 2.0]
@test NeuroAnalyzer._lap_suffix("plap", true, true) == "-wplapmed"
@test NeuroAnalyzer._lap_name("planar", 4, true, true) ==
      "weighted planar Laplacian (median, 4)"
lo = NeuroAnalyzer._locs_ch_order(eb, [3, 1])
@test lo[!, :label] == ["F3", "Fp1"]
d = NeuroAnalyzer._locs_gc_dist(lo)
@test d[1, 1] ≈ 0 atol = 1e-7
@test d[1, 2] ≈ acos(
    clamp(
        sum(
            collect(lo[1, [:loc_x, :loc_y, :loc_z]]) .*
            collect(lo[2, [:loc_x, :loc_y, :loc_z]]),
        ) /
        (
            NeuroAnalyzer.norm(collect(lo[1, [:loc_x, :loc_y, :loc_z]])) *
            NeuroAnalyzer.norm(collect(lo[2, [:loc_x, :loc_y, :loc_z]]))
        ), -1, 1),
)
px, py = NeuroAnalyzer._locs_polar_xy(lo)
@test px ≈ lo.loc_radius .* cosd.(lo.loc_theta) && py ≈ lo.loc_radius .* sind.(lo.loc_theta)
ch_l, s_l = NeuroAnalyzer._laplacian_reference(
    eb;
    nn = 3,
    weighted = false,
    med = false,
    spherical = true,
)
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
@test r.history[end] == "resample(obj; new_sr=128)"

@info "Test: detect_powerline()"
x =
    0.2 .* randn(NeuroAnalyzer.Random.MersenneTwister(3), 2560) .+
    sin.(2π * 50 .* (0:2559) ./ 256)
@test NeuroAnalyzer._detect_powerline(x; fs = 256) == 50.0
e_tmp = keep_epoch(eb; ep = 1)
e_tmp.data[1, :, 1] = x
@test detect_powerline(keep_channel(e_tmp; ch = "Fp1")) == 50.0
# only 50 or 60 Hz: a stronger 10 Hz rhythm is not power-line noise; 0 if there is none
tt_pl = (0:2559) ./ 256
x_pl =
    5 .* sin.(2π * 10 .* tt_pl) .+
    0.5 .* randn(NeuroAnalyzer.Random.MersenneTwister(4), 2560)
@test NeuroAnalyzer._detect_powerline(x_pl; fs = 256) == 0.0
@test NeuroAnalyzer._detect_powerline(x_pl .+ 0.5 .* sin.(2π * 60 .* tt_pl); fs = 256) ==
      60.0
@test NeuroAnalyzer._detect_powerline(x_pl .+ 0.5 .* sin.(2π * 50 .* tt_pl); fs = 256) ==
      50.0
@test NeuroAnalyzer._detect_powerline(
    x_pl .+ 0.5 .* sin.(2π * 50 .* tt_pl);
    fs = 256,
    pr = 100,
) == 0.0
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
@test_throws ArgumentError mlinterpolate_channel(
    eb;
    ch = "Cz",
    ep = 1,
    ep_ref = 1:3,
    model = cr,
)
@test_throws ArgumentError mlinterpolate_channel(eb; ch = "ECG", ep = 1, model = cr)
e_tmp = deepcopy(eb)
e_tmp.header.recording[:bad_channel][get_channel(eb; ch = "Cz")[1]] = true
@test !mlinterpolate_channel(
    e_tmp;
    ch = "Cz",
    ep = 1,
    ep_ref = 2:4,
    model = cr,
    reset_bads = true,
).header.recording[:bad_channel][get_channel(eb; ch = "Cz")[1]]
mlinterpolate_channel!(
    e_tmp;
    ch = "Cz",
    ep = 1,
    ep_ref = 2:4,
    model = cr,
    reset_bads = true,
)
@test !e_tmp.header.recording[:bad_channel][get_channel(eb; ch = "Cz")[1]]

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
        @test edit_montage(joinpath(NeuroAnalyzer.PATH, "montages", "bip_long.mnt")) ===
              nothing
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
    e_k = zeros(3)
    e_k[k] = h
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
U = NeuroAnalyzer._tsss_projector(
    0.01 .* randn(NeuroAnalyzer.Random.MersenneTwister(1), 5, 200) .+
    [1.0, 2, 3, 4, 5] * tk_',
    0.01 .* randn(NeuroAnalyzer.Random.MersenneTwister(2), 5, 200) .+
    [5.0, 4, 3, 2, 1] * tk_', 0.98)
@test size(U, 2) >= 1
@test abs(NeuroAnalyzer.norm(U' * tk_) / NeuroAnalyzer.norm(tk_)) > 0.9
# coil integration points: weights of a magnetometer sum to its area normal
meg_raw = import_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fif"))
mr = meg_raw.header.recording
pts, wn, cidx = NeuroAnalyzer._sss_coil_points(
    mr[:coil_loc][mr[:magnetometers][1:2], :],
    mr[:coil_id][mr[:magnetometers][1:2]],
    NeuroAnalyzer.Matrix{Float64}(NeuroAnalyzer.I(4)),
)
@test size(pts, 2) == 3 && size(wn) == size(pts) && sort(unique(cidx)) == [1, 2]
@test size(NeuroAnalyzer._vv_mag_acc(0.0258), 2) == 4
@test size(NeuroAnalyzer._vv_mag_normal(0.0258), 2) == 4
cp_ = NeuroAnalyzer._chpi_coil_points(
    mr[:coil_loc][mr[:magnetometers][1:30], :],
    mr[:coil_id][mr[:magnetometers][1:30]],
)
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
pts, wn, cidx = NeuroAnalyzer._sss_coil_points(
    mr[:coil_loc][sel, :],
    mr[:coil_id][sel],
    NeuroAnalyzer.Matrix{Float64}(NeuroAnalyzer.I(4)),
)
S_in, S_out = NeuroAnalyzer._sss_basis(pts, wn, cidx, length(sel), [0.0, 0.0, 0.0], 3, 2)
@test size(S_in) == (40, NeuroAnalyzer._sss_ncomp(3))
@test size(S_out) == (40, NeuroAnalyzer._sss_ncomp(2))
b =
    NeuroAnalyzer._chpi_mag_dipole_fwd([0.0, 0.0, 0.0], pts, wn, cidx, length(sel)) *
    [1e-8, -2e-8, 0.5e-8]
B1 = S_in[:, 1:3]
@test NeuroAnalyzer.norm(b .- B1 * (B1 \ b)) / NeuroAnalyzer.norm(b) < 1e-6
# an off-origin dipole needs higher orders
b2 =
    NeuroAnalyzer._chpi_mag_dipole_fwd([0.0, 0.02, 0.04], pts, wn, cidx, length(sel)) *
    [1e-8, 0.0, 0.0]
@test NeuroAnalyzer.norm(b2 .- B1 * (B1 \ b2)) / NeuroAnalyzer.norm(b2) > 1e-3
# information-based regularization returns a subset of the internal components
keep = NeuroAnalyzer._sss_regularize_in(
    hcat(S_in, S_out),
    size(S_in, 2),
    3,
    collect(1:size(S_out, 2)),
)
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
@test NeuroAnalyzer._sss_avg_trans([0], [T1], [[0.1, 0.0, 0.0, 0.01, 0.0, 0.0]], 10, 20) ==
      T1
Ta = NeuroAnalyzer._sss_avg_trans(
    [0, 15],
    [T1, T1],
    [[0.1, 0.0, 0.0, 0.01, 0.0, 0.0], [0.1, 0.0, 0.0, -0.01, 0.0, 0.0]],
    10,
    20,
)
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
    (denoise_wien, denoise_wien!, eb, (ch = ["Fp1", "Fp2", "F3"],)),
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
    (reference_rest, reference_rest!, eb, (ndip = 500,)),
    (reference_custom, reference_custom!, eb, ()),
    (resample, resample!, eb, (new_sr = 128,)),
    (downsample, downsample!, eb, (new_sr = 128,)),
    (upsample, upsample!, eb, (new_sr = 512,)),
    (lrinterpolate_channel, lrinterpolate_channel!, eb, (ch = "Cz", ep = 1)),
    (plinterpolate_channel, plinterpolate_channel!, eb, (ch = "Cz", ep = 1)),
    (
        artrem_cwd,
        artrem_cwd!,
        eb,
        (ch = "Fp1", ep = 1, tseg = (0.2, 0.4), fseg = (10, 12.5)),
    ),
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
@test e_tmp.data ≈
      ica_remove(eeg; ch = "all", ic_idx = 1:2, ic = ica_d.ic, ic_mw = ica_d.ic_mw).data
e_tmp = keep_epoch(eb; ep = 1)
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
@test NeuroAnalyzer._detrend(x; type = :linear) ≈
      NeuroAnalyzer._detrend(sin.(2π .* (1:100) ./ 20); type = :linear) atol = 1e-10
@test sum(NeuroAnalyzer._remove_dc(x)) ≈ 0 atol = 1e-9
# taper: elementwise product
@test NeuroAnalyzer._taper([1.0, 2.0, 3.0]; t = [0.0, 0.5, 1.0]) == [0.0, 1.0, 3.0]
# moving average of a constant, Savitzky-Golay of a polynomial of its order
@test NeuroAnalyzer._filter_sg(collect(1.0:50) .^ 2; order = 2, window = 11)[10:40] ≈
      (collect(10.0:40)) .^ 2 rtol = 1e-8
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
e_pl.data[1, :, 1] =
    sin.(2π * 10 .* tp_) .+ sin.(2π * 50 .* tp_) .+
    0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(1), length(tp_))
e_out, df_pl = remove_powerline(e_pl; ch = "Fp1", pl_frq = 50)
amp_f(s, f) = (
    mid = 1001:(length(s) - 1000);
    A = [sin.(2π * f .* tp_[mid]) cos.(2π * f .* tp_[mid])];
    hypot((A \ s[mid])...)
)
@test amp_f(e_out.data[1, :, 1], 50) < 0.05
@test amp_f(e_out.data[1, :, 1], 10) > 0.95
@test df_pl[1, :channel] == "Fp1" && df_pl[1, :frequency] == 50 && df_pl[1, :bandwidth] > 0
# the removed frequencies are recorded in :notch
@test e_out.header.recording[:notch][1] == Float64.(df_pl[df_pl.channel .== "Fp1", :frequency])
@test isempty(e_out.header.recording[:notch][2])
# only channels with a peak are filtered; no harmonic in this signal
@test nrow(df_pl) == 1
# other channels are not modified
@test e_out.data[2, :, 1] == e_pl.data[2, :, 1]
# harmonics: 50 Hz + 100 Hz removed; a channel without line noise is not modified
e_pl.data[1, :, 1] =
    sin.(2π * 10 .* tp_) .+ sin.(2π * 50 .* tp_) .+ 0.3 .* sin.(2π * 100 .* tp_) .+
    0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(1), length(tp_))
e_pl.data[2, :, 1] =
    sin.(2π * 10 .* tp_) .+
    0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(2), length(tp_))
e_out, df_pl = remove_powerline(e_pl; ch = ["Fp1", "Fp2"], pl_frq = 50)
@test df_pl[!, :channel] == ["Fp1", "Fp1"] && df_pl[!, :frequency] == [50, 100]
@test amp_f(e_out.data[1, :, 1], 50) < 0.05 && amp_f(e_out.data[1, :, 1], 100) < 0.03
@test amp_f(e_out.data[1, :, 1], 10) > 0.95
@test e_out.data[2, :, 1] == e_pl.data[2, :, 1]
# 60 Hz
e_pl.data[1, :, 1] =
    sin.(2π * 10 .* tp_) .+ sin.(2π * 60 .* tp_) .+
    0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(1), length(tp_))
e_out, df_pl = remove_powerline(e_pl; ch = "Fp1", pl_frq = 60)
@test df_pl[!, :frequency] == [60] && amp_f(e_out.data[1, :, 1], 60) < 0.05
@test_throws ArgumentError remove_powerline(e_pl; ch = "Fp1", pl_frq = 0)
@test_throws ArgumentError remove_powerline(e_pl; ch = "Fp1", pl_frq = 50, d = 0)
@test_throws ArgumentError remove_powerline(epoch(eeg; ep_len = 10); ch = "Fp1")

@info "Test: remove_powerline(method = :zapline)"
# 19 EEG channels, 120 s: low-pass "brain" activity mixed over the channels + line noise
# (50 Hz and 100 Hz) with one spatial pattern
e_zl = keep_epoch(epoch(eeg; ep_len = 120); ep = 1)
ch_zl = get_channel(e_zl; ch = "eeg")
rng_zl = NeuroAnalyzer.Random.MersenneTwister(5)
n_zl = epoch_len(e_zl)
tz = (0:(n_zl - 1)) ./ sr(e_zl)
brain = NeuroAnalyzer._zl_smooth(randn(rng_zl, n_zl, 8), 8) * randn(rng_zl, 8, length(ch_zl))
line = (sin.(2π * 50 .* tz) .+ 0.3 .* sin.(2π * 100 .* tz .+ 1)) * randn(rng_zl, 1, length(ch_zl))
e_zl.data[ch_zl, :, 1] = permutedims(brain .+ 3 .* line)
rel_err(e) = LinearAlgebra.norm(permutedims(e.data[ch_zl, :, 1]) .- brain) / LinearAlgebra.norm(brain)
@test rel_err(e_zl) > 1
e_out, df_zl = remove_powerline(e_zl; ch = "eeg", pl_frq = 50, method = :zapline, chunk = :none)
@test rel_err(e_out) < 0.05
@test nrow(df_zl) == 1 && df_zl[1, :nremove] == 1
@test df_zl[1, :reduction] > 20 && df_zl[1, :prominence] < 2
@test df_zl[1, :start] == 0 && df_zl[1, :stop] == 120
@test occursin("method=:zapline", e_out.history[end])
# other channels are not modified
ch_other = setdiff(1:nchannels(e_zl), ch_zl)
@test e_out.data[ch_other, :, 1] == e_zl.data[ch_other, :, 1]
# fixed chunks; a fixed number of components
_, df_zl = remove_powerline(e_zl; ch = "eeg", method = :zapline, chunk = :fixed, chunk_len = 40)
@test df_zl.start == [0, 40, 80] && all(df_zl.nremove .== 1)
_, df_zl = remove_powerline(e_zl; ch = "eeg", method = :zapline, chunk = :none, nremove = 2)
@test df_zl[1, :nremove] == 2
# adaptive chunks cover the signal
_, df_zl = remove_powerline(e_zl; ch = "eeg", method = :zapline, min_chunk = 20)
@test df_zl[1, :start] == 0 && df_zl[end, :stop] == 120
@test all(df_zl.stop[1:(end - 1)] .== df_zl.start[2:end])
# the in-place version
e_tmp = deepcopy(e_zl)
remove_powerline!(e_tmp; ch = "eeg", method = :zapline, chunk = :none)
@test rel_err(e_tmp) < 0.05
# heavily contaminated: line noise in 6 components (amplitude-modulated sources); the default
# max_prop = 0.2 allows 3 of 19 components, a larger max_prop removes more of the noise
line6 = reduce(hcat, [(1 .+ 0.5 .* sin.(2π * 0.05 * k .* tz)) .* sin.(2π * 50 .* tz .+ k) for k in 1:6]) *
        randn(rng_zl, 6, length(ch_zl))
e_tmp.data[ch_zl, :, 1] = permutedims(brain .+ 3 .* line6)
e_out, df_zl = remove_powerline(e_tmp; ch = "eeg", method = :zapline, chunk = :none)
@test df_zl[1, :nremove] == 3 && df_zl[1, :prominence] > 10
e_out2, df_zl2 = remove_powerline(e_tmp; ch = "eeg", method = :zapline, chunk = :none, max_prop = 0.35)
@test df_zl2[1, :nremove] > 3 && df_zl2[1, :prominence] < df_zl[1, :prominence]
@test_throws ArgumentError remove_powerline(e_zl; ch = "Fp1", method = :zapline)
@test_throws ArgumentError remove_powerline(e_zl; ch = "eeg", method = :zapline, chunk = :x)
@test_throws ArgumentError remove_powerline(e_zl; ch = "eeg", method = :zapline, nremove = 19)
@test_throws ArgumentError remove_powerline(e_zl; ch = "eeg", method = :zapline, nremove = :x)
@test_throws ArgumentError remove_powerline(e_zl; ch = "eeg", method = :x)

@info "Test: remove_powerline(method = :cleanline)"
# one channel: 10 Hz signal + line noise at 50 Hz and 100 Hz (+ a little noise)
e_cl = deepcopy(e_pl)
e_cl.data[1, :, 1] =
    sin.(2π * 10 .* tp_) .+ sin.(2π * 50 .* tp_ .+ 0.4) .+ 0.3 .* sin.(2π * 100 .* tp_) .+
    0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(1), length(tp_))
e_cl.data[2, :, 1] =
    sin.(2π * 10 .* tp_) .+
    0.05 .* randn(NeuroAnalyzer.Random.MersenneTwister(2), length(tp_))
e_out, df_cl = remove_powerline(e_cl; ch = "Fp1", pl_frq = 50, method = :cleanline)
@test amp_f(e_out.data[1, :, 1], 50) < 0.05 && amp_f(e_out.data[1, :, 1], 100) < 0.03
@test amp_f(e_out.data[1, :, 1], 10) > 0.95
@test df_cl.channel == ["Fp1", "Fp1"] && df_cl.frequency == [50, 100]
@test all(df_cl.windows .== 1) && all(df_cl.reduction .> 20)
@test occursin("method=:cleanline", e_out.history[end])
# other channels are not modified
@test e_out.data[2:end, :, 1] == e_cl.data[2:end, :, 1]
# a channel without line noise is (almost) not changed
e_out, df_cl = remove_powerline(e_cl; ch = "Fp2", pl_frq = 50, method = :cleanline)
@test maximum(abs, e_out.data[2, :, 1] .- e_cl.data[2, :, 1]) < 0.05
@test all(df_cl.windows .< 0.2)
# in-place
e_tmp = deepcopy(e_cl)
remove_powerline!(e_tmp; ch = "Fp1", pl_frq = 50, method = :cleanline, wlen = 2, wstep = 0.5)
@test amp_f(e_tmp.data[1, :, 1], 50) < 0.05
@test_throws ArgumentError remove_powerline(e_cl; ch = "Fp1", method = :cleanline, p = 0)
@test_throws ArgumentError remove_powerline(e_cl; ch = "Fp1", method = :cleanline, wstep = 5)
@test_throws ArgumentError remove_powerline(e_cl; ch = "Fp1", method = :cleanline, wlen = 1)
@test_throws ArgumentError remove_powerline(e_cl; ch = "Fp1", method = :cleanline, wlen = 100)
@test_throws ArgumentError remove_powerline(e_cl; ch = "Fp1", method = :cleanline, maxiter = 0)

@info "Test: in-place and out-of-place versions"
e_io = epoch(eeg; ep_len = 10)
keep_epoch!(e_io; ep = 1:5)
load_locs!(
    e_io;
    file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"),
)
# the in-place version modifies the object, the out-of-place version returns the same result
# and does not modify its input
for (f!, f, kw) in (
    (reference_avg!, reference_avg, (;)),
    (reference_ce!, reference_ce, (; ch = "Cz")),
    (reference_a!, reference_a, (; type = :i)),
    (reference_plap!, reference_plap, (;)),
    (remove_dc!, remove_dc, (; ch = "eeg")),
    (
        NeuroAnalyzer.filter!,
        NeuroAnalyzer.filter,
        (; ch = "eeg", fprototype = :fir, ftype = :lp, cutoff = 30, bw = 5),
    ),
    (keep_channel!, keep_channel, (; ch = ["Fp1", "Cz"])),
    (delete_epoch!, delete_epoch, (; ep = 2)),
    (resample!, resample, (; new_sr = 128)),
    (plinterpolate_channel!, plinterpolate_channel, (; ch = "Fp1", ep = 1)),
    (sphinterpolate_channel!, sphinterpolate_channel, (; ch = "Fp1", ep = 1)),
    (edit_channel!, edit_channel, (; ch = "Fp1", field = :label, value = "FP1")),
)
    e_in = deepcopy(e_io)
    local d0
    d0, l0, h0 = copy(e_io.data), copy(labels(e_io)), copy(e_io.history)
    e_oop = f(e_io; kw...)
    @test e_io.data == d0 && labels(e_io) == l0 && e_io.history == h0
    @test f!(e_in; kw...) === nothing
    @test e_in.data == e_oop.data
    @test labels(e_in) == labels(e_oop)
    @test e_in.history == e_oop.history
    @test e_in.markers == e_oop.markers
end
# reference_ce!: a reference channel that is also re-referenced (copied before subtraction)
e_ce = reference_ce(e_io; ch = "Cz")
@test all(iszero, e_ce.data[get_channel(e_ce; ch = "Cz-Cz"), :, :])
@test e_ce.header.recording[:reference] == "common (Cz)"

true
