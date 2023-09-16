using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@test size(covm(zeros(2))) == (2, 2)
@test size(covm(zeros(2, 3, 2))) == (2, 2, 3, 2)
@test covm(ones(2), zeros(2)) == zeros(2, 2)

@test size(corm(rand(2))) == (2, 2)
@test size(corm(rand(2, 3, 2))) == (2, 2, 3, 2)
@test round.(corm(rand(2), rand(2))) == ones(2, 2)

@test linspace(1, 10, 10) == 1.0:10.0
@test logspace(0, 1, 3) == [1.0, 3.1622776601683795, 10.0]
@test length(pad0(ones(3), 1)) == 4
@test size(pad0(ones(3, 3), 1)) == (3, 4)
@test length(pad2(ones(10))) == 16
@test size(pad2(ones(3, 3))) == (3, 4)
@test size(m_pad0(ones(3, 4))) == (4, 4)
@test vsearch([1, 2], [1, 2, 3, 4]) == [1, 2]
@test cart2pol(1, 1) == (1.414, 45.0)
@test pol2cart(1.414, 45.0) == (1.0, 1.0)
@test sph2cart(1, 45, 45) == (0.5, 0.5, 0.707)
@test cart2sph(0.5, 0.5, 0.707) == (1.0, 45.0, 44.991)
@test generate_window(:hann, 3) == [0.0, 1.0, 0.0]
@test length(fft0(ones(10), 10)) == 20
@test length(ifft0(ones(20), 10)) == 10
@test nextpow2(5) == 8
@test vsplit(ones(4), 2) == [[1.0, 1.0], [1.0, 1.0]]
@test length(rms(ones(10))) == 1
@test length(generate_sine(2, collect(1:10))) == 10
@test length(generate_csine(2, collect(1:10))) == 10
@test length(generate_square(collect(1:10), 0.5)) == 10
@test length(generate_triangle(collect(1:10))) == 10
@test length(freqs(rand(10))) == 2
@test m_sortperm([1 2; 4 6]) == [1 1; 2 2]
@test m_sort([1 2; 4 6], [2, 1]) == [4 6; 1 2]
@test hz2rads(1) == 6.283185307179586
@test rads2hz(1) == 1.5707963267948966
@test cmax([1.0im, 2.0im]) == 0.0 + 2.0im
@test cmin([1.0im, 2.0im]) == 0.0 + 1.0im
@test length(generate_sinc(-2:2)) == 5
@test length(generate_morlet(10, 10, ncyc=2)) == 21
@test length(generate_gaussian(10, 10)) == 21
@test tuple_order((2, 1)) == (1, 2)
@test rmse(ones(10), ones(10)) == 0.0
@test size(m_norm(ones(4, 4, 1))) == (4, 4, 1)
@test dft(ones(4), fs=10) == (s_fft = ComplexF64[4.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im], s_sf = [0.0, 0.025, -0.05, -0.025])
@test msci95(ones(4)) == (s_m = 1.0, s_s = 0.0, s_u = 1.0, s_l = 1.0)
@test msci95(ones(4), zeros(4)) == (s_m = 1.0, s_s = 0.0, s_u = 1.0, s_l = 1.0)
@test length(difference(ones(4), zeros(4))) == 3
@test acov(ones(4)) == (acov = [3.0, 4.0, 3.0], lags = [-1, 0, 1])
@test xcov(ones(4), ones(4)) == (xcov = [3.0, 4.0, 3.0], lags = [-1, 0, 1])
@test spectrum(ones(4)) == (s_fft = ComplexF64[4.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im], s_amp = [1.0, 0.0], s_pow = [1.0, 0.0], s_pha = [0.0, 0.0, 0.0, 0.0])
@test total_power(ones(4), fs=10) == 0.0
@test band_power(ones(4), fs=10, f=(1,2)) == 0.0
@test taper(ones(10), t=zeros(10)) == zeros(10)
@test detrend(ones(10)) == zeros(10)
@test remove_dc(ones(10)) == zeros(10)
@test normalize_zscore([1, 2, 3]) == [-1.0, 0.0, 1.0]
@test normalize_minmax([1, 2, 3]) == [-1.0, 0.0, 1.0]
@test normalize_log([0, 0, 0]) == [0.0, 0.0, 0.0]
@test length(add_noise(ones(10), randn(10))) == 10

s, t = resample(ones(10), t=1:10, new_sr=20)
@test t == 1.0:0.05:10.0

@test derivative(ones(10)) == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test NeuroAnalyzer.filter(ones(10), fprototype=:mavg, order=2) == [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]

p, f = psd(ones(100), fs=10)
@test p[1] == 0.0
@test f[end] == 5.0

@test round.(stationarity_hilbert(ones(10))) == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test stationarity_mean(ones(10), window=1) == [1.0]
@test stationarity_var(ones(10), window=1) == [0.0]
@test trim(ones(10), segment=(1,5)) == ones(5)
@test mutual_information(ones(10), ones(10)) == 0.0
@test entropy([1.0, 2.0, 3.0]) == (ent = 1.5849625007211552, sent = 0.8304717124362917, leent = 4.333653050389665)
@test negentropy([1, 2, 3]) == -0.16602396751648252
@test average(ones(10, 10, 1)) == ones(1, 10, 1)
@test average(ones(5, 5, 1), zeros(5, 5, 1)) == [0.5; 0.5; 0.5; 0.5; 0.5;;;]
@test tcoherence([1, 2], [3, 4]) == (c = [5.25, 0.25], msc = [27.5625, 0.0625], ic = [-0.0, 0.0])

pc, pc_w, pc_m, pc_model = pca(ones(2, 10, 1), n=1)
@test pc == [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;;;]

s = pca_reconstruct(ones(2, 10, 1), pc=pc, pc_model=pc_model)
@test s == ones(2, 10, 1)

@test tconv(ones(5), kernel=[1.0, 2.0]) == [1.0, 3.0, 3.0, 3.0, 3.0]
@test round.(fconv(ones(5), kernel=[1.0, 2.0])) == [0.0, 1.0, 1.0, 1.0, 1.0]

i, m = ica([1.0 2.0; 3.0 4.0;;;], n=1)
@test size(i) == (1, 2, 1)
@test ica_reconstruct([1.0 2.0; 3.0 4.0;;;], ic=i, ic_mw=m, ic_idx=[1]) == zeros(2, 2, 1)

p, f, t = spectrogram(ones(100), fs=10)
@test size(p) == (51, 46)

@test detect_channel_flat(ones(100)) == 0.9583333333333334
@test snr(ones(10)) == Inf
@test snr2(ones(10)) == 0.0
@test findpeaks(repeat([0, 1], 100)) == [6, 38, 70, 102, 134, 166, 198]
@test length(denoise_wavelet(rand(100), wt=wavelet(WT.haar))) == 100
@test ispc([1.0, 1.0, 1.0], [0.0, 0.0, 0.0]) == (ispc_value = 1.0, ispc_angle = 0.0, signal_diff = [-1.0, -1.0, -1.0], phase_diff = [0.0, 0.0, 0.0], s1_phase = [0.0, 0.0, 0.0], s2_phase = [0.0, 0.0, 0.0])
@test itpc(ones(1, 10, 10), t=1) == (itpc_value = 1.0, itpcz = 10.0, itpc_angle = 0.0, itpc_phases = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
@test pli([1.0, 1.0, 1.0], [0.0, 0.0, 0.0]) == (pli_value = 0.0, signal_diff = [-1.0, -1.0, -1.0], phase_diff = [0.0, 0.0, 0.0], s1_phase = [0.0, 0.0, 0.0], s2_phase = [0.0, 0.0, 0.0])
@test length(ged(ones(10, 10), zeros(10, 10))) == 3
@test round.(frqinst(ones(10), fs=10)) == zeros(10)

s, _ = fftdenoise(rand(10))
@test length(s) == 10

@test length(ghspectrogram(rand(100), fs=10, frq_lim=(1, 5), frq_n=10)) == 3
@test tkeo(ones(5)) == [1.0, 0.0, 0.0, 0.0, 1.0]
@test length(mwspectrogram(rand(100), fs=10, frq_lim=(1, 5), frq_n=10)) == 4
@test length(mwpsd(rand(100), fs=10, frq_lim=(1, 5), frq_n=10)) == 2
@test length(perm_cmp(ones(10,10,10), zeros(10,10,10))) == 2
@test length(fcoherence(ones(2, 10), fs=1)) == 3
@test length(fcoherence(ones(10), ones(10), fs=1)) == 3
@test l1(ones(10,10,10), zeros(10,10,10)) == 1000.0
@test l2(ones(10,10,10), zeros(10,10,10)) == 31.622776601683793
@test cums(zeros(2,2,2)) == zeros(2,2,2)
@test gfp(ones(10)) == 1.0
@test gfp_norm(ones(10)) == ones(10)
@test diss(ones(10), ones(10)) == (glob_diss = 0.0, c = 1.0)
@test length(generate_morlet_fwhm(10, 10)) == 21
@test f_nearest([(1.0, 1.0) (0.0, 0.0); (0.0, 0.0) (0.0, 0.0)], (1.0, 0.0)) == (1, 1)

mb, mf, mxb = band_mpower(ones(100), f=(1,2), fs=10)
@test round(mb) == 0.0
@test mf == 1.0
@test round(mxb) == 0.0

p, f = psd_rel(ones(10), fs=10)
@test f == [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0]

@test wbp(ones(10), fs=10, frq=3) == zeros(10)
@test round.(normalize_gauss([1, 0, 1]), digits=2) == [-0.26, -0.26, 0.55]
@test cbp(ones(100), fs=10, frq=1) == zeros(100)

sp, sf, st = spectrogram(rand(2560), fs=256)
segp, segs, tidx, fidx = specseg(sp, st, sf, t=(0.5,2), f=(10,20))
@test size(segp) == (101, 6)

@test size(denoise_wien(ones(2, 4, 1))) == (2, 4, 1)

p, _, _ = cps(zeros(100), ones(100), fs=10)
@test p == zeros(65)

@test phdiff(ones(10), zeros(10)) == zeros(10)
@test round.(normalize_log10([1, 2, 3]), digits=2) == [0.48, 0.6, 0.7]
@test round.(normalize_neglog([1, 2, 3]), digits=2) == [0.0, -0.69, -1.1]
@test round.(normalize_neglog10([1, 2, 3]), digits=2) == [0.0, -0.3, -0.48]
@test normalize_neg([1, 2, 3]) == [-2, -1, 0]
@test normalize_pos([1, 2, 3]) == [2, 3, 4]
@test normalize_perc([1, 2, 3]) == [0.0, 0.5, 1.0]
@test normalize_n([1, 2, 3], 2) == [0.0, 1.0, 2.0]
@test normalize([1, 2, 3], method=:zscore) == normalize_zscore([1, 2, 3])
@test phases(ones(ComplexF64, 10)) == zeros(10)
@test length(cwtspectrogram(rand(100), wt=wavelet(Morlet(2π), β=2), fs=10, frq_lim=(0, 5))) == 2
@test size(dw_trans(rand(100), type=:sdwt, wt=wavelet(WT.haar))) == (3, 100)
@test length(idw_trans(dw_trans(rand(100), type=:sdwt, wt=wavelet(WT.haar)), type=:sdwt, wt=wavelet(WT.haar))) == 100
@test round.(normalize_invroot([1, 2, 3]), digits=2) == [0.71, 0.58, 0.5]
@test size(cw_trans(rand(100), wt=wavelet(Morlet(2π), β=2))) == (14, 100)
@test length(icw_trans(cw_trans(rand(100), wt=wavelet(Morlet(2π), β=2)), wt=wavelet(Morlet(2π), β=2), type=:pd)) == 100
@test t2s(1, 256) == 256
@test s2t(256, 256) == 1.0
@test length(generate_noise(256, 10)) == 256

true