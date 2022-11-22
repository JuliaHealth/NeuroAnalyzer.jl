using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@test linspace(1, 10, 10) == 1.0:10.0
@test logspace(0, 1, 3) == [1.0, 3.1622776601683795, 10.0]
@test length(pad0(ones(3), 1)) == 4
@test size(pad0(ones(3, 3), 1)) == (3, 4)
@test length(pad2(ones(10))) == 16
@test size(pad2(ones(3, 3))) == (3, 4)
@test size(m_pad0(ones(3, 4))) == (4, 4)
@test vsearch([1, 2], [1, 2, 3, 4]) == [1, 2]
@test cart2pol(1, 1) == (1.41, 45.0)
@test pol2cart(1, 1) == (1.0, 0.02)
@test sph2cart(1, 45, 45) == (0.5, 0.5, 0.71)
@test cart2sph(1, 1, 1) == (1.73, 45.0, 54.69)
@test generate_window(:hann, 3) == [0.0, 1.0, 0.0]
@test length(fft0(ones(10), 10)) == 20
@test length(ifft0(ones(10), 10)) == 20
@test nextpow2(5) == 8
@test vsplit(ones(4), 2) == [[1.0, 1.0], [1.0, 1.0]]
@test length(s_rms(ones(10))) == 1
@test length(generate_sine(2, collect(1:10))) == 10
@test length(s_freqs(rand(10))) == 2
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
@test s2_rmse(ones(10), ones(10)) == 0.0
@test size(m_norm(ones(4, 4, 1))) == (4, 4, 1)
@test size(s_cov(zeros(2))) == (2, 2)
@test s2_cov(ones(2), zeros(2)) == zeros(2, 2)
@test s_dft(ones(4), fs=10) == (s_fft = ComplexF64[4.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im], s_sf = [0.0, 0.025, -0.05, -0.025])
@test s_msci95(ones(4)) == (1.0, 0.0, 1.0, 1.0)
@test s2_mean(ones(4), zeros(4)) == (1.0, 0.0, 1.0, 1.0)
@test length(s2_difference(ones(4), zeros(4))) == 3
@test s_acov(ones(4)) == ([3.0, 4.0, 3.0], [-1, 0, 1])
@test s2_xcov(ones(4), ones(4)) == ([3.0, 4.0, 3.0], [-1, 0, 1])
@test s_spectrum(ones(4)) == (s_fft = ComplexF64[4.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im], s_amp = [1.0, 0.0], s_pow = [1.0, 0.0], s_pha = [0.0, 0.0, 0.0, 0.0])
@test s_total_power(ones(4), fs=10) == 0.0
@test s_band_power(ones(4), fs=10, f=(1,2)) == 0.0
@test s_taper(ones(10), taper=zeros(10)) == zeros(10)
@test s_detrend(ones(10)) == zeros(10)
@test s_demean(ones(10)) == zeros(10)
@test s_normalize_zscore([1, 2, 3]) == [-1.0, 0.0, 1.0]
@test s_normalize_minmax([1, 2, 3]) == [-1.0, 0.0, 1.0]
@test s_normalize_max([1, 2, 1]) == [0.5, 1.0, 0.5]
@test s_normalize_log([0, 0, 0]) == [0.0, 0.0, 0.0]
@test length(s_add_noise(ones(10))) == 10

s, t = s_resample(ones(10), t=1:10, new_sr=20)
@test t == 1.0:0.05:10.0

@test s_derivative(ones(10)) == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test s_tconv(ones(10), kernel=[1.0, 1.0]) == [1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
@test s_filter(ones(10), fs=1, fprototype=:mavg) == ones(10)

p, f = s_psd(ones(10), fs=10)
@test f == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

@test s_stationarity_hilbert(ones(10)) == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test s_stationarity_mean(ones(10), window=1) == [1.0;;]
@test s_stationarity_var(ones(10), window=1) == [0.0;;]
@test s_trim(ones(10), segment=(1,5)) == ones(5)
@test s2_mi(ones(10), ones(10)) == 0.0
@test s_entropy([1.0, 2.0, 3.0]) == (ent = 1.5849625007211552, sent = 0.8304717124362917, leent = 4.333653050389665)
@test s_negentropy([1, 2, 3]) == -0.16602396751648252
@test s_average(ones(10, 10, 1)) == ones(1, 10, 1)
@test s2_average(ones(5, 5, 1), zeros(5, 5, 1)) == [0.5; 0.5; 0.5; 0.5; 0.5;;;]
@test s2_tcoherence([1, 2], [3, 4]) == (c = [0.04166666666666667, 0.027777777777777776], msc = [0.006944444444444444, 0.0007716049382716049], ic = [0.07216878364870322, -0.0])

p, w, m, pca = s_pca(ones(2, 10, 1), n=1)
@test p == [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;;;]

s = s_pca_reconstruct(ones(2, 10, 1), pc=p, pca=pca)
@test s == ones(2, 10, 1)

@test round.(s_fconv(ones(10), kernel=[1.0, 2.0])) == [1.00 + 0.0im, 3.0 + 0.0im, 3.0 + 0.0im, 3.0 + 0.0im, 3.0 + 0.0im, 3.0 + 0.0im, 3.0 + 0.0im, 3.0 + 0.0im, 3.0 + 0.0im, 3.0 + 0.0im]

i, m = s_ica([1.0 2.0; 3.0 4.0;;;], n=1)
@test size(i) == (1, 2, 1)
@test s_ica_reconstruct([1.0 2.0; 3.0 4.0;;;], ica=i, ica_mw=m, ic=[1]) == zeros(2, 2, 1)

p, f, t = s_spectrogram(ones(100), fs=10)
@test size(p) == (51, 46)

@test s_detect_channel_flat(ones(100)) == 0.9583333333333334
@test s_snr(ones(10)) == Inf
@test s_findpeaks(repeat([0, 1], 100)) == [6, 38, 70, 102, 134, 166, 198]
@test length(s_wdenoise(rand(100), wt=wavelet(WT.haar))) == 100
@test s2_ispc([1.0, 1.0, 1.0], [0.0, 0.0, 0.0]) == (ispc = 1.0, ispc_angle = 0.0, signal_diff = [-1.0, -1.0, -1.0], phase_diff = [0.0, 0.0, 0.0], s1_phase = [0.0, 0.0, 0.0], s2_phase = [0.0, 0.0, 0.0])
@test s_itpc(ones(1, 10, 10), t=1) == (itpc = 1.0, itpcz = 10.0, itpc_angle = 0.0, itpc_phases = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
@test s2_pli([1.0, 1.0, 1.0], [0.0, 0.0, 0.0]) == (pli = 0.0, signal_diff = [-1.0, -1.0, -1.0], phase_diff = [0.0, 0.0, 0.0], s1_phase = [0.0, 0.0, 0.0], s2_phase = [0.0, 0.0, 0.0])
@test length(s2_ged(ones(10, 10), zeros(10, 10))) == 3
@test s_frqinst(ones(10), fs=10) == zeros(10)

s, _ = s_fftdenoise(rand(10))
@test length(s) == 10

@test length(s_ghspectrogram(rand(100), fs=10, frq_lim=(1, 5), frq_n=10)) == 3
@test s_tkeo(ones(5)) == [1.0, 0.0, 0.0, 0.0, 1.0]
@test length(s_wspectrogram(rand(100), fs=10, frq_lim=(1, 5), frq_n=10)) == 4
@test length(s_mwpsd(rand(100), fs=10, frq_lim=(1, 5), frq_n=10)) == 2
@test length(a2_cmp(ones(10,10,10), zeros(10,10,10))) == 2
@test length(s_fcoherence(ones(2, 10), fs=1)) == 3
@test length(s2_fcoherence(ones(10), ones(10), fs=1)) == 3
@test a2_l1(ones(10,10,10), zeros(10,10,10)) == 1000.0
@test a2_l2(ones(10,10,10), zeros(10,10,10)) == 31.622776601683793
@test s_cums(ones(10)) == 1:10
@test s_cums(zeros(2,2,2)) == zeros(2,2,2)
@test s_gfp(ones(10)) == 1.0
@test s_gfp_norm(ones(10)) == ones(10)
@test s2_diss(ones(10), ones(10)) == (diss = 0.0, c = 1.0)
@test length(generate_morlet_fwhm(10, 10)) == 21
@test f_nearest([(1.0, 1.0) (0.0, 0.0); (0.0, 0.0) (0.0, 0.0)], (1.0, 0.0)) == (1, 1)
@test s_band_mpower(ones(100), f=(1,2), fs=10) == (mbp = 0.0, maxfrq = 1.0, maxbp = 0.0)

p, f = s_rel_psd(ones(10), fs=10)
@test f == [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0]

@test s_wbp(ones(10), fs=10, frq=3) == zeros(10)
@test round.(s_normalize_gauss([1, 0, 1]), digits=2) == [-0.26, -0.26, 0.55]
@test s_cbp(ones(100), fs=10, frq=1) == zeros(100)

sp, sf, st = s_spectrogram(rand(2560), fs=256)
segp, segs, tidx, fidx = s_specseg(sp, st, sf, t=(0.5,2), f=(10,20))
@test size(segp) == (101, 10)

@test size(s_denoise_wien(ones(2, 4, 1))) == (2, 4, 1)

p, _, _ = s2_cps(zeros(100), ones(100), fs=10)
@test p == zeros(65)

@test s2_phdiff(ones(100), zeros(100)) == zeros(100)
@test round.(s_normalize_log10([1, 2, 3]), digits=2) == [0.48, 0.6, 0.7]
@test round.(s_normalize_neglog([1, 2, 3]), digits=2) == [0.0, -0.69, -1.1]
@test round.(s_normalize_neglog10([1, 2, 3]), digits=2) == [0.0, -0.3, -0.48]
@test s_normalize_neg([1, 2, 3]) == [-2, -1, 0]
@test s_normalize_pos([1, 2, 3]) == [2, 3, 4]
@test s_normalize_perc([1, 2, 3]) == [0.0, 0.5, 1.0]
@test s_normalize([1, 2, 3], method=:zscore) == s_normalize_zscore([1, 2, 3])
@test s_phases(ones(ComplexF64, 10)) == zeros(10)
@test length(s_cwtspectrogram(rand(100), wt=wavelet(Morlet(π), β=2), fs=10, frq_lim=(0, 5))) == 2
@test size(s_dwt(rand(100), type=:sdwt, wt=wavelet(WT.haar))) == (3, 100)
@test length(s_idwt(s_dwt(rand(100), type=:sdwt, wt=wavelet(WT.haar)), type=:sdwt, wt=wavelet(WT.haar))) == 100
@test round.(s_normalize_invroot([1, 2, 3]), digits=2) == [0.71, 0.58, 0.5]
@test size(s_cwt(rand(100), wt=wavelet(Morlet(π), β=2))) == (14, 100)
@test length(s_icwt(s_cwt(rand(100), wt=wavelet(Morlet(π), β=2)), wt=wavelet(Morlet(π), β=2), type=:pd)) == 100
@test t2s(1, 256) == 256
@test s2t(256, 256) == 1.0

true