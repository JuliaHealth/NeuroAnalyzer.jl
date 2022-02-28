using NeuroJ
using Test

fs = 10
t = collect(0:1/fs:10)
signal_v = generate_sine(1, t)
signal_v1 = generate_sine(2, t)
signal_v2 = generate_sine(4, t)
signal_a = Array{Float64}(undef, 2, length(signal_v), 2)
signal_m = Matrix{Float64}(undef, 2, length(signal_v))
signal_m[1, :] = signal_v1
signal_m[2, :] = signal_v2
signal_a[1, :, 1] = signal_v1
signal_a[2, :, 1] = signal_v2
signal_a[1, :, 2] = signal_v1
signal_a[2, :, 2] = signal_v2
signal_m1 = signal_m .* 0.1
signal_m2 = signal_m .* 0.2
signal_a1 = signal_a .* 0.1
signal_a2 = signal_a .* 0.2

@test size(signal_derivative(signal_v)) == (101, )

@test typeof(signal_total_power(signal_v, fs=fs)) == Float64
@test size(signal_total_power(signal_a, fs=fs)) == (2, 2)

@test typeof(signal_band_power(signal_v, fs=fs, f=(2, 4))) == Float64
@test size(signal_band_power(signal_a, fs=fs, f=(2, 4))) == (2, 2)

signal_fft, signal_sf = signal_make_spectrum(signal_v, fs=fs)
@test size(signal_sf) == (101, )

signal_fft, signal_sf = signal_make_spectrum(signal_a, fs=fs)
@test size(abs.(signal_fft)) == (2, 101, 2)
@test size(signal_sf) == (2, 101, 2)

@test typeof(signal_detrend(signal_v, type=:linear)) == Vector{Float64}
@test size(signal_detrend(signal_v, type=:constant)) == (101, )
@test typeof(signal_detrend(signal_a, type=:linear)) == Array{Float64, 3}
@test size(signal_detrend(signal_a, type=:constant)) == (2, 101, 2)

m, s, u, l = signal_ci95(signal_v)
@test typeof(m) == Float64
@test typeof(s) == Float64
@test typeof(u) == Float64
@test typeof(l) == Float64

m, s, u, l = signal_ci95(signal_m)
@test typeof(m) == Vector{Float64}
@test typeof(s) == Vector{Float64}
@test typeof(u) == Vector{Float64}
@test typeof(l) == Vector{Float64}

m, s, u, l = signal_ci95(signal_a)
@test size(m) == (2, 101)
@test size(s) == (2, 101)
@test size(u) == (2, 101)
@test size(l) == (2, 101)

@test typeof(signal_mean(signal_v1, signal_v2)) == NTuple{4, Float64}
m, s, u, l = signal_mean(signal_a1, signal_a2)
@test typeof(m) == Matrix{Float64}
@test typeof(s) == Matrix{Float64}
@test typeof(u) == Matrix{Float64}
@test typeof(l) == Matrix{Float64}

signals_statistic, signals_statistic_single, p = signal_difference(signal_m1, signal_m2)
@test p == 1.0
signals_statistic, signals_statistic_single, p = signal_difference(signal_a1, signal_a2)
@test size(signals_statistic) == (2, 6)
@test size(signals_statistic_single) == (2, )

acov, lags = signal_autocov(signal_v, lag=10)
@test size(acov) == (21, )
acov, lags = signal_autocov(signal_a)
@test size(acov) == (2, 3, 2)

ccov, lags = signal_crosscov(signal_v1, signal_v2, lag=10)
@test size(ccov) == (21, )
ccov, lags = signal_crosscov(signal_a)
@test size(ccov) == (4, 3, 2)
ccov, lags = signal_crosscov(signal_a1, signal_a2)
@test size(ccov) == (2, 3, 2)

signal_fft, signal_amplitudes, signal_powers, signal_phases = signal_spectrum(signal_v, pad=9)
@test size(signal_fft) == (110, )
@test size(signal_amplitudes) == (110, )
@test size(signal_powers) == (110, )
@test size(signal_phases) == (110, )

signal_fft, signal_amplitudes, signal_powers, signal_phases = signal_spectrum(signal_a, pad=9)
@test size(signal_fft) == (2, 110, 2)
@test size(signal_amplitudes) == (2, 110, 2)
@test size(signal_powers) == (2, 110, 2)
@test size(signal_phases) == (2, 110, 2)

epochs = signal_epochs(signal_v, epoch_n=10)
@test size(epochs) == (10, 10)
epochs = signal_epochs(signal_m, epoch_n=10, average=true)
@test size(epochs) == (2, 10)

signal = signal_delete_channel(signal_m, channel=1)
@test size(signal) == (1, 101)

signal = signal_delete_channel(signal_a, channel=1)
@test size(signal) == (1, 101, 2)

signal_ref = signal_reference_channel(signal_m, channel=1)
@test size(signal_ref) == (2, 101)

signal_ref = signal_reference_channel(signal_a, channel=1)
@test size(signal_ref) == (2, 101, 2)

signal_ref = signal_reference_car(signal_m)
@test size(signal_ref) == (2, 101)

signal_ref = signal_reference_car(signal_a)
@test size(signal_ref) == (2, 101, 2)

signal_tap = signal_taper(signal_v, taper=signal_v)
@test size(signal_tap) == (101, )

signal_tap = signal_taper(signal_a, taper=signal_v)
@test size(signal_tap) == (2, 101, 2)

signal_dem = signal_demean(signal_v)
@test size(signal_dem) == (101, )

signal_dem = signal_demean(signal_a)
@test size(signal_dem) == (2, 101, 2)

signal_norm = signal_normalize_zscore(signal_v)
@test size(signal_norm) == (101, )

signal_norm = signal_normalize_zscore(signal_a)
@test size(signal_norm) == (2, 101, 2)

signal_norm = signal_normalize_minmax(signal_v)
@test size(signal_norm) == (101, )

signal_norm = signal_normalize_minmax(signal_a)
@test size(signal_norm) == (2, 101, 2)

cov_mat = signal_cov(signal_a)
@test size(cov_mat) == (2, 2, 2)

cor_mat = signal_cor(signal_a)
@test size(cov_mat) == (2, 2, 2)

sn = signal_add_noise(signal_v)
@test length(sn) == 101
sn = signal_add_noise(signal_a)
@test size(sn) == (2, 101, 2)

s_up, t_up = signal_upsample(signal_v, t=0:1/fs:10, new_sr=15)
@test size(s_up) == (151,)
s_up, t_up = signal_upsample(signal_a, t=0:1/fs:10, new_sr=15)
@test size(s_up) == (2, 151, 2)

s_tcov = signal_tconv(signal_v, kernel=[1, 0, 1])
@test size(s_tcov) == (101, )
s_tcov = signal_tconv(signal_a, kernel=[1, 0, 1])
@test size(s_tcov) == (2, 101, 2)

s_filt = signal_filter(signal_v, fprototype=:butterworth, ftype=:lp, cutoff=2, fs=fs, order=8)
@test size(s_filt) == (101, )
s_filt = signal_filter(signal_a, fprototype=:butterworth, ftype=:lp, cutoff=2, fs=fs, order=8)
@test size(s_filt) == (2, 101, 2)
s_filt = signal_filter(signal_v, fprototype=:mavg, d=10)
@test size(s_filt) == (101,)
s_filt = signal_filter(signal_add_noise(signal_v), fprototype=:mavg, window=generate_gaussian(fs, 32, 0.01))
@test size(s_filt) == (101, )
s_filt = signal_filter(signal_a, fprototype=:mmed, t=2.2)
@test size(s_filt) == (2, 101, 2)
s_filt = signal_filter(signal_v, fprototype=:poly, fs=fs, order=40)
@test size(s_filt) == (101, )

s_down, t_down = signal_downsample(signal_v, t=0:1/fs:10, new_sr=2)
@test typeof(s_down) == Vector{Float64}
s_down, t_down = signal_downsample(signal_a, t=0:1/fs:10, new_sr=2)
@test typeof(s_down) == Array{Float64, 3}

p, f = signal_psd(signal_v, fs=fs)
@test length(p) == 21
@test size(f) == (21, )

p, f = signal_psd(signal_m, fs=fs)
@test length(p) == 42
@test size(f) == (2, 21)

p, f = signal_psd(signal_a, fs=fs)
@test length(p) == 84
@test size(f) == (2, 21, 2)

s = signal_stationarity_hilbert(signal_v)
@test size(s) == (100, )
s = signal_stationarity_mean(signal_v, window=10)
@test size(s) == (1, 10)
s = signal_stationarity_var(signal_v, window=10)
@test size(s) == (1, 10)

p = signal_stationarity(signal_a, method=:hilbert)
@test size(p) == (2, 100, 2)
p = signal_stationarity(signal_a, window=100, method=:euclid)
@test size(p) == (2, 2)
p = signal_stationarity(signal_a, window=5, method=:mean)
@test size(p) == (2, 5, 2)
p = signal_stationarity(signal_a, window=10, method=:var)
@test size(p) == (2, 10, 2)

s = signal_trim(signal_v, len=10, offset=80, from=:start)
@test length(s) == 91
s = signal_trim(signal_v, len=11, from=:end)
@test length(s) == 90
s = signal_trim(signal_a, len=10, offset=30, from=:start)
@test size(s) == (2, 91, 2) 
s = signal_trim(signal_a, len=11, from=:end)
@test size(s) == (2, 90, 2)

m = signal_mi(signal_v1, signal_v2)
@test m == 2.41069306603075
m = signal_mi(signal_a)
@test size(m) == (2, 2, 2)
m = signal_mi(signal_a1, signal_a2)
@test size(m) == (2, 2, 2)

e = signal_entropy(signal_v)
@test typeof(e) == Float64
e = signal_entropy(signal_a)
@test size(e) == (2, 2)

s = signal_average(signal_v1, signal_v2)
@test size(s) == (101, 1)
s = signal_average(signal_a1, signal_a2)
@test size(s) == (2, 101, 2)

s = signal_coherence(signal_v1, signal_v2)
@test size(s) == (101, )
s = signal_coherence(signal_a1, signal_a2)
@test size(s) == (2, 101, 2)

p, v = signal_pca(signal_a, n=2)
@test size(p) == (2, 101, 2)
@test size(v) == (2, 2)

s = signal_fconv(signal_v1, kernel=[1, 2, 3])
@test size(s) == (101, )
s = signal_fconv(signal_a1, kernel=[1, 2, 3])
@test size(s) == (2, 101, 2)

m, md, s, v, k = signal_epochs_stats(signal_a)
@test size(s) == (2, )

p, f, t = signal_spectrogram(signal_v, fs=fs)
@test size(p) == (51, 46)
p, f, t = signal_spectrogram(signal_a, fs=fs)
@test size(p) == (51, 46, 2, 2)

i, mw = signal_ica(signal_a, n=1, tol=1.0)
@test size(i) == (1, 101, 2)
s = signal_ica_reconstruct(signal_a, ic_activations=i, ic_mw=mw, ic_v=1)
@test size(s) == (2, 101, 2)

c, b = signal_detect_flat(signal_a, len=2)
@test typeof(c) == Set{Int64}

true