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

@test floor(sum(signal_derivative(signal_v)), digits=2) == 0.58

@test floor(signal_total_power(signal_v, fs=fs), digits=2) == 0.33
@test floor(sum(signal_total_power(signal_a, fs=fs)), digits=2) == 1.33

@test floor(signal_band_power(signal_v, fs=fs, f1=2, f2=4), digits=2) == 0.0
@test floor(sum(signal_band_power(signal_a, fs=fs, f1=2, f2=4)), digits=2) == 0.66

signal_fft, signal_sf = signal_make_spectrum(signal_v, fs)
@test floor(sum(signal_sf), digits=2) == 0.0

signal_fft, signal_sf = signal_make_spectrum(signal_a, fs)
@test floor(sum(abs.(signal_fft)), digits=2) == 1072.35
@test floor(sum(signal_sf), digits=2) == -0.01

@test floor(sum(signal_detrend(signal_v, type=:linear)), digits=2) == 0.0
@test floor(sum(signal_detrend(signal_v, type=:constant)), digits=2) == 0.0
@test floor(sum(signal_detrend(signal_a, type=:linear)), digits=2) == 0.0
@test floor(sum(signal_detrend(signal_a, type=:constant)), digits=2) == 0.0

m, s, u, l = signal_ci95(signal_v)
@test floor(m, digits=2) == 0.0
@test floor(s, digits=2) == 0.07
@test floor(u, digits=2) == 0.13
@test floor(l, digits=2) == -0.14

m, s, u, l = signal_ci95(signal_m)
@test floor(sum(m), digits=2) == -0.01
@test floor(sum(s), digits=2) == 38.04
@test floor(sum(u), digits=2) == 74.56
@test floor(sum(l), digits=2) == -74.57

m, s, u, l = signal_ci95(signal_a)
@test floor(sum(m), digits=2) == -0.01
@test floor(sum(s), digits=2) == 76.08
@test floor(sum(u), digits=2) == 149.12
@test floor(sum(l), digits=2) == -149.13

@test floor(sum(signal_mean(signal_v1, signal_v2)), digits=2) == 0.09
m, s, u, l = signal_mean(signal_a1, signal_a2)
@test floor(sum(m), digits=2) == 0.0
@test floor(sum(s), digits=2) == 2.39
@test floor(sum(u), digits=2) == 4.69
@test floor(sum(l), digits=2) == -4.7

signals_statistic, signals_statistic_single, p = signal_difference(signal_m1, signal_m2)
@test p == 1.0
signals_statistic, signals_statistic_single, p = signal_difference(signal_a1, signal_a2)
@test size(signals_statistic) == (2, 6)
@test floor(sum(signals_statistic_single), digits=2) == 0.15

acov, lags = signal_autocov(signal_v, lag=10)
@test floor(sum(acov), digits=2) == 44.99
acov, lags = signal_autocov(signal_a)
@test floor(sum(acov), digits=2) == 99.99

ccov, lags = signal_crosscov(signal_v1, signal_v2, lag=10)
@test floor(sum(ccov), digits=2) == 0.0
ccov, lags = signal_crosscov(signal_a)
@test floor(sum(ccov), digits=2) == 100.0
ccov, lags = signal_crosscov(signal_a1, signal_a2)
@test floor(sum(ccov), digits=2) == 1.99

signal_fft, signal_amplitudes, signal_powers, signal_phases = signal_spectrum(signal_v, pad=9)
@test round(sum(signal_fft), digits=2) == -0.0 - 0.0im
@test floor(sum(signal_amplitudes), digits=2) == 4.88
@test floor(sum(signal_powers), digits=2) == 2.15
@test floor(sum(signal_phases), digits=2) == 3.14

signal_fft, signal_amplitudes, signal_powers, signal_phases = signal_spectrum(signal_a, pad=9)
@test round(sum(signal_fft), digits=2) == -0.0 - 0.0im
@test floor(sum(signal_amplitudes), digits=2) == 20.38
@test floor(sum(signal_powers), digits=2) == 8.62
@test floor(sum(signal_phases), digits=2) == 18.84

epochs = signal_epochs(signal_v, epochs_no=10)
@test floor(sum(epochs), digits=2) == 0.0
epochs = signal_epochs(signal_m, epochs_no=10, average=true)
@test floor(sum(epochs), digits=2) == -0.01

signal = signal_drop_channel(signal_m, 1)
@test floor(sum(signal), digits=2) == -0.01

signal = signal_drop_channel(signal_a, 1)
@test floor(sum(signal), digits=2) == -0.01

signal_ref = signal_reference_channel(signal_m, 1)
@test floor(sum(signal_ref), digits=2) == -0.01

signal_ref = signal_reference_channel(signal_a, 1)
@test floor(sum(signal_ref), digits=2) == -0.01

signal_ref = signal_reference_car(signal_m)
@test floor(sum(signal_ref), digits=2) == -0.01

signal_ref = signal_reference_car(signal_a)
@test floor(sum(signal_ref), digits=2) == -0.01

signal_tap = signal_taper(signal_v, signal_v)
@test floor(sum(signal_tap), digits=2) == 49.99

signal_tap = signal_taper(signal_a, signal_v)
@test floor(sum(signal_tap), digits=2) == 0.0

signal_dem = signal_demean(signal_v)
@test floor(sum(signal_dem), digits=2) == 0.0

signal_dem = signal_demean(signal_a)
@test floor(sum(signal_dem), digits=2) == 0.0

signal_norm = signal_normalize_zscore(signal_v)
@test floor(sum(signal_norm), digits=2) == 0.0

signal_norm = signal_normalize_zscore(signal_a)
@test floor(sum(signal_norm), digits=2) == 0.0

signal_norm = signal_normalize_minmax(signal_v)
@test floor(sum(signal_norm), digits=2) == 50.49

signal_norm = signal_normalize_minmax(signal_a)
@test floor(sum(signal_norm), digits=2) == 202.0

cov_mat = signal_cov(signal_v1, signal_v2)
@test floor(sum(cov_mat), digits=2) == 0.0

cov_mat = signal_cov(signal_a)
@test floor(sum(cov_mat), digits=2) == 1.0

cor_mat = signal_cor(signal_a)
@test floor(sum(cov_mat), digits=2) == 1.0

sn = signal_add_noise(signal_v)
@test length(sn) == 101
sn = signal_add_noise(signal_a)
@test size(sn) == (2, 101, 2)

s_up, t_up = signal_upsample(signal_v, t=0:1/fs:10, new_sr=15)
@test floor(sum(s_up), digits=2) == 0.0
s_up, t_up = signal_upsample(signal_a, t=0:1/fs:10, new_sr=15)
@test floor(sum(s_up), digits=2) == -0.01

s_tcov = signal_tconv(signal_v, [1, 0, 1])
@test floor(sum(s_tcov), digits=2) == 0.0
s_tcov = signal_tconv(signal_a, [1, 0, 1])
@test floor(sum(s_tcov), digits=2) == -0.01

s_filt = signal_filter(signal_v, fprototype=:butterworth, ftype=:lp, cutoff=2, fs=fs, order=8)
@test floor(sum(s_filt), digits=2) == -0.01
s_filt = signal_filter(signal_a, fprototype=:butterworth, ftype=:lp, cutoff=2, fs=fs, order=8)
@test floor(sum(s_filt), digits=2) == -0.01
s_filt = signal_filter(signal_v, fprototype=:mavg, d=10)
@test length(s_filt) == 101
s_filt = signal_filter(signal_add_noise(signal_v), fprototype=:mavg, window=generate_gaussian(fs, 32, 0.01))
@test length(s_filt) == 101
s_filt = signal_filter(signal_a, fprototype=:mmed, t=2.2)
@test length(s_filt) == 404

s_down, t_down = signal_downsample(signal_v, t=0:1/fs:10, new_sr=2)
@test typeof(s_down) == Vector{Float64}
s_down, t_down = signal_downsample(signal_a, t=0:1/fs:10, new_sr=2)
@test typeof(s_down) == Array{Float64, 3}

p, f = signal_psd(signal_v, fs=fs)
@test length(p) == 21
@test floor(sum(f), digits=2) == 52.5

p, f = signal_psd(signal_m, fs=fs)
@test length(p) == 42
@test floor(sum(f), digits=2) == 105.0

p, f = signal_psd(signal_a, fs=fs)
@test length(p) == 84
@test floor(sum(f), digits=2) == 210.0

s = signal_stationarity_hilbert(signal_v)
@test floor(sum(s), digits=2) == 43.98
s = signal_stationarity_mean(signal_v, window=10)
@test floor(sum(s), digits=2) == -0.14
s = signal_stationarity_var(signal_v, window=10)
@test floor(sum(s), digits=2) == 1.12

p = signal_stationarity(signal_a, method=:hilbert)
@test floor(sum(p), digits=2) == 753.98
p = signal_stationarity(signal_a, window=100, method=:euclid)
@test floor(sum(p), digits=2) == 0.08
p = signal_stationarity(signal_a, window=5, method=:mean)
@test floor(sum(p), digits=2) == -0.01
p = signal_stationarity(signal_a, window=10, method=:var)
@test floor(sum(p), digits=2) == 22.22

s = signal_trim(signal_v, trim_len=10, offset=80, from=:start)
@test length(s) == 91 
s = signal_trim(signal_v, trim_len=11, from=:end)
@test length(s) == 90
s = signal_trim(signal_a, trim_len=10, offset=30, from=:start)
@test size(s) == (2, 91, 2) 
s = signal_trim(signal_a, trim_len=11, from=:end)
@test size(s) == (2, 90, 2)

m = signal_mi(signal_v1, signal_v2)
@test m == 2.41069306603075
m = signal_mi(signal_a)
@test floor(sum(m), digits=2) == 19.65
m = signal_mi(signal_a1, signal_a2)
@test floor(sum(m), digits=2) == 19.65

e = signal_entropy(signal_v)
@test floor(sum(e), digits=2) == 2.28
e = signal_entropy(signal_a)
@test floor(sum(e), digits=2) == 6.61

s = signal_average(signal_v1, signal_v2)
@test floor(sum(s), digits=2) == -0.01
s = signal_average(signal_a1, signal_a2)
@test floor(sum(s), digits=2) == -0.01

true