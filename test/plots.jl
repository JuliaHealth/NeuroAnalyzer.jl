using NeuroJ
using Plots
using Test

edf = eeg_import_edf("eeg-test-edf.edf")
eeg_load_electrodes!(edf, file_name="../locs/standard-10-20-cap19-elmiko.ced")
signal_v = rand(1000)
signal_m = rand(10, 1000)
t = linspace(0, 10, 1000)
isfile("test.png") && rm("test.png")

p = signal_plot(t, signal_v)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = signal_plot(t, signal_m)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot(edf)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
eeg_plot_save(p, file_name="test.png")
@test isfile("test.png") == true
isfile("test.png") && rm("test.png")

p = eeg_plot_filter_response(edf, fprototype=:butterworth, ftype=:hp, cutoff=10, order=8)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = signal_plot_avg(t, signal_m)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_avg(edf)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = signal_plot_butterfly(t, signal_m)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_butterfly(edf, head=true)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

signal_pow, signal_frq = signal_psd(signal_v, fs=100, norm=true)
p = signal_plot_psd(signal_pow, signal_frq)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = signal_plot_psd(signal_v, fs=100)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = signal_plot_psd(signal_m, fs=100)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_psd(edf, norm=true, average=true, head=true)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_electrodes(edf, head=true)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

edf_cor = eeg_cor(edf)
p = eeg_plot_matrix(edf, edf_cor)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

ac, lags = eeg_autocov(edf, lag=5, norm=false)
p = eeg_plot_covmatrix(edf, ac, lags)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

cc, lags = eeg_crosscov(edf, lag=5, norm=false)
p = eeg_plot_covmatrix(edf, cc, lags)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_spectrogram(edf, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_spectrogram(edf, channel=1:10, len=1024)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = signal_plot_histogram(signal_v)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = signal_plot_histogram(signal_m, type=:kd)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_histogram(edf, channel=1)

p = signal_plot_band(signal_v, fs=100, type=:abs)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_band(edf, channel=1, type=:abs)

p = eeg_plot_topo(edf, offset=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

true