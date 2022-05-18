using NeuroJ
using Plots
using Test

edf = eeg_import_edf("eeg-test-edf.edf")
ecg = eeg_extract_channel(edf, channel=24)
eeg_delete_channel!(edf, channel=24)
eog2 = eeg_extract_channel(edf, channel=23)
eeg_delete_channel!(edf, channel=23)
eog1 = eeg_extract_channel(edf, channel=22)
eeg_delete_channel!(edf, channel=22)
a2 = eeg_extract_channel(edf, channel=18)
eeg_delete_channel!(edf, channel=18)
a1 = eeg_extract_channel(edf, channel=17)
eeg_delete_channel!(edf, channel=17)
eeg_load_electrodes!(edf, file_name="../locs/standard-10-20-cap19-elmiko.ced")
isfile("test.png") && rm("test.png")
e10 = eeg_epochs(edf, epoch_n=10)

p = eeg_plot_filter_response(edf, fprototype=:butterworth, ftype=:hp, cutoff=10, order=8)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
eeg_plot_save(p, file_name="test.png")
@test isfile("test.png") == true
isfile("test.png") && rm("test.png")

p = eeg_plot_electrodes(edf, head=true)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_signal(edf)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_signal_avg(edf)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_signal_butterfly(edf, head=true)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_signal_details(edf, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_signal_avg_details(edf)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_signal_butterfly_details(edf)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_signal_psd(edf, norm=true, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_signal_psd_avg(edf, norm=true)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_signal_psd_butterfly(edf, norm=true)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_signal_spectrogram(edf, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_signal_spectrogram(edf, channel=1:10, len=1024)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_signal_spectrogram_avg(edf, channel=1:10)
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

p = eeg_plot_histogram(edf, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_bands(edf, channel=1, type=:abs)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_signal_topo(edf, offset=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
e10 = eeg_epochs(edf, epoch_len=10*256)
_, a, pow, ph = eeg_spectrum(e10)
p = eeg_plot_acomponent_topo(e10, c=ph, epoch=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_weights_topo(edf, weights=rand(19), epoch=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_mcomponent_topo(edf, c=rand(19, 1), epoch=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p1 = eeg_plot_signal(e10, epoch=1)
p2 = eeg_plot_signal(e10, epoch=2)
pp = [p1, p2]
l = (2, 1)
p = eeg_plot_compose(pp, layout=l)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

snr = eeg_snr(edf)
p = eeg_plot_channels(edf, c=snr, epoch=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

e = eeg_epochs_stats(e10)
p = eeg_plot_epochs(e10, c=e[4])
@test typeof(p) == Plots.Plot{Plots.GRBackend}
eeg_add_component!(e10, c=:epochs_var, v=e[4])
p = eeg_plot_epochs(e10, c=:epochs_var)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

_, a, pow, ph = eeg_spectrum(e10)
p = eeg_plot_component(e10, c=pow, epoch=10, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_avg(e10, c=pow, epoch=10, channel=1:4)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_butterfly(e10, c=pow, epoch=10, channel=1:4)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_component_idx(e10, c=pow, epoch=10, c_idx=1:5)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_idx_avg(e10, c=pow, epoch=10, c_idx=1:5)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_idx_butterfly(e10, c=pow, epoch=10, c_idx=1:5)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

s = eeg_tconv(e10, kernel=generate_window(:hann, 128))
p = eeg_plot_component_psd(e10, c=s, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_psd_avg(e10, c=s, epoch=1, channel=1:10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_psd_butterfly(e10, c=s, epoch=1, channel=1:10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_component_idx_psd(e10, c=pow, epoch=10, c_idx=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_idx_psd_avg(e10, c=pow, epoch=10, c_idx=1:5)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_idx_psd_butterfly(e10, c=pow, epoch=10, c_idx=1:5)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_component_idx_spectrogram(e10, c=pow, epoch=10, c_idx=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_idx_spectrogram(e10, c=pow, epoch=10, c_idx=1:5)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_idx_spectrogram_avg(e10, c=pow, epoch=10, c_idx=1:5)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_component_spectrogram(e10, c=s, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_spectrogram(e10, c=s, epoch=1, channel=1:10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_component_spectrogram_avg(e10, c=s, epoch=1, channel=1:10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

e10=eeg_epochs(edf, epoch_len=10*256)
ic, icm = eeg_ica(e10, n=16, tol=0.99)
eeg_add_component!(e10, c=:ica, v=ic)
eeg_add_component!(e10, c=:ica_mw, v=icm)
p = eeg_plot_ica_topo(e10, epoch=1, offset=0, len=10, ic=1:10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_env(e10, type=:amp, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_env(e10, type=:amp, average=:mean, dims=3, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_env(e10, type=:pow, average=:median, dims=2, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_env(e10, type=:spec, average=:mean, dims=1, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_ispc(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=2)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_itpc(e10, channel=1, t=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=2)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_signal_psd_3d(edf, channel=1:5, offset=25600, mt=true)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

true