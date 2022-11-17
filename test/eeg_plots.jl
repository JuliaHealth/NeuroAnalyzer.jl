using NeuroAnalyzer
using Plots
using GLMakie
using Test

edf = eeg_import_edf("eeg-test-edf.edf")
eeg_load_electrodes!(edf, file_name="../locs/standard-10-20-cap19-elmiko-correct.ced")
isfile("test.png") && rm("test.png")
e10 = eeg_epochs(edf, epoch_n=10)

p = plot_filter_response(fs=eeg_sr(edf), fprototype=:butterworth, ftype=:hp, cutoff=10, order=8)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
eeg_plot_save(p, file_name="test.png")
@test isfile("test.png") == true
isfile("test.png") && rm("test.png")

p = eeg_plot(e10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot(e10, type=:mean)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot(e10, type=:butterfly)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_psd(e10, norm=true, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_psd(e10, norm=true, epoch=1, channel=1, method=:mw)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_psd(e10, norm=true, epoch=1, channel=1, method=:mt)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_psd(e10, norm=true, epoch=1, channel=1, ref=:delta)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_psd(e10, norm=true, epoch=1, channel=1:10, type=:mean)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_psd(e10, norm=true, epoch=1, channel=1:10, type=:butterfly)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_psd(e10, norm=true, epoch=1, channel=1:10, type=:w3d)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_psd(e10, norm=true, epoch=1, channel=1:10, type=:topo)
@test typeof(p) == Makie.Figure

p = eeg_plot_spectrogram(e10, norm=true, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_spectrogram(e10, norm=true, epoch=1, channel=1, method=:stft)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_spectrogram(e10, norm=true, epoch=1, channel=1, method=:mt)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_spectrogram(e10, norm=true, epoch=1, channel=1, method=:mw)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_spectrogram(e10, norm=true, epoch=1, channel=1:10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_electrodes(e10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_electrodes(e10, selected=1:4)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_electrodes(e10, threed=true)
@test typeof(p) == Makie.Figure

edf_cor = eeg_cor(e10)
channels = eeg_channel_idx(e10, type=Symbol(e10.eeg_header[:signal_type]))
p = plot_matrix(edf_cor[:, :, 1], xlabels=eeg_labels(e10)[channels], ylabels=eeg_labels(e10)[channels])
@test typeof(p) == Plots.Plot{Plots.GRBackend}
ac, lags = eeg_acov(edf, lag=5, norm=false)
p = plot_covmatrix(ac[1, :, 1], lags)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
cc, lags = eeg_xcov(edf, lag=5, norm=false)
p = plot_covmatrix(cc[1, :, 1], lags)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_weights(e10, weights=rand(19), channel=1:19)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_connections(e10, connections=rand(19, 19), channel=1:19, threshold=0.5)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = eeg_plot_topo(e10, segment=(1, 2))
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_topo(e10, segment=(1, 2), amethod=:median)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_topo(e10, segment=(1, 2), amethod=:median, imethod=:nn)

stats = rand(19, 10)
p = eeg_plot_stats(e10, stats, plot_by=:channels, channel=1:10, epoch=1, type=:hist)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_stats(e10, stats, plot_by=:channels, channel=1:10, epoch=1, type=:kd)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_stats(e10, stats, plot_by=:channels, channel=1:10, epoch=1, type=:bar)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = eeg_plot_stats(e10, stats, plot_by=:epochs, channel=1:10, epoch=1:10, type=:line)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p1 = eeg_plot(e10, epoch=1, xlabel="")
p2 = plot_empty()
pp = [p1, p2]
l = (2, 1)
p = eeg_plot_compose(pp, layout=l)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

#####

#=

e10=eeg_epochs(edf, epoch_len=10*256)
ic, icm = eeg_ica(e10, n=16, tol=0.99)
eeg_add_component!(e10, c=:ica, v=ic)
eeg_add_component!(e10, c=:ica_mw, v=icm)
p = eeg_plot_ica_topo(e10, epoch=1, offset=0, len=10, ic=1:10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

=#

true