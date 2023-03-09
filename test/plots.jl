using NeuroAnalyzer
using Plots
using GLMakie
using Test
 
eeg = import_edf("eeg-test-edf.edf")
load_locs!(eeg, file_name="../locs/standard-10-20-cap19-elmiko.ced")
isfile("test.png") && rm("test.png")
e10 = epoch(eeg, ep_n=10)

p = plot_filter_response(fs=sr(eeg), fprototype=:butterworth, ftype=:hp, cutoff=10, order=8)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
plot_save(p, file_name="test.png")
@test isfile("test.png") == true
isfile("test.png") && rm("test.png")

p = plot(e10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot(e10, type=:mean)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot(e10, type=:butterfly)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = plot_psd(e10, norm=true, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_psd(e10, norm=true, epoch=1, channel=1, method=:mw)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_psd(e10, norm=true, epoch=1, channel=1, method=:mt)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_psd(e10, norm=true, epoch=1, channel=1, ref=:delta)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_psd(e10, norm=true, epoch=1, channel=1:10, type=:mean)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_psd(e10, norm=true, epoch=1, channel=1:10, type=:butterfly)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_psd(e10, norm=true, epoch=1, channel=1:10, type=:w3d)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_psd(e10, norm=true, epoch=1, channel=1:10, type=:topo)
@test typeof(p) == Makie.Figure

p = plot_spectrogram(e10, norm=true, epoch=1, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_spectrogram(e10, norm=true, epoch=1, channel=1, method=:stft)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_spectrogram(e10, norm=true, epoch=1, channel=1, method=:mt)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_spectrogram(e10, norm=true, epoch=1, channel=1, method=:mw)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_spectrogram(e10, norm=true, epoch=1, channel=1:10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = plot_electrodes(e10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_electrodes(e10, selected=1:4)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_electrodes(e10, threed=true)
@test typeof(p) == Makie.Figure

c = cor(e10)
channels = get_channel_bytype(e10, type=Symbol(e10.header[:signal_type]))
p = plot_matrix(c[:, :, 1], xlabels=labels(e10)[channels], ylabels=labels(e10)[channels])
@test typeof(p) == Plots.Plot{Plots.GRBackend}
ac, lags = acov(eeg, lag=5, norm=false)
p = plot_covmatrix(ac[1, :, 1], lags)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
cc, lags = xcov(eeg, lag=5, norm=false)
p = plot_covmatrix(cc[1, :, 1], lags)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = plot_weights(e10, weights=rand(19), channel=1:19)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_connections(e10, connections=rand(19, 19), channel=1:19, threshold=0.5)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = plot_topo(e10, segment=(1, 2))
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_topo(e10, segment=(1, 2), amethod=:median)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_topo(e10, segment=(1, 2), amethod=:median, imethod=:nn)

stats = rand(2, 10)
p = plot_histogram(stats[1, :])
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_bar(stats[:, 1], labels=["1", "2"])
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_line(stats[:, 1], labels=["1", "2"])
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_dots([stats[1, :], stats[2, :]], labels=["1", "2"])
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_paired([stats[1, :], stats[2, :]], labels=["1", "2"])
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_polar(stats')
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p1 = plot(e10, epoch=1, xlabel="")
p2 = plot_empty()
pp = [p1, p2]
l = (2, 1)
p = plot_compose(pp, layout=l)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

p = plot_erp(e10, channel=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_erp(e10, channel=1, type=:mean)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_erp(e10, channel=1, type=:butterfly)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_erp(e10, channel=1, type=:stack)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

c = e10.signals .^ 2
p = plot_erp(e10, c, c_idx=1)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_erp(e10, c, c_idx=1, type=:mean)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_erp(e10, c, c_idx=1, type=:butterfly)
@test typeof(p) == Plots.Plot{Plots.GRBackend}
p = plot_erp(e10, c, c_idx=1, type=:stack)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

#####

#=

e10=epoch(eeg, epoch_len=10*256)
ic, icm = ica(e10, n=16, tol=0.99)
add_component!(e10, c=:ica, v=ic)
add_component!(e10, c=:ica_mw, v=icm)
p = plot_ica_topo(e10, epoch=1, offset=0, len=10, ic=1:10)
@test typeof(p) == Plots.Plot{Plots.GRBackend}

=#

true