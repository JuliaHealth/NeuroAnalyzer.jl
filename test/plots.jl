using NeuroAnalyzer
using Plots
using GLMakie
using Test
 
@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
load_locs!(e10, file_name=joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))

isfile("test.png") && rm("test.png")

@info "test 1/21: plot_compose()"
p1 = NeuroAnalyzer.plot(e10, ep=1, xlabel="")
p2 = plot_empty()
pp = [p1, p2]
l = (2, 1)
p = NeuroAnalyzer.plot_compose(pp, layout=l)
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 2/21: plot_connections()"
p = NeuroAnalyzer.plot_connections(e10, connections=rand(19, 19), ch=1:19, threshold=0.5)
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 3/21: plot_erp()"
e10_erp = erp(e10)
p = NeuroAnalyzer.plot_erp(e10_erp, ch=1)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=1:10, type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=1, type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=1:10, type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=1, type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=1:10, type=:stack)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=1, type=:stack)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=1:10, type=:topo)
@test p isa GLMakie.Figure

@info "test 4/21: plot_filter_response()"
p = NeuroAnalyzer.plot_filter_response(fs=sr(eeg), fprototype=:butterworth, ftype=:hp, cutoff=10, order=8)
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 5/21: plot_locs()"
p = NeuroAnalyzer.plot_locs(e10)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_locs(e10, selected=1:4)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_locs(e10, threed=true)
@test p isa Makie.Figure

@info "test 6/21: plot_psd()"
p = NeuroAnalyzer.plot_psd(e10, norm=true, ep=1, ch=1)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, norm=true, ep=1, ch=1, method=:mw)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, norm=true, ep=1, ch=1, method=:mt)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, norm=true, ep=1, ch=1, ref=:delta)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, norm=true, ep=1, ch=1:10, type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, norm=true, ep=1, ch=1:10, type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, norm=true, ep=1, ch=1:10, type=:w3d)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, norm=true, ep=1, ch=1:10, type=:topo)
@test p isa Makie.Figure

@info "test 7/21: plot_save()"
p = NeuroAnalyzer.plot(e10)
plot_save(p, file_name="test.png")
@test isfile("test.png") == true
isfile("test.png") && rm("test.png")

@info "test 8/21: plot()"
p = NeuroAnalyzer.plot(e10)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot(e10, ch=1:19, type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot(e10, ch=1:19, type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 9/21: plot_spectrogram()"
p = NeuroAnalyzer.plot_spectrogram(e10, norm=true, ep=1, ch=1)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, norm=true, ep=1, ch=1, method=:stft)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, norm=true, ep=1, ch=1, method=:mt)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, norm=true, ep=1, ch=1, method=:mw)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, norm=true, ep=1, ch=1:10)
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 10/21: plot_topo()"
p = NeuroAnalyzer.plot_topo(e10, seg=(0, 1))
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_topo(e10, seg=(0, 1), amethod=:median)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_topo(e10, seg=(0, 1), amethod=:median, imethod=:nn)

@info "test 11/21: plot_matrix()"
c = corm(e10)
channels = signal_channels(e10)
p = NeuroAnalyzer.plot_matrix(c[:, :, 1, 1], xlabels=labels(e10)[channels], ylabels=labels(e10)[channels])
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 12/21: plot_covmatrix()"
ac, lags = acov(e10, lag=5, norm=false)
p = NeuroAnalyzer.plot_covmatrix(ac[1, :, 1], lags)
@test p isa Plots.Plot{Plots.GRBackend}
cc, lags = xcov(eeg, lag=5, norm=false)
p = NeuroAnalyzer.plot_covmatrix(cc[1, :, 1], lags)
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 13/21: plot_histogram()"
stats = rand(2, 10)
p = NeuroAnalyzer.plot_histogram(stats[1, :])
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 14/21: plot_bar()"
p = NeuroAnalyzer.plot_bar(stats[:, 1], xlabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 15/21: plot_line()"
p = NeuroAnalyzer.plot_line(stats[:, 1], xlabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 16/21: plot_dots()"
p = NeuroAnalyzer.plot_dots([stats[1, :], stats[2, :]], glabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 17/21: plot_paired()"
p = NeuroAnalyzer.plot_paired([stats[1, :], stats[2, :]], glabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 18/21: plot_polar()"
p = NeuroAnalyzer.plot_polar(stats')
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 19/21: plot_weights()"
p = NeuroAnalyzer.plot_weights(e10, weights=rand(19), channel=1:19)
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 20/21: plot_dipole2d()"
d = NeuroAnalyzer.DIPOLE((0, 0, 0), (1, 1, 1))
p = NeuroAnalyzer.plot_dipole2d(d)
@test p isa Plots.Plot{Plots.GRBackend}

@info "test 21/21: plot_dipole3d()"
d = NeuroAnalyzer.DIPOLE((0, 0, 1), (1, 1, 1))
p = NeuroAnalyzer.plot_dipole3d(d)
@test p isa Makie.Figure

true