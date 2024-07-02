using NeuroAnalyzer
using Plots
using Test

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
load_locs!(e10, file_name=joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))

isfile("test.png") && rm("test.png")

@info "Test 1/34: plot_compose()"
p1 = NeuroAnalyzer.plot(e10, ep=1, xlabel="")
p2 = plot_empty()
pp = [p1, p2]
l = (2, 1)
p = NeuroAnalyzer.plot_compose(pp, layout=l)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 2/34: plot_connections()"
p = NeuroAnalyzer.plot_connections(e10, connections=rand(19, 19), ch=1:19, threshold=0.5)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 3/34: plot_erp()"
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
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 4/34: plot_filter_response()"
p = NeuroAnalyzer.plot_filter_response(fs=sr(eeg), fprototype=:butterworth, ftype=:hp, cutoff=10, order=8)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 5/34: plot_locs()"
p = NeuroAnalyzer.plot_locs(e10)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_locs(e10, selected=1:4)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_locs(e10, threed=true, interactive=false)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 6/34: plot_psd()"
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1, method=:mw)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1, method=:mt)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1, method=:stft)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1, method=:fft)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1, ref=:delta)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1:10, type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1:10, type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1:10, type=:w3d)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=1:10, type=:topo)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 7/34: plot_save()"
p = NeuroAnalyzer.plot(e10)
plot_save(p, file_name="test.png")
@test isfile("test.png")
isfile("test.png") && rm("test.png")

@info "Test 8/34: plot()"
p = NeuroAnalyzer.plot(e10)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot(e10, ch=1:19, type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot(e10, ch=1:19, type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 9/34: plot_spectrogram()"
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch=1)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch=1, method=:stft)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch=1, method=:mt)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch=1, method=:mw)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch=1, method=:gh)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch=1, method=:cwt)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch=1:10)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 10/34: plot_topo()"
p = NeuroAnalyzer.plot_topo(e10, seg=(0, 1))
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_topo(e10, seg=(0, 1), amethod=:median)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_topo(e10, seg=(0, 1), amethod=:median, imethod=:nn)

@info "Test 11/34: plot_matrix()"
channels = get_channel_bytype(e10, type=datatype(e10))
c = corm(e10, ch=channels)
p = NeuroAnalyzer.plot_matrix(c[:, :, 1, 1], xlabels=labels(e10)[channels], ylabels=labels(e10)[channels])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 12/34: plot_xac()"
ac, lags = acov(e10)
p = NeuroAnalyzer.plot_xac(ac[1, :, 1], lags)
@test p isa Plots.Plot{Plots.GRBackend}
xc, lags = xcov(e10, e10, ch1=1, ch2=2)
p = NeuroAnalyzer.plot_xac(xc[1, :, 1], lags)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 13/34: plot_histogram()"
stats = rand(2, 10)
p = NeuroAnalyzer.plot_histogram(stats[1, :])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 14/34: plot_bar()"
p = NeuroAnalyzer.plot_bar(stats[:, 1], xlabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 15/34: plot_line()"
p = NeuroAnalyzer.plot_line(stats[:, 1], xlabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 16/34: plot_dots()"
p = NeuroAnalyzer.plot_dots([stats[1, :], stats[2, :]], glabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 17/34: plot_paired()"
p = NeuroAnalyzer.plot_paired([stats[1, :], stats[2, :]], glabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 18/34: plot_polar()"
p = NeuroAnalyzer.plot_polar(stats')
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 19/34: plot_weights()"
p = NeuroAnalyzer.plot_weights(e10, weights=rand(19), ch=1:19)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 20/34: plot_dipole2d()"
d = NeuroAnalyzer.DIPOLE((0, 0, 0), (1, 1, 1))
p = NeuroAnalyzer.plot_dipole2d(d)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 21/34: plot_dipole3d()"
d = NeuroAnalyzer.DIPOLE((0, 0, 1), (1, 1, 1))
p = NeuroAnalyzer.plot_dipole3d(d)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 22/34: plot_eros()"
s, f, t = eros(e10, ch=1)
p = NeuroAnalyzer.plot_eros(s, f, t)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 23/34: plot_erop()"
p, f = erop(e10, ch=1)
p = NeuroAnalyzer.plot_erop(p, f)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 24/34: plot(obj1, obj2)"
p = NeuroAnalyzer.plot(e10, e10, ch=1)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 25/34: plot_icatopo()"
eeg_new = keep_epoch(e10, ep=1)
ic, ic_mw, ic_var = ica_decompose(eeg_new, iter=10)
p = plot_icatopo(eeg_new, ic, ic_mw)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 26/34: add_locs()"
c = add_locs(NeuroAnalyzer.plot(eeg), NeuroAnalyzer.plot_locs(eeg), view=false, file_name="")
@test c isa Cairo.CairoSurfaceBase{UInt32}

@info "Test 27/34: plot2canvas()"
@test plot2canvas(p) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test 28/34: resize_canvas()"
c = plot2canvas(p)
@test resize_canvas(c, r=2) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test 29/34: add_topmargin_canvas()"
c = plot2canvas(p)
@test add_topmargin_canvas(c, c) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test 30/34: add_to_canvas()"
c = plot2canvas(p)
@test add_to_canvas(c, c, x=0, y=0, view=false) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test 31/34: plot_mep()"
mep = import_duomag(joinpath(testfiles_path, "mep-duomag.m"))
p = plot_mep(mep, ch=1)
@test p isa Plots.Plot{Plots.GRBackend}
p = plot_mep(mep, ch=1:10, type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = plot_mep(mep, ch=1:10, type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = plot_mep(mep, ch=1:10, type=:stack)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 32/34: plot_ci()"
s = eeg.data[1, 1:100, :]
t = eeg.epoch_time[1:100]
s_avg, s_l, s_u = bootstrap_ci(s, ci=0.95)
p = plot_ci(s_avg, s_l, s_u, t)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 33/34: plot_phsd()"
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch=1)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch=1:10, type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch=1:10, type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch=1:10, type=:w3d)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch=1:10, type=:topo)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test 34/34: plot_coherence()"
coh, mscoh, f = coherence(e10, e10, ch1=1:3, ch2=2:4, ep1=1, ep2=1, frq_lim=(45, 55))
p = NeuroAnalyzer.plot_coherence(coh[1, :, 1], f)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_coherence(coh[:, :, 1], f)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_coherence_avg(coh[:, :, 1], f)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_coherence_butterfly(coh[:, :, 1], f)
@test p isa Plots.Plot{Plots.GRBackend}

true