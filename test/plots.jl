using NeuroAnalyzer
using Plots
using Cairo
using Test

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
load_locs!(e10, file_name=joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))

isfile("test.png") && rm("test.png")

@info "Test: plot_compose()"
p1 = NeuroAnalyzer.plot(e10, ch="Fp1", ep=1, xlabel="")
p2 = plot_empty()
pp = [p1, p2]
l = (2, 1)
p = NeuroAnalyzer.plot_compose(pp, layout=l)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_connections()"
p = NeuroAnalyzer.plot_locs(e10, ch="eeg", connections=rand(19, 19), threshold=0.5)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_erp()"
e10_erp = average_epochs(e10)
p = NeuroAnalyzer.plot_erp(e10_erp, ch="Fp1")
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=["Fp1", "Fp2"], type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch="Fp1", type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=["Fp1", "Fp2"], type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch="Fp1", type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=["Fp1", "Fp2"], type=:stack)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch="Fp1", type=:stack)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_erp(e10_erp, ch=["Fp1", "Fp2"], type=:topo)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_filter_response()"
p = NeuroAnalyzer.plot_filter_response(fs=sr(eeg), fprototype=:butterworth, ftype=:hp, cutoff=10, order=8)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_locs()"
p = NeuroAnalyzer.plot_locs(e10, ch="eeg")
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_locs(e10, ch="eeg", selected=["Fp1", "Fp2"])
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_locs(e10, ch="eeg", d=3, interactive=false)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_psd()"
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch="Fp1")
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch="Fp1", method=:mw)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch="Fp1", method=:mt)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch="Fp1", method=:stft)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch="Fp1", method=:fft)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch="Fp1", method=:gh)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch="Fp1", method=:cwt)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch="Fp1", ref=:delta)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], type=:w3d)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], type=:topo)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_save()"
p = NeuroAnalyzer.plot(e10, ch="Fp1")
plot_save(p, file_name="test.png")
@test isfile("test.png")
isfile("test.png") && rm("test.png")

@info "Test: plot()"
p = NeuroAnalyzer.plot(e10, ch="Fp1")
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot(e10, ch=["Fp1", "Fp2"], type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot(e10, ch=["Fp1", "Fp2"], type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_spectrogram()"
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch="Fp1")
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:stft)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:mt)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:mw)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:gh)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:cwt)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_spectrogram(e10, db=true, ep=1, ch=["Fp1", "Fp2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_topo()"
p = NeuroAnalyzer.plot_topo(e10, ch="eeg", seg=(0, 1))
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_topo(e10, ch="eeg", seg=(0, 1), amethod=:median)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_topo(e10, ch="eeg", seg=(0, 1), amethod=:median, imethod=:nn)

@info "Test: plot_matrix()"
channels = get_channel(e10, ch=get_channel(e10, type=datatype(e10)))
c = corm(e10, ch="eeg")
p = NeuroAnalyzer.plot_matrix(c[:, :, 1, 1], xlabels=labels(e10)[channels], ylabels=labels(e10)[channels])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_xac()"
ac, lags = acov(e10, ch="eeg")
p = NeuroAnalyzer.plot_xac(ac[1, :, 1], lags)
@test p isa Plots.Plot{Plots.GRBackend}
xc, lags = xcov(e10, e10, ch1="Fp1", ch2="Fp2")
p = NeuroAnalyzer.plot_xac(xc[1, :, 1], lags)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_histogram()"
stats = rand(2, 10)
p = NeuroAnalyzer.plot_histogram(stats[1, :])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_bar()"
p = NeuroAnalyzer.plot_bar(stats[:, 1], xlabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_line()"
p = NeuroAnalyzer.plot_line(stats[:, 1], xlabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_dots()"
p = NeuroAnalyzer.plot_dots([stats[1, :], stats[2, :]], glabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_paired()"
p = NeuroAnalyzer.plot_paired([stats[1, :], stats[2, :]], glabels=["1", "2"])
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_polar()"
p = NeuroAnalyzer.plot_polar(stats')
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_weights()"
p = NeuroAnalyzer.plot_locs(e10, weights=rand(19), ch="eeg")
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_dipole2d()"
d = NeuroAnalyzer.DIPOLE((0, 0, 0), (1, 1, 1))
p = NeuroAnalyzer.plot_dipole2d(d)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_dipole3d()"
d = NeuroAnalyzer.DIPOLE((0, 0, 1), (1, 1, 1))
p = NeuroAnalyzer.plot_dipole3d(d)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_eros()"
s, f, t = eros(e10, ch="Fp1")
p = NeuroAnalyzer.plot_eros(s, f, t)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_erop()"
p, f = erop(e10, ch="Fp1")
p = NeuroAnalyzer.plot_erop(p, f)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot(obj1, obj2)"
p = NeuroAnalyzer.plot(e10, e10, ch="all")
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_icatopo()"
eeg_new = keep_epoch(e10, ep=1)
ic, ic_mw, ic_var = ica_decompose(eeg_new, ch="eeg", iter=10)
p = plot_icatopo(eeg_new, ch="eeg", ic, ic_mw)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: add_plot_locs()"
c = add_plot_locs(NeuroAnalyzer.plot(e10, ch="Fp1"), NeuroAnalyzer.plot_locs(e10, ch="Fp1", large=false), view=false, file_name="")
@test c isa Cairo.CairoSurfaceBase{UInt32}

@info "Test: plot2canvas()"
@test plot2canvas(p) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test: resize_canvas()"
c = plot2canvas(p)
@test resize_canvas(c, r=2) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test: add_topmargin_canvas()"
c = plot2canvas(p)
@test add_topmargin_canvas(c, c) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test: add_to_canvas()"
c = plot2canvas(p)
@test add_to_canvas(c, c, x=0, y=0, view=false) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test: plot_mep()"
mep = import_duomag(joinpath(testfiles_path, "mep-duomag.m"))
p = plot_mep(mep, ch="MEP1")
@test p isa Plots.Plot{Plots.GRBackend}
p = plot_mep(mep, ch=["MEP1", "MEP2"], type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = plot_mep(mep, ch=["MEP1", "MEP2"], type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = plot_mep(mep, ch=["MEP1", "MEP2"], type=:stack)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_ci()"
s = eeg.data[1, 1:100, :]
t = eeg.epoch_time[1:100]
s_avg, s_l, s_u = NeuroStats.bootstrap_ci(s, cl=0.95)
p = plot_ci(s_avg, s_l, s_u, t)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_phsd()"
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch="Fp1")
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch=["Fp1", "Fp2"], type=:mean)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch=["Fp1", "Fp2"], type=:butterfly)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch=["Fp1", "Fp2"], type=:w3d)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_phsd(e10, ep=1, ch=["Fp1", "Fp2"], type=:topo)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_coherence()"
coh, mscoh, f = coherence(e10, e10, ch1=["Fp1", "Fp2"], ch2=["Fp1", "Fp2"], ep1=1, ep2=1, frq_lim=(45, 55))
p = NeuroAnalyzer.plot_coherence(coh[1, :, 1], f)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_coherence(coh[:, :, 1], f)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_coherence_avg(coh[:, :, 1], f)
@test p isa Plots.Plot{Plots.GRBackend}
p = NeuroAnalyzer.plot_coherence_butterfly(coh[:, :, 1], f)
@test p isa Plots.Plot{Plots.GRBackend}

@info "Test: plot_heatmap()"
m = rand(nchannels(e10), epoch_len(e10))
p = NeuroAnalyzer.plot_heatmap(m, x=e10.epoch_time, y=1:nchannels(e10))
@test p isa Plots.Plot{Plots.GRBackend}

true
