using NeuroAnalyzer
using Cairo
using Test
using GLMakie

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
load_locs!(e10, file_name=joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))

isfile("test.png") && rm("test.png")

@info "Test: plot_compose()"
p1 = NeuroAnalyzer.plot(e10, ch="Fp1", ep=1, n_epochs=1, xlabel="", gui=false);
p2 = NeuroAnalyzer.plot(e10, ch="Fp2", ep=1, n_epochs=1, xlabel="", gui=false);
p3 = NeuroAnalyzer.plot(e10, ch="Cz", ep=1, n_epochs=1, xlabel="", gui=false);
pp = [p1, p2, p3]
l = (2, 2)
p = NeuroAnalyzer.plot_compose(pp, layout=l)
@test p isa GLMakie.Figure

@info "Test: plot_connections()"
p = plot_locs(e10, ch="eeg", connections=rand(19, 19), threshold=0.5, gui=false)
@test p isa GLMakie.Figure
p = plot_locs(e10, ch="eeg", connections=rand(-10:1:10, 19, 19), threshold=12, threshold_type=:g, gui=false)
@test p isa GLMakie.Figure

@info "Test: plot_erp()"
e10_erp = average_epochs(e10)
p = plot_erp(e10_erp, ch="Fp1")
@test p isa GLMakie.Figure
p = plot_erp(e10_erp, ch=["Fp1", "Fp2"], avg=true)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp, ch=["Fp1", "Fp2"], avg=false)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp, ch=["Fp1", "Fp2"], ci95=true)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp, ch=["Fp1", "Fp2"], type=:stack)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp, ch="Fp1", type=:stack)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp, ch=["Fp1", "Fp2"], type=:topo)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp, ch=["Fp1", "Fp2"], type=:gfp)
@test p isa GLMakie.Figure

@info "Test: plot_filter()"
p = plot_filter(fs=256, fprototype=:butterworth, ftype=:hp, cutoff=10, order=8, gui=false)
@test p isa GLMakie.Figure
p = plot_filter(fs=256, fprototype=:butterworth, ftype=:lp, cutoff=10, order=8, gui=false)
@test p isa GLMakie.Figure
p = plot_filter(fs=256, fprototype=:butterworth, ftype=:bp, cutoff=(10, 12), order=8, gui=false)
@test p isa GLMakie.Figure
p = plot_filter(fs=256, fprototype=:butterworth, ftype=:bs, cutoff=(10, 12), order=8, gui=false)
@test p isa GLMakie.Figure

@info "Test: plot_locs()"
p = plot_locs(e10, ch="eeg", gui=false)
@test p isa GLMakie.Figure
p = plot_locs(e10, ch="eeg", sch=["Fp1", "Fp2"], gui=false)
@test p isa GLMakie.Figure

@info "Test: plot_psd()"
p = plot_psd(e10, db=true, ep=1, ch="Fp1")
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch="Fp1", method=:mw)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch="Fp1", method=:mt)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch="Fp1", method=:stft)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch="Fp1", method=:fft)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch="Fp1", method=:gh)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch="Fp1", ref=:delta)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], avg=true)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], avg=false)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], ci95=true)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], type=:w3d)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], type=:s3d)
@test p isa GLMakie.Figure
p = plot_psd(e10, db=true, ep=1, ch=["Fp1", "Fp2"], type=:topo)
@test p isa GLMakie.Figure

@info "Test: plot_save()"
p = NeuroAnalyzer.plot(e10, ch="Fp1")
NeuroAnalyzer.plot_save(p, file_name="test.png")
@test isfile("test.png")
isfile("test.png") && rm("test.png")

@info "Test: plot()"
t = 0:0.1:1
s = sin.(t)
p = NeuroAnalyzer.plot(t, s)
@test p isa GLMakie.Figure
p = NeuroAnalyzer.plot(e10, ch="Fp1", gui=false)
@test p isa GLMakie.Figure
p = NeuroAnalyzer.plot(e10, ch=["Fp1", "Fp2"], type=:butterfly)
@test p isa GLMakie.Figure
p = NeuroAnalyzer.plot(e10, ch=["Fp1", "Fp2"], type=:butterfly, avg=false)
@test p isa GLMakie.Figure
p = NeuroAnalyzer.plot(e10, ch=["Fp1", "Fp2"], type=:butterfly, ci95=true)
@test p isa GLMakie.Figure

@info "Test: plot_spectrogram()"
p = plot_spectrogram(e10, db=true, ep=1, ch="Fp1")
@test p isa GLMakie.Figure
p = plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:stft)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:mt)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:mw)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:gh)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:cwt)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10, db=true, ep=1, ch="Fp1", method=:hht)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10, db=true, ep=1, ch=["Fp1", "Fp2"])
@test p isa GLMakie.Figure
p = plot_spectrogram(e10, db=true, ep=1, ch="eeg", type=:topo)
@test p isa GLMakie.Figure

@info "Test: plot_topo()"
p = plot_topo(e10, ch="eeg", tpos=0)
@test p isa GLMakie.Figure
p = plot_topo(e10, ch="eeg", tpos=0, imethod=:nn)
@test p isa GLMakie.Figure

@info "Test: plot_matrix()"
channels = get_channel(e10, ch=get_channel(e10, type=datatype(e10)))
c = corm(e10, ch="eeg")
p = plot_matrix(c[:, :, 1, 1], xlabels=labels(e10)[channels], ylabels=labels(e10)[channels])
@test p isa GLMakie.Figure

@info "Test: plot_xac()"
ac, lags = acov(e10, ch="eeg")
p = plot_xac(ac[1, :, 1], lags)
@test p isa GLMakie.Figure
xc, lags = xcov(e10, e10, ch1="Fp1", ch2="Fp2")
p = plot_xac(xc[1, :, 1], lags)
@test p isa GLMakie.Figure

@info "Test: plot_histogram()"
stats = rand(2, 100)
p = plot_histogram(stats[1, :], 0.8)
@test p isa GLMakie.Figure

@info "Test: plot_bar()"
p = plot_bar([1, 2, 3, 4, 5], xlabels=["G1", "G2", "G3", "G4", "G5"])
@test p isa GLMakie.Figure

@info "Test: plot_line()"
p = plot_line([1, 2, 1, 4.1, 1.5], xlabels=["G1", "G2", "G3", "G4", "G5"])
@test p isa GLMakie.Figure
p = plot_line(rand(2, 5), rlabels=["pre", "post"], xlabels=["G1", "G2", "G3", "G4", "G5"])
@test p isa GLMakie.Figure

@info "Test: plot_box()"
s = rand(3, 10)
s[1, :] .*= 2
p = plot_box(s, xlabels=["G1", "G2", "G3"])
@test p isa GLMakie.Figure

@info "Test: plot_violin()"
s = rand(3, 10)
s[1, :] .*= 2
p = plot_violin(s, xlabels=["G1", "G2", "G3"])
@test p isa GLMakie.Figure

@info "Test: plot_dots()"
s = rand(2, 10)
p = plot_dots(s, xlabels=["G1", "G2"])
@test p isa GLMakie.Figure

@info "Test: plot_paired()"
s = rand(3, 10)
p = plot_paired(s, xlabels=["V1", "V2", "V3"])
@test p isa GLMakie.Figure

@info "Test: plot_polar()"
s = rand(0:0.1:10, 10, 2)
p = plot_polar(s)
@test p isa GLMakie.Figure
s = rand(10)
p = plot_polar(s)
@test p isa GLMakie.Figure
p = plot_polar(s, m=(1, 1))
@test p isa GLMakie.Figure
s = rand(0:0.1:10, 10, 2)
p = plot_polar(s, m=(1, 1))
@test p isa GLMakie.Figure

@info "Test: plot_weights()"
p = plot_locs(e10, weights=rand(19), ch="eeg", gui=false)
@test p isa GLMakie.Figure

@info "Test: plot_dipole2d()"
d = NeuroAnalyzer.DIPOLE((0, 0, 0), (1, 1, 1))
p = plot_dipole2d(d)
@test p isa GLMakie.Figure

@info "Test: plot_dipole3d()"
d = NeuroAnalyzer.DIPOLE((0, 0, 1), (1, 1, 1))
p = plot_dipole3d(d)
@test p isa GLMakie.Figure

@info "Test: plot_eros()"
s, f, t = eros(e10, ch="Fp1")
p = plot_eros(s, f, t, tm=1000)
@test p isa GLMakie.Figure
e10_erp = average_epochs(e10)
s, f, t = eros(e10_erp, ch="Fp1")
p = plot_eros(s, f, t)
@test p isa GLMakie.Figure

@info "Test: plot_erop()"
sp, sf = erop(e10, ch="Fp1")
p = plot_erop(sp, sf)
@test p isa GLMakie.Figure
e10_erp = average_epochs(e10)
sp, sf = erop(e10_erp, ch="Fp1")
p = plot_erop(sp, sf)
@test p isa GLMakie.Figure

@info "Test: plot(obj1, obj2)"
p = NeuroAnalyzer.plot(eeg, eeg, ch="all", gui=false)
@test p isa GLMakie.Figure

@info "Test: plot_icatopo()"
eeg_new = keep_epoch(e10, ep=1)
ic, ic_mw, ic_var = ica_decompose(eeg_new, ch="eeg", iter=10)
p = plot_icatopo(eeg_new, ch="eeg", ic=ic, ic_mw=ic_mw, ic_idx=1:3, tpos=0)
@test p isa GLMakie.Figure

@info "Test: add_pl()"
p = NeuroAnalyzer.plot(e10, ep=1)
pl = plot_locs(eeg, ch="eeg", sch="eeg", ps=:s, gui=false)
pp = add_pl(p, pl)
@test pp isa GLMakie.Figure

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
@test p isa GLMakie.Figure
p = plot_mep(mep, ch=["MEP1", "MEP2"], avg=true)
@test p isa GLMakie.Figure
p = plot_mep(mep, ch=["MEP1", "MEP2"], avg=false)
@test p isa GLMakie.Figure
p = plot_mep(mep, ch=["MEP1", "MEP2"], ci95=true)
@test p isa GLMakie.Figure
p = plot_mep(mep, ch="all", type=:stack)
@test p isa GLMakie.Figure

@info "Test: plot_ci()"
s = eeg.data[1, 1:100, :]
t = eeg.epoch_time[1:100]
s_avg, s_l, s_u = NeuroAnalyzer.bootstrap_ci(s, cl=0.95)
p = plot_ci(s_avg, s_l, s_u, t)
@test p isa GLMakie.Figure

@info "Test: plot_phsd()"
p = plot_phsd(e10, ep=1, ch="Fp1")
@test p isa GLMakie.Figure
p = plot_phsd(e10, ep=1, ch="eeg", avg=true)
@test p isa GLMakie.Figure
p = plot_phsd(e10, ep=1, ch="eeg", avg=false)
@test p isa GLMakie.Figure
p = plot_phsd(e10, ep=1, ch="eeg", ci95=true)
@test p isa GLMakie.Figure
p = plot_phsd(e10, ep=1, ch="eeg", type=:w3d)
@test p isa GLMakie.Figure
p = plot_phsd(e10, ep=1, ch="eeg", type=:topo)
@test p isa GLMakie.Figure

@info "Test: plot_coherence()"
coh, imcoh, mscoh, f = coherence(e10, e10, ch1=["Fp1", "Fp2"], ch2=["Fp1", "Fp2"], ep1=1, ep2=1, flim=(15, 25))
p = plot_coherence(abs.(coh[1, :, 1]), f)
@test p isa GLMakie.Figure
p = plot_coherence(abs.(coh[:, :, 1]), f)
@test p isa GLMakie.Figure
p = plot_coherence(abs.(coh[:, :, 1]), f, avg=true)
@test p isa GLMakie.Figure
p = plot_coherence(abs.(coh[:, :, 1]), f, avg=false)
@test p isa GLMakie.Figure
p = plot_coherence(abs.(coh[:, :, 1]), f, ci95=true)
@test p isa GLMakie.Figure

@info "Test: plot_heatmap()"
m = rand(nchannels(e10), epoch_len(e10))
p = plot_heatmap(m, x=e10.epoch_time, y=1:nchannels(e10))
@test p isa GLMakie.Figure

@info "Test: plot_connectivity_circle()"
l = get_channel(eeg, type="eeg")
m = rand(-10:0.1:10, length(l), length(l))
p = plot_connectivity_circle(m, clabels=l)
@test p isa GLMakie.Figure

@info "Test: plot_imf()"
imf = rand(5, 100)
t = collect(1:100)
p = plot_imf(imf, t=t)
@test p isa GLMakie.Figure

@info "Test: plot_gridlocs()"
p = plot_gridlocs()
@test p isa GLMakie.Figure

@info "Test: plot_hs()"
hms, t = hmspectrum(e10, ch="Fp1")
hms = vec(hms[:, :, 1])
p = plot_hs(hms, t)
@test p isa GLMakie.Figure

@info "Test: plot_fi()"
t = e10.epoch_time
fi = frqinst(e10, ch="Fp1")
p = plot_fi(fi[1, :, 1], t)
@test p isa GLMakie.Figure

@info "Test: plot_phase()"
s = e10.data[1, :, 1]
X = ftransform(s)
f, _ = freqs(s, sr(e10))
p = plot_phase(rad2deg.(X.ph[1:100]), f[1:100], unit=:rad, type=:stem)
@test p isa GLMakie.Figure
p = plot_phase(rad2deg.(X.ph[1:100]), f[1:100], unit=:deg, type=:stem)
@test p isa GLMakie.Figure
p = plot_phase(rad2deg.(X.ph[1:100]), f[1:100], unit=:rad, type=:line)
@test p isa GLMakie.Figure
p = plot_phase(rad2deg.(X.ph[1:100]), f[1:100], unit=:deg, type=:line)
@test p isa GLMakie.Figure

@info "Test: plot_locs3d()"
p = plot_locs3d(e10, ch="eeg", mesh_type=:disabled, gui=false)
@test p isa GLMakie.Figure
p = plot_locs3d(e10, ch="eeg", mesh_type=:brain, gui=false)
@test p isa GLMakie.Figure
p = plot_locs3d(e10, ch="eeg", mesh_type=:head, gui=false)
@test p isa GLMakie.Figure

@info "Test: plot_polezero()"
p = plot_polezero(rand(ComplexF64, 2), rand(ComplexF64, 2))
@test p isa GLMakie.Figure

@info "Test: plot_dwc()"
dc = rand(16, 100)
t = collect(1:100)
p = plot_dwc(dc, t=t, n=4)
@test p isa GLMakie.Figure

true
