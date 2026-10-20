using NeuroAnalyzer
using Cairo
using Test
using GLMakie

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)
load_locs!(e10; file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"))

isfile("test.png") && rm("test.png")

@info "Test: plot_compose()"
p1 = NeuroAnalyzer.plot(e10; ch = "Fp1", ep = 1, n_epochs = 1, xlabel = "", gui = false);
p2 = NeuroAnalyzer.plot(e10; ch = "Fp2", ep = 1, n_epochs = 1, xlabel = "", gui = false);
p3 = NeuroAnalyzer.plot(e10; ch = "Cz", ep = 1, n_epochs = 1, xlabel = "", gui = false);
pp = [p1, p2, p3]
l = (2, 2)
p = NeuroAnalyzer.plot_compose(pp; layout = l)
@test p isa GLMakie.Figure

@info "Test: plot_connections()"
p = plot_locs(e10; ch = "eeg", connections = rand(19, 19), threshold = 0.5, gui = false)
@test p isa GLMakie.Figure
p = plot_locs(
    e10;
    ch = "eeg",
    connections = rand(-10:1:10, 19, 19),
    threshold = 12,
    threshold_type = :g,
    gui = false,
)
@test p isa GLMakie.Figure

@info "Test: plot_erp()"
e10_erp = average_epochs(e10)
p = plot_erp(e10_erp; ch = "Fp1", gui = false)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp; ch = ["Fp1", "Fp2"], avg = true, gui = false)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp; ch = ["Fp1", "Fp2"], avg = false, gui = false)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp; ch = ["Fp1", "Fp2"], ci95 = true, gui = false)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp; ch = ["Fp1", "Fp2"], type = :stack, gui = false)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp; ch = "Fp1", type = :stack, gui = false)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp; ch = ["Fp1", "Fp2"], type = :topo, gui = false)
@test p isa GLMakie.Figure
p = plot_erp(e10_erp; ch = ["Fp1", "Fp2"], type = :gfp, gui = false)
@test p isa GLMakie.Figure

@info "Test: plot_filter()"
p = plot_filter(;
    fs = 256,
    fprototype = :butterworth,
    ftype = :hp,
    cutoff = 10,
    order = 8,
    gui = false,
)
@test p isa GLMakie.Figure
p = plot_filter(;
    fs = 256,
    fprototype = :chebyshev1,
    ftype = :hp,
    cutoff = 10,
    order = 8,
    gui = false,
)
@test p isa GLMakie.Figure
p = plot_filter(;
    fs = 256,
    fprototype = :chebyshev2,
    ftype = :hp,
    cutoff = 10,
    order = 8,
    gui = false,
)
@test p isa GLMakie.Figure
p = plot_filter(;
    fs = 256,
    fprototype = :elliptic,
    ftype = :hp,
    cutoff = 10,
    order = 8,
    rs = 4,
    rp = 2,
    gui = false,
)
@test p isa GLMakie.Figure
p = plot_filter(;
    fs = 256,
    fprototype = :iirnotch,
    ftype = :lp,
    cutoff = 10,
    bw = 0.5,
    gui = false,
)
@test p isa GLMakie.Figure
p = plot_filter(;
    fs = 256,
    fprototype = :fir,
    ftype = :lp,
    cutoff = 10,
    order = 11,
    bw = 0.5,
    gui = false,
)
@test p isa GLMakie.Figure
p = plot_filter(;
    fs = 256,
    fprototype = :fir,
    ftype = :hp,
    cutoff = 10,
    order = 11,
    bw = 0.5,
    gui = false,
)
@test p isa GLMakie.Figure
p = plot_filter(;
    fs = 256,
    fprototype = :firls,
    ftype = :bp,
    cutoff = (10, 12),
    order = 11,
    bw = 0.5,
    gui = false,
)
@test p isa GLMakie.Figure
p = plot_filter(;
    fs = 256,
    fprototype = :remez,
    ftype = :bs,
    cutoff = (10, 12),
    order = 11,
    bw = 0.5,
    gui = false,
)
@test p isa GLMakie.Figure

@info "Test: plot_locs()"
p = plot_locs(e10; ch = "eeg", gui = false)
@test p isa GLMakie.Figure
p = plot_locs(e10; ch = "eeg", sch = ["Fp1", "Fp2"], gui = false)
@test p isa GLMakie.Figure

@info "Test: plot_psd()"
p = plot_psd(e10; db = true, ep = 1, ch = "Fp1", gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = "Fp1", method = :mw, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = "Fp1", method = :mt, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = "Fp1", method = :stft, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = "Fp1", method = :fft, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = "Fp1", method = :gh, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = "Fp1", ref = :delta, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = ["Fp1", "Fp2"], avg = true, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = ["Fp1", "Fp2"], avg = false, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = ["Fp1", "Fp2"], ci95 = true, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = ["Fp1", "Fp2"], type = :w3d, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = ["Fp1", "Fp2"], type = :s3d, gui = false)
@test p isa GLMakie.Figure
p = plot_psd(e10; db = true, ep = 1, ch = ["Fp1", "Fp2"], type = :topo, gui = false)
@test p isa GLMakie.Figure

@info "Test: plot_save()"
p = NeuroAnalyzer.plot(e10; ch = "Fp1", gui = false)
NeuroAnalyzer.plot_save(p; file_name = "test.png")
@test isfile("test.png")
isfile("test.png") && rm("test.png")

@info "Test: plot_signal() / plot()"
t = 0:0.1:1
s = sin.(t)
p = plot_signal(t, s)
@test p isa GLMakie.Figure
p = NeuroAnalyzer.plot(e10; ch = "Fp1", gui = false)
@test p isa GLMakie.Figure
p = NeuroAnalyzer.plot(e10; ch = ["Fp1", "Fp2"], type = :butterfly, gui = false)
@test p isa GLMakie.Figure
p = NeuroAnalyzer.plot(
    e10;
    ch = ["Fp1", "Fp2"],
    type = :butterfly,
    avg = false,
    gui = false,
)
@test p isa GLMakie.Figure
p = NeuroAnalyzer.plot(
    e10;
    ch = ["Fp1", "Fp2"],
    type = :butterfly,
    avg = false,
    ci95 = true,
    gui = false,
)
@test p isa GLMakie.Figure

@info "Test: plot_spectrogram()"
p = plot_spectrogram(e10; db = true, ep = 1, ch = "Fp1", gui = false)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10; db = true, ep = 1, ch = "Fp1", method = :stft, gui = false)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10; db = true, ep = 1, ch = "Fp1", method = :mt, gui = false)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10; db = true, ep = 1, ch = "Fp1", method = :mw, gui = false)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10; db = true, ep = 1, ch = "Fp1", method = :gh, gui = false)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10; db = true, ep = 1, ch = "Fp1", method = :cwt, gui = false)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10; db = true, ep = 1, ch = "Fp1", method = :hht, gui = false)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10; db = true, ep = 1, ch = ["Fp1", "Fp2"], gui = false)
@test p isa GLMakie.Figure
p = plot_spectrogram(e10; db = true, ep = 1, ch = "eeg", type = :topo, gui = false)
@test p isa GLMakie.Figure

@info "Test: plot_topo()"
p = plot_topo(e10; ch = "eeg", tpos = 0)
@test p isa GLMakie.Figure
p = plot_topo(e10; ch = "eeg", tpos = 0, imethod = :nn)
@test p isa GLMakie.Figure

@info "Test: plot_matrix()"
channels = get_channel(e10; ch = get_channel(e10; type = datatype(e10)))
c = corm(e10; ch = "eeg")
p = plot_matrix(
    c[:, :, 1, 1];
    xlabels = labels(e10)[channels],
    ylabels = labels(e10)[channels],
)
@test p isa GLMakie.Figure

@info "Test: plot_xac()"
ac, lags = acov(e10; ch = "eeg")
p = plot_xac(ac[1, :, 1], lags)
@test p isa GLMakie.Figure
xc, lags = xcov(e10, e10; ch1 = "Fp1", ch2 = "Fp2")
p = plot_xac(xc[1, :, 1], lags)
@test p isa GLMakie.Figure

@info "Test: plot_histogram()"
stats = rand(2, 100)
p = plot_histogram(stats[1, :], 0.8)
@test p isa GLMakie.Figure

@info "Test: plot_bar()"
p = plot_bar([1, 2, 3, 4, 5]; glabels = ["G1", "G2", "G3", "G4", "G5"])
@test p isa GLMakie.Figure

@info "Test: plot_line()"
p = plot_line([1, 2, 1, 4.1, 1.5]; glabels = ["G1", "G2", "G3", "G4", "G5"])
@test p isa GLMakie.Figure
p = plot_line(
    rand(2, 5);
    rlabels = ["pre", "post"],
    glabels = ["G1", "G2", "G3", "G4", "G5"],
)
@test p isa GLMakie.Figure

@info "Test: plot_box()"
s = rand(3, 10)
s[1, :] .*= 2
p = plot_box(s; glabels = ["G1", "G2", "G3"])
@test p isa GLMakie.Figure

@info "Test: plot_violin()"
s = rand(3, 10)
s[1, :] .*= 2
p = plot_violin(s; glabels = ["G1", "G2", "G3"])
@test p isa GLMakie.Figure

@info "Test: plot_dots()"
s = rand(2, 10)
p = plot_dots(s; glabels = ["G1", "G2"])
@test p isa GLMakie.Figure

@info "Test: plot_paired()"
s = rand(3, 10)
p = plot_paired(s; glabels = ["V1", "V2", "V3"])
@test p isa GLMakie.Figure

@info "Test: plot_polar()"
s = rand(0:0.1:10, 10, 2)
p = plot_polar(s)
@test p isa GLMakie.Figure
s = rand(10)
p = plot_polar(s)
@test p isa GLMakie.Figure
p = plot_polar(s; m = (1, 1))
@test p isa GLMakie.Figure
s = rand(0:0.1:10, 10, 2)
p = plot_polar(s; m = (1, 1))
@test p isa GLMakie.Figure

@info "Test: plot_weights()"
p = plot_locs(e10; weights = rand(19), ch = "eeg", gui = false)
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
eros_data = eros(e10; ch = "Fp1")
s, f, t = eros_data.s, eros_data.f, eros_data.t
p = plot_eros(s, f, t; tm = 1000)
@test p isa GLMakie.Figure
e10_erp = average_epochs(e10)
s, f, t = eros(e10_erp; ch = "Fp1")
p = plot_eros(s, f, t)
@test p isa GLMakie.Figure

@info "Test: plot_erop()"
sp, sf = erop(e10; ch = "Fp1")
p = plot_erop(sp, sf)
@test p isa GLMakie.Figure
e10_erp = average_epochs(e10)
sp, sf = erop(e10_erp; ch = "Fp1")
p = plot_erop(sp, sf)
@test p isa GLMakie.Figure

@info "Test: plot(obj1, obj2)"
p = NeuroAnalyzer.plot(eeg, eeg; ch = "all", gui = false)
@test p isa GLMakie.Figure

@info "Test: plot_icatopo()"
eeg_new = keep_epoch(e10; ep = 1)
ica_data = ica_decompose(eeg_new; ch = "eeg", iter = 10)
p = plot_icatopo(
    eeg_new;
    ch = "eeg",
    ic = ica_data.ic,
    ic_mw = ica_data.ic_mw,
    ic_idx = 2:3,
    tpos = 8.01,
)
@test p isa GLMakie.Figure

@info "Test: add_inset()"
p = NeuroAnalyzer.plot(e10; ep = 1, gui = false)
pl = plot_locs(eeg; ch = "eeg", sch = "eeg", ps = :s, gui = false)
pp = add_inset(p, pl)
@test pp isa GLMakie.Figure

@info "Test: plot2canvas()"
@test plot2canvas(p) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test: resize_canvas()"
c = plot2canvas(p)
@test resize_canvas(c, r = 2) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test: add_topmargin_canvas()"
c = plot2canvas(p)
@test add_topmargin_canvas(c, c) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test: add_to_canvas()"
c = plot2canvas(p)
@test add_to_canvas(c, c, x = 0, y = 0) isa Cairo.CairoSurfaceBase{UInt32}

@info "Test: plot_mep()"
mep = import_duomag(joinpath(testfiles_path, "mep", "mep-duomag.m"))
p = plot_mep(mep; ch = "MEP1")
@test p isa GLMakie.Figure
p = plot_mep(mep; ch = ["MEP1", "MEP2"], avg = true)
@test p isa GLMakie.Figure
p = plot_mep(mep; ch = ["MEP1", "MEP2"], avg = false)
@test p isa GLMakie.Figure
p = plot_mep(mep; ch = ["MEP1", "MEP2"], avg = false, ci95 = true)
@test p isa GLMakie.Figure
p = plot_mep(mep; ch = "all", type = :stack)
@test p isa GLMakie.Figure

@info "Test: plot_ci()"
s = eeg.data[1, 1:100, :]
t = eeg.epoch_time[1:100]
s_avg, s_l, s_u = NeuroAnalyzer.bootstrap_ci(s; cl = 0.95)
p = plot_ci(s_avg, s_l, s_u, t)
@test p isa GLMakie.Figure

@info "Test: plot_phs()"
p = plot_phs(e10; ep = 1, ch = "Fp1", gui = false)
@test p isa GLMakie.Figure
p = plot_phs(e10; ep = 1, ch = "eeg", avg = true, gui = false)
@test p isa GLMakie.Figure
p = plot_phs(e10; ep = 1, ch = "eeg", avg = false, gui = false)
@test p isa GLMakie.Figure
p = plot_phs(e10; ep = 1, ch = "eeg", ci95 = true, gui = false)
@test p isa GLMakie.Figure
p = plot_phs(e10; ep = 1, ch = "eeg", type = :w3d, gui = false)
@test p isa GLMakie.Figure
p = plot_phs(e10; ep = 1, ch = "eeg", type = :topo, gui = false)
@test p isa GLMakie.Figure

@info "Test: plot_coherence()"
coh, imcoh, mscoh, f = NeuroAnalyzer.coherence(
    e10,
    e10;
    ch1 = ["Fp1", "Fp2"],
    ch2 = ["Fp1", "Fp2"],
    flim = (15, 25),
)
p = plot_coherence(abs.(coh[1, :, 1]), f)
@test p isa GLMakie.Figure
p = plot_coherence(abs.(coh[:, :, 1]), f)
@test p isa GLMakie.Figure
p = plot_coherence(abs.(coh[:, :, 1]), f; avg = true)
@test p isa GLMakie.Figure
p = plot_coherence(abs.(coh[:, :, 1]), f; avg = false)
@test p isa GLMakie.Figure
p = plot_coherence(abs.(coh[:, :, 1]), f; ci95 = true)
@test p isa GLMakie.Figure

@info "Test: plot_heatmap()"
m = rand(nchannels(e10), epoch_len(e10))
p = plot_heatmap(m; x = e10.epoch_time, y = 1:nchannels(e10))
@test p isa GLMakie.Figure

@info "Test: plot_connectivity_circle()"
l = get_channel(eeg; type = "eeg")
m = rand(-10:0.1:10, length(l), length(l))
p = plot_connectivity_circle(m; clabels = l)
@test p isa GLMakie.Figure

@info "Test: plot_imf()"
imf = rand(5, 100)
t = collect(1:100)
p = plot_imf(imf; t = t)
@test p isa GLMakie.Figure

@info "Test: plot_gridlocs()"
p = plot_gridlocs()
@test p isa GLMakie.Figure

@info "Test: plot_fi()"
t = e10.epoch_time
fi = frqinst(e10; ch = "Fp1")
p = plot_fi(fi[1, :, 1], t)
@test p isa GLMakie.Figure

@info "Test: plot_phase()"
s = e10.data[1, :, 1]
X = ftransform(s)
f, _ = NeuroAnalyzer._freqs(s, sr(e10))
p = plot_phase(rad2deg.(X.ph[1:100]), f[1:100]; unit = :rad, type = :stem)
@test p isa GLMakie.Figure
p = plot_phase(rad2deg.(X.ph[1:100]), f[1:100]; unit = :deg, type = :stem)
@test p isa GLMakie.Figure
p = plot_phase(rad2deg.(X.ph[1:100]), f[1:100]; unit = :rad, type = :line)
@test p isa GLMakie.Figure
p = plot_phase(rad2deg.(X.ph[1:100]), f[1:100]; unit = :deg, type = :line)
@test p isa GLMakie.Figure

@info "Test: plot_locs3d()"
p = plot_locs3d(e10; ch = "eeg", mesh_type = :disabled, gui = false)
@test p isa GLMakie.Figure
p = plot_locs3d(e10; ch = "eeg", mesh_type = :brain, gui = false)
@test p isa GLMakie.Figure
p = plot_locs3d(e10; ch = "eeg", mesh_type = :head, gui = false)
@test p isa GLMakie.Figure

@info "Test: plot_polezero()"
p = plot_polezero(rand(ComplexF64, 2), rand(ComplexF64, 2))
@test p isa GLMakie.Figure

@info "Test: plot_dwc()"
dc = rand(16, 100)
t = collect(1:100)
p = plot_dwc(dc; t = t, n = 4)
@test p isa GLMakie.Figure

@info "Test: plot regressions"
# frequency range / log axis in 3-D and topographic PSD and phase plots
@test plot_psd(e10; ep = 1, ch = ["Fp1", "Fp2"], type = :w3d, frq = :log, gui = false) isa GLMakie.Figure
@test plot_psd(e10; ep = 1, ch = "eeg", type = :topo, frq = :log, gui = false) isa GLMakie.Figure
@test plot_phs(e10; ep = 1, ch = ["Fp1", "Fp2"], type = :s3d, frq = :log, gui = false) isa GLMakie.Figure
# relative PSD keeps the plotted range; Morlet wavelets with a band starting below 1 Hz
@test plot_psd(e10; ep = 1, ch = "Fp1", method = :mw, ref = :delta, gui = false) isa GLMakie.Figure
@test plot_spectrogram(e10; ep = 1, ch = "eeg", type = :topo, method = :gh, gui = false) isa GLMakie.Figure
@test NeuroAnalyzer.plot(e10_erp; ch = ["Fp1", "Fp2"], type = :butterfly, gui = false) isa GLMakie.Figure
@test plot_mep(mep; ch = labels(mep)[1], type = :stack) isa GLMakie.Figure
@test plot_topo(e10; ch = "eeg", tpos = 1, sch = r"^Fp") isa GLMakie.Figure
# vector method with the default channels (all rows of locs)
@test plot_topo(rand(nrow(e10.locs)); locs = e10.locs) isa GLMakie.Figure
@test plot_locs(e10; ch = "eeg", head_labels = true, gui = false) isa GLMakie.Figure
@test plot_efield2d([1, -1], [-1.0, 1.0], [0.0, 0.0]) isa GLMakie.Figure
f = collect(0.0:0.5:50.0)
@test plot_coherence(rand(length(f)), f; frq = :log) isa GLMakie.Figure
# logarithmic frequency axes are base 10
@test first(c for c in plot_coherence(rand(length(f)), f; frq = :log).content if c isa GLMakie.Axis).xscale[] === log10
@test plot_eros(rand(41, 21, 2), collect(0.0:40.0), collect(0.0:0.1:2.0); frq = :log) isa GLMakie.Figure
m = rand(19, 19)
@test plot_connectivity_circle(m; clabels = labels(e10)[1:19], threshold = (0.2, 0.8), threshold_type = :in) isa GLMakie.Figure
v = [plot_empty(), plot_empty()]
plot_compose(v; layout = (2, 2))
@test length(v) == 2
for fn in ("test.pdf", "test.PNG")
    isfile(fn) && rm(fn)
    plot_save(plot_empty(); file_name = fn)
    @test isfile(fn)
    rm(fn)
end
@info "Test: ch_locs head inset"
@test plot_psd(e10; ep = 1, ch = ["Fp1", "Cz", "O2"], ch_locs = true, gui = false) isa GLMakie.Figure
@test plot_psd(e10; ep = 1, ch = "T3", ch_locs = true, gui = false) isa GLMakie.Figure
@test plot_phs(e10; ep = 1, ch = ["Fp1", "Cz"], ch_locs = true, gui = false) isa GLMakie.Figure
@test plot_erp(e10_erp; ch = ["F3", "F4", "P3"], ch_locs = true, gui = false) isa GLMakie.Figure
@test plot_erp(e10_erp; ch = "Cz", ch_locs = true, gui = false) isa GLMakie.Figure
@test NeuroAnalyzer.plot(e10_erp; ch = ["F3", "F4"], ci95 = true, avg = false, ch_locs = true, gui = false) isa GLMakie.Figure
# the inset is the second layout cell, next to the plot
p = plot_psd(e10; ep = 1, ch = ["Fp1", "Cz"], ch_locs = true, gui = false)
@test count(c -> c isa GLMakie.Axis, p.content) == 2
@test count(c -> c isa GLMakie.Legend, p.content) == 1
# inside the plot: the inset shares the plot's cell; the legend stays outside
for pos in (:rt, :lt, :rb, :lb)
    local p = plot_psd(e10; ep = 1, ch = ["Fp1", "Cz"], ch_locs = true, ch_locs_pos = pos, gui = false)
    @test count(c -> c isa GLMakie.Axis, p.content) == 2
    @test count(c -> c isa GLMakie.Legend, p.content) == 1
end
@test plot_erp(e10_erp; ch = ["F3", "F4"], ch_locs = true, ch_locs_pos = :lb, gui = false) isa GLMakie.Figure
@test plot_phs(e10; ep = 1, ch = "Cz", ch_locs = true, ch_locs_pos = :lt, gui = false) isa GLMakie.Figure
@test_throws ArgumentError plot_psd(e10; ep = 1, ch = "Cz", ch_locs = true, ch_locs_pos = :middle, gui = false)
@info "Test: legends outside the plot"
nleg(p) = count(c -> c isa GLMakie.Legend, p.content)
f = collect(0.0:1.0:50.0)
p = plot_psd(f, rand(3, length(f)); clabels = ["a", "b", "c"])
@test nleg(p) == 1
# the legend is in the side column, not in the plot
@test any(c -> c isa GLMakie.GridLayout, GLMakie.contents(p[1, 2]))
@test nleg(plot_psd(f, rand(3, length(f)); leg = false)) == 0
# too many channels: no legend
@test nleg(plot_psd(f, rand(80, length(f)))) == 0
# fits in 2 columns
p = plot_psd(f, rand(25, length(f)))
@test nleg(p) == 1 && only(c for c in p.content if c isa GLMakie.Legend).nbanks[] == 2
@test nleg(plot_phs(f, rand(3, length(f)))) == 1
@test nleg(plot_coherence(rand(3, length(f)), f)) == 1
@test nleg(plot_erp(e10_erp; ch = ["F3", "F4"], gui = false)) == 1
@test nleg(plot_mep(mep.epoch_time, mep.data[:, :, 1])) == 1
# legend below an outside head inset
p = plot_psd(e10; ep = 1, ch = ["Fp1", "Cz", "O2"], ch_locs = true, gui = false)
@test nleg(p) == 1 && count(c -> c isa GLMakie.Axis, p.content) == 2

@info "Test: add_inset()"
for pos in (:outside, :rt, :lt, :rb, :lb)
    local p = plot_psd(e10; ep = 1, ch = "Cz", gui = false)
    @test add_inset(p, plot_locs(e10; ch = "eeg", ps = :s, gui = false); pos = pos) === p
    @test count(c -> c isa GLMakie.Axis, p.content) == 2
end
p = plot_psd(e10; ep = 1, ch = "Cz", gui = false)
add_inset(p, plot_locs(e10; ch = "eeg", ps = :m, gui = false); pos = :rb, width = 150)
@test any(c -> c isa GLMakie.Axis && c.width[] == 150 && c.height[] == 150, p.content)
@test_throws ArgumentError add_inset(plot_empty(), plot_empty(); pos = :middle)
@test_throws ArgumentError add_inset(plot_empty(), plot_empty(); width = 0)
# no locations: the plot is drawn without the inset
e_noloc = epoch(eeg; ep_len = 10)
e_noloc.locs = e_noloc.locs[1:0, :]
p = plot_psd(e_noloc; ep = 1, ch = ["Fp1", "Cz"], ch_locs = true, gui = false)
@test count(c -> c isa GLMakie.Axis, p.content) == 1
@test plot_locs(e10; ch = "Cz", gui = false) isa GLMakie.Figure
nirs = import_nirs(joinpath(testfiles_path, "nirs", "fnirs-test-nirs.nirs"))
r = nirs.header.recording
@test plot_locs_nirs(nirs.locs, r[:optode_pairs], length(r[:src_labels]), length(r[:det_labels])) isa GLMakie.Figure


##
## untested plot functions and helpers; content checks (labels, plotted data)
##

axes_of(fig) = [c for c in fig.content if c isa GLMakie.Axis]
# the (x, y) data of every Lines plot in an axis
line_data(ax) = [p[1][] for p in ax.scene.plots if p isa GLMakie.Lines]

@info "Test: plot_cont() / plot_ep()"
p = plot_cont(eeg; ch = ["Fp1", "Fp2"], seg = (0, 5), gui = false)
@test p isa GLMakie.Figure
p = plot_cont(eeg; ch = "Fp1", seg = (0, 5), type = :butterfly, gui = false)
@test p isa GLMakie.Figure
p = plot_cont(eeg, eeg; ch = ["Fp1", "Fp2"], seg = (0, 5), gui = false)
@test p isa GLMakie.Figure
@test_throws ArgumentError plot_cont(eeg; ch = "Fp1", res = 0, gui = false)
@test_throws ArgumentError plot_cont(eeg, e10; ch = "Fp1", gui = false)
p = plot_ep(e10; ch = ["Fp1", "Fp2"], ep = 2, gui = false)
@test p isa GLMakie.Figure
@test_throws ArgumentError plot_ep(e10; ch = "Fp1", ep = 11, gui = false)

@info "Test: plot_gfp() / plot_erp_stack() / plot_erp_topo()"
tg = collect(0.0:0.01:1.0)
g = sin.(2π .* tg) .^ 2
p = plot_gfp(tg, g; title = "GFP test", xlabel = "Time [s]", ylabel = "GFP")
ax = axes_of(p)[1]
@test ax.title[] == "GFP test" && ax.xlabel[] == "Time [s]"
@test any(d -> [q[2] for q in d] ≈ g, line_data(ax))
p = plot_erp_stack(tg, rand(3, length(tg)); clabels = ["a", "b", "c"])
@test p isa GLMakie.Figure
e_erp = average_epochs(e10)
ch_eeg = get_channel(e_erp; ch = get_channel(e_erp; type = "eeg"))
locs_e = e_erp.locs[NeuroAnalyzer._find_bylabel(e_erp.locs, labels(e_erp)[ch_eeg]), :]
p = plot_erp_topo(locs_e, e_erp.epoch_time, e_erp.data[ch_eeg, :, 1]; clabels = labels(e_erp)[ch_eeg], gui = false)
# one trace per channel, placed at the channel locations
@test p isa GLMakie.Figure

@info "Test: plot_mep_stack()"
mep = import_duomag(joinpath(testfiles_path, "mep", "mep-duomag.ascii"))
p = plot_mep_stack(mep.epoch_time, mep.data[:, :, 1])
@test p isa GLMakie.Figure

@info "Test: plot_psd_3d() / plot_psd_topo() / plot_phs_3d() / plot_phs_topo() / plot_spectrogram_topo()"
pw, fw = NeuroAnalyzer._psd(e10.data[ch_eeg, :, 1]; fs = sr(e10))
p = plot_psd_3d(fw, pw; clabels = labels(e10)[ch_eeg])
@test p isa GLMakie.Figure
p = plot_psd_3d(fw, pw; clabels = labels(e10)[ch_eeg], db = false, flim = (0, 40))
@test p isa GLMakie.Figure
p = plot_psd_topo(locs_e, fw, pw; gui = false)
@test p isa GLMakie.Figure
ph, fp = NeuroAnalyzer._phs(e10.data[ch_eeg, :, 1]; fs = sr(e10))
p = plot_phs_3d(fp, ph; clabels = labels(e10)[ch_eeg])
@test p isa GLMakie.Figure
p = plot_phs_topo(locs_e, fp, ph; gui = false)
@test p isa GLMakie.Figure
sp = NeuroAnalyzer._spectrogram(e10.data[ch_eeg, :, 1]; fs = sr(e10), method = :stft, db = false)
p = plot_spectrogram_topo(locs_e, sp.t, sp.f, sp.p; gui = false)
@test p isa GLMakie.Figure

@info "Test: plot_filter() helpers"
f = filter_create(; fprototype = :butterworth, ftype = :lp, cutoff = 30, fs = 256, order = 4)
fr = collect(0.0:1.0:128.0)
H = NeuroAnalyzer._filter_freqz(f, fr, 256)
@test H ≈ NeuroAnalyzer.DSP.freqresp(f, 2π .* fr ./ 256)
@test abs(H[1]) ≈ 1
g1, ph1, gd1 = NeuroAnalyzer._filter_response(f, fr, 256, :onepass)
g2, ph2, gd2 = NeuroAnalyzer._filter_response(f, fr, 256, :twopass)
# two-pass: squared magnitude (dB × 2), zero phase
@test g2 ≈ 2 .* g1
@test all(iszero, ph2)
@test g1[31] ≈ -3.01 atol = 0.02
# FIR (symmetric, N taps): constant group delay (N - 1) / 2
fir = filter_create(; fprototype = :fir, ftype = :lp, cutoff = 30, fs = 256, order = 51)
Hf = NeuroAnalyzer._filter_freqz(fir, fr, 256)
gd = NeuroAnalyzer._group_delay(Hf, fr, 256)
@test all(x -> isnan(x) || isapprox(x, 25; atol = 1e-3), gd[2:20])
@test NeuroAnalyzer._bw_effective(fir, :fir, 256, :hamming, false, nothing) ≈ 3.3 * 256 / 51
@test isnothing(NeuroAnalyzer._bw_effective(f, :butterworth, 256, :hamming, false, nothing))
@test NeuroAnalyzer._bw_effective(fir, :remez, 256, :hamming, false, 4) == 4
@test NeuroAnalyzer._fmt(1.234) == "1.23"
@test NeuroAnalyzer._fmt((1.0, 2.555)) == "(1.0, 2.56)"
@test NeuroAnalyzer._fmt(Float64[]) == "not reached"
ttl = NeuroAnalyzer._filter_title(f; fprototype = :butterworth, ftype = :lp, cutoff = 30, order = 4, bw = nothing,
                                  bw_eff = nothing, rp = nothing, rs = nothing, fs = 256, dir = :twopass,
                                  window = :hamming, custom_w = false, status = "")
@test occursin("Butterworth", ttl) && occursin("LP", ttl) && occursin("order: 4", ttl)
fig = GLMakie.Figure()
ax = NeuroAnalyzer._filter_axis(fig, (1, 1), "t", "Gain [dB]", (0, 128))
@test ax isa GLMakie.Axis && ax.ylabel[] == "Gain [dB]" && ax.xlabel[] == "Frequency [Hz]"

@info "Test: plot helpers (varia, locs)"
@test NeuroAnalyzer._ylims_padded([1.0, 4.0]) == (0, 6)
@test NeuroAnalyzer._ylims_padded([-2.0, 4.0]) == (-3, 6)
@test NeuroAnalyzer._passes_threshold(5, 3, :g)
@test !NeuroAnalyzer._passes_threshold(5, 3, :l)
@test NeuroAnalyzer._passes_threshold(3, 3, :eq)
@test NeuroAnalyzer._passes_threshold(3, 4, :neq)
@test NeuroAnalyzer._passes_threshold(3, 3, :geq) && NeuroAnalyzer._passes_threshold(3, 3, :leq)
@test NeuroAnalyzer._passes_threshold(3, (1, 5), :in) && !NeuroAnalyzer._passes_threshold(5, (1, 5), :bin)
fig = GLMakie.Figure()
ax = GLMakie.Axis(fig[1, 1])
GLMakie.current_axis!(ax)
n0 = length(ax.scene.plots)
NeuroAnalyzer._draw_connection!([0.0, 1.0], [0.0, 1.0], 1, 2, 0.5, 2.0, false, true)
@test length(ax.scene.plots) > n0
n0 = length(ax.scene.plots)
NeuroAnalyzer._draw_connection_label!([0.0, 1.0], [0.0, 1.0], 1, 2, 0.5, 8, false)
@test length(ax.scene.plots) > n0


@info "Test: plot_env()"
r = tenv(e10; ch = "Fp1")
p = plot_env(r.t, r.e[1, :, 1]; title = "env")
@test p isa GLMakie.Figure
ax = axes_of(p)[1]
@test ax.title[] == "env" && ax.xlabel[] == "Time [s]"
@test any(d -> [q[2] for q in d] ≈ r.e[1, :, 1], line_data(ax))
r = tenv_mean(e10; ch = "Fp1", dims = 3)
p = plot_env(r.t, r.em[:, 1]; el = r.el[:, 1], eu = r.eu[:, 1])
# the mean and its 95% CI band
@test count(x -> x isa GLMakie.Band, axes_of(p)[1].scene.plots) == 1
@test_throws ArgumentError plot_env(r.t, r.em[1:10, 1])
@test_throws ArgumentError plot_env(r.t, r.em[:, 1]; el = r.el[:, 1])

@info "Test: plot_ispc() / plot_pli()"
r = ispc(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
p = plot_ispc(e10.epoch_time, e10.data[1, :, 1], e10.data[2, :, 1]; ispcv = r.ispcv[1], ispca = r.ispca[1],
              sd = r.sd[1, :, 1], phd = r.phd[1, :, 1], s1ph = r.s1ph[1, :, 1], s2ph = r.s2ph[1, :, 1], unit = "μV")
@test p isa GLMakie.Figure
# 4 time-course axes (signals, difference, phases, phase difference) + 2 polar axes
@test length(axes_of(p)) == 4
@test count(c -> c isa GLMakie.PolarAxis, p.content) == 2
@test axes_of(p)[1].ylabel[] == "Amplitude [μV]"
@test any(c -> c isa GLMakie.PolarAxis && occursin("ISPC = $(round(r.ispcv[1]; digits = 2))", c.title[]), p.content)
r = pli(e10, e10; ch1 = "Fp1", ch2 = "Fp2", ep1 = 1, ep2 = 1)
p = plot_pli(e10.epoch_time, e10.data[1, :, 1], e10.data[2, :, 1]; pv = r.pv[1], sd = r.sd[1, :, 1],
             phd = r.phd[1, :, 1], s1ph = r.s1ph[1, :, 1], s2ph = r.s2ph[1, :, 1])
@test any(c -> c isa GLMakie.PolarAxis && occursin("PLI = $(round(r.pv[1]; digits = 2))", c.title[]), p.content)
@test_throws ArgumentError plot_pli(e10.epoch_time, e10.data[1, 1:10, 1], e10.data[2, :, 1]; pv = r.pv[1],
                                    sd = r.sd[1, :, 1], phd = r.phd[1, :, 1], s1ph = r.s1ph[1, :, 1], s2ph = r.s2ph[1, :, 1])

@info "Test: plot_itpc() / plot_itpc_spec() / plot_itpc_frq()"
r = itpc(e10; ch = "Fp1", t = 256)
p = plot_itpc(r.itpcph[1, :]; itpcv = r.itpcv[1], itpca = r.itpca[1])
@test p isa GLMakie.Figure
@test any(c -> c isa GLMakie.PolarAxis && occursin("ITPC = $(round(r.itpcv[1]; digits = 2))", c.title[]), p.content)
p = plot_itpc(r.itpcph[1, :]; itpcv = r.itpcz[1], itpca = r.itpca[1], label = "ITPCz")
@test any(c -> c isa GLMakie.PolarAxis && occursin("ITPCz", c.title[]), p.content)
@test_throws ArgumentError plot_itpc(Float64[]; itpcv = 0.5, itpca = 0.0)
r = itpc_spec(e10; ch = "Fp1", flim = (2, 20), nfrq = 10)
p = plot_itpc_spec(e10.epoch_time, r.f, r.itpcs)
@test p isa GLMakie.Figure
p = plot_itpc_frq(e10.epoch_time, r.itpcs[5, :]; title = "ITPC at $(r.f[5]) Hz")
ax = axes_of(p)[1]
@test ax.title[] == "ITPC at $(r.f[5]) Hz"
@test any(d -> [q[2] for q in d] ≈ r.itpcs[5, :], line_data(ax))
@test_throws ArgumentError plot_itpc_frq(e10.epoch_time, r.itpcs[5, 1:10])

true
