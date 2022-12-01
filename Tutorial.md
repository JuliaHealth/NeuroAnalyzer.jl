# Neuro.jl Tutorial

View markers:
```julia
eeg_view_markers(edf)
```

### EEG Analyze

Channels stats:
```julia
eeg_channels_stats(edf)
```

Calculate signal total power:
```julia
eeg_total_power(edf)
```

Calculate band power:
```julia
eeg_band_power(edf, f=(8, 12.5))
```

Calculate mean and maximum band power and frequency of maximum band power:
```julia
_, mfrq, _ = eeg_band_mpower(edf, f=eeg_band(edf, band=:alpha))
eeg_plot_channels(edf, c=mfrq, epoch=1, title="Maximum α band frequency\n[epoch: 1]")
```

Calculate covariance matrix:
```julia
eeg_cov(edf)
```

Calculate correlation matrix
```julia
eeg_cor(edf)
```

Calculate auto-covariance:
```julia
eeg_acov(edf, lag=20, norm=false)
```

Calculate cross-covariance:
```julia
cc, lags = eeg_xcov(edf, lag=20, demean=true)
# channel by channel, all combinations
plot(lags, cc[1, :])

edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=8)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=12)
cc, lags = eeg_xcov(edf1, edf2, channel1=1, channel2=1, epoch1=1, epoch2=1, lag=20, demean=true, norm=true)
plot(lags, cc[1, :])
```

Stationarity:
```julia
p = eeg_stationarity(edf, method=:mean)
p = eeg_stationarity(edf, method=:var)
plot(p[1, :, :], legend=false)
eeg_stationarity!(edf, method=:var)

p = eeg_stationarity(edf, method=:hilbert)
signal_mi(p[1, :, 1], p[2, :, 1])
m = eeg_mi(edf)
eeg_plot_matrix(edf, m)

plot(p[1, :, :], ylims=(-10, 10), legend=false)
p = eeg_stationarity(edf, window=100, method=:euclid)
plot(p[10:end])
```

Entropy:
```julia
e = eeg_entropy(edf)
plot(eeg_labels(edf), e, seriestype=:bar)
```

Coherence over time:
```julia
c, msc, ic = eeg_tcoherence(edf1, edf2)
plot(c[2, 1:2560, 1])
plot(ic[2, 1:2560, 1])

edf_alpha = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(edf, band=:alpha), order=8)
eeg_labels(edf_alpha)
# O1 vs O2
c, ic = eeg_tcoherence(edf_alpha, edf_alpha, channel1=9, channel2=10, epoch1=1, epoch2=1)
plot(c[1:2560])
plot(ic[1:2560])
```

Coherence over frequencies:
```julia
c, msc, f = eeg_fcoherence(edf, edf, channel1=[1, 2], channel2=3:4, epoch1=1, epoch2=1)
plot(f[1, :], c[1, :, 1])

# O1 vs O2, alpha range
c, msc, f = eeg_fcoherence(edf_alpha, edf_alpha, channel1=9, channel2=10, epoch1=1, epoch2=1, frq_lim=eeg_band(edf, band=:alpha))
plot(f[1, :], c[1, :, 1])
```

Comparing two signals:
```julia
edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(edf, band=:delta), order=8)
edf1 = eeg_epochs(edf1, epoch_len=10*eeg_sr(edf1), average=true)
edf1 = eeg_keep_channel(edf1, 4)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(edf, band=:beta), order=8)
edf2 = eeg_epochs(edf2, epoch_len=10*eeg_sr(edf2), average=true)
edf2 = eeg_keep_channel(edf2, 4)
s, ss, p = eeg_difference(edf1, edf2, n=10, method=:absdiff)
s, ss, p = eeg_difference(edf1, edf2, n=10, method=:diff2int)
# the distribution of bootstrapped signal statistic
histogram(s)
# statistic for signal1 and signal2
vline!([ss])
# how many % of bootstrapped signal statistic are > statistic for signal1 and signal2
ss
# if p < alpha → signals are different
p
```

Signal stats:
```julia
m, s, v = eeg_epochs_stats(edf)
bar(v)
```

ISPC:
```julia
e10 = eeg_epochs(edf, epoch_len=10*256)
i, _, _, _, _, _ = eeg_ispc(e10, e10, channel1=1:5, channel2=6:10, epoch1=1, epoch2=1)
```

PLI:
```julia
e10 = eeg_epochs(edf, epoch_len=10*256)
p, _, _, _, _ = eeg_pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)
```

Amplitude Envelope Correlation:
```julia
aec, aec_p = eeg_aec(e10, e10, channel1=3, channel2=4, epoch1=10, epoch2=10)
```

Perform discrete wavelet transform (DWT) and continuous wavelet transform (CWT):
```julia
dwt_c = eeg_dwt(e10, wt=wavelet(WT.haar), type=:sdwt)
cwt_c = eeg_cwt(e10, wt=wavelet(Morlet(π), β=2))
```

Calculate PSD slope of the alpha band:
```julia
f, psd_slope, frq = eeg_psdslope(eeg, f=(8, 14), norm=true, mt=false)
```

### EEG Plots

Plot filter response:
```julia
p = plot_filter_response(fs=eeg_sr(edf), fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=8)
eeg_plot_save(p, file_name="images/butter_bs_45-55_8.png")
```

![BS Butterworth (45, 55), 8](images/butter_bs_45-55_8.png)

Plot band power:
```julia
using DSP
bands = [:delta, :theta, :alpha, :beta]
frq = Vector{Tuple{Float64, Float64}}()
for idx in bands
    push!(frq, eeg_band(edf, band=idx))
end
bp = Vector{Float64}()
for idx in frq
    push!(bp, pow2db.(eeg_band_power(edf, f=idx))[1])
end
p = eeg_plot_stats(edf, bp, epoch=1, channel=1, plot_by=:labels, labels=string.(bands), type=:bar, xlabel="", ylabel="Power [dB]", title="Band powers\n[epoch: 1, channel: 1 ($(eeg_labels(edf)[1]))]")
eeg_plot_save(p, file_name="images/edf_bands.png")
```

![edf bands](images/edf_bands.png)


Plot PSD relative to alpha band power:
```julia
p = eeg_plot_psd(edf, epoch=1, channel=1, ref=:alpha)
eeg_plot_save(p, file_name="images/edf_rel_psd.png")
```
![](images/edf_rel_psd.png)

Plot phase of the convoluted signal:
```julia
mw = generate_morlet(256, 10, 32, complex=true)
e10_tconv = eeg_tconv(e10, kernel=mw)
pt = s_phases(e10_tconv)
p = eeg_plot(e10, pt, epoch=1, c_idx=1:4, scale=false, emarkers=false)
eeg_plot_save(p, file_name="images/e10_tconv_phases.png")
```
![](images/e10_tconv_phases.png)

Topographical plot:
```julia
p = eeg_plot_topo(edf, segment=(1, 2560))
eeg_plot_save(p, file_name="images/edf_amp.png")
```

![edf topo :amp](images/edf_amp.png)

Topographical plots:
```
p1 = eeg_plot_topo(edf, segment=(1, 2560), title="0:1s", cb=false)
p2 = eeg_plot_topo(edf, segment=(1 * 256 + 1, 2 * 2560), title="1:2s", cb=false)
p3 = eeg_plot_topo(edf, segment=(2 * 256 + 1, 3 * 2560), title="2:3s", cb=false)
p4 = eeg_plot_topo(edf, segment=(3 * 256 + 1, 4 * 2560), title="3:4s", cb=false)
p = plot(p1, p2, p3, p4, layout=(2, 2))
eeg_plot_save(p, file_name="images/edf_topos.png")
```

![](images/edf_topos.png)

Topographical map of PSD
```julia
p = eeg_plot_psd(edf, channel=1:19, epoch=1, frq_lim=(0, 8), mt=true, variant=:topo)
eeg_plot_save(p, file_name="images/edf_psd_topo.png")
```
![](images/edf_psd_topo.png)

Plot PCA components:
```julia
pc, pc_m, pc_var = eeg_pca(edf, n=10)
p = eeg_plot(edf, channel=1:5, segment=(10*eeg_sr(edf), 20*eeg_sr(edf)))
eeg_plot_save(p, file_name="images/edf_1_5.png")
p = eeg_plot(edf, pc, c_idx=1:5, segment=(10*eeg_sr(edf), 20*eeg_sr(edf)))
eeg_plot_save(p, file_name="images/edf_pca_1_5.png")
bar(vec(pc_var))
```

![pca](images/edf_1_5.png)
![pca](images/edf_pca_1_5.png)

Plot ICA components:
```julia
ic, icm = eeg_ica(edf, n=16, tol=0.99)
p = eeg_plot(edf, ic, c_idx=1:10, epoch=1)
eeg_plot_save(p, file_name="images/edf_ica_1_10.png")
```

![edf amplitude :ica](images/edf_ica_1_10.png)

```julia
ic, icm = eeg_ica(edf, n=16, tol=0.99)
eeg_add_component!(edf, c=:ica, v=ic)
eeg_add_component!(edf, c=:ica_mw, v=icm)
# reconstruct signal using ICA 1:10
s_reconstructed = s_ica_reconstruct(edf.eeg_signals, ic=ic, ic_mw=icm, ic_v=1:10)
p1 = eeg_plot_topo(edf, epoch=1, title="Original signal", cb=false)
p2 = eeg_plot_topo(edf, s_reconstructed, epoch=1, title="Signal reconstructed from ICA 1:10", cb=false)
p = plot(p1, p2)
eeg_plot_save(p, file_name="images/edf_ica_1_10.png")
```

![edf topo :ica](images/edf_ica_1_10.png)

Plot alpha band power:
```julia
alpha_power = eeg_band_power(edf, f=eeg_band(edf, band=:alpha))
p = eeg_plot_topo(edf, alpha_power)
eeg_plot_save(p, file_name="images/edf_alpha_topo.png")
```

![edf topo :power](images/edf_alpha_topo.png)

Plot phase difference at time = 1s (sample = 256 as fs = 256 Hz)
```julia
pdiff = eeg_phdiff(edf)
p = eeg_plot_topo(edf, pdiff[:, 256, :])
eeg_plot_save(p, file_name="images/edf_phdiff_topo.png")
```
![](images/edf_phdiff_topo.png)

Plot amp difference at time = 1s (sample = 256 as fs = 256 Hz)
```julia
ampdiff = eeg_ampdiff(edf)
p = eeg_plot_topo(edf, ampdiff[:, 256, :], epoch=1)
eeg_plot_save(p, file_name="images/edf_ampdiff_topo.png")
```
![](images/edf_ampdiff_topo.png)

Plot weights:
```julia
w = (1:19) * 0.05
p = plot_weights(edf.eeg_locs, channel=1:19, weights=round.(w, digits=2), head_labels=false)
eeg_plot_save(p, file_name="images/edf_weights.png")
```

![](images/edf_weights.png)

Plot covariance matrix:
```julia
edf_cov = eeg_cov(edf)
p = plot_matrix(edf_cov[:, :, 1], title="Covariance matrix", labels_x=eeg_labels(edf)[1:19], labels_y=eeg_labels(edf)[1:19])
eeg_plot_save(p, file_name="images/edf_cov.png")
```

![edf cov](images/edf_cov.png)

Plot auto-covariance matrix:
```julia
ac, lags = eeg_acov(edf, lag=5, norm=false)
p = plot_covmatrix(ac[1, :, 1], lags, xlabel="", title=eeg_labels(edf)[eeg_channel_idx(edf, type=Symbol(edf.eeg_header[:signal_type]))][1])
eeg_plot_save(p, file_name="images/edf_autocov.png")
```
![edf autocov](images/edf_autocov.png)

Plot cross-covariance matrix:
```julia
xc, lags = eeg_xcov(edf, lag=5, norm=false)
p = plot_covmatrix(xc[2, :, 1], lags, xlabel="", title="Fp1-Fp2")
eeg_plot_save(p, file_name="images/edf_xcov.png")
```

![edf xcov](images/edf_xcov.png)

Plot channels stats:
```julia
c = eeg_channels_stats(edf)
e = eeg_epochs_stats(edf)
eeg_add_component!(edf, c=:channels_var, v=c[4])
eeg_add_component!(edf, c=:epochs_var, v=e[4])
p = eeg_plot_stats(edf, :channels_var, epoch=1, title="Channels variance\n[epoch: 1]", plot_by=:channels, type=:line)
eeg_plot_save(p, file_name="images/edf_channels.png")
p = eeg_plot_stats(edf, :epochs_var, epoch=1:10, title="Epochs 1:10 variance", plot_by=:epochs, type=:line)
eeg_plot_save(p, file_name="images/e10_epochs.png")
```

![e10 channels](images/edf_channels.png)

![e10 epochs](images/e10_epochs.png)

Envelopes:
```julia
p = eeg_plot_signal_psd(e10, epoch=1, channel=1)
eeg_plot_save(p, file_name="images/e10_psd.png")

p = eeg_plot_env(e10, type=:pow, average=:mean, dims=1, epoch=1, channel=1)
eeg_plot_save(p, file_name="images/e10_penv.png")

p = eeg_plot_signal_spectrogram(e10, epoch=1, channel=1, frq_lim=(0,10))
eeg_plot_save(p, file_name="images/e10_spec.png")

p = eeg_plot_env(e10, type=:spec, average=:median, epoch=1, channel=1, dims=3, frq_lim=(0,10))
eeg_plot_save(p, file_name="images/e10_senv.png")

p = eeg_plot_env(e10, type=:hamp, average=:median, epoch=1, channel=1, dims=3)
eeg_plot_save(p, file_name="images/e10_henv.png")
```

![e10 power](images/e10_psd.png)

![e10 power envelope mean](images/e10_penv.png)

![e10 spectrogram](images/e10_spec.png)

![e10 spectrogram envelope mean](images/e10_senv.png)

![e10 Hilbert spectrum amplitude median](images/e10_henv.png)

ISPC:
```julia
p = eeg_plot_ispc(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)
eeg_plot_save(p, file_name="images/e10_ispc.png")
m = eeg_ispc(e10)
p = plot_matrix(m[:, :, 1], labels_x=eeg_labels(e10), labels_y=eeg_labels(e10))
eeg_plot_save(p, file_name="images/e10_ispc_m.png")
p = eeg_plot_connections(e10, m=m[:, :, 1], threshold=0.90, threshold_type=:geq)
eeg_plot_save(p, file_name="images/e10_ispc_connections.png")
```

![e10 ISPC](images/e10_ispc.png)

![e10 ISPC matrix](images/e10_ispc_m.png)

![e10 ISPC matrix](images/e10_ispc_connections.png)

ITPC:
```julia
p = eeg_plot_itpc(e10, channel=1, t=256)
eeg_plot_save(p, file_name="images/e10_itpc.png")
p = eeg_plot_itpc_s(e10, channel=1, frq_lim=(1, 20), frq_n=20)
eeg_plot_save(p, file_name="images/e10_itpc_s.png")
i, f = eeg_itpc_s(e10, channel=1, frq_lim=(1, 20), frq_n=20)
# plot ITCP at 4 Hz frequency over epoch time
eeg_plot_itpc_f(e10, channel=1, frq_lim=(0, 10), frq_n=10, f=4, frq=:lin)
```

![e10 ITPC](images/e10_itpc.png)

![e10 ITPC spectrogram](images/e10_itpc_s.png)

PLI:
```julia
p = eeg_pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1)
p1 = eeg_plot(e10, segment=(1, 256), channel=1:2, variant=:butterfly, title="Signals")
p2 = eeg_plot(e10, p.signal_diff, segment=(1, 256), title="Signals difference [μV]", scale=false)
p3 = eeg_plot(e10, [p.s1_phase; p.s2_phase], segment=(1, 256), channel=1:2, variant=:butterfly, title="Phases [rad]")
p4 = eeg_plot(e10, p.phase_dif, segment=(1, 256), channel=1:2, title="Phases difference [rad]", scale=false)
p5 = eeg_plot_stats(e10, hcat(p.s1_phase[1, :, 1], ones(length(p.s1_phase[:, :, 1]))), plot_by=:channels, channel=1:2, type=:polar, title="Phases")
p6 = eeg_plot_stats(e10, hcat(p.phase_dif[1, :, 1], ones(length(p.s1_phase[:, :, 1]))), plot_by=:channels, channel=1:2, type=:polar, title="Phases difference and PLI = $(p.pli[1])")
p = plot(p1, p2, p3, p4, p5, p6, layout=(3, 2))
eeg_plot_save(p, file_name="images/e10_pli.png")

m = eeg_pli(e10)
p = plot_matrix(m[:, :, 1], labels_x=eeg_labels(e10), labels_y=eeg_labels(e10))
eeg_plot_save(p, file_name="images/e10_pli_m.png")
```

![e10 PLI](images/e10_pli.png)

![e10 PLI matrix](images/e10_pli_m.png)

Connections based on Hilbert transform amplitude envelope:
```julia
h, t = eeg_henv(edf)
m = zeros(size(h, 1), size(h, 1))
for idx1 in 1:size(h, 1)
    for idx2 in 1:size(h, 1)
        c = s2_cor(h[idx1, :, 1], h[idx2, :, 1])
        m[idx1, idx2] = c.r
    end
end
p = eeg_plot_connections(edf, m=m, threshold=0.2, threshold_type=:geq)
eeg_plot_save(p, file_name="images/h_connections.png")
```
![](images/h_connections.png)

Animate:
```julia
anim = @animate for i ∈ 1::10:2*eeg_sr(edf)
    eeg_plot_topo(edf, segment=(2560 + i, 2560 + i + 1), channel=1:19)
end
gif(anim, "/tmp/anim_fps15.gif", fps = 15)
```

Using kwargs:
```julia
c = rand(-pi:0.01:pi, 100)
p = eeg_plot_stats(edf, c, plot_by=:epochs, epoch=1, type=:hist)
p = plot!(title="Phases in radians", xticks=([-pi, pi], ["-π", "π"]))
eeg_plot_save(p, file_name="images/kwargs.png")
```

![](images/kwargs.png)

### Statistics

Generate spectrogram segments:
```julia
sp, sf, st = eeg_spectrogram(e10) 
segp1, segs1, tidx1, fidx1 = s_specseg(sp, st, sf, channel=1, t=(1.0, 4.0), f=(10.0, 20.0))
segp2, segs2, tidx2, fidx2 = s_specseg(sp, st, sf, channel=1, t=(5.0, 8.0), f=(45.0, 55.0))

p = eeg_plot_spectrogram(e10, channel=1, epoch=1)
p = plot!(segs1, lc=:black, fill=nothing, label=false)
p = plot!(segs2, lc=:white, fill=nothing, label=false)
eeg_plot_save(p, file_name="images/spec_seg.png")

segp1 = seg_mean(segp1)
segp2 = seg_mean(segp2)
tt, t, c, df, p = s2_cmp(segp1, segp2, paired=true, type=:p);
println("segment 1: mean $(round(mean(segp1), digits=2)), sd $(round(std(segp1), digits=2))")
println("segment 2: mean $(round(mean(segp2), digits=2)), sd $(round(std(segp2), digits=2))")
println("test statistic $(t[2]): $(t[1]) (df = $df), p: $p")
p = boxplot([segp1, segp1], xticks=([1, 2], ["segment 1", "segment 2"]), legend=false, outliers=false)
segp = _labeled_matrix2dict(["segment 1", "segment 2"], [segp1, segp2])
p = eeg_plot_stats(e10, segp, plot_by=:labels, type=:box, title="Mean power [dB]", outliers=false)
eeg_plot_save(p, file_name="images/spec_seg_box.png")
```
![](images/spec_seg.png)
![](images/spec_seg_box.png)

### Pipelines

Pipelines are stored in `~/NeuroAnalyzer/pipelines`.

```julia
include("$(homedir())/NeuroAnalyzer/pipelines/test_pipeline.jl")
edf2 = test_pipeline(edf)
```

### Study

Create study object:
```julia
edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=4)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=8)
edf3 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=12)
my_study = eeg_study_create([edf1, edf2, edf3], [:g1, :g2, :g3])
my_study.study_eeg[1].eeg_signals
```

### EEG Misc

hz, nyq = eeg_freqs(eeg)