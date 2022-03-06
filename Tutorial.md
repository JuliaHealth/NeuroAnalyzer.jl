# Neuro.jl Tutorial

## EEG

For interactive processing it's best to use Pluto:
```julia
import Pluto
Pluto.run()
```

Load package:
```julia
using NeuroJ
```

Show version (for reproducibility):
```julia
neuroj_version()
```

Get help:
```julia
?eeg_plot
```

### EEG IO

Load EDF file:
```julia
edf = eeg_import_edf("test/eeg-test-edf.edf")
```

Load electrode positions:
```julia
eeg_load_electrodes!(edf, file_name="locs/standard-10-20-cap19-elmiko.ced")
p = eeg_plot_electrodes(edf, labels=true, head=true, size=(300, 300))
eeg_plot_save(p, file_name="images/edf_electrodes.png")
```

![EDF electrodes](images/edf_electrodes.png)


Save EEG object:
```julia
eeg_save(edf, file_name="test.bin")
eeg_save(edf, file_name="test.bin", overwrite=true)
eeg_info(edf)
```

Load EEG object:
```julia
edf = eeg_load("test.bin")
```

Export EEG data and header to .csv:
```julia
eeg_export_csv(edf, file_name="edf.csv", header=true)
```

### EEG Edit

Show EEG object properties:
```julia
eeg_info(edf)
eeg_sr(edf)
eeg_channel_n(edf)
eeg_epoch_n(edf)
eeg_signal_len(edf)
eeg_epoch_len(edf)
edf.eeg_header[:eeg_duration_seconds]
```

Edit EEG header
```julia
eeg_show_header(edf)
eeg_edit_header!(edf, field=:patient, value="N.N.")
```

Show labels:
```julia
eeg_labels(edf)
```

Get channel index:
```julia
eeg_get_channel(edf, channel="Cz")
eeg_get_channel(edf, channel=18)
```

Rename channel:
```julia
edf = eeg_rename_channel(edf, channel="Cz", new_name="CZ")
eeg_rename_channel!(edf, channel=18, new_name="Cz")
```

Show processing history:
```julia
eeg_history(edf)
```

### EEG Process

Show components (e.g. ICA, PCA):
```julia
eeg_list_components(edf)
```

Any action that changes EEG signal data (e.g. channel removal, filtering) resets embedded components.

Get component:
```julia
eeg_extract_component(edf, c=:epochs_mean)
```

Delete component:
```julia
edf = eeg_delete_component(edf, c=:ica)
eeg_delete_component!(edf, c=:ica)
```

Upsample:
```julia
eeg_sr(edf)
eeg_info(edf)
edf_512 = eeg_upsample(edf, new_sr=512)
eeg_info(edf_512)
eeg_upsample!(edf, new_sr=1024)
```

Downsample:
```julia
eeg_sr(edf)
edf_128 = eeg_downsample(edf, new_sr=128)
eeg_info(edf_128)
eeg_downsample!(edf, new_sr=64)
```

Remove parts of the signal:
```julia
edf = eeg_trim(edf, len=(10*eeg_sr(edf)), from=:start)
edf = eeg_trim(edf, len=(10*eeg_sr(edf)), offset=(10*eeg_sr(edf)), from=:start)
eeg_trim!(edf, len=(10*eeg_sr(edf)), from=:start)
```

Split into 10-second epochs:
```julia
e10 = eeg_epochs(edf, epoch_len=10*eeg_sr(edf))
eeg_info(e10)
eeg_plot(e10)
eeg_plot(e10, len=10*eeg_sr(edf))
eeg_plot(e10, len=10*eeg_sr(edf), offset=12*256)
eeg_plot(edf)
eeg_epochs!(edf, epoch_len=10*eeg_sr(edf))
```

Trim 1 second from each epoch:
```julia
e9 = eeg_trim(e10, trim_len=(1*eeg_sr(e10)), from=:start)
eeg_info(e9)
eeg_plot(e9, len=11*eeg_sr(e9))
eeg_plot(e9, len=2*eeg_sr(e9), offset=5*256)
eeg_plot(e9, len=60*eeg_sr(e9), offset=0)
```

Get 1st epoch:
```julia
e10e1 = eeg_extract_epoch(e10, epoch=1)
eeg_info(e10e1)
```

Remove epochs:
```julia
e = eeg_delete_epoch(e10, epoch=8:10)
eeg_delete_epoch!(e10, epoch=8:10)
```

Keep epochs:
```julia
e1 = eeg_keep_epoch(e, epoch=[1, 3, 5, 9])
eeg_keep_epoch!(e, epoch=[1, 3, 5, 9])
```

Split into 5-second averaged epoch
```julia
e2avg = eeg_epochs(edf, epoch_len=5*eeg_sr(edf), average=true)
eeg_info(e2avg)
```

Re-reference
```julia
edf = eeg_reference_channel(edf, channel=[1, 2])
edf = eeg_reference_channel(edf, channel=2:4)
edf = eeg_reference_channel(edf, channel=18)
edf = eeg_reference_car(edf)
eeg_reference_car!(edf)
```

Remove channel:
```julia
edf = eeg_delete_channel(edf, channel=1)
edf = eeg_delete_channel(edf, channel=10:18)
eeg_delete_channel!(edf, channel=10:18)
```

Keep channel:
```julia
edf = eeg_keep_channel(edf, channel=1:4)
eeg_keep_channel!(edf, channel=1:4)
```

Filter response:
```julia
p = eeg_plot_filter_response(edf, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=8)
eeg_plot_save(p, file_name="images/butter_bs_45-55_8.png")
```

![BS Butterworth (45, 55), 8](images/butter_bs_45-55_8.png)

FIR filtering:
```julia
eeg_filter!(edf, fprototype=:fir, ftype=:bs, cutoff=(45, 55), order=8, window=generate_hanning(128))
eeg_filter!(edf, fprototype=:fir, ftype=:lp, cutoff=45.0, order=8, window=generate_hanning(128))
eeg_filter!(edf, fprototype=:fir, ftype=:hp, cutoff=0.1, order=8, window=generate_hanning(128))
```

IIR filtering:
```julia
eeg_filter!(edf, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8)
eeg_filter!(edf, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)
eeg_filter!(edf, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=8)
eeg_filter!(edf, fprototype=:elliptic, ftype=:bs, cutoff=(45, 55), rs=10, rp=1, order=12)
eeg_filter!(edf, fprototype=:chebyshev1, ftype=:bs, cutoff=(45, 55), rs=1, order=8)
eeg_filter!(edf, fprototype=:chebyshev2, ftype=:bs, cutoff=(45, 55), rp=1, order=8)
```

Plot:
```julia
p = eeg_plot(edf)
eeg_plot_save(p, file_name="images/edf_channels.png")
```

![edf_channels](images/edf_channels.png)

```julia
p = eeg_plot(edf, channel=1, frq_lim=(0, 20))
eeg_plot_save(p, file_name="images/edf_channel_1.png")
```

![edf_channel1](images/edf_channel_1.png)

```julia
eeg_plot(edf, len=5*256)
eeg_plot(edf, channel=1:4)
eeg_plot(edf, offset=20*eeg_sr(edf), len=20*eeg_sr(edf))
eeg_plot(edf, norm=false)
p = eeg_plot(e9, head=true)
eeg_save_plot(p, file_name="/tmp/e9.png")
```

```juia
p = eeg_plot_avg(edf, frq_lim=(0, 20), channel=1:4)
eeg_plot_save(p, file_name="images/edf_avg.png")
```

![edf_avg](images/edf_avg.png)

```julia
eeg_plot_avg(edf, channel=1:4)
eeg_plot_avg(edf)
eeg_plot_avg(e9)
eeg_plot_avg(e10, len=125*eeg_sr(edf))
eeg_plot_avg(e10, epoch=1:5)
eeg_plot_avg(e9, len=5*eeg_sr(edf), offset=6*eeg_sr(e9))

eeg_plot_butterfly(edf)
eeg_plot_butterfly(edf, offset=20*256, len=120*eeg_sr(edf), channel=1:4, norm=true)
eeg_plot_butterfly(edf, channel=1:4, norm=true)
eeg_plot_butterfly(e10, epoch=1:5, channel=1:4, norm=true)
p = eeg_plot(edf)
eeg_save_plot(p, file_name="/tmp/edf.pdf")

e9 = eeg_load_electrodes(e9, file_name="locs/standard-10-20-cap19-elmiko.ced")
eeg_plot_butterfly(e9, len=9*eeg_sr(edf))
eeg_plot_butterfly(e9, len=55*eeg_sr(edf))
eeg_plot_butterfly(e9, len=55*eeg_sr(edf), head=true)
eeg_plot_butterfly(e10, len=55*eeg_sr(edf), head=true)
eeg_plot_butterfly(e10)
eeg_plot_butterfly(e10, len=11*eeg_sr(edf), head=true)

# use kwargs
p1 = eeg_plot(edf, title="edf1")
p2 = eeg_plot(edf)
# or modify plots
p2 = plot!(title="edf2")
plot(p1, p2, layout=(1, 2))
```

Calculate covariance matrix:
```julia
edf_cov = eeg_cov(edf)
p = eeg_plot_matrix(edf, edf_cov)
eeg_plot_save(p, file_name="images/edf_cov.png")
eeg_cov!(edf)
```

![edf cov](images/edf_cov.png)

Calculate correlation matrix
```julia
edf_cor = eeg_cor(edf)
eeg_plot_matrix(edf, edf_cor)
eeg_cor!(edf)
```

Calculate auto-covariance:
```julia
ac, lags = eeg_autocov(edf, lag=20, norm=false)
eeg_plot_covmatrix(edf, ac, lags)
eeg_autocov!(edf, lag=20, norm=false)
```

Calculate cross-covariance:
```julia
cc, lags = eeg_crosscov(edf, lag=20, demean=true)
# channel by channel, all combinations
plot(lags, cc[1, :])

edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=8)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=12)
cc, lags = eeg_crosscov(edf1, edf2, lag=20, demean=true, norm=true)
plot(lags, cc[1, :])
```

Normalize:
```julia
eeg_norm_zscore(edf)
eeg_norm_minmax(edf)
eeg_norm_zscore!(edf)
eeg_norm_minmax!(edf)
```

Detect bad epochs:
```julia
bad_epochs = eeg_detect_bad_epochs(e10, method=[:flat, :cor], r=0.7)
eeg_check_bad_epochs(e10, bad_epochs)
eeg_delete_epoch!(e10, epoch=bad_epochs)
```

Remove DC:
```julia
eeg_demean(edf)
eeg_demean!(edf)
```

Taper:
```julia
h = hann(edf.eeg_header[:epoch_duration_samples])
eeg_taper(edf, taper=h)
eeg_taper!(edf, taper=h)
```

Calculate signal derivative:
```julia
eeg_derivative(edf)
eeg_derivative!(edf)
```

Detrend:
```julia
eeg_detrend(edf, type=:linear)
eeg_detrend!(edf, type=:constant)
```

Calculate signal total power:
```julia
tbp = eeg_total_power(edf)
bar(eeg_labels(edf), tbp, xticks=(1:length(eeg_labels(edf)), eeg_labels(edf)))
eeg_total_power!(edf)
```

Calculate signal band power:
```julia
abp = eeg_band_power(edf, f=(8, 12.5))
p = eeg_plot_bands(edf, channel=1, type=:abs)
eeg_plot_save(p, file_name="images/edf_bands.png")
```

![edf bands](images/edf_bands.png)

Get channel index/name:
```julia
f3 = eeg_get_channel(edf, channel="F3")
f4 = eeg_get_channel(edf, channel=4)
f3_f = signal_filter(f3, fprototype=:butterworth, ftype=:hp, cutoff=0.1, fs=eeg_sr(edf), order=8)
```

Time-domain convolution:
```julia
mw = generate_morlet(256, 1, 32, complex=true)
eeg_tconv(e10, kernel=mw)
eeg_tconv!(e10, kernel=mw)
```

Frequency-domain convolution:
```julia
mw = generate_morlet(256, 1, 32, complex=true)
eeg_fconv(e10, kernel=mw)
eeg_fconv!(e10, kernel=mw)
```

### EEG Plots

Spectral analysis:
```julia
edf_pow, edf_frq = eeg_psd(edf, norm=true)
plot(edf_frq[10, :], edf_pow[10, :])
eeg_plot_psd(edf, frq_lim=20.0, channel=1:4)
eeg_plot_psd(edf, norm=true, average=false, frq_lim=50)
eeg_plot_psd(edf, norm=true, average=true, frq_lim=20)
eeg_plot_psd(edf, channel=1:4, average=true)
f3 = eeg_extract_channel(edf, channel="F3")
f4 = eeg_extract_channel(edf, channel=4)
signal_psd(f4, fs=256)
signal_plot_psd(f3, fs=256)
signal_plot_psd(f4, fs=256)

eeg_psd!(edf, norm=true)

eeg_plot_spectrogram(edf, channel=9, norm=true)
eeg_plot_spectrogram(e10, channel=9, norm=true, len=110*eeg_sr(edf))
eeg_plot_spectrogram(e9, channel=9, norm=true, ylim=80, len=75*eeg_sr(edf))
eeg_plot_spectrogram(e9, channel=9, norm=true, ylim=40, len=80*eeg_sr(edf), offset=18*256)
eeg_plot_spectrogram(e9, channel=9, norm=true, ylim=40)
eeg_plot_psd(e10)
eeg_plot_psd(e10, len=120*eeg_sr(edf))
eeg_plot_psd(edf, len=240*eeg_sr(edf), channel=4, frq_lim=20)
```

Topographical plots:
```julia
p = eeg_plot_topo(edf, offset=1, len=2560, frq_lim=(0, 20))
eeg_plot_save(p, file_name="images/edf_amp.png")
```

![edf_topo :amp](images/edf_amp.png)

```julia
eeg_plot_topo(edf, offset=2560, c=:ica, c_idx=1:8)
eeg_plot_topo(edf, offset=2560, len=2560, c=:pca)
eeg_plot_topo(edf, offset=2560, len=2560, c=:power, c_idx=eeg_band(edf, band=:alpha))
```

# Stationarity

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

# Entropy

```julia
e = eeg_entropy(edf)
plot(eeg_labels(edf), e, seriestype=:bar)
eeg_entropy!(edf)
```

# Coherence

```julia
m = eeg_coherence(edf1, edf2)
hz, nyq = eeg_freqs(edf1)
plot(hz, abs.(m[1, 1:length(hz)]), xlims=(0, 40))

edf_alpha = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(edf, band=:alpha), order=8)
eeg_labels(edf_alpha)
# O1 vs O2
m = eeg_coherence(edf_alpha, channel1=9, channel2=10)
plot(hz, abs.(m[1:length(hz)]), xlims=(0, 20))
# Fp1 vs Fp2
m = eeg_coherence(edf_alpha, channel1=1, channel2=2)
plot!(hz, abs.(m[1:length(hz)]), xlims=(0, 20))
```

### EEG PCA

```julia
edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(edf, band=:beta), order=8)
edf1 = eeg_epochs(edf1, epoch_len=10*eeg_sr(edf1), average=true)
edf1 = eeg_keep_channel(edf1, 3)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(edf, band=:beta), order=8)
edf2 = eeg_epochs(edf2, epoch_len=10*eeg_sr(edf2), average=true)
edf2 = eeg_keep_channel(edf2, 4)
pc, pc_var = eeg_pca(edf1, n=4)
plot(pc[1, :, 1, 1])
bar(vec(pc_var))

pc, pc_var = eeg_pca(edf, n=4)
plot(pc[1, :, 1, 1])
bar(vec(pc_var))
eeg_pca!(edf, n=4)
```

### EEG ICA

Generate ICAs:
```julia
i = eeg_ica(edf, n=15, tol=1.0)
eeg_ica!(edf, n=5, tol=1.0, iter=1000, f=:gauss)
```

Plot ICAs:
```julia
eeg_plot_ica(edf, ica=1:5)
eeg_plot_ica(edf, ica=1)
```

Remove ICA #001 component from the signal:
```julia
eeg_ica_reconstruct!(edf, ica=1)
```

Remove ICA #001-007 component from the signal:
```julia
eeg_ica_reconstruct!(edf, ica=1:7)
```

Remove ICA #001, 003 and 007 component from the signal:
```julia
eeg_ica_reconstruct!(edf, ica=[1, 3, 7])
```

### EEG Comparing two signals

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

### EEG Misc

```julia
eeg_band(edf, band=:alpha)
edf = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(edf, band=:alpha), order=8)
hz, nyq = eeg_freqs(edf)

e = eeg_pick(edf, pick=:left)
eeg_labels(edf)[e]
e = eeg_pick(edf, pick=[:l, :f, :t])
eeg_labels(edf)[e]
```

Convert samples to seconds/seconds to samples:
```julia
eeg_s2t(edf, t=1234)
eeg_t2s(edf, t=10)
```

Mutators store results as embedded components in the EEG object:
```julia
eeg_autocov!(edf)
eeg_cor!(edf)
eeg_cov!(edf)
eeg_entropy!(edf)
eeg_freqs!(edf)
eeg_mi!(edf)
eeg_pca!(edf, n=2)
eeg_psd!(edf)
eeg_stationarity!(edf)
eeg_total_power!(edf)

eeg_ica!(edf, n=4)
eeg_component(edf, c=:ica)

eeg_epochs_stats!(edf)
eeg_component(edf, c=:epochs_mean)

eeg_info(edf)
```

### EEG Benchmarking

```julia
edf = eeg_import_edf("test/eeg-test-edf.edf");
function eeg_benchmark(n::Int64)
    for idx in 1:n
        edf_new = eeg_reference_car(edf)
        e10 = eeg_epochs(edf_new, epoch_len=10*eeg_sr(edf_new))
        e10 = eeg_filter(e10, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8)
        e10 = eeg_filter(e10, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)
        e10 = eeg_filter(e10, fprototype=:butterworth, ftype=:bs, cutoff=(45, 55), order=8)
        tbp = eeg_total_power(e10)
        ac = eeg_autocov(e10, norm=false)
        cc = eeg_crosscov(e10, lag=10, demean=true)
        mconv = eeg_tconv(e10, kernel=generate_morlet(256, 1, 32, complex=true))
    end
end
@time eeg_benchmark(10)

# workstation: 24.408695 seconds (55.97 M allocations: 274.946 GiB, 27.82% gc time)
# laptop: 40.043600 seconds (55.26 M allocations: 274.877 GiB, 17.69% gc time, 0.04% compilation time)
```