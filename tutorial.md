# Neuro.jl

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

Load EDF file:
```julia
edf = eeg_import_edf("test/eeg-test-edf.edf")
```

Show EEG object properties:
```julia
eeg_info(edf)
edf.eeg_header[:sampling_rate][1]
edf.eeg_header[:eeg_duration_seconds]
edf.eeg_header[:eeg_duration_samples]
eeg_samplingrate(edf)
```

Show labels:
```julia
eeg_labels(edf)
```

Save EEG object:
```julia
eeg_save(edf, "test.bin")
```

Load EEG object:
```julia
edf = eeg_load("test.bin")
```

Get channel index:
```julia
eeg_get_channel_idx(edf, "Cz")
eeg_get_channel_name(edf, 18)
```

Rename channel:
```julia
edf = eeg_rename_channel(edf, "Cz", "CZ")
edf = eeg_rename_channel(edf, 18, "Cz")
```

Show processing history:
```julia
eeg_history(edf)
```

Upsample:
```julia
eeg_samplingrate(edf)
eeg_info(edf)
edf_512 = eeg_upsample(edf, new_sr=512)
eeg_info(edf_512)
```

Downsample:
```julia
eeg_samplingrate(edf)
edf_128 = eeg_downsample(edf, new_sr=128)
eeg_info(edf_128)
```

Remove parts of the signal:
```julia
edf = eeg_trim(edf, trim_len=(20 * eeg_samplingrate(edf)), from=:start)
edf = eeg_trim(edf, trim_len=(10 * eeg_samplingrate(edf)), offset=(10 * eeg_samplingrate(edf)), from=:start)
```

Split into 10-second epochs:
```julia
e10 = eeg_epochs(edf, epoch_len=10*eeg_samplingrate(edf))
eeg_info(e10)
eeg_plot(e10)
eeg_plot(e10, len=10)
eeg_plot(e10, len=10, offset=12*256)
eeg_plot(edf)
```

Trim 1 second from each epoch:
```julia
e9 = eeg_trim(e10, trim_len=(1 * eeg_samplingrate(edf)), from=:start)
eeg_info(e9)
eeg_plot(e9, len=11)
eeg_plot(e9, len=2, offset=5*256)
eeg_plot(e9, len=60, offset=0)
```

Get 1st epoch:
```julia
e10e1 = eeg_extract_epoch(e10, 1)
eeg_info(e10e1)
```

Split into 5-second averaged epoch
```julia
e2avg = eeg_epochs(edf, epoch_len=5*eeg_samplingrate(edf), average=true)
eeg_info(e2avg)
```

Re-reference
```julia
edf = eeg_reference_channel(edf, [1, 2])
edf = eeg_reference_channel(edf, 2:4)
edf = eeg_reference_channel(edf, 18)
edf = eeg_reference_car(edf)
```

Remove channel:
```julia
eeg_delete_channel(edf, 10:18)
```

Keep channel:
```julia
eeg_keep_channel(edf, 1:4)
```

Show filter response:
```julia
filter_response(fprototype=:butterworth, ftype=:hp, cutoff=10, fs=eeg_samplingrate(edf), order=8)
filter_response(fprototype=:chebyshev1, ftype=:bp, cutoff=[40, 50], fs=eeg_samplingrate(edf), rs=1, order=12)
filter_response(fprototype=:chebyshev2, ftype=:bs, cutoff=[45, 55], fs=eeg_samplingrate(edf), rp=10, order=4)
filter_response(fprototype=:elliptic, ftype=:bs, cutoff=[45, 55], fs=eeg_samplingrate(edf), rs=1, rp=10, order=4)
```

FIR filtering:
```julia
edf = eeg_filter(edf, fprototype=:fir, ftype=:bs, cutoff=[45.0, 55.0], order=8, window=hann(128))
edf = eeg_filter(edf, fprototype=:fir, ftype=:lp, cutoff=45.0, order=8)
edf = eeg_filter(edf, fprototype=:fir, ftype=:hp, cutoff=0.1, order=8)
```

IIR filtering:
```julia
edf = eeg_filter(edf, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8)
edf = eeg_filter(edf, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)
edf = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=8)
edf = eeg_filter(edf, fprototype=:elliptic, ftype=:bs, cutoff=[45.0, 55.0], rs=10, rp=1, order=12)
edf = eeg_filter(edf, fprototype=:chebyshev1, ftype=:bs, cutoff=[45.0, 55.0], rs=1, order=8)
edf = eeg_filter(edf, fprototype=:chebyshev2, ftype=:bs, cutoff=[45.0, 55.0], rp=1, order=8)
```

Plot:
```julia
eeg_plot(edf)
eeg_plot(edf, channel=1:4)
eeg_plot(edf, offset=20*eeg_samplingrate(edf), len=20)
eeg_plot(edf, normalize=false)
eeg_plot(e9, head=true, figure="/tmp/1.png")

eeg_plot_avg(edf, channel=1:4)
eeg_plot_avg(edf)
eeg_plot_avg(e9)
eeg_plot_avg(e10, len=125)
eeg_plot_avg(e9, len=5, offset=6 * eeg_samplingrate(e9))

eeg_plot_butterfly(edf)
eeg_plot_butterfly(edf, offset=20*256, len=120, channel=1:4, normalize=true)
eeg_plot_butterfly(edf, channel=1:4, normalize=true)
eeg_plot(edf, figure="/tmp/test.png")

e9 = eeg_load_electrode_positions(e9, "locs/standard-10-20-cap19.ced")
eeg_plot_butterfly(e9, len=9)
eeg_plot_butterfly(e9, len=55)
eeg_plot_butterfly(e9, len=55, head=true)
eeg_plot_butterfly(e10, len=55, head=true)
eeg_plot_butterfly(e10)
eeg_plot_butterfly(e10, len=11, head=true)

p1 = eeg_plot(edf)
p1 = plot!(title="edf1")
p2 = eeg_plot(edf)
p2 = plot!(title="edf2")
plot(p1, p2, layout=(1, 2))
```

Calculate covariance matrix:
```julia
edf_cov = eeg_cov(edf)
```

Calculate correlation matrix
```julia
edf_cor = eeg_cor(edf)
eeg_plot_matrix(edf, edf_cor)
```

Calculate auto-covariance:
```julia
ac, lags = eeg_autocov(edf, lag=20, normalize=false)
eeg_plot_covmatrix(edf, ac, lags)
```

Calculate cross-covariance:
```julia
cc, lags = eeg_crosscov(edf, lag=20, demean=true)
# channel by channel, all combinations
plot(lags, cc[1, :])

edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=8)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=12)
cc, lags = eeg_crosscov(edf1, edf2, lag=20, demean=true, norm=true)
plot(lags, cc[1, :])
```

Normalize:
```julia
eeg_normalize_zscore(edf)
eeg_normalize_minmax(edf)
```

Remove DC:
```julia
eeg_demean(edf)
eeg_demean(edf)
```

Taper:
```julia
h = hann(edf.eeg_header[:epoch_duration_samples])
eeg_taper(edf, h)
```

Calculate signal derivative:
```julia
eeg_derivative(edf)
```

Detrend:
```julia
eeg_detrend(edf, type=:linear)
eeg_detrend(edf, type=:constant)
```

Calculate signal total power:
```julia
tbp = eeg_total_power(edf)
bar(eeg_labels(edf), tbp, xticks=(1:length(eeg_labels(edf)), eeg_labels(edf)))
```

Calculate signal band power:
```julia
abp = eeg_band_power(edf, f1=8.0, f2=12.0)
bar(eeg_labels(edf), abp, xticks=(1:length(eeg_labels(edf)), eeg_labels(edf)))
```

Get separate channels:
```julia
f3 = eeg_get_channel(edf, "F3")
f4 = eeg_get_channel(edf, 4)
f3_f = signal_filter(f3, fprototype=:butterworth, ftype=:hp, cutoff=0.1, fs=eeg_samplingrate(edf), order=8)
```

Time-domain convolution:
```julia
mw = generate_morlet(256, 1, 32, complex=true)
eeg_tconv(e10, kernel=mw)
```

Frequency-domain convolution:
```julia
mw = generate_morlet(256, 1, 32, complex=true)
eeg_fconv(e10, kernel=mw)
```

Spectral analysis:
```julia
edf_pow, edf_frq = eeg_psd(edf, normalize=true)
plot(edf_frq[10, :], edf_pow[10, :])
eeg_plot_psd(edf, frq_lim=20.0, channel=1:4)
eeg_plot_psd(edf, normalize=true, average=false, frq_lim=50)
eeg_plot_psd(edf, normalize=true, average=true, frq_lim=20)
eeg_plot_psd(edf, channel=1:4, average=true)
f3 = eeg_extract_channel(edf, "F3")
f4 = eeg_extract_channel(edf, 4)
signal_psd(f4, fs=256)
signal_plot_psd(f3, fs=256)
signal_plot_psd(f4, fs=256)

eeg_plot_spectrogram(edf, channel=9, normalize=true)
eeg_plot_spectrogram(e10, channel=9, normalize=true, len=110)
eeg_plot_spectrogram(e9, channel=9, normalize=true, ylim=80, len=75)
eeg_plot_spectrogram(e9, channel=9, normalize=true, ylim=40, len=80, offset=18*256)
eeg_plot_spectrogram(e9, channel=9, normalize=true, ylim=40)
eeg_plot_psd(e10)
eeg_plot_psd(e10, len=120)
eeg_plot_psd(edf, len=240, channel=4, frq_lim=20)
```

Electrode positioning:
```julia
eeg_info(edf)
using NeuroJ
edf = eeg_import_edf("test/eeg-test-edf.edf")
edf = eeg_load_electrode_positions(edf, "locs/standard-10-20-cap19.ced")
eeg_plot_electrodes(edf, labels=true, head=true)
eeg_plot(edf, channel=1:10)
eeg_plot(edf, channel=1:10, head=false)
eeg_plot_butterfly(e10, channel=1:10, head=true, len=55)
eeg_plot_psd(edf, normalize=true, channel=1:19, head=true, figure="/tmp/1.pdf", frq_lim=40)
eeg_plot_electrodes(edf, labels=true, selected=1:, small=false)
eeg_plot_electrodes(edf, labels=true, selected=1:15, small=true)
```

Stationarity:
```julia
p = eeg_stationarity(edf, method=:mean)
p = eeg_stationarity(edf, method=:var)
plot(p[1, :, :], legend=false)

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

Coherence:
```julia
m = eeg_coherence(edf1, edf2)
hz, nyq = eeg_freqs(edf1)
plot(hz, abs.(m[1, 1:length(hz)]), xlims=(0, 40))

edf_alpha = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(:alpha), order=8)
eeg_labels(edf_alpha)
# O1 vs O2
m = eeg_coherence(edf_alpha, channel1=9, channel2=10)
plot(hz, abs.(m[1:length(hz)]), xlims=(0, 20))
# Fp1 vs Fp2
m = eeg_coherence(edf_alpha, channel1=1, channel2=2)
plot!(hz, abs.(m[1:length(hz)]), xlims=(0, 20))
```

PCA:
```julia
edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(:beta), order=8)
edf1 = eeg_epochs(edf1, epoch_len=10*eeg_samplingrate(edf1), average=true)
edf1 = eeg_keep_channel(edf1, 3)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(:beta), order=8)
edf2 = eeg_epochs(edf2, epoch_len=10*eeg_samplingrate(edf1), average=true)
edf2 = eeg_keep_channel(edf2, 4)
pc, pc_var = eeg_pca(edf1, edf2, n=4)
plot(pc[1, :, 1, 1])
bar(vec(pc_var))
```

Comparing two signals:
```julia
edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(:delta), order=8)
edf1 = eeg_epochs(edf1, epoch_len=10*eeg_samplingrate(edf1), average=true)
edf1 = eeg_keep_channel(edf1, 4)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(:beta), order=8)
edf2 = eeg_epochs(edf2, epoch_len=10*eeg_samplingrate(edf1), average=true)
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

Misc:
```julia
eeg_band(:alpha)
alpha_band = eeg_filter(edf, fprototype=:butterworth, ftype=:bp, cutoff=eeg_band(:alpha), order=8)
hz, nyq = eeg_freqs(edf)
```

Benchmarking:
```julia
edf = eeg_import_edf("test/eeg-test-edf.edf")
function eeg_benchmark(n::Int64)
    for idx in 1:n
        edf_new = eeg_reference_car(edf)
        e10 = eeg_epochs(edf_new, epoch_len=10*eeg_samplingrate(edf))
        e10 = eeg_filter(e10, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8)
        e10 = eeg_filter(e10, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)
        e10 = eeg_filter(e10, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=8)
        tbp = eeg_total_power(e10)
        ac = eeg_autocov(e10, norm=false)
        cc = eeg_crosscov(e10, lag=10, demean=true)
        mconv = eeg_tconv(e10, kernel=generate_morlet(256, 1, 32, complex=true))
    end
end
@time eeg_benchmark(10)

# workstation: 71.059454 seconds (116.09 M allocations: 644.196 GiB, 18.44% gc time)
# workstation: 54.460665 seconds (126.95 M allocations: 622.011 GiB, 21.19% gc time, 0.04% compilation time)

# laptop: 92.209796 seconds (114.60 M allocations: 642.148 GiB, 10.88% gc time)
```