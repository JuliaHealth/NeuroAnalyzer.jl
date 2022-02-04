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
NeuroJ_version()
```

Load EDF file:
```julia
edf = eeg_import_edf("eeg-test.edf")
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

Split into 10-second epochs:
```julia
e10 = eeg_epochs(edf, epochs_len=10*eeg_samplingrate(edf))
eeg_info(e10)
```

Get 1st epoch:
```julia
e10e1 = eeg_get_epoch(e10, 1)
eeg_info(e10e1)
```

Split into 5-second averaged epoch
```julia
e2avg = eeg_epochs(edf, epochs_len=5*eeg_samplingrate(edf), average=true)
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
edf = eeg_filter(edf, fprototype=:elliptic, ftype=:bs, cutoff=[45.0, 55.0], rs=1, rp=10, order=12)
edf = eeg_filter(edf, fprototype=:chebyshev1, ftype=:bs, cutoff=[45.0, 55.0], rs=1, order=8)
edf = eeg_filter(edf, fprototype=:chebyshev2, ftype=:bs, cutoff=[45.0, 55.0], rp=1, order=8)
```

Plot:
```julia
eeg_plot(edf)
eeg_plot(edf, channels=1:4)
eeg_plot(edf, offset=20*eeg_samplingrate(edf), len=20)
eeg_plot(edf, normalize=false)
eeg_plot_avg(edf, channels=1:4, offset=20)
eeg_plot_butterfly(edf)
eeg_plot_butterfly(edf, offset=20*256, len=120, channels=1:4, normalize=true)
eeg_plot(edf, figure="/tmp/test.png")

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

edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=8)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=12)
cc, lags = eeg_crosscov(edf1, edf2)
```

Normalize:
```julia
eeg_normalize_mean(edf)
eeg_normalize_minmax(edf)
```

Remove DC:
```julia
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
mw = morlet(256, 1, 32, complex=true)
eeg_tconv(e10, kernel=mw)
```

Spectral analysis:
```julia
edf_pow, edf_frq = eeg_psd(edf, normalize=true)
plot(edf_frq[10, :], edf_pow[10, :])
eeg_plot_psd(edf, frq_lim=20.0, channels=1:4)
eeg_plot_psd(edf, normalize=true, average=false, frq_lim=50)
eeg_plot_psd(edf, normalize=true, average=true, frq_lim=20)
eeg_plot_psd(edf, channels=1:4, average=true)
f3 = eeg_get_channel(edf, "F3")
f4 = eeg_get_channel(edf, 4)
signal_psd(f4, fs=256)
signal_plot_psd(f3, fs=256)
signal_plot_psd(f4, fs=256)
```

Electrode positioning:
```julia
eeg_info(edf)
edf = eeg_load_electrode_positions(edf, in_file)
eeg_plot_electrodes(edf, labels=false)
p = eeg_plot(edf, channels=1:10)
h = eeg_plot_electrodes(edf, labels=true, selected=1:10)
plot(p, h, layout=(1, 2))
```

Benchmarking:
```julia
using BenchmarkTools
function eeg_benchmark(n::Int64)
    for idx in 1:n
        edf_new = eeg_reference_car(edf)
        e10 = eeg_epochs(edf, epochs_len=10*eeg_samplingrate(edf))
        e10 = eeg_filter(e10, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8)
        e10 = eeg_filter(e10, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)
        e10 = eeg_filter(e10, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=8)
        tbp = eeg_total_power(e10)
        ac = eeg_autocov(e10, normalize=false)
        cc = eeg_crosscov(e10, lag=10, demean=true)
        mconv = eeg_tconv(e10, kernel=morlet(256, 1, 32, complex=true))
    end
end
@time eeg_benchmark(10)
# 71.059454 seconds (116.09 M allocations: 644.196 GiB, 18.44% gc time)
```