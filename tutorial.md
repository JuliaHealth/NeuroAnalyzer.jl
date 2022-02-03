# Neuro.jl

## EEG

```julia
using Pkg
Pkg.update()
using NeuroJ

# load EDF file
edf = eeg_load_edf("eeg-test.edf")

# show properties
eeg_info(edf)
eeg_info(e2avg)
edf.eeg_header[:sampling_rate][1]
edf.eeg_header[:eeg_duration_seconds]
edf.eeg_header[:eeg_duration_samples]
eeg_samplingrate(edf)

# show labels
eeg_labels(edf)

# save
eeg_save(edf, "test.bin")

# load
edf = eeg_load("test.bin")

# get channel index
eeg_get_channel_idx(edf, "Cz")
eeg_get_channel_name(edf, 18)

# channel rename
edf = eeg_rename_channel(edf, "Cz", "CZ")
edf = eeg_rename_channel(edf, 18, "Cz")

# show processing history
eeg_history(edf)

# upsample
eeg_samplingrate(edf)
eeg_info(edf)
edf_512 = eeg_upsample(edf, new_sr=512)
eeg_info(edf_512)

# downsample
eeg_samplingrate(edf)
edf_128 = eeg_downsample(edf, new_sr=128)
eeg_info(edf_128)

# split into 10 epochs
e10 = eeg_epochs(edf, epochs_no=10)
eeg_info(e10)

# split into 10-second epochs
e10 = eeg_epochs(edf, epochs_len=10*eeg_samplingrate(edf))
eeg_info(e10)

# get 1st epoch
e10e1 = eeg_get_epoch(e10, 1)
eeg_info(e10e1)

# split into 5-second epochs, average
e2avg = eeg_epochs(edf, epochs_len=5*eeg_samplingrate(edf), average=true)
eeg_info(e2avg)

# re-reference
edf = eeg_reference_channel(edf, [1, 2])
edf = eeg_reference_channel(edf, 18)
edf = eeg_reference_car(edf)

# filtering
filter_response(fprototype=:butterworth, ftype=:hp, cutoff=10, fs=eeg_samplingrate(edf), order=8)
filter_response(fprototype=:chebyshev1, ftype=:bp, cutoff=[40, 50], fs=eeg_samplingrate(edf), rs=1, order=12)
filter_response(fprototype=:chebyshev2, ftype=:bs, cutoff=[45, 55], fs=eeg_samplingrate(edf), rp=10, order=4)
## FIR
edf_fir = eeg_filter(edf, fprototype=:fir, ftype=:bs, cutoff=[45.0, 55.0], order=8, window=hann(128))
edf_fir = eeg_filter(edf_fir, fprototype=:fir, ftype=:lp, cutoff=45.0, order=8)
edf_fir = eeg_filter(edf_fir, fprototype=:fir, ftype=:hp, cutoff=0.1, order=8)
## IIR
edf = eeg_filter(edf, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8)
edf = eeg_filter(edf, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)
edf = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=8)
edf_cheb1 = eeg_filter(edf, fprototype=:chebyshev1, ftype=:bs, cutoff=[45.0, 55.0], rs=1, order=8)
edf_cheb2 = eeg_filter(edf, fprototype=:chebyshev2, ftype=:bs, cutoff=[45.0, 55.0], rp=1, order=8)

# plot
eeg_plot(edf)
eeg_plot(e10e1)

# remove channel
edf19 = eeg_delete_channel(edf, 10:18)

# keep channel
edf14 = eeg_keep_channel(edf, 1:4)

# plot channels
eeg_plot(edf)
eeg_plot(edf, channels=1:4)
eeg_plot(edf, offset=20, len=20)
eeg_plot(edf, normalize=false)
eeg_plot_avg(edf, channels=1:4, offset=20)
eeg_plot_avg(edf)
eeg_plot_butterfly(edf)
eeg_plot_butterfly(edf, offset=20*256, len=120, channels=1:4, normalize=true)
eeg_plot_butterfly(e2avg)
eeg_plot_avg(e2avg)
eeg_plot(edf, figure="/test.png")
eeg_plot(edf, figure="/tmp/test.png")
eeg_plot(edf, offset=60, figure="/tmp/test.png")
p1 = eeg_plot(e2avg)
p1 = plot!(title="e2avg")
p2 = eeg_plot(edf_but)
p2 = plot!(title="butterworth")
p3 = eeg_plot(edf_fir)
p3 = plot!(title="FIR")
plot(p1, p2, p3, layout=(1, 3))
eeg_plot(edf, offset=60)
eeg_plot(e10, epoch=10)
eeg_plot(e2avg)
eeg_plot(edf, figure="figure1.pdf")

# covariance
edf_cov = eeg_cov(edf)
heatmap(edf_cov)
edf_cov = eeg_cov(e10)
heatmap(edf_cov[:, :, 1])

# correlation
edf_cor = eeg_cor(edf)
heatmap(edf_cor)

# autocovariance
ac = eeg_autocov(e10, normalize=false)
heatmap(ac)
cc = eeg_crosscov(e10, lag=100, demean=true)
heatmap(cc)
edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=8)
edf2 = eeg_filter(edf, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=12)
cc = eeg_crosscov(edf1, edf2, lag=100)
heatmap(cc)

# normalize
eeg_normalize_mean(edf)
eeg_normalize_minmax(edf)

# remove DC
eeg_demean(edf)

# taper with
h = hann(e10.eeg_header[:epoch_duration_samples])
e10_t = eeg_taper(e10, h)
eeg_plot(e10_t)

# derivative
eeg_derivative(e10)

# detrend
edf1 = eeg_detrend(edf, type=:linear)
e2avg = eeg_detrend(e2avg, type=:linear)

# total power
tbp = eeg_total_power(e10)
bar(eeg_labels(e10),
    tbp,
    xticks=(1:length(eeg_labels(e10)), eeg_labels(e10)))

# alpha power
abp = eeg_band_power(e10, f1=8.0, f2=12.0)
bar(eeg_labels(e10),
    abp,
    xticks=(1:length(eeg_labels(e10)), eeg_labels(e10)))

# get separate channels
f3 = eeg_get_channel(edf, "F3")
f4 = eeg_get_channel(edf, 4)
f3_f = signal_filter(f3, fprototype=:butterworth, ftype=:hp, cutoff=0.1, fs=eeg_samplingrate(edf), order=8)

# time-domain convolution
mw = morlet(256, 1, 32, complex=true)
eeg_tconv(e10, kernel=mw)

# PSD
using Pkg
Pkg.update()
using NeuroJ
edf = eeg_load_edf("eeg-test.edf")
edf_pow, edf_frq = eeg_psd(edf, normalize=true)
plot(edf_frq[10, :], edf_pow[10, :])
eeg_plot_psd(edf, frq_lim=20.0)
eeg_plot_psd(edf, normalize=true, average=false, frq_lim=50)
eeg_plot_psd(edf, normalize=true, average=true, frq_lim=20)
eeg_plot_psd(edf, channels=1:4, average=true)
f3 = eeg_get_channel(edf, "F3")
t=collect(0:1/fs:length(f3))
signal_plot(t, f3)
f4 = eeg_get_channel(edf, 4)
signal_psd(f4, fs=256)
signal_plot_psd(f3, fs=256)
signal_plot_psd(f4, fs=256)

# benchmarking
using BenchmarkTools
function eeg_benchmark(n::Int64)
    for idx in 1:n
        edf_new = eeg_reference_car(edf)
        edf10 = eeg_epochs(edf_new, epochs_len=10*eeg_samplingrate(edf))
        edf10 = eeg_filter(edf10, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8)
        edf10 = eeg_filter(edf10, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)
        edf10 = eeg_filter(edf10, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=8)
        tbp = eeg_total_power(edf10)
        ac = eeg_autocov(edf10, normalize=false)
        cc = eeg_crosscov(edf10, lag=10, demean=true)
        mconv = eeg_tconv(e10, kernel=morlet(256, 1, 32, complex=true))
    end
end
@time eeg_benchmark(10)
```