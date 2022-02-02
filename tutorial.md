# Neuro.jl

## EEG

```julia
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
eeg_save(edf, "test.bin", overwrite=true)

# load
edf = eeg_load("test.bin")

# split into 10 epochs
e10 = eeg_epochs(edf, epochs_no=10)
eeg_info(e10)

# split into 5-second epochs, average
e2avg = eeg_epochs(edf, epochs_len=10*256, average=true)
eeg_info(e2avg)

# get 1st epoch
e10e1 = eeg_get_epoch(e10, 1)

# get channel index
eeg_get_channel_idx(edf, "Cz")
eeg_get_channel_name(edf, 18)

# channel rename
edf = eeg_rename_channel(edf, "Cz", "CZ")
edf = eeg_rename_channel(edf, 18, "Cz")

# re-reference
edf = eeg_reference_channel(edf, [1, 2])
edf = eeg_reference_channel(edf, 18)
edf = eeg_reference_car(edf)
eeg_reference_car(e10)

# filtering
filter_response(fprototype=:butterworth, ftype=:hp, cutoff=0.1, fs=eeg_samplingrate(edf), order=8, response=true)
## FIR
edf_fir = eeg_filter(e2avg, fprototype=:fir, ftype=:bs, cutoff=[45.0, 55.0], order=8)
edf_fir = eeg_filter(edf_fir, fprototype=:fir, ftype=:lp, cutoff=45.0, order=8)
edf_fir = eeg_filter(edf_fir, fprototype=:fir, ftype=:hp, cutoff=0.1, order=8)

## IIR
edf_but = eeg_filter(e2avg, fprototype=:butterworth, ftype=:bs, cutoff=[45.0, 55.0], order=8)
edf_but = eeg_filter(edf_but, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8)
edf_but = eeg_filter(edf_but, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)

# remove channel
edf1 = eeg_drop_channel(edf, 10)

# show processing history
eeg_history(edf)
eeg_history(e10)
eeg_history(e2avg)

# upsample
edf_512 = eeg_upsample(edf, new_sr=512)

# plot channels
p1 = eeg_plot(e2avg)
p2 = eeg_plot(edf_but)
p3 = eeg_plot(edf_fir)
plot(p1, p2, p3, layout=(1, 3))

eeg_plot(edf)
eeg_plot(edf_512, offset=60*256)
eeg_plot(edf1, offset=100*256)
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

# normalize
eeg_normalize_mean(edf)
eeg_normalize_minmax(e10)

# remove DC
eeg_demean(e10)

# taper with
h = hann(edf.eeg_header[:epoch_duration_samples])
edf_t = eeg_taper(edf, h);
eeg_plot(e10_t)

# derivative
eeg_derivative(e10)

# detrend
edf1 = eeg_detrend(edf, type=:linear)
e2avg = eeg_detrend(e2avg, type=:linear)

# total power
tbp = eeg_total_power(e2avg)
bar(eeg_labels(e2avg),
    tbp,
    xticks=(1:length(eeg_labels(e2avg)), eeg_labels(e2avg)))

# alpha power
abp = eeg_band_power(e2avg, f1=8.0, f2=12.0)
bar(eeg_labels(e2avg),
    abp,
    xticks=(1:length(eeg_labels(e2avg)), eeg_labels(e2avg)))

# get separate channels
f3 = eeg_get_channel(edf, "F3")
f4 = eeg_get_channel(edf, 4)
f3_f = signal_filter(f3, fprototype=:butterworth, ftype=:hp, cutoff=0.1, fs=eeg_samplingrate(edf), order=8)

# time-domain convolution
mw = morlet(256, 1, 10, complex=false)
eeg_tconv(edf, mw)

```
