# Neuro.jl

## EEG

```julia
using NeuroJ

# load EDF file
edf = eeg_load_edf("eeg-test.edf")

# show properties
eeg_info(edf)
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
edf = eeg_reference_car(edf)
eeg_reference_car(e10)

# filtering
edf = eeg_filter_butter(edf, filter_type=:bs, cutoff=[45.0, 55.0], poles=8);
edf = eeg_filter_butter(edf, filter_type=:lp, cutoff=45.0, poles=8);
edf = eeg_filter_butter(edf, filter_type=:hp, cutoff=1.0, poles=8);

# remove channel
edf1 = eeg_drop_channel(edf, 10)

# show processing history
eeg_show_processing_history(edf)
eeg_show_processing_history(e10)
eeg_show_processing_history(e2avg)

# upsample
edf_512 = eeg_upsample(edf, new_sr=512)

# plot channels
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
```
