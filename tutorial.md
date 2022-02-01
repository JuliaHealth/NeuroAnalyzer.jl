# Neuro.jl

## EEG

```
using NeuroJ

edf = eeg_load_edf("eeg-test.edf")

# save
eeg_save(edf, "test.bin")

# load
edf = eeg_load("test.bin")

# show info
eeg_info(edf)

# split into epochs
eeg_epochs(edf, epochs_no=10)
# 2-second epochs
eeg_epochs(edf, epochs_len=2 * eeg_signal_header[:sampling_rate][1], average=true)

# show labels
edf.eeg_signal_header[:labels]

# show properties
edf.eeg_signal_header[:sampling_rate][1]
edf.eeg_object_header[:eeg_duration_seconds]
edf.eeg_object_header[:eeg_duration_samples]

# channel index
eeg_get_channel_idx(edf, "Cz")
eeg_get_channel_name(edf, 18)

# channel rename
edf = eeg_rename_channel(edf, "Cz", "CZ")
edf = eeg_rename_channel(edf, 18, "Cz")

# re-reference
eeg = eeg_reference_channel(edf, [1, 2])
eeg = eeg_reference_car(edf)

# filtering
edf = eeg_filter_butter(edf, filter_type=:bs, cutoff=[45.0, 55.0], poles=8)
edf = eeg_filter_butter(edf, filter_type=:lp, cutoff=45.0, poles=8)
edf = eeg_filter_butter(edf, filter_type=:hp, cutoff=1.0, poles=8)

# remove channel
edf = eeg_drop_channel(edf, 10)

# show processing history
eeg_show_processing_history(edf)

# show sampling rate
edf.eeg_signal_header[:sampling_rate][1]

# upsample
eeg_upsample(edf, new_sr=512)

# plot channels
eeg_plot(edf)
eeg_plot(edf, figure="figure1.pdf")

# covariance
edf_cov = eeg_cov(edf)
heatmap(edf_cov)

# correlation
edf_cor = eeg_cor(edf)
heatmap(edf_cor)

# normalize
eeg_normalize_mean(edf)
eeg_normalize_minmax(edf)

# remove DC
eeg_demean(edf)

# derivative
eeg_derivative(edf)

# detrend
eeg_detrend(edf)

# total power
tbp = eeg_total_power(edf)
bar(edf.eeg_signal_header[:labels],
    tbp,
    xticks=(1:length(edf.eeg_signal_header[:labels]),
    edf.eeg_signal_header[:labels]))

# alpha power
abp = eeg_band_power(edf, f1=8.0, f2=12.0)
bar(edf.eeg_signal_header[:labels],
    abp,
    xticks=(1:length(edf.eeg_signal_header[:labels]),
    edf.eeg_signal_header[:labels]))

# get separate channels
f3 = eeg_get_channel(edf, "F3")
f4 = eeg_get_channel(edf, 4)
```

