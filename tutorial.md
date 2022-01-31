# Neuro.jl

## EEG

```
edf = eeg_load_edf("test.edf")

# save
eeg_save(edf, "test.bin")

# load
edf = eeg_load("test.bin")

# channel index
eeg_get_channel_idx(edf, "Cz")
eeg_get_channel_name(edf, 20)

# channel rename
edf = eeg_rename_channel(edf, "Cz, "CZ")
edf = eeg_rename_channel(edf, 20, "Cz")

# re-reference
eeg = eeg_reference_channel(edf, [1, 2])
eeg = eeg_reference_car(edf)

# filtering
edf = eeg_filter_butter(edf, filter_type=:bs, cutoff=[45, 55], poles=8)
edf = eeg_filter_butter(edf, filter_type=:lp, cutoff=45.0, poles=8)
edf = eeg_filter_butter(edf, filter_type=:hp, cutoff=1.0, poles=8)

# remove channel
edf = eeg_drop_channel(edf, 22)

# show labels
edf.eeg_signal_header[:labels]

# show processing history
edf.eeg_object_header[:history]

# show sampling rate
edf.eeg_signal_header[:sampling_rate][1]

# upsample
eeg_upsample(edf, new_sr=512)

# plot channels
eeg_plot(edf)

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
eeg_total_power(edf)

# alpha power
eeg_band_power(eeg, 8, 12)

# get separate channels
f3 = eeg_get_channel(edf, "F3")
f4 = eeg_get_channel(edf, 4)
```

