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

# show labels
eeg_show_labels(edf)

# save
eeg_save(edf, "test.bin")
eeg_save(edf, "test.bin", overwrite=true)

# load
edf = eeg_load("test.bin")

# split into 10 epochs
e10 = eeg_epochs(edf, epochs_no=10)

# split into 2-second epochs, average
e2avg = eeg_epochs(edf, epochs_len=512, average=true)

# get channel index
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

# upsample
eeg_upsample(edf, new_sr=512)

# plot channels
eeg_plot(edf)
eeg_plot(e2avg)
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
tbp = eeg_total_power(e2avg)
bar(eeg_labels(e2avg, tbp, xticks=(1:length(eeg_labels(e2avg)), eeg_labels(e2avg)))

# alpha power
abp = eeg_band_power(edf, f1=8.0, f2=12.0)
bar(edf.eeg_header[:labels],
    abp,
    xticks=(1:length(edf.eeg_header[:labels]),
    edf.eeg_header[:labels]))

# get separate channels
f3 = eeg_get_channel(edf, "F3")
f4 = eeg_get_channel(edf, 4)
```

