# Neuro.jl

## EEG

```
edf = eeg_load_edf("test.edf")

eeg = eeg_rereference_channel(edf, 1)
eeg = eeg_rereference_car(edf)

edf = eeg_filter_butter(edf, filter_type=:bs, cutoff=[45, 55], poles=8)
edf = eeg_filter_butter(edf, filter_type=:lp, cutoff=45.0, poles=8)
edf = eeg_filter_butter(edf, filter_type=:hp, cutoff=1.0, poles=8)

# plot channels
eeg_plot(edf)

# remove channel
edf = eeg_drop_channel(edf, 22)

# show labels
edf.eeg_signal_header[:labels]

# show processing history
edf.eeg_object_header[:history]

# show sampling rate
edf.eeg_signal_header[:sampling_rate][1]

# save
eeg_save(edf, "test.bin")

# load
edf = eeg_load("test.bin")
```