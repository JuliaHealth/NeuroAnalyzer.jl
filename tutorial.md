# Neuro.jl

## EEG

```
edf = eeg_load_edf("test.edf")

eeg = eeg_rereference_channel(edf, 1)

edf = eeg_filter_butter(edf, filter_type=:bs, cutoff=[45, 55], poles=8)
edf = eeg_filter_butter(edf, filter_type=:lp, cutoff=45.0, poles=8)
edf = eeg_filter_butter(edf, filter_type=:hp, cutoff=1.0, poles=8)

eeg_plot(edf)
```