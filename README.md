# NeuroJ.jl

Welcome fellow researcher!

## Installation

```
using Pkg
Pkg.add(url="https://notabug.org/AdamWysokinski/NeuroJ.jl")
```

## Documentation

Single-channel signals are column vectors.

Multi-channel or multi-trial signals are matrices of channels/trials by signals.

Function name prefix:
- signal_  :: functions taking signal or signals as argument
- eeg_     :: functions taking EEG object as argument

EEG object (headers + data) is stored in the EEG structure:
```
struct EEG
    eeg_file_header::Dict
    eeg_signal_header::Dict
    eeg_signals::Matrix
end
```

### signal.jl

```
signal_derivative(signal)
```
Returns the derivative of the `signal` vector with length same as the signal.

```
signal_derivative(signal)
```
Returns the derivative of each the `signal` matrix channels with length same as the signal.

```
signal_total_power(signal)
```
Calculates total power for the `signal` vector.

```
signal_total_power(signal)
```
Calculates total power for each the `signal` matrix channels.

```
signal_band_power(signal, fs, f1, f2)
```
Calculates absolute band power between frequencies `f1` and `f2` for the `signal` vector.

```
signal_band_power(signal, fs, f1, f2)
```
Calculates absolute band power between frequencies `f1` and `f2` for each the `signal` matrix channels.

```
make_spectrum(signal, fs)
```
Returns FFT and DFT sample frequencies for a DFT for the `signal` vector.

```
make_spectrum(signal, fs)
```
Returns FFT and DFT sample frequencies for a DFT for each the `signal` matrix channels.

```
signal_detrend(signal, type=:linear)
```
Removes linear trend from the `signal` vector.

```
signal_detrend(signal, type=:linear)
```
Removes linear trend for each the `signal` matrix channels.

```
signal_ci95(signal, n=3; method=:normal)
```
Calculates mean, std and 95% confidence interval for each the `signal` matrix channels.

```
signal_mean(signal1, signal2)
```
Calculates mean and 95% confidence interval for 2 signals.

```
signal_difference(signal1::Matrix, signal2::Matrix, n=3; method=:absdiff)
```
Calculates mean difference and 95% confidence interval for 2 signals.

```
signal_autocov(signal, lag=1, remove_dc=false, normalize=false)
```
Calculates autocovariance of the `signal` vector.

```
signal_autocov(signal, lag=1, remove_dc=false, normalize=false)
```
Calculates autocovariance of each the `signal` matrix channels.

```
signal_crosscov(signal1, signal2, lag=1, remove_dc=false, normalize=false)
```
Calculates cross-covariance between `signal1` and `signal2` vectors.

```
signal_crosscov(signal1, signal2, lag=1, remove_dc=false, normalize=false)
```
Calculates cross-covariance between same channels in `signal1` and `signal2` matrices.

```
signal_crosscov(signal, lag=1, remove_dc=false, normalize=false)
```
Calculates cross-covariance for all channels in the `signal` matrix

```
signal_spectrum(signal, pad=0, remove_dc=false, detrend=false, taper=nothing)
```
Calculates FFT, amplitudes, powers and phases of the `signal` vector.

```
signal_spectrum(signal, pad=0, remove_dc=false, detrend=false, taper=nothing)
```
Calculates FFT, amplitudes, powers and phases for each channel of the `signal` matrix.

```
signal_epoch(signal; epoch_no, epoch_len, average=true, remove_dc=false, detrend=false, derivative=false, taper=nothing)
```
Splits `signal` vector into epochs.

```
signal_epoch(signal, n; average=true, remove_dc=false, detrend=false, derivative=false, taper=nothing)
```
Splits `signal` matrix into epochs.

```
signal_filter_butter(signal; filter_type, cutoff, fs, poles=8)
```
Filters `signal` vector using Butterworth filter.

```
signal_filter_butter(signal; filter_type, cutoff, fs, poles=8)
```
Filters `signal` matrix using Butterworth filter.

```
signal_plot(t, signal; offset=0, labels=[], normalize=false, xlabel="Time [s]", ylabel="Amplitude [μV]", figure="")
```
Plots `signal` against time vector `t`.

```
signal_plot(t, signal; channels=[], labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels", figure="")
```
Plots `signal` matrix against time vector `t`.

```
signal_drop_channel(signal, channels)
```
Removes `channels` from the `signal` matrix.

```
signal_rereference_channel(signal, reference)
```
Re-references channels of the `signal` matrix to specific signal channel.

### eeg.jl

Functions for importing, processing, analyzing and displaying EEG data go here.

```
eeg_load_edf(in_file, read_annotations=true, header_only=false, clean_labels=true)
```
Loads EDF/EDFPlus file.

```
eeg_plot(eeg; t=nothing, offset=0, channels=[], labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels")
```
Plots `eeg` signals.

```
eeg_drop_channel(eeg, channels)
```
Removes `channels` from the `eeg` set.

```
eeg_filter_butter(eeg; filter_type, cutoff, fs, poles=8)
```
Filters `eeg` channels using Butterworth filter.

```
eeg_derivative(eeg)
```
Returns the derivative of each the `eeg` channels with length same as the signal.

```
eeg_total_power(eeg)
```
Calculates total power for each the `eeg` signal channels.

```
eeg_band_power(eeg, f1, f2)
```
Calculates absolute band power between frequencies `f1` and `f2` for each the `eeg` signal channels.

```
eeg_make_spectrum(eeg)
```
Returns FFT and DFT sample frequencies for a DFT for each the `eeg` signal channels.

```
eeg_detrend(eeg, type=:linear)
```
Removes linear trend for each the `eeg` signal channels.

```
eeg_draw_head(p, loc_x, loc_y)
```
Draws head over a topographical plot `p`.

```
eeg_rereference_channel(eeg, reference)
```
Re-references the `eeg` signal channels to specific signal channel.

```
eeg_save(eeg, file_name)
```
Saves the `eeg` object to `file_name` file.

```
eeg_load(file_name)
```
Loads the `eeg` object from `file_name` file.

### mri.jl

### nirs.jl

### nstim.jl

- neuromodulation modeling

### misc.jl

Various, miscellaneous, general or not belonging to above categories, functions go here.

```
linspace(start, stop, length)
```
Generates `length`-long sequence of evenly spaced numbers between `start` and `stop`.

```
logspace(start, stop, length)
```
Generates `length`-long sequence of log10-spaced numbers between `start` and `stop`.

```
zero_pad(m)
```
Pads the matrix `m` with zeros to make it square.

```
vsearch(x, y; return_distance=false)
```
Returns the positions of the `y` value in the vector `x`.

```
vsearch(x, y)
```
Returns the positions of the `y` vector in the vector `x`.

```
cart2pol(x, y)
```
Converts cartographic coordinates `x` and `y` to polar.

```
pol2cart(theta, rho)
```
Converts polar coordinates `theta` and `rho` to cartographic.

```
cvangle(x)
```
Returns the phase angles, in radians, of the vector `x` with complex elements.

```
hann(n)
```
Returns the `n`-point long symmetric Hanning window.

```
hildebrand_rule(x)
```
Calculates Hildebrand rule for vector `x`.
If H < 0.2 then the vector `x` is symmetrical.

```
jaccard_similarity(x, y)
```
Calculates Jaccard similarity between two vectors `x` and `y`.

```
fft0(x, n)
```
Calculates FFT for the vector `x` padded with `n` zeros at the end.

```
ifft0(x, n)
```
Calculates IFFT for the vector `x` padded with `n` zeros at the end.

```
nexpow2(x)
```
Returns the next power of 2 for given number `x`.

```
vsplit(x, n)
```
Splits the vector `x` into `n`-long pieces.

```
rms(x)
```
Calculates Root Mean Square of the vector `x`.

```
db(x)
```
Converts values of the vector `x` to dB.

```
sine(f, t, a, p)
```
Generates sine wave of `f` frequency over `t` time; optional arguments are: `a` amplitude and  `p` phase.

```
frequencies(t)
```
Returns vector of frequencies and Nyquist frequency for given time vector `t`.

```
matrix_sortperm(m; dims=1)
```
Generates sorting index for matrix `m` by columns (`dims` = 1) or by rows (`dims` = 2).

```
matrix_sort(m, m_idx; dims=1)
```
Sorts matrix `m` using sorting index `m_idx` by columns (`dims` = 1) or by rows (`dims` = 2).

```
pad0(x, n)
```
Pads the vector `x` with `n` zeros at the beginning and at the end.

```
hz2rads(f)
```
Converts frequency `f` in Hz to rad/s.

```
rads2hz(f)
```
Converts frequency `f` in rad/s to Hz.

```
z_score(x)
```
Calculates Z-scores for each value of the vector `x`.

```
k(n)
```
Calculates number of categories for a given sample size `n`.

```
demean(signal)
```
Demean `signal` vector.

```
normalize_mean(signal)
```
Normalize (scales around the mean) `signal` vector.

```
normalize_minmax(signal)
```
Normalize (to 0…1) `signal` vector.

```
mag2db(x)
```
Converts magnitude of the vector or matrix `x` into dB values. Maximum value of `x` is 0 dB.

```
db2mag(x)
```
Converts dB values of the vector or matrix `x` into magnitude.

```
pow2db(x)
```
Converts power values of the vector or matrix `x` into dB values. Maximum value of `x` is 0 dB.

```
db2pow(x)
```
Converts dB values of the vector or matrix `x` into corresponding power value.

## TO DO

EEG:
- re-reference
- generate time vector and add to eeg_signal_header
- import channel location files
- channel locations data to eeg_signal_header

## Contributing

Please feel free to contribute documentation, tutorials, paid and free books.

## Contributors

If you've contributed, add your name below!

[Adam Wysokiński](adam.wysokinski@umed.lodz.pl)

## License

The program is licensed under GPL-2.0-only.