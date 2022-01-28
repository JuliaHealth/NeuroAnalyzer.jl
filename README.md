# NeuroJ.jl

Welcome fellow researcher!

## Installation

```
using Pkg
Pkg.add(url="https://notabug.org/AdamWysokinski/NeuroJ.jl")
```

## Conventions

Single-channel signals are column vectors.

Multi-channel or multi-trial signals are matrices of channels/trials by signals.

Function name prefix:
- signal_  :: functions taking signal or signals as argument
- eeg_     :: functions taking EEG object as argument

## Sub-modules

### signal.jl

### eeg.jl

### mri.jl

### nirs.jl

### nstim.jl

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
    vsearch(x::Vector, y::Number; return_distance=false)
```
Returns the positions of the `y` value in the vector `x`.

```
    vsearch(x::Vector, y::Vector)
```
Returns the positions of the `y` vector in the vector `x`.

```
    cart2pol(x, y)

Converts cartographic coordinates `x` and `y` to polar.
```

```
    pol2cart(theta, rho)

Converts polar coordinates `theta` and `rho` to cartographic.
```

```
    cvangle(x)

Returns the phase angles, in radians, of the vector `x` with complex elements.
```

```
    hann(n)

Returns the `n`-point long symmetric Hanning window.
```

```
    hildebrand_rule(x)

Calculates Hildebrand rule for vector `x`.
If H < 0.2 then the vector `x` is symmetrical.
```

```
    jaccard_similarity(x, y)

Calculates Jaccard similarity between two vectors `x` and `y`.
```

```
    fft0(x, n)

Calculates FFT for the vector `x` padded with `n` zeros at the end.
```

```
    ifft0(x, n)

Calculates IFFT for the vector `x` padded with `n` zeros at the end.
```

```
    nexpow2(x)

Returns the next power of 2 for given number `x`.
```

```
    vsplit(x, n)

Splits the vector `x` into `n`-long pieces.
```

```
    rms(x)

Calculates Root Mean Square of the vector `x`.
```

```
    db(x)

Converts values of the vector `x` to dB.
```

```
    sine(f, t, a, p)

Generates sine wave of `f` frequency over `t` time; optional arguments are: `a` amplitude and  `p` phase.
```

```
    frequencies(t)

Returns vector of frequencies and Nyquist frequency for given time vector `t`.
```

```
    matrix_sortperm(m::AbstractMatrix; dims=1)

Generates sorting index for matrix `m` by columns (`dims` = 1) or by rows (`dims` = 2).
```

```
    matrix_sort(m::Matrix, m_idx::Vector{Int}; dims=1)

Sorts matrix `m` using sorting index `m_idx` by columns (`dims` = 1) or by rows (`dims` = 2).
```

```
    pad0(x, n)

Pads the vector `x` with `n` zeros at the beginning and at the end.
```

```
    hz2rads(f)

Converts frequency `f` in Hz to rad/s.
```

```
    rads2hz(f)

Converts frequency `f` in rad/s to Hz.
```

```
    z_score(x)

Calculates Z-scores for each value of the vector `x`.
```

```
    k(n)

Calculates number of categories for a given sample size `n`.
```

```
    demean(signal)

Demean `signal` vector.
```

```
    normalize_mean(signal)

Normalize (scales around the mean) `signal` vector.
```

```
    normalize_minmax(signal)

Normalize (to 0…1) `signal` vector.
```

## Documentation

```
struct EEG
    eeg_file_header::Dict
    eeg_signal_header::Dict
    eeg_signals::Matrix
end
```

## TO DO

- save figures (PDF, PNG) for eeg_plot(), signal_plot()

## Contributing

Please feel free to contribute documentation, tutorials, paid and free books.

## Contributors

If you've contributed, add your name below!

[Adam Wysokiński](adam.wysokinski@umed.lodz.pl)

## License

The program is licensed under GPL-2.0-only.