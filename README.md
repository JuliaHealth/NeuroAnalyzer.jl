![neuroj](images/neuroj.png)

Welcome fellow researcher!

NeuroJ.jl is a [Julia](https://julialang.org) package for analyzing of EEG data. Future versions will also process NIRS data and use MRI data for source localization techniques. Also, various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will be included.

NeuroJ.jl contains a set of separate (high-level) functions, it does not have a graphical user interface (although one could built it upon these). NeuroJ.jl functions can be combined into an analysis pipeline, i.e. a Julia script containing all steps of your analysis. This combined with processing power of Julia language and easiness of distributing calculations across computing cluster, will make NeuroJ.jl particularly useful for processing large amounts of research data.

NeuroJ.jl is a non-commercial project, developed for researchers in psychiatry, neurology and neuroscience.

Initially NeuroJ.jl will be focused on resting-state EEG analysis, but ERP and other type of analyses will be developed in future versions. The goal is to make a powerful, expandable and elastic environment for EEG/MRI/NIRS/NIBS analyses.

Every contribution (bug reports, fixes, new ideas, feature requests or additions, documentation improvements, etc.) to the project is highly welcomed.

## Installation

```
using Pkg
Pkg.add(url="https://codeberg.org/AdamWysokinski/NeuroJ.jl")
```

## Requirements

Julia version ≥ 1.0 is required. Julia [current stable version](https://julialang.org/downloads/#current_stable_release) is recommended, as NeuroJ.jl is only tested against it.

The following packages are required:
- CSV
- DataFrames
- Distances
- DSP
- FFTW
- HypothesisTests
- InformationMeasures
- Interpolations
- JLD2
- LinearAlgebra
- Loess
- MultivariateStats
- Pkg
- Plots
- Polynomials
- ScatteredInterpolation
- [Simpson](https://codeberg.org/AdamWysokinski/Simpson.jl)
- StatsKit
- StatsPlots

NeuroJ.jl will be 100% Julia based.

## General remarks

NeuroJ.jl functions operate on NeuroJ objects (NeuroJ.EEG: EEG metadata header + time + epoched signals + components).

EEG signal is `Array{Float64, 3}` (channels × signals × epochs). If epochs are not defined, the whole signal is an epoch, i.e. there is always at least one epoch.

Functions name prefix:
- `eeg_` functions taking EEG object as an argument
- `eeg_plot_` plotting functions

The majority of `eeg_` functions will process all channels and epochs of the input EEG object. To process individual channels/epochs, you need to extract them from the EEG object first (`eeg_keep_epoch()`, `eeg_keep_channel()` to process as NeuroJ.EEG object or `eeg_extract_channel()`, `eeg_extract_epoch()` to process as multi-channel array)

`eeg_` functions use named arguments for all arguments other than input signal(s), e.g. `eeg_delete_epoch!(my_eeg, epoch=12)`.

EEG object (headers + time + components + EEG signal) is stored in the EEG structure:
```julia
mutable struct EEG
    eeg_header::Dict
    eeg_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
    eeg_components::Vector{Any}
end
```

Many `eeg_` functions have a mutator variant (e.g. `eeg_delete_epoch!()`). These functions modifies the input EEG object, e.g. `eeg_delete_channel!(my_eeg, channel=1)` instead of `my_eeg = eeg_delete_channel(my_eeg, channel=1)`.

## Documentation

NeuroJ.jl documentation is available [here](https://codeberg.org/AdamWysokinski/NeuroJ.jl/src/master/Documentation.md).

Tutorial introducing NeuroJ.jl functions is [here](https://codeberg.org/AdamWysokinski/NeuroJ.jl/src/master/Tutorial.md).

## Plugins (extensions)

Rudimentary plugin architecture is already available, just put .jl plugin scripts into `~/Documents/NeuroJ/plugins`.

Run `neuroj_reload_plugins()` to refresh plugins.

```julia
neuroj_reload_plugins()
neuroj_plugin_demo()
```

## Known bugs

- ignore non-eeg channels for processing, analysis and plotting; currently NeuroJ does not analyze/process/plot EEG containing non-eeg channels, you have to manually extract these to another EEG object 

## To do

.. so much to do ..

The lists below are not in any particular order.

General:
- performance optimization
- CUDA/AMD ROCm acceleration

EEG:
- new plotting architecture:
    - eeg_plot() meta function
    - update tutorial.md
    - check for wrong epoch number
    - raw component vs reconstructed signal (ica)
    - tile plot
- complex kernel convolution: plot magnitude and phase of the convoluted signal
- multi-trial data
- split signal into frequency bands
- continuous wavelet transform (using ContinuousWavelets.jl)
- spectrogram: mean temporal and spectral envelope
- conectomes graph
- coherence spectrum (y: relative amplitude, x: frequencies)
- spectrogram: extract area of specific frequencies
- spectral envelope: peak frequencies
- CDR: current density reconstruction (GCDR, CDR spectrum), activity within specified band
- EEG bands: medial vs left vs right channels within each band
- MEG data (fT instead of μV)
- merge EEG objects
- reports in .md format
- brain topography
- events markers; epoch by event markers; rewrite epoching
- eeg_keep_eeg_channels → keep_channels(type)
- export channel locs to .CED
- preview 2d/3d channel locs
- swap channel locs axes
- insert channel
- modify channel data
- virtual channels (e.g. F3 + 2.5 × Fp1 - 3 × Cz / 4)
- 3d headplots
- small plots (amplitude, PSD, spectrogram) at electrode locations
- phase-amplitude cross-frequency coupling (PAC)
- more re-referencing methods: Laplacian, REST
- io: import from CSV
- create EEG object
- dipoles
- automated DC line cleaning
- automated channel rejection
- automated cleaning of artifacts
- bad channel marking / rejection
- bad epoch marking / rejection
- time-frequency analysis
- signals/spectra comparison
- more channel location formats
- io: import from EDF+, BDF and other formats
- channel interpolation: manual, automated
- source localization
- ERPs
- ML/DL analysis

NIRS
- import and process data

MRI
- import and process data for EEG source localization

NSTIM
- TES modelling

## Contributors

If you've contributed, add your name below!

[Adam Wysokiński](mailto:adam.wysokinski@umed.lodz.pl)

![umed](images/umed.jpg)

## License

The program is licensed under GPL-2.0-only.
