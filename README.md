# NeuroJ.jl

Welcome fellow researcher! NeuroJ.jl is a [Julia](https://julialang.org) package for analyzing of EEG data. Future versions will also process NIRS data and use MRI data for source localization techniques. Also, various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will be included.

This is a non-commercial projected, aimed for researchers in psychiatry, neurology and neuroscience.

Initially NeuroJ.jl will be focused on resting-state EEG analysis, but ERP analysis will developed in future versions.

Every contribution (bug reports, fixes, new ideas, feature requests or additions, documentation improvements, etc.) to the project is highly welcomed.

## Installation

```
using Pkg
Pkg.add(url="https://notabug.org/AdamWysokinski/NeuroJ.jl")
```

### Requirements

- CSV
- DataFrames
- Distances
- DSP
- FFTW
- Interpolations
- JLD2
- LinearAlgebra
- MultivariateStats
- Pkg
- Plots
- [Simpson](https://notabug.org/AdamWysokinski/Simpson.jl)
- StatsKit

NeuroJ will be 100% Julia based.

## General remarks

NeuroJ.jl processes both EEG objects (EEG epoched signals + header) and signals (single-channel or multi-channel/multi-trial). These have no header, therefore some functions will be limited.

The following conventions are used:

- single-channel signals and time      :: `Vector{Float64}`
- multi-channel or multi-trial signals :: `Array{Float64, 3}` channels/trials × signals × epochs

If epochs are not defined, the whole signal is an epoch, i.e. there is always at least one epoch.

Functions name prefix:

- `signal_`  :: functions taking single-/multi-channel signals as an argument
- `eeg_`     :: functions taking EEG object as an argument

All `eeg_` functions will process all channels and epochs of the input EEG object. To process individual channels/epochs, you need to extract them from the EEG object first (`eeg_extract_channel()`, `eeg_extract_epoch()`)

`eeg_` and `signal_` functions use named arguments for all arguments other than input signal(s).

For `eeg_plot_*()` channels and epochs may be specified.

EEG object (headers + data) is stored in the EEG structure:
```
mutable struct EEG
    eeg_header::Dict
    eeg_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
    eeg_components::Vector{Any}
end
```

Many `eeg_` functions have a mutator variant (e.g. `eeg_delete_epoch!()`). These functions modifies the input object 

## Roadmap

### Bugs

### To do

.. so much to do ..

General:
- performance optimization

EEG:

- plot spectrogram/psd/ica: use embedded spectrogram/psd/ica
- topographical plots
- more re-referencing methods
- improve: compose eeg_plots with electrode plots
- io: import from CSV
- create EEG object
- automated DC line cleaning
- automated channel rejection
- automated cleaning of artifacts
- bad channel marking / rejection
- bad epoch marking / rejection
- time-frequency analysis
- signals comparison
- more channel location formats
- io: import from EDF+, BDF and other formats
- channel interpolation
- source localization
- ERPs
- CUDA/AMD ROCm acceleration

NIRS
- import and process data

MRI
- import and process data for EEG source localization

NSTIM
- TES modelling

## Contributors

If you've contributed, add your name below!

[Adam Wysokiński](adam.wysokinski@umed.lodz.pl)

## License

The program is licensed under GPL-2.0-only.