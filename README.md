# NeuroJ.jl

Welcome fellow researcher! NeuroJ.jl is a [Julia](https://julialang.org) package for analyzing of EEG data. Future versions will also process NIRS data and use MRI data for source localization techniques. Also, various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will be included.

This is a non-commercial projected, aimed for researchers in psychiatry, neurology and neuroscience.

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

All eeg_* functions will process all channels and epochs of the input EEG object. To process individual channels/epochs, you need to extract them from the EEG object first (`eeg_get_channel()`, `eeg_get_epoch()`)

For `eeg_plot_*()` channels and epochs may be specified.

EEG object (headers + data) is stored in the EEG structure:
```
mutable struct EEG
    eeg_header::Dict
    eeg_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
end
```

## Roadmap

### Bugs

- plots: time ticks
- eeg_crosscov(edf1, edf2)

### To do

.. so much to do ..

General:
- verbose flag
- performance optimization

EEG:
- eeg_picks: left, right, midline, dlpfc, frontal, temporal, etc.
- create EEG object
- clean DC line
- automated channel rejection
- channel interpolation
- automated clean artifacts
- channel locations data to eeg_signal_header
- time-frequency analysis
- spectrogram
- topoplots
- signals comparison
- kwargs for plots
- add re-referencing methods
- import channel location files
- source localization
- edit eeg header
- compose eeg_plots with electrode plots
- io: import EDF+, BDF, other formats
- io: export to CSV
- ERPs

NIRS
- import and process data

MRI
- import and process data for EEG source localization

NSTIM
- modelling

## Contributors

If you've contributed, add your name below!

[Adam Wysokiński](adam.wysokinski@umed.lodz.pl)

## License

The program is licensed under GPL-2.0-only.