# NeuroJ.jl

Welcome fellow researcher! NeuroJ.jl is a [Julia](https://julialang.org) package for analysis of EEG data. Future versions will also process NIRS data and use MRI data for source localization techniques. Also, various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will be included.

This is a non-commercial projected, aimed for researchers in psychiatry, neurology and neuroscience.

Anyone willing to participate is highly welcomed, open issues, make fork and pull requests.

## Installation

```
using Pkg
Pkg.add(url="https://notabug.org/AdamWysokinski/NeuroJ.jl")
```

## Requirements

- CSV
- DataFrames
- DSP
- FFTW
- Interpolations
- JLD2
- LinearAlgebra
- Pkg
- Plots
- [Simpson](https://notabug.org/AdamWysokinski/Simpson.jl)
- StatsKit

## Roadmap

### To do

.. so much to do ..

General:
- performance optimization
- describe return in __doc__

EEG:
- channel locations data to eeg_signal_header
- time-frequency analysis
- spectrogram
- topoplots
- signals comparison
- kwargs for plots
- plots: fix: time ticks
- fix: eeg_crosscov(edf1, edf2)
- add re-referencing methods
- import channel location files
- source localization
- edit eeg header
- compose eeg_plots with electrode plots
- io: import formats
- io: export to CSV

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