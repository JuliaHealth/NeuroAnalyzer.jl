![neuroj](images/neuroj.png)

Welcome fellow researcher!

NeuroJ.jl is a [Julia](https://julialang.org) package for analyzing of EEG data. Future versions will also process NIRS data and use MRI data for source localization techniques. Also, various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will be included.

This is a non-commercial projected, aimed for researchers in psychiatry, neurology and neuroscience.

Initially NeuroJ.jl will be focused on resting-state EEG analysis, but ERP and other type of analyses will be developed in future versions.

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
- Interpolations
- JLD2
- LinearAlgebra
- MultivariateStats
- Pkg
- Plots
- [Simpson](https://codeberg.org/AdamWysokinski/Simpson.jl)
- StatsKit

NeuroJ.jl will be 100% Julia based.

## General remarks

NeuroJ.jl will process both NeuroJ.EEG objects (EEG metadata header + time + epoched signals + components) and signals (single-channel or multi-channel/multi-trial). The latter have no header, therefore some functions will be limited.

The following conventions are used:

- single-channel signals and time `Vector{Float64}`
- multi-channel signals `Array{Float64, 3}` (channels × signals × epochs)

If epochs are not defined, the whole signal is an epoch, i.e. there is always at least one epoch.

Functions name prefix:

- `signal_` functions taking single-/multi-channel signals as an argument
- `eeg_` functions taking EEG object as an argument

The majority of `eeg_` functions will process all channels and epochs of the input EEG object. To process individual channels/epochs, you need to extract them from the EEG object first (`eeg_keep_epoch()`, `eeg_keep_channel()` to process as NeuroJ.EEG object or `eeg_extract_channel()`, `eeg_extract_epoch()` to process as multi-channel array)

`eeg_` and `signal_` functions use named arguments for all arguments other than input signal(s).

EEG object (headers + data) is stored in the EEG structure:
```
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

- epoch markers drawn too low in eeg_plot_ica() for few ICs
- channel labels should be displayed as strings not vector
- eeg_plot_topo() check minimum len value for frequency to analyze
- ignore non-eeg channels for processing, analysis and plotting; currently NeuroJ does not analyze/process/plot EEG containing non-eeg channels

## To do

.. so much to do ..

General:
- performance optimization: CUDA/AMD ROCm acceleration

EEG:
- remove embedded components that are not useful
- rewrite plotting functions to be more modular
- piecewise detrending
- store any calculations (e.g. median delta power) as a component for topo plots
- export channel locs to .CED
- preview 2d/3d channel locs
- swap channel locs axes
- insert channel
- modify channel data
- view channel parameters
- virtual channels (e.g. F3 + 2.5 × Fp1 - 3 × Cz / 4)
- 3d headplots
- amplitude plots at electrode locations
- negative index channel to invert polarity
- phase-amplitude cross-frequency coupling (PAC)
- plot spectrogram/psd: use embedded spectrogram/psd
- more re-referencing methods: Laplacian, REST
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