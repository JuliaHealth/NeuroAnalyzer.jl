![](images/neuroanalyzer.png)

Welcome fellow researcher!

NeuroAnalyzer is a [Julia](https://julialang.org) package for analyzing of EEG data. Future versions will also process MEG, ECoG, depth electrodes and NIRS data. Also, it will use MRI data for source localization techniques. Various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will also be included.

NeuroAnalyzer contains a set of separate (high-level) functions, it does not have a graphical user interface (although one could built it upon these). NeuroAnalyzer functions can be combined into an analysis pipeline, i.e. a Julia script containing all steps of your analysis. This combined with processing power of Julia language and easiness of distributing calculations across computing cluster, will make NeuroAnalyzer particularly useful for processing large amounts of research data.

NeuroAnalyzer is a collaborative, non-commercial project, developed for researchers in psychiatry, neurology and neuroscience.

Currently NeuroAnalyzer is focused on resting-state EEG analysis. ERP and other type of analyses will be developed in future versions. The goal is to make a powerful, expandable and flexible environment for EEG/MEG/NIRS/NIBS processing workflows.

Every contribution (bug reports, fixes, new ideas, feature requests or additions, documentation improvements, etc.) to the project is highly welcomed.

NeuroAnalyzer website is located at [https://neuroanalyzer.org](https://neuroanalyzer.org)

## Installation

First, download [Julia](https://julialang.org/downloads/) 1.7.0 or later. 

There are two branches of NeuroAnalyzer:
- [stable](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/branch/master): released once per month, recommended for research tasks
- [devel](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/branch/devel): a rolling release for NeuroAnalyzer developers, not for production use

You can add NeuroAnalyzer using Julia package manager, by typing:

```Julia
using Pkg
Pkg.update()
# for master branch:
Pkg.add(url="https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl#master")
# for development branch:
Pkg.add(url="https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl#devel")
# activate the package
using NeuroAnalyzer
# check if correctly installed
na_info()
```

Another option is to initialize a new Julia environment for the package:
```shell
git clone https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl
cd NeuroAnalyzer.jl
julia --project
```

Next, in Julia REPL do the following:
```Julia
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
# if necessary:
Pkg.resolve()
Pkg.update()
# activate the package
using NeuroAnalyzer
# check if NeuroAnalyzer has been correctly installed
na_info()
```

## Requirements

See [https://neuroanalyzer.org/requirements.html](https://neuroanalyzer.org/requirements.html) for more details.

## General remarks

NeuroAnalyzer functions operate on NeuroAnalyzer objects.

For EEG this is NeuroAnalyzer.EEG (EEG metadata header + time + epoched signals + components + annotations). EEG signal is `Array{Float64, 3}` (channels × signals × epochs). If epochs are not defined, the whole signal is an epoch, i.e. there is always at least one epoch.

Functions name prefix:
- `eeg_` functions taking EEG object as an argument
- `eeg_plot_` plotting functions

There are also low level functions operating on single-/multi-channel signal vector/matrix/array (`s_` and `s2_`).

The majority of `eeg_` functions will process all channels and epochs of the input EEG object. To process individual channels/epochs, you need to extract them from the EEG object first (`eeg_keep_epoch()`, `eeg_keep_channel()` to process as NeuroAnalyzer.EEG object or `eeg_extract_channel()`, `eeg_extract_epoch()` to process as multi-channel array).

`eeg_` functions use named arguments for all arguments other than input signal(s), e.g. `eeg_delete_epoch!(my_eeg, epoch=12)`.

EEG object (headers + time + epochs time + EEG signal + (optional) components + annotations) is stored in the EEG structure:
```julia
mutable struct EEG
    eeg_header::Dict
    eeg_time::Vector{Float64}
    eeg_epochs_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
    eeg_components::Vector{Any}
    eeg_annotations::DataFrame
end
```

Study object is stored in the STUDY structure:
```julia
mutable struct STUDY
    study_header::Dict{Symbol, Any}
    study_eeg::Vector{NeuroAnalyzer.EEG}
    study_group::Vector{Symbol}
end
```
Certain `eeg_` functions use multiple dispatch mechanism to analyze STUDY object inter- and intra- groups.

Many `eeg_` functions have a mutator variant (e.g. `eeg_delete_epoch!()`). These functions modifies the input EEG object in-place, e.g. you may use `eeg_delete_channel!(my_eeg, channel=1)` instead of `my_eeg = eeg_delete_channel(my_eeg, channel=1)`.

## Preferences

To modify the plugins path, use `na_set_plugins_path("/new/path/to/plugins")` to set the variable `plugins_path`.

CUDA is disabled by default. For some low-level operations (e.g. FFT and IFFT) CUDA acceleration is may be used if compatible NVIDIA card and drivers are installed. To enable CUDA, use `na_set_cuda(true)` to set the variable `use_cuda`.

Most NeuroAnalyzer function take less than one second to complete. Functions that take more than few seconds (e.g. `eeg_itpc_s()`) may optionally display progress bar. To disable the progress bar, use `na_set_progress_bar(false)` to set the variable `progress_bar`.

Preferences are saved to `LocalPreferences.toml` in the local NeuroAnalyzer folder or in `~/.julia/environments/v*/LocalPreferences.toml`. Julia session needs to be restarted to use new preferences.

## Documentation

Complete NeuroAnalyzer documentation is available in [Markdown](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/master/Documentation.md) and [HTML](https://neuroanalyzer.org/docs/index.html).

Tutorial introducing NeuroAnalyzer functions is [here](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/master/Tutorial.md). More tutorials and FAQ are available at [https://neuroanalyzer.org](https://neuroanalyzer.org#tutorials)

Changelog and commit details are available at [https://neuroanalyzer.org/changelog.htm](https://neuroanalyzer.org/changelog.html).

## What's next

This [roadmap](https://neuroanalyzer.org/roadmap.html) of the future developments of NeuroAnalyzer is not complete and not in any particular order.

## Performance

Due to the nature of Julia JIT (just-in-time compiler), first run of functions is slow. Also, due to modules precompilation, NeuroAnalyzer startup may take some time.

For testing performance between individual machines, a complete set of benchmarks is available [here](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/master/Benchmarking.md).

## Plugins (extensions)

See [https://neuroanalyzer.org/plugins.html](https://neuroanalyzer.org/plugins.html) for more details.

## Known bugs

List of reported bugs is [here](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/issues?labels=70759).

## Contributors

If you have contributed, please add your name below.

[Adam Wysokiński](mailto:adam.wysokinski@umed.lodz.pl)

![umed](images/umed.png)

## License

The program is licensed under [GPL-2.0-only](LICENSE).