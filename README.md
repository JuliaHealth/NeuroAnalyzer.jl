![NeuroAnalyzer.jl](images/neuroanalyzer.png)

[![DOI: 10.5281/zenodo.7686075](https://zenodo.org/badge/DOI/10.5281/zenodo.7686075.svg)](https://doi.org/10.5281/zenodo.7372648) [![status-badge](https://ci.codeberg.org/api/badges/AdamWysokinski/NeuroAnalyzer.jl/status.svg)](https://ci.codeberg.org/AdamWysokinski/NeuroAnalyzer.jl) [![docs-badge](https://img.shields.io/badge/docs-stable-blue.svg)](https://neuroanalyzer.org/docs-stable/) [![docs-badge](https://img.shields.io/badge/docs-devel-blue.svg)](https://neuroanalyzer.org/docs-devel/)

Welcome fellow researcher!

NeuroAnalyzer is a [Julia](https://julialang.org) toolbox for analyzing neurophysiological data. Currently it allows importing, processing and analyzing EEG and NIRS data. Future versions will also support MEG, ECoG and depth electrodes recordings. Also, it will use MRI data for source localization techniques. Various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will also be included.

NeuroAnalyzer contains a set of separate (high-level) functions, it does not have a graphical user interface (although one could built it upon these). NeuroAnalyzer functions can be combined into an analysis pipeline, i.e. a Julia script containing all steps of your analysis. This combined with processing power of Julia language and easiness of distributing calculations across computing cluster, will make NeuroAnalyzer particularly useful for processing large amounts of research data.

NeuroAnalyzer is a collaborative, non-commercial project, developed for researchers in psychiatry, neurology and neuroscience.

Every contribution (bug reports, fixes, new ideas, feature requests or additions, documentation improvements, etc.) to the project is highly welcomed.

NeuroAnalyzer website is located at [https://neuroanalyzer.org](https://neuroanalyzer.org).

Note: this toolbox is under active development and is subject to change.

## Documentation

1. [Stable](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/branch/stable) branch:
    - [Markdown](https://codeberg.org/AdamWysokinski/NeuroAnalyzer-docs/src/branch/stable/Documentation-stable.md)
    - [HTML](https://neuroanalyzer.org/docs-stable) 
2. [Devel](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/branch/devel) branch:
    - [Markdown](https://codeberg.org/AdamWysokinski/NeuroAnalyzer-docs/src/branch/stable/Documentation-devel.md)
    - [HTML](https://neuroanalyzer.org/docs-devel)

Changelog and commit details are available at [https://neuroanalyzer.org/changelog.html](https://neuroanalyzer.org/changelog.html).

## Tutorials

NeuroAnalyzer tutorials are available at [https://neuroanalyzer.org#tutorials](https://neuroanalyzer.org#tutorials).

## Requirements

See [https://neuroanalyzer.org/requirements.html](https://neuroanalyzer.org/requirements.html) for more details.

## What's next

This [roadmap](https://neuroanalyzer.org/roadmap.html) of the future developments of NeuroAnalyzer is neither complete, nor in any particular order.

## Performance

For testing performance between individual machines, a [complete set of benchmarks](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/master/Benchmarking.md) is available.

## Plugins (extensions)

See [https://neuroanalyzer.org/plugins.html](https://neuroanalyzer.org/plugins.html) for more details.

## Contributors

If you have contributed, please add your name below.

[Adam Wysokiński](mailto:adam.wysokinski@umed.lodz.pl)

![Medical University of Lodz](images/umed.png)

## License

The program is licensed under [GPL-2.0-only](LICENSE).