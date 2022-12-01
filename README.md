![](images/neuroanalyzer.png)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7372648.svg)](https://doi.org/10.5281/zenodo.7372648) 

Welcome fellow researcher!

NeuroAnalyzer is a [Julia](https://julialang.org) package for analyzing of EEG data. Future versions will also process MEG, ECoG, depth electrodes and NIRS data. Also, it will use MRI data for source localization techniques. Various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will also be included.

NeuroAnalyzer contains a set of separate (high-level) functions, it does not have a graphical user interface (although one could built it upon these). NeuroAnalyzer functions can be combined into an analysis pipeline, i.e. a Julia script containing all steps of your analysis. This combined with processing power of Julia language and easiness of distributing calculations across computing cluster, will make NeuroAnalyzer particularly useful for processing large amounts of research data.

NeuroAnalyzer is a collaborative, non-commercial project, developed for researchers in psychiatry, neurology and neuroscience.

Currently NeuroAnalyzer is focused on resting-state EEG analysis. ERP and other type of analyses will be developed in future versions. The goal is to make a powerful, expandable and flexible environment for EEG/MEG/NIRS/NIBS processing workflows.

Every contribution (bug reports, fixes, new ideas, feature requests or additions, documentation improvements, etc.) to the project is highly welcomed.

NeuroAnalyzer website is located at [https://neuroanalyzer.org](https://neuroanalyzer.org).

## Documentation

Complete NeuroAnalyzer documentation is available in [Markdown](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/master/Documentation.md) and [HTML](https://neuroanalyzer.org/docs/index.html) formats.

Changelog and commit details are available at [https://neuroanalyzer.org](https://neuroanalyzer.org/changelog.html).

## Tutorials

NeuroAnalyzer tutorials are available at [https://neuroanalyzer.org](https://neuroanalyzer.org#tutorials).

## Requirements

See [https://neuroanalyzer.org](https://neuroanalyzer.org/requirements.html) for more details.

## What's next

This [roadmap](https://neuroanalyzer.org/roadmap.html) of the future developments of NeuroAnalyzer is neither complete, nor in any particular order.

## Performance

For testing performance between individual machines, a complete set of benchmarks is available [here](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/master/Benchmarking.md).

## Plugins (extensions)

See [https://neuroanalyzer.org](https://neuroanalyzer.org/plugins.html) for more details.

## Known bugs

List of reported bugs is [here](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/issues?labels=70759).

## Contributors

If you have contributed, please add your name below.

[Adam Wysokiński](mailto:adam.wysokinski@umed.lodz.pl)

![umed](images/umed.png)

## License

The program is licensed under [GPL-2.0-only](LICENSE).