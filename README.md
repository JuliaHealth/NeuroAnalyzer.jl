![NeuroAnalyzer.jl](images/neuroanalyzer.png)

[![DOI: 10.5281/zenodo.7372648](https://zenodo.org/badge/DOI/10.5281/zenodo.7372648.svg)](https://doi.org/10.5281/zenodo.7372648) [![status-badge](https://ci.codeberg.org/api/badges/AdamWysokinski/NeuroAnalyzer.jl/status.svg)](https://ci.codeberg.org/AdamWysokinski/NeuroAnalyzer.jl)

Welcome fellow researcher!

NeuroAnalyzer is a [Julia](https://julialang.org) package for analyzing of EEG data. Future versions will also process MEG, ECoG, depth electrodes and NIRS data. Also, MRI data will be used for source localization techniques. Stimuli presentation module and various methods for modelling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS) will also be included.

NeuroAnalyzer contains a set of separate (high-level) functions, it does not have a graphical user interface (although one could built it upon these). NeuroAnalyzer functions can be combined into an analysis pipeline, i.e. a Julia script containing all steps of your analysis. This combined with processing power of Julia language and easiness of distributing calculations across computing cluster, will make NeuroAnalyzer particularly useful for processing large amounts of research data.

NeuroAnalyzer is a collaborative, non-commercial project, developed for researchers in psychiatry, neurology and neuroscience.

Every contribution (bug reports, fixes, new ideas, feature requests or additions, documentation improvements, etc.) to the project is highly welcomed.

NeuroAnalyzer website is located at [https://neuroanalyzer.org](https://neuroanalyzer.org).

## Documentation

Complete NeuroAnalyzer (stable branch) documentation is available in [Markdown](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/master/Documentation.md) and [HTML](https://neuroanalyzer.org/docs/index.html) formats.

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

## Known bugs

List of reported bugs is [here](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/issues?labels=70759).

## Contributors

If you have contributed, please add your name below.

[Adam Wysokiński](mailto:adam.wysokinski@umed.lodz.pl)

![Medical University of Lodz](images/umed.png)

## License

The program is licensed under [GPL-2.0-only](LICENSE).