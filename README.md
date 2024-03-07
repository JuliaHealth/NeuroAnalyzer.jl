![NeuroAnalyzer.jl](images/neuroanalyzer.png)

[![DOI: 10.5281/zenodo.7372648](images/doi.png)](https://doi.org/10.5281/zenodo.10718848) [![status-badge](https://ci.codeberg.org/api/badges/AdamWysokinski/NeuroAnalyzer.jl/status.svg)](https://ci.codeberg.org/AdamWysokinski/NeuroAnalyzer.jl) [![docs-badge](https://img.shields.io/badge/documentation-blue.svg)](https://neuroanalyzer.org/docs/) [![tuts-badge](https://img.shields.io/badge/tutorials-blue.svg)](https://neuroanalyzer.org#tutorials) [![license-badge](https://img.shields.io/badge/licence-BSD_2C-blue.svg)](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/branch/main/LICENSE) 

NeuroAnalyzer is a [Julia](https://julialang.org) toolbox for analyzing neurophysiological data. Currently it allows importing, editing, processing and analyzing EEG, MEP and NIRS data. Preliminary functionality is also available for EDA, ECoG, SEEG and iEEG recordings. Future versions will also support MEG recordings and source localization techniques.

Various methods for modeling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS/TUS/INS) will also be implemented ([NeuroStim](https://codeberg.org/AdamWysokinski/NeuroStim.jl) submodule). Another submodule, [NeuroTester](https://codeberg.org/AdamWysokinski/NeuroTester.jl), will allow designing and running psychological studies. Certain neurophysiological data can be recorded using [NeuroRecorder](https://codeberg.org/AdamWysokinski/NeuroRecorder.jl) submodule.

NeuroAnalyzer contains a set of separate (high- and low-level) functions. Some interactive graphical user interface (GUI) functions are also available. NeuroAnalyzer functions can be combined into an analysis pipeline, i.e. a Julia script containing all steps of your analysis. This, combined with processing power of Julia language and easiness of distributing calculations across computing cluster, will make NeuroAnalyzer particularly useful for processing large amounts of neurophysiological data.

NeuroAnalyzer is a collaborative, non-commercial project, developed for researchers in psychiatry, neurology and neuroscience.

Every contribution (bug reports, fixes, new ideas, feature requests or additions, documentation improvements, etc.) to the project is highly welcomed.

NeuroAnalyzer website is located at [https://neuroanalyzer.org](https://neuroanalyzer.org).

Note: this toolbox is under active development and is subject to change without prior notice.

## Quickstart

Add NeuroAnalyzer from the [Pkg REPL](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode): `pkg> add NeuroAnalyzer`.

## Documentation

Documentation is available in the following formats:

- [HTML](https://neuroanalyzer.org/docs)
- [Markdown](https://codeberg.org/AdamWysokinski/NeuroAnalyzer-docs/src/branch/main/Documentation.md)

Changelog and commit details are available at [https://neuroanalyzer.org/changelog.html](https://neuroanalyzer.org/changelog.html).

## Tutorials

NeuroAnalyzer tutorials are available at [https://neuroanalyzer.org#tutorials](https://neuroanalyzer.org#tutorials).

## Requirements

See [https://neuroanalyzer.org/requirements.html](https://neuroanalyzer.org/requirements.html) for more details.

## What's next

This [roadmap](https://neuroanalyzer.org/roadmap.html) of the future developments of NeuroAnalyzer is neither complete, nor in any particular order.

## Performance

For testing performance between individual machines, a [complete set of benchmarks](https://codeberg.org/AdamWysokinski/NeuroAnalyzer-benchmarks) is available.

## Plugins (extensions)

See [https://neuroanalyzer.org/plugins.html](https://neuroanalyzer.org/plugins.html) for more details.

## License

This software is licensed under [The 2-Clause BSD License](LICENSE).

## Financial support

If you would like to support the project financially, we have the Liberapay account:
<a href="https://liberapay.com/~1829183/donate"><img alt="Donate using Liberapay" src="https://liberapay.com/assets/widgets/donate.svg"></a>

## Contributors

Below is the list of contributors and their affiliations.

[Adam Wysokiński](mailto:adam.wysokinski@neuroanalyzer.org) [![ORCID](images/orcid.png)](https://orcid.org/0000-0002-6159-6579)

[![Medical University of Lodz](images/umed.png)](https://en.umed.pl)
