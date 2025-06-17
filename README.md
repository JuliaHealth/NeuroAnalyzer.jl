![NeuroAnalyzer.jl](images/neuroanalyzer.png)

[![DOI: 10.5281/zenodo.7372648](images/doi.png)](https://doi.org/10.5281/zenodo.7372648) [![status-badge](https://ci.codeberg.org/api/badges/AdamWysokinski/NeuroAnalyzer.jl/status.svg)](https://ci.codeberg.org/AdamWysokinski/NeuroAnalyzer.jl) [![docs-badge](https://img.shields.io/badge/documentation-blue.svg)](https://neuroanalyzer.org/docs/) [![tuts-badge](https://img.shields.io/badge/tutorials-blue.svg)](https://neuroanalyzer.org#tutorials) [![license-badge](https://img.shields.io/badge/licence-BSD_2C-blue.svg)](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/branch/main/LICENSE) [![JOSS-status](https://joss.theoj.org/papers/cefb2c86708619ddd786b49ce47c98df/status.svg)](https://joss.theoj.org/papers/cefb2c86708619ddd786b49ce47c98df)

NeuroAnalyzer is a [Julia](https://julialang.org) toolbox for analyzing neurophysiological data. Currently it covers importing, editing, processing, visualizing, and analyzing EEG, MEP and EDA data. Preliminary functionality is also available for MEG, NIRS, ECoG, SEEG and iEEG recordings.

Various methods for modeling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS/TUS/INS) will also be implemented (NeuroStim submodule). Another submodule, NeuroTester, will allow designing and running psychological studies. Certain neurophysiological data can be recorded using NeuroRecorder submodule.

NeuroAnalyzer contains a set of separate (high- and low-level) functions. Some interactive graphical user interface (GUI) functions are also available. NeuroAnalyzer functions can be combined into an analysis pipeline, i.e. a Julia script containing all steps of your analysis. This, combined with processing power of Julia language and easiness of distributing calculations across computing cluster, will make NeuroAnalyzer particularly useful for processing large amounts of neurophysiological data.

NeuroAnalyzer is a collaborative, non-commercial project, developed for researchers in psychiatry, neurology and neuroscience.

Every contribution (bug reports, fixes, new ideas, feature requests or additions, documentation improvements, etc.) to the project is highly welcomed.

NeuroAnalyzer website is located at [https://neuroanalyzer.org](https://neuroanalyzer.org).

You may also follow NeuroAnalyzer on [Mastodon](https://fediscience.org/web/tags/neuroanalyzer).

Note: this toolbox is under active development and is subject to change without prior notice.

## Quickstart

Install FIRLSFilterDesign.jl from https://codeberg.org/AdamWysokinski/FIRLSFilterDesign.jl:
```julia
using Pkg
Pkg.add("https://codeberg.org/AdamWysokinski/FIRLSFilterDesign.jl")
```

Start NeuroAnalyzer:
```julia
using NeuroAnalyzer
```

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

## How to Cite

If you use this toolbox, please acknowledge us by citing our [paper](https://neuroanalyzer.org#how-to-cite).

## Contributing

Every contribution (bug reports, fixes, new ideas, feature requests or additions, speed optimization, better code, documentation improvements, typos, etc.) to the project is highly welcomed.

You are very welcome to raise issues and start pull requests. Bugs, suggestions and questions should be reported using the Codeberg [AdamWysokinski/NeuroAnalyzer.jl](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/issues) (preferred method) or Github [JuliaHealth/NeuroAnalyzer.jl](https://github.com/JuliaHealth/NeuroAnalyzer.jl/issues) issues page.

If you notice any bugs, such as crashing code, incorrect results or speed issues, please raise a Codeberg/GitHub issue. Before filing an issue please:

- check that there are no similar existing issues already
- check that your versions are up to date

If you want to report a bug, include your version and system information, and all relevant information. If possible, condense your bug into the shortest example possible that the maintainers can replicate, a so called "minimal working example" or MWE.

If you want to suggest a new feature, for example functionality that other plotting packages offer already, include supplementary material such as example images if possible, so it's clear what you are asking for.

When opening a pull request, please add a short but meaningful description of the changes/features you implemented. Moreover, please add tests (where appropriate) to ensure that your code is working as expected.

For each feature you want to contribute, please file a separate PR to keep the complexity down and time to merge short. Add PRs in draft mode if you want to discuss your approach first.

## Contributors

Below is the list of contributors and their affiliations.

[Adam Wysokiński](mailto:adam.wysokinski@neuroanalyzer.org) [![ORCID](images/orcid.png)](https://orcid.org/0000-0002-6159-6579)

[![Medical University of Lodz](images/umed.png)](https://en.umed.pl)
