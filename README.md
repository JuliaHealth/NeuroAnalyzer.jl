![NeuroAnalyzer.jl](images/neuroanalyzer.png)

[![DOI: 10.5281/zenodo.7372648](images/doi.png)](https://doi.org/10.5281/zenodo.7372648) [![status-badge](https://ci.codeberg.org/api/badges/AdamWysokinski/NeuroAnalyzer.jl/status.svg)](https://ci.codeberg.org/AdamWysokinski/NeuroAnalyzer.jl) [![docs-badge](https://img.shields.io/badge/documentation-blue.svg)](https://neuroanalyzer.org/docs/) [![tuts-badge](https://img.shields.io/badge/tutorials-blue.svg)](https://neuroanalyzer.org#tutorials) [![license-badge](https://img.shields.io/badge/licence-BSD_2C-blue.svg)](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/src/branch/main/LICENSE) [![JOSS-status](https://joss.theoj.org/papers/cefb2c86708619ddd786b49ce47c98df/status.svg)](https://joss.theoj.org/papers/cefb2c86708619ddd786b49ce47c98df)

NeuroAnalyzer is a high-performance [Julia](https://julialang.org) toolbox for exploring, visualizing, and analyzing neurophysiological data - including EEG, MEG, ECoG, SEEG, iEEG, NIRS, MEP, and EDA recordings.

The toolbox provides both high-level and low-level functions, along with interactive GUI components for exploratory work. Individual functions can be chained into reproducible analysis pipelines - Julia scripts that capture every step of your workflow. Combined with Julia's computational performance and native support for distributed computing, this makes NeuroAnalyzer especially well-suited for processing large neurophysiological datasets across computing clusters.

NeuroAnalyzer is a free, non-commercial, collaborative project built for researchers in neuroscience, neurology, and psychiatry.

NeuroAnalyzer website is located at [https://neuroanalyzer.org](https://neuroanalyzer.org).

You may also follow NeuroAnalyzer on [Mastodon](https://fediscience.org/web/tags/neuroanalyzer).

Note: this toolbox is under active development and its API is subject to changes.

## Quickstart

Start NeuroAnalyzer:

```julia
using NeuroAnalyzer
```

## Documentation

Complete [API reference](https://neuroanalyzer.org/docs) is available.

Changelog and commit details are available at [https://neuroanalyzer.org/changelog.html](https://neuroanalyzer.org/changelog.html).

## Tutorials

NeuroAnalyzer tutorials are available at [https://neuroanalyzer.org#tutorials](https://neuroanalyzer.org#tutorials).

## Requirements

See [Requirements](https://neuroanalyzer.org/requirements.html) for more details.

## What's next

This [roadmap](roadmap.html) gives an overview of where NeuroAnalyzer is headed. It is not exhaustive, and the order of items does not reflect implementation priority.

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

If NeuroAnalyzer contributed to your research, please cite the [JOSS article →](https://doi.org/10.21105/joss.07734):

```{bibtex}
@article{Wysokiński_2025,
    doi = {10.21105/joss.07734},
    url = {https://doi.org/10.21105/joss.07734},
    year = {2025},
    publisher = {The Open Journal},
    volume = {10},
    number = {107},
    pages = {7734},
    author = {Adam Wysokiński},
    title = {NeuroAnalyzer: Julia toolbox for analyzing neurophysiological data},
    journal = {Journal of Open Source Software}
}
```

## Contributing

Every contribution (bug reports, fixes, new ideas, feature requests or additions, speed optimization, better code, documentation improvements, typos, etc.) to the project is highly welcomed.

You are very welcome to raise issues and start pull requests. Bugs, suggestions and questions should be reported using the Codeberg [AdamWysokinski/NeuroAnalyzer.jl](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/issues) (preferred method) or Github [JuliaHealth/NeuroAnalyzer.jl](https://github.com/JuliaHealth/NeuroAnalyzer.jl/issues) issues page.

If you notice any bugs, such as crashing code, incorrect results or speed issues, please raise a Codeberg/GitHub issue. Before filing an issue please:

- check that there are no similar existing issues already
- check that your version is up to date

If you want to report a bug, include your version and system information, and all relevant information. If possible, condense your bug into the shortest example possible that the maintainers can replicate, a so called "minimal working example" or MWE.

If you want to suggest a new feature, for example functionality that other plotting packages offer already, include supplementary material such as example images if possible, so it's clear what you are asking for.

When opening a pull request, please add a short but meaningful description of the changes/features you implemented. Moreover, please add tests (where appropriate) to ensure that your code is working as expected.

For each feature you want to contribute, please file a separate PR to keep the complexity down and time to merge short. Add PRs in draft mode if you want to discuss your approach first.

## Contributors

Below is the list of contributors and their affiliations.

[Adam Wysokiński](mailto:adam.wysokinski@neuroanalyzer.org) [![ORCID](images/orcid.png)](https://orcid.org/0000-0002-6159-6579)

[![Medical University of Lodz](images/umed.png)](https://en.umed.pl)
