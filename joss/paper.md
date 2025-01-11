---
title: 'NeuroAnalyzer: Julia toolbox for analyzing neurophysiological data'

tags:
- Julia
- eeg
- ecog
- ieeg
- seeg
- meg
- nirs
- julia
- neuroscience

authors:
- name: Adam Wysokiński
  orcid: 0000-0002-6159-6579
  affiliation: 1

affiliations:
- name: Medical University of Lodz, Poland
  index: 1

date: 01 November 2024
bibliography: paper.bib
---

# Summary

NeuroAnalyzer is a collaborative, non-commercial Julia toolbox developed for researchers in psychiatry, neurology and neuroscience. It's repository is located at [https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl](https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl) and mirrored at JuliaHealth [https://github.com/JuliaHealth/NeuroAnalyzer.jl](https://github.com/JuliaHealth/NeuroAnalyzer.jl). This software is licensed under The 2-Clause BSD License.

# Background

NeuroAnalyzer is a Julia [@bezanson_2017] toolbox for analyzing neurophysiological (such as EEG, MEG, MEP and EDA) data. To my knowledge, there is no complete Julia solution for neurophysiology data analysis. There are some available packages (some are not actively developed) to perform individual tasks, such as importing EDF/BDF data ([EDF.jl](https://beacon-biosignals.github.io/EDF.jl/stable/), [EDFPlus.jl](https://github.com/wherrera10/EDFPlus.jl), [BDF.jl](https://github.com/sam81/BDF.jl)), processing EEG data ([EEG.jl](https://github.com/JuliaPackageMirrors/EEG.jl), abandoned), ClusterDepth multiple comparison ([ClusterDepth.jl](https://github.com/s-ccs/ClusterDepth.jl)), performing linear/GAM/hierarchical/deconvolution regression on biological signals ([Unfold.jl](https://github.com/unfoldtoolbox/Unfold.jl)) or plotting topomaps ([TopoPlots.jl](https://github.com/MakieOrg/TopoPlots.jl)). Therefore, the development of `NeuroAnalyzer` was started in 2022 in order to provide researchers in the field of neuroscience, psychiatry or neurology a complete, complex (but still easy to use) solution for importing, editing, processing, analyzing and visualizing neurophysiological data.

# Statement of need

There are many excellent MATLAB and Python based EEG/MEG/NIRS applications (e.g. EEGLAB [@delorme_2011], Fieldtrip [@oostenveld_2011], Brainstorm or MNE [@gramfort_2013]). They have been in development for many years and are well established in the scientific community. Many state-of-the-art papers were published using data prepared using these programs. However, compared with Python and MATLAB, there are many advantages of Julia, which underlie my decision to start developing such a toolbox in Julia. I believe that Julia is the future of scientific computing and scientific data analysis [@selvaraj_2022]. Major advantages of Julia are listed in Julia documentation:

- Julia is fast [@bezanson_2018]. In many situations Julia is considerably faster than Python (without having to use numba/cython) and MATLAB. Moreover, Julia provides unlimited scalability. Julia programs can easily be ran on a large cluster or across distributed computers.
- Julia is open-source and free. Increasing MATLAB licensing costs are prohibitive to individual researchers and many research institutions.
- From its very beginning Julia is being focused on scientific computations. Currently only Julia, C, C++ and Fortran belong to the HPC (High Performance Computing) Petaflop Club. Julia is designed for distributed and parallel computations, making it great for distributed analyzes of large data sets.
- Most of the Julia packages are written in pure Julia. It’s easier to understand and modify their code if you already know Julia.
- Julia is beautifully designed, making programming in Julia a pure pleasure. This elegant design makes Julia easy to read and write.

# Three years of development

Currently NeuroAnalyzer includes functions for importing, editing, processing, visualizing, and analyzing EEG, MEP and EDA data. Preliminary functionality is also available for MEG, NIRS, ECoG, SEEG and iEEG recordings.

Various methods for modeling non-invasive brain stimulation protocols (tDCS/tACS/tRNS/tPCS/TMS/TUS/INS) are also being implemented (NeuroStim submodule). Another submodule, NeuroTester, will allow designing and running psychological studies. Certain neurophysiological data can be recorded using NeuroRecorder submodule.

NeuroAnalyzer contains a set of separate (high- and low-level) functions. Some interactive graphical user interface (GUI) functions are also available. NeuroAnalyzer functions can be combined into an analysis pipeline, i.e. a Julia script containing all steps of your analysis. This, combined with processing power of Julia language and easiness of distributing calculations across computing cluster, will make NeuroAnalyzer particularly useful for processing large amounts of neurophysiological data.

Currently NeuroAnalyzer is focused on resting-state analysis. Some ERP functions are already available, while other type of analyses will be developed in future versions. The goal is to make a powerful, expandable and flexible environment for processing and analysis of various types of neurophysiological data.

The following list of already implemented functionalities is presented below, with many more to [come in the future releases](https://neuroanalyzer.org/roadmap.html).

1. Load neurophysiological recordings:
   - EEG (EDF, EDF+, BDF, BDF+, GDF, Alice4, DigiTrack, BrainVision, CSV, EEGLAB, NPY, Thymatron, NCS, CNT, XDF)
   - MEG (FIFF)
   - NIRS (SNIRF, NIRS, NIRX)
   - MEP (DuoMAG)
   - body sensors: acceleration, magnetic field, angular velocity and orientation
   - electrode positions (CED, LOCS, ELC, TSV, SFP, CSD, GEO, MAT, TXT, DAT, ASC)
2. Edit:
   - edit channel data (unit, type, label)
   - edit electrode locations
   - trim (remove part of the signal)
   - resample (up/down)
   - divide into epochs (fixed and by event markers)
   - delete channels/epochs
   - auto-detect bad channels/epochs
   - interpolate channels (planar interpolation/linear regression)
3. Process:
   - reference (common/averaged/auricular/mastoid/Laplacian/custom montage)
   - filter (FIR/IIR/Remez/moving average/moving median/polynomial filters), all types (HP, LP, BP, BS); with preview of filter response
   - remove power line noise
   - auto-detect and remove electrode pops
   - ICA/PCA decompose/reconstruct
   - convolution (in time and frequency domain)
   - create ERP (event-related potentials)
   - NIRS: convert raw light intensity to optical density and HbO/HbR/HbT concentrations
4. Analyze:
   - signal comparison
   - stationarity
   - frequency analysis: total power, band power (absolute and relative)
   - auto- and cross- covariance and correlation (biased and unbiased)
   - time-frequency analysis: various spectrogram methods (FFT-based, short-time Fourier transform, multi-tapered periodogram, Morlet wavelet convolution, Gaussian and Hilbert transform, continuous wavelet transformation)
   - coherence and magnitude-squared coherence
   - mutual information
   - entropy, negentropy
   - envelopes (amplitude, power, spectrogram)
   - power spectrum slope
   - PLI/ISPC/ITPC
   - ERP: detect peaks, analyze amplitude and average amplitude
   - EROs (event-related oscillations): spectrogram, power spectrum
   - HRV (heart rate variability): time-domain analysis (MENN, MDNN, VNN, SDNN, RMSSD, SDSD, NN50, pNN50, NN20, pNN20)
   - MEPs: detect peaks, analyze amplitude and average amplitude
5. Plotting:
   - signal (single-/multi-channel)
   - power spectrum (single-/multi-channel, 2D/3D)
   - spectrogram (single-/multi-channel)
   - topographical maps (various methods of interpolation)
   - weights at channel locations
   - weighted inter-channel connections
   - matrices (channel × channel)
   - heatmaps
   - channel/epoch data (histogram/bar plot/dot plot/box plot/violin plot/polar plot/paired data)
   - ERPs: amplitude, topographical distribution
   - EROs: spectrogram, power spectrum
   - MEPs: amplitude

Interactive (Gtk3-based) plotting (signal, power spectrum, spectrograms, topomaps, ICA components) and editing (signal, electrode positions) are also implemented.

All computations are performed using the double-precision 64-bit floating point format. NeuroAnalyzer data is stored using standard Julia Array and can be easily exported as DataFrame. Thus, external processing of those data using Julia packages is readily available.

NeuroAnalyzer also includes NeuroRecorder, a set of functions for recording various neurophysiological signals:

1. Finger Tapping Test (FTT) -- using computer keyboard or external panel attached to Raspberry Pi
2. Electrodermal Activity (EDA) = Galvanic Skin Response (GSR) -- via Raspberry Pi
3. Two-point Pinch Test (TPT) -- using finger-worn accelerator attached to Raspberry Pi (in development)
4. Angular Velocity Sensors (AVS) -- via Raspberry Pi (in development)

Extensive documentation and tutorials are available at [https://www.neuroanalyzer.org/docs](https://www.neuroanalyzer.org/docs) and [https://neuroanalyzer.org/#tutorials](https://neuroanalyzer.org/#tutorials), respectively.

NeuroAnalyzer functionality can be easily expanded using [plugins](https://neuroanalyzer.org/tut-plugins.html) (written in Julia).

For common tasks (importing, filtering, referencing) NeuroAnalyzer performance is comparable with MNE and ~4× higher than EEGLAB. Benchmarks are available at [https://neuroanalyzer.org/benchmarks.html](https://neuroanalyzer.org/benchmarks.html).

# Research projects using the software

NeuroAnalyzer has been used in the preparation of the following publications:

1. Wysokiński A, Szczepocka E, Szczakowska A. Improved cognitive performance, increased theta, alpha, beta and decreased delta powers after cognitive rehabilitation augmented with tDCS in a patient with post-COVID-19 cognitive impairment (brain-fog). Psychiatry Research Case Reports. 2023 DOI: 10.1016/j.psycr.2023.100164
2. Sochal M. et al. The relationship between sleep quality measured by polysomnography and selected neurotrophic factors. Journal of Clinical Medicine. 2024 DOI: 10.3390/jcm13030893
3. Sochal M. et al. Circadian Rhythm Genes and Their Association with Sleep and Sleep Restriction. International Journal of Molecular Sciences. 2024 DOI: 10.3390/ijms251910445
4. Dataseris G, Zelko J. Physiological signal analysis and open science using the Julia language and associated software. Frontiers in Network Physiology. 2024 DOI: 10.3389/fnetp.2024.1478280
5. Wysokiński A, Pazdrak M. Complete resolution of auditory verbal hallucinations (AVH) in a patient with schizophrenia after transcranial direct current stimulation (tDCS) therapy. Psychiatria i Psycholologia Kliniczna. 2024 (in press)

In the Department of Old Age Psychiatry and Psychotic Disorders (Medical University of Lodz, Poland) we are using NeuroAnalyzer for analyzing EEG data to analyze EEG data from our research projects on functional (EEG-based) connectivity in schizophrenia and EEG markers of tDCS stimulation.

# References
