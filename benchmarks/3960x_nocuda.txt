# System info

OS: Linux shogun 6.0.0-2-amd64 #1 SMP PREEMPT_DYNAMIC Debian 6.0.3-1 (2022-10-21) x86_64 GNU/Linux
CPU: AMD Ryzen Threadripper 3960X 24-Core Processor 24 cores (3.66 GHz)
RAM: 113.6 GB free / 125.7 GB

# Loading NeuroAnalyzer.jl

170.769677 seconds (196.65 M allocations: 11.715 GiB, 2.44% gc time, 8.21% compilation time: 46% of which was recompilation)

# NeuroAnalyzer info

    NeuroAnalyzer: 0.22.10
            Julia: 1.8.2
             CUDA: 11.7.0 (use_cuda = false)
     Plugins path: /home/eb/NeuroAnalyzer/plugins/
Show progress bar: true
          Verbose: true
          Threads: 24 [set using `JULIA_NUM_THREADS` environment variable or Julia --threads command-line option]

Imported packages:
            ColorSchemes 3.19.0
                     CSV 0.10.7
            CubicSplines 0.2.1
                    CUDA 3.12.0
              DataFrames 1.4.2
           Deconvolution 1.1.1
               Distances 0.10.7
                     DSP 0.7.7
                    FFTW 1.5.0
                  FileIO 1.16.0
             FindPeaks1D 0.1.7
                     Git 1.2.1
                     GLM 1.8.1
                 GLMakie 0.7.1
         HypothesisTests 0.10.11
     InformationMeasures 0.3.1
          Interpolations 0.13.6
                    JLD2 0.4.25
                   Loess 0.5.4
       MultivariateStats 0.10.0
                   Plots 1.35.5
             Polynomials 3.2.0
             Preferences 1.3.0
           ProgressMeter 1.7.2
  ScatteredInterpolation 0.3.6
                 Simpson 1.0.1
               StatsFuns 1.0.1
                StatsKit 0.3.1
             StatsModels 0.6.33
              StatsPlots 0.15.4
                Wavelets 0.9.5
             WaveletsExt 0.2.1
      ContinuousWavelets 1.1.1

# IO

Import EDF+                       0.369811 seconds (3.12 M allocations: 501.760 MiB, 13.55% gc time)
Import BDF+                       0.070117 seconds (20.57 k allocations: 108.687 MiB, 51.86% gc time)
Import Digitrack                  0.327285 seconds (18.82 M allocations: 1.268 GiB, 41.35% gc time)
Import EDF                        0.101252 seconds (42.02 k allocations: 152.768 MiB)

# EDIT

Delete channel                    0.063947 seconds (200 allocations: 87.479 MiB, 74.06% gc time)
Keep channel                      0.009379 seconds (575 allocations: 47.316 MiB)
Delete epoch                      0.103334 seconds (139 allocations: 89.472 MiB, 84.39% gc time)
Keep epoch                        0.007875 seconds (171 allocations: 45.316 MiB)
Virtual channel                   0.425623 seconds (2.50 M allocations: 385.071 MiB, 0.02% compilation time)

# PROCESS

A referencing                     0.127524 seconds (26.09 k allocations: 349.228 MiB, 23.64% gc time)
M referencing                     0.097635 seconds (35.41 k allocations: 349.686 MiB)
CAR referencing                   0.134995 seconds (106.22 k allocations: 953.537 MiB)
Channel referencing               0.074392 seconds (27.53 k allocations: 47.241 MiB)
Filter notch                      0.165784 seconds (76.34 k allocations: 588.841 MiB)
Filter LP                         0.070062 seconds (1.05 M allocations: 646.583 MiB)
Filter HP                         0.076694 seconds (1.05 M allocations: 646.863 MiB)

# ANALYZE

Total power                       0.148286 seconds (120.05 k allocations: 76.679 MiB)
Band power                        0.142437 seconds (154.72 k allocations: 88.520 MiB)
Covariance matrix                 0.004314 seconds (3.20 k allocations: 46.115 MiB)
Correlation matrix                0.005779 seconds (3.32 k allocations: 46.139 MiB)
Auto-covariance                   0.047684 seconds (57.15 k allocations: 138.443 MiB)
Cross-covariance 1                3.351507 seconds (3.04 M allocations: 10.212 GiB, 11.84% gc time)
Cross-covariance 2                0.000209 seconds (317 allocations: 24.406 KiB)
PSD 1                             0.117417 seconds (113.61 k allocations: 116.203 MiB)
PSD 2                             0.795581 seconds (170.87 k allocations: 2.226 GiB)
Stationarity: mean                0.018017 seconds (47.83 k allocations: 48.687 MiB)
Stationarity: var                 0.018650 seconds (50.11 k allocations: 49.002 MiB)
Stationarity: euclid              0.808090 seconds (4.09 M allocations: 3.847 GiB, 22.49% gc time)
Stationarity: hilbert             0.222616 seconds (201.24 k allocations: 280.800 MiB)
Stationarity: adf                 1.378489 seconds (411.01 k allocations: 2.516 GiB)
Mutual information 1              0.627359 seconds (542.66 k allocations: 3.149 GiB)
Mutual information 2              2.098574 seconds (979.07 k allocations: 5.979 GiB)
Entropy                           0.064699 seconds (45.92 k allocations: 138.544 MiB)
Negentropy                        0.079850 seconds (50.38 k allocations: 183.472 MiB)
Time coherence                    0.000931 seconds (320 allocations: 1.193 MiB)
Signal difference 1               0.337141 seconds (251.53 k allocations: 5.727 GiB)
Signal difference 2               1.008554 seconds (335.77 k allocations: 6.405 GiB, 58.19% gc time)
Epoch stats                       0.273187 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1                     0.495469 seconds (127.64 k allocations: 4.183 GiB, 8.04% gc time)
Spectrogram 2                     0.222780 seconds (169.79 k allocations: 546.828 MiB)
Spectrum 1                        0.107416 seconds (134.06 k allocations: 524.109 MiB)
Spectrum 2                        0.257573 seconds (229.72 k allocations: 506.472 MiB)
Channel stats                     0.035285 seconds (36.76 k allocations: 182.076 MiB)
SNR                               0.073549 seconds (18.33 k allocations: 1.665 MiB)
Standardize                       0.036289 seconds (3.77 k allocations: 134.982 MiB)
Frequency convolution             0.388068 seconds (272.99 k allocations: 3.795 GiB)
Time convolution                  0.410715 seconds (208.35 k allocations: 1.856 GiB)
DFT                               0.088422 seconds (102.85 k allocations: 320.574 MiB)
MSCI95                            0.003764 seconds (5.75 k allocations: 38.099 MiB)
Mean                              0.004494 seconds (10.13 k allocations: 47.866 MiB)
Difference                        0.367220 seconds (251.52 k allocations: 5.727 GiB)
Temporal envelope 1               2.279270 seconds (36.32 M allocations: 5.704 GiB, 73.80% gc time)
Temporal envelope 2               1.337660 seconds (36.32 M allocations: 5.771 GiB, 70.63% gc time)
Temporal envelope 3               0.851759 seconds (41.43 M allocations: 6.110 GiB, 22.09% gc time)
Power envelope 1                  0.260661 seconds (8.04 M allocations: 4.995 GiB)
Power envelope 2                  0.167681 seconds (534.57 k allocations: 400.789 MiB)
Power envelope 3                  0.213667 seconds (1.04 M allocations: 413.790 MiB)
Spectral envelope 1               0.579758 seconds (437.73 k allocations: 5.574 GiB, 24.68% gc time, 1.27% compilation time)
Spectral envelope 2               0.504421 seconds (431.51 k allocations: 5.575 GiB, 23.62% gc time)
Spectral envelope 3               0.337595 seconds (493.85 k allocations: 5.577 GiB)
Hilbert amplitude envelope 1      0.387936 seconds (36.32 M allocations: 5.704 GiB)
Hilbert amplitude envelope 2      1.396390 seconds (36.32 M allocations: 5.771 GiB, 71.64% gc time)
Hilbert amplitude envelope 3      0.860412 seconds (41.43 M allocations: 6.110 GiB, 19.65% gc time)
ISPC 1                            4.087702 seconds (4.20 M allocations: 8.982 GiB, 6.41% gc time, 0.10% compilation time)
ISPC 2                            0.000572 seconds (376 allocations: 506.562 KiB)
ITPC                              0.019466 seconds (8.26 k allocations: 19.418 MiB)
ITPC spectrogram                269.721754 seconds (282.58 M allocations: 309.151 GiB, 4.70% gc time)
PLI 1                             3.524386 seconds (4.11 M allocations: 7.665 GiB)
PLI 2                             0.000537 seconds (372 allocations: 446.016 KiB)
AEC                               0.018209 seconds (32.60 k allocations: 95.798 MiB)
GED                               0.013713 seconds (8.43 k allocations: 183.243 MiB)
Instant frequency                 0.181573 seconds (236.18 k allocations: 550.326 MiB)
Wavelet spectrogram               2.834780 seconds (3.42 M allocations: 40.486 GiB, 4.66% gc time)
TKEO                              0.042847 seconds (22.70 k allocations: 91.550 MiB)
Wavelet spectrum                  1.695696 seconds (3.42 M allocations: 10.703 GiB)
Frequency coherence               0.011933 seconds (441 allocations: 3.425 MiB)
F-test                            0.053591 seconds (71.88 k allocations: 7.028 MiB)
F-test                            0.023307 seconds (58.76 k allocations: 4.755 MiB)
Band power                        0.122535 seconds (124.24 k allocations: 93.180 MiB)
Relative PSD                      0.233425 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split              1.614749 seconds (19.36 M allocations: 6.861 GiB, 46.93% gc time)
Channel difference                0.012887 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1           24.256289 seconds (2.20 M allocations: 40.411 GiB, 38.67% gc time)
Cross power spectrum 2            0.006265 seconds (130 allocations: 1.741 MiB)
Amplitude difference              0.051309 seconds (100.63 k allocations: 143.207 MiB)
DWT                               0.296247 seconds (90.04 k allocations: 1.714 GiB)
CWT                               5.060684 seconds (3.51 M allocations: 57.340 GiB, 17.05% gc time, 0.33% compilation time)
PSD slope                         0.143714 seconds (585.86 k allocations: 290.779 MiB)
