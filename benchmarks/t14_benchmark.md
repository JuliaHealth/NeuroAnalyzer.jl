# System info

OS: Linux samurai 5.19.0-2.1-liquorix-amd64 #1 ZEN SMP PREEMPT_DYNAMIC liquorix 5.19-2.1~bookworm (2022-08-1 x86_64 GNU/Linux
CPU: AMD Ryzen 5 PRO 4650U with Radeon Graphics 6 cores (2.59 GHz)
RAM: 12.2 GB free / 22.7 GB

# Loading NeuroAnalyzer.jl

167.338015 seconds (158.41 M allocations: 9.456 GiB, 1.95% gc time, 4.98% compilation time: 72% of which was recompilation)

# NeuroAnalyzer info

NeuroAnalyzer: 0.22.9
        Julia: 1.8.1
         CUDA: CUDA not available (use_cuda = false)
 Plugins path: ~/Documents/NeuroAnalyzer/plugins/
      Threads: 12 [set using using the `JULIA_NUM_THREADS` environment variable]

Imported packages:
            ColorSchemes 3.19.0
                     CSV 0.10.4
            CubicSplines 0.2.1
                    CUDA 3.12.0
              DataFrames 1.3.4
           Deconvolution 1.1.1
               Distances 0.10.7
                     DSP 0.7.7
                    FFTW 1.5.0
                  FileIO 1.15.0
             FindPeaks1D 0.1.7
                     Git 1.2.1
                     GLM 1.8.0
                 GLMakie 0.6.13
         HypothesisTests 0.10.10
     InformationMeasures 0.3.1
          Interpolations 0.14.4
                    JLD2 0.4.22
                   Loess 0.5.4
       MultivariateStats 0.10.0
                   Plots 1.31.7
             Polynomials 3.2.0
             Preferences 1.3.0
  ScatteredInterpolation 0.3.6
                 Simpson 1.0.1
               StatsFuns 1.0.1
                StatsKit 0.3.1
             StatsModels 0.6.31
              StatsPlots 0.15.2
                Wavelets 0.9.4

# IO

Import EDF+               0.442573 seconds (3.12 M allocations: 501.761 MiB, 17.49% gc time)
Import BDF+               0.083266 seconds (20.57 k allocations: 108.687 MiB, 45.81% gc time)
Import EDF                0.117755 seconds (42.02 k allocations: 152.768 MiB)

# PROCESS

A referencing             0.228555 seconds (17.15 k allocations: 348.457 MiB, 20.21% gc time)
M referencing             0.651851 seconds (26.43 k allocations: 348.913 MiB, 68.83% gc time)
CAR referencing           0.413617 seconds (96.21 k allocations: 933.965 MiB)
Channel referencing       0.042605 seconds (18.64 k allocations: 46.471 MiB)
Filter notch              0.108656 seconds (62.46 k allocations: 587.633 MiB)
Filter LP                 0.157521 seconds (1.04 M allocations: 645.382 MiB)
Filter HP                 0.525180 seconds (1.04 M allocations: 645.663 MiB, 40.87% gc time)

# ANALYZE

Total power               0.227678 seconds (107.76 k allocations: 75.739 MiB)
Band power                0.245401 seconds (142.38 k allocations: 87.556 MiB)
Covariance matrix         0.007630 seconds (2.77 k allocations: 46.094 MiB)
Correlation matrix        0.008355 seconds (2.89 k allocations: 46.118 MiB)
Auto-covariance           0.028092 seconds (48.11 k allocations: 137.602 MiB)
Cross-covariance 1        2.669378 seconds (3.02 M allocations: 10.210 GiB, 6.46% gc time)
Cross-covariance 2        0.000141 seconds (244 allocations: 16.781 KiB)
PSD 1                     0.230457 seconds (103.12 k allocations: 115.317 MiB)
PSD 2                     2.141663 seconds (153.63 k allocations: 2.225 GiB, 2.56% gc time, 0.40% compilation time)
Stationarity: mean        0.015900 seconds (38.90 k allocations: 47.849 MiB)
Stationarity: var         0.017273 seconds (41.22 k allocations: 48.166 MiB)
Stationarity: euclid      1.007061 seconds (4.08 M allocations: 3.845 GiB, 14.31% gc time)
Stationarity: hilbert     0.399091 seconds (179.58 k allocations: 279.596 MiB)
Stationarity: adf         3.042321 seconds (418.98 k allocations: 2.516 GiB, 50.37% gc time, 1.36% compilation time)
Mutual information 1      2.634189 seconds (523.78 k allocations: 3.147 GiB, 28.90% gc time)
Mutual information 2      2.793190 seconds (970.13 k allocations: 5.978 GiB, 24.02% gc time)
Entropy                   0.078501 seconds (27.47 k allocations: 47.624 MiB)
Negentropy                0.081822 seconds (32.07 k allocations: 92.636 MiB)
Time coherence            0.001863 seconds (244 allocations: 1.186 MiB)
Signal difference 1       1.053608 seconds (242.78 k allocations: 5.726 GiB, 2.46% gc time)
Signal difference 2       1.177142 seconds (327.00 k allocations: 6.404 GiB, 3.85% gc time)
Epoch stats               0.330338 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             1.890032 seconds (113.67 k allocations: 4.182 GiB, 1.33% gc time)
Spectrogram 2             0.637818 seconds (151.92 k allocations: 545.717 MiB)
Spectrum 1                0.576259 seconds (97.15 k allocations: 544.470 MiB)
Spectrum 2                0.481780 seconds (187.71 k allocations: 504.747 MiB)
Channel stats             0.061660 seconds (27.49 k allocations: 181.139 MiB)
SNR                       0.041975 seconds (9.03 k allocations: 858.195 KiB)
Standardize               0.051013 seconds (3.76 k allocations: 134.981 MiB)
Frequency convolution     1.343188 seconds (241.95 k allocations: 3.794 GiB)
Time covolution           0.818917 seconds (212.29 k allocations: 1.855 GiB, 9.76% gc time, 3.43% compilation time)
DFT                       0.129535 seconds (90.92 k allocations: 319.790 MiB)
MSCI95                    0.013778 seconds (5.67 k allocations: 38.092 MiB)
Mean                      0.014698 seconds (10.05 k allocations: 47.859 MiB)
Difference                0.981083 seconds (242.77 k allocations: 5.726 GiB, 6.21% gc time)
Temporal envelope 1       1.428575 seconds (36.31 M allocations: 5.703 GiB, 28.22% gc time)
Temporal envelope 2       1.477108 seconds (36.32 M allocations: 5.770 GiB, 29.50% gc time)
Temporal envelope 3       1.504397 seconds (41.43 M allocations: 6.110 GiB, 15.48% gc time)
Power envelope 1          0.973321 seconds (8.03 M allocations: 4.994 GiB, 19.55% gc time)
Power envelope 2          0.466470 seconds (541.71 k allocations: 400.267 MiB, 13.96% gc time)
Power envelope 3          0.361907 seconds (1.03 M allocations: 413.006 MiB, 6.14% gc time)
Spectral envelope 1       1.229856 seconds (414.16 k allocations: 5.572 GiB, 10.14% gc time)
Spectral envelope 2       1.160051 seconds (422.61 k allocations: 5.574 GiB, 4.70% gc time)
Spectral envelope 3       1.160387 seconds (481.14 k allocations: 5.576 GiB, 4.53% gc time)
ISPC 1                    6.137096 seconds (3.90 M allocations: 8.094 GiB, 2.94% gc time)
ISPC 2                    0.000681 seconds (299 allocations: 459.031 KiB)
ITPC                      0.019340 seconds (8.02 k allocations: 17.049 MiB)
ITPC spectrogram        459.215772 seconds (275.76 M allocations: 276.149 GiB, 1.81% gc time)
PLI 1                     6.499148 seconds (3.84 M allocations: 6.777 GiB, 1.95% gc time)
PLI 2                     0.000632 seconds (295 allocations: 398.672 KiB)
AEC                       0.032786 seconds (32.45 k allocations: 95.783 MiB)
GED                       0.049309 seconds (7.63 k allocations: 183.206 MiB)
Instant frequency         0.478529 seconds (212.02 k allocations: 504.578 MiB, 22.46% gc time)
Wavelet spectrogram      15.362512 seconds (2.96 M allocations: 40.472 GiB, 10.49% gc time)
TKEO                      0.071680 seconds (13.46 k allocations: 90.725 MiB)
Wavelet spectrum          3.951890 seconds (3.13 M allocations: 10.694 GiB, 2.84% gc time)
Frequency coherence       0.012607 seconds (368 allocations: 3.418 MiB)
F-test                    0.082113 seconds (53.00 k allocations: 5.388 MiB)
F-test                    0.170110 seconds (58.68 k allocations: 4.747 MiB)
Band power                0.258858 seconds (112.56 k allocations: 92.214 MiB)
Relative PSD              0.589568 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split      3.165563 seconds (19.31 M allocations: 6.858 GiB, 10.74% gc time)
Channel difference        0.009096 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   30.460245 seconds (1.99 M allocations: 40.396 GiB, 3.42% gc time)
Cross power spectrum 2    0.014931 seconds (180 allocations: 1.744 MiB, 58.28% compilation time)
Amplitude difference      0.051570 seconds (91.80 k allocations: 142.373 MiB)

# PLOT

