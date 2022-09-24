# System info

OS: Linux shogun 5.19.0-1-amd64 #1 SMP PREEMPT_DYNAMIC Debian 5.19.6-1 (2022-09-01) x86_64 GNU/Linux
CPU: AMD Ryzen Threadripper 3960X 24-Core Processor 24 cores (3.71 GHz)
RAM: 97.2 GB free / 125.7 GB

# Loading NeuroAnalyzer.jl

455.710348 seconds (156.27 M allocations: 9.312 GiB, 0.75% gc time, 1.62% compilation time: 78% of which was recompilation)

# NeuroAnalyzer info

NeuroAnalyzer: 0.22.9
        Julia: 1.8.1
         CUDA: CUDA not available (use_cuda = false)
 Plugins path: ~/Documents/NeuroAnalyzer/plugins/
      Threads: 48 [set using using the `JULIA_NUM_THREADS` environment variable]

Imported packages:
            ColorSchemes 3.19.0
                     CSV 0.10.4
            CubicSplines 0.2.1
                    CUDA 3.12.0
              DataFrames 1.3.6
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
          Interpolations 0.14.5
                    JLD2 0.4.23
                   Loess 0.5.4
       MultivariateStats 0.10.0
                   Plots 1.33.0
             Polynomials 3.2.0
             Preferences 1.3.0
  ScatteredInterpolation 0.3.6
                 Simpson 1.0.1
               StatsFuns 1.0.1
                StatsKit 0.3.1
             StatsModels 0.6.32
              StatsPlots 0.15.3
                Wavelets 0.9.4

# IO

Import EDF+               0.359889 seconds (3.12 M allocations: 501.761 MiB, 12.13% gc time)
Import BDF+               0.047690 seconds (20.57 k allocations: 108.687 MiB)
Import EDF                0.084607 seconds (42.02 k allocations: 152.768 MiB)

# PROCESS

A referencing             0.474144 seconds (43.05 k allocations: 350.744 MiB, 8.00% gc time)
M referencing             0.459176 seconds (52.37 k allocations: 351.202 MiB)
CAR referencing           0.376692 seconds (121.17 k allocations: 898.276 MiB)
Channel referencing       0.814932 seconds (44.58 k allocations: 48.760 MiB, 52.64% gc time)
Filter notch              0.379074 seconds (88.43 k allocations: 590.255 MiB)
Filter LP                 0.394065 seconds (1.06 M allocations: 648.003 MiB)
Filter HP                 0.369885 seconds (1.06 M allocations: 648.283 MiB)

# ANALYZE

Total power               0.360574 seconds (121.12 k allocations: 77.844 MiB)
Band power                0.264759 seconds (159.51 k allocations: 89.842 MiB)
Covariance matrix         0.005686 seconds (2.99 k allocations: 46.114 MiB)
Correlation matrix        0.003748 seconds (3.11 k allocations: 46.138 MiB)
Auto-covariance           0.374827 seconds (74.08 k allocations: 140.091 MiB)
Cross-covariance 1        3.465999 seconds (3.10 M allocations: 10.217 GiB, 16.96% gc time)
Cross-covariance 2        0.003112 seconds (460 allocations: 39.578 KiB)
PSD 1                     0.282909 seconds (119.06 k allocations: 117.501 MiB)
PSD 2                     0.882011 seconds (178.35 k allocations: 2.227 GiB)
Stationarity: mean        0.385148 seconds (64.88 k allocations: 50.339 MiB)
Stationarity: var         0.375719 seconds (67.18 k allocations: 50.655 MiB)
Stationarity: euclid      1.187244 seconds (4.13 M allocations: 3.850 GiB, 15.41% gc time)
Stationarity: hilbert     0.348452 seconds (211.02 k allocations: 282.186 MiB)
Stationarity: adf         1.238768 seconds (430.81 k allocations: 2.518 GiB)
Mutual information 1      1.458894 seconds (577.28 k allocations: 3.152 GiB)
Mutual information 2      6.259236 seconds (1.01 M allocations: 5.981 GiB, 84.74% gc time, 0.33% compilation time)
Entropy                   0.369914 seconds (53.39 k allocations: 50.045 MiB)
Negentropy                0.369190 seconds (57.98 k allocations: 95.057 MiB)
Time coherence            0.003386 seconds (460 allocations: 1.207 MiB)
Signal difference 1       0.562314 seconds (269.00 k allocations: 5.728 GiB)
Signal difference 2       0.435702 seconds (353.43 k allocations: 6.406 GiB)
Epoch stats               0.285989 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             0.532784 seconds (138.61 k allocations: 4.185 GiB, 5.74% gc time)
Spectrogram 2             0.278019 seconds (179.44 k allocations: 548.254 MiB)
Spectrum 1                0.391572 seconds (118.30 k allocations: 547.010 MiB)
Spectrum 2                0.377482 seconds (221.99 k allocations: 507.622 MiB)
Channel stats             0.380254 seconds (53.41 k allocations: 183.893 MiB)
SNR                       0.371309 seconds (34.99 k allocations: 3.261 MiB)
Standardize               0.036099 seconds (3.76 k allocations: 134.981 MiB)
Frequency convolution     0.474532 seconds (285.13 k allocations: 3.797 GiB)
Time convolution          0.444950 seconds (217.87 k allocations: 1.857 GiB)
DFT                       0.394111 seconds (109.11 k allocations: 322.107 MiB)
MSCI95                    0.004416 seconds (5.89 k allocations: 38.114 MiB)
Mean                      0.005585 seconds (10.27 k allocations: 47.881 MiB)
Difference                0.409443 seconds (269.09 k allocations: 5.728 GiB)
Temporal envelope 1       1.856054 seconds (36.33 M allocations: 5.705 GiB, 72.55% gc time)
Temporal envelope 2       2.396082 seconds (36.34 M allocations: 5.773 GiB, 83.23% gc time)
Temporal envelope 3       0.864627 seconds (41.46 M allocations: 6.112 GiB, 22.41% gc time)
Power envelope 1          0.408266 seconds (8.05 M allocations: 4.997 GiB)
Power envelope 2          0.361772 seconds (544.75 k allocations: 402.169 MiB)
Power envelope 3          0.456843 seconds (1.05 M allocations: 415.128 MiB, 6.25% gc time)
Spectral envelope 1       0.681840 seconds (445.74 k allocations: 5.575 GiB, 28.59% gc time, 1.23% compilation time)
Spectral envelope 2       0.577879 seconds (443.90 k allocations: 5.577 GiB, 23.86% gc time)
Spectral envelope 3       0.612498 seconds (518.71 k allocations: 5.579 GiB, 26.47% gc time)
ISPC 1                    6.486481 seconds (4.26 M allocations: 8.109 GiB, 5.57% gc time, 0.06% compilation time)
ISPC 2                    0.003025 seconds (515 allocations: 481.266 KiB)
ITPC                      0.017794 seconds (8.02 k allocations: 17.049 MiB)
ITPC spectrogram        409.198798 seconds (273.65 M allocations: 276.091 GiB, 4.07% gc time)
PLI 1                     6.445921 seconds (4.34 M allocations: 6.793 GiB, 3.12% gc time)
PLI 2                     0.003091 seconds (511 allocations: 420.344 KiB)
AEC                       0.018252 seconds (32.88 k allocations: 95.824 MiB)
GED                       0.014344 seconds (7.85 k allocations: 183.229 MiB)
Instant frequency         0.896746 seconds (445.99 k allocations: 510.117 MiB, 38.33% gc time)
Wavelet spectrogram       3.842331 seconds (3.23 M allocations: 40.482 GiB, 5.46% gc time)
TKEO                      0.378701 seconds (39.59 k allocations: 93.153 MiB)
Wavelet spectrum          3.490768 seconds (3.41 M allocations: 10.702 GiB, 11.04% gc time)
Frequency coherence       0.011943 seconds (585 allocations: 3.439 MiB)
F-test                    0.765208 seconds (104.99 k allocations: 10.169 MiB)
F-test                    0.014930 seconds (58.91 k allocations: 4.770 MiB)
Band power                0.347035 seconds (126.26 k allocations: 94.462 MiB)
Relative PSD              0.235329 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split      1.751132 seconds (19.33 M allocations: 6.860 GiB, 45.36% gc time)
Channel difference        0.012384 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   27.135606 seconds (2.48 M allocations: 40.437 GiB, 18.63% gc time)
Cross power spectrum 2    0.011985 seconds (180 allocations: 1.744 MiB, 48.38% compilation time)
Amplitude difference      0.375460 seconds (117.76 k allocations: 144.862 MiB)

# PLOT

