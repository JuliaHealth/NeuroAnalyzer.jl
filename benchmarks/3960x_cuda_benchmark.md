# NeuroAnalyzer benchmarks

OS: OS: Linux shogun 5.18.0-4-amd64 #1 SMP PREEMPT_DYNAMIC Debian 5.18.16-1 (2022-08-10) x86_64 GNU/Linux
CPU: AMD Ryzen Threadripper 3960X 24-Core Processor 24 cores (3.71 GHz)
RAM: 91.5 GB free /  125.7 GB

NeuroAnalyzer: 0.22.9
        Julia: 1.8.0
         CUDA: 11.7.0 (use_cuda = true)
 Plugins path: ~/Documents/NeuroAnalyzer/plugins/
      Threads: 48 [set using using the `JULIA_NUM_THREADS` environment variable]

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

## IO

Import EDF+               0.691072 seconds (3.12 M allocations: 498.707 MiB, 60.58% gc time)
Import BDF+               0.045072 seconds (20.57 k allocations: 107.749 MiB)
Import EDF                0.096466 seconds (42.01 k allocations: 150.405 MiB)

## PROCESS

CAR referencing           0.379549 seconds (121.19 k allocations: 903.472 MiB)
Channel referencing       0.402788 seconds (44.59 k allocations: 48.779 MiB)
Filter notch              0.394333 seconds (88.43 k allocations: 590.275 MiB)
Filter LP                 0.409550 seconds (1.06 M allocations: 648.022 MiB)
Filter HP                 0.373810 seconds (1.06 M allocations: 648.303 MiB)

## ANALYZE

Total power               0.351514 seconds (121.75 k allocations: 77.863 MiB)
Band power                0.350810 seconds (155.94 k allocations: 89.733 MiB)
Covariance matrix         0.003812 seconds (2.98 k allocations: 46.114 MiB)
Correlation matrix        0.004015 seconds (3.10 k allocations: 46.138 MiB)
Auto-covariance           0.367425 seconds (74.08 k allocations: 140.091 MiB)
Cross-covariance 1        1.796288 seconds (3.09 M allocations: 10.217 GiB, 15.92% gc time)
Cross-covariance 2        0.003028 seconds (460 allocations: 39.656 KiB)
PSD 1                     0.256914 seconds (120.73 k allocations: 117.552 MiB)
PSD 2                     1.204085 seconds (204.80 k allocations: 2.228 GiB, 15.60% gc time)
Stationarity: mean        0.382429 seconds (64.88 k allocations: 50.339 MiB)
Stationarity: var         0.382483 seconds (67.18 k allocations: 50.655 MiB)
Stationarity: euclid      1.185207 seconds (4.13 M allocations: 3.850 GiB, 15.04% gc time)
Stationarity: hilbert     0.357371 seconds (212.16 k allocations: 282.221 MiB)
Stationarity: adf         6.025229 seconds (445.05 k allocations: 2.518 GiB, 77.15% gc time, 0.42% compilation time)
Mutual information 1      1.498560 seconds (577.29 k allocations: 3.152 GiB)
Mutual information 2      4.985739 seconds (996.79 k allocations: 5.980 GiB, 71.27% gc time)
Entropy                   0.374307 seconds (53.39 k allocations: 50.045 MiB)
Negentropy                0.372636 seconds (57.98 k allocations: 95.057 MiB)
Time coherence            0.005045 seconds (639 allocations: 1002.398 KiB)
Signal difference 1       0.556682 seconds (269.27 k allocations: 5.728 GiB, 17.35% gc time)
Signal difference 2       0.512975 seconds (353.40 k allocations: 6.406 GiB, 8.71% gc time)
Epoch stats               0.283592 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             0.586758 seconds (142.40 k allocations: 4.185 GiB, 5.95% gc time)
Spectrogram 2             0.288011 seconds (179.22 k allocations: 548.247 MiB)
Spectrum 1                0.680585 seconds (342.60 k allocations: 510.163 MiB)
Spectrum 2                0.833932 seconds (282.59 k allocations: 509.331 MiB, 20.14% gc time)
Channel stats             0.373977 seconds (53.41 k allocations: 183.893 MiB)
SNR                       0.376106 seconds (34.99 k allocations: 3.261 MiB)
Standardize               0.035420 seconds (3.77 k allocations: 135.001 MiB)
Frequency convolution   101.087896 seconds (702.80 k allocations: 3.299 GiB, 81.28% gc time)
Time covolution           0.484145 seconds (217.83 k allocations: 1.857 GiB)
DFT                       0.658777 seconds (377.65 k allocations: 286.824 MiB, 9.33% gc time, 2.58% compilation time)
MSCI95                    0.004643 seconds (5.89 k allocations: 38.113 MiB)
Mean                      0.005605 seconds (10.28 k allocations: 47.881 MiB)
Difference                0.724110 seconds (299.22 k allocations: 5.729 GiB, 25.57% gc time)
Temporal envelope 1       1.474432 seconds (36.33 M allocations: 5.903 GiB, 69.92% gc time)
Temporal envelope 2       1.402232 seconds (36.34 M allocations: 5.970 GiB, 69.56% gc time)
Temporal envelope 3       0.909740 seconds (41.45 M allocations: 6.349 GiB, 22.71% gc time)
Power envelope 1          0.505883 seconds (8.05 M allocations: 4.997 GiB, 19.38% gc time)
Power envelope 2          0.489712 seconds (545.52 k allocations: 402.192 MiB, 19.15% gc time)
Power envelope 3          0.403497 seconds (1.05 M allocations: 415.219 MiB, 10.34% gc time)
Spectral envelope 1       1.116070 seconds (449.55 k allocations: 5.575 GiB, 34.75% gc time)
Spectral envelope 2       0.625019 seconds (446.81 k allocations: 5.577 GiB)
Spectral envelope 3       0.463734 seconds (510.96 k allocations: 5.579 GiB)
ISPC 1                    5.815651 seconds (3.99 M allocations: 8.102 GiB)
ISPC 2                    0.003038 seconds (515 allocations: 481.266 KiB)
ITPC                      0.018102 seconds (7.89 k allocations: 19.375 MiB)
ITPC spectrogram        404.225747 seconds (285.33 M allocations: 534.661 GiB, 4.25% gc time)
PLI 1                     6.733993 seconds (4.42 M allocations: 6.795 GiB, 7.06% gc time)
PLI 2                     0.003164 seconds (511 allocations: 420.344 KiB)
AEC                       0.022458 seconds (32.88 k allocations: 96.039 MiB)
GED                       0.014864 seconds (7.85 k allocations: 183.229 MiB)
Instant frequency         0.341650 seconds (257.39 k allocations: 506.116 MiB)
Wavelet spectrogram      34.520191 seconds (6.88 M allocations: 43.256 GiB, 45.19% gc time)
TKEO                      0.395439 seconds (39.59 k allocations: 93.153 MiB)
Wavelet spectrum        474.296258 seconds (6.72 M allocations: 9.608 GiB, 93.51% gc time)
Frequency coherence       0.012308 seconds (585 allocations: 3.439 MiB)
F-test                    0.785218 seconds (104.99 k allocations: 10.169 MiB)
F-test                    0.017667 seconds (58.91 k allocations: 4.770 MiB)
Band power                0.362437 seconds (125.80 k allocations: 94.448 MiB)
Relative PSD              0.366889 seconds (211.63 k allocations: 155.908 MiB, 23.20% gc time)
Frequency band split      2.621714 seconds (19.33 M allocations: 6.860 GiB, 63.36% gc time)
Channel difference        0.015943 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   25.281561 seconds (2.47 M allocations: 40.437 GiB, 15.84% gc time)
Cross power spectrum 1    0.011687 seconds (180 allocations: 1.744 MiB, 49.28% compilation time)
Amplitude difference      0.378285 seconds (117.76 k allocations: 144.862 MiB)

## PLOT

