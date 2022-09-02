# System info

OS: Linux shogun 5.18.0-4-amd64 #1 SMP PREEMPT_DYNAMIC Debian 5.18.16-1 (2022-08-10) x86_64 GNU/Linux
CPU: AMD Ryzen Threadripper 3960X 24-Core Processor 24 cores (3.71 GHz)
RAM: 87.6 GB free / 125.7 GB

# Loading NeuroAnalyzer.jl

144.631469 seconds (165.13 M allocations: 9.882 GiB, 2.40% gc time, 5.35% compilation time: 76% of which was recompilation)

# NeuroAnalyzer info

NeuroAnalyzer: 0.22.9
        Julia: 1.8.0
         CUDA: 11.7.0 (use_cuda = false)
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

# IO

Import EDF+               0.712520 seconds (3.12 M allocations: 501.761 MiB, 59.49% gc time)
Import BDF+               0.046347 seconds (20.57 k allocations: 108.687 MiB)
Import EDF                0.082218 seconds (42.02 k allocations: 152.768 MiB)

# PROCESS

CAR referencing           0.454603 seconds (120.99 k allocations: 817.832 MiB, 17.76% gc time)
Channel referencing       0.390714 seconds (44.58 k allocations: 48.760 MiB)
Filter notch              0.391235 seconds (88.43 k allocations: 590.255 MiB)
Filter LP                 0.385744 seconds (1.06 M allocations: 648.003 MiB)
Filter HP                 0.357251 seconds (1.06 M allocations: 648.283 MiB)

# ANALYZE

Total power               0.348468 seconds (121.81 k allocations: 77.865 MiB)
Band power                0.315500 seconds (157.72 k allocations: 89.787 MiB)
Covariance matrix         0.004023 seconds (2.99 k allocations: 46.114 MiB)
Correlation matrix        0.003820 seconds (3.11 k allocations: 46.138 MiB)
Auto-covariance           0.385383 seconds (74.08 k allocations: 140.091 MiB)
Cross-covariance 1        2.433943 seconds (3.09 M allocations: 10.217 GiB, 12.54% gc time)
Cross-covariance 2        0.003206 seconds (460 allocations: 39.578 KiB)
PSD 1                     0.318135 seconds (117.93 k allocations: 117.466 MiB)
PSD 2                     1.083755 seconds (179.51 k allocations: 2.227 GiB, 11.37% gc time)
Stationarity: mean        0.378989 seconds (64.88 k allocations: 50.339 MiB)
Stationarity: var         0.383299 seconds (67.18 k allocations: 50.655 MiB)
Stationarity: euclid      1.215452 seconds (4.13 M allocations: 3.850 GiB, 14.63% gc time)
Stationarity: hilbert     0.358362 seconds (211.97 k allocations: 282.215 MiB)
Stationarity: adf         1.220848 seconds (430.74 k allocations: 2.518 GiB)
Mutual information 1      1.481287 seconds (577.26 k allocations: 3.152 GiB)
Mutual information 2      6.457149 seconds (996.78 k allocations: 5.980 GiB, 72.80% gc time)
Entropy                   0.361813 seconds (53.39 k allocations: 50.045 MiB)
Negentropy                0.372333 seconds (57.98 k allocations: 95.057 MiB)
Time coherence            0.003511 seconds (463 allocations: 1.207 MiB)
Signal difference 1       0.551868 seconds (269.15 k allocations: 5.728 GiB)
Signal difference 2       0.495973 seconds (353.44 k allocations: 6.406 GiB, 9.38% gc time)
Epoch stats               0.273794 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             0.547814 seconds (139.77 k allocations: 4.185 GiB, 5.88% gc time)
Spectrogram 2             0.292407 seconds (180.35 k allocations: 548.282 MiB)
Spectrum 1                0.405266 seconds (120.68 k allocations: 546.872 MiB)
Spectrum 2                0.598339 seconds (220.60 k allocations: 507.580 MiB, 25.55% gc time)
Channel stats             0.361656 seconds (53.41 k allocations: 183.893 MiB)
SNR                       0.373506 seconds (34.99 k allocations: 3.261 MiB)
Standardize               0.034394 seconds (3.76 k allocations: 134.981 MiB)
Frequency convolution     0.443899 seconds (299.67 k allocations: 3.797 GiB)
Time covolution           1.199848 seconds (245.79 k allocations: 1.858 GiB, 50.01% gc time, 3.08% compilation time)
DFT                       0.380401 seconds (111.42 k allocations: 321.967 MiB)
MSCI95                    0.004548 seconds (5.89 k allocations: 38.114 MiB)
Mean                      0.005526 seconds (10.28 k allocations: 47.881 MiB)
Difference                0.476137 seconds (269.33 k allocations: 5.728 GiB)
Temporal envelope 1       0.386171 seconds (36.33 M allocations: 5.705 GiB)
Temporal envelope 2       2.000115 seconds (36.34 M allocations: 5.773 GiB, 80.15% gc time)
Temporal envelope 3       0.861306 seconds (41.46 M allocations: 6.112 GiB, 22.81% gc time)
Power envelope 1          0.415071 seconds (8.05 M allocations: 4.997 GiB)
Power envelope 2          0.453371 seconds (543.72 k allocations: 402.137 MiB, 10.49% gc time)
Power envelope 3          0.470543 seconds (1.04 M allocations: 415.115 MiB, 5.15% gc time)
Spectral envelope 1       0.950116 seconds (449.58 k allocations: 5.575 GiB, 27.84% gc time, 0.78% compilation time)
Spectral envelope 2       0.574481 seconds (453.67 k allocations: 5.577 GiB, 22.11% gc time)
Spectral envelope 3       0.396067 seconds (506.12 k allocations: 5.579 GiB)
ISPC 1                    7.204499 seconds (4.34 M allocations: 8.110 GiB, 9.52% gc time)
ISPC 2                    0.003183 seconds (515 allocations: 481.266 KiB)
ITPC                      0.017457 seconds (7.89 k allocations: 19.375 MiB)
ITPC spectrogram        415.943914 seconds (283.63 M allocations: 534.613 GiB, 4.46% gc time)
PLI 1                     5.987076 seconds (3.93 M allocations: 6.785 GiB)
PLI 2                     0.003064 seconds (511 allocations: 420.344 KiB)
AEC                       0.013930 seconds (32.88 k allocations: 95.825 MiB)
GED                       0.016087 seconds (7.85 k allocations: 183.229 MiB)
Instant frequency         0.366457 seconds (258.18 k allocations: 506.140 MiB)
Wavelet spectrogram       3.753382 seconds (3.34 M allocations: 40.485 GiB, 4.88% gc time)
TKEO                      0.401001 seconds (39.59 k allocations: 93.153 MiB)
Wavelet spectrum          2.518878 seconds (3.24 M allocations: 10.699 GiB)
Frequency coherence       0.012020 seconds (585 allocations: 3.439 MiB)
F-test                    0.755464 seconds (104.99 k allocations: 10.169 MiB)
F-test                    0.014928 seconds (58.91 k allocations: 4.770 MiB)
Band power                0.353072 seconds (126.13 k allocations: 94.458 MiB)
Relative PSD              0.231676 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split     17.382625 seconds (19.33 M allocations: 6.860 GiB, 93.51% gc time)
Channel difference        0.015592 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   21.795510 seconds (2.47 M allocations: 40.437 GiB, 3.92% gc time)
Cross power spectrum 1    0.011690 seconds (180 allocations: 1.744 MiB, 49.36% compilation time)
Amplitude difference      0.374314 seconds (117.76 k allocations: 144.862 MiB)

# PLOT

