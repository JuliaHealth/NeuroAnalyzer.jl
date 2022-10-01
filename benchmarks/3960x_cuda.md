# System info

OS: Linux shogun 5.19.0-1-amd64 #1 SMP PREEMPT_DYNAMIC Debian 5.19.6-1 (2022-09-01) x86_64 GNU/Linux
CPU: AMD Ryzen Threadripper 3960X 24-Core Processor 24 cores (3.71 GHz)
RAM: 116.1 GB free / 125.7 GB

# Loading NeuroAnalyzer.jl

140.651180 seconds (156.27 M allocations: 9.311 GiB, 2.39% gc time, 5.15% compilation time: 78% of which was recompilation)

# NeuroAnalyzer info

NeuroAnalyzer: 0.22.9
        Julia: 1.8.1
         CUDA: 11.7.0 (use_cuda = true)
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

Import EDF+               0.723181 seconds (3.12 M allocations: 501.761 MiB, 58.90% gc time)
Import BDF+               0.045713 seconds (20.57 k allocations: 108.687 MiB)
Import EDF                0.093580 seconds (42.02 k allocations: 152.768 MiB)

# PROCESS

A referencing             0.458169 seconds (43.05 k allocations: 350.745 MiB, 8.27% gc time)
M referencing             0.418491 seconds (52.37 k allocations: 351.202 MiB)
CAR referencing           0.367360 seconds (121.14 k allocations: 897.219 MiB)
Channel referencing       0.377850 seconds (44.58 k allocations: 48.760 MiB)
Filter notch              0.386256 seconds (88.43 k allocations: 590.255 MiB)
Filter LP                 0.380783 seconds (1.06 M allocations: 648.003 MiB)
Filter HP                 0.365044 seconds (1.06 M allocations: 648.284 MiB)

# ANALYZE

Total power               0.326462 seconds (122.03 k allocations: 77.871 MiB)
Band power                0.266497 seconds (160.24 k allocations: 89.864 MiB)
Covariance matrix         0.004962 seconds (2.99 k allocations: 46.114 MiB)
Correlation matrix        0.004896 seconds (3.11 k allocations: 46.138 MiB)
Auto-covariance           0.362837 seconds (74.08 k allocations: 140.091 MiB)
Cross-covariance 1        2.672607 seconds (3.09 M allocations: 10.217 GiB, 9.31% gc time)
Cross-covariance 2        0.003081 seconds (460 allocations: 39.578 KiB)
PSD 1                     0.285844 seconds (119.68 k allocations: 117.520 MiB)
PSD 2                     0.947419 seconds (184.60 k allocations: 2.228 GiB)
Stationarity: mean        0.363487 seconds (64.88 k allocations: 50.339 MiB)
Stationarity: var         0.367373 seconds (67.18 k allocations: 50.655 MiB)
Stationarity: euclid      1.183409 seconds (4.13 M allocations: 3.850 GiB, 15.02% gc time)
Stationarity: hilbert     0.363666 seconds (212.77 k allocations: 282.239 MiB)
Stationarity: adf         1.284909 seconds (430.69 k allocations: 2.518 GiB)
Mutual information 1      1.496367 seconds (577.38 k allocations: 3.152 GiB)
Mutual information 2      0.887429 seconds (996.80 k allocations: 5.980 GiB)
Entropy                   0.370385 seconds (53.39 k allocations: 50.045 MiB)
Negentropy                0.361351 seconds (57.98 k allocations: 95.057 MiB)
Time coherence            0.014253 seconds (639 allocations: 1002.398 KiB)
Signal difference 1       0.448002 seconds (269.45 k allocations: 5.728 GiB)
Signal difference 2       0.592219 seconds (353.40 k allocations: 6.406 GiB, 26.09% gc time)
Epoch stats               0.270329 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             0.584559 seconds (142.89 k allocations: 4.185 GiB, 6.15% gc time)
Spectrogram 2             0.296749 seconds (180.87 k allocations: 548.297 MiB)
Spectrum 1                0.684177 seconds (340.35 k allocations: 510.024 MiB)
Spectrum 2                0.363555 seconds (222.36 k allocations: 507.634 MiB)
Channel stats             0.379407 seconds (53.41 k allocations: 183.893 MiB)
SNR                       0.380808 seconds (34.99 k allocations: 3.261 MiB)
Standardize               0.036614 seconds (3.76 k allocations: 134.981 MiB)
Frequency convolution   147.062383 seconds (706.30 k allocations: 3.299 GiB, 87.98% gc time)
Time convolution          0.456221 seconds (219.99 k allocations: 1.857 GiB)
DFT                       0.565075 seconds (331.19 k allocations: 285.122 MiB)
MSCI95                    0.006122 seconds (5.89 k allocations: 38.113 MiB)
Mean                      0.006610 seconds (10.28 k allocations: 47.881 MiB)
Difference                0.503199 seconds (269.16 k allocations: 5.728 GiB, 10.62% gc time)
Temporal envelope 1       0.988339 seconds (36.34 M allocations: 5.705 GiB, 60.16% gc time)
Temporal envelope 2       0.399678 seconds (36.34 M allocations: 5.773 GiB)
Temporal envelope 3       0.817910 seconds (41.46 M allocations: 6.112 GiB, 18.75% gc time)
Power envelope 1          0.429682 seconds (8.05 M allocations: 4.997 GiB)
Power envelope 2          0.539263 seconds (549.54 k allocations: 402.315 MiB, 31.29% gc time)
Power envelope 3          0.369625 seconds (1.05 M allocations: 415.187 MiB)
Spectral envelope 1       0.591232 seconds (439.80 k allocations: 5.575 GiB, 22.54% gc time)
Spectral envelope 2       0.451352 seconds (448.49 k allocations: 5.577 GiB)
Spectral envelope 3       0.445133 seconds (510.98 k allocations: 5.579 GiB)
ISPC 1                    6.447711 seconds (4.27 M allocations: 8.109 GiB, 5.65% gc time, 0.06% compilation time)
ISPC 2                    0.002934 seconds (515 allocations: 481.266 KiB)
ITPC                      0.017479 seconds (8.02 k allocations: 17.049 MiB)
ITPC spectrogram        403.328044 seconds (272.27 M allocations: 276.050 GiB, 4.11% gc time)
PLI 1                     5.834424 seconds (3.95 M allocations: 6.786 GiB)
PLI 2                     0.003109 seconds (511 allocations: 420.344 KiB)
AEC                       0.022355 seconds (32.89 k allocations: 95.825 MiB)
GED                       0.014516 seconds (7.85 k allocations: 183.229 MiB)
Instant frequency         0.340658 seconds (265.10 k allocations: 506.351 MiB)
Wavelet spectrogram      34.432400 seconds (6.78 M allocations: 37.463 GiB, 44.78% gc time)
TKEO                      0.372365 seconds (39.59 k allocations: 93.153 MiB)
Wavelet spectrum        549.406397 seconds (6.72 M allocations: 9.608 GiB, 94.38% gc time)
Frequency coherence       0.011956 seconds (585 allocations: 3.439 MiB)
F-test                    0.748960 seconds (104.99 k allocations: 10.169 MiB)
F-test                    0.016166 seconds (58.90 k allocations: 4.770 MiB)
Band power                0.356742 seconds (125.78 k allocations: 94.448 MiB)
Relative PSD              0.238782 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split      2.218011 seconds (19.33 M allocations: 6.860 GiB, 56.89% gc time)
Channel difference        0.012933 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   23.313774 seconds (2.47 M allocations: 40.437 GiB, 8.79% gc time)
Cross power spectrum 2    0.011757 seconds (180 allocations: 1.744 MiB, 49.08% compilation time)
Amplitude difference      0.366064 seconds (117.76 k allocations: 144.862 MiB)

# PLOT

