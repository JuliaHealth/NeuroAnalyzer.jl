# System info

OS: Linux samurai 5.19.0-12.1-liquorix-amd64 #1 ZEN SMP PREEMPT_DYNAMIC liquorix 5.19-14.1~bookworm (2022-09- x86_64 GNU/Linux
CPU: AMD Ryzen 5 PRO 4650U with Radeon Graphics 6 cores (3.9 GHz)
RAM: 16.5 GB free / 22.7 GB

# Loading NeuroAnalyzer.jl

168.603299 seconds (160.27 M allocations: 9.584 GiB, 2.12% gc time, 4.50% compilation time: 73% of which was recompilation)

# NeuroAnalyzer info

NeuroAnalyzer: 0.22.10
        Julia: 1.8.2
         CUDA: not available (use_cuda = false)
 Plugins path: /home/eb/NeuroAnalyzer/plugins/
      Threads: 8 [set using `JULIA_NUM_THREADS` environment variable or Julia --threads command-line option]

Imported packages:
            ColorSchemes 3.19.0
                     CSV 0.10.4
            CubicSplines 0.2.1
                    CUDA 3.12.0
              DataFrames 1.4.1
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
          Interpolations 0.14.6
                    JLD2 0.4.25
                   Loess 0.5.4
       MultivariateStats 0.10.0
                   Plots 1.35.3
             Polynomials 3.2.0
             Preferences 1.3.0
  ScatteredInterpolation 0.3.6
                 Simpson 1.0.1
               StatsFuns 1.0.1
                StatsKit 0.3.1
             StatsModels 0.6.33
              StatsPlots 0.15.4
                Wavelets 0.9.4

# IO

Import EDF+               0.873421 seconds (3.12 M allocations: 501.760 MiB, 53.24% gc time)
Import BDF+               0.050935 seconds (20.57 k allocations: 108.687 MiB)
Import Digitrack          1.075113 seconds (18.82 M allocations: 1.268 GiB, 46.71% gc time)
Import EDF                0.221062 seconds (42.02 k allocations: 152.768 MiB, 41.36% gc time)

# PROCESS

A referencing             0.212459 seconds (14.18 k allocations: 348.200 MiB, 21.39% gc time)
M referencing             0.160856 seconds (23.48 k allocations: 348.658 MiB)
CAR referencing           0.300446 seconds (92.71 k allocations: 918.103 MiB, 30.78% gc time)
Channel referencing       0.043787 seconds (15.72 k allocations: 46.216 MiB)
Filter notch              0.174409 seconds (59.52 k allocations: 587.340 MiB, 34.82% gc time)
Filter LP                 0.305698 seconds (1.03 M allocations: 645.087 MiB, 46.87% gc time)
Filter HP                 0.323890 seconds (1.03 M allocations: 645.370 MiB)

# ANALYZE

Total power               0.311979 seconds (101.50 k allocations: 75.360 MiB)
Band power                0.217968 seconds (136.45 k allocations: 87.179 MiB)
Covariance matrix         0.010494 seconds (3.11 k allocations: 46.106 MiB)
Correlation matrix        0.008414 seconds (3.23 k allocations: 46.130 MiB)
Auto-covariance           0.029747 seconds (45.19 k allocations: 137.325 MiB)
Cross-covariance 1        1.531581 seconds (3.01 M allocations: 10.209 GiB, 10.23% gc time)
Cross-covariance 2        0.000157 seconds (221 allocations: 14.281 KiB)
PSD 1                     0.292055 seconds (97.19 k allocations: 114.948 MiB)
PSD 2                     2.777568 seconds (148.44 k allocations: 2.225 GiB, 2.09% gc time)
Stationarity: mean        0.015504 seconds (36.00 k allocations: 47.573 MiB)
Stationarity: var         0.015599 seconds (38.30 k allocations: 47.888 MiB)
Stationarity: euclid      0.972215 seconds (4.07 M allocations: 3.844 GiB, 16.03% gc time)
Stationarity: hilbert     0.431842 seconds (168.55 k allocations: 279.079 MiB)
Stationarity: adf         2.219382 seconds (399.29 k allocations: 2.515 GiB)
Mutual information 1      1.588979 seconds (518.31 k allocations: 3.147 GiB)
Mutual information 2      2.441301 seconds (967.12 k allocations: 5.977 GiB, 3.99% gc time)
Entropy                   0.094082 seconds (24.48 k allocations: 47.352 MiB)
Negentropy                0.107810 seconds (29.07 k allocations: 92.364 MiB)
Time coherence            0.001168 seconds (220 allocations: 1.184 MiB)
Signal difference 1       0.916701 seconds (239.79 k allocations: 5.726 GiB, 5.56% gc time)
Signal difference 2       1.056316 seconds (323.97 k allocations: 6.404 GiB, 3.71% gc time)
Epoch stats               0.313705 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             1.308073 seconds (107.04 k allocations: 4.182 GiB, 7.74% gc time)
Spectrogram 2             0.465847 seconds (146.07 k allocations: 545.350 MiB)
Spectrum 1                1.449910 seconds (91.90 k allocations: 544.070 MiB, 78.68% gc time)
Spectrum 2                0.444249 seconds (175.78 k allocations: 504.179 MiB)
Channel stats             0.069019 seconds (24.51 k allocations: 180.830 MiB)
SNR                       0.041822 seconds (6.11 k allocations: 581.602 KiB)
Standardize               0.050617 seconds (3.76 k allocations: 134.981 MiB)
Frequency convolution     0.971736 seconds (253.06 k allocations: 3.794 GiB, 6.17% gc time, 0.82% compilation time)
Time convolution          0.660629 seconds (189.60 k allocations: 1.854 GiB, 7.26% gc time, 1.21% compilation time)
DFT                       0.134220 seconds (84.62 k allocations: 319.373 MiB)
MSCI95                    0.010472 seconds (5.65 k allocations: 38.090 MiB)
Mean                      0.017112 seconds (10.03 k allocations: 47.856 MiB)
Difference                0.856256 seconds (239.81 k allocations: 5.726 GiB, 3.22% gc time)
Temporal envelope 1       1.077442 seconds (36.31 M allocations: 5.703 GiB, 19.65% gc time)
Temporal envelope 2       1.178868 seconds (36.31 M allocations: 5.770 GiB, 21.67% gc time)
Temporal envelope 3       1.487479 seconds (41.43 M allocations: 6.109 GiB, 19.13% gc time)
Power envelope 1          0.903409 seconds (8.02 M allocations: 4.994 GiB, 14.04% gc time)
Power envelope 2          0.364430 seconds (523.97 k allocations: 399.650 MiB, 9.15% gc time)
Power envelope 3          0.334374 seconds (1.02 M allocations: 412.618 MiB)
Spectral envelope 1       1.142767 seconds (411.18 k allocations: 5.572 GiB, 6.31% gc time)
Spectral envelope 2       1.174905 seconds (414.01 k allocations: 5.574 GiB, 6.35% gc time)
Spectral envelope 3       1.206400 seconds (481.30 k allocations: 5.576 GiB, 5.98% gc time)
ISPC 1                    6.393993 seconds (3.58 M allocations: 8.086 GiB, 1.90% gc time)
ISPC 2                    0.000642 seconds (275 allocations: 456.562 KiB)
ITPC                      0.019503 seconds (8.02 k allocations: 17.049 MiB)
ITPC spectrogram        453.656068 seconds (254.16 M allocations: 275.508 GiB, 1.69% gc time)
PLI 1                     6.845823 seconds (3.52 M allocations: 6.769 GiB, 1.92% gc time)
PLI 2                     0.000546 seconds (271 allocations: 396.266 KiB)
AEC                       0.033947 seconds (32.41 k allocations: 95.779 MiB)
GED                       0.357195 seconds (8.33 k allocations: 183.233 MiB, 53.36% gc time)
Instant frequency         0.354247 seconds (192.92 k allocations: 504.132 MiB)
Wavelet spectrogram      14.505681 seconds (2.93 M allocations: 40.471 GiB, 10.97% gc time)
TKEO                      0.040383 seconds (10.69 k allocations: 90.460 MiB)
Wavelet spectrum          3.862125 seconds (3.05 M allocations: 10.691 GiB, 2.79% gc time)
Frequency coherence       0.012637 seconds (344 allocations: 3.416 MiB)
F-test                    0.087902 seconds (47.42 k allocations: 4.864 MiB)
F-test                    0.102835 seconds (58.66 k allocations: 4.745 MiB)
Band power                0.277465 seconds (106.90 k allocations: 91.839 MiB)
Relative PSD              0.248962 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split      2.706756 seconds (19.30 M allocations: 6.857 GiB, 14.76% gc time)
Channel difference        0.010074 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   39.698694 seconds (1.92 M allocations: 40.392 GiB, 8.18% gc time)
Cross power spectrum 2    0.006457 seconds (130 allocations: 1.741 MiB)
Amplitude difference      0.063089 seconds (88.83 k allocations: 142.094 MiB)
Virtual channel           2.490425 seconds (2.50 M allocations: 385.070 MiB, 5.37% gc time, 0.01% compilation time)

# PLOT

