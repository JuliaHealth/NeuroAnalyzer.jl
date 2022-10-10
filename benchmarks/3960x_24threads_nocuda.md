# System info

OS: Linux shogun 5.19.0-2-amd64 #1 SMP PREEMPT_DYNAMIC Debian 5.19.11-1 (2022-09-24) x86_64 GNU/Linux
CPU: AMD Ryzen Threadripper 3960X 24-Core Processor 24 cores (3.84 GHz)
RAM: 114.0 GB free / 125.7 GB

# Loading NeuroAnalyzer.jl

249.524751 seconds (158.81 M allocations: 9.550 GiB, 1.40% gc time, 2.65% compilation time: 74% of which was recompilation)

# NeuroAnalyzer info

NeuroAnalyzer: 0.22.10
        Julia: 1.8.2
         CUDA: 11.7.0 (use_cuda = false)
 Plugins path: /home/eb/NeuroAnalyzer/plugins/
      Threads: 24 [set using using the `JULIA_NUM_THREADS` environment variable]

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

Import EDF+               0.749330 seconds (3.12 M allocations: 501.760 MiB, 57.98% gc time)
Import BDF+               0.047110 seconds (20.57 k allocations: 108.687 MiB)
Import Digitrack          0.214298 seconds (18.82 M allocations: 1.268 GiB)
Import EDF                0.087217 seconds (42.02 k allocations: 152.768 MiB)

# PROCESS

A referencing             0.138334 seconds (26.05 k allocations: 349.227 MiB, 21.02% gc time)
M referencing             0.129990 seconds (35.43 k allocations: 349.687 MiB, 20.02% gc time)
CAR referencing           0.145535 seconds (106.28 k allocations: 953.539 MiB)
Channel referencing       0.071164 seconds (27.64 k allocations: 47.244 MiB)
Filter notch              0.161225 seconds (71.66 k allocations: 588.523 MiB)
Filter LP                 0.069480 seconds (1.04 M allocations: 646.272 MiB)
Filter HP                 0.069940 seconds (1.04 M allocations: 646.549 MiB)

# ANALYZE

Total power               0.131550 seconds (119.57 k allocations: 76.664 MiB)
Band power                0.130800 seconds (154.22 k allocations: 88.504 MiB)
Covariance matrix         0.004094 seconds (3.20 k allocations: 46.115 MiB)
Correlation matrix        0.004086 seconds (3.32 k allocations: 46.139 MiB)
Auto-covariance           0.020928 seconds (56.97 k allocations: 138.437 MiB)
Cross-covariance 1        2.806623 seconds (3.04 M allocations: 10.212 GiB, 7.11% gc time)
Cross-covariance 2        0.000307 seconds (317 allocations: 24.406 KiB)
PSD 1                     0.132368 seconds (112.69 k allocations: 116.174 MiB)
PSD 2                     0.809550 seconds (172.30 k allocations: 2.226 GiB)
Stationarity: mean        0.017737 seconds (47.78 k allocations: 48.685 MiB)
Stationarity: var         0.018200 seconds (50.08 k allocations: 49.001 MiB)
Stationarity: euclid      0.810801 seconds (4.09 M allocations: 3.847 GiB, 22.05% gc time)
Stationarity: hilbert     0.223693 seconds (201.83 k allocations: 280.818 MiB)
Stationarity: adf         1.151738 seconds (410.96 k allocations: 2.516 GiB)
Mutual information 1      0.698661 seconds (542.88 k allocations: 3.149 GiB)
Mutual information 2      2.062983 seconds (979.04 k allocations: 5.979 GiB)
Entropy                   0.046454 seconds (36.78 k allocations: 48.451 MiB)
Negentropy                0.058215 seconds (41.35 k allocations: 93.462 MiB)
Time coherence            0.001050 seconds (320 allocations: 1.193 MiB)
Signal difference 1       0.390221 seconds (251.51 k allocations: 5.727 GiB)
Signal difference 2       0.477417 seconds (335.70 k allocations: 6.405 GiB, 21.89% gc time)
Epoch stats               0.270231 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             0.505774 seconds (129.12 k allocations: 4.183 GiB, 6.06% gc time)
Spectrogram 2             0.226029 seconds (171.26 k allocations: 546.872 MiB)
Spectrum 1                0.127421 seconds (112.49 k allocations: 545.402 MiB)
Spectrum 2                0.283749 seconds (214.08 k allocations: 506.161 MiB)
Channel stats             0.030575 seconds (36.75 k allocations: 182.076 MiB)
SNR                       0.081113 seconds (18.33 k allocations: 1.665 MiB)
Standardize               0.053689 seconds (3.77 k allocations: 134.982 MiB)
Frequency convolution     0.381845 seconds (273.85 k allocations: 3.796 GiB)
Time convolution          0.414824 seconds (206.98 k allocations: 1.855 GiB)
DFT                       0.830628 seconds (131.46 k allocations: 322.012 MiB, 66.49% gc time, 4.31% compilation time)
MSCI95                    0.005856 seconds (5.75 k allocations: 38.099 MiB)
Mean                      0.007696 seconds (10.13 k allocations: 47.866 MiB)
Difference                0.445199 seconds (251.44 k allocations: 5.727 GiB)
Temporal envelope 1       0.387171 seconds (36.32 M allocations: 5.704 GiB)
Temporal envelope 2       1.355878 seconds (36.33 M allocations: 5.771 GiB, 70.32% gc time)
Temporal envelope 3       0.852701 seconds (41.44 M allocations: 6.110 GiB, 22.84% gc time)
Power envelope 1          0.653930 seconds (8.05 M allocations: 4.996 GiB, 51.39% gc time)
Power envelope 2          0.322114 seconds (541.88 k allocations: 400.949 MiB, 23.89% gc time)
Power envelope 3          0.219162 seconds (1.04 M allocations: 413.925 MiB)
Spectral envelope 1       0.589598 seconds (429.80 k allocations: 5.573 GiB, 24.65% gc time)
Spectral envelope 2       0.569718 seconds (448.10 k allocations: 5.576 GiB, 21.34% gc time)
Spectral envelope 3       0.524671 seconds (510.55 k allocations: 5.577 GiB, 18.07% gc time)
ISPC 1                    4.136480 seconds (4.35 M allocations: 8.110 GiB, 5.99% gc time, 0.10% compilation time)
ISPC 2                    0.000481 seconds (372 allocations: 466.469 KiB)
ITPC                      0.017495 seconds (8.02 k allocations: 17.049 MiB)
ITPC spectrogram        280.309665 seconds (283.99 M allocations: 276.394 GiB, 4.63% gc time)
PLI 1                     4.249946 seconds (4.36 M allocations: 6.792 GiB, 6.34% gc time)
PLI 2                     0.000534 seconds (367 allocations: 405.891 KiB)
AEC                       0.010169 seconds (32.60 k allocations: 95.798 MiB)
GED                       0.014304 seconds (8.43 k allocations: 183.243 MiB)
Instant frequency         0.804568 seconds (481.39 k allocations: 510.419 MiB, 42.69% gc time)
Wavelet spectrogram       3.886606 seconds (3.36 M allocations: 40.484 GiB, 4.57% gc time)
TKEO                      0.045751 seconds (22.56 k allocations: 91.546 MiB)
Wavelet spectrum          2.523181 seconds (3.82 M allocations: 10.711 GiB, 13.12% gc time)
Frequency coherence       0.011898 seconds (441 allocations: 3.425 MiB)
F-test                    0.054428 seconds (71.84 k allocations: 7.027 MiB)
F-test                    0.023345 seconds (58.76 k allocations: 4.755 MiB)
Band power                0.132292 seconds (125.19 k allocations: 93.209 MiB)
Relative PSD              0.663174 seconds (211.63 k allocations: 155.908 MiB, 42.81% gc time)
Frequency band split      1.662616 seconds (19.32 M allocations: 6.859 GiB, 43.32% gc time, 0.31% compilation time)
Channel difference        0.012661 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   15.503398 seconds (2.12 M allocations: 40.409 GiB, 4.56% gc time)
Cross power spectrum 2    0.005916 seconds (130 allocations: 1.741 MiB)
Amplitude difference      0.038938 seconds (100.91 k allocations: 143.216 MiB)
Virtual channel           0.431798 seconds (2.50 M allocations: 385.078 MiB, 0.02% compilation time)

# PLOT

