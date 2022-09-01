# NeuroAnalyzer benchmarks

OS: Linux shogun 5.18.0-4-amd64 #1 SMP PREEMPT_DYNAMIC Debian 5.18.16-1 (2022-08-10) x86_64 GNU/Linux
CPU: AMD Ryzen Threadripper 3960X 24-Core Processor 24 cores (3.71 GHz)
RAM: 90.8 GB free /  125.7 GB

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

## IO

Import EDF+               0.417975 seconds (3.12 M allocations: 498.707 MiB, 23.23% gc time)
Import BDF+               0.076875 seconds (20.57 k allocations: 107.749 MiB, 45.05% gc time)
Import EDF                0.099286 seconds (42.01 k allocations: 150.405 MiB)

## PROCESS

CAR referencing           0.453473 seconds (120.96 k allocations: 808.391 MiB, 16.30% gc time)
Channel referencing       0.389278 seconds (44.59 k allocations: 48.779 MiB)
Filter notch              0.389280 seconds (88.43 k allocations: 590.275 MiB)
Filter LP                 0.388503 seconds (1.06 M allocations: 648.022 MiB)
Filter HP                 0.366686 seconds (1.06 M allocations: 648.303 MiB)

## ANALYZE

Total power               0.370409 seconds (121.12 k allocations: 77.843 MiB)
Band power                0.364485 seconds (155.72 k allocations: 89.726 MiB)
Covariance matrix         0.004073 seconds (2.98 k allocations: 46.114 MiB)
Correlation matrix        0.005016 seconds (3.10 k allocations: 46.138 MiB)
Auto-covariance           0.391305 seconds (74.08 k allocations: 140.091 MiB)
Cross-covariance 1        2.659942 seconds (3.09 M allocations: 10.217 GiB, 16.10% gc time)
Cross-covariance 2        0.003165 seconds (460 allocations: 39.656 KiB)
PSD 1                     0.266184 seconds (119.79 k allocations: 117.523 MiB)
PSD 2                     1.197982 seconds (204.03 k allocations: 2.228 GiB, 15.65% gc time, 0.67% compilation time)
Stationarity: mean        0.380869 seconds (64.88 k allocations: 50.339 MiB)
Stationarity: var         0.374964 seconds (67.18 k allocations: 50.655 MiB)
Stationarity: euclid      1.185207 seconds (4.13 M allocations: 3.850 GiB, 15.02% gc time)
Stationarity: hilbert     0.354309 seconds (211.14 k allocations: 282.190 MiB)
Stationarity: adf         1.339939 seconds (430.67 k allocations: 2.518 GiB)
Mutual information 1      1.408313 seconds (577.05 k allocations: 3.152 GiB)
Mutual information 2      2.623180 seconds (996.75 k allocations: 5.980 GiB, 31.12% gc time)
Entropy                   0.370121 seconds (53.39 k allocations: 50.045 MiB)
Negentropy                0.381115 seconds (57.98 k allocations: 95.057 MiB)
Time coherence            0.003029 seconds (463 allocations: 1.207 MiB)
Signal difference 1       1.411311 seconds (283.58 k allocations: 5.729 GiB, 67.42% gc time, 1.36% compilation time)
Signal difference 2       0.458847 seconds (353.46 k allocations: 6.406 GiB, 9.57% gc time)
Epoch stats               0.284052 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             0.535448 seconds (139.12 k allocations: 4.185 GiB, 5.27% gc time)
Spectrogram 2             0.292904 seconds (180.53 k allocations: 548.287 MiB)
Spectrum 1                0.398053 seconds (122.98 k allocations: 547.012 MiB)
Spectrum 2                1.927965 seconds (219.15 k allocations: 507.536 MiB, 76.68% gc time)
Channel stats             0.375078 seconds (53.41 k allocations: 183.893 MiB)
SNR                       0.372479 seconds (34.99 k allocations: 3.261 MiB)
Standardize               0.046490 seconds (3.76 k allocations: 135.000 MiB)
Frequency convolution     0.461261 seconds (286.08 k allocations: 3.797 GiB)
Time covolution           1.189952 seconds (245.53 k allocations: 1.858 GiB, 49.56% gc time, 2.97% compilation time)
DFT                       0.355312 seconds (113.78 k allocations: 322.109 MiB)
MSCI95                    0.006997 seconds (5.89 k allocations: 38.114 MiB)
Mean                      0.007518 seconds (10.28 k allocations: 47.881 MiB)
Difference                0.554696 seconds (269.19 k allocations: 5.728 GiB)
Temporal envelope 1       2.624528 seconds (36.33 M allocations: 5.903 GiB, 73.24% gc time)
Temporal envelope 2       2.444439 seconds (36.34 M allocations: 5.970 GiB, 83.00% gc time)
Temporal envelope 3       0.925418 seconds (41.45 M allocations: 6.349 GiB, 26.23% gc time)
Power envelope 1          0.422567 seconds (8.05 M allocations: 4.997 GiB)
Power envelope 2          0.337673 seconds (550.28 k allocations: 402.338 MiB, 12.38% gc time)
Power envelope 3          0.369334 seconds (1.05 M allocations: 415.267 MiB, 6.53% gc time)
Spectral envelope 1       0.777644 seconds (444.01 k allocations: 5.575 GiB, 31.44% gc time)
Spectral envelope 2       0.427056 seconds (445.02 k allocations: 5.577 GiB)
Spectral envelope 3       0.472207 seconds (507.80 k allocations: 5.579 GiB)
ISPC 1                    6.941082 seconds (4.42 M allocations: 8.113 GiB, 6.91% gc time)
ISPC 2                    0.002966 seconds (515 allocations: 481.266 KiB)
ITPC                      0.017637 seconds (7.89 k allocations: 19.375 MiB)
ITPC spectrogram        410.317339 seconds (286.04 M allocations: 534.682 GiB, 4.70% gc time)
PLI 1                     5.849859 seconds (3.95 M allocations: 6.786 GiB)
PLI 2                     0.003145 seconds (511 allocations: 420.344 KiB)
AEC                       0.015322 seconds (32.88 k allocations: 96.039 MiB)
GED                       0.015318 seconds (7.85 k allocations: 183.229 MiB)
Instant frequency         0.347715 seconds (263.40 k allocations: 506.299 MiB)
Wavelet spectrogram       8.302736 seconds (3.38 M allocations: 46.275 GiB, 49.19% gc time)
TKEO                      0.382308 seconds (39.59 k allocations: 93.153 MiB)
Wavelet spectrum          3.960182 seconds (3.61 M allocations: 10.707 GiB, 24.70% gc time)
Frequency coherence       0.012239 seconds (585 allocations: 3.439 MiB)
F-test                    0.742024 seconds (104.99 k allocations: 10.169 MiB)
F-test                    0.016787 seconds (58.91 k allocations: 4.770 MiB)
Band power                0.385288 seconds (126.18 k allocations: 94.460 MiB)
Relative PSD              0.264666 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split      1.467082 seconds (19.33 M allocations: 6.860 GiB, 35.10% gc time)
Channel difference        0.015910 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   20.977132 seconds (2.42 M allocations: 40.436 GiB, 0.79% gc time)
Cross power spectrum 1    0.011941 seconds (180 allocations: 1.744 MiB, 48.17% compilation time)
Amplitude difference      0.382161 seconds (117.76 k allocations: 144.862 MiB)

## PLOT

