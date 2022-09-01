# NeuroAnalyzer benchmarks

CPU: AMD Ryzen 5 PRO 4650U with Radeon Graphics 6 cores (3.9 GHz)
RAM: 10.5 GB free /  22.7 GB

NeuroAnalyzer: 0.22.8
        Julia: 1.8.0
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

## IO

Import EDF+               0.445492 seconds (3.12 M allocations: 498.707 MiB, 17.95% gc time)
Import BDF+               0.088119 seconds (20.57 k allocations: 107.749 MiB, 45.79% gc time)
Import EDF                0.118767 seconds (42.01 k allocations: 150.405 MiB)

## PROCESS

CAR referencing           0.374914 seconds (96.22 k allocations: 934.004 MiB)
Channel referencing       0.040107 seconds (18.65 k allocations: 46.490 MiB)
Filter notch              0.133047 seconds (62.49 k allocations: 587.653 MiB)
Filter LP                 0.748905 seconds (1.04 M allocations: 645.401 MiB, 78.63% gc time)
Filter HP                 0.406025 seconds (1.04 M allocations: 645.683 MiB)

## ANALYZE

Total power               0.252331 seconds (106.86 k allocations: 75.712 MiB)
Band power                0.303711 seconds (140.62 k allocations: 87.502 MiB)
Covariance matrix         0.017796 seconds (2.77 k allocations: 46.094 MiB)
Correlation matrix        0.019414 seconds (2.89 k allocations: 46.118 MiB)
Auto-covariance           0.184151 seconds (48.20 k allocations: 137.605 MiB)
Cross-covariance 1        1.850023 seconds (3.02 M allocations: 10.210 GiB, 10.09% gc time)
Cross-covariance 2        0.000134 seconds (244 allocations: 16.859 KiB)
PSD 1                     0.260424 seconds (103.19 k allocations: 115.320 MiB)
PSD 2                     2.063021 seconds (156.19 k allocations: 2.225 GiB, 3.06% gc time)
Stationarity: mean        0.020875 seconds (38.91 k allocations: 47.850 MiB)
Stationarity: var         0.021289 seconds (41.21 k allocations: 48.165 MiB)
Stationarity: euclid      1.213733 seconds (4.08 M allocations: 3.845 GiB, 31.02% gc time)
Stationarity: hilbert     0.395343 seconds (178.66 k allocations: 279.568 MiB)
Stationarity: adf         2.187837 seconds (404.51 k allocations: 2.515 GiB)
Mutual information 1      1.297275 seconds (523.77 k allocations: 3.147 GiB, 8.33% gc time)
Mutual information 2      2.103725 seconds (970.16 k allocations: 5.978 GiB, 8.07% gc time)
Entropy                   0.076283 seconds (27.45 k allocations: 47.623 MiB)
Negentropy                0.081226 seconds (32.06 k allocations: 92.636 MiB)
Time coherence            0.000993 seconds (244 allocations: 1.186 MiB)
Signal difference 1       1.016393 seconds (242.78 k allocations: 5.726 GiB, 7.09% gc time)
Signal difference 2       1.068266 seconds (326.99 k allocations: 6.404 GiB, 4.48% gc time)
Epoch stats               0.356498 seconds (2.94 k allocations: 179.983 MiB, 8.96% gc time)
Spectrogram 1             1.300768 seconds (113.89 k allocations: 4.182 GiB, 2.69% gc time)
Spectrogram 2             0.428097 seconds (153.27 k allocations: 545.758 MiB)
Spectrum 1                1.572107 seconds (128.77 k allocations: 545.482 MiB, 73.33% gc time, 0.27% compilation time)
Spectrum 2                0.486036 seconds (187.78 k allocations: 504.749 MiB)
Channel stats             0.063892 seconds (27.47 k allocations: 181.139 MiB)
SNR                       0.044087 seconds (9.04 k allocations: 858.477 KiB)
Standardize               0.629874 seconds (3.76 k allocations: 135.000 MiB, 89.18% gc time)
Frequency convolution     1.983204 seconds (258.37 k allocations: 3.795 GiB, 29.73% gc time, 1.16% compilation time)
Time covolution           0.757865 seconds (186.26 k allocations: 1.854 GiB)
DFT                       0.154282 seconds (94.48 k allocations: 319.969 MiB)
MSCI95                    0.014999 seconds (5.67 k allocations: 38.092 MiB)
Mean                      0.019080 seconds (10.05 k allocations: 47.859 MiB)
Difference                0.979963 seconds (242.74 k allocations: 5.726 GiB, 5.81% gc time)
Temporal envelope 1       1.213916 seconds (36.30 M allocations: 5.900 GiB, 16.69% gc time)
Temporal envelope 2       1.452973 seconds (36.31 M allocations: 5.967 GiB, 28.01% gc time)
Temporal envelope 3       1.763259 seconds (41.43 M allocations: 6.347 GiB, 23.78% gc time)
Power envelope 1          0.968455 seconds (8.03 M allocations: 4.994 GiB, 8.32% gc time)
Power envelope 2          0.358842 seconds (530.09 k allocations: 400.025 MiB, 8.16% gc time)
Power envelope 3          0.413006 seconds (1.03 M allocations: 412.995 MiB, 5.40% gc time)
Spectral envelope 1       1.097702 seconds (414.17 k allocations: 5.572 GiB, 5.92% gc time)
Spectral envelope 2       1.126373 seconds (422.85 k allocations: 5.574 GiB, 4.87% gc time)
Spectral envelope 3       1.158876 seconds (485.58 k allocations: 5.576 GiB, 4.89% gc time)
ISPC 1                    6.047754 seconds (3.85 M allocations: 8.093 GiB, 2.63% gc time)
ISPC 2                    0.001413 seconds (300 allocations: 459.062 KiB)
ITPC                      0.019973 seconds (7.89 k allocations: 19.375 MiB)
ITPC spectrogram        413.099822 seconds (286.50 M allocations: 534.696 GiB, 2.43% gc time)
PLI 1                     5.697892 seconds (3.85 M allocations: 6.778 GiB, 1.45% gc time)
PLI 2                     0.000623 seconds (296 allocations: 398.703 KiB)
AEC                       0.036550 seconds (32.45 k allocations: 95.998 MiB)
GED                       0.057057 seconds (7.63 k allocations: 183.206 MiB)
Instant frequency         0.483032 seconds (214.97 k allocations: 504.661 MiB, 21.65% gc time)
Wavelet spectrogram      16.329059 seconds (3.09 M allocations: 46.266 GiB, 17.56% gc time)
TKEO                      0.038149 seconds (13.63 k allocations: 90.730 MiB)
Wavelet spectrum          3.513731 seconds (3.26 M allocations: 10.696 GiB, 4.27% gc time)
Frequency coherence       0.013228 seconds (369 allocations: 3.418 MiB)
F-test                    0.068460 seconds (52.96 k allocations: 5.387 MiB)
F-test                    0.113173 seconds (58.68 k allocations: 4.747 MiB)
Band power                0.226510 seconds (112.40 k allocations: 92.209 MiB)
Relative PSD              0.273646 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split      1.798697 seconds (19.30 M allocations: 6.858 GiB, 6.84% gc time)
Channel difference        0.011360 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   27.111315 seconds (2.00 M allocations: 40.397 GiB, 3.38% gc time)
Cross power spectrum 1    0.015419 seconds (180 allocations: 1.744 MiB, 58.30% compilation time)
Amplitude difference      0.074223 seconds (91.78 k allocations: 142.372 MiB)

## PLOT

