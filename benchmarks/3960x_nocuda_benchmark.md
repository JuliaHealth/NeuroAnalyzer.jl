# System info

OS: Linux shogun 5.19.0-1-amd64 #1 SMP PREEMPT_DYNAMIC Debian 5.19.6-1 (2022-09-01) x86_64 GNU/Linux
CPU: AMD Ryzen Threadripper 3960X 24-Core Processor 24 cores (3.95 GHz)
RAM: 112.6 GB free / 125.7 GB

# Loading NeuroAnalyzer.jl

257.589180 seconds (154.54 M allocations: 9.258 GiB, 1.21% gc time, 2.81% compilation time: 73% of which was recompilation)

# NeuroAnalyzer info

NeuroAnalyzer: 0.22.9
        Julia: 1.8.1
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

Import EDF+               0.376361 seconds (3.12 M allocations: 501.761 MiB, 17.27% gc time)
Import BDF+               0.071652 seconds (20.57 k allocations: 108.687 MiB, 42.96% gc time)
Import EDF                0.098781 seconds (42.02 k allocations: 152.768 MiB)

# PROCESS

A referencing             0.843407 seconds (43.05 k allocations: 350.744 MiB, 49.65% gc time)
M referencing             0.456846 seconds (52.37 k allocations: 351.202 MiB, 5.37% gc time)
CAR referencing           0.371279 seconds (121.09 k allocations: 895.928 MiB)
Channel referencing       0.385689 seconds (44.58 k allocations: 48.760 MiB)
Filter notch              0.376770 seconds (88.43 k allocations: 590.255 MiB)
Filter LP                 0.383216 seconds (1.06 M allocations: 648.003 MiB)
Filter HP                 0.368220 seconds (1.06 M allocations: 648.283 MiB)

# ANALYZE

Total power               0.311118 seconds (124.25 k allocations: 77.939 MiB)
Band power                0.299298 seconds (159.01 k allocations: 89.827 MiB)
Covariance matrix         0.005030 seconds (2.99 k allocations: 46.114 MiB)
Correlation matrix        0.004855 seconds (3.11 k allocations: 46.138 MiB)
Auto-covariance           0.366878 seconds (74.08 k allocations: 140.091 MiB)
Cross-covariance 1        3.188031 seconds (3.10 M allocations: 10.217 GiB, 8.52% gc time)
Cross-covariance 2        0.002991 seconds (460 allocations: 39.578 KiB)
PSD 1                     0.306439 seconds (118.80 k allocations: 117.493 MiB)
PSD 2                     0.898634 seconds (179.53 k allocations: 2.227 GiB)
Stationarity: mean        0.373936 seconds (64.88 k allocations: 50.339 MiB)
Stationarity: var         0.373473 seconds (67.18 k allocations: 50.655 MiB)
Stationarity: euclid      1.104919 seconds (4.13 M allocations: 3.850 GiB, 16.48% gc time)
Stationarity: hilbert     0.363830 seconds (211.69 k allocations: 282.206 MiB)
Stationarity: adf         1.074087 seconds (430.84 k allocations: 2.518 GiB)
Mutual information 1      1.236480 seconds (576.80 k allocations: 3.152 GiB)
Mutual information 2      2.036735 seconds (996.76 k allocations: 5.980 GiB)
Entropy                   0.367610 seconds (53.39 k allocations: 50.045 MiB)
Negentropy                0.363420 seconds (57.98 k allocations: 95.057 MiB)
Time coherence            0.003416 seconds (463 allocations: 1.207 MiB)
Signal difference 1       0.867096 seconds (269.27 k allocations: 5.728 GiB)
Signal difference 2      12.715415 seconds (353.37 k allocations: 6.406 GiB, 96.35% gc time)
Epoch stats               0.299070 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             0.576527 seconds (139.76 k allocations: 4.185 GiB, 6.27% gc time)
Spectrogram 2             0.310033 seconds (180.90 k allocations: 548.298 MiB)
Spectrum 1                0.393986 seconds (120.65 k allocations: 546.871 MiB)
Spectrum 2                0.398755 seconds (221.47 k allocations: 507.606 MiB)
Channel stats             0.361307 seconds (53.41 k allocations: 183.893 MiB)
SNR                       0.367982 seconds (34.99 k allocations: 3.261 MiB)
Standardize               0.079581 seconds (3.76 k allocations: 134.981 MiB, 55.04% gc time)
Frequency convolution     0.450615 seconds (289.45 k allocations: 3.797 GiB)
Time covolution           0.427285 seconds (216.37 k allocations: 1.857 GiB)
DFT                       0.375894 seconds (111.44 k allocations: 321.968 MiB)
MSCI95                    0.005764 seconds (5.89 k allocations: 38.113 MiB)
Mean                      0.006792 seconds (10.28 k allocations: 47.881 MiB)
Difference                0.531477 seconds (269.13 k allocations: 5.728 GiB, 14.00% gc time)
Temporal envelope 1       0.435527 seconds (36.33 M allocations: 5.705 GiB)
Temporal envelope 2       1.380238 seconds (36.34 M allocations: 5.773 GiB, 70.02% gc time)
Temporal envelope 3       0.849163 seconds (41.46 M allocations: 6.112 GiB, 20.73% gc time)
Power envelope 1          0.455275 seconds (8.05 M allocations: 4.997 GiB)
Power envelope 2          0.365147 seconds (545.89 k allocations: 402.204 MiB)
Power envelope 3          0.482299 seconds (1.05 M allocations: 415.158 MiB, 7.25% gc time)
Spectral envelope 1       0.594044 seconds (436.99 k allocations: 5.575 GiB, 21.97% gc time)
Spectral envelope 2       0.415190 seconds (444.43 k allocations: 5.577 GiB)
Spectral envelope 3       0.434286 seconds (507.96 k allocations: 5.579 GiB)
ISPC 1                    6.869175 seconds (4.41 M allocations: 8.112 GiB, 6.75% gc time)
ISPC 2                    0.003090 seconds (515 allocations: 481.266 KiB)
ITPC                      0.017721 seconds (7.89 k allocations: 19.375 MiB)
ITPC spectrogram        412.699010 seconds (285.78 M allocations: 534.674 GiB, 4.64% gc time)
PLI 1                     5.920231 seconds (3.94 M allocations: 6.785 GiB)
PLI 2                     0.002987 seconds (511 allocations: 420.344 KiB)
AEC                       0.022202 seconds (32.88 k allocations: 95.825 MiB)
GED                       0.016494 seconds (7.85 k allocations: 183.229 MiB)
Instant frequency         0.363098 seconds (268.00 k allocations: 506.439 MiB)
Wavelet spectrogram       3.913532 seconds (3.40 M allocations: 40.486 GiB, 5.73% gc time)
TKEO                      0.376796 seconds (39.59 k allocations: 93.153 MiB)
Wavelet spectrum          3.277540 seconds (3.55 M allocations: 10.705 GiB, 11.85% gc time)
Frequency coherence       0.011956 seconds (585 allocations: 3.439 MiB)
F-test                    0.719889 seconds (104.99 k allocations: 10.169 MiB)
F-test                    0.015536 seconds (58.91 k allocations: 4.770 MiB)
Band power                0.354766 seconds (126.63 k allocations: 94.473 MiB)
Relative PSD              0.240613 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split      1.447256 seconds (19.33 M allocations: 6.860 GiB, 34.97% gc time)
Channel difference        0.013176 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   22.593389 seconds (2.48 M allocations: 40.437 GiB, 6.86% gc time)
Cross power spectrum 2    0.011719 seconds (180 allocations: 1.744 MiB, 49.03% compilation time)
Amplitude difference      0.373758 seconds (117.76 k allocations: 144.862 MiB)

# PLOT

