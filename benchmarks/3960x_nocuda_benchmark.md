# System info

OS: Linux shogun 5.18.0-4-amd64 #1 SMP PREEMPT_DYNAMIC Debian 5.18.16-1 (2022-08-10) x86_64 GNU/Linux
CPU: AMD Ryzen Threadripper 3960X 24-Core Processor 24 cores (3.71 GHz)
RAM: 79.0 GB free / 125.7 GB

# Loading NeuroAnalyzer.jl

238.367735 seconds (156.70 M allocations: 9.381 GiB, 1.42% gc time, 3.06% compilation time: 74% of which was recompilation)

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

Import EDF+               0.700253 seconds (3.12 M allocations: 501.761 MiB, 60.08% gc time)
Import BDF+               0.045826 seconds (20.57 k allocations: 108.687 MiB)
Import EDF                0.099344 seconds (42.02 k allocations: 152.768 MiB)

# PROCESS

A referencing             0.474713 seconds (43.05 k allocations: 350.744 MiB, 8.08% gc time)
M referencing             0.457425 seconds (52.37 k allocations: 351.202 MiB, 4.90% gc time)
CAR referencing           0.379882 seconds (121.19 k allocations: 900.054 MiB)
Channel referencing       0.390979 seconds (44.58 k allocations: 48.760 MiB)
Filter notch              0.396645 seconds (88.43 k allocations: 590.255 MiB)
Filter LP                 0.401306 seconds (1.06 M allocations: 648.003 MiB)
Filter HP                 0.396505 seconds (1.06 M allocations: 648.284 MiB)

# ANALYZE

Total power               0.347172 seconds (122.47 k allocations: 77.885 MiB)
Band power                0.336086 seconds (157.43 k allocations: 89.778 MiB)
Covariance matrix         0.005451 seconds (2.99 k allocations: 46.114 MiB)
Correlation matrix        0.005146 seconds (3.11 k allocations: 46.138 MiB)
Auto-covariance           0.375944 seconds (74.08 k allocations: 140.091 MiB)
Cross-covariance 1        2.087867 seconds (3.09 M allocations: 10.217 GiB, 9.97% gc time)
Cross-covariance 2        0.003135 seconds (460 allocations: 39.578 KiB)
PSD 1                     0.269087 seconds (119.75 k allocations: 117.522 MiB)
PSD 2                     0.931599 seconds (182.79 k allocations: 2.228 GiB)
Stationarity: mean        0.385498 seconds (64.88 k allocations: 50.339 MiB)
Stationarity: var         0.380308 seconds (67.18 k allocations: 50.655 MiB)
Stationarity: euclid      1.196470 seconds (4.13 M allocations: 3.850 GiB, 15.00% gc time)
Stationarity: hilbert     0.366667 seconds (212.77 k allocations: 282.239 MiB)
Stationarity: adf         1.242435 seconds (430.65 k allocations: 2.518 GiB)
Mutual information 1      1.382742 seconds (577.03 k allocations: 3.152 GiB)
Mutual information 2      5.013229 seconds (1.01 M allocations: 5.981 GiB, 75.75% gc time, 0.42% compilation time)
Entropy                   0.379258 seconds (53.39 k allocations: 50.045 MiB)
Negentropy                0.372340 seconds (57.98 k allocations: 95.057 MiB)
Time coherence            0.003147 seconds (463 allocations: 1.207 MiB)
Signal difference 1       0.538418 seconds (268.95 k allocations: 5.728 GiB)
Signal difference 2       0.477394 seconds (353.37 k allocations: 6.406 GiB)
Epoch stats               0.270850 seconds (2.94 k allocations: 179.983 MiB)
Spectrogram 1             0.625156 seconds (143.33 k allocations: 4.185 GiB, 11.04% gc time)
Spectrogram 2             0.297530 seconds (181.17 k allocations: 548.307 MiB)
Spectrum 1                0.390226 seconds (120.63 k allocations: 546.871 MiB)
Spectrum 2                0.405891 seconds (221.30 k allocations: 507.601 MiB)
Channel stats             0.654827 seconds (53.41 k allocations: 183.893 MiB, 27.14% gc time)
SNR                       0.368331 seconds (34.99 k allocations: 3.261 MiB)
Standardize               0.033194 seconds (3.76 k allocations: 134.981 MiB)
Frequency convolution     0.453365 seconds (304.02 k allocations: 3.798 GiB)
Time covolution           0.486932 seconds (217.96 k allocations: 1.857 GiB)
DFT                       0.386935 seconds (111.41 k allocations: 321.967 MiB)
MSCI95                    0.004481 seconds (5.89 k allocations: 38.114 MiB)
Mean                      0.005559 seconds (10.28 k allocations: 47.881 MiB)
Difference                0.490660 seconds (269.06 k allocations: 5.728 GiB)
Temporal envelope 1       0.412724 seconds (36.33 M allocations: 5.705 GiB)
Temporal envelope 2       0.971417 seconds (36.34 M allocations: 5.773 GiB, 58.59% gc time)
Temporal envelope 3       0.837607 seconds (41.46 M allocations: 6.112 GiB, 20.33% gc time)
Power envelope 1          0.414785 seconds (8.05 M allocations: 4.997 GiB)
Power envelope 2          0.404974 seconds (547.03 k allocations: 402.239 MiB, 11.55% gc time)
Power envelope 3          0.376400 seconds (1.05 M allocations: 415.195 MiB)
Spectral envelope 1       0.616452 seconds (439.14 k allocations: 5.575 GiB, 22.58% gc time)
Spectral envelope 2       0.619064 seconds (458.92 k allocations: 5.577 GiB, 20.40% gc time, 1.21% compilation time)
Spectral envelope 3       0.659322 seconds (521.74 k allocations: 5.579 GiB, 24.83% gc time)
ISPC 1                    6.679597 seconds (4.33 M allocations: 8.110 GiB, 6.08% gc time, 0.06% compilation time)
ISPC 2                    0.003100 seconds (515 allocations: 481.266 KiB)
ITPC                      0.018076 seconds (7.89 k allocations: 19.375 MiB)
ITPC spectrogram        411.502770 seconds (283.43 M allocations: 534.609 GiB, 4.23% gc time)
PLI 1                     6.537737 seconds (4.26 M allocations: 6.792 GiB, 4.82% gc time)
PLI 2                     0.003140 seconds (511 allocations: 420.344 KiB)
AEC                       0.014798 seconds (32.88 k allocations: 95.825 MiB)
GED                       0.013897 seconds (7.85 k allocations: 183.229 MiB)
Instant frequency         1.055440 seconds (526.21 k allocations: 511.744 MiB, 37.74% gc time)
Wavelet spectrogram       4.358414 seconds (3.48 M allocations: 40.488 GiB, 5.99% gc time)
TKEO                      0.380039 seconds (39.59 k allocations: 93.153 MiB)
Wavelet spectrum          8.698889 seconds (3.73 M allocations: 10.709 GiB, 62.27% gc time)
Frequency coherence       0.012694 seconds (585 allocations: 3.439 MiB)
F-test                    0.759772 seconds (104.99 k allocations: 10.169 MiB)
F-test                    0.017456 seconds (58.90 k allocations: 4.770 MiB)
Band power                0.367762 seconds (126.08 k allocations: 94.457 MiB)
Relative PSD              0.275298 seconds (211.63 k allocations: 155.908 MiB)
Frequency band split      2.670813 seconds (19.33 M allocations: 6.860 GiB, 60.76% gc time)
Channel difference        0.016033 seconds (2.27 k allocations: 11.894 MiB)
Cross power spectrum 1   22.200708 seconds (2.48 M allocations: 40.437 GiB, 4.62% gc time)
Cross power spectrum 1    0.011755 seconds (180 allocations: 1.744 MiB, 49.59% compilation time)
Amplitude difference      0.380257 seconds (117.76 k allocations: 144.862 MiB)

# PLOT

