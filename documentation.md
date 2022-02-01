# NeuroJ.jl

Welcome to the NeuroJ.jl documentation!

## Installation

```
using Pkg
Pkg.add(url="https://notabug.org/AdamWysokinski/NeuroJ.jl")
```

## General remarks

Single-channel signals are column vectors.

Multi-channel or multi-trial signals are matrices of channels/trials by signals.

Function name prefix:
- signal_  :: functions taking signal or signals as argument
- eeg_     :: functions taking EEG object as argument

EEG object (headers + data) is stored in the EEG structure:
```
struct EEG
    eeg_object_header::Dict
    eeg_signal_header::Dict
    eeg_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
end
```