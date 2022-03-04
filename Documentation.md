
<a id='NeuroJ.jl-Documentation'></a>

<a id='NeuroJ.jl-Documentation-1'></a>

# NeuroJ.jl Documentation


This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).


<a id='EEG-plots'></a>

<a id='EEG-plots-1'></a>

## EEG plots

<a id='NeuroJ.signal_plot-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Vector{Float64}}' href='#NeuroJ.signal_plot-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Vector{Float64}}'>#</a>
**`NeuroJ.signal_plot`** &mdash; *Method*.



```julia
signal_plot(t, signal; <keyword arguments>)
```

Plot single-channel `signal`.

**Arguments**

  * `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
  * `signal::Vector{Float64}`
  * `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, AbstractArray}' href='#NeuroJ.signal_plot-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, AbstractArray}'>#</a>
**`NeuroJ.signal_plot`** &mdash; *Method*.



```julia
signal_plot(t, signal; <keyword arguments>)
```

Plot multi-channel `signal`.

**Arguments**

  * `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
  * `signal::AbstractArray`
  * `labels::Vector{String}=[""]`: labels vector
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot`** &mdash; *Method*.



```julia
eeg_plot(eeg; <keyword arguments>)
```

Plot `eeg` channels. If signal is multi-channel, only channel amplitudes are plotted. For single-channel signal, the histogram, amplitude, power density and spectrogram are plotted.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epochs to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `labels::Vector{String}=[""]`: channel labels vector
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `head::Bool=true`: add head with electrodes
  * `hist::Symbol[:hist, :kd]=:hist`: histogram type
  * `norm::Bool=true`: convert power to dB
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: frequency limit for PSD and spectrogram
  * `kwargs`: other arguments for plot() function; <keyword arguments>

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_avg-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Matrix{Float64}}' href='#NeuroJ.signal_plot_avg-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Matrix{Float64}}'>#</a>
**`NeuroJ.signal_plot_avg`** &mdash; *Method*.



```julia
signal_plot_avg(t, signal; <keyword arguments>)
```

Plot averaged `signal` channels.

**Arguments**

  * `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange`
  * `signal::Matrix{Float64}`
  * `norm::Bool=true`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_avg`** &mdash; *Method*.



```julia
eeg_plot_avg(eeg; <keyword arguments>)
```

Plot averaged `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epoch number to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: frequency limit for PSD and spectrogram
  * `hist::Symbol=:hist`: histogram type: :hist, :kd
  * `head::Bool=true`: add head plot
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_butterfly-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Matrix{Float64}}' href='#NeuroJ.signal_plot_butterfly-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Matrix{Float64}}'>#</a>
**`NeuroJ.signal_plot_butterfly`** &mdash; *Method*.



```julia
signal_plot_butterfly(t, signal; <keyword arguments>)
```

Butterfly plot of `signal` channels.

**Arguments**

  * `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange`
  * `signal::Matrix{Float64}`
  * `labels::Vector{String}=[""]`: channel labels vector
  * `norm::Bool=true`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple`: y-axis limits, default (0, 0)
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_butterfly-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_butterfly-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_butterfly`** &mdash; *Method*.



```julia
eeg_plot_butterfly(eeg; <keyword arguments>)
```

Butterfly plot of `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epoch number to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `labels::Vector{String}=[""]`: channel labels vector
  * `norm::Bool=false`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
  * `head::Bool=true`: add head with electrodes
  * `hist::Bool=true`: add histograms
  * `average::Bool=true`: plot averaged signal
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_psd-Tuple{Vector{Float64}, Vector{Float64}}' href='#NeuroJ.signal_plot_psd-Tuple{Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroJ.signal_plot_psd`** &mdash; *Method*.



```julia
signal_plot_psd(s_powers, s_freqs; <keyword arguments>)
```

Plot power spectrum density.

**Arguments**

  * `s_powers::Vector{Float64}`: signal powers
  * `s_freqs::Vector{Float64}`: signal frequencies
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String="Power [μV^2/Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_psd-Tuple{Vector{Float64}}' href='#NeuroJ.signal_plot_psd-Tuple{Vector{Float64}}'>#</a>
**`NeuroJ.signal_plot_psd`** &mdash; *Method*.



```julia
signal_plot_psd(signal; <keyword arguments>)
```

Plot `signal` channel power spectrum density.

**Arguments**

  * `signal::Vector{Float64}`
  * `fs::Int64`: sampling frequency
  * `norm::Bool=false`: converts power to dB
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String="Power [μV^2/Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_psd-Tuple{Matrix{Float64}}' href='#NeuroJ.signal_plot_psd-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.signal_plot_psd`** &mdash; *Method*.



```julia
signal_plot_psd(signal; <keyword arguments>)
```

Plot `signal` channels power spectrum density.

**Arguments**

  * `signal::Matrix{Float64}`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: power in dB
  * `average::Bool=false`: plots average power and 95%CI for all channels
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: x-axis limit
  * `labels::Vector{String}=[""]`: channel labels vector
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String="Power [μV^2/Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_psd-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_psd-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_psd`** &mdash; *Method*.



```julia
eeg_plot_psd(eeg; <keyword arguments>)
```

Plot `eeg` channels power spectrum density.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epoch number to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `labels::Vector{String}=[""]`: channel labels vector
  * `norm::Bool=false`: power in dB
  * `average::Bool=false`: plots average power and 95%CI for all channels
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String="Power [μV^2/Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `head::Bool=false`: add head with electrodes
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_spectrogram-Tuple{Vector{Float64}}' href='#NeuroJ.signal_plot_spectrogram-Tuple{Vector{Float64}}'>#</a>
**`NeuroJ.signal_plot_spectrogram`** &mdash; *Method*.



```julia
signal_plot_spectrogram(signal; <keyword arguments>)
```

Plot spectrogram of `signal`.

**Arguments**

  * `signal::Vector{Float64}`
  * `fs::Int64`: sampling frequency
  * `offset::Int64=0`: displayed segment offset in samples
  * `norm::Bool=true`: normalize powers to dB
  * `demean::Bool=true`: demean signal prior to analysis
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_spectrogram-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_spectrogram-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_spectrogram`** &mdash; *Method*.



```julia
eeg_plot_spectrogram(eeg; <keyword arguments>)
```

Plots spectrogram of `eeg` channel.

**Arguments**

  * `eeg:EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epoch to plot
  * `channel::Int64`: channel to plot
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize powers to dB
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_histogram-Tuple{Vector{Float64}}' href='#NeuroJ.signal_plot_histogram-Tuple{Vector{Float64}}'>#</a>
**`NeuroJ.signal_plot_histogram`** &mdash; *Method*.



```julia
signal_plot_histogram(signal; <keyword arguments>)
```

Plot histogram of `signal`.

**Arguments**

  * `signal::Vector{Float64}`
  * `type::Symbol`: type of histogram: regular `:hist` or kernel density `:kd`
  * `label::String=""`: channel label
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_histogram-Tuple{VecOrMat{Float64}}' href='#NeuroJ.signal_plot_histogram-Tuple{VecOrMat{Float64}}'>#</a>
**`NeuroJ.signal_plot_histogram`** &mdash; *Method*.



```julia
signal_plot_histogram(signal; <keyword arguments>)
```

Plot histogram of `signal`.

**Arguments**

  * `signal::Matrix{Float64}`
  * `type::Symbol`: type of histogram: :hist or :kd
  * `labels::Vector{String}=[""]`
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_histogram-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_histogram-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_histogram`** &mdash; *Method*.



```julia
eeg_plot_histogram(eeg; <keyword arguments>)
```

Plot `eeg` channel histograms.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `type::Symbol: type of histogram: :hist or :kd
  * `epoch::Int64=1`: epoch number to display
  * `channel::Int64`: channel to display
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default 1 epoch or 20 seconds
  * `label::String=""`: channel label
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_matrix-Tuple{NeuroJ.EEG, Union{Array{Float64, 3}, Matrix{Float64}}}' href='#NeuroJ.eeg_plot_matrix-Tuple{NeuroJ.EEG, Union{Array{Float64, 3}, Matrix{Float64}}}'>#</a>
**`NeuroJ.eeg_plot_matrix`** &mdash; *Method*.



```julia
eeg_plot_matrix(eeg, m; <keyword arguments>)
```

Plot matrix `m` of `eeg` channels.

**Arguments**

  * `eeg:EEG`
  * `m::Union{Matrix{Float64}, Array{Float64, 3}}`: channels by channels matrix
  * `epoch::Int64=1`: epoch number to display
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_covmatrix-Tuple{NeuroJ.EEG, Union{Array{Float64, 3}, Matrix{Float64}}, Union{Vector{Float64}, Vector{Int64}}}' href='#NeuroJ.eeg_plot_covmatrix-Tuple{NeuroJ.EEG, Union{Array{Float64, 3}, Matrix{Float64}}, Union{Vector{Float64}, Vector{Int64}}}'>#</a>
**`NeuroJ.eeg_plot_covmatrix`** &mdash; *Method*.



```julia
eeg_plot_matrix(eeg, cov_m, lags; <keyword arguments>)
```

Plot covariance matrix `m` of `eeg` channels.

**Arguments**

  * `eeg:EEG`
  * `cov_m::Union{Matrix{Float64}, Array{Float64, 3}}`: covariance matrix
  * `lags::Union{Vector{Int64}, Vector{Float64}}`: covariance lags
  * `channel::Union{Int64, Vector{Int64}, UnitRange{Int64}, Nothing}`: channel to display
  * `epoch::Int64=1`: epoch number to display
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_electrodes-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_electrodes-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_electrodes`** &mdash; *Method*.



```julia
eeg_plot_electrodes(eeg; <keyword arguments>)
```

Plot `eeg` electrodes.

**Arguments**

  * `eeg:EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: which channel should be highlighted, default is all channels
  * `labels::Bool=true`: plot electrode labels
  * `head::Bool`=true: plot head
  * `head_labels::Bool=false`: plot head labels
  * `small::Bool=false`: draws small plot
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_draw_head-Tuple{Plots.Plot{Plots.GRBackend}, Vector{Float64}, Vector{Float64}}' href='#NeuroJ.eeg_draw_head-Tuple{Plots.Plot{Plots.GRBackend}, Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroJ.eeg_draw_head`** &mdash; *Method*.



```julia
eeg_draw_head(p, loc_x, loc_y; head_labels, kwargs)
```

Draw head over a topographical plot `p`.

**Arguments**

  * `p::Plots.Plot{Plots.GRBackend}`: electrodes plot
  * `loc_x::Vector{Float64}`: vector of x electrode position
  * `loc_y::Vector{Float64}`: vector of y electrode position
  * `head_labels::Bool=true`: add text labels to the plot
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_filter_response-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_filter_response-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_filter_response`** &mdash; *Method*.



```julia
eeg_plot_filter_response(eeg; <keyword arguments>)
```

Plot filter response.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `fprototype::Symbol`: filter class: :butterworth, :chebyshev1, :chebyshev2, :elliptic
  * `ftype::Symbol`: filter type: :lp, :hp, :bp, :bs
  * `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64`: filter order
  * `rp::Union{Int64, Float64}`: dB ripple in the passband
  * `rs::Union{Int64, Float64}`: dB attenuation in the stopband
  * `window::window::Union{Vector{Float64}, Nothing}`: window, required for FIR filter
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

