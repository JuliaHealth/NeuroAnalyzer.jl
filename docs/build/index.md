
<a id='NeuroJ.jl-Documentation'></a>

<a id='NeuroJ.jl-Documentation-1'></a>

# NeuroJ.jl Documentation


This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).


<a id='NeuroJ'></a>

<a id='NeuroJ-1'></a>

## NeuroJ

<a id='NeuroJ.neuroj_version-Tuple{}' href='#NeuroJ.neuroj_version-Tuple{}'>#</a>
**`NeuroJ.neuroj_version`** &mdash; *Method*.



```julia
neuroj_version()
```

Shows NeuroJ and imported packages versions.

<a id='NeuroJ.neuroj_reload_plugins-Tuple{}' href='#NeuroJ.neuroj_reload_plugins-Tuple{}'>#</a>
**`NeuroJ.neuroj_reload_plugins`** &mdash; *Method*.



```julia
neuroj_reload_plugins()
```

Reload NeuroJ plugins. Plugins path is: ~/Documents/NeuroJ/plugins/


<a id='EEG-io'></a>

<a id='EEG-io-1'></a>

## EEG io

<a id='NeuroJ.eeg_import_edf-Tuple{String}' href='#NeuroJ.eeg_import_edf-Tuple{String}'>#</a>
**`NeuroJ.eeg_import_edf`** &mdash; *Method*.



```julia
eeg_edit(file_name; read_annotations=true, header_only=false, clean_labels=true)
```

Loads EDF/EDFPlus file and returns EEG object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `read_annotations::Bool`: read annotations from EDF+ file (currently not implemented)
  * `clean_labels::Bool`: only keep channel names in channel labels

**Returns**

  * `eeg:EEG`

**Notes**

sampling*rate = n.samples / data.record.duration gain = (physical*maximum - physical*minimum) / (digital*maximum - digital*minimum) value = (value - digital*minimum ) * gain + physical_minimum

**Source**

Kemp B, Värri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3. 

<a id='NeuroJ.eeg_import_ced-Tuple{String}' href='#NeuroJ.eeg_import_ced-Tuple{String}'>#</a>
**`NeuroJ.eeg_import_ced`** &mdash; *Method*.



```julia
eeg_import_ced(file_name)
```

Loads electrode positions from CED file.

**Arguments**

  * `file_name::String`

**Returns**

  * `sensors::DataFrame`

<a id='NeuroJ.eeg_import_locs-Tuple{String}' href='#NeuroJ.eeg_import_locs-Tuple{String}'>#</a>
**`NeuroJ.eeg_import_locs`** &mdash; *Method*.



```julia
eeg_import_locs(file_name)
```

Loads electrode positions from LOCS file.

**Arguments**

  * `file_name::String`

**Returns**

  * `sensors::DataFrame`

<a id='NeuroJ.eeg_import_elc-Tuple{String}' href='#NeuroJ.eeg_import_elc-Tuple{String}'>#</a>
**`NeuroJ.eeg_import_elc`** &mdash; *Method*.



```julia
eeg_import_elc(file_name)
```

Loads electrode positions from ELC file.

**Arguments**

  * `file_name::String`

**Returns**

  * `sensors::DataFrame`

<a id='NeuroJ.eeg_load_electrodes-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_load_electrodes-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_load_electrodes`** &mdash; *Method*.



```julia
eeg_load_electrodes(eeg; file_name)
```

Loads electrode positions from `file_name`. Accepted formats:

  * CED
  * LOCS
  * ELC

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `file_name::String`

**Returns**

  * `eeg:EEG`

<a id='NeuroJ.eeg_load_electrodes!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_load_electrodes!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_load_electrodes!`** &mdash; *Method*.



```julia
eeg_load_electrodes!(eeg; file_name)
```

Loads electrode positions from `file_name`. Accepted formats:

  * CED
  * LOCS
  * ELC

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `file_name::String`

<a id='NeuroJ.eeg_load-Tuple{String}' href='#NeuroJ.eeg_load-Tuple{String}'>#</a>
**`NeuroJ.eeg_load`** &mdash; *Method*.



```julia
eeg_load(file_name)
```

Loads the `eeg` from `file_name` file (HDF5-based).

**Arguments**

  * `file_name::String`: file name

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_save-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_save-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_save`** &mdash; *Method*.



```julia
eeg_save(eeg; file_name, overwrite=false)
```

Saves the `eeg` to `file_name` file (HDF5-based).

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `file_name::String`: file name
  * `overwrite::Bool=false`

**Returns**

  * `success::Bool`

<a id='NeuroJ.eeg_export_csv-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_export_csv-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_export_csv`** &mdash; *Method*.



```julia
eeg_export_csv(eeg, file_name, header, overwrite)
```

Exports EEG data as CSV.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `file_name::String`
  * `header::Bool=false`: export header
  * `components::Bool=false`: export components
  * `overwrite::Bool=false`

**Returns**

  * `success::Bool`


<a id='EEG-edit'></a>

<a id='EEG-edit-1'></a>

## EEG edit


<a id='EEG-process'></a>

<a id='EEG-process-1'></a>

## EEG process

<a id='NeuroJ.eeg_reference_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_channel`** &mdash; *Method*.



```julia
eeg_reference_channel(eeg; channel)
```

Reference the `eeg` to specific channel `channel`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_reference_channel!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_channel!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_channel!`** &mdash; *Method*.



```julia
eeg_reference_channel!(eeg; channel)
```

Reference the `eeg` to specific channel `channel`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference

<a id='NeuroJ.eeg_reference_car-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_car-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_car`** &mdash; *Method*.



```julia
eeg_reference_car(eeg)
```

Reference the `eeg` to common average reference.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_reference_car!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_car!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_car!`** &mdash; *Method*.



```julia
eeg_reference_car!(eeg)
```

Reference the `eeg` to common average reference.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_derivative-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_derivative-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_derivative`** &mdash; *Method*.



```julia
eeg_derivative(eeg)
```

Return the derivative of the `eeg` with length same as the signal.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_derivative!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_derivative!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_derivative!`** &mdash; *Method*.



```julia
eeg_derivative!(eeg)
```

Return the derivative of the `eeg` with length same as the signal.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_detrend-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_detrend-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_detrend`** &mdash; *Method*.



```julia
eeg_detrend(eeg; type)
```

Remove linear trend from the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Symbol=:linear`, optional

      * `:linear`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `:constant`: the mean of `signal` is subtracted

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_detrend!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_detrend!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_detrend!`** &mdash; *Method*.



```julia
eeg_detrend!(eeg; type)
```

Remove linear trend from the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Symbol=:linear`, optional

      * `:linear`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `:constant`: the mean of `signal` is subtracted

<a id='NeuroJ.eeg_taper-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_taper-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_taper`** &mdash; *Method*.



```julia
eeg_taper(eeg; taper)
```

Taper `eeg` with `taper`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `taper::Vector`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_taper!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_taper!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_taper!`** &mdash; *Method*.



```julia
eeg_taper!(eeg; taper)
```

Taper `eeg` with `taper`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `taper::Vector`

<a id='NeuroJ.eeg_demean-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_demean-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_demean`** &mdash; *Method*.



```julia
eeg_demean(eeg)
```

Remove mean value (DC offset).

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_demean!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_demean!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_demean!`** &mdash; *Method*.



```julia
eeg_demean!(eeg)
```

Remove mean value (DC offset).

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_normalize_zscore-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_normalize_zscore-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_normalize_zscore`** &mdash; *Method*.



```julia
eeg_normalize_zscore(eeg)
```

Normalize by z-score.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_normalize_zscore!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_normalize_zscore!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_normalize_zscore!`** &mdash; *Method*.



```julia
eeg_normalize_zscore!(eeg)
```

Normalize by z-score.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_normalize_minmax-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_normalize_minmax-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_normalize_minmax`** &mdash; *Method*.



```julia
eeg_normalize_minmax(eeg)
```

Normalize to 0...1

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_normalize_minmax!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_normalize_minmax!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_normalize_minmax!`** &mdash; *Method*.



```julia
eeg_normalize_minmax!(eeg)
```

Normalize to 0...1

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_average-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_average-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_average`** &mdash; *Method*.



```julia
eeg_average(eeg)
```

Returns the average signal of all `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`


!!! warning "Missing docstring."
    Missing docstring for `eeg_average!(eeg)`. Check Documenter's build log for details.


<a id='NeuroJ.eeg_resample-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_resample-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_resample`** &mdash; *Method*.



```julia
eeg_resample(eeg; new_sr)
```

Resample all channels of `eeg` to `new_sr` sampling frequency.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_resample!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_resample!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_resample!`** &mdash; *Method*.



```julia
eeg_resample!(eeg; new_sr)
```

Resample all channels of `eeg` to `new_sr` sampling frequency.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `new_sr::Int64`: new sampling rate

<a id='NeuroJ.eeg_upsample-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_upsample-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_upsample`** &mdash; *Method*.



```julia
eeg_upsample(eeg; new_sr)
```

Upsample all channels of `eeg` to `new_sr` sampling frequency.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_upsample!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_upsample!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_upsample!`** &mdash; *Method*.



```julia
eeg_upsample!(eeg; new_sr)
```

Upsample all channels of `eeg` to `new_sr` sampling frequency.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_downsample-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_downsample-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_downsample`** &mdash; *Method*.



```julia
eeg_downsample(eeg; new_sr)
```

Downsamples all channels of `eeg` to `new_sr` sampling frequency.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_downsample!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_downsample!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_downsample!`** &mdash; *Method*.



```julia
eeg_downsample!(eeg; new_sr)
```

Downsamples all channels of `eeg` to `new_sr` sampling frequency.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `new_sr::Int64`: new sampling rate

<a id='NeuroJ.eeg_filter-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_filter-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_filter`** &mdash; *Method*.



```julia
eeg_filter(eeg; <keyword arguments>)
```

Filter `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `fprototype::Symbol`: filter prototype:

      * `:butterworth`
      * `:chebyshev1`
      * `:chebyshev2`
      * `:elliptic`
      * `:fir`
      * `:mavg`: moving average (with threshold and/or weight window)
      * `:mmed`: moving median (with threshold and/or weight window)
      * `:poly`: polynomial of `order` order
  * `ftype::Symbol`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Int64, Float64, Tuple}=0`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64`=0: filter order
  * `rp::Union{Int64, Float64}=-1`: dB ripple in the passband
  * `rs::Union{Int64, Float64}=-1`: dB attentuation in the stopband
  * `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass)
  * `d::Int64=1`: window length for mean average and median average filter
  * `t::Union{Int64, Float64}=0`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
  * `window::Union{Vector{Float64}, Nothing}=nothing`: window, required for FIR filter

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_filter!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_filter!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_filter!`** &mdash; *Method*.



```julia
eeg_filter!(eeg; <keyword arguments>)
```

Filter `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `fprototype::Symbol`: filter prototype:

      * `:butterworth`
      * `:chebyshev1`
      * `:chebyshev2`
      * `:elliptic`
      * `:fir`
      * `:mavg`: moving average (with threshold and/or weight window)
      * `:mmed`: moving median (with threshold and/or weight window)
      * `:poly`: polynomial of `order` order
  * `ftype::Symbol`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Int64, Float64, Tuple}=0`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64`=0: filter order
  * `rp::Union{Int64, Float64}=-1`: dB ripple in the passband
  * `rs::Union{Int64, Float64}=-1`: dB attentuation in the stopband
  * `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass)
  * `d::Int64=1`: window length for mean average and median average filter
  * `t::Union{Int64, Float64}=0`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
  * `window::Union{Vector{Float64}, Nothing}=nothing`: window, required for FIR filter

<a id='NeuroJ.eeg_tconv-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_tconv-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_tconv`** &mdash; *Method*.



```julia
eeg_tconv(eeg; kernel)
```

Perform convolution in the time domain.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel used for convolution

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_tconv!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_tconv!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_tconv!`** &mdash; *Method*.



```julia
eeg_tconv!(eeg; kernel)
```

Perform convolution in the time domain.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel used for convolution

<a id='NeuroJ.eeg_fconv-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_fconv-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_fconv`** &mdash; *Method*.



```julia
eeg_fconv(eeg, kernel)
```

Performs convolution in the time domain.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel for convolution

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_fconv!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_fconv!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_fconv!`** &mdash; *Method*.



```julia
eeg_fconv!(eeg, kernel)
```

Performs convolution in the time domain.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel for convolution

<a id='NeuroJ.eeg_pca-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_pca-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_pca`** &mdash; *Method*.



```julia
eeg_pca(eeg; n)
```

Calculates `n` first PCs for `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `n::Int64`: number of PCs

**Returns**

  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pc_var::Matrix{Float64}`: PC*VAR(1)..PC*VAR(n) × epoch

<a id='NeuroJ.eeg_pca!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_pca!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_pca!`** &mdash; *Method*.



```julia
eeg_pca!(eeg; n)
```

Calculates `n` first PCs for `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `n::Int64`: number of PCs

<a id='NeuroJ.eeg_ica-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_ica-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_ica`** &mdash; *Method*.



```julia
eeg_ica(eeg; <keyword arguments>)
```

Calculates `n` first ICs for `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `n::Int64`: number of ICs
  * `tol::Float64=1.0e-6`: tolerance for ICA
  * `iter::Int64=100`: maximum number of iterations
  * `f::Symbol=:tanh`: neg-entropy functor: :tanh, :gaus

**Returns**

  * `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
  * `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)

<a id='NeuroJ.eeg_ica!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_ica!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_ica!`** &mdash; *Method*.



```julia
eeg_ica!(eeg; <keyword arguments>)
```

Calculates `n` first ICs for `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `n::Int64`: number of ICs
  * `tol::Float64=1.0e-6`: tolerance for ICA
  * `iter::Int64=100`: maximum number of iterations
  * `f::Symbol=:tanh`: neg-entropy functor: :tanh, :gaus

<a id='NeuroJ.eeg_ica_reconstruct-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_ica_reconstruct-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_ica_reconstruct`** &mdash; *Method*.



```julia
eeg_ica_reconstruct(eeg; ica)
```

Reconstructs `eeg` signals using removal of `ica` ICA components.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `ica::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_ica_reconstruct!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_ica_reconstruct!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_ica_reconstruct!`** &mdash; *Method*.



```julia
eeg_ica_reconstruct!(eeg; ica)
```

Reconstructs `eeg` signals using removal of `ica` ICA components.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `ica::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove


<a id='EEG-analyze'></a>

<a id='EEG-analyze-1'></a>

## EEG analyze


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

<a id='NeuroJ.signal_plot_ica-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Vector{Float64}}' href='#NeuroJ.signal_plot_ica-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Vector{Float64}}'>#</a>
**`NeuroJ.signal_plot_ica`** &mdash; *Method*.



```julia
signal_plot_ica(t, ica; <keyword arguments>)
```

Plot `ica` components against time vector `t`.

**Arguments**

  * `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`: the time vector
  * `ica::Vector{Float64}`
  * `label::String=""`: channel label
  * `norm::Bool=true`: normalize the `ica` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits (-ylim:ylim)
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_ica-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Matrix{Float64}}' href='#NeuroJ.signal_plot_ica-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Matrix{Float64}}'>#</a>
**`NeuroJ.signal_plot_ica`** &mdash; *Method*.



```julia
signal_plot_ica(t, ica; <keyword arguments>)
```

Plots `ica` components.

**Arguments**

  * `t::Union{Vector{Float64}, Vector{Int64}, AbstractRange}`
  * `ica::Matrix{Float64}`
  * `labels::Vector{String}=[""]`: labels vector
  * `norm::Bool=true`: normalize the ICs prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_ica-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_ica-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_ica`** &mdash; *Method*.



```julia
eeg_plot_ica(eeg; <keyword arguments>)
```

Plots embedded ICs components.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Int64=1`: epoch number to display
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default 1 epoch or 20 seconds
  * `ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing`: which IC to plot, default all
  * `norm::Bool=true`: normalize the ICs prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_topo-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_topo-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_topo`** &mdash; *Method*.



```julia
eeg_plot_topo(eeg; <keyword arguments>)
```

Plot topographical view of `eeg` component.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `offset::Int64`: time (in samples) at which to plot
  * `len::Int64=0`: interpolation window, default 100 ms
  * `m::Symbol=:shepard`: interpolation method :shepard (Shepard), :mq (Multiquadratic), :tp (ThinPlate)
  * `c::Symbol=:amp`: component name (:ica, :pca, :amp, :power)
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange, Tuple, Nothing}=nothing`: component index, e.g. ICA number or frequency range
  * `norm::Bool=true`: convert power as dB
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: frequency limit for PSD and spectrogram
  * `head_labels::Bool=false`: plot head labels
  * `cb::Bool=false`: add color bars to plots
  * `cb_label::String=""`: color bar label
  * `average::Bool=true`: plot averaged signal and PSD
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_bands-Tuple{Vector{Float64}}' href='#NeuroJ.signal_plot_bands-Tuple{Vector{Float64}}'>#</a>
**`NeuroJ.signal_plot_bands`** &mdash; *Method*.



```julia
signal_plot_bands(signal; <keyword arguments>)
```

Plot absolute/relative band powers of `signal`.

**Arguments**

  * `signal::Vector{Float64}`
  * `fs::Int64`: sampling rate
  * `band:Vector{Symbols}=:all`: band name, e.g. :delta (see `eeg_band()`)
  * `type::Symbol`: plots absolute (:abs) or relative power (:rel)
  * `norm::Bool=true`: convert power to dB if true
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: other arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_bands-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_bands-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_bands`** &mdash; *Method*.



```julia
eeg_plot_bands(eeg; <keyword arguments>)
```

Plots `eeg` channels. If signal is multichannel, only channel amplitudes are plotted. For single-channel signal, the histogram, amplitude, power density and spectrogram are plotted.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}=1`: epochs to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channels to display
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `band:Vector{Symbols}=:all`: band name, e.g. :delta (see `eeg_band()`)
  * `type::Symbol`: plots absolute (:abs) or relative power (:rel)
  * `norm::Bool=true`: convert power to dB if true
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
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

<a id='NeuroJ.eeg_plot_save-Tuple{Plots.Plot{Plots.GRBackend}}' href='#NeuroJ.eeg_plot_save-Tuple{Plots.Plot{Plots.GRBackend}}'>#</a>
**`NeuroJ.eeg_plot_save`** &mdash; *Method*.



```julia
eeg_plot_save(p; file_name::String)
```

Saves plot as file (PDF/PNG/TIFF). File format is determined using `file_name` extension.

**Arguments**

  * `p::Plots.Plot{Plots.GRBackend}`
  * `file_name::String`


<a id='Signal'></a>

<a id='Signal-1'></a>

## Signal


<a id='Misc'></a>

<a id='Misc-1'></a>

## Misc


<a id='NSTIM'></a>

<a id='NSTIM-1'></a>

## NSTIM

<a id='NeuroJ.tes_dose-Tuple{Union{Float64, Int64}, Union{Float64, Int64}, Int64}' href='#NeuroJ.tes_dose-Tuple{Union{Float64, Int64}, Union{Float64, Int64}, Int64}'>#</a>
**`NeuroJ.tes_dose`** &mdash; *Method*.



```julia
tes_dose(current, pad_area, duration)
```

Converts `current`, `pad_area` and stimulation `duration` into `charge`, `current_density` and `charge_ density`.

**Arguments**

  * `current::Union{Int64, Float64}`: stimulation current [mA]
  * `pad_area::Union{Int64, Float64}`: electrode pad area [cm^2]
  * `duration::Int64`: stimulation duration [s]

**Returns**

  * `charge::Float64`: charge [C]
  * `current_density::Float64`: current density [A/m^2]
  * `charge_density::Float64`: delibvered charge density [kC/m^2]

**Source**

1. Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.

