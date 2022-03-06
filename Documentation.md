
<a id='NeuroJ.jl-Documentation'></a>

<a id='NeuroJ.jl-Documentation-1'></a>

# NeuroJ.jl Documentation


This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).


<a id='TOC'></a>

<a id='TOC-1'></a>

## TOC

- [NeuroJ.jl Documentation](index.md#NeuroJ.jl-Documentation)
    - [TOC](index.md#TOC)
    - [NeuroJ](index.md#NeuroJ)
    - [EEG io](index.md#EEG-io)
    - [EEG edit](index.md#EEG-edit)
    - [EEG process](index.md#EEG-process)
    - [EEG analyze](index.md#EEG-analyze)
    - [EEG plots](index.md#EEG-plots)
    - [Signal](index.md#Signal)
    - [Misc](index.md#Misc)
    - [NSTIM](index.md#NSTIM)
    - [Index](index.md#Index)


<a id='NeuroJ'></a>

<a id='NeuroJ-1'></a>

## NeuroJ

<a id='NeuroJ.neuroj_version-Tuple{}' href='#NeuroJ.neuroj_version-Tuple{}'>#</a>
**`NeuroJ.neuroj_version`** &mdash; *Method*.



```julia
neuroj_version()
```

Show NeuroJ and imported packages versions.

<a id='NeuroJ.neuroj_reload_plugins-Tuple{}' href='#NeuroJ.neuroj_reload_plugins-Tuple{}'>#</a>
**`NeuroJ.neuroj_reload_plugins`** &mdash; *Method*.



```julia
neuroj_reload_plugins()
```

Reload NeuroJ plugins. Plugins path is: `~/Documents/NeuroJ/plugins/`.


<a id='EEG-io'></a>

<a id='EEG-io-1'></a>

## EEG io

<a id='NeuroJ.eeg_import_edf-Tuple{String}' href='#NeuroJ.eeg_import_edf-Tuple{String}'>#</a>
**`NeuroJ.eeg_import_edf`** &mdash; *Method*.



```julia
eeg_import_edf(file_name; read_annotations, clean_labels)
```

Load EDF/EDFPlus file and return and `NeuroJ.EEG` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `read_annotations::Bool=true`: read annotations from EDF+ file (currently not implemented)
  * `clean_labels::Bool=true`: only keep channel names in channel labels

**Returns**

  * `eeg:EEG`

**Notes**

  * sampling_rate = n.samples / data.record.duration
  * gain = (physical*maximum - physical*minimum) / (digital*maximum - digital*minimum)
  * value = (value - digital*minimum ) * gain + physical*minimum

**Source**

Kemp B, Värri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3. 

<a id='NeuroJ.eeg_import_ced-Tuple{String}' href='#NeuroJ.eeg_import_ced-Tuple{String}'>#</a>
**`NeuroJ.eeg_import_ced`** &mdash; *Method*.



```julia
eeg_import_ced(file_name)
```

Load electrode positions from CED file.

**Arguments**

  * `file_name::String`

**Returns**

  * `sensors::DataFrame`

<a id='NeuroJ.eeg_import_locs-Tuple{String}' href='#NeuroJ.eeg_import_locs-Tuple{String}'>#</a>
**`NeuroJ.eeg_import_locs`** &mdash; *Method*.



```julia
eeg_import_locs(file_name)
```

Load electrode positions from LOCS file.

**Arguments**

  * `file_name::String`

**Returns**

  * `sensors::DataFrame`

<a id='NeuroJ.eeg_import_elc-Tuple{String}' href='#NeuroJ.eeg_import_elc-Tuple{String}'>#</a>
**`NeuroJ.eeg_import_elc`** &mdash; *Method*.



```julia
eeg_import_elc(file_name)
```

Load electrode positions from ELC file.

**Arguments**

  * `file_name::String`

**Returns**

  * `sensors::DataFrame`

<a id='NeuroJ.eeg_load_electrodes-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_load_electrodes-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_load_electrodes`** &mdash; *Method*.



```julia
eeg_load_electrodes(eeg; file_name)
```

Load electrode positions from `file_name` and return `NeuroJ.EEG` object with metadata: `:channel_locations`, `:loc_theta`, `:loc_radius`, `:loc_x`, `:loc_x`, `:loc_y`, `:loc_radius_sph`, `:loc_theta_sph`, `:loc_phi_sph`. Accepted formats:

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

Load electrode positions from `file_name` and set `eeg` metadata: `:channel_locations`, `:loc_theta`, `:loc_radius`, `:loc_x`, `:loc_x`, `:loc_y`, `:loc_radius_sph`, `:loc_theta_sph`, `:loc_phi_sph`. Accepted formats:

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

Load the `eeg` from `file_name` file (HDF5-based).

**Arguments**

  * `file_name::String`: file name

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_save-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_save-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_save`** &mdash; *Method*.



```julia
eeg_save(eeg; file_name, overwrite)
```

Save the `eeg` to `file_name` file (HDF5-based).

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

Export EEG data as CSV.

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

<a id='NeuroJ.eeg_delete_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_channel`** &mdash; *Method*.



```julia
eeg_delete_channel(eeg; channel)
```

Remove `channel` from the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to be removed, vector of numbers or range

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_delete_channel!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_channel!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_channel!`** &mdash; *Method*.



```julia
eeg_delete_channel!(eeg; channel)
```

Remove `channel` from the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to be removed

<a id='NeuroJ.eeg_keep_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_keep_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_keep_channel`** &mdash; *Method*.



```julia
eeg_keep_channel(eeg; channel)
```

Keep `channels` in the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to keep

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_keep_channel!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_keep_channel!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_keep_channel!`** &mdash; *Method*.



```julia
eeg_keep_channel!(eeg; channel)
```

Keep `channels` in the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to keep

<a id='NeuroJ.eeg_get_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_get_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_get_channel`** &mdash; *Method*.



```julia
eeg_get_channel(eeg; channel)
```

Returnthe `channel` index / name.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, String}`: channel name

**Returns**

  * `channel_idx::Int64`

<a id='NeuroJ.eeg_rename_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_rename_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_rename_channel`** &mdash; *Method*.



```julia
eeg_rename_channel(eeg; channel, new_name)
```

Renames the `eeg` `channel`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, String}`
  * `new_name::String`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_rename_channel!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_rename_channel!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_rename_channel!`** &mdash; *Method*.



```julia
eeg_rename_channel!(eeg; channel, new_name)
```

Renames the `eeg` `channel`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, String}`
  * `new_name::String`

<a id='NeuroJ.eeg_extract_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_extract_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_extract_channel`** &mdash; *Method*.



```julia
eeg_extract_channel(eeg; channel)
```

Extract `channel` number or name.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, String}`

**Returns**

  * `channel::Vector{Float64}`

<a id='NeuroJ.eeg_history-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_history-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_history`** &mdash; *Method*.



```julia
eeg_history(eeg)
```

Show processing history.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_labels-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_labels-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_labels`** &mdash; *Method*.



```julia
eeg_labels(eeg)
```

Returnlabels.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_sr-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_sr-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_sr`** &mdash; *Method*.



```julia
eeg_sr(eeg)
```

Returnsampling rate.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_channel_n-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_channel_n-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_channel_n`** &mdash; *Method*.



```julia
eeg_channel_n(eeg; type=:eeg)
```

Returnnumber of `eeg` channels of `type`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Vector{Symbol}[:all, :eeg, :ecg, :eog, :emg]`

**Returns**

  * `channel_n::Int64`

<a id='NeuroJ.eeg_epoch_n-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epoch_n-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epoch_n`** &mdash; *Method*.



```julia
eeg_epoch_n(eeg)
```

Returnnumber of `eeg` epochs.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `epoch_n::Int64`

<a id='NeuroJ.eeg_signal_len-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_signal_len-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_signal_len`** &mdash; *Method*.



```julia
eeg_signal_len(eeg)
```

Returnlength of `eeg` signal.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `signal_len::Int64`

<a id='NeuroJ.eeg_epoch_len-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epoch_len-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epoch_len`** &mdash; *Method*.



```julia
eeg_epoch_len(eeg)
```

Returnlength of `eeg` signal.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `epoch_len::Int64`

<a id='NeuroJ.eeg_info-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_info-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_info`** &mdash; *Method*.



```julia
eeg_info(eeg)
```

Show info.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_epochs-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epochs-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epochs`** &mdash; *Method*.



```julia
eeg_epochs(eeg; epoch_n=nothing, epoch_len=nothing, average=false)
```

Splits `eeg` into epochs.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch_n::Union{Int64, Nothing}`: number of epochs
  * `epoch_len::Union{Int64, Nothing}`: epoch length in samples
  * `average::Bool`: average all epochs, returnone averaged epoch; if false than returnarray of epochs, each row is one epoch

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_epochs!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epochs!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epochs!`** &mdash; *Method*.



```julia
eeg_epochs!(eeg; epoch_n=nothing, epoch_len=nothing, average=false)
```

Splits `eeg` into epochs.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch_n::Union{Int64, Nothing}`: number of epochs
  * `epoch_len::Union{Int64, Nothing}`: epoch length in samples
  * `average::Bool`: average all epochs, returnone averaged epoch; if false than returnarray of epochs, each row is one epoch

<a id='NeuroJ.eeg_extract_epoch-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_extract_epoch-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_extract_epoch`** &mdash; *Method*.



```julia
eeg_extract_epoch(eeg; epoch)
```

Extract the `epoch` epoch.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Int64`: epoch index

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_trim-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_trim-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_trim`** &mdash; *Method*.



```julia
eeg_trim(eeg:EEG; len, offset=0, from=:start, keep_epochs::Bool=true)
```

Remove `len` samples from the beginning + `offset` (`from` = :start, default) or end (`from` = :end) of the `eeg`.

**Arguments**

  * `eeg:EEG`
  * `len::Int64`: number of samples to remove
  * `offset::Int64`: offset from which trimming starts, only works for `from` = :start
  * `from::Symbol[:start, :end]`: trims from the signal start (default) or end
  * `keep_epochs::Bool`: remove epochs containing signal to trim (keep_epochs=true) or remove signal and remove epoching

**Returns**

  * `eeg:EEG`

<a id='NeuroJ.eeg_trim!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_trim!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_trim!`** &mdash; *Method*.



```julia
eeg_trim!(eeg:EEG; len, offset=0, from=:start, keep_epochs::Bool=true)
```

Remove `len` samples from the beginning + `offset` (`from` = :start, default) or end (`from` = :end) of the `eeg`.

**Arguments**

  * `eeg:EEG`
  * `len::Int64`: number of samples to remove
  * `offset::Int64`: offset from which trimming starts, only works for `from` = :start
  * `from::Symbol[:start, :end]`: trims from the signal start (default) or end
  * `keep_epochs::Bool`: remove epochs containing signal to trim (keep_epochs=true) or remove signal and remove epoching

<a id='NeuroJ.eeg_edit_header-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_edit_header-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_edit_header`** &mdash; *Method*.



```julia
eeg_edit_header(eeg; field, value)
```

Change value of `eeg` `field` to `value`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `field::Symbol`
  * `value::Any`

**Returns**

  * `eeg:EEG`

<a id='NeuroJ.eeg_edit_header!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_edit_header!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_edit_header!`** &mdash; *Method*.



```julia
eeg_edit_header!(eeg; field, value)
```

Change value of `eeg` `field` to `value`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `field::Symbol`
  * `value::Any`

**Returns**

  * `eeg:EEG`

<a id='NeuroJ.eeg_show_header-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_show_header-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_show_header`** &mdash; *Method*.



```julia
eeg_show_header(eeg)
```

Show keys and values of `eeg` header.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_delete_epoch-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_epoch-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_epoch`** &mdash; *Method*.



```julia
eeg_delete_epoch(eeg; epoch)
```

Remove `epoch` from the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch index to be removed, vector of numbers or range

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_delete_epoch!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_epoch!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_epoch!`** &mdash; *Method*.



```julia
eeg_delete_epoch!(eeg; epoch)
```

Remove `epoch` from the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch index to be removed, vector of numbers or range

<a id='NeuroJ.eeg_keep_epoch-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_keep_epoch-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_keep_epoch`** &mdash; *Method*.



```julia
eeg_keep_epoch(eeg; epoch)
```

Keep `epoch` in the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch index to keep, vector of numbers or range

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_keep_epoch!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_keep_epoch!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_keep_epoch!`** &mdash; *Method*.



```julia
eeg_keep_epoch!(eeg; epoch)
```

Keep `epoch` in the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch index to keep, vector of numbers or range

<a id='NeuroJ.eeg_detect_bad_epochs-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_detect_bad_epochs-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_detect_bad_epochs`** &mdash; *Method*.



```julia
eeg_detect_bad_epochs(eeg; method=[:flat, :rmse, :rmsd, :euclid, :p2p], ch_t)
```

Detect bad `eeg` epochs based on:

  * flat channel(s)
  * RMSE
  * RMSD
  * Euclidean distance
  * peak-to-peak amplitude

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p]`
  * `ch_t::Float64`: percentage of bad channels to mark the epoch as bad

**Returns**

  * `bad_epochs_idx::Vector{Int64}`

<a id='NeuroJ.eeg_delete_bad_epochs-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_bad_epochs-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_bad_epochs`** &mdash; *Method*.



```julia
eeg_check_bad_epochs(eeg; bad_epochs, confirm=true)
```

Delete bad `eeg` epochs.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `bad_epochs_idx::Vector{Int64}`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_delete_bad_epochs!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_bad_epochs!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_bad_epochs!`** &mdash; *Method*.



```julia
eeg_delete_bad_epochs!(eeg; bad_epochs, confirm=true)
```

Delete bad `eeg` epochs.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `bad_epochs_idx::Vector{Int64}`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_add_labels-Tuple{NeuroJ.EEG, Vector{String}}' href='#NeuroJ.eeg_add_labels-Tuple{NeuroJ.EEG, Vector{String}}'>#</a>
**`NeuroJ.eeg_add_labels`** &mdash; *Method*.



```julia
eeg_add_labels(eeg::NeuroJ.EEG, labels::Vector{String})
```

Add `labels` to `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `labels::Vector{String}`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_add_labels!-Tuple{NeuroJ.EEG, Vector{String}}' href='#NeuroJ.eeg_add_labels!-Tuple{NeuroJ.EEG, Vector{String}}'>#</a>
**`NeuroJ.eeg_add_labels!`** &mdash; *Method*.



```julia
eeg_add_labels!(eeg::NeuroJ.EEG, labels::Vector{String})
```

Add `labels` to `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `labels::Vector{String}`

<a id='NeuroJ.eeg_edit_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_edit_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_edit_channel`** &mdash; *Method*.



```julia
eeg_edit_channel(eeg; channel, field, value)
```

Edits `eeg` `channel` properties.

**Arguments**

  * `eeg:EEG`
  * `channel::Int64`
  * `field::Any`
  * `value::Any`

**Returns**

  * `eeg_new::NeuroJ.EEG`

<a id='NeuroJ.eeg_edit_channel!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_edit_channel!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_edit_channel!`** &mdash; *Method*.



```julia
eeg_edit_channel!(eeg; channel, field, value)
```

Edit `eeg` `channel` properties.

**Arguments**

  * `eeg:EEG`
  * `channel::Int64`
  * `field::Any`
  * `value::Any`

<a id='NeuroJ.eeg_keep_eeg_channels-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_keep_eeg_channels-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_keep_eeg_channels`** &mdash; *Method*.



```julia
eeg_keep_eeg_channels(eeg::NeuroJ.EEG)
```

Keep only EEG channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_keep_eeg_channels!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_keep_eeg_channels!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_keep_eeg_channels!`** &mdash; *Method*.



```julia
eeg_keep_eeg_channels!(eeg::NeuroJ.EEG)
```

Keep only EEG channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_list_components-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_list_components-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_list_components`** &mdash; *Method*.



```julia
eeg_list_components(eeg)
```

List `eeg` components.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `components::Vector{Symbol}`

<a id='NeuroJ.eeg_extract_component-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_extract_component-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_extract_component`** &mdash; *Method*.



```julia
eeg_extract_component(eeg, c)
```

Extract `component` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `component::Symbol`

**Returns**

  * `component::Any`

<a id='NeuroJ.eeg_delete_component-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_component-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_component`** &mdash; *Method*.



```julia
eeg_delete_component(eeg, c)
```

Delete `component` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `component::Symbol`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_delete_component!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_component!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_component!`** &mdash; *Method*.



```julia
eeg_delete_component!(eeg, c)
```

Delete `component` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `component::Symbol`

<a id='NeuroJ.eeg_reset_components-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reset_components-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reset_components`** &mdash; *Method*.



```julia
eeg_reset_components(eeg)
```

Reset `eeg` components.

**Arguments**

  * `eeg:EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_reset_components!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reset_components!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reset_components!`** &mdash; *Method*.



```julia
eeg_reset_components!(eeg)
```

Reset `eeg` components.

**Arguments**

  * `eeg:EEG`


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

<a id='NeuroJ.eeg_average!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_average!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_average!`** &mdash; *Method*.



```julia
eeg_average!(eeg)
```

Returns the average signal of all `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`

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

Downsample all channels of `eeg` to `new_sr` sampling frequency.

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

Downsample all channels of `eeg` to `new_sr` sampling frequency.

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

Reconstruct `eeg` signals using removal of `ica` ICA components.

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

Reconstruct `eeg` signals using removal of `ica` ICA components.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `ica::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove


<a id='EEG-analyze'></a>

<a id='EEG-analyze-1'></a>

## EEG analyze

<a id='NeuroJ.eeg_total_power-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_total_power-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_total_power`** &mdash; *Method*.



```julia
eeg_total_power(eeg)
```

Calculate total power of the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `stp::Vector{Float64}`

<a id='NeuroJ.eeg_total_power!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_total_power!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_total_power!`** &mdash; *Method*.



```julia
eeg_total_power!(eeg)
```

Calculate total power of the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_band_power-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_band_power-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_band_power`** &mdash; *Method*.



```julia
eeg_band_power(eeg; f)
```

Calculate absolute band power between frequencies `f[1]` and `f[2]` of the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `f::Tuple(Union(Int64, Float64}, Union(Int64, Float64}}`: lower and upper frequency bounds

**Returns**

  * `sbp::Vector{Float64}`

<a id='NeuroJ.eeg_cov-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_cov-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_cov`** &mdash; *Method*.



```julia
eeg_cov(eeg; norm=true)
```

Calculate covariance between all channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `norm::Bool`: normalize covariance

**Returns**

  * `cov_mat::Array{Float64, 3}`

<a id='NeuroJ.eeg_cov!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_cov!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_cov!`** &mdash; *Method*.



```julia
eeg_cov!(eeg; norm)
```

Calculate covariance between all channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `norm::Bool=true`: normalize covariance

<a id='NeuroJ.eeg_cor-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_cor-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_cor`** &mdash; *Method*.



```julia
eeg_cor(eeg)
```

Calculate correlation coefficients between all channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `cov_mat::Array{Float64, 3}`

<a id='NeuroJ.eeg_cor!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_cor!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_cor!`** &mdash; *Method*.



```julia
eeg_cor!(eeg)
```

Calculate correlation coefficients between all channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_autocov-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_autocov-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_autocov`** &mdash; *Method*.



```julia
eeg_autocov(eeg; lag, demean, norm)
```

Calculate autocovariance of each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `lag::Int64=1`: lags range is `-lag:lag`
  * `demean::Bool=false`: demean signal prior to analysis
  * `norm::Bool=false`: normalize autocovariance

**Returns**

  * `acov::Matrix{Float64}`
  * `lags::Vector{Float64}

<a id='NeuroJ.eeg_autocov!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_autocov!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_autocov!`** &mdash; *Method*.



```julia
eeg_autocov!(eeg; lag, demean, norm)
```

Calculate autocovariance of each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `lag::Int64=1`: lags range is `-lag:lag`
  * `demean::Bool=false`: demean signal prior to analysis
  * `norm::Bool=false`: normalize autocovariance

<a id='NeuroJ.eeg_crosscov-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_crosscov-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_crosscov`** &mdash; *Method*.



```julia
eeg_crosscov(eeg; lag, demean, norm)
```

Calculate cross-covariance of each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `lag::Int64=1`: lags range is `-lag:lag`
  * `demean::Bool=false`: demean signal prior to analysis
  * `norm::Bool=false`: normalize cross-covariance

**Returns**

  * `ccov::Matrix{Float64}`
  * `lags::Vector{Float64}

<a id='NeuroJ.eeg_crosscov!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_crosscov!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_crosscov!`** &mdash; *Method*.



```julia
eeg_crosscov!(eeg; lag, demean, norm)
```

Calculate cross-covariance of each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `lag::Int64=1`: lags range is `-lag:lag`
  * `demean::Bool=false`: demean signal prior to analysis
  * `norm::Bool=false`: normalize cross-covariance

<a id='NeuroJ.eeg_crosscov-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_crosscov-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_crosscov`** &mdash; *Method*.



```julia
eeg_crosscov(eeg1, eeg2; lag, demean, norm)
```

Calculate cross-covariance between `eeg1` and `eeg2` channels.

**Arguments**

  * `eeg1::NeuroJ.EEG`
  * `eeg2::NeuroJ.EEG`
  * `lag::Int64=1`: lags range is `-lag:lag`
  * `demean::Bool=false`: demean signal prior to analysis
  * `norm::Bool=false`: normalize cross-covariance

**Returns**

  * `ccov::Matrix{Float64}`
  * `lags::Vector{Float64}

<a id='NeuroJ.eeg_psd-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_psd-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_psd`** &mdash; *Method*.



```julia
eeg_psd(eeg; norm)
```

Calculate total power for each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `norm::Bool=false`: normalize do dB

**Returns**

  * `powers::Array{Float64, 3}`
  * `frequencies::Array{Float64, 3}`

<a id='NeuroJ.eeg_psd!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_psd!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_psd!`** &mdash; *Method*.



```julia
eeg_psd!(eeg; norm)
```

Calculate total power for each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `norm::Bool=false`: normalize do dB

<a id='NeuroJ.eeg_stationarity-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_stationarity-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_stationarity`** &mdash; *Method*.



```julia
eeg_stationarity(eeg; window, method)
```

Calculate stationarity.

**Arguments**

  * `eeg:EEG`
  * `window::Int64=10`: time window in samples
  * `method::Symbol=:euclid`: stationarity method: :mean, :var, :euclid, :hilbert

**Returns**

  * `stationarity::Union{Matrix{Float64}, Array{Float64, 3}}`

<a id='NeuroJ.eeg_stationarity!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_stationarity!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_stationarity!`** &mdash; *Method*.



```julia
eeg_stationarity!(eeg; window, method)
```

Calculate stationarity.

**Arguments**

  * `eeg:EEG`
  * `window::Int64=10`: time window in samples
  * `method::Symbol=:euclid`: stationarity method: :mean, :var, :euclid, :hilbert

<a id='NeuroJ.eeg_mi-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_mi-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_mi`** &mdash; *Method*.



```julia
eeg_mi(eeg)
```

Calculate mutual information between all channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `mi::Array{Float64, 3}`

<a id='NeuroJ.eeg_mi!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_mi!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_mi!`** &mdash; *Method*.



```julia
eeg_mi!(eeg)
```

Calculate mutual information between all channels of `eeg` and store into :mi component.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_mi-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_mi-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_mi`** &mdash; *Method*.



```julia
eeg_mi(eeg1, eeg2)
```

Calculate mutual information between all channels of `eeg1` and `eeg2`.

**Arguments**

  * `eeg1::NeuroJ.EEG`
  * `eeg2::NeuroJ.EEG`

**Returns**

  * `mi::Array{Float64, 3}`

<a id='NeuroJ.eeg_entropy-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_entropy-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_entropy`** &mdash; *Method*.



```julia
eeg_entropy(eeg)
```

Calculate entropy of all channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `entropy::Matrix{Float64}`

<a id='NeuroJ.eeg_entropy!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_entropy!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_entropy!`** &mdash; *Method*.



```julia
eeg_entropy!(eeg)
```

Calculate entropy of all channels of `eeg` and store into :entropy component.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_band-Tuple{Any}' href='#NeuroJ.eeg_band-Tuple{Any}'>#</a>
**`NeuroJ.eeg_band`** &mdash; *Method*.



```julia
eeg_band(eeg, band)
```

Return frequency limits for a `band` range.

**Arguments**

  * `eeg:EEG`
  * `band::Symbol`: name of band range: :delta, :theta, :alpha, :beta, :beta*high, :gamma, :gamma*1, :gamma*2, :gamma*lower, :gamma_higher. If lower or upper band frequency limit exceeds Nyquist frequency of `eeg`, than bound is truncated to `eeg` range.

**Returns**

  * `band_frequency::Tuple{Float64, Float64}`

<a id='NeuroJ.eeg_coherence-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_coherence-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_coherence`** &mdash; *Method*.



```julia
eeg_coherence(eeg1, eeg2)
```

Calculate coherence between all channels of `eeg1` and `eeg2`.

**Arguments**

  * `eeg1::NeuroJ.EEG`
  * `eeg2::NeuroJ.EEG`

**Returns**

  * `coherence::Union{Matrix{Float64}, Array{ComplexF64, 3}}`

<a id='NeuroJ.eeg_coherence-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_coherence-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_coherence`** &mdash; *Method*.



```julia
eeg_coherence(eeg; channel1, channel2, epoch1, epoch2)
```

Calculate coherence between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel1::Int64`
  * `channel2::Int64`
  * `epoch1::Int64`
  * `epoch2::Int64`

**Returns**

  * `coherence::Vector{ComplexF64}`

<a id='NeuroJ.eeg_freqs-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_freqs-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_freqs`** &mdash; *Method*.



```julia
eeg_freqs(eeg)
```

Return vector of frequencies and Nyquist frequency for `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `hz::Vector{Float64}`
  * `nyquist::Float64`

<a id='NeuroJ.eeg_freqs!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_freqs!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_freqs!`** &mdash; *Method*.



```julia
eeg_freqs!(eeg)
```

Return vector of frequencies and Nyquist frequency for `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_difference-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_difference-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_difference`** &mdash; *Method*.



```julia
eeg_difference(eeg1, eeg2; n, method)
```

Calculate mean difference and its 95% CI between `eeg1` and `eeg2`.

**Arguments**

  * `eeg1::NeuroJ.EEG`
  * `eeg2::NeuroJ.EEG`
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:absdiff`

      * `:absdiff`: maximum difference
      * `:diff2int`: integrated area of the squared difference

**Returns**

  * `signals_statistic::Matrix{Float64}`
  * `signals_statistic_single::Vector{Float64}`
  * `p::Vector{Float64}`

<a id='NeuroJ.eeg_pick-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_pick-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_pick`** &mdash; *Method*.



```julia
eeg_picks(eeg; pick)
```

Return `pick` of electrodes for `eeg` electrodes.

**Arguments**

  * `pick::Vector{Symbol}`

**Returns**

  * `channels::Vector{Int64}`

<a id='NeuroJ.eeg_epochs_stats-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epochs_stats-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epochs_stats`** &mdash; *Method*.



```julia
eeg_epochs_stats(eeg)
```

Calculate mean, median, standard deviation, variance and kurtosis of each `eeg` epoch.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `mean::Vector{Float64}`
  * `median::Vector{Float64}`
  * `sd::Vector{Float64}`
  * `var::Vector{Float64}`
  * `kurtosis::Vector{Float64}`

<a id='NeuroJ.eeg_epochs_stats!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epochs_stats!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epochs_stats!`** &mdash; *Method*.



```julia
eeg_epochs_stats!(eeg)
```

Calculate mean, median, standard deviation, variance and kurtosis of each `eeg` epoch and store these in `eeg` components: `:epochs_mean`, `:epochs_median`, `:epochs_sd`, `:epochs_var`, `:epochs_kurtosis`. 

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_spectrogram-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_spectrogram-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_spectrogram`** &mdash; *Method*.



```julia
eeg_spectrogram(eeg; norm, demean)
```

Return spectrogram of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `norm::Bool`=true: normalize powers to dB
  * `demean::Bool`=true: demean signal prior to analysis

**Returns**

  * `spec.power::Array{Float64, 3}`
  * `spec.freq::Matrix{Float64}`
  * `spec.time::Matrix{Float64}`

<a id='NeuroJ.eeg_spectrogram!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_spectrogram!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_spectrogram!`** &mdash; *Method*.



```julia
eeg_spectrogram!(eeg, norm, demean)
```

Calculate spectrogram of `eeg`. and sore in `eeg` components: `:spectrogram_pow`, `:spectrogram_frq`, `:spectrogram_t`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `norm::Bool=true`: normalize powers to dB
  * `demean::Bool=true`: demean signal prior to analysis

<a id='NeuroJ.eeg_spectrum-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_spectrum-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_spectrum`** &mdash; *Method*.



```julia
eeg_spectrum(eeg; pad)
```

Calculate FFT, amplitudes, powers and phases for each channel of the `eeg`. For `pad` > 0 channels are padded with 0s.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `pad::Int64=0`: number of 0s to pad

**Returns**

  * `fft::Array{ComplexF64, 3}`
  * `amplitudes::Array{Float64, 3}`
  * `powers::Array{Float64, 3}`
  * `phases::Array{Float64, 3}

<a id='NeuroJ.eeg_spectrum!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_spectrum!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_spectrum!`** &mdash; *Method*.



```julia
eeg_spectrum!(eeg; pad)
```

Calculate FFT, amplitudes, powers and phases for each channel of the `eeg` and store in `eeg` components: `:spectrum_fft`, `:spectrum_amp`, `:spectrum_pow`, `:spectrum_phase`. For `pad` > 0 channels are padded with 0s.

**Arguments**

  * `eeg::NeuroJ.EEG`: the signal
  * `pad::Int64`: pad channels `pad` zeros

<a id='NeuroJ.eeg_s2t-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_s2t-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_s2t`** &mdash; *Method*.



```julia
eeg_s2t(eeg; t)
```

Convert time `t` in samples to seconds using `eeg` sampling rate.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `t::Int64`: time in samples

**Returns**

  * `t_s::Float64`: time in seconds

<a id='NeuroJ.eeg_t2s-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_t2s-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_t2s`** &mdash; *Method*.



```julia
eeg_t2s(eeg; t)
```

Convert time `t` in seconds to samples using `eeg` sampling rate.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `t::Union{Int64, Float64}`: time in seconds

**Returns**

  * `t_s::Float64`: time in samples


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

<a id='NeuroJ.signal_derivative-Tuple{AbstractArray}' href='#NeuroJ.signal_derivative-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_derivative`** &mdash; *Method*.



```julia
signal_derivative(signal)
```

Returns the derivative of the `signal` with length same as the signal.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_der::Union{Vector{Int64}, Vector{Float64}}`

<a id='NeuroJ.signal_derivative-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_derivative-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_derivative`** &mdash; *Method*.



```julia
signal_derivative(signal)
```

Returns the derivative of each the `signal` channels with length same as the signal.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `s_der::Array{Float64, 3}`

<a id='NeuroJ.signal_total_power-Tuple{AbstractArray}' href='#NeuroJ.signal_total_power-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_total_power`** &mdash; *Method*.



```julia
signal_total_power(signal; fs)
```

Calculates total power for the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate

**Returns**

  * `stp::Float64`

<a id='NeuroJ.signal_total_power-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_total_power-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_total_power`** &mdash; *Method*.



```julia
signal_total_power(signal; fs)
```

Calculates total power for each the `signal` channels.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `fs::Int64`: sampling rate

**Returns**

  * `stp::Matrix{Float64}`

<a id='NeuroJ.signal_band_power-Tuple{AbstractArray}' href='#NeuroJ.signal_band_power-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_band_power`** &mdash; *Method*.



```julia
signal_band_power(signal; fs, f)
```

Calculates absolute band power between frequencies `f[1]` and `f[2]` for the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate of the signal
  * `f::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}`: lower and upper frequency bound

**Returns**

  * `sbp::Float64`

<a id='NeuroJ.signal_band_power-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_band_power-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_band_power`** &mdash; *Method*.



```julia
signal_band_power(signal; fs, f)
```

Calculates absolute band power between frequencies `f[1]` and `f[2]` for the `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}`: lower and upper frequency bound

**Returns**

  * `sbp::Matrix{Float64}`

<a id='NeuroJ.signal_make_spectrum-Tuple{AbstractArray}' href='#NeuroJ.signal_make_spectrum-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_make_spectrum`** &mdash; *Method*.



```julia
signal_make_spectrum(signal; fs)
```

Returns FFT and DFT sample frequencies for a DFT for the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate

**Returns**

  * `s_fft::Vector{ComplexF64}`
  * `s_sf::Vector{Float64}`

<a id='NeuroJ.signal_make_spectrum-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_make_spectrum-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_make_spectrum`** &mdash; *Method*.



```julia
signal_make_spectrum(signal; fs)
```

Returns FFT and DFT sample frequencies for a DFT for each the `signal` channels.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `fs::Int64`: sampling rate

**Returns**

  * `s_fft::Array{ComplexF64, 3}`
  * `s_sf::Array{Float64, 3}`

<a id='NeuroJ.signal_detrend-Tuple{AbstractArray}' href='#NeuroJ.signal_detrend-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_detrend`** &mdash; *Method*.



```julia
signal_detrend(signal; type=:linear)
```

Removes linear trend from the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `type::Symbol[:linear, :constant]`, optional

      * `linear`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `constant`: the mean of `signal` is subtracted

**Returns**

  * `s_detrended::Vector{Float64}`

<a id='NeuroJ.signal_detrend-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_detrend-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_detrend`** &mdash; *Method*.



```julia
signal_detrend(signal; type=:linear)
```

Removes linear trend for each the `signal` channels.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `type::Symbol[:linear, :constant]`, optional

      * `linear`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `constant`: the mean of `signal` is subtracted

**Returns**

  * `s_detrended::Array{Float64, 3}`

<a id='NeuroJ.signal_ci95-Tuple{Vector{Float64}}' href='#NeuroJ.signal_ci95-Tuple{Vector{Float64}}'>#</a>
**`NeuroJ.signal_ci95`** &mdash; *Method*.



```julia
signal_ci95(signal; n=3, method=:normal)
```

Calculates mean, std and 95% confidence interval for `signal`.

**Arguments**

  * `signal::Vector{Float64}`
  * `n::Int64`: number of bootstraps
  * `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

**Returns**

  * `s_m::Float64`: mean
  * `s_s::Float64`: standard deviation
  * `s_u::Float64`: upper 95% CI
  * `s_l::Float64`: lower 95% CI

<a id='NeuroJ.signal_ci95-Tuple{AbstractArray}' href='#NeuroJ.signal_ci95-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_ci95`** &mdash; *Method*.



```julia
signal_ci95(signal; n=3, method=:normal)
```

Calculates mean, std and 95% confidence interval for each the `signal` channels.

**Arguments**

  * `signal::AbstractArray`
  * `n::Int64`: number of bootstraps
  * `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

**Returns**

  * `s_m::Vector{Float64}`: mean
  * `s_s::Vector{Float64}`: standard deviation
  * `s_u::Vector{Float64}`: upper 95% CI
  * `s_l::Vector{Float64}`: lower 95% CI

<a id='NeuroJ.signal_ci95-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_ci95-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_ci95`** &mdash; *Method*.



```julia
signal_ci95(signal; n::=3, method=:normal)
```

Calculates mean, std and 95% confidence interval for each the `signal` channels.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `n::Int64`: number of bootstraps
  * `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

**Returns**

  * `s_m::Matrix{Float64}`: mean
  * `s_s::Matrix{Float64}`: standard deviation
  * `s_u::Matrix{Float64}`: upper 95% CI
  * `s_l::Matrix{Float64}`: lower 95% CI

<a id='NeuroJ.signal_mean-Tuple{Vector{Float64}, Vector{Float64}}' href='#NeuroJ.signal_mean-Tuple{Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroJ.signal_mean`** &mdash; *Method*.



```julia
signal_mean(signal1, signal2)
```

Calculates mean and 95% confidence interval for 2 signals.

**Arguments**

  * `signal1::Vector{Float64}`
  * `signal2:Vector{Float64}`

**Returns**

  * `s_m::Float64`: mean
  * `s_s::Float64`: standard deviation
  * `s_u::Float64`: upper 95% CI
  * `s_l::Float64`: lower 95% CI

<a id='NeuroJ.signal_mean-Tuple{Array{Float64, 3}, Array{Float64, 3}}' href='#NeuroJ.signal_mean-Tuple{Array{Float64, 3}, Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_mean`** &mdash; *Method*.



```julia
signal_mean(signal1, signal2)
```

Calculates mean and 95% confidence interval for 2 signals.

**Arguments**

  * `signal1::Array{Float64, 3}`
  * `signal2:Array{Float64, 3}`

**Returns**

  * `s_m::Matrix{Float64}`
  * `s_s::Matrix{Float64}`
  * `s_u::Matrix{Float64}`
  * `s_l::Matrix{Float64}`

<a id='NeuroJ.signal_difference-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.signal_difference-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.signal_difference`** &mdash; *Method*.



```julia
signal_difference(signal1, signal2; n=3, method=:absdiff)
```

Calculates mean difference and 95% confidence interval for 2 signals.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`
  * `n::Int64`: number of bootstraps
  * `method::Symbol[:absdiff, :diff2int]`

      * `:absdiff`: maximum difference
      * `:diff2int`: integrated area of the squared difference

**Returns**

  * `s_statistic::Vector{Float64}`
  * `s_statistic_single::Float64`
  * `p::Float64`

<a id='NeuroJ.signal_difference-Tuple{Array{Float64, 3}, Array{Float64, 3}}' href='#NeuroJ.signal_difference-Tuple{Array{Float64, 3}, Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_difference`** &mdash; *Method*.



```julia
signal_difference(signal1, signal2; n=3, method=:absdiff)
```

Calculates mean difference and 95% confidence interval for 2 signals.

**Arguments**

  * `signal1::Array{Float64, 3}`
  * `signal2:Array{Float64, 3}`
  * `n::Int64`: number of bootstraps
  * `method::Symbol[:absdiff, :diff2int]`

      * `:absdiff`: maximum difference
      * `:diff2int`: integrated area of the squared difference

**Returns**

  * `s_statistic::Matrix{Float64}`
  * `s_statistic_single::Vector{Float64}`
  * `p::Vector{Float64}`

<a id='NeuroJ.signal_autocov-Tuple{AbstractArray}' href='#NeuroJ.signal_autocov-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_autocov`** &mdash; *Method*.



signal_autocov(signal; lag=1, demean=false, norm=false)

Calculates autocovariance of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `lag::Int64`: lags range is `-lag:lag`
  * `demean::Bool`: demean `signal` prior to calculations
  * `norm::Bool`: normalize autocovariance

**Returns**

  * `acov::Vector{Float64}`
  * `lags::Vector{Int64}`

<a id='NeuroJ.signal_autocov-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_autocov-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_autocov`** &mdash; *Method*.



signal_autocov(signal; lag=1, demean=false, norm=false)

Calculates autocovariance of each the `signal` channels.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `lag::Int64`: lags range is `-lag:lag`
  * `demean::Bool`: demean signal prior to analysis
  * `norm::Bool`: normalize autocovariance

**Returns**

  * `acov::Matrix{Float64}`
  * `lags::Vector{Int64}`

<a id='NeuroJ.signal_crosscov-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.signal_crosscov-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.signal_crosscov`** &mdash; *Method*.



signal_crosscov(signal1, signal2; lag=1, demean=false, norm=false)

Calculates cross-covariance between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`
  * `lag::Int64`: lags range is `-lag:lag`
  * `demean::Bool`: demean signal prior to analysis
  * `norm::Bool`: normalize cross-covariance

**Returns**

  * `ccov::Vector{Float64}`
  * `lags::Vector{Int64}`

<a id='NeuroJ.signal_crosscov-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_crosscov-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_crosscov`** &mdash; *Method*.



signal_crosscov(signal; lag=1, demean=false, norm=false)

Calculates cross-covariance between all channels in the `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`: the signal
  * `lag::Int64`: lags range is `-lag:lag`
  * `demean::Bool`: demean `signal` prior to analysis
  * `norm::Bool`: normalize cross-covariance

**Returns**

  * `ccov::Array{Float64, 3}`
  * `lags::Vector{Int64}`

<a id='NeuroJ.signal_crosscov-Tuple{Array{Float64, 3}, Array{Float64, 3}}' href='#NeuroJ.signal_crosscov-Tuple{Array{Float64, 3}, Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_crosscov`** &mdash; *Method*.



signal_crosscov(signal1, signal2; lag=1, demean=false, norm=false)

Calculates cross-covariance between same channels in `signal1` and `signal2`.

**Arguments**

  * `signal1::Array{Float64, 3}`
  * `signal2::Array{Float64, 3}`
  * `lag::Int64`: lags range is `-lag:lag`
  * `demean::Bool`: demean signal prior to analysis
  * `norm::Bool`: normalize cross-covariance

**Returns**

  * `ccov::Array{Float64, 3}`
  * `lags::Vector{Int64}`

<a id='NeuroJ.signal_spectrum-Tuple{AbstractArray}' href='#NeuroJ.signal_spectrum-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_spectrum`** &mdash; *Method*.



```julia
signal_spectrum(signal; pad=0)
```

Calculates FFT, amplitudes, powers and phases of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `pad::Int64`: pad the `signal` with `pad` zeros

**Returns**

  * `fft::Vector(ComplexF64}`
  * `amplitudes::Vector{Float64}`
  * `powers::Vector{Float64}`
  * `phases::Vector{Float64}

<a id='NeuroJ.signal_spectrum-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_spectrum-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_spectrum`** &mdash; *Method*.



```julia
signal_spectrum(signal; pad=0)
```

Calculates FFT, amplitudes, powers and phases for each channel of the `signal` matrix.

**Arguments**

  * `signal::Array{Float64, 3}`: the signal
  * `pad::Int64`: pad the `signal` with `pad` zeros

**Returns**

  * `fft::Array{ComplexF64, 3}`
  * `amplitudes::Array{Float64, 3}`
  * `powers::Array{Float64, 3}`
  * `phases::Array{Float64, 3}

<a id='NeuroJ.signal_epochs-Tuple{Vector{Float64}}' href='#NeuroJ.signal_epochs-Tuple{Vector{Float64}}'>#</a>
**`NeuroJ.signal_epochs`** &mdash; *Method*.



```julia
signal_epochs(signal; epoch_n, epoch_len, average=true)
```

Splits `signal` into epochs.

**Arguments**

  * `signal::Vector{Float64}`
  * `epoch_n::Union{Int64, Nothing}`: number of epochs
  * `epoch_len::Union{Int64, Nothing}`: epoch length in samples
  * `average::Bool`: average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

**Returns**

  * `epochs::Matrix{Float64}`

<a id='NeuroJ.signal_epochs-Tuple{Matrix{Float64}}' href='#NeuroJ.signal_epochs-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.signal_epochs`** &mdash; *Method*.



```julia
signal_epochs(signal; epoch_n=nothing, epoch_len=nothing, average=true)
```

Splits `signal` into epochs.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `epoch_n::Union{Int64, Nothing}`: number of epochs
  * `epoch_len::Union{Int64, Nothing}`: epoch length in samples
  * `average::Bool`: average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

**Returns**

  * `epochs::Array{Float64, 3}`

<a id='NeuroJ.signal_delete_channel-Tuple{Matrix{Float64}}' href='#NeuroJ.signal_delete_channel-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.signal_delete_channel`** &mdash; *Method*.



```julia
signal_delete_channel(signal; channel)
```

Removes `channel` from the `signal`.

**Arguments**

  * `signal::Matrix{Float64}`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel to be removed, vector of numbers or range

**Returns**

  * `signal::Matrix{Float64}`

<a id='NeuroJ.signal_delete_channel-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_delete_channel-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_delete_channel`** &mdash; *Method*.



```julia
signal_delete_channel(signal; channel)
```

Removes `channel` from the `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel to be removed, vector of numbers or range

**Returns**

  * `signal::Matrix{Float64}`

<a id='NeuroJ.signal_reference_channel-Tuple{Matrix{Float64}}' href='#NeuroJ.signal_reference_channel-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.signal_reference_channel`** &mdash; *Method*.



```julia
signal_reference_channel(signal, reference)
```

Re-references channels of the `signal` to specific signal channel.

**Arguments**

  * `signal::Matrix{Float64}`
  * `reference::Union{Int64, Vector{Int64}, AbstractRange}}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference

**Returns**

  * `s_referenced::Matrix{Float64}`

<a id='NeuroJ.signal_reference_channel-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_reference_channel-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_reference_channel`** &mdash; *Method*.



```julia
signal_reference_channel(signal, channel)
```

Re-references channels of the `signal` to specific signal channel.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference

**Returns**

  * `s_referenced::Matrix{Float64}`

<a id='NeuroJ.signal_reference_car-Tuple{Matrix{Float64}}' href='#NeuroJ.signal_reference_car-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.signal_reference_car`** &mdash; *Method*.



```julia
signal_reference_car(signal)
```

Re-references channels of the `signal` to common average reference.

**Arguments**

  * `signal::Matrix{Float64}`

**Returns**

  * `s_referenced::Matrix{Float64}`

<a id='NeuroJ.signal_reference_car-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_reference_car-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_reference_car`** &mdash; *Method*.



```julia
signal_reference_car(signal)
```

Re-references channels of the `signal` to common average reference.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `s_referenced::Array{Float64, 3}`

<a id='NeuroJ.signal_taper-Tuple{AbstractArray}' href='#NeuroJ.signal_taper-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_taper`** &mdash; *Method*.



```julia
signal_taper(signal; taper)
```

Taper the `signal` with `taper`.

**Arguments**

  * `signal::AbstractArray`
  * `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

**Returns**

  * `s_tapered::Vector{Union{Float64, ComplexF64}}`

<a id='NeuroJ.signal_taper-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_taper-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_taper`** &mdash; *Method*.



```julia
signal_taper(signal; taper)
```

Tapers channels of the `signal` with `taper`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

**Returns**

  * `s_tapered::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`

<a id='NeuroJ.signal_demean-Tuple{AbstractArray}' href='#NeuroJ.signal_demean-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_demean`** &mdash; *Method*.



```julia
signal_demean(signal)
```

Removes mean value (DC offset) from the `signal`.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_demeaned::Vector{Float64}`

<a id='NeuroJ.signal_demean-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_demean-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_demean`** &mdash; *Method*.



```julia
signal_demean(signal)
```

Removes mean value (DC offset) for each the `signal` channels.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `s_demeaned::Array{Float64, 3}`

<a id='NeuroJ.signal_normalize_zscore-Tuple{AbstractArray}' href='#NeuroJ.signal_normalize_zscore-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_normalize_zscore`** &mdash; *Method*.



```julia
signal_normalize_zscore(signal)
```

Normalize (by z-score) `signal`.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroJ.signal_normalize_zscore-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_normalize_zscore-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_normalize_zscore`** &mdash; *Method*.



```julia
signal_normalize_zscore(signal)
```

Normalize (by z-score) each the `signal` channel.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `s_normalized::Array{Float64, 3}`

<a id='NeuroJ.signal_normalize_minmax-Tuple{AbstractArray}' href='#NeuroJ.signal_normalize_minmax-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_normalize_minmax`** &mdash; *Method*.



```julia
signal_normalize_minmax(signal)
```

Normalize `signal` in [-1, +1].

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroJ.signal_normalize_minmax-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_normalize_minmax-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_normalize_minmax`** &mdash; *Method*.



```julia
signal_normalize_minmax(signal)
```

Normalize each the `signal` channel in [-1, +1].

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `s_normalized::Array{Float64, 3}`

<a id='NeuroJ.signal_cov-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.signal_cov-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.signal_cov`** &mdash; *Method*.



signal_cov(signal1, signal2; norm=true)

Calculates covariance between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`
  * `norm::Bool`: normalize covariance

**Returns**

  * `cov_mat::Matrix{Float64}`

<a id='NeuroJ.signal_cov-Tuple{AbstractArray}' href='#NeuroJ.signal_cov-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_cov`** &mdash; *Method*.



signal_cov(signal; norm=true)

Calculates covariance between all channels of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `norm::Bool`: normalize covariance

**Returns**

  * `cov_mat::Matrix{Float64}`

<a id='NeuroJ.signal_cov-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_cov-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_cov`** &mdash; *Method*.



signal_cov(signal; norm=true)

Calculates covariance between all channels of the `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `norm::Bool`: normalize covariance

**Returns**

  * `cov_mat::Array{Float64, 3}`

<a id='NeuroJ.signal_cor-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_cor-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_cor`** &mdash; *Method*.



signal_cor(signal)

Calculates correlation coefficients between all channels of the `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `cor_mat::Array{Float64, 3}`

<a id='NeuroJ.signal_add_noise-Tuple{AbstractArray}' href='#NeuroJ.signal_add_noise-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_add_noise`** &mdash; *Method*.



```julia
signal_add_noise(signal)
```

Adds random noise to the `signal`.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_noisy::Vector{Float64}`

<a id='NeuroJ.signal_add_noise-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_add_noise-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_add_noise`** &mdash; *Method*.



```julia
signal_add_noise(signal)
```

Add random noise to the `signal` channels.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `s_noisy::Array{Float64, 3}`

<a id='NeuroJ.signal_upsample-Tuple{AbstractArray}' href='#NeuroJ.signal_upsample-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_upsample`** &mdash; *Method*.



```julia
signal_upsample(signal; t, new_sr)
```

Upsamples `signal` to `new_sr` sampling frequency.

**Arguments**

  * `signal::AbstractArray`
  * `t::AbstractRange`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `s_upsampled::Vector{Float64}`
  * `t_upsampled::AbstractRange`

<a id='NeuroJ.signal_upsample-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_upsample-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_upsample`** &mdash; *Method*.



```julia
signal_upsample(signal; t, new_sr)
```

Upsamples all channels of `signal` to `new_sr` sampling frequency.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `t::AbstractRange`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `s_upsampled::Array{Float64, 3}`
  * `t_upsampled::AbstractRange`

<a id='NeuroJ.signal_tconv-Tuple{AbstractArray}' href='#NeuroJ.signal_tconv-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_tconv`** &mdash; *Method*.



```julia
signal_tconv(signal; kernel)
```

Performs convolution in the time domain between `signal` and `kernel`.

**Arguments**

  * `signal::AbstractArray`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

**Returns**

  * `s_conv::Union{Vector{Float64}, Vector{ComplexF64}}`

<a id='NeuroJ.signal_tconv-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_tconv-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_tconv`** &mdash; *Method*.



```julia
signal_tconv(signal; kernel)
```

Performs convolution in the time domain between `signal` and `kernel`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

**Returns**

  * `s_conv::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`

<a id='NeuroJ.signal_filter-Tuple{AbstractArray}' href='#NeuroJ.signal_filter-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_filter`** &mdash; *Method*.



```julia
signal_filter(signal; fprototype, ftype=nothing, cutoff, fs, order, rp, rs, dir=:twopass, d=1, window=nothing)
```

Filters `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fprototype::Symbol[:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir]`: filter prototype:

      * `:mavg`: moving average (with threshold and/or weight window)
      * `:mmed`: moving median (with threshold and/or weight window)
      * `:poly`: polynomial of `order` order
  * `ftype::Union{Symbol[:lp, :hp, :bp, :bs], Nothing}`: filter type
  * `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64`: filter order
  * `rp::Union{Int64, Float64}`: dB ripple in the passband
  * `rs::Union{Int64, Float64}`: dB attentuation in the stopband
  * `dir:Symbol[:onepass, :onepass_reverse, :twopass]`: filter direction
  * `d::Int64`: window length for mean average and median average filter
  * `t::Union{Int64, Float64}`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
  * `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

**Returns**

  * `s_filtered::Vector{Float64}`

<a id='NeuroJ.signal_filter-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_filter-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_filter`** &mdash; *Method*.



```julia
signal_filter(signal; fprototype, ftype=nothing, cutoff, fs, order, rp, rs, dir=:twopass, d=1, window=nothing)
```

Filters `signal` using zero phase distortion filter.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `fprototype::Symbol[:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir]`: filter prototype:

      * `:mavg`: moving average (with threshold and/or weight window)
      * `:mmed`: moving median (with threshold and/or weight window)
      * `:poly`: polynomial of `order` order
  * `ftype::Union{Symbol[:lp, :hp, :bp, :bs], Nothing}`: filter type
  * `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64`: filter order
  * `rp::Union{Int64, Float64}`: dB ripple in the passband
  * `rs::Union{Int64, Float64}`: dB attentuation in the stopband
  * `dir:Symbol[:onepass, :onepass_reverse, :twopass]`: filter direction
  * `d::Int64`: window length for mean average and median average filter
  * `t::Union{Int64, Float64}`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
  * `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

**Returns**

  * `s_filtered::Array{Float64, 3}`

<a id='NeuroJ.signal_downsample-Tuple{AbstractArray}' href='#NeuroJ.signal_downsample-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_downsample`** &mdash; *Method*.



```julia
signal_downsample(signal; t, new_sr)
```

Downsamples the`signal` to `new_sr` sampling frequency.

**Arguments**

  * `signal::AbstractArray`
  * `t::AbstractRange`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `s_downsampled::Vector{Float64}`
  * `t_downsampled::Vector{Float64}`

<a id='NeuroJ.signal_downsample-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_downsample-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_downsample`** &mdash; *Method*.



```julia
signal_downsample(signal; t, new_sr)
```

Downsamples all channels of the`signal` to `new_sr` sampling frequency.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `new_sr::Int64`: new sampling rate
  * `t::AbstractRange`

**Returns**

  * `s_downsampled::Array{Float64, 3}`
  * `t_downsampled::AbstractRange`

<a id='NeuroJ.signal_psd-Tuple{AbstractArray}' href='#NeuroJ.signal_psd-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_psd`** &mdash; *Method*.



```julia
signal_psd(signal; fs, norm=false)
```

Calculates power spectrum density of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate
  * `norm::Bool`: normalize do dB

**Returns**

  * `psd_pow::Vector{Float64}`
  * `psd_frq::Vector{Float64}`

<a id='NeuroJ.signal_psd-Tuple{Matrix{Float64}}' href='#NeuroJ.signal_psd-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.signal_psd`** &mdash; *Method*.



```julia
signal_psd(signal; fs, norm=false)
```

Calculates power spectrum density for each the `signal` channels.

**Arguments**

  * `signal::Matrix{Float64}`
  * `fs::Int64`: sampling rate
  * `norm::Bool`: normalize do dB

**Returns**

  * `psd_pow::Matrix{Float64}`
  * `psd_frq::Matrix{Float64}`

<a id='NeuroJ.signal_psd-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_psd-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_psd`** &mdash; *Method*.



```julia
signal_psd(signal; fs, norm=false)
```

Calculates power spectrum density for each the `signal` channels.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `fs::Int64` sampling rate
  * `norm::Bool`: normalize do dB

**Returns**

  * `psd_pow::Array{Float64, 3}`
  * `psd_frq::Array{Float64, 3}`

<a id='NeuroJ.signal_stationarity_hilbert-Tuple{AbstractArray}' href='#NeuroJ.signal_stationarity_hilbert-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_stationarity_hilbert`** &mdash; *Method*.



```julia
signal_stationarity_hilbert(signal::Vector{Float64})
```

Calculates phase stationarity using Hilbert transformation.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `phase_stationarity::Vector{Float64}`

<a id='NeuroJ.signal_stationarity_mean-Tuple{AbstractArray}' href='#NeuroJ.signal_stationarity_mean-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_stationarity_mean`** &mdash; *Method*.



```julia
signal_stationarity_mean(signal)
```

Calculates mean stationarity.

**Arguments**

  * `signal::AbstractArray`
  * `window::Int64`: time window in samples

**Returns**

  * `mean_stationarity::Vector{Float64}`

<a id='NeuroJ.signal_stationarity_var-Tuple{AbstractArray}' href='#NeuroJ.signal_stationarity_var-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_stationarity_var`** &mdash; *Method*.



```julia
signal_stationarity_var(signal)
```

Calculates variance stationarity.

**Arguments**

  * `signal::AbstractArray`
  * `window::Int64`: time window in samples

**Returns**

  * `var_stationarity::Vector{Float64}`

<a id='NeuroJ.signal_stationarity-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_stationarity-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_stationarity`** &mdash; *Method*.



```julia
signal_stationarity(signal; window=10, method=:euclid)
```

Calculates stationarity.

**Arguments**

  * `signal:Array{Float64, 3}`
  * `window::Int64`: time window in samples
  * `method::Symbol[:mean, :var, :euclid, :hilbert]

**Returns**

  * `stationarity::Union{Matrix{Float64}, Array{Float64, 3}}`

<a id='NeuroJ.signal_trim-Tuple{AbstractArray}' href='#NeuroJ.signal_trim-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_trim`** &mdash; *Method*.



```julia
signal_trim(signal; len::Int64)
```

Removes `len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `len::Int64`: trimming length in samples
  * `offset::Int64`: offset from which trimming starts, only works for `from` = :start
  * `from::Symbol[:start, :end]

**Returns**

  * `s_trimmed::Vector{Float64}`

<a id='NeuroJ.signal_trim-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_trim-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_trim`** &mdash; *Method*.



```julia
signal_trim(signal; len, offset=0, from=:start)
```

Removes `len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `len::Int64`: number of samples to remove
  * `offset::Int64`: offset from which trimming starts, only works for `from` = :start
  * `from::Symbol[:start, :end]`

**Returns**

  * `s_trimmed::Array{Float64, 3}`

<a id='NeuroJ.signal_mi-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.signal_mi-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.signal_mi`** &mdash; *Method*.



```julia
signal_mi(signal1::Vector{Float64}, signal2::Vector{Float64})
```

Calculates mutual information between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

  * `mi::Float64`

<a id='NeuroJ.signal_mi-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_mi-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_mi`** &mdash; *Method*.



```julia
signal_mi(signal)
```

Calculates mutual information between each the `signal` channels.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `mi::Array{Float64, 3}`

<a id='NeuroJ.signal_mi-Tuple{Array{Float64, 3}, Array{Float64, 3}}' href='#NeuroJ.signal_mi-Tuple{Array{Float64, 3}, Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_mi`** &mdash; *Method*.



```julia
signal_mi(signal1, signal2)
```

Calculates mutual information between each the `signal1` and `signal2` channels.

**Arguments**

  * `signal1::Array{Float64, 3}`
  * `signal2::Array{Float64, 3}`

**Returns**

  * `mi::Array{Float64, 3}`

<a id='NeuroJ.signal_entropy-Tuple{AbstractArray}' href='#NeuroJ.signal_entropy-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_entropy`** &mdash; *Method*.



```julia
signal_entropy(signal)
```

Calculates entropy of `signal`.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `ent::Float64`

<a id='NeuroJ.signal_entropy-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_entropy-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_entropy`** &mdash; *Method*.



```julia
signal_entropy(signal)
```

Calculates mutual information between each the `signal1` and `signal2` channels.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `s_entropy::Matrix{Float64}`

<a id='NeuroJ.signal_average-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.signal_average-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.signal_average`** &mdash; *Method*.



```julia
signal_average(signal1, signal2)
```

Averages `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

  * `s_averaged::Vector{Float64}`

<a id='NeuroJ.signal_average-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_average-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_average`** &mdash; *Method*.



```julia
signal_average(signal)
```

Averages all channels of `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `s_averaged::Array{Float64, 3}`

<a id='NeuroJ.signal_average-Tuple{Array{Float64, 3}, Array{Float64, 3}}' href='#NeuroJ.signal_average-Tuple{Array{Float64, 3}, Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_average`** &mdash; *Method*.



```julia
signal_average(signal1, signal2)
```

Averages `signal1` and `signal2`.

**Arguments**

  * `signal1::Array{Float64, 3}`
  * `signal2::Array{Float64, 3}`

**Returns**

  * `s_averaged::Array{Float64, 3}`

<a id='NeuroJ.signal_coherence-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.signal_coherence-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.signal_coherence`** &mdash; *Method*.



```julia
signal_coherence(signal1, signal2)
```

Calculates coherence between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

  * `coherence::Vector{ComplexF64}`

<a id='NeuroJ.signal_coherence-Tuple{Matrix{Float64}, Matrix{Float64}}' href='#NeuroJ.signal_coherence-Tuple{Matrix{Float64}, Matrix{Float64}}'>#</a>
**`NeuroJ.signal_coherence`** &mdash; *Method*.



```julia
signal_coherence(signal1, signal2)
```

Calculates coherence between `signal1` and `signal2`.

**Arguments**

  * `signal1::Matrix{Float64}`
  * `signal2::Matrix{Float64}`

**Returns**

  * `coherence::Matrix{Float64}`

<a id='NeuroJ.signal_coherence-Tuple{Array{Float64, 3}, Array{Float64, 3}}' href='#NeuroJ.signal_coherence-Tuple{Array{Float64, 3}, Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_coherence`** &mdash; *Method*.



```julia
signal_coherence(signal1, signal2)
```

Calculates coherence between `signal1` and `signal2`.

**Arguments**

  * `signal1::Array{Float64, 3}`
  * `signal2::Array{Float64, 3}`

**Returns**

  * `coherence::Array{ComplexF64, 3}`

<a id='NeuroJ.signal_pca-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_pca-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_pca`** &mdash; *Method*.



```julia
signal_pca(signal, n)
```

Calculates `n` first PCs for `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `n::Int64`: number of PCs

**Returns**

  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pc_var::Matrix{Float64}`: PC*VAR(1)..PC*VAR(n) × epoch

<a id='NeuroJ.signal_fconv-Tuple{AbstractArray}' href='#NeuroJ.signal_fconv-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_fconv`** &mdash; *Method*.



```julia
signal_fconv(signal; kernel)
```

Performs convolution in the frequency domain between `signal` and `kernel`.

**Arguments**

  * `signal::AbstractArray`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

**Returns**

  * `s_conv::Vector{ComplexF64}`

<a id='NeuroJ.signal_fconv-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_fconv-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_fconv`** &mdash; *Method*.



```julia
signal_fconv(signal; kernel)
```

Performs convolution in the frequency domain between `signal` and `kernel`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

**Returns**

  * `s_conv::Array{ComplexF64, 3}`

<a id='NeuroJ.signal_ica-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_ica-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_ica`** &mdash; *Method*.



```julia
signal_ica(signal, n, tol=1.0e-6, iter=100, f=:tanh)
```

Calculates `n` first ICs for `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `n::Int64`: number of PCs
  * `tol::Float64`: tolerance for ICA
  * `iter::Int64`: maximum number of iterations
  * `f::Symbol[:tanh, :gaus]`: neg-entropy functor

**Returns**

  * `ic::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
  * `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch

<a id='NeuroJ.signal_ica_reconstruct-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_ica_reconstruct-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_ica_reconstruct`** &mdash; *Method*.



```julia
signal_ica_reconstruct(signal, ic_activations, ic_mw, ic_v)
```

Reconstructs `signal` using removal of `ic_v` ICA components.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `ic_activation::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
  * `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
  * `ic_v::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

**Returns**

  * `signal_reconstructed::Array{Float64, 3}`

<a id='NeuroJ.signal_epochs_stats-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_epochs_stats-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_epochs_stats`** &mdash; *Method*.



```julia
signal_epochs_var(signal)
```

Calculates variance for all `signal` epochs.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `var::Vector{Float64}`

<a id='NeuroJ.signal_spectrogram-Tuple{AbstractArray}' href='#NeuroJ.signal_spectrogram-Tuple{AbstractArray}'>#</a>
**`NeuroJ.signal_spectrogram`** &mdash; *Method*.



```julia
signal_spectrogram(signal; fs, norm=true, demean=true)
```

Calculates spectrogram of `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling frequency
  * `norm::Bool`: normalize powers to dB
  * `demean::Bool`: demean signal prior to analysis

**Returns**

  * `spec.power::Matrix{Float64}`
  * `spec.freq::Vector{Float64}`
  * `spec.time::Vector{Float64}`

<a id='NeuroJ.signal_spectrogram-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_spectrogram-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_spectrogram`** &mdash; *Method*.



```julia
signal_spectrogram(signal; fs, norm=true, demean=true)
```

Calculates spectrogram of `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `fs::Int64`: sampling frequency
  * `norm::Bool`: normalize powers to dB
  * `demean::Bool`: demean signal prior to analysis

**Returns**

  * `spec.power::Array{Float64, 4}`
  * `spec.freq::Matrix{Float64}`
  * `spec.time::Matrix{Float64}`

<a id='NeuroJ.signal_band-Tuple{Union{Float64, Int64}, Symbol}' href='#NeuroJ.signal_band-Tuple{Union{Float64, Int64}, Symbol}'>#</a>
**`NeuroJ.signal_band`** &mdash; *Method*.



```julia
signal_band(fs, band)
```

Returns EEG band frequency limits.

**Arguments**

  * `fs::Int64`
  * `band::Symbol`: EEG band name (:delta, :theta, :alpha, :beta, :beta*high, :gamma, :gamma*1, :gamma*2, :gamma*lower, :gamma_higher)

**Returns**

  * `band_frequency::Tuple{Float64, Float64}`: lower and upper bounds

<a id='NeuroJ.signal_detect_epoch_flat-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_detect_epoch_flat-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_detect_epoch_flat`** &mdash; *Method*.



```julia
signal_detect_epoch_flat(signal::Array{Float64, 3}, threshold=0.1)
```

Detect bad `signal` epochs based on: flat channel(s)

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.signal_detect_epoch_rmse-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_detect_epoch_rmse-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_detect_epoch_rmse`** &mdash; *Method*.



```julia
signal_detect_epoch_rmse(signal::Array{Float64, 3})
```

Detect bad `signal` epochs based on: RMSE vs average channel > 95%CI.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.signal_detect_epoch_rmsd-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_detect_epoch_rmsd-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_detect_epoch_rmsd`** &mdash; *Method*.



```julia
signal_detect_epoch_rmsd(signal::Array{Float64, 3})
```

Detect bad `signal` epochs based on: RMSD vs average channel > 95%CI.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.signal_detect_epoch_euclid-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_detect_epoch_euclid-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_detect_epoch_euclid`** &mdash; *Method*.



```julia
signal_detect_epoch_euclid(signal::Array{Float64, 3})
```

Detect bad `signal` epochs based on: Euclidean distance vs median channel > 95% CI.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.signal_detect_epoch_p2p-Tuple{Array{Float64, 3}}' href='#NeuroJ.signal_detect_epoch_p2p-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.signal_detect_epoch_p2p`** &mdash; *Method*.



```julia
signal_detect_epoch_p2p(signal::Array{Float64, 3})
```

Detect bad `signal` epochs based on: p2p amplitude > upper 95% CI p2p amplitude.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch


<a id='Misc'></a>

<a id='Misc-1'></a>

## Misc

<a id='NeuroJ.zero_pad-Tuple{Union{Matrix{ComplexF64}, Matrix{Float64}, Matrix{Int64}}}' href='#NeuroJ.zero_pad-Tuple{Union{Matrix{ComplexF64}, Matrix{Float64}, Matrix{Int64}}}'>#</a>
**`NeuroJ.zero_pad`** &mdash; *Method*.



```julia
zero_pad(m)
```

Pads the matrix `m` with zeros to make it square.

**Arguments**

  * `m::Union{Matrix{Int64}, Matrix{Float64}, Matrix{ComplexF64}}`

**Returns**

  * `m::Union{Matrix{Int64}, Matrix{Float64}, Matrix{ComplexF64}}`

<a id='NeuroJ.vsearch-Tuple{Union{Float64, Int64}, Union{Vector{Float64}, Vector{Int64}}}' href='#NeuroJ.vsearch-Tuple{Union{Float64, Int64}, Union{Vector{Float64}, Vector{Int64}}}'>#</a>
**`NeuroJ.vsearch`** &mdash; *Method*.



```julia
vsearch(y, x; return_distance=false)
```

Return the positions of the `y` value in the vector `x` and the difference between `y` and `x[vsearch(x, y)].

**Arguments**

  * `y::Union{Int64, Float64}`
  * `x::Union{Vector{Int64}, Vector{Float64}}
  * `return_distance::Bool`

**Returns**

  * `y_idx::Int64`

-`y_dist::Union{Int64, Float64}`

<a id='NeuroJ.vsearch-Tuple{Union{Vector{Float64}, Vector{Int64}}, Union{Vector{Float64}, Vector{Int64}}}' href='#NeuroJ.vsearch-Tuple{Union{Vector{Float64}, Vector{Int64}}, Union{Vector{Float64}, Vector{Int64}}}'>#</a>
**`NeuroJ.vsearch`** &mdash; *Method*.



```julia
vsearch(y, x; return_distance=false)
```

Return the positions of the `y` vector in the vector `x`.

**Arguments**

  * `x::Union{Vector{Int64}, Vector{Float64}}`
  * `y::Union{Vector{Int64}, Vector{Float64}}
  * `return_distance::Bool`

**Returns**

  * `y_idx::Int64`
  * `y_dist::Union{Int64, Float64}`

<a id='NeuroJ.jaccard_similarity-Tuple{Union{Vector{Float64}, Vector{Int64}}, Union{Vector{Float64}, Vector{Int64}}}' href='#NeuroJ.jaccard_similarity-Tuple{Union{Vector{Float64}, Vector{Int64}}, Union{Vector{Float64}, Vector{Int64}}}'>#</a>
**`NeuroJ.jaccard_similarity`** &mdash; *Method*.



```julia
jaccard_similarity(x, y)
```

Calculate Jaccard similarity between two vectors `x` and `y`.

**Arguments**

  * `n::Int64`

**Returns**

  * `j::Float64`

<a id='NeuroJ.fft0-Tuple{AbstractArray, Int64}' href='#NeuroJ.fft0-Tuple{AbstractArray, Int64}'>#</a>
**`NeuroJ.fft0`** &mdash; *Method*.



```julia
fft0(x, n)
```

Calculate FFT for the vector `x` padded with `n` or `n - length(x)` zeros at the end.

**Arguments**

  * `x::AbstractArray`
  * `n::Int64`

**Returns**

  * `fft0::Vector{ComplexF64}`

<a id='NeuroJ.ifft0-Tuple{AbstractArray, Int64}' href='#NeuroJ.ifft0-Tuple{AbstractArray, Int64}'>#</a>
**`NeuroJ.ifft0`** &mdash; *Method*.



```julia
ifft0(x, n)
```

Calculate IFFT for the vector `x` padded with `n` or `n - length(x)` zeros at the end.

**Arguments**

  * `x::AbstractArray`
  * `n::Int64`

**Returns**

  * `ifft0::Vector{ComplexF64}`

<a id='NeuroJ.nextpow2-Tuple{Int64}' href='#NeuroJ.nextpow2-Tuple{Int64}'>#</a>
**`NeuroJ.nextpow2`** &mdash; *Method*.



```julia
nextpow2(x)
```

Return the next power of 2 for given number `x`.

**Argument**

  * `x::Int64`

**Returns**

  * `nextpow::Int64`

<a id='NeuroJ.vsplit' href='#NeuroJ.vsplit'>#</a>
**`NeuroJ.vsplit`** &mdash; *Function*.



```julia
vsplit(x, n)
```

Splits the vector `x` into `n`-long pieces.

**Argument**

  * `x::Union{Vector{Int64}, Vector{Float64}}`
  * `n::Int64`

**Returns**

  * `x::Vector{Union{Vector{Int64}, Vector{Float64}}}`

<a id='NeuroJ.freqs-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}}' href='#NeuroJ.freqs-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}}'>#</a>
**`NeuroJ.freqs`** &mdash; *Method*.



```julia
freqs(t)
```

Return vector of frequencies and Nyquist frequency for given time vector `t`.

**Arguments**

  * `t::Union{Vector{Int64}, Vector{Float64}, AbstractRange}`

**Returns**

  * `hz::Vector{Float64}
  * `nyquist_freq::Float64`

<a id='NeuroJ.freqs-Tuple{Vector{Float64}, Union{Float64, Int64}}' href='#NeuroJ.freqs-Tuple{Vector{Float64}, Union{Float64, Int64}}'>#</a>
**`NeuroJ.freqs`** &mdash; *Method*.



```julia
freqs(signal, fs)
```

Return vector of frequencies and Nyquist frequency for given `signal` and `fs`.

**Arguments**

  * `signal::Vector{Float64}`
  * `fs::Union{Int64, Float64}`

**Returns**

  * `hz::Vector{Float64}
  * `nyquist_freq::Float64`

<a id='NeuroJ.matrix_sortperm-Tuple{Matrix}' href='#NeuroJ.matrix_sortperm-Tuple{Matrix}'>#</a>
**`NeuroJ.matrix_sortperm`** &mdash; *Method*.



```julia
matrix_sortperm(m; rev=false, dims=1)
```

Generates sorting index for matrix `m` by columns (`dims` = 1) or by rows (`dims` = 2).

**Arguments**

  * `m::Matrix`
  * `rev::Bool`
  * `dims::Int64`

**Returns**

  * `m_idx::Matrix{Int64}`

<a id='NeuroJ.matrix_sort-Tuple{Matrix, Vector{Int64}}' href='#NeuroJ.matrix_sort-Tuple{Matrix, Vector{Int64}}'>#</a>
**`NeuroJ.matrix_sort`** &mdash; *Method*.



```julia
matrix_sort(m, m_idx; rev=false, dims=1)
```

Sorts matrix `m` using sorting index `m_idx` by columns (`dims` = 1) or by rows (`dims` = 2).

**Arguments**

  * `m::Matrix`
  * `m_idx::Vector{Int64}`
  * `rev::Bool`
  * `dims::Int64`

**Returns**

  * `m_sorted::Matrix`

<a id='NeuroJ.pad0-Tuple{Union{Vector{Float64}, Vector{Int64}}, Any}' href='#NeuroJ.pad0-Tuple{Union{Vector{Float64}, Vector{Int64}}, Any}'>#</a>
**`NeuroJ.pad0`** &mdash; *Method*.



```julia
pad0(x, n)
```

Pads the vector `x` with `n` zeros at the beginning and at the end.

To do: check if x is numeric vector

**Arguments**

  * `x::Union{Vector{Int64}, Vector{Float64}}`
  * `n::Int64`

**Returns**

  * `v_pad::Union{Vector{Int64}, Vector{Float64}}`

<a id='NeuroJ.generate_sinc' href='#NeuroJ.generate_sinc'>#</a>
**`NeuroJ.generate_sinc`** &mdash; *Function*.



```julia
generate_sinc(t, f, peak)
```

Generate sinc function.

**Arguments**

  * `t::AbstractRange=-2:0.01:2`: time
  * `f::Union{Int64, Float64}=10.0`: frequency
  * `peak::Union{Int64, Float64}=0`: sinc peak time

**Returns**

  * `sinc::Vector{Float64}

<a id='NeuroJ.generate_morlet-Tuple{Int64, Union{Float64, Int64}, Union{Float64, Int64}}' href='#NeuroJ.generate_morlet-Tuple{Int64, Union{Float64, Int64}, Union{Float64, Int64}}'>#</a>
**`NeuroJ.generate_morlet`** &mdash; *Method*.



```julia
generate_morlet(fs, wt, wf)
```

Generate Morlet wavelet.

**Arguments**

  * `fs::Int64`: sampling rate
  * `wt::Union{Int64, Float64}`: length = -wt:1/fs:wt
  * `wf::Union{Int64, Float64}`: frequency
  * `ncyc::Int64=5`: number of cycles
  * `complex::Bool=false`: generate complex Morlet

**Returns**

  * `morlet::Union{Vector{Float64}, Vector{ComplexF64}}`

<a id='NeuroJ.generate_gaussian-Tuple{Int64, Union{Float64, Int64}, Union{Float64, Int64}}' href='#NeuroJ.generate_gaussian-Tuple{Int64, Union{Float64, Int64}, Union{Float64, Int64}}'>#</a>
**`NeuroJ.generate_gaussian`** &mdash; *Method*.



```julia
generate_gaussian(fs, wt, wf)
```

Generate Gaussian wave.

**Arguments**

  * `fs::Int64`: sampling rate
  * `gt::Union{Int64, Float64}`: length = -wt:1/fs:wt
  * `gw::Union{Int64, Float64}`: width

**Returns**

  * `gaussian::Vector{Float64}`

<a id='NeuroJ.tuple_order' href='#NeuroJ.tuple_order'>#</a>
**`NeuroJ.tuple_order`** &mdash; *Function*.



```julia
tuple_order(t, rev)
```

Order tuple elements in ascending or descending (rev=true) order.

**Arguments**

  * `t::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}`
  * `rev::Bool=false`

**Returns**

  * `t::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}`

<a id='NeuroJ.rmse-Tuple{Vector{Float64}, Vector{Float64}}' href='#NeuroJ.rmse-Tuple{Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroJ.rmse`** &mdash; *Method*.



```julia
rmse(signal1, signal2)
```

Calculate RMSE between `signal1` and `signal2`.

**Arguments**

  * `signal1::Vector{Float64}`
  * `signal2::Vector{Float64}`

**Returns**

  * `r::Float64`


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


<a id='Index'></a>

<a id='Index-1'></a>

## Index

- [`NeuroJ.eeg_add_labels`](index.md#NeuroJ.eeg_add_labels-Tuple{NeuroJ.EEG, Vector{String}})
- [`NeuroJ.eeg_add_labels!`](index.md#NeuroJ.eeg_add_labels!-Tuple{NeuroJ.EEG, Vector{String}})
- [`NeuroJ.eeg_autocov`](index.md#NeuroJ.eeg_autocov-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_autocov!`](index.md#NeuroJ.eeg_autocov!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_average`](index.md#NeuroJ.eeg_average-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_average!`](index.md#NeuroJ.eeg_average!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_band`](index.md#NeuroJ.eeg_band-Tuple{Any})
- [`NeuroJ.eeg_band_power`](index.md#NeuroJ.eeg_band_power-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_channel_n`](index.md#NeuroJ.eeg_channel_n-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_coherence`](index.md#NeuroJ.eeg_coherence-Tuple{NeuroJ.EEG, NeuroJ.EEG})
- [`NeuroJ.eeg_coherence`](index.md#NeuroJ.eeg_coherence-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_cor`](index.md#NeuroJ.eeg_cor-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_cor!`](index.md#NeuroJ.eeg_cor!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_cov`](index.md#NeuroJ.eeg_cov-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_cov!`](index.md#NeuroJ.eeg_cov!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_crosscov`](index.md#NeuroJ.eeg_crosscov-Tuple{NeuroJ.EEG, NeuroJ.EEG})
- [`NeuroJ.eeg_crosscov`](index.md#NeuroJ.eeg_crosscov-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_crosscov!`](index.md#NeuroJ.eeg_crosscov!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_delete_bad_epochs`](index.md#NeuroJ.eeg_delete_bad_epochs-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_delete_bad_epochs!`](index.md#NeuroJ.eeg_delete_bad_epochs!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_delete_channel`](index.md#NeuroJ.eeg_delete_channel-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_delete_channel!`](index.md#NeuroJ.eeg_delete_channel!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_delete_component`](index.md#NeuroJ.eeg_delete_component-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_delete_component!`](index.md#NeuroJ.eeg_delete_component!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_delete_epoch`](index.md#NeuroJ.eeg_delete_epoch-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_delete_epoch!`](index.md#NeuroJ.eeg_delete_epoch!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_demean`](index.md#NeuroJ.eeg_demean-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_demean!`](index.md#NeuroJ.eeg_demean!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_derivative`](index.md#NeuroJ.eeg_derivative-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_derivative!`](index.md#NeuroJ.eeg_derivative!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_detect_bad_epochs`](index.md#NeuroJ.eeg_detect_bad_epochs-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_detrend`](index.md#NeuroJ.eeg_detrend-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_detrend!`](index.md#NeuroJ.eeg_detrend!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_difference`](index.md#NeuroJ.eeg_difference-Tuple{NeuroJ.EEG, NeuroJ.EEG})
- [`NeuroJ.eeg_downsample`](index.md#NeuroJ.eeg_downsample-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_downsample!`](index.md#NeuroJ.eeg_downsample!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_draw_head`](index.md#NeuroJ.eeg_draw_head-Tuple{Plots.Plot{Plots.GRBackend}, Vector{Float64}, Vector{Float64}})
- [`NeuroJ.eeg_edit_channel`](index.md#NeuroJ.eeg_edit_channel-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_edit_channel!`](index.md#NeuroJ.eeg_edit_channel!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_edit_header`](index.md#NeuroJ.eeg_edit_header-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_edit_header!`](index.md#NeuroJ.eeg_edit_header!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_entropy`](index.md#NeuroJ.eeg_entropy-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_entropy!`](index.md#NeuroJ.eeg_entropy!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_epoch_len`](index.md#NeuroJ.eeg_epoch_len-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_epoch_n`](index.md#NeuroJ.eeg_epoch_n-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_epochs`](index.md#NeuroJ.eeg_epochs-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_epochs!`](index.md#NeuroJ.eeg_epochs!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_epochs_stats`](index.md#NeuroJ.eeg_epochs_stats-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_epochs_stats!`](index.md#NeuroJ.eeg_epochs_stats!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_export_csv`](index.md#NeuroJ.eeg_export_csv-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_extract_channel`](index.md#NeuroJ.eeg_extract_channel-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_extract_component`](index.md#NeuroJ.eeg_extract_component-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_extract_epoch`](index.md#NeuroJ.eeg_extract_epoch-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_fconv`](index.md#NeuroJ.eeg_fconv-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_fconv!`](index.md#NeuroJ.eeg_fconv!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_filter`](index.md#NeuroJ.eeg_filter-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_filter!`](index.md#NeuroJ.eeg_filter!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_freqs`](index.md#NeuroJ.eeg_freqs-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_freqs!`](index.md#NeuroJ.eeg_freqs!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_get_channel`](index.md#NeuroJ.eeg_get_channel-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_history`](index.md#NeuroJ.eeg_history-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_ica`](index.md#NeuroJ.eeg_ica-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_ica!`](index.md#NeuroJ.eeg_ica!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_ica_reconstruct`](index.md#NeuroJ.eeg_ica_reconstruct-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_ica_reconstruct!`](index.md#NeuroJ.eeg_ica_reconstruct!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_import_ced`](index.md#NeuroJ.eeg_import_ced-Tuple{String})
- [`NeuroJ.eeg_import_edf`](index.md#NeuroJ.eeg_import_edf-Tuple{String})
- [`NeuroJ.eeg_import_elc`](index.md#NeuroJ.eeg_import_elc-Tuple{String})
- [`NeuroJ.eeg_import_locs`](index.md#NeuroJ.eeg_import_locs-Tuple{String})
- [`NeuroJ.eeg_info`](index.md#NeuroJ.eeg_info-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_keep_channel`](index.md#NeuroJ.eeg_keep_channel-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_keep_channel!`](index.md#NeuroJ.eeg_keep_channel!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_keep_eeg_channels`](index.md#NeuroJ.eeg_keep_eeg_channels-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_keep_eeg_channels!`](index.md#NeuroJ.eeg_keep_eeg_channels!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_keep_epoch`](index.md#NeuroJ.eeg_keep_epoch-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_keep_epoch!`](index.md#NeuroJ.eeg_keep_epoch!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_labels`](index.md#NeuroJ.eeg_labels-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_list_components`](index.md#NeuroJ.eeg_list_components-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_load`](index.md#NeuroJ.eeg_load-Tuple{String})
- [`NeuroJ.eeg_load_electrodes`](index.md#NeuroJ.eeg_load_electrodes-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_load_electrodes!`](index.md#NeuroJ.eeg_load_electrodes!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_mi`](index.md#NeuroJ.eeg_mi-Tuple{NeuroJ.EEG, NeuroJ.EEG})
- [`NeuroJ.eeg_mi`](index.md#NeuroJ.eeg_mi-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_mi!`](index.md#NeuroJ.eeg_mi!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_normalize_minmax`](index.md#NeuroJ.eeg_normalize_minmax-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_normalize_minmax!`](index.md#NeuroJ.eeg_normalize_minmax!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_normalize_zscore`](index.md#NeuroJ.eeg_normalize_zscore-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_normalize_zscore!`](index.md#NeuroJ.eeg_normalize_zscore!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_pca`](index.md#NeuroJ.eeg_pca-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_pca!`](index.md#NeuroJ.eeg_pca!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_pick`](index.md#NeuroJ.eeg_pick-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot`](index.md#NeuroJ.eeg_plot-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_avg`](index.md#NeuroJ.eeg_plot_avg-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_bands`](index.md#NeuroJ.eeg_plot_bands-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_butterfly`](index.md#NeuroJ.eeg_plot_butterfly-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_covmatrix`](index.md#NeuroJ.eeg_plot_covmatrix-Tuple{NeuroJ.EEG, Union{Array{Float64, 3}, Matrix{Float64}}, Union{Vector{Float64}, Vector{Int64}}})
- [`NeuroJ.eeg_plot_electrodes`](index.md#NeuroJ.eeg_plot_electrodes-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_filter_response`](index.md#NeuroJ.eeg_plot_filter_response-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_histogram`](index.md#NeuroJ.eeg_plot_histogram-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_ica`](index.md#NeuroJ.eeg_plot_ica-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_matrix`](index.md#NeuroJ.eeg_plot_matrix-Tuple{NeuroJ.EEG, Union{Array{Float64, 3}, Matrix{Float64}}})
- [`NeuroJ.eeg_plot_psd`](index.md#NeuroJ.eeg_plot_psd-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_save`](index.md#NeuroJ.eeg_plot_save-Tuple{Plots.Plot{Plots.GRBackend}})
- [`NeuroJ.eeg_plot_spectrogram`](index.md#NeuroJ.eeg_plot_spectrogram-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_plot_topo`](index.md#NeuroJ.eeg_plot_topo-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_psd`](index.md#NeuroJ.eeg_psd-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_psd!`](index.md#NeuroJ.eeg_psd!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_reference_car`](index.md#NeuroJ.eeg_reference_car-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_reference_car!`](index.md#NeuroJ.eeg_reference_car!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_reference_channel`](index.md#NeuroJ.eeg_reference_channel-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_reference_channel!`](index.md#NeuroJ.eeg_reference_channel!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_rename_channel`](index.md#NeuroJ.eeg_rename_channel-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_rename_channel!`](index.md#NeuroJ.eeg_rename_channel!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_resample`](index.md#NeuroJ.eeg_resample-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_resample!`](index.md#NeuroJ.eeg_resample!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_reset_components`](index.md#NeuroJ.eeg_reset_components-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_reset_components!`](index.md#NeuroJ.eeg_reset_components!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_s2t`](index.md#NeuroJ.eeg_s2t-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_save`](index.md#NeuroJ.eeg_save-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_show_header`](index.md#NeuroJ.eeg_show_header-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_signal_len`](index.md#NeuroJ.eeg_signal_len-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_spectrogram`](index.md#NeuroJ.eeg_spectrogram-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_spectrogram!`](index.md#NeuroJ.eeg_spectrogram!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_spectrum`](index.md#NeuroJ.eeg_spectrum-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_spectrum!`](index.md#NeuroJ.eeg_spectrum!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_sr`](index.md#NeuroJ.eeg_sr-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_stationarity`](index.md#NeuroJ.eeg_stationarity-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_stationarity!`](index.md#NeuroJ.eeg_stationarity!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_t2s`](index.md#NeuroJ.eeg_t2s-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_taper`](index.md#NeuroJ.eeg_taper-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_taper!`](index.md#NeuroJ.eeg_taper!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_tconv`](index.md#NeuroJ.eeg_tconv-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_tconv!`](index.md#NeuroJ.eeg_tconv!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_total_power`](index.md#NeuroJ.eeg_total_power-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_total_power!`](index.md#NeuroJ.eeg_total_power!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_trim`](index.md#NeuroJ.eeg_trim-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_trim!`](index.md#NeuroJ.eeg_trim!-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_upsample`](index.md#NeuroJ.eeg_upsample-Tuple{NeuroJ.EEG})
- [`NeuroJ.eeg_upsample!`](index.md#NeuroJ.eeg_upsample!-Tuple{NeuroJ.EEG})
- [`NeuroJ.fft0`](index.md#NeuroJ.fft0-Tuple{AbstractArray, Int64})
- [`NeuroJ.freqs`](index.md#NeuroJ.freqs-Tuple{Vector{Float64}, Union{Float64, Int64}})
- [`NeuroJ.freqs`](index.md#NeuroJ.freqs-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}})
- [`NeuroJ.generate_gaussian`](index.md#NeuroJ.generate_gaussian-Tuple{Int64, Union{Float64, Int64}, Union{Float64, Int64}})
- [`NeuroJ.generate_morlet`](index.md#NeuroJ.generate_morlet-Tuple{Int64, Union{Float64, Int64}, Union{Float64, Int64}})
- [`NeuroJ.generate_sinc`](index.md#NeuroJ.generate_sinc)
- [`NeuroJ.ifft0`](index.md#NeuroJ.ifft0-Tuple{AbstractArray, Int64})
- [`NeuroJ.jaccard_similarity`](index.md#NeuroJ.jaccard_similarity-Tuple{Union{Vector{Float64}, Vector{Int64}}, Union{Vector{Float64}, Vector{Int64}}})
- [`NeuroJ.matrix_sort`](index.md#NeuroJ.matrix_sort-Tuple{Matrix, Vector{Int64}})
- [`NeuroJ.matrix_sortperm`](index.md#NeuroJ.matrix_sortperm-Tuple{Matrix})
- [`NeuroJ.neuroj_reload_plugins`](index.md#NeuroJ.neuroj_reload_plugins-Tuple{})
- [`NeuroJ.neuroj_version`](index.md#NeuroJ.neuroj_version-Tuple{})
- [`NeuroJ.nextpow2`](index.md#NeuroJ.nextpow2-Tuple{Int64})
- [`NeuroJ.pad0`](index.md#NeuroJ.pad0-Tuple{Union{Vector{Float64}, Vector{Int64}}, Any})
- [`NeuroJ.rmse`](index.md#NeuroJ.rmse-Tuple{Vector{Float64}, Vector{Float64}})
- [`NeuroJ.signal_add_noise`](index.md#NeuroJ.signal_add_noise-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_add_noise`](index.md#NeuroJ.signal_add_noise-Tuple{AbstractArray})
- [`NeuroJ.signal_autocov`](index.md#NeuroJ.signal_autocov-Tuple{AbstractArray})
- [`NeuroJ.signal_autocov`](index.md#NeuroJ.signal_autocov-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_average`](index.md#NeuroJ.signal_average-Tuple{AbstractArray, AbstractArray})
- [`NeuroJ.signal_average`](index.md#NeuroJ.signal_average-Tuple{Array{Float64, 3}, Array{Float64, 3}})
- [`NeuroJ.signal_average`](index.md#NeuroJ.signal_average-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_band`](index.md#NeuroJ.signal_band-Tuple{Union{Float64, Int64}, Symbol})
- [`NeuroJ.signal_band_power`](index.md#NeuroJ.signal_band_power-Tuple{AbstractArray})
- [`NeuroJ.signal_band_power`](index.md#NeuroJ.signal_band_power-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_ci95`](index.md#NeuroJ.signal_ci95-Tuple{Vector{Float64}})
- [`NeuroJ.signal_ci95`](index.md#NeuroJ.signal_ci95-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_ci95`](index.md#NeuroJ.signal_ci95-Tuple{AbstractArray})
- [`NeuroJ.signal_coherence`](index.md#NeuroJ.signal_coherence-Tuple{Matrix{Float64}, Matrix{Float64}})
- [`NeuroJ.signal_coherence`](index.md#NeuroJ.signal_coherence-Tuple{AbstractArray, AbstractArray})
- [`NeuroJ.signal_coherence`](index.md#NeuroJ.signal_coherence-Tuple{Array{Float64, 3}, Array{Float64, 3}})
- [`NeuroJ.signal_cor`](index.md#NeuroJ.signal_cor-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_cov`](index.md#NeuroJ.signal_cov-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_cov`](index.md#NeuroJ.signal_cov-Tuple{AbstractArray, AbstractArray})
- [`NeuroJ.signal_cov`](index.md#NeuroJ.signal_cov-Tuple{AbstractArray})
- [`NeuroJ.signal_crosscov`](index.md#NeuroJ.signal_crosscov-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_crosscov`](index.md#NeuroJ.signal_crosscov-Tuple{AbstractArray, AbstractArray})
- [`NeuroJ.signal_crosscov`](index.md#NeuroJ.signal_crosscov-Tuple{Array{Float64, 3}, Array{Float64, 3}})
- [`NeuroJ.signal_delete_channel`](index.md#NeuroJ.signal_delete_channel-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_delete_channel`](index.md#NeuroJ.signal_delete_channel-Tuple{Matrix{Float64}})
- [`NeuroJ.signal_demean`](index.md#NeuroJ.signal_demean-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_demean`](index.md#NeuroJ.signal_demean-Tuple{AbstractArray})
- [`NeuroJ.signal_derivative`](index.md#NeuroJ.signal_derivative-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_derivative`](index.md#NeuroJ.signal_derivative-Tuple{AbstractArray})
- [`NeuroJ.signal_detect_epoch_euclid`](index.md#NeuroJ.signal_detect_epoch_euclid-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_detect_epoch_flat`](index.md#NeuroJ.signal_detect_epoch_flat-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_detect_epoch_p2p`](index.md#NeuroJ.signal_detect_epoch_p2p-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_detect_epoch_rmsd`](index.md#NeuroJ.signal_detect_epoch_rmsd-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_detect_epoch_rmse`](index.md#NeuroJ.signal_detect_epoch_rmse-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_detrend`](index.md#NeuroJ.signal_detrend-Tuple{AbstractArray})
- [`NeuroJ.signal_detrend`](index.md#NeuroJ.signal_detrend-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_difference`](index.md#NeuroJ.signal_difference-Tuple{AbstractArray, AbstractArray})
- [`NeuroJ.signal_difference`](index.md#NeuroJ.signal_difference-Tuple{Array{Float64, 3}, Array{Float64, 3}})
- [`NeuroJ.signal_downsample`](index.md#NeuroJ.signal_downsample-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_downsample`](index.md#NeuroJ.signal_downsample-Tuple{AbstractArray})
- [`NeuroJ.signal_entropy`](index.md#NeuroJ.signal_entropy-Tuple{AbstractArray})
- [`NeuroJ.signal_entropy`](index.md#NeuroJ.signal_entropy-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_epochs`](index.md#NeuroJ.signal_epochs-Tuple{Vector{Float64}})
- [`NeuroJ.signal_epochs`](index.md#NeuroJ.signal_epochs-Tuple{Matrix{Float64}})
- [`NeuroJ.signal_epochs_stats`](index.md#NeuroJ.signal_epochs_stats-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_fconv`](index.md#NeuroJ.signal_fconv-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_fconv`](index.md#NeuroJ.signal_fconv-Tuple{AbstractArray})
- [`NeuroJ.signal_filter`](index.md#NeuroJ.signal_filter-Tuple{AbstractArray})
- [`NeuroJ.signal_filter`](index.md#NeuroJ.signal_filter-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_ica`](index.md#NeuroJ.signal_ica-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_ica_reconstruct`](index.md#NeuroJ.signal_ica_reconstruct-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_make_spectrum`](index.md#NeuroJ.signal_make_spectrum-Tuple{AbstractArray})
- [`NeuroJ.signal_make_spectrum`](index.md#NeuroJ.signal_make_spectrum-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_mean`](index.md#NeuroJ.signal_mean-Tuple{Vector{Float64}, Vector{Float64}})
- [`NeuroJ.signal_mean`](index.md#NeuroJ.signal_mean-Tuple{Array{Float64, 3}, Array{Float64, 3}})
- [`NeuroJ.signal_mi`](index.md#NeuroJ.signal_mi-Tuple{AbstractArray, AbstractArray})
- [`NeuroJ.signal_mi`](index.md#NeuroJ.signal_mi-Tuple{Array{Float64, 3}, Array{Float64, 3}})
- [`NeuroJ.signal_mi`](index.md#NeuroJ.signal_mi-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_normalize_minmax`](index.md#NeuroJ.signal_normalize_minmax-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_normalize_minmax`](index.md#NeuroJ.signal_normalize_minmax-Tuple{AbstractArray})
- [`NeuroJ.signal_normalize_zscore`](index.md#NeuroJ.signal_normalize_zscore-Tuple{AbstractArray})
- [`NeuroJ.signal_normalize_zscore`](index.md#NeuroJ.signal_normalize_zscore-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_pca`](index.md#NeuroJ.signal_pca-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_plot`](index.md#NeuroJ.signal_plot-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Vector{Float64}})
- [`NeuroJ.signal_plot`](index.md#NeuroJ.signal_plot-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, AbstractArray})
- [`NeuroJ.signal_plot_avg`](index.md#NeuroJ.signal_plot_avg-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Matrix{Float64}})
- [`NeuroJ.signal_plot_bands`](index.md#NeuroJ.signal_plot_bands-Tuple{Vector{Float64}})
- [`NeuroJ.signal_plot_butterfly`](index.md#NeuroJ.signal_plot_butterfly-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Matrix{Float64}})
- [`NeuroJ.signal_plot_histogram`](index.md#NeuroJ.signal_plot_histogram-Tuple{VecOrMat{Float64}})
- [`NeuroJ.signal_plot_histogram`](index.md#NeuroJ.signal_plot_histogram-Tuple{Vector{Float64}})
- [`NeuroJ.signal_plot_ica`](index.md#NeuroJ.signal_plot_ica-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Matrix{Float64}})
- [`NeuroJ.signal_plot_ica`](index.md#NeuroJ.signal_plot_ica-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}, Vector{Float64}})
- [`NeuroJ.signal_plot_psd`](index.md#NeuroJ.signal_plot_psd-Tuple{Vector{Float64}, Vector{Float64}})
- [`NeuroJ.signal_plot_psd`](index.md#NeuroJ.signal_plot_psd-Tuple{Vector{Float64}})
- [`NeuroJ.signal_plot_psd`](index.md#NeuroJ.signal_plot_psd-Tuple{Matrix{Float64}})
- [`NeuroJ.signal_plot_spectrogram`](index.md#NeuroJ.signal_plot_spectrogram-Tuple{Vector{Float64}})
- [`NeuroJ.signal_psd`](index.md#NeuroJ.signal_psd-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_psd`](index.md#NeuroJ.signal_psd-Tuple{AbstractArray})
- [`NeuroJ.signal_psd`](index.md#NeuroJ.signal_psd-Tuple{Matrix{Float64}})
- [`NeuroJ.signal_reference_car`](index.md#NeuroJ.signal_reference_car-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_reference_car`](index.md#NeuroJ.signal_reference_car-Tuple{Matrix{Float64}})
- [`NeuroJ.signal_reference_channel`](index.md#NeuroJ.signal_reference_channel-Tuple{Matrix{Float64}})
- [`NeuroJ.signal_reference_channel`](index.md#NeuroJ.signal_reference_channel-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_spectrogram`](index.md#NeuroJ.signal_spectrogram-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_spectrogram`](index.md#NeuroJ.signal_spectrogram-Tuple{AbstractArray})
- [`NeuroJ.signal_spectrum`](index.md#NeuroJ.signal_spectrum-Tuple{AbstractArray})
- [`NeuroJ.signal_spectrum`](index.md#NeuroJ.signal_spectrum-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_stationarity`](index.md#NeuroJ.signal_stationarity-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_stationarity_hilbert`](index.md#NeuroJ.signal_stationarity_hilbert-Tuple{AbstractArray})
- [`NeuroJ.signal_stationarity_mean`](index.md#NeuroJ.signal_stationarity_mean-Tuple{AbstractArray})
- [`NeuroJ.signal_stationarity_var`](index.md#NeuroJ.signal_stationarity_var-Tuple{AbstractArray})
- [`NeuroJ.signal_taper`](index.md#NeuroJ.signal_taper-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_taper`](index.md#NeuroJ.signal_taper-Tuple{AbstractArray})
- [`NeuroJ.signal_tconv`](index.md#NeuroJ.signal_tconv-Tuple{AbstractArray})
- [`NeuroJ.signal_tconv`](index.md#NeuroJ.signal_tconv-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_total_power`](index.md#NeuroJ.signal_total_power-Tuple{AbstractArray})
- [`NeuroJ.signal_total_power`](index.md#NeuroJ.signal_total_power-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_trim`](index.md#NeuroJ.signal_trim-Tuple{Array{Float64, 3}})
- [`NeuroJ.signal_trim`](index.md#NeuroJ.signal_trim-Tuple{AbstractArray})
- [`NeuroJ.signal_upsample`](index.md#NeuroJ.signal_upsample-Tuple{AbstractArray})
- [`NeuroJ.signal_upsample`](index.md#NeuroJ.signal_upsample-Tuple{Array{Float64, 3}})
- [`NeuroJ.tes_dose`](index.md#NeuroJ.tes_dose-Tuple{Union{Float64, Int64}, Union{Float64, Int64}, Int64})
- [`NeuroJ.tuple_order`](index.md#NeuroJ.tuple_order)
- [`NeuroJ.vsearch`](index.md#NeuroJ.vsearch-Tuple{Union{Float64, Int64}, Union{Vector{Float64}, Vector{Int64}}})
- [`NeuroJ.vsearch`](index.md#NeuroJ.vsearch-Tuple{Union{Vector{Float64}, Vector{Int64}}, Union{Vector{Float64}, Vector{Int64}}})
- [`NeuroJ.vsplit`](index.md#NeuroJ.vsplit)
- [`NeuroJ.zero_pad`](index.md#NeuroJ.zero_pad-Tuple{Union{Matrix{ComplexF64}, Matrix{Float64}, Matrix{Int64}}})

