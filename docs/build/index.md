
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
eeg_edit(file_name; read_annotations, clean_labels)
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

Removes `channel` from the `eeg`.

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

Removes `channel` from the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to be removed

<a id='NeuroJ.eeg_keep_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_keep_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_keep_channel`** &mdash; *Method*.



```julia
eeg_keep_channel(eeg; channel)
```

Keeps `channels` in the `eeg`.

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

Keeps `channels` in the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to keep

<a id='NeuroJ.eeg_get_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_get_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_get_channel`** &mdash; *Method*.



```julia
eeg_get_channel(eeg; channel)
```

Returns the `channel` index / name.

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

Extracts `channel` number or name.

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

Shows processing history.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_labels-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_labels-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_labels`** &mdash; *Method*.



```julia
eeg_labels(eeg)
```

Returns labels.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_sr-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_sr-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_sr`** &mdash; *Method*.



```julia
eeg_sr(eeg)
```

Returns sampling rate.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_channel_n-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_channel_n-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_channel_n`** &mdash; *Method*.



```julia
eeg_channel_n(eeg; type=:eeg)
```

Returns number of `eeg` channels of `type`.

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

Returns number of `eeg` epochs.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `epoch_n::Int64`

<a id='NeuroJ.eeg_signal_len-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_signal_len-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_signal_len`** &mdash; *Method*.



```julia
eeg_signal_len(eeg)
```

Returns length of `eeg` signal.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `signal_len::Int64`

<a id='NeuroJ.eeg_epoch_len-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epoch_len-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epoch_len`** &mdash; *Method*.



```julia
eeg_epoch_len(eeg)
```

Returns length of `eeg` signal.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `epoch_len::Int64`

<a id='NeuroJ.eeg_info-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_info-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_info`** &mdash; *Method*.



```julia
eeg_info(eeg)
```

Shows info.

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
  * `average::Bool`: average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

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
  * `average::Bool`: average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

<a id='NeuroJ.eeg_extract_epoch-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_extract_epoch-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_extract_epoch`** &mdash; *Method*.



```julia
eeg_extract_epoch(eeg; epoch)
```

Extracts the `epoch` epoch.

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

Removes `len` samples from the beginning + `offset` (`from` = :start, default) or end (`from` = :end) of the `eeg`.

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

Removes `len` samples from the beginning + `offset` (`from` = :start, default) or end (`from` = :end) of the `eeg`.

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

Changes value of `eeg` `field` to `value`.

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

Changes value of `eeg` `field` to `value`.

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

Shows keys and values of `eeg` header.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_delete_epoch-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_epoch-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_epoch`** &mdash; *Method*.



```julia
eeg_delete_epoch(eeg; epoch)
```

Removes `epoch` from the `eeg`.

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

Removes `epoch` from the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch index to be removed, vector of numbers or range

<a id='NeuroJ.eeg_keep_epoch-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_keep_epoch-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_keep_epoch`** &mdash; *Method*.



```julia
eeg_keep_epoch(eeg; epoch)
```

Keeps `epoch` in the `eeg`.

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

Keeps `epoch` in the `eeg`.

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

Deletes bad `eeg` epochs.

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

Deletes bad `eeg` epochs.

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

Adds `labels` to `eeg` channels.

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

Adds `labels` to `eeg` channels.

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

Edits `eeg` `channel` properties.

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

Keeps only EEG channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_keep_eeg_channels!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_keep_eeg_channels!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_keep_eeg_channels!`** &mdash; *Method*.



```julia
eeg_keep_eeg_channels!(eeg::NeuroJ.EEG)
```

Keeps only EEG channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_list_components-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_list_components-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_list_components`** &mdash; *Method*.



```julia
eeg_list_components(eeg)
```

Lists `eeg` components.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `components::Vector{Symbol}`

<a id='NeuroJ.eeg_extract_component-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_extract_component-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_extract_component`** &mdash; *Method*.



```julia
eeg_extract_component(eeg, c)
```

Extracts `component` of `eeg`.

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

Deletes `component` of `eeg`.

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

Deletes `component` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `component::Symbol`

<a id='NeuroJ.eeg_reset_components-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reset_components-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reset_components`** &mdash; *Method*.



```julia
eeg_reset_components(eeg)
```

Resets `eeg` components.

**Arguments**

  * `eeg:EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_reset_components!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reset_components!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reset_components!`** &mdash; *Method*.



```julia
eeg_reset_components!(eeg)
```

Resets `eeg` components.

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

