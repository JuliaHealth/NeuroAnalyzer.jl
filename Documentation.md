
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

<a id='NeuroJ.eeg_copy-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_copy-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_copy`** &mdash; *Method*.



```julia
eeg_copy(eeg::NeuroJ.EEG)
```

Make copy of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_add_component-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_add_component-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_add_component`** &mdash; *Method*.



```julia
eeg_add_component(eeg; c, v)
```

Add component name `c` of value `v` to `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `c::Symbol`: component name
  * `v::Any`: component value

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_add_component!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_add_component!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_add_component!`** &mdash; *Method*.



```julia
eeg_add_component!(eeg; c, v)
```

Add component name `c` of value `v` to `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `c::Symbol`: component name
  * `v::Any`: component value

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

Extract component `c` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `c::Symbol`: component name

**Returns**

  * `component::Any`

<a id='NeuroJ.eeg_delete_component-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_component-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_component`** &mdash; *Method*.



```julia
eeg_delete_component(eeg; c)
```

Delete component `c` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `c::Symbol`: component name

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_delete_component!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_delete_component!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_delete_component!`** &mdash; *Method*.



```julia
eeg_delete_component!(eeg; c)
```

Delete component `c` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `c::Symbol`: component name

<a id='NeuroJ.eeg_reset_components-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reset_components-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reset_components`** &mdash; *Method*.



```julia
eeg_reset_components(eeg)
```

Remove all `eeg` components.

**Arguments**

  * `eeg:EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_reset_components!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reset_components!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reset_components!`** &mdash; *Method*.



```julia
eeg_reset_components!(eeg)
```

Remove all `eeg` components.

**Arguments**

  * `eeg:EEG`

<a id='NeuroJ.eeg_component_idx-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_component_idx-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_component_idx`** &mdash; *Method*.



```julia
eeg_component_idx(eeg, c)
```

Return index of `eeg` component.

**Arguments**

  * `eeg:EEG`
  * `c::Symbol`: component name

**Return**

  * `c_idx::Int64`

<a id='NeuroJ.eeg_component_type-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_component_type-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_component_type`** &mdash; *Method*.



```julia
eeg_component_type(eeg, c)
```

Return type of `eeg` components.

**Arguments**

  * `eeg:EEG`
  * `c::Symbol`: component name

**Return**

  * `c_type::DataType`

<a id='NeuroJ.eeg_rename_component-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_rename_component-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_rename_component`** &mdash; *Method*.



```julia
eeg_rename_component(eeg, c_old, c_new)
```

Return type of `eeg` components.

**Arguments**

  * `eeg:EEG`
  * `c_old::Symbol`: old component name
  * `c_new::Symbol`: new component name

**Return**

  * `eeg_new:EEG`

<a id='NeuroJ.eeg_rename_component!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_rename_component!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_rename_component!`** &mdash; *Method*.



```julia
eeg_rename_component(eeg, c_old, c_new)
```

Return type of `eeg` components.

**Arguments**

  * `eeg:EEG`
  * `c_old::Symbol`: old component name
  * `c_new::Symbol`: new component name

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

Return the `channel` index / name.

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

Return `eeg` labels.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `labels::Vector{String}`

<a id='NeuroJ.eeg_sr-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_sr-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_sr`** &mdash; *Method*.



```julia
eeg_sr(eeg)
```

Return `eeg` sampling rate.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `sr::Int64`

<a id='NeuroJ.eeg_channel_n-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_channel_n-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_channel_n`** &mdash; *Method*.



```julia
eeg_channel_n(eeg; type=:eeg)
```

Return number of `eeg` channels of `type`.

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

Return number of `eeg` epochs.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `epoch_n::Int64`

<a id='NeuroJ.eeg_signal_len-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_signal_len-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_signal_len`** &mdash; *Method*.



```julia
eeg_signal_len(eeg)
```

Return length of `eeg` signal.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `signal_len::Int64`

<a id='NeuroJ.eeg_epoch_len-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epoch_len-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epoch_len`** &mdash; *Method*.



```julia
eeg_epoch_len(eeg)
```

Return length of `eeg` signal.

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
  * `average::Bool`: average all epochs, return one averaged epoch

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

<a id='NeuroJ.eeg_comment-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_comment-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_comment`** &mdash; *Method*.



```julia
eeg_comment(eeg)
```

Return `eeg` comment.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_epochs_time-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epochs_time-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epochs_time`** &mdash; *Method*.



```julia
eeg_epochs_time(eeg; ts)
```

Edit `eeg` epochs time start.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `ts::Union{Int64, Float64}`: time start in seconds

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_epochs_time!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epochs_time!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epochs_time!`** &mdash; *Method*.



```julia
eeg_epochs_time!(eeg; ts)
```

Edit `eeg` epochs time start.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `ts::Union{Int64, Float64}`: time start in seconds

**Returns**

  * `eeg::NeuroJ.EEG`


<a id='EEG-process'></a>

<a id='EEG-process-1'></a>

## EEG process

<a id='NeuroJ.eeg_reference_channel-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_channel-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_channel`** &mdash; *Method*.



```julia
eeg_reference_channel(eeg; channel)
```

Reference the `eeg` to specific `channel`.

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

Perform piecewise detrending of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Symbol`, optional

      * `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `:linear`: linear trend is subtracted from `signal`
      * `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
      * `:poly`: polynomial of `order` is subtracted
      * `:loess`: fit and subtract loess approximation
  * `offset::Union{Int64, Float64}=0`: constant for :constant detrending
  * `order::Int64=1`: polynomial fitting order
  * `span::Float64=0.5`: smoothing of loess

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
  * `type::Symbol`, optional

      * `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `:linear`: linear trend is subtracted from `signal`
      * `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
      * `:poly`: polynomial of `order` order is subtracted
      * `:loess`: fit and subtract loess approximation
  * `offset::Union{Int64, Float64}=0`: constant for :constant detrending
  * `order::Int64=1`: polynomial fitting order
  * `span::Float64`: smoothing of loess

<a id='NeuroJ.eeg_taper-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_taper-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_taper`** &mdash; *Method*.



```julia
eeg_taper(eeg; taper)
```

Taper `eeg` with `taper`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}``

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
  * `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}``

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

Normalize each `eeg` channel by z-score.

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

<a id='NeuroJ.eeg_add_noise-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_add_noise-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_add_noise`** &mdash; *Method*.



```julia
eeg_add_noise(eeg)
```

Add random noise to the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_add_noise!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_add_noise!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_add_noise!`** &mdash; *Method*.



```julia
eeg_add_noise!(eeg)
```

Add random noise to the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_filter-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_filter-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_filter`** &mdash; *Method*.



```julia
eeg_filter(eeg; <keyword arguments>)
```

Filter `eeg` channels.

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
  * `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64=2`: filter order
  * `rp::Union{Int64, Float64}=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
  * `rs::Union{Int64, Float64}=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
  * `dir:Symbol=:twopass`: filter direction (`:onepass`, `:onepass_reverse`, `:twopass`), for causal filter use `:onepass`
  * `d::Int64=1`: window length for mean average and median average filter
  * `t::Union{Int64, Float64}`: threshold for `:mavg` and `:mmed` filters; threshold = threshold * std(signal) + mean(signal) for `:mavg` or threshold = threshold * std(signal) + median(signal) for `:mmed` filter
  * `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

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
  * `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64=2`: filter order
  * `rp::Union{Int64, Float64}=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
  * `rs::Union{Int64, Float64}=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
  * `dir:Symbol=:twopass`: filter direction (`:onepass`, `:onepass_reverse`, `:twopass`), for causal filter use `:onepass`
  * `d::Int64=1`: window length for mean average and median average filter
  * `t::Union{Int64, Float64}`: threshold for `:mavg` and `:mmed` filters; threshold = threshold * std(signal) + mean(signal) for `:mavg` or threshold = threshold * std(signal) + median(signal) for `:mmed` filter
  * `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

<a id='NeuroJ.eeg_pca-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_pca-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_pca`** &mdash; *Method*.



```julia
eeg_pca(eeg; n)
```

Calculate `n` first PCs for `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `n::Int64`: number of PCs

**Returns**

Named tuple containing:

  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pc_var::Matrix{Float64}`: PCVAR(1)..PCVAR(n) × epoch
  * `pc_m::Matrix{Float64}`: PC mean

<a id='NeuroJ.eeg_ica-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_ica-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_ica`** &mdash; *Method*.



```julia
eeg_ica(eeg; <keyword arguments>)
```

Calculate `n` first ICs for `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `n::Int64`: number of ICs
  * `tol::Float64=1.0e-6`: tolerance for ICA
  * `iter::Int64=100`: maximum number of iterations
  * `f::Symbol=:tanh`: neg-entropy functor: :tanh, :gaus

**Returns**

Named tuple containing:

  * `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
  * `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)

<a id='NeuroJ.eeg_average-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_average-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_average`** &mdash; *Method*.



```julia
eeg_average(eeg)
```

Return the average signal of all `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_average!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_average!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_average!`** &mdash; *Method*.



```julia
eeg_average!(eeg)
```

Return the average signal of all `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_average-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_average-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_average`** &mdash; *Method*.



```julia
eeg_average(eeg1, eeg2)
```

Return the average signal of all `eeg1` and `eeg2` channels.

**Arguments**

  * `eeg1::NeuroJ.EEG`
  * `eeg2::NeuroJ.EEG`

**Returns**

  * `eeg_new::NeuroJ.EEG`

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

<a id='NeuroJ.eeg_invert_polarity-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_invert_polarity-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_invert_polarity`** &mdash; *Method*.



```julia
eeg_invert_polarity(eeg; channel)
```

Invert polarity of `channel` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Int64`: channel to invert

**Returns**

  * `eeg_new::NeuroJ.EEG`

<a id='NeuroJ.eeg_invert_polarity!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_invert_polarity!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_invert_polarity!`** &mdash; *Method*.



```julia
eeg_invert_polarity!(eeg; channel)
```

Invert polarity of `channel` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to invert

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

  * `stp::Matrix{Float64}`: total power for each channel per epoch

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

  * `sbp::Matrix{Float64}`: band power for each channel per epoch

<a id='NeuroJ.eeg_cov-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_cov-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_cov`** &mdash; *Method*.



```julia
eeg_cov(eeg; norm=true)
```

Calculate covariance matrix for all channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `norm::Bool`: normalize covariance

**Returns**

  * `cov_mat::Array{Float64, 3}`: covariance matrix for each epoch

<a id='NeuroJ.eeg_cor-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_cor-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_cor`** &mdash; *Method*.



```julia
eeg_cor(eeg)
```

Calculate correlation coefficients between all channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `cov_mat::Array{Float64, 3}`: correlation matrix for each epoch

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

Named tuple containing:

  * `ccov::Matrix{Float64}`
  * `lags::Vector{Float64}

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

Named tuple containing:

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

Named tuple containing:

  * `psd_pow::Array{Float64, 3}`:powers
  * `psd_frq::Array{Float64, 3}`: frequencies

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

<a id='NeuroJ.eeg_band-Tuple{Any}' href='#NeuroJ.eeg_band-Tuple{Any}'>#</a>
**`NeuroJ.eeg_band`** &mdash; *Method*.



```julia
eeg_band(eeg, band)
```

Return frequency limits for a `band` range.

**Arguments**

  * `eeg:EEG`
  * `band::Symbol`: name of band range: :total, :delta, :theta, :alpha, :beta, :beta*high, :gamma, :gamma*1, :gamma*2, :gamma*lower, :gamma_higher. If lower or upper band frequency limit exceeds Nyquist frequency of `eeg`, than bound is truncated to `eeg` range.

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

  * `coh::Vector{ComplexF64}`

<a id='NeuroJ.eeg_freqs-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_freqs-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_freqs`** &mdash; *Method*.



```julia
eeg_freqs(eeg)
```

Return vector of frequencies and Nyquist frequency for `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

Named tuple containing:

  * `hz::Vector{Float64}`
  * `nyquist::Float64`

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

Named tuple containing:

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

Calculate `eeg` epochs statistics.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

Named tuple containing:

  * `e_mean::Vector(Float64)`: mean
  * `e_median::Vector(Float64)`: median
  * `e_std::Vector(Float64)`: standard deviation
  * `e_var::Vector(Float64)`: variance
  * `e_kurt::Vector(Float64)`: kurtosis
  * `e_mean_diff::Vector(Float64)`: mean diff value
  * `e_median_diff::Vector(Float64)`: median diff value
  * `e_max_dif::Vector(Float64)`: max difference
  * `e_dev_mean::Vector(Float64)`: deviation from channel mean

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

Named tuple containing:

  * `s_pow::Array{Float64, 3}`
  * `s_frq::Matrix{Float64}`
  * `s_t::Matrix{Float64}`

<a id='NeuroJ.eeg_spectrum-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_spectrum-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_spectrum`** &mdash; *Method*.



```julia
eeg_spectrum(eeg; pad)
```

Calculate FFT, amplitudes, powers and phases for each channel of the `eeg`. For `pad` > 0 channels are padded with 0s.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `pad::Int64=0`: pad with `pad` zeros

**Returns**

Named tuple containing:

  * `fft::Array{ComplexF64, 3}`
  * `amp::Array{Float64, 3}`
  * `pow::Array{Float64, 3}`
  * `phase::Array{Float64, 3}

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

<a id='NeuroJ.eeg_channels_stats-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_channels_stats-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_channels_stats`** &mdash; *Method*.



```julia
eeg_channels_stats(eeg)
```

Calculate `eeg` channels statistics.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

Named tuple containing:

  * `c_mean::Matrix(Float64)`: mean
  * `c_median::Matrix(Float64)`: median
  * `c_std::Matrix(Float64)`: standard deviation
  * `c_var::Matrix(Float64)`: variance
  * `c_kurt::Matrix(Float64)`: kurtosis
  * `c_mean_diff::Matrix(Float64)`: mean diff value
  * `c_median_diff::Matrix(Float64)`: median diff value
  * `c_max_dif::Matrix(Float64)`: max difference
  * `c_dev_mean::Matrix(Float64)`: deviation from channel mean

<a id='NeuroJ.eeg_snr-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_snr-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_snr`** &mdash; *Method*.



```julia
eeg_snr(eeg)
```

Calculate SNR of `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `snr::Matrix(Float64)`: SNR for each channel per epoch

**Source**

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278

<a id='NeuroJ.eeg_standardize-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_standardize-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_standardize`** &mdash; *Method*.



```julia
eeg_standardize(eeg)
```

Standardize `eeg` channels for ML.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg_new::NeuroJ.EEG`: standardized EEG
  * `scaler::Matrix{Float64}`: standardized EEG

<a id='NeuroJ.eeg_standardize!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_standardize!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_standardize!`** &mdash; *Method*.



```julia
eeg_standardize!(eeg)
```

Standardize `eeg` channels for ML.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `scaler::Matrix{Float64}`: standardized EEG

<a id='NeuroJ.eeg_fconv-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_fconv-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_fconv`** &mdash; *Method*.



```julia
eeg_fconv(eeg, kernel)
```

Perform convolution of all `eeg` channels in the frequency domain using `kernel`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel for convolution

**Returns**

  * `s_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal

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

  * `s_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal

<a id='NeuroJ.eeg_dft-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_dft-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_dft`** &mdash; *Method*.



```julia
eeg_make_spectrum(eeg)
```

Returns FFT and DFT sample frequencies for a DFT for each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

Named tuple containing:

  * `fft::Array{ComplexF64, 3}`: FFT
  * `sf::Array{Float64, 3}`: sample frequencies

<a id='NeuroJ.eeg_msci95-Tuple{Array{Float64, 3}}' href='#NeuroJ.eeg_msci95-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.eeg_msci95`** &mdash; *Method*.



```julia
eeg_msci95(eeg; n::=3, method=:normal)
```

Calculates mean, std and 95% confidence interval for each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `n::Int64`: number of bootstraps
  * `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

**Returns**

Named tuple containing:

  * `s_m::Matrix{Float64}`: mean
  * `s_s::Matrix{Float64}`: standard deviation
  * `s_u::Matrix{Float64}`: upper 95% CI
  * `s_l::Matrix{Float64}`: lower 95% CI

<a id='NeuroJ.eeg_mean-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_mean-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_mean`** &mdash; *Method*.



```julia
eeg_mean(eeg1, eeg2)
```

Calculates mean and 95% confidence interval for `eeg1` and `eeg2` channels.

**Arguments**

  * `eeg1::NeuroJ.EEG`
  * `eeg2:NeuroJ.EEG`

**Returns**

Named tuple containing:

  * `s_m::Matrix{Float64}`: mean by epochs
  * `s_s::Matrix{Float64}`: std by epochs
  * `s_u::Matrix{Float64}`: upper 95% CI bound by epochs
  * `s_l::Matrix{Float64}`: lower 95% CI bound by epochs

<a id='NeuroJ.eeg_difference-Tuple{Array{Float64, 3}, Array{Float64, 3}}' href='#NeuroJ.eeg_difference-Tuple{Array{Float64, 3}, Array{Float64, 3}}'>#</a>
**`NeuroJ.eeg_difference`** &mdash; *Method*.



```julia
eeg_difference(eeg1, eeg2; n=3, method=:absdiff)
```

Calculates mean difference and 95% confidence interval for `eeg1` and `eeg2`.

**Arguments**

  * `eeg1::Array{Float64, 3}`
  * `eeg2:Array{Float64, 3}`
  * `n::Int64`: number of bootstraps
  * `method::Symbol[:absdiff, :diff2int]`

      * `:absdiff`: maximum difference
      * `:diff2int`: integrated area of the squared difference

**Returns**

Named tuple containing:

  * `statistic::Matrix{Float64}`
  * `statistic_single::Vector{Float64}`
  * `p::Vector{Float64}`

<a id='NeuroJ.eeg_autocov-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_autocov-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_autocov`** &mdash; *Method*.



eeg_autocov(eeg; lag=1, demean=false, norm=false)

Calculate autocovariance of each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `lag::Int64`: lags range is `-lag:lag`
  * `demean::Bool`: demean eeg prior to analysis
  * `norm::Bool`: normalize autocovariance

**Returns**

Named tuple containing:

  * `acov::Matrix{Float64}`
  * `lags::Vector{Float64}`


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
  * `kwargs`: optional arguments for plot() function

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
  * `ylabel::String="Channels"`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal`** &mdash; *Method*.



```julia
eeg_plot_signal(eeg; <keyword arguments>)
```

Plot `eeg` channel or channels.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=0`: epochs to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component`** &mdash; *Method*.



```julia
eeg_plot_component(eeg; <keyword arguments>)
```

Plot `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `v::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component `v`
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_avg`** &mdash; *Method*.



```julia
eeg_plot_signal_avg(eeg; <keyword arguments>)
```

Plot averaged `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=1`: epoch number to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_avg`** &mdash; *Method*.



```julia
eeg_plot_component_avg(eeg; <keyword arguments>)
```

Plot averaged `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `v::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component `v`
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_butterfly-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_butterfly-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_butterfly`** &mdash; *Method*.



```julia
eeg_plot_signal_butterfly(eeg; <keyword arguments>)
```

Butterfly plot of `eeg` channels.

**Arguments**

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=1`: epoch number to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_butterfly-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_butterfly-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_butterfly`** &mdash; *Method*.



```julia
eeg_plot_component_butterfly(eeg; <keyword arguments>)
```

Butterfly plot of `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `v::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component `v`
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `norm::Bool=true`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_psd-Tuple{Vector{Float64}, Vector{Float64}}' href='#NeuroJ.signal_plot_psd-Tuple{Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroJ.signal_plot_psd`** &mdash; *Method*.



```julia
signal_plot_psd(s_pow, s_frq; <keyword arguments>)
```

Plot power spectrum density.

**Arguments**

  * `s_pow::Vector{Float64}`: signal powers
  * `s_frq::Vector{Float64}`: signal frequencies
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String="Power [μV^2/Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `epoch::Union{Int64, AbstractRange}=1`: epoch number to display
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
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_spectrogram-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_spectrogram-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_spectrogram`** &mdash; *Method*.



```julia
eeg_plot_spectrogram(eeg; <keyword arguments>)
```

Plots spectrogram of `eeg` channel(s).

**Arguments**

  * `eeg:EEG`
  * `epoch::Union{Int64, AbstractRange}=1`: epoch to plot
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize powers to dB
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: y-axis limits
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `epoch::Int64=0`: epoch number to display
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default 1 epoch or 20 seconds
  * `ic::Union{Int64, Vector{Int64}, AbstractRange, Nothing}=nothing`: which IC to plot, default all
  * `norm::Bool=true`: normalize the ICs prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

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
  * `epoch::Union{Int64, AbstractRange}=1`: epochs to display
  * `offset::Int64=1`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 second
  * `m::Symbol=:shepard`: interpolation method `:shepard` (Shepard), `:mq` (Multiquadratic), `:tp` (ThinPlate)
  * `c::Symbol=:amp`: component name (:ica, :pca, :amp, :power)
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange, Tuple, Nothing}=nothing`: component index, e.g. ICA number or frequency range
  * `norm::Bool=true`: convert power as dB
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: frequency limit for PSD and spectrogram
  * `head_labels::Bool=false`: plot head labels
  * `cb::Bool=false`: add color bars to plots
  * `cb_label::String=""`: color bar label
  * `average::Bool=true`: plot averaged signal and PSD
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.signal_plot_bands-Tuple{Vector{Float64}}' href='#NeuroJ.signal_plot_bands-Tuple{Vector{Float64}}'>#</a>
**`NeuroJ.signal_plot_bands`** &mdash; *Method*.



```julia
signal_plot_bands(signal; <keyword arguments>)
```

Plot absolute/relative bands powers of a single-channel `signal`.

**Arguments**

  * `signal::Vector{Float64}`
  * `fs::Int64`: sampling rate
  * `band::Vector{Symbol}=[:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]`: band names, e.g. [:delta, alpha](see `eeg_band()`)
  * `band_frq::Vector{Tuple{Float64, Float64}}`: vector of band frequencies
  * `type::Symbol`: plots absolute (:abs) or relative power (:rel)
  * `norm::Bool=true`: convert power to dB if true
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

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
  * `epoch::Union{Int64, AbstractRange}=1`: epochs to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channels to display
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `band:Vector{Symbols}=:all`: band name, e.g. :delta (see `eeg_band()`)
  * `type::Symbol`: plots absolute (:abs) or relative power (:rel)
  * `norm::Bool=true`: convert power to dB if true
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_channels-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_channels-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_channels`** &mdash; *Method*.



```julia
eeg_plot_channels(eeg; <keyword arguments>)
```

Plot values of `v` for selected channels of `eeg`.

**Arguments**

  * `eeg:NeuroJ.EEG`
  * `v::Union{Matrix{Int64}, Matrix{Float64}, Symbol}`: values to plot; if symbol, than use embedded component `v`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: list of channels to plot
  * `epoch::Int64`: number of epoch for which `v` should be plotted
  * `xlabel::String="Channels"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_epochs-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_epochs-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_epochs`** &mdash; *Method*.



```julia
eeg_plot_epochs(eeg; <keyword arguments>)
```

Plot values of `v` for selected epoch of `eeg`.

**Arguments**

  * `eeg:NeuroJ.EEG`
  * `v::Union{Vector{Int64}, Vector{Float64}, Symbol}`: values to plot; if symbol, than use embedded component `v`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: list of epochs to plot
  * `xlabel::String="Epochs"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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
  * `kwargs`: optional arguments for plot() function

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

<a id='NeuroJ.eeg_plot_compose-Tuple{Vector{Plots.Plot{Plots.GRBackend}}}' href='#NeuroJ.eeg_plot_compose-Tuple{Vector{Plots.Plot{Plots.GRBackend}}}'>#</a>
**`NeuroJ.eeg_plot_compose`** &mdash; *Method*.



```julia
eeg_plot_compose(p, l; <keyword arguments>)
```

Compose a complex plot of various plots contained in vector `p` using layout `layout`. Layout scheme is:

  * `(2, 2)`: 2 × 2 plots, regular layout
  * `@layout [a{0.2w} b{0.8w};_ c{0.6}]`: complex layout using Plots.jl `@layout` macro

**Arguments**

  * `p::Vector{Plots.Plot{Plots.GRBackend}}`: vector of plots
  * `layout::Union(Matrix{Any}, Tuple{Int64, Int64}}`: layout
  * `kwargs`: optional arguments for `p` vector plots

**Returns**

  * `pc::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_details-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_details-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_details`** &mdash; *Method*.



```julia
eeg_plot_signal_details(eeg; <keyword arguments>)
```

Plot details of `eeg` channels: amplitude, histogram, power density, phase histogram and spectrogram.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Union{Int64, AbstractRange}=0`: epoch number to display
  * `channel::Int64`: channel to display
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `head::Bool=true`: add head plot
  * `hist::Symbol=:hist`: histogram type: :hist, :kd
  * `frq_lim::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}=(0, 0)`: frequency limit for PSD and spectrogram
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `pc::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_avg_details-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_avg_details-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_avg_details`** &mdash; *Method*.



```julia
eeg_plot_avg_details(eeg; <keyword arguments>)
```

Plot details of averaged `eeg` channels: amplitude, histogram, power density, phase histogram and spectrogram.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Int64=1`: epoch number to display
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
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_butterfly_details-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_butterfly_details-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_butterfly_details`** &mdash; *Method*.



```julia
eeg_plot_signal_butterfly_details(eeg; <keyword arguments>)
```

Plot details butterfly plot of `eeg` channels: amplitude, histogram, power density, phase histogram and spectrogram.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Int64=1`: epoch number to display
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
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


<a id='Low-level-functions'></a>

<a id='Low-level-functions-1'></a>

## Low level functions

<a id='NeuroJ.linspace-Tuple{Union{Float64, Int64}, Union{Float64, Int64}, Int64}' href='#NeuroJ.linspace-Tuple{Union{Float64, Int64}, Union{Float64, Int64}, Int64}'>#</a>
**`NeuroJ.linspace`** &mdash; *Method*.



```julia
linspace(start, stop, length)
```

Generates `length`-long sequence of evenly spaced numbers between `start` and `stop`.

**Arguments**

  * `start::Union{Int64, Float64}`
  * `stop::Union{Int64, Float64}`
  * `length::Int64`

**Returns**

  * `range::Number`

<a id='NeuroJ.logspace-Tuple{Union{Float64, Int64}, Union{Float64, Int64}, Int64}' href='#NeuroJ.logspace-Tuple{Union{Float64, Int64}, Union{Float64, Int64}, Int64}'>#</a>
**`NeuroJ.logspace`** &mdash; *Method*.



```julia
logspace(start, stop, length)
```

Generates `length`-long sequence of log10-spaced numbers between `start` and `stop`.

**Arguments**

  * `start::Union{Int64, Float64}`
  * `stop::Union{Int64, Float64}`
  * `length::Int64`

**Returns**

  * `range::Number`

<a id='NeuroJ.m_pad0-Tuple{Union{Matrix{ComplexF64}, Matrix{Float64}, Matrix{Int64}}}' href='#NeuroJ.m_pad0-Tuple{Union{Matrix{ComplexF64}, Matrix{Float64}, Matrix{Int64}}}'>#</a>
**`NeuroJ.m_pad0`** &mdash; *Method*.



```julia
m_pad0(m)
```

Pad the matrix `m` with zeros to make it square.

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

<a id='NeuroJ.cart2pol-Tuple{Union{Float64, Int64}, Union{Float64, Int64}}' href='#NeuroJ.cart2pol-Tuple{Union{Float64, Int64}, Union{Float64, Int64}}'>#</a>
**`NeuroJ.cart2pol`** &mdash; *Method*.



```julia
cart2pol(x, y)
```

Convert cartographic coordinates `x` and `y` to polar.

**Arguments**

  * `x::Union{Int64, Float64}`
  * `y::Union{Int64, Float64}`

**Returns**

  * `phi::Float64`
  * `theta::Float64`

<a id='NeuroJ.pol2cart-Tuple{Union{Float64, Int64}, Union{Float64, Int64}}' href='#NeuroJ.pol2cart-Tuple{Union{Float64, Int64}, Union{Float64, Int64}}'>#</a>
**`NeuroJ.pol2cart`** &mdash; *Method*.



```julia
pol2cart(theta, phi)
```

Convert polar coordinates `theta` and `phi` to cartographic.

**Arguments**

  * `phi::Union{Float64, Int64}`
  * `theta::Union{Float64, Int64}`

**Returns**

  * `x::Float64`
  * `y::Float64`

<a id='NeuroJ.sph2cart' href='#NeuroJ.sph2cart'>#</a>
**`NeuroJ.sph2cart`** &mdash; *Function*.



```julia
pol2cart_sph(rho, theta, phi=0)
```

Convert spherical coordinates `theta` and `phi` and `rho` to cartographic.

**Arguments**

  * `phi::Union{Float64, Int64}`: the angle with respect to the z-axis (elevation)
  * `theta::Union{Float64, Int64}`: the angle in the xy plane with respect to the x-axis (azimuth)
  * `rho::Union{Float64, Int64}`: the distance from the origin to the point

**Returns**

  * `x::Float64`
  * `y::Float64`
  * `z::Float64`

<a id='NeuroJ.generate_window-Tuple{Symbol, Int64}' href='#NeuroJ.generate_window-Tuple{Symbol, Int64}'>#</a>
**`NeuroJ.generate_window`** &mdash; *Method*.



```julia
generate_window(type, n; even)
```

Return the `n`-point long symmetric window `type`.

**Arguments**

  * `type::Symbol`: window type:

      * `:hann`: Hann
      * `:bh`: Blackman-Harris
      * `:bohman`: Bohman
      * `:flat`: Flat-top window
      * `:bn`: Blackman-Nuttall
      * `:nutall`: Nuttall
  * `n::Int64`: window length
  * `even::Bool=false`: if true, make the window of even length (+1 for odd n)

**Returns**

  * `w::Vector{Float64}`:: generated window

<a id='NeuroJ.hildebrand_rule-Tuple{Union{Vector{Float64}, Vector{Int64}}}' href='#NeuroJ.hildebrand_rule-Tuple{Union{Vector{Float64}, Vector{Int64}}}'>#</a>
**`NeuroJ.hildebrand_rule`** &mdash; *Method*.



```julia
hildebrand_rule(x)
```

Calculate Hildebrand rule for vector `x`. If H < 0.2 then the vector `x` is symmetrical.

**Arguments**

  * `x::Union{Vector{Int64}, Vector{Float64}}`

**Returns**

  * `h::Float64`

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

<a id='NeuroJ.s_rms-Tuple{Union{Vector{Float64}, Vector{Int64}}}' href='#NeuroJ.s_rms-Tuple{Union{Vector{Float64}, Vector{Int64}}}'>#</a>
**`NeuroJ.s_rms`** &mdash; *Method*.



```julia
s_rms(signal)
```

Calculate Root Mean Square of `signal`.

**Arguments**

  * `signal::Union{Vector{Int64}, Vector{Float64}}`

**Returns**

  * rms::Float64`

<a id='NeuroJ.generate_sine' href='#NeuroJ.generate_sine'>#</a>
**`NeuroJ.generate_sine`** &mdash; *Function*.



```julia
generate_sine(f, t, a, p)
```

Generates sine wave of `f` frequency over `t` time; optional arguments are: `a` amplitude and  `p` phase.

**Arguments**

  * `f::Union{Int64, Float64}`
  * `t::Union{Vector{Int64}, Vector{Float64}}`
  * `a::Union{Int64, Float64}`
  * `p::Union{Int64, Float64}`

**Returns**

  * sine::Vector{Float64}`

<a id='NeuroJ.s_freqs-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}}' href='#NeuroJ.s_freqs-Tuple{Union{Vector{Float64}, Vector{Int64}, AbstractRange}}'>#</a>
**`NeuroJ.s_freqs`** &mdash; *Method*.



```julia
s_freqs(t)
```

Return vector of frequencies and Nyquist frequency for given time vector `t`.

**Arguments**

  * `t::Union{Vector{Int64}, Vector{Float64}, AbstractRange}`

**Returns**

  * `hz::Vector{Float64}
  * `nyquist_freq::Float64`

<a id='NeuroJ.s_freqs-Tuple{Vector{Float64}, Union{Float64, Int64}}' href='#NeuroJ.s_freqs-Tuple{Vector{Float64}, Union{Float64, Int64}}'>#</a>
**`NeuroJ.s_freqs`** &mdash; *Method*.



```julia
s_freqs(signal, fs)
```

Return vector of frequencies and Nyquist frequency for given `signal` and `fs`.

**Arguments**

  * `signal::Vector{Float64}`
  * `fs::Union{Int64, Float64}`

**Returns**

  * `hz::Vector{Float64}
  * `nyquist_freq::Float64`

<a id='NeuroJ.m_sortperm-Tuple{Matrix}' href='#NeuroJ.m_sortperm-Tuple{Matrix}'>#</a>
**`NeuroJ.m_sortperm`** &mdash; *Method*.



```julia
m_sortperm(m; rev=false, dims=1)
```

Generates sorting index for matrix `m` by columns (`dims` = 1) or by rows (`dims` = 2).

**Arguments**

  * `m::Matrix`
  * `rev::Bool`
  * `dims::Int64`

**Returns**

  * `m_idx::Matrix{Int64}`

<a id='NeuroJ.m_sort-Tuple{Matrix, Vector{Int64}}' href='#NeuroJ.m_sort-Tuple{Matrix, Vector{Int64}}'>#</a>
**`NeuroJ.m_sort`** &mdash; *Method*.



```julia
m_sort(m, m_idx; rev=false, dims=1)
```

Sorts matrix `m` using sorting index `m_idx` by columns (`dims` = 1) or by rows (`dims` = 2).

**Arguments**

  * `m::Matrix`
  * `m_idx::Vector{Int64}`
  * `rev::Bool`
  * `dims::Int64`

**Returns**

  * `m_sorted::Matrix`

<a id='NeuroJ.pad0' href='#NeuroJ.pad0'>#</a>
**`NeuroJ.pad0`** &mdash; *Function*.



```julia
pad0(x, n, sym)
```

Pad the vector `x` with `n` zeros.

**Arguments**

  * `x::Union{Vector{Int64}, Vector{Float64}}`
  * `n::Int64`
  * `sym::Bool=false`: if true, than pad at the beginning and at the end, otherwise only at the end.

**Returns**

  * `v_pad::Union{Vector{Int64}, Vector{Float64}}`

<a id='NeuroJ.hz2rads-Tuple{Union{Float64, Int64}}' href='#NeuroJ.hz2rads-Tuple{Union{Float64, Int64}}'>#</a>
**`NeuroJ.hz2rads`** &mdash; *Method*.



```julia
hz2rads(f)
```

Convert frequency `f` in Hz to rad/s.

**Arguments**

  * `f::Union{Int64, Float64}`

**Returns**

  * `f_rads::Float64`

<a id='NeuroJ.rads2hz-Tuple{Union{Float64, Int64}}' href='#NeuroJ.rads2hz-Tuple{Union{Float64, Int64}}'>#</a>
**`NeuroJ.rads2hz`** &mdash; *Method*.



```julia
rads2hz(f)
```

Convert frequency `f` in rad/s to Hz.

**Arguments**

  * `f::Union{Int64, Float64}`

**Returns**

  * `f_rads::Float64`

<a id='NeuroJ.z_score-Tuple{Union{Vector{Float64}, Vector{Int64}}}' href='#NeuroJ.z_score-Tuple{Union{Vector{Float64}, Vector{Int64}}}'>#</a>
**`NeuroJ.z_score`** &mdash; *Method*.



```julia
z_score(x)
```

Calculate Z-scores for each value of the vector `x`.

**Arguments**

  * `x::Union{Vector{Int64}, Vector{Float64}}`

**Returns**

  * `z_score::Vector{Float64}`

<a id='NeuroJ.k_categories-Tuple{Int64}' href='#NeuroJ.k_categories-Tuple{Int64}'>#</a>
**`NeuroJ.k_categories`** &mdash; *Method*.



```julia
k_categories(n)
```

Calculate number of categories for a given sample size `n`.

**Arguments**

  * `n::Int64`

**Returns**

  * `k::Float64`

<a id='NeuroJ.cmax-Tuple{Vector{ComplexF64}}' href='#NeuroJ.cmax-Tuple{Vector{ComplexF64}}'>#</a>
**`NeuroJ.cmax`** &mdash; *Method*.



```julia
cmax(x)
```

Return maximum value of the complex vector`x`.

**Arguments**

  * `x::Vector{ComplexF64}`

**Returns**

  * `cmax::ComplexF64`

<a id='NeuroJ.cmin-Tuple{Vector{ComplexF64}}' href='#NeuroJ.cmin-Tuple{Vector{ComplexF64}}'>#</a>
**`NeuroJ.cmin`** &mdash; *Method*.



```julia
cmin(x)
```

Return minimum value of the complex vector`x`.

**Arguments**

  * `x::Vector{ComplexF64}`

**Returns**

  * `cmin::ComplexF64`

<a id='NeuroJ.generate_sinc' href='#NeuroJ.generate_sinc'>#</a>
**`NeuroJ.generate_sinc`** &mdash; *Function*.



```julia
generate_sinc(t, f, peak, norm)
```

Generate normalized or unnormalized sinc function.

**Arguments**

  * `t::AbstractRange=-2:0.01:2`: time
  * `f::Union{Int64, Float64}=10.0`: frequency
  * `peak::Union{Int64, Float64}=0`: sinc peak time
  * `norm::Bool=true`: generate normalzied function

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

<a id='NeuroJ.generate_gaussian' href='#NeuroJ.generate_gaussian'>#</a>
**`NeuroJ.generate_gaussian`** &mdash; *Function*.



```julia
generate_gaussian(fs, gt, gw, pt, pa)
```

Generate Gaussian wave.

**Arguments**

  * `fs::Int64`: sampling rate
  * `gt::Union{Int64, Float64}`: length = 0:1/fs:gt
  * `gw::Union{Int64, Float64}=1`: width
  * `pt::Union{Int64, Float64}=0`: peak time
  * `pa::Union{Int64, Float64}=1`: peak amp
  * 

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

<a id='NeuroJ.s2_rmse-Tuple{Vector{Float64}, Vector{Float64}}' href='#NeuroJ.s2_rmse-Tuple{Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroJ.s2_rmse`** &mdash; *Method*.



```julia
s2_rmse(signal1, signal2)
```

Calculate RMSE between `signal1` and `signal2`.

**Arguments**

  * `signal1::Vector{Float64}`
  * `signal2::Vector{Float64}`

**Returns**

  * `r::Float64`

<a id='NeuroJ.m_norm-Tuple{Array{Float64, 3}}' href='#NeuroJ.m_norm-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.m_norm`** &mdash; *Method*.



```julia
m_norm(m)
```

Normalize matrix `m`.

**Arguments**

  * `m::Matrix{Float64}`

**Returns**

  * `m_norm::Matrix{Float64}`

<a id='NeuroJ.s_cov-Tuple{AbstractArray}' href='#NeuroJ.s_cov-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_cov`** &mdash; *Method*.



s_cov(signal; norm=true)

Calculates covariance between all channels of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `norm::Bool`: normalize covariance

**Returns**

  * `cov_mat::Matrix{Float64}`

<a id='$-Tuple{Any, Any}' href='#$-Tuple{Any, Any}'>#</a>
**`$`** &mdash; *Method*.



```julia
$
```

Interpolation operator for interpolating into e.g. [strings](@ref string-interpolation) and [expressions](@ref man-expression-interpolation).

**Examples**

```julia-repl
julia> name = "Joe"
"Joe"

julia> "My name is $name."
"My name is Joe."
```


<a target='_blank' href='https://github.com/JuliaLang/julia/blob/bf534986350a991e4a1b29126de0342ffd76205e/base/docs/basedocs.jl#L559-L573' class='documenter-source'>source</a><br>

<a id='NeuroJ.s_msci95-Tuple{Vector{Float64}}' href='#NeuroJ.s_msci95-Tuple{Vector{Float64}}'>#</a>
**`NeuroJ.s_msci95`** &mdash; *Method*.



```julia
s_msci95(signal)
```

Calculate mean, std and 95% confidence interval for `signal`.

**Arguments**

  * `signal::Vector{Float64}`

**Returns**

  * `s_m::Float64`: mean
  * `s_s::Float64`: standard deviation
  * `s_u::Float64`: upper 95% CI
  * `s_l::Float64`: lower 95% CI

<a id='NeuroJ.s_msci95-Tuple{AbstractArray}' href='#NeuroJ.s_msci95-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_msci95`** &mdash; *Method*.



```julia
s_msci95(signal; n=3, method=:normal)
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

<a id='NeuroJ.s2_mean-Tuple{Vector{Float64}, Vector{Float64}}' href='#NeuroJ.s2_mean-Tuple{Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroJ.s2_mean`** &mdash; *Method*.



```julia
s2_mean(signal1, signal2)
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

<a id='NeuroJ.s2_difference-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s2_difference-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s2_difference`** &mdash; *Method*.



```julia
s2_difference(signal1, signal2; n=3, method=:absdiff)
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

<a id='NeuroJ.s_acov-Tuple{AbstractArray}' href='#NeuroJ.s_acov-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_acov`** &mdash; *Method*.



s_acov(signal; lag=1, demean=false, norm=false)

Calculate autocovariance of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `lag::Int64`: lags range is `-lag:lag`
  * `demean::Bool`: demean `signal` prior to calculations
  * `norm::Bool`: normalize autocovariance

**Returns**

  * `acov::Vector{Float64}`
  * `lags::Vector{Int64}`

<a id='NeuroJ.s_xcov-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s_xcov-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s_xcov`** &mdash; *Method*.



s_xcov(signal1, signal2; lag=1, demean=false, norm=false)

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

<a id='NeuroJ.s_spectrum-Tuple{AbstractArray}' href='#NeuroJ.s_spectrum-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_spectrum`** &mdash; *Method*.



```julia
s_spectrum(signal; pad=0)
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

<a id='NeuroJ.s_total_power-Tuple{AbstractArray}' href='#NeuroJ.s_total_power-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_total_power`** &mdash; *Method*.



```julia
s_total_power(signal; fs)
```

Calculates `signal` total power.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate

**Returns**

  * `stp::Float64`: signal total power

<a id='NeuroJ.s_band_power-Tuple{AbstractArray}' href='#NeuroJ.s_band_power-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_band_power`** &mdash; *Method*.



```julia
s_band_power(signal; fs, f)
```

Calculates `signal` power between `f[1]` and `f[2]`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate
  * `f::Tuple`: lower and upper frequency bounds

**Returns**

  * `stp::Float64`: signal total power

<a id='NeuroJ.s_taper-Tuple{AbstractArray}' href='#NeuroJ.s_taper-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_taper`** &mdash; *Method*.



```julia
s_taper(signal; taper)
```

Taper the `signal` with `taper`.

**Arguments**

  * `signal::AbstractArray`
  * `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

**Returns**

  * `s_tapered::Vector{Union{Float64, ComplexF64}}`

<a id='NeuroJ.s_detrend-Tuple{AbstractArray}' href='#NeuroJ.s_detrend-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_detrend`** &mdash; *Method*.



```julia
s_detrend(signal; type, offset, order, span)
```

Perform piecewise detrending of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Symbol`, optional

      * `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `:linear`: linear trend is subtracted from `signal`
      * `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
      * `:poly`: polynomial of `order` is subtracted
      * `:loess`: fit and subtract loess approximation
  * `offset::Union{Int64, Float64}=0`: constant for :constant detrending
  * `order::Int64=1`: polynomial fitting order
  * `span::Float64=0.5`: smoothing of loess

**Returns**

  * `s_det::Vector{Float64}`

<a id='NeuroJ.s_demean-Tuple{AbstractArray}' href='#NeuroJ.s_demean-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_demean`** &mdash; *Method*.



```julia
s_demean(signal)
```

Remove mean value (DC offset) from the `signal`.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_demeaned::Vector{Float64}`

<a id='NeuroJ.s_normalize_zscore-Tuple{AbstractArray}' href='#NeuroJ.s_normalize_zscore-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_normalize_zscore`** &mdash; *Method*.



```julia
s_normalize_zscore(signal)
```

Normalize (by z-score) `signal`.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroJ.s_normalize_minmax-Tuple{AbstractArray}' href='#NeuroJ.s_normalize_minmax-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_normalize_minmax`** &mdash; *Method*.



```julia
s_normalize_minmax(signal)
```

Normalize `signal` in [-1, +1].

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroJ.s_add_noise-Tuple{AbstractArray}' href='#NeuroJ.s_add_noise-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_add_noise`** &mdash; *Method*.



```julia
s_add_noise(signal)
```

Adds random noise to the `signal`.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_noisy::Vector{Float64}`

<a id='NeuroJ.s_resample-Tuple{AbstractArray}' href='#NeuroJ.s_resample-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_resample`** &mdash; *Method*.



```julia
s_resample(signal; t, new_sr)
```

Resample `signal` to `new_sr` sampling frequency.

**Arguments**

  * `signal::AbstractArray`
  * `t::AbstractRange`: time
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `s_resampled::Vector{Float64}`
  * `t_resampled::AbstractRange`

<a id='NeuroJ.s_resample-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_resample-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_resample`** &mdash; *Method*.



```julia
s_resample(signal; t, new_sr)
```

Resamples all channels of the`signal` and time vector `t` to `new_sr` sampling frequency.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `t::AbstractRange`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `s_downsampled::Array{Float64, 3}`
  * `t_downsampled::AbstractRange`

<a id='NeuroJ.s_derivative-Tuple{AbstractArray}' href='#NeuroJ.s_derivative-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_derivative`** &mdash; *Method*.



```julia
s_derivative(signal)
```

Return derivative of `signal` of the same length.

<a id='NeuroJ.s_tconv-Tuple{AbstractArray}' href='#NeuroJ.s_tconv-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_tconv`** &mdash; *Method*.



```julia
s_tconv(signal; kernel)
```

Performs convolution in the time domain between `signal` and `kernel`.

**Arguments**

  * `signal::AbstractArray`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

**Returns**

  * `s_conv::Union{Vector{Float64}, Vector{ComplexF64}}`

<a id='NeuroJ.s_filter-Tuple{AbstractArray}' href='#NeuroJ.s_filter-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_filter`** &mdash; *Method*.



```julia
s_filter(signal; <keyword arguments>)
```

Filter `signal`.

**Arguments**

  * `signal::AbstractArray`
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
  * `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64=8`: filter order
  * `rp::Union{Int64, Float64}=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
  * `rs::Union{Int64, Float64}=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
  * `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass), for causal filter use :onepass
  * `d::Int64=1`: window length for mean average and median average filter
  * `t::Union{Int64, Float64}`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
  * `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

**Returns**

  * `s_filtered::Vector{Float64}`

<a id='NeuroJ.s_psd-Tuple{AbstractArray}' href='#NeuroJ.s_psd-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_psd`** &mdash; *Method*.



```julia
s_psd(signal; fs, norm=false)
```

Calculate power spectrum density of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate
  * `norm::Bool`: normalize do dB

**Returns**

  * `psd_pow::Vector{Float64}`
  * `psd_frq::Vector{Float64}`

<a id='NeuroJ.s_psd-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_psd-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_psd`** &mdash; *Method*.



```julia
s_psd(signal; fs, norm=false)
```

Calculate power spectrum density of the `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `fs::Int64`: sampling rate
  * `norm::Bool`: normalize do dB

**Returns**

  * `psd_pow::Array{Float64, 3}`
  * `psd_frq::Array{Float64, 3}`

<a id='NeuroJ.s_stationarity_hilbert-Tuple{AbstractArray}' href='#NeuroJ.s_stationarity_hilbert-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_stationarity_hilbert`** &mdash; *Method*.



```julia
s_stationarity_hilbert(signal::Vector{Float64})
```

Calculate phase stationarity using Hilbert transformation.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `phase_stationarity::Vector{Float64}`

<a id='NeuroJ.s_stationarity_mean-Tuple{AbstractArray}' href='#NeuroJ.s_stationarity_mean-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_stationarity_mean`** &mdash; *Method*.



```julia
s_stationarity_mean(signal)
```

Calculate mean stationarity.

**Arguments**

  * `signal::AbstractArray`
  * `window::Int64`: time window in samples

**Returns**

  * `mean_stationarity::Vector{Float64}`

<a id='NeuroJ.s_stationarity_var-Tuple{AbstractArray}' href='#NeuroJ.s_stationarity_var-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_stationarity_var`** &mdash; *Method*.



```julia
s_stationarity_var(signal)
```

Calculate variance stationarity.

**Arguments**

  * `signal::AbstractArray`
  * `window::Int64`: time window in samples

**Returns**

  * `var_stationarity::Vector{Float64}`

<a id='NeuroJ.s_trim-Tuple{AbstractArray}' href='#NeuroJ.s_trim-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_trim`** &mdash; *Method*.



```julia
s_trim(signal; len::Int64)
```

Remove `len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `len::Int64`: trimming length in samples
  * `offset::Int64`: offset from which trimming starts, only works for `from` = :start
  * `from::Symbol[:start, :end]

**Returns**

  * `s_trimmed::Vector{Float64}`

<a id='NeuroJ.s2_mi-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s2_mi-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s2_mi`** &mdash; *Method*.



```julia
s2_mi(signal1::AbstractArray, signal2::AbstractArray)
```

Calculate mutual information between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

  * `mi::Float64`

<a id='NeuroJ.s_entropy-Tuple{AbstractArray}' href='#NeuroJ.s_entropy-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_entropy`** &mdash; *Method*.



```julia
s_entropy(signal)
```

Calculate entropy of `signal`.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `ent::Float64`

<a id='NeuroJ.s_average-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_average-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_average`** &mdash; *Method*.



```julia
s_average(signal)
```

Average all channels of `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `s_averaged::Array{Float64, 3}`

<a id='NeuroJ.s2_average-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s2_average-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s2_average`** &mdash; *Method*.



```julia
s2_average(signal1, signal2)
```

Averages `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

  * `s_averaged::Vector{Float64}`

<a id='NeuroJ.s2_coherence-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s2_coherence-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s2_coherence`** &mdash; *Method*.



```julia
s2_coherence(signal1, signal2)
```

Calculate coherence between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

  * `coherence::Vector{ComplexF64}`

<a id='NeuroJ.s_pca-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_pca-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_pca`** &mdash; *Method*.



```julia
s_pca(signal, n)
```

Calculates `n` first PCs for `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `n::Int64`: number of PCs

**Returns**

  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pc_var::Matrix{Float64}`: PC*VAR(1)..PC*VAR(n) × epoch

<a id='NeuroJ.s_fconv-Tuple{AbstractArray}' href='#NeuroJ.s_fconv-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_fconv`** &mdash; *Method*.



```julia
s_fconv(signal; kernel)
```

Perform convolution in the frequency domain between `signal` and `kernel`.

**Arguments**

  * `signal::AbstractArray`
  * `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

**Returns**

  * `s_conv::Vector{ComplexF64}`

<a id='NeuroJ.s_ica-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_ica-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_ica`** &mdash; *Method*.



```julia
s_ica(signal, n, tol=1.0e-6, iter=100, f=:tanh)
```

Calculate `n` first ICs for `signal`.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `n::Int64`: number of PCs
  * `tol::Float64`: tolerance for ICA
  * `iter::Int64`: maximum number of iterations
  * `f::Symbol[:tanh, :gaus]`: neg-entropy functor

**Returns**

  * `ic::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
  * `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch

<a id='NeuroJ.s_ica_reconstruct-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_ica_reconstruct-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_ica_reconstruct`** &mdash; *Method*.



```julia
s_ica_reconstruct(signal, ic, ic_mw, ic_v)
```

Reconstructs `signal` using removal of `ic_v` ICA components.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `ic::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
  * `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
  * `ic_v::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

**Returns**

  * `s_reconstructed::Array{Float64, 3}`

<a id='NeuroJ.s_spectrogram-Tuple{AbstractArray}' href='#NeuroJ.s_spectrogram-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_spectrogram`** &mdash; *Method*.



```julia
s_spectrogram(signal; fs, norm=true, demean=true)
```

Calculate spectrogram of `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling frequency
  * `norm::Bool`: normalize powers to dB
  * `demean::Bool`: demean signal prior to analysis

**Returns**

  * `s_pow::Matrix{Float64}`: powers
  * `s_frq::Vector{Float64}`: frequencies
  * `s_t::Vector{Float64}`: time

<a id='NeuroJ.s_detect_epoch_flat-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_flat-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_flat`** &mdash; *Method*.



```julia
s_detect_epoch_flat(signal::Array{Float64, 3}, threshold=0.1)
```

Detect bad `signal` epochs based on: flat channel(s)

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.s_detect_epoch_rmse-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_rmse-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_rmse`** &mdash; *Method*.



```julia
s_detect_epoch_rmse(signal::Array{Float64, 3})
```

Detect bad `signal` epochs based on: RMSE vs average channel > 95%CI.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.s_detect_epoch_rmsd-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_rmsd-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_rmsd`** &mdash; *Method*.



```julia
detect_epoch_rmsd(signal::Array{Float64, 3})
```

Detect bad `signal` epochs based on: RMSD vs average channel > 95%CI.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.s_detect_epoch_euclid-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_euclid-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_euclid`** &mdash; *Method*.



```julia
s_detect_epoch_euclid(signal::Array{Float64, 3})
```

Detect bad `signal` epochs based on: Euclidean distance vs median channel > 95% CI.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.s_detect_epoch_p2p-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_p2p-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_p2p`** &mdash; *Method*.



```julia
s_detect_epoch_p2p(signal::Array{Float64, 3})
```

Detect bad `signal` epochs based on: p2p amplitude > upper 95% CI p2p amplitude.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.s_snr-Tuple{AbstractArray}' href='#NeuroJ.s_snr-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_snr`** &mdash; *Method*.



```julia
s_snr(signal)
```

Calculate SNR of `signal`.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `snr::Float64`: SNR

**Source**

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278


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

Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.

