
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
  * `average::Bool`: average all epochs, return one averaged epoch; if false than return array of epochs, each row is one epoch

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

<a id='NeuroJ.eeg_epochs_time-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_epochs_time-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_epochs_time`** &mdash; *Method*.



```julia
eeg_epochs_time(eeg; ts)
```

Edit `eeg` epochs time start.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `ts::Real`: time start in seconds

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
  * `ts::Real`: time start in seconds

**Returns**

  * `eeg::NeuroJ.EEG`


<a id='EEG-process'></a>

<a id='EEG-process-1'></a>

## EEG process

<a id='NeuroJ.eeg_reference_ch-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_ch-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_ch`** &mdash; *Method*.



```julia
eeg_reference_ch(eeg; channel)
```

Reference the `eeg` to specific `channel`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_reference_ch!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_ch!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_ch!`** &mdash; *Method*.



```julia
eeg_reference_ch!(eeg; channel)
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
  * `exclude_fpo::Bool=true`: exclude Fp1, Fp2, O1, O2 from CAR mean calculation
  * `exclude_current::Bool=true`: exclude current electrode from CAR mean calculation

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

<a id='NeuroJ.eeg_reference_a-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_a-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_a`** &mdash; *Method*.



```julia
eeg_reference_a(eeg; type)
```

Reference the `eeg` to auricular channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Symbol=:link`: :l (linked, average of A1 and A2), :i (ipsilateral, A1 for left channels) or :c (contraletral, A1 for right channels)

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_reference_a!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_a!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_a!`** &mdash; *Method*.



```julia
eeg_reference_a!(eeg; type)
```

Reference the `eeg` to auricular channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Symbol=:link`: :l (linked, average of A1 and A2), :i (ipsilateral, A1 for left channels) or :c (contraletral, A1 for right channels)

<a id='NeuroJ.eeg_reference_m-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_m-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_m`** &mdash; *Method*.



```julia
eeg_reference_m(eeg; type)
```

Reference the `eeg` to mastoid channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Symbol=:link`: :l (linked, average of M1 and M2), :i (ipsilateral, M1 for left channels) or :c (contraletral, M1 for right channels)

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_reference_m!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_reference_m!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_reference_m!`** &mdash; *Method*.



```julia
eeg_reference_m!(eeg; type)
```

Reference the `eeg` to mastoid channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Symbol=:link`: :l (linked, average of M1 and M2), :i (ipsilateral, M1 for left channels) or :c (contraletral, M1 for right channels)

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
  * `offset::Real=0`: constant for :constant detrending
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
  * `offset::Real=0`: constant for :constant detrending
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
  * `taper::Union{Vector{Real, Vector{ComplexF64}}``

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
  * `taper::Union{Vector{<:Real}, Vector{ComplexF64}}``

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
      * `:iirnotch`
      * `:remez`
      * `:mavg`: moving average (with threshold)
      * `:mmed`: moving median (with threshold)
      * `:poly`: polynomial of `order` order
  * `ftype::Symbol`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64=8`: filter order, number of taps for :remez filter, k-value for :mavg and :mmed (window length = 2 × k + 1)
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
  * `bw::Real=-1`: bandwidth for :iirnotch and :remez filters
  * `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass), for causal filter use :onepass
  * `t::Real`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
  * `window::Union{Vector{<:Real}, Nothing} - window, required for FIR filter

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
      * `:iirnotch`
      * `:remez`
      * `:mavg`: moving average (with threshold)
      * `:mmed`: moving median (with threshold)
      * `:poly`: polynomial of `order` order
  * `ftype::Symbol`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64=8`: filter order, number of taps for :remez filter, k-value for :mavg and :mmed (window length = 2 × k + 1)
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
  * `bw::Real=-1`: bandwidth for :iirnotch and :remez filters
  * `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass), for causal filter use :onepass
  * `t::Real`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
  * `window::Union{Vector{<:Real}, Nothing} - window, required for FIR filter

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
  * `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
  * `pc_m::PCA{Float64}`: PC mean

<a id='NeuroJ.eeg_pca_reconstruct-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_pca_reconstruct-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_pca_reconstruct`** &mdash; *Method*.



```julia
eeg_pca_reconstruct(eeg)
```

Reconstruct `eeg` signals using PCA components.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `eeg::NeuroJ.EEG`

<a id='NeuroJ.eeg_pca_reconstruct!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_pca_reconstruct!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_pca_reconstruct!`** &mdash; *Method*.



```julia
eeg_pca_reconstruct!(eeg)
```

Reconstruct `eeg` signals using PCA components.

**Arguments**

  * `eeg::NeuroJ.EEG`

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

<a id='NeuroJ.eeg_wdenoise-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_wdenoise-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_wdenoise`** &mdash; *Method*.



```julia
eeg_wdenoise(eeg; wt)
```

Perform wavelet denoising.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `wt::Symbol=:db4`: wavelet type: :db2, :db4, :db8, :db10, :haar, :coif2, :coif4, :coif8

**Returns**

  * `eeg_new::NeuroJ.EEG`

<a id='NeuroJ.eeg_wdenoise!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_wdenoise!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_wdenoise!`** &mdash; *Method*.



```julia
eeg_wdenoise!(eeg; wt)
```

Perform wavelet denoising.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `wt::Symbol=:db4`: wavelet type: db2, db4, db8, db10, haar

<a id='NeuroJ.eeg_fftdenoise-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_fftdenoise-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_fftdenoise`** &mdash; *Method*.



```julia
eeg_fftdenoise(eeg; pad, threshold)
```

Perform wavelet denoising.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `pad::Int64=0`: pad signal with `pad` zeros
  * `threshold::Int64=100`: PSD threshold for keeping frequency components

**Returns**

  * `eeg_new::NeuroJ.EEG`

<a id='NeuroJ.eeg_fftdenoise!-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_fftdenoise!-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_fftdenoise!`** &mdash; *Method*.



```julia
eeg_fftdenoise!(eeg; pad, threshold)
```

Perform wavelet denoising.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `pad::Int64=0`: pad signal with `pad` zeros
  * `threshold::Int64=100`: PSD threshold for keeping frequency components


<a id='EEG-analyze'></a>

<a id='EEG-analyze-1'></a>

## EEG analyze

<a id='NeuroJ.eeg_total_power' href='#NeuroJ.eeg_total_power'>#</a>
**`NeuroJ.eeg_total_power`** &mdash; *Function*.



```julia
eeg_total_power(eeg, mt)
```

Calculate total power of the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

  * `stp::Matrix{Float64}`: total power for each channel per epoch

<a id='NeuroJ.eeg_band_power-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_band_power-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_band_power`** &mdash; *Method*.



```julia
eeg_band_power(eeg; f, mt)
```

Calculate absolute band power between frequencies `f[1]` and `f[2]` of the `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `mt::Bool=false`: if true use multi-tapered periodogram

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
  * `lags::Vector{Float64}`

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
  * `lags::Vector{Float64}`

<a id='NeuroJ.eeg_psd-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_psd-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_psd`** &mdash; *Method*.



```julia
eeg_psd(eeg; norm)
```

Calculate total power for each the `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `norm::Bool=false`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram

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

<a id='NeuroJ.eeg_negentropy-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_negentropy-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_negentropy`** &mdash; *Method*.



```julia
eeg_negentropy(eeg)
```

Calculate negentropy of all channels of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `ne::Matrix{Float64}`

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

  * `band_frequency::Tuple{Real, Real}`

<a id='NeuroJ.eeg_tcoherence-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_tcoherence-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_tcoherence`** &mdash; *Method*.



```julia
eeg_tcoherence(eeg1, eeg2)
```

Calculate coherence (mean over time) and MSC (magnitude-squared coherence) between all channels of `eeg1` and `eeg2`.

**Arguments**

  * `eeg1::NeuroJ.EEG`
  * `eeg2::NeuroJ.EEG`

**Returns**

Named tuple containing:

  * `c::Array{Float64, 3}`: coherence
  * `msc::Array{Float64, 3}`: MSC
  * `ic::Array{Float64, 3}`: imaginary part of coherence

<a id='NeuroJ.eeg_tcoherence-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_tcoherence-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_tcoherence`** &mdash; *Method*.



```julia
eeg_tcoherence(eeg; channel1, channel2, epoch1, epoch2)
```

Calculate coherence (mean over time), IC (imaginary coherence) and MSC (magnitude-squared coherence) between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel1::Int64`
  * `channel2::Int64`
  * `epoch1::Int64`
  * `epoch2::Int64`

**Returns**

Named tuple containing:

  * `c::Vector{Float64}`: coherence
  * `msc::Vector{Float64}`: MSC
  * `ic::Vector{Float64}`: imaginary part of coherence

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
  * `e_skew::Vector(Float64)`: skewness
  * `e_mean_diff::Vector(Float64)`: mean diff value
  * `e_median_diff::Vector(Float64)`: median diff value
  * `e_max_dif::Vector(Float64)`: max difference
  * `e_dev_mean::Vector(Float64)`: deviation from channel mean

<a id='NeuroJ.eeg_spectrogram-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_spectrogram-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_spectrogram`** &mdash; *Method*.



```julia
eeg_spectrogram(eeg; norm, mt, demean)
```

Return spectrogram of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `norm::Bool`=true: normalize powers to dB
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `demean::Bool`=true: demean signal prior to analysis

**Returns**

Named tuple containing:

  * `s_pow::Array{Float64, 3}`
  * `s_frq::Matrix{Float64}`
  * `s_t::Matrix{Float64}`

<a id='NeuroJ.eeg_spectrum-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_spectrum-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_spectrum`** &mdash; *Method*.



```julia
eeg_spectrum(eeg; pad, h)
```

Calculate FFT, amplitudes, powers and phases for each channel of the `eeg`. For `pad` > 0 channels are padded with 0s.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `pad::Int64=0`: pad with `pad` zeros
  * `h::Bool=false`: use Hilbert transform for calculations instead of FFT

**Returns**

Named tuple containing:

  * `fft::Array{ComplexF64, 3}`: Fourier or Hilbert components
  * `amp::Array{Float64, 3}`: amplitudes
  * `pow::Array{Float64, 3}`: powers
  * `phase::Array{Float64, 3}: phase angles

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
  * `t::Real`: time in seconds

**Returns**

  * `t_s::Int64`: time in samples

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
  * `c_skew::Matrix(Float64)`: skewness
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
eeg_fconv(eeg, kernel, norm)
```

Perform convolution of all `eeg` channels in the frequency domain using `kernel`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel for convolution
  * `norm::Bool=false`: normalize kernel

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
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

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

<a id='NeuroJ.eeg_tenv-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_tenv-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_tenv`** &mdash; *Method*.



```julia
eeg_tenv(eeg; d)
```

Calculate temporal envelope of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `t_env::Array{Float64, 3}`: temporal envelope
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroJ.eeg_tenv_mean-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_tenv_mean-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_tenv_mean`** &mdash; *Method*.



```julia
eeg_tenv_mean(eeg; dims, d)
```

Calculate temporal envelope of `eeg`: mean and 95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
  * `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
  * `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroJ.eeg_tenv_median-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_tenv_median-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_tenv_median`** &mdash; *Method*.



```julia
eeg_tenv_median(eeg; dims, d)
```

Calculate temporal envelope of `eeg`: median and 95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
  * `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
  * `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroJ.eeg_penv-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_penv-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_penv`** &mdash; *Method*.



```julia
eeg_penv(eeg; d)
```

Calculate power (in dB) envelope of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `p_env::Array{Float64, 3}`: power spectrum envelope
  * `p_env_frq::Vector{Float64}`: frequencies for each envelope

<a id='NeuroJ.eeg_penv_mean-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_penv_mean-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_penv_mean`** &mdash; *Method*.



```julia
eeg_penv_mean(eeg; dims, d)
```

Calculate power (in dB) envelope of `eeg`: mean and 95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `p_env_m::Array{Float64, 3}`: power spectrum envelope: mean
  * `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
  * `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
  * `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)

<a id='NeuroJ.eeg_penv_median-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_penv_median-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_penv_median`** &mdash; *Method*.



```julia
eeg_penv_median(eeg; dims, d)
```

Calculate power (in dB) envelope of `eeg`: median and 95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
  * `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `p_env_m::Array{Float64, 3}`: power spectrum envelope: median
  * `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
  * `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
  * `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)

<a id='NeuroJ.eeg_senv-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_senv-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_senv`** &mdash; *Method*.



```julia
eeg_senv(eeg; d, mt)
```

Calculate spectral (in dB) envelope of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered spectrogram

**Returns**

Named tuple containing:

  * `s_env::Array{Float64, 3}`: spectral envelope
  * `s_env_t::Vector{Float64}`: spectrogram time

<a id='NeuroJ.eeg_senv_mean-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_senv_mean-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_senv_mean`** &mdash; *Method*.



```julia
eeg_senv_mean(eeg; dims, d, mt)
```

Calculate spectral (in dB) envelope of `eeg`: mean and 95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered spectrogram

**Returns**

Named tuple containing:

  * `s_env_m::Array{Float64, 3}`: spectral envelope: mean
  * `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
  * `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
  * `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)

<a id='NeuroJ.eeg_senv_median-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_senv_median-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_senv_median`** &mdash; *Method*.



```julia
eeg_senv_median(eeg; dims, d, mt)
```

Calculate spectral (in dB) envelope of `eeg`: median and 95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `dims::Int64`: mean over chan (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered spectrogram

**Returns**

Named tuple containing:

  * `s_env_m::Array{Float64, 3}`: spectral envelope: median
  * `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
  * `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
  * `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)

<a id='NeuroJ.eeg_ispc-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_ispc-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_ispc`** &mdash; *Method*.



```julia
eeg_ispc(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculate ISPC (Inter-Site-Phase Clustering) between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

**Arguments**

  * `eeg1::NeuroJ.EEG`
  * `eeg2::NeuroJ.EEG`
  * `channel1::Int64`
  * `channel2::Int64`
  * `epoch1::Int64`
  * `epoch2::Int64`

**Returns**

Named tuple containing:

  * `ispc::Float64`: ISPC value
  * `ispc_angle::Float64`: ISPC angle
  * `signal_diff::Vector{Float64}`: signal difference (signal2 - signal1)
  * `phase_diff::Vector{Float64}`: phase difference (signal2 - signal1)
  * `s1_phase::Vector{Float64}`: signal 1 phase
  * `s2_phase::Vector{Float64}`: signal 2 phase

<a id='NeuroJ.eeg_itpc-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_itpc-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_itpc`** &mdash; *Method*.



```julia
eeg_itpc(eeg; channel)
```

Calculate ITPC (Inter-Trial-Phase Clustering) at time `t` over epochs/trials of `channel` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Int64`
  * `t::Int64`: time point (sample number) at which ITPC is calculated
  * `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc::Float64`: ITPC or wITPC value
  * `itpcz::Float64`: Rayleigh's ITPC Z value
  * `itpc_angle::Float64`: ITPC angle
  * `phase_diff::Array{Float64, 3}`: phase difference (channel2 - channel1)

<a id='NeuroJ.eeg_pli-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_pli-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_pli`** &mdash; *Method*.



```julia
eeg_pli(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculate PLI (Phase Lag Index) between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel1::Int64`
  * `channel2::Int64`
  * `epoch1::Int64`
  * `epoch2::Int64`

**Returns**

Named tuple containing:

  * `pli::Float64`: PLI value
  * `signal_diff::Vector{Float64}`: signal difference (signal2 - signal1)
  * `phase_diff::Vector{Float64}`: phase difference (signal2 - signal1)
  * `s1_phase::Vector{Float64}`: signal 1 phase
  * `s2_phase::Vector{Float64}`: signal 2 phase

<a id='NeuroJ.eeg_pli_m-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_pli_m-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_pli_m`** &mdash; *Method*.



```julia
eeg_pli_m(eeg; epoch)
```

Calculate matrix of PLIs (Phase Lag Index) between all channels of `eeg` at `epoch`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch1::Int64`

**Returns**

  * `pli_m::Matrix{Float64}`: PLI values matrix

<a id='NeuroJ.eeg_ispc_m-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_ispc_m-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_ispc_m`** &mdash; *Method*.



```julia
eeg_ispc_m(eeg; epoch)
```

Calculate matrix of ISPCs (Inter-Site-Phase Clustering) between all channels of `eeg` at `epoch`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch1::Int64`

**Returns**

  * `ispc_m::Matrix{Float64}`: ISPC values matrix

<a id='NeuroJ.eeg_aec-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_aec-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_aec`** &mdash; *Method*.



```julia
eeg_aec(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculate amplitude envelope correlation between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

**Arguments**

  * `eeg1::NeuroJ.EEG`
  * `eeg2::NeuroJ.EEG`
  * `channel1::Int64`
  * `channel2::Int64`
  * `epoch1::Int64`
  * `epoch2::Int64`

**Returns**

Named tuple containing:

  * `aec::Float64`: power correlation value
  * `aec_p::Float64`: power correlation p-value

<a id='NeuroJ.eeg_ged-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_ged-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_ged`** &mdash; *Method*.



```julia
eeg_ged(eeg1, eeg2)
```

Perform generalized eigendecomposition between `eeg1` and `eeg2`.

**Arguments**

  * `eeg1::NeuroJ.EEG`: signal data to be analyzed
  * `eeg2::NeuroJ.EEG`: original signal data

**Returns**

  * `sged::Array{Float64, 3}`
  * `ress::Matrix{Float64}`
  * `ress_normalized::Matrix{Float64}`

<a id='NeuroJ.eeg_frqinst-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_frqinst-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_frqinst`** &mdash; *Method*.



```julia
eeg_frqinst(eeg)
```

Calculate instantaneous frequency of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `frqinst::Array{Float64, 3}`

<a id='NeuroJ.eeg_itpc_s-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_itpc_s-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_itpc_s`** &mdash; *Method*.



```julia
eeg_itpc_s(eeg; <keyword arguments>)
```

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering) for `channel` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Int64`
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc_s::Array{Float64, 3}`: spectrogram of ITPC values
  * `itpc_z_s::Array{Float64, 3}`: spectrogram ITPCz values
  * `itpc_frq::Vector{Float64}`: frequencies list

<a id='NeuroJ.eeg_wspectrogram-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_wspectrogram-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_wspectrogram`** &mdash; *Method*.



```julia
eeg_wspectrogram(eeg; norm, mt, demean)
```

Return spectrogram of `eeg` using Morlet wavelet convolution.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `norm::Bool`=true: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `fs::Int64`: sampling rate
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet
  * `demean::Bool`=true: demean signal prior to analysis

**Returns**

Named tuple containing:

  * `w_pow::Array{Float64, 4}`
  * `w_frq::Matrix{Float64}`
  * `w_t::Matrix{Float64}`

<a id='NeuroJ.eeg_tkeo-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_tkeo-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_tkeo`** &mdash; *Method*.



```julia
eeg_tkeo(eeg)
```

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1)x(t+1)

**Arguments**

  * `eeg::NeuroJ.EEG`

**Returns**

  * `tkeo::Array{Float64, 3}`

<a id='NeuroJ.eeg_wspectrum-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_wspectrum-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_wspectrum`** &mdash; *Method*.



```julia
eeg_wspectrum(eeg; norm, mt, demean)
```

Return power spectrogrum of `eeg` using Morlet wavelet convolution.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `norm::Bool`=true: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `fs::Int64`: sampling rate
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet

**Returns**

Named tuple containing:

  * `w_pow::Array{Float64, 4}`
  * `w_frq::Matrix{Float64}`


<a id='EEG-plots'></a>

<a id='EEG-plots-1'></a>

## EEG plots

<a id='NeuroJ.plot_signal_scaled-Tuple{Union{AbstractRange, Vector{<:Real}}, AbstractArray}' href='#NeuroJ.plot_signal_scaled-Tuple{Union{AbstractRange, Vector{<:Real}}, AbstractArray}'>#</a>
**`NeuroJ.plot_signal_scaled`** &mdash; *Method*.



```julia
plot_signal_scaled(t, signal; <keyword arguments>)
```

Plot scaled multi-channel `signal`.

**Arguments**

  * `t::Union{Vector{<:Real}, AbstractRange}`
  * `signal::AbstractArray`
  * `labels::Vector{String}=[""]`: labels vector
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Channel"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_signal-Tuple{Union{AbstractRange, Vector{<:Real}}, Vector{<:Real}}' href='#NeuroJ.plot_signal-Tuple{Union{AbstractRange, Vector{<:Real}}, Vector{<:Real}}'>#</a>
**`NeuroJ.plot_signal`** &mdash; *Method*.



```julia
plot_signal(t, signal; <keyword arguments>)
```

Plot single-channel `signal`.

**Arguments**

  * `t::Union{Vector{<:Real}, AbstractRange}`
  * `signal::Vector{<:Real}`
  * `ylim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_signal-Tuple{Union{AbstractRange, Vector{<:Real}}, AbstractArray}' href='#NeuroJ.plot_signal-Tuple{Union{AbstractRange, Vector{<:Real}}, AbstractArray}'>#</a>
**`NeuroJ.plot_signal`** &mdash; *Method*.



```julia
plot_signal(t, signal; <keyword arguments>)
```

Plot multi-channel `signal`.

**Arguments**

  * `t::Union{Vector{<:Real}, AbstractRange}`
  * `signal::AbstractArray`
  * `labels::Vector{String}=[""]`: labels vector
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Channel"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
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
  * `scaled::Bool=false`: if true than scale signals before plotting so all signals will fit the plot
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

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
  * `mt::Bool=false`: if true use multi-tapered periodogram/spectrogram
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `head::Bool=true`: add head plot
  * `hist::Symbol=:hist`: histogram type: :hist, :kd
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for PSD and spectrogram
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `pc::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component`** &mdash; *Method*.



```julia
eeg_plot_component(eeg; <keyword arguments>)
```

Plot `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_idx-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_idx-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_idx`** &mdash; *Method*.



```julia
eeg_plot_component_idx(eeg; <keyword arguments>)
```

Plot indexed `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component index to display, default is all components
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_idx_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_idx_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_idx_avg`** &mdash; *Method*.



```julia
eeg_plot_component_idx_avg(eeg; <keyword arguments>)
```

Plot indexed `eeg` external or embedded component: mean and ±95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component index to display, default is all components
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_idx_butterfly-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_idx_butterfly-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_idx_butterfly`** &mdash; *Method*.



```julia
eeg_plot_component_idx_butterfly(eeg; <keyword arguments>)
```

Butterfly plot of indexed `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component index to display, default is all components
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_idx_psd-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_idx_psd-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_idx_psd`** &mdash; *Method*.



```julia
eeg_plot_component_idx_psd(eeg; <keyword arguments>)
```

Plot PSD of indexed `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `c_idx::Int64`: component index to display, default is all components
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_idx_psd_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_idx_psd_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_idx_psd_avg`** &mdash; *Method*.



```julia
eeg_plot_component_idx_psd_avg(eeg; <keyword arguments>)
```

Plot PSD of indexed `eeg` external or embedded component: mean ± 95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component index to display, default is all components
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_idx_psd_butterfly-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_idx_psd_butterfly-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_idx_psd_butterfly`** &mdash; *Method*.



```julia
eeg_plot_component_idx_psd_butterfly(eeg; <keyword arguments>)
```

Plot PSD of indexed `eeg` external or embedded component: mean ± 95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component index to display, default is all components
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_signal_avg-Tuple{Union{AbstractRange, Vector{<:Real}}, Matrix{Float64}}' href='#NeuroJ.plot_signal_avg-Tuple{Union{AbstractRange, Vector{<:Real}}, Matrix{Float64}}'>#</a>
**`NeuroJ.plot_signal_avg`** &mdash; *Method*.



```julia
plot_signal_avg(t, signal; <keyword arguments>)
```

Plot `signal` channels: mean and ±95% CI.

**Arguments**

  * `t::Union{Vector{<:Real}, AbstractRange}`
  * `signal::Matrix{<:Real}`
  * `norm::Bool=false`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_avg`** &mdash; *Method*.



```julia
eeg_plot_signal_avg(eeg; <keyword arguments>)
```

Plot `eeg` channels: mean and ±95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=0`: epoch number to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=false`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_avg_details-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_avg_details-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_avg_details`** &mdash; *Method*.



```julia
eeg_plot_avg_details(eeg; <keyword arguments>)
```

Plot details of averaged `eeg` channels: amplitude, histogram, power density, phase histogram and spectrogram.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=0`: epoch number to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=false`: normalize the `signal` prior to calculations
  * `mt::Bool=false`: if true use multi-tapered periodogram/spectrogram
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for PSD and spectrogram
  * `hist::Symbol=:hist`: histogram type: :hist, :kd
  * `head::Bool=true`: add head plot
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_avg`** &mdash; *Method*.



```julia
eeg_plot_component_avg(eeg; <keyword arguments>)
```

Plot `eeg` external or embedded component: mean and ±95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_signal_butterfly-Tuple{Union{AbstractRange, Vector{<:Real}}, Matrix{Float64}}' href='#NeuroJ.plot_signal_butterfly-Tuple{Union{AbstractRange, Vector{<:Real}}, Matrix{Float64}}'>#</a>
**`NeuroJ.plot_signal_butterfly`** &mdash; *Method*.



```julia
plot_signal_butterfly(t, signal; <keyword arguments>)
```

Butterfly plot of `signal` channels.

**Arguments**

  * `t::Union{Vector{<:Real}, AbstractRange}`
  * `signal::Matrix{<:Real}`
  * `labels::Vector{String}=[""]`: channel labels vector
  * `norm::Bool=false`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple`: y-axis limits, default (0, 0)
  * `mono::Bool=false`: use color or grey palette
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
  * `norm::Bool=false`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `mono::Bool=false`: use color or grey palette
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
  * `norm::Bool=false`: normalize the `signal` prior to calculations
  * `mt::Bool=false`: if true use multi-tapered periodogram/spectrogram
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for PSD and spectrogram
  * `hist::Symbol=:hist`: histogram type: :hist, :kd
  * `head::Bool=true`: add head plot
  * `mono::Bool=false`: use color or grey palette
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
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `norm::Bool=false`: normalize the `signal` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_psd-Tuple{Vector{<:Real}}' href='#NeuroJ.plot_psd-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.plot_psd`** &mdash; *Method*.



```julia
plot_psd(signal; <keyword arguments>)
```

Plot `signal` channel power spectrum density.

**Arguments**

  * `signal::Vector{<:Real}`
  * `fs::Int64`: sampling frequency
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered periodogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_psd_avg-Tuple{Matrix{Float64}}' href='#NeuroJ.plot_psd_avg-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.plot_psd_avg`** &mdash; *Method*.



```julia
plot_psd_avg(signal; <keyword arguments>)
```

Plot `signal` channels power spectrum density: mean and ±95% CI.

**Arguments**

  * `signal::Matrix{<:Real}`
  * `fs::Int64`: sampling rate
  * `norm::Bool=true`: normalize powers to dB
  * `mt::Bool=false`: if true use multi-tapered periodogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `labels::Vector{String}=[""]`: channel labels vector
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_psd_butterfly-Tuple{Matrix{Float64}}' href='#NeuroJ.plot_psd_butterfly-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.plot_psd_butterfly`** &mdash; *Method*.



```julia
plot_psd_butterfly(signal; <keyword arguments>)
```

Butterfly plot of `signal` channels power spectrum density.

**Arguments**

  * `signal::Matrix{<:Real}`
  * `fs::Int64`: sampling rate
  * `norm::Bool=true`: normalize powers to dB
  * `mt::Bool=false`: if true use multi-tapered periodogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `labels::Vector{String}=[""]`: channel labels vector
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_psd-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_psd-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_psd`** &mdash; *Method*.



```julia
eeg_plot_signal_psd(eeg; <keyword arguments>)
```

Plot `eeg` channels power spectrum density.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=0`: epoch number to display
  * `channel::Int64`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_psd_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_psd_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_psd_avg`** &mdash; *Method*.



```julia
eeg_plot_signal_psd_avg(eeg; <keyword arguments>)
```

Plot `eeg` channels power spectrum density: mean and ±95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=0`: epoch number to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `labels::Vector{String}=[""]`: channel labels vector
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_psd_butterfly-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_psd_butterfly-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_psd_butterfly`** &mdash; *Method*.



```julia
eeg_plot_signal_psd_butterfly(eeg; <keyword arguments>)
```

Plot `eeg` channels power spectrum density: mean and ±95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=0`: epoch number to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `labels::Vector{String}=[""]`: channel labels vector
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_psd-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_psd-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_psd`** &mdash; *Method*.



```julia
eeg_plot_component_psd(eeg; <keyword arguments>)
```

Plot PSD of `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `channel::Int64`: channel to display
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_psd_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_psd_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_psd_avg`** &mdash; *Method*.



```julia
eeg_plot_component_psd_avg(eeg; <keyword arguments>)
```

Plot PSD of `eeg` external or embedded component: mean and ±95% CI.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_psd_butterfly-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_psd_butterfly-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_psd_butterfly`** &mdash; *Method*.



```julia
eeg_plot_component_psd_butterfly(eeg; <keyword arguments>)
```

Butterfly plot PSD of `eeg` external or embedded component:.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channels to display, default is all channels
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_spectrogram-Tuple{Vector{<:Real}}' href='#NeuroJ.plot_spectrogram-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.plot_spectrogram`** &mdash; *Method*.



```julia
plot_spectrogram(signal; <keyword arguments>)
```

Plot spectrogram of `signal`.

**Arguments**

  * `signal::Vector{<:Real}`
  * `fs::Int64`: sampling frequency
  * `offset::Real`: displayed segment offset in seconds
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_spectrogram-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_spectrogram-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_spectrogram`** &mdash; *Method*.



```julia
eeg_plot_signal_spectrogram(eeg; <keyword arguments>)
```

Plots spectrogram of `eeg` channel(s).

**Arguments**

  * `eeg:EEG`
  * `epoch::Union{Int64, AbstractRange}=1`: epoch to plot
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_spectrogram_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_spectrogram_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_spectrogram_avg`** &mdash; *Method*.



```julia
eeg_plot_signal_spectrogram_avg(eeg; <keyword arguments>)
```

Plots spectrogram of `eeg` channel(s).

**Arguments**

  * `eeg:EEG`
  * `epoch::Union{Int64, AbstractRange}=1`: epoch to plot
  * `channel::Union{Vector{Int64}, AbstractRange}`: channels to plot
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_spectrogram-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_spectrogram-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_spectrogram`** &mdash; *Method*.



```julia
eeg_plot_component_spectrogram(eeg; <keyword arguments>)
```

Plots spectrogram of `eeg` external or embedded component.

**Arguments**

  * `eeg:EEG`
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `channel::Int64`: channel to display
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_spectrogram_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_spectrogram_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_spectrogram_avg`** &mdash; *Method*.



```julia
eeg_plot_component_spectrogram_avg(eeg; <keyword arguments>)
```

Plots spectrogram of `eeg` channel(s).

**Arguments**

  * `eeg:EEG`
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Union{Int64, AbstractRange}=1`: epoch to plot
  * `channel::Union{Vector{Int64}, AbstractRange}`: channels to plot
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_idx_spectrogram-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_idx_spectrogram-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_idx_spectrogram`** &mdash; *Method*.



```julia
eeg_plot_component_idx_spectrogram(eeg; <keyword arguments>)
```

Plot spectrogram of indexed `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `c_idx::Int64`: component index to display, default is all components
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Times [s]`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_component_idx_spectrogram_avg-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_component_idx_spectrogram_avg-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_component_idx_spectrogram_avg`** &mdash; *Method*.



```julia
eeg_plot_component_idx_spectrogram_avg(eeg; <keyword arguments>)
```

Plot spectrogram of averaged indexed `eeg` external or embedded component.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `c::Union{Array{Float64, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component index to display, default is all components
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
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
  * `mono::Bool=false`: use color or grey palette
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
  * `m::Union{Matrix{<:Real}, Array{Float64, 3}}`: channels by channels matrix
  * `epoch::Int64=1`: epoch number to display
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_covmatrix-Tuple{NeuroJ.EEG, Union{Array{Float64, 3}, Matrix{Float64}}, Vector{<:Real}}' href='#NeuroJ.eeg_plot_covmatrix-Tuple{NeuroJ.EEG, Union{Array{Float64, 3}, Matrix{Float64}}, Vector{<:Real}}'>#</a>
**`NeuroJ.eeg_plot_covmatrix`** &mdash; *Method*.



```julia
eeg_plot_matrix(eeg, cov_m, lags; <keyword arguments>)
```

Plot covariance matrix `m` of `eeg` channels.

**Arguments**

  * `eeg:EEG`
  * `cov_m::Union{Matrix{<:Real}, Array{Float64, 3}}`: covariance matrix
  * `lags::Vector{<:Real}`: covariance lags
  * `channel::Union{Int64, Vector{Int64}, AbstractRange, Nothing}`: channel to display
  * `epoch::Int64=1`: epoch number to display
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_histogram-Tuple{Vector{<:Real}}' href='#NeuroJ.plot_histogram-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.plot_histogram`** &mdash; *Method*.



```julia
plot_histogram(signal; <keyword arguments>)
```

Plot histogram of `signal`.

**Arguments**

  * `signal::Vector{<:Real}`
  * `type::Symbol`: type of histogram: regular `:hist` or kernel density `:kd`
  * `label::String=""`: channel label
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_histogram-Tuple{Matrix{Float64}}' href='#NeuroJ.plot_histogram-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.plot_histogram`** &mdash; *Method*.



```julia
plot_histogram(signal; <keyword arguments>)
```

Plot histogram of `signal`.

**Arguments**

  * `signal::Matrix{<:Real}`
  * `type::Symbol`: type of histogram: :hist or :kd
  * `labels::Vector{String}=[""]`
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
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
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_ica-Tuple{Union{AbstractRange, Vector{<:Real}}, Vector{Float64}}' href='#NeuroJ.plot_ica-Tuple{Union{AbstractRange, Vector{<:Real}}, Vector{Float64}}'>#</a>
**`NeuroJ.plot_ica`** &mdash; *Method*.



```julia
plot_ica(t, ica; <keyword arguments>)
```

Plot `ica` components against time vector `t`.

**Arguments**

  * `t::Union{Vector{<:Real}, AbstractRange}`: the time vector
  * `ica::Vector{Float64}`
  * `label::String=""`: channel label
  * `norm::Bool=true`: normalize the `ica` prior to calculations
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Amplitude [μV]"`: y-axis label
  * `title::String=""`: plot title
  * `ylim::Tuple{Real, Real}=(0, 0)`: y-axis limits (-ylim:ylim)
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_topo-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_topo-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_topo`** &mdash; *Method*.



```julia
eeg_plot_signal_topo(eeg; <keyword arguments>)
```

Plot topographical view of `eeg` signal.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Union{Int64, AbstractRange}=1`: epochs to display
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 second
  * `m::Symbol=:shepard`: interpolation method `:shepard` (Shepard), `:mq` (Multiquadratic), `:tp` (ThinPlate)
  * `cb::Bool=true`: draw color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_acomponent_topo-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_acomponent_topo-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_acomponent_topo`** &mdash; *Method*.



```julia
eeg_plot_acomponent_topo(eeg; <keyword arguments>)
```

Plot topographical view of `eeg` external or embedded component (array type: many values per channel per epoch).

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `c::Union{Array{<:Real, 3}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Int64`: epoch to display
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 second
  * `m::Symbol=:shepard`: interpolation method `:shepard` (Shepard), `:mq` (Multiquadratic), `:tp` (ThinPlate)
  * `cb_label::String="[A.U.]"`: color bar label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_weights_topo-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_weights_topo-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_weights_topo`** &mdash; *Method*.



```julia
eeg_plot_weights_topo(eeg; <keyword arguments>)
```

Topographical plot `eeg` of weights values at electrodes locations.

**Arguments**

  * `eeg:EEG`
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display, default is all channels
  * `weights=Matrix{<:Real}`: weights to plot
  * `head::Bool`=true: plot head
  * `small::Bool=false`: draws small plot
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_mcomponent_topo-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_mcomponent_topo-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_mcomponent_topo`** &mdash; *Method*.



```julia
eeg_plot_mcomponent_topo(eeg; <keyword arguments>)
```

Plot topographical view of `eeg` external or embedded component (matrix type: 1 value per channel per epoch).

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `epoch::Int64`: epoch to display
  * `c::Union{Matrix{<:Real}, Symbol}`: values to plot; if symbol, than use embedded component
  * `m::Symbol=:shepard`: interpolation method `:shepard` (Shepard), `:mq` (Multiquadratic), `:tp` (ThinPlate)
  * `cb::Bool=false`: draw color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_ica_topo-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_ica_topo-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_ica_topo`** &mdash; *Method*.



```julia
eeg_plot_ica_topo(eeg; <keyword arguments>)
```

Plot topographical view of `eeg` ICAs (each plot is signal reconstructed from this ICA).

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Int64`: epoch to display
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 second
  * `ic::Union{Vector{Int64}, AbstractRange}=0`: list of ICAs plot, default is all ICAs
  * `m::Symbol=:shepard`: interpolation method `:shepard` (Shepard), `:mq` (Multiquadratic), `:tp` (ThinPlate)
  * `cb::Bool=false`: draw color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_tile' href='#NeuroJ.eeg_plot_tile'>#</a>
**`NeuroJ.eeg_plot_tile`** &mdash; *Function*.



```julia
eeg_plot_tile(p)
```

Plot vector of plots `p` as tiles.

**Arguments**

  * `p::Vector{Any}`: vector of plots
  * `w::Int64=800`: single plot width (px)
  * `h::Int64=800`: single plot height (px)
  * `rows::Int64=2`: number of rows; if number of plots > 10 then number of rows = rows × 2
  * `mono::Bool=false`: use color or grey palette

**Returns**

  * `p_tiled::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_bands-Tuple{Vector{<:Real}}' href='#NeuroJ.plot_bands-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.plot_bands`** &mdash; *Method*.



```julia
plot_bands(signal; <keyword arguments>)
```

Plot absolute/relative bands powers of a single-channel `signal`.

**Arguments**

  * `signal::Vector{<:Real}`
  * `fs::Int64`: sampling rate
  * `band::Vector{Symbol}=[:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]`: band names, e.g. [:delta, alpha](see `eeg_band()`)
  * `band_frq::Vector{Tuple{Real, Real}}`: vector of band frequencies
  * `type::Symbol`: plots absolute (:abs) or relative power (:rel)
  * `norm::Bool=true`: normalize powers to dB
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
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
  * `norm::Bool=true`: normalize powers to dB
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
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

<a id='NeuroJ.eeg_plot_channels-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_channels-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_channels`** &mdash; *Method*.



```julia
eeg_plot_channels(eeg; <keyword arguments>)
```

Plot values of `c` for selected channels of `eeg`.

**Arguments**

  * `eeg:NeuroJ.EEG`
  * `c::Union{Matrix{Int64}, Matrix{<:Real}, Symbol}`: values to plot; if symbol, than use embedded component
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: list of channels to plot
  * `epoch::Int64`: number of epoch for which `c` should be plotted
  * `xlabel::String="Channel"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_epochs-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_epochs-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_epochs`** &mdash; *Method*.



```julia
eeg_plot_epochs(eeg; <keyword arguments>)
```

Plot values of `c` for selected epoch of `eeg`.

**Arguments**

  * `eeg:NeuroJ.EEG`
  * `c::Union{Vector{<:Real}, Symbol}`: values to plot; if symbol, than use embedded component
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: list of epochs to plot
  * `xlabel::String="Epochs"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
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
  * `fprototype::Symbol`: filter class: :fir, :butterworth, :chebyshev1, :chebyshev2, :elliptic
  * `ftype::Symbol`: filter type: :lp, :hp, :bp, :bs
  * `cutoff::Union{Real, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64`: filter order
  * `rp::Real`: dB ripple in the passband
  * `rs::Real`: dB attenuation in the stopband
  * `window::window::Union{Vector{Float64}, Nothing}`: window, required for FIR filter
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_compose-Tuple{Vector{Plots.Plot{Plots.GRBackend}}}' href='#NeuroJ.eeg_plot_compose-Tuple{Vector{Plots.Plot{Plots.GRBackend}}}'>#</a>
**`NeuroJ.eeg_plot_compose`** &mdash; *Method*.



```julia
eeg_plot_compose(p; <keyword arguments>)
```

Compose a complex plot of various plots contained in vector `p` using layout `layout`. Layout scheme is:

  * `(2, 2)`: 2 × 2 plots, regular layout
  * `@layout [a{0.2w} b{0.8w};_ c{0.6}]`: complex layout using Plots.jl `@layout` macro

**Arguments**

  * `p::Vector{Plots.Plot{Plots.GRBackend}}`: vector of plots
  * `layout::Union(Matrix{Any}, Tuple{Int64, Int64}}`: layout
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for `p` vector plots

**Returns**

  * `pc::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_env-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_env-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_env`** &mdash; *Method*.



```julia
eeg_plot_env(eeg; <keyword arguments>)
```

Plot envelope of `eeg` channels.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `type::Symbol`: envelope type: :amp (amplitude over time), :pow (power over frequencies), :spec (frequencies over time)
  * `average::Symbol`: averaging method: :no, :mean or :median
  * `dims::Union{Int64, Nothing}=nothing`: average over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `epoch::Int64`: epoch number to display
  * `channel::Int64`: channel to display
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `y_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for PSD and spectrogram
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_ispc-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_plot_ispc-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_ispc`** &mdash; *Method*.



```julia
eeg_plot_ispc(eeg1, eeg2; <keyword arguments>)
```

Plot ISPC `eeg1` and `eeg2` channels/epochs.

**Arguments**

  * `eeg1:NeuroJ.EEG`
  * `eeg2:NeuroJ.EEG`
  * `channel1::Int64`: epoch to plot
  * `channel2::Int64`: epoch to plot
  * `epoch1::Int64`: epoch to plot
  * `epoch2::Int64`: epoch to plot
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_itpc-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_itpc-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_itpc`** &mdash; *Method*.



```julia
eeg_plot_itpc(eeg; <keyword arguments>)
```

Plot ITPC (Inter-Trial-Phase Clustering) at time `t` over epochs/trials of `channel` of `eeg`.

**Arguments**

  * `eeg:NeuroJ.EEG`
  * `channel::Int64`: channel to plot
  * `t::Int64`: time point to plot
  * `z::Bool=false`: plot ITPCz instead of ITPC
  * `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_pli-Tuple{NeuroJ.EEG, NeuroJ.EEG}' href='#NeuroJ.eeg_plot_pli-Tuple{NeuroJ.EEG, NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_pli`** &mdash; *Method*.



```julia
eeg_plot_pli(eeg1, eeg2; <keyword arguments>)
```

Plot pli `eeg1` and `eeg2` channels/epochs.

**Arguments**

  * `eeg1:NeuroJ.EEG`
  * `eeg2:NeuroJ.EEG`
  * `channel1::Int64`: epoch to plot
  * `channel2::Int64`: epoch to plot
  * `epoch1::Int64`: epoch to plot
  * `epoch2::Int64`: epoch to plot
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_itpc_s-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_itpc_s-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_itpc_s`** &mdash; *Method*.



```julia
eeg_plot_itpc_s(eeg; <keyword arguments>)
```

Plot spectrogram of ITPC (Inter-Trial-Phase Clustering) for `channel` of `eeg`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Int64`
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:lin`: linear (:lin) or logarithmic (:log) frequencies
  * `z::Bool=false`: plot ITPCz instead of ITPC
  * `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String="ITPC spectrogram"`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_connections-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_connections-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_connections`** &mdash; *Method*.



```julia
eeg_plot_connections(eeg; <keyword arguments>)
```

Plot connections between `eeg` electrodes.

**Arguments**

  * `eeg:EEG`
  * `m::Matrix{<:Real}`: matrix of connections weights
  * `threshold::Float64`: plot all connection above threshold
  * `threshold_type::Symbol=:g`: rule for thresholding: :eq =, :geq ≥, :leq ≤, :g >, :l <
  * `labels::Bool=false`: plot electrode labels
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_itpc_f-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_itpc_f-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_itpc_f`** &mdash; *Method*.



```julia
eeg_plot_itpc_f(eeg; <keyword arguments>)
```

Plot time-frequency plot of ITPC (Inter-Trial-Phase Clustering) for `channel` of `eeg` for frequency `f`.

**Arguments**

  * `eeg::NeuroJ.EEG`
  * `channel::Int64`
  * `f::Int64`: frequency to plot
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:lin`: linear (:lin) or logarithmic (:log) frequencies
  * `z::Bool=false`: plot ITPCz instead of ITPC
  * `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation
  * `xlabel::String="Time [s]"`: x-axis label
  * `ylabel::String="Frequency [Hz]"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_psd_3dw-Tuple{Matrix{Float64}}' href='#NeuroJ.plot_psd_3dw-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.plot_psd_3dw`** &mdash; *Method*.



```julia
plot_psd_3dw(signal; <keyword arguments>)
```

Plot 3-d waterfall plot of `signal` channels power spectrum density.

**Arguments**

  * `signal::Matrix{Float64}`
  * `fs::Int64`: sampling frequency
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered periodogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel="Channel"`: y-axis label
  * `zlabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.plot_psd_3ds-Tuple{Matrix{Float64}}' href='#NeuroJ.plot_psd_3ds-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.plot_psd_3ds`** &mdash; *Method*.



```julia
plot_psd_3ds(signal; <keyword arguments>)
```

Plot 3-d surface plot of `signal` channels power spectrum density.

**Arguments**

  * `signal::Matrix{Float64}`
  * `fs::Int64`: sampling frequency
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered periodogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel="Channel"`: y-axis label
  * `zlabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroJ.eeg_plot_signal_psd_3d-Tuple{NeuroJ.EEG}' href='#NeuroJ.eeg_plot_signal_psd_3d-Tuple{NeuroJ.EEG}'>#</a>
**`NeuroJ.eeg_plot_signal_psd_3d`** &mdash; *Method*.



```julia
eeg_plot_signal_psd_3d(eeg; <keyword arguments>)
```

Plot 3-d waterfall plot of `eeg` channels power spectrum density.

**Arguments**

  * `eeg::NeuroJ.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=0`: epoch number to display
  * `channel::Int64`: channel to display, default is all channels
  * `offset::Int64=0`: displayed segment offset in samples
  * `len::Int64=0`: displayed segment length in samples, default is 1 epoch or 20 seconds
  * `type::Symbol=:w`: plot type: :w waterfall, :s surface
  * `norm::Bool=true`: normalize powers to dB
  * `mw::Bool=false`: if true use Morlet wavelet convolution
  * `mt::Bool=false`: if true use multi-tapered periodogram
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `xlabel::String="Frequency [Hz]`: x-axis label
  * `ylabel="Channel"`: y-axis label
  * `zlabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


<a id='Low-level-functions'></a>

<a id='Low-level-functions-1'></a>

## Low level functions

<a id='NeuroJ.linspace-Tuple{Real, Real, Int64}' href='#NeuroJ.linspace-Tuple{Real, Real, Int64}'>#</a>
**`NeuroJ.linspace`** &mdash; *Method*.



```julia
linspace(start, stop, length)
```

Generates `length`-long sequence of evenly spaced numbers between `start` and `stop`.

**Arguments**

  * `start::Real`
  * `stop::Real`
  * `length::Int64`

**Returns**

  * `range::Number`

<a id='NeuroJ.logspace-Tuple{Real, Real, Int64}' href='#NeuroJ.logspace-Tuple{Real, Real, Int64}'>#</a>
**`NeuroJ.logspace`** &mdash; *Method*.



```julia
logspace(start, stop, length)
```

Generates `length`-long sequence of log10-spaced numbers between `start` and `stop`.

**Arguments**

  * `start::Real`
  * `stop::Real`
  * `length::Int64`

**Returns**

  * `range::Number`

<a id='NeuroJ.m_pad0-Tuple{Matrix{<:Number}}' href='#NeuroJ.m_pad0-Tuple{Matrix{<:Number}}'>#</a>
**`NeuroJ.m_pad0`** &mdash; *Method*.



```julia
m_pad0(m)
```

Pad the matrix `m` with zeros to make it square.

**Arguments**

  * `m::Matrix{<:Number}`

**Returns**

  * `m::Matrix{Number}`

<a id='NeuroJ.vsearch-Tuple{Real, Vector{<:Real}}' href='#NeuroJ.vsearch-Tuple{Real, Vector{<:Real}}'>#</a>
**`NeuroJ.vsearch`** &mdash; *Method*.



```julia
vsearch(y, x; return_distance=false)
```

Return the positions of the `y` value in the vector `x` and the difference between `y` and `x[vsearch(x, y)].

**Arguments**

  * `y::Real`
  * `x::Vector{<:Real}`
  * `return_distance::Bool`

**Returns**

  * `y_idx::Int64`

-`y_dist::Real`

<a id='NeuroJ.vsearch-Tuple{Vector{<:Real}, Vector{<:Real}}' href='#NeuroJ.vsearch-Tuple{Vector{<:Real}, Vector{<:Real}}'>#</a>
**`NeuroJ.vsearch`** &mdash; *Method*.



```julia
vsearch(y, x; return_distance=false)
```

Return the positions of the `y` vector in the vector `x`.

**Arguments**

  * `x::Vector{<:Real}`
  * `y::Vector{<:Real}`
  * `return_distance::Bool`

**Returns**

  * `y_idx::Int64`
  * `y_dist::Real`

<a id='NeuroJ.cart2pol-Tuple{Real, Real}' href='#NeuroJ.cart2pol-Tuple{Real, Real}'>#</a>
**`NeuroJ.cart2pol`** &mdash; *Method*.



```julia
cart2pol(x, y)
```

Convert cartographic coordinates `x` and `y` to polar.

**Arguments**

  * `x::Real`
  * `y::Real`

**Returns**

  * `phi::Float64`
  * `theta::Float64`

<a id='NeuroJ.pol2cart-Tuple{Real, Real}' href='#NeuroJ.pol2cart-Tuple{Real, Real}'>#</a>
**`NeuroJ.pol2cart`** &mdash; *Method*.



```julia
pol2cart(theta, phi)
```

Convert polar coordinates `theta` and `phi` to cartographic.

**Arguments**

  * `phi::Real`
  * `theta::Real`

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

  * `phi::Real`: the angle with respect to the z-axis (elevation)
  * `theta::Real`: the angle in the xy plane with respect to the x-axis (azimuth)
  * `rho::Real`: the distance from the origin to the point

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
      * `:triangle`: symmetric triangle (left half ↑, right half ↓)
      * `:exp`: symmetric exponential (left half ↑, right half ↓)
  * `n::Int64`: window length
  * `even::Bool=false`: if true, make the window of even length (+1 for odd n)

**Returns**

  * `w::Vector{Float64}`:: generated window

<a id='NeuroJ.hildebrand_rule-Tuple{Vector{<:Real}}' href='#NeuroJ.hildebrand_rule-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.hildebrand_rule`** &mdash; *Method*.



```julia
hildebrand_rule(x)
```

Calculate Hildebrand rule for vector `x`. If H < 0.2 then the vector `x` is symmetrical.

**Arguments**

  * `x::Vector{<:Real}`

**Returns**

  * `h::Float64`

<a id='NeuroJ.jaccard_similarity-Tuple{Vector{<:Real}, Vector{<:Real}}' href='#NeuroJ.jaccard_similarity-Tuple{Vector{<:Real}, Vector{<:Real}}'>#</a>
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

  * `x::Vector{<:Real}`
  * `n::Int64`

**Returns**

  * `x::Vector{Vector{<:Real}}`

<a id='NeuroJ.s_rms-Tuple{Vector{<:Real}}' href='#NeuroJ.s_rms-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.s_rms`** &mdash; *Method*.



```julia
s_rms(signal)
```

Calculate Root Mean Square of `signal`.

**Arguments**

  * `signal::Vector{<:Real}`

**Returns**

  * rms::Float64`

<a id='NeuroJ.generate_sine' href='#NeuroJ.generate_sine'>#</a>
**`NeuroJ.generate_sine`** &mdash; *Function*.



```julia
generate_sine(f, t, a, p)
```

Generates sine wave of `f` frequency over `t` time; optional arguments are: `a` amplitude and  `p` phase.

**Arguments**

  * `f::Real`: frequency
  * `t::Union{Vector{<:Real}, AbstractRange}`: time vector
  * `a::Real`: amplitude
  * `p::Real`: initial phase

**Returns**

  * sine::Vector{Float64}`

<a id='NeuroJ.s_freqs-Tuple{Union{AbstractRange, Vector{<:Real}}}' href='#NeuroJ.s_freqs-Tuple{Union{AbstractRange, Vector{<:Real}}}'>#</a>
**`NeuroJ.s_freqs`** &mdash; *Method*.



```julia
s_freqs(t)
```

Return vector of frequencies and Nyquist frequency for given time vector `t`.

**Arguments**

  * `t::Vector{<:Real}, AbstractRange}`

**Returns**

  * `hz::Vector{Float64}`
  * `nyquist_freq::Float64`

<a id='NeuroJ.s_freqs-Tuple{Vector{Float64}, Real}' href='#NeuroJ.s_freqs-Tuple{Vector{Float64}, Real}'>#</a>
**`NeuroJ.s_freqs`** &mdash; *Method*.



```julia
s_freqs(signal, fs)
```

Return vector of frequencies and Nyquist frequency for given `signal` and `fs`.

**Arguments**

  * `signal::Vector{Float64}`
  * `fs::Real`

**Returns**

  * `hz::Vector{Float64`
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

  * `x::Vector{<:Real}`
  * `n::Int64`
  * `sym::Bool=false`: if true, than pad at the beginning and at the end, otherwise only at the end.

**Returns**

  * `v_pad::Vector{<:Real}`

<a id='NeuroJ.hz2rads-Tuple{Real}' href='#NeuroJ.hz2rads-Tuple{Real}'>#</a>
**`NeuroJ.hz2rads`** &mdash; *Method*.



```julia
hz2rads(f)
```

Convert frequency `f` in Hz to rad/s.

**Arguments**

  * `f::Real`

**Returns**

  * `f_rads::Float64`

<a id='NeuroJ.rads2hz-Tuple{Real}' href='#NeuroJ.rads2hz-Tuple{Real}'>#</a>
**`NeuroJ.rads2hz`** &mdash; *Method*.



```julia
rads2hz(f)
```

Convert frequency `f` in rad/s to Hz.

**Arguments**

  * `f::Real`

**Returns**

  * `f_rads::Float64`

<a id='NeuroJ.z_score-Tuple{Vector{<:Real}}' href='#NeuroJ.z_score-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.z_score`** &mdash; *Method*.



```julia
z_score(x)
```

Calculate Z-scores for each value of the vector `x`.

**Arguments**

  * `x::Vector{<:Real}`

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
  * `f::Real=10.0`: frequency
  * `peak::Real=0`: sinc peak time
  * `norm::Bool=true`: generate normalized function

**Returns**

  * `sinc::Vector{Float64}`

<a id='NeuroJ.generate_morlet' href='#NeuroJ.generate_morlet'>#</a>
**`NeuroJ.generate_morlet`** &mdash; *Function*.



```julia
generate_morlet(fs, f, t; ncyc, complex)
```

Generate Morlet wavelet.

**Arguments**

  * `fs::Int64`: sampling rate
  * `f::Real`: frequency
  * `t::Real=1`: length = -t:1/fs:t
  * `ncyc::Int64=5`: number of cycles
  * `complex::Bool=false`: generate complex Morlet

**Returns**

  * `morlet::Union{Vector{Float64}, Vector{ComplexF64}}`

<a id='NeuroJ.generate_gaussian' href='#NeuroJ.generate_gaussian'>#</a>
**`NeuroJ.generate_gaussian`** &mdash; *Function*.



```julia
generate_gaussian(fs, f, t; ncyc, a)
```

Generate Gaussian wave.

**Arguments**

  * `fs::Int64`: sampling rate
  * `f::Real`: frequency
  * `t::Real=1`: length = -t:1/fs:t
  * `ncyc::Int64`: : number of cycles, width, SD of the Gaussian
  * `a::Real=1`: peak amp

**Returns**

  * `gaussian::Vector{Float64}`

<a id='NeuroJ.tuple_order' href='#NeuroJ.tuple_order'>#</a>
**`NeuroJ.tuple_order`** &mdash; *Function*.



```julia
tuple_order(t, rev)
```

Order tuple elements in ascending or descending (rev=true) order.

**Arguments**

  * `t::Tuple{Real, Real}`
  * `rev::Bool=false`

**Returns**

  * `t::Tuple{Real, Real}`

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

<a id='NeuroJ.s2_cov-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s2_cov-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s2_cov`** &mdash; *Method*.



s2_cov(signal1, signal2; norm=true)

Calculates covariance between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`
  * `norm::Bool`: normalize covariance

**Returns**

  * `cov_mat::Matrix{Float64}`

<a id='NeuroJ.s_dft-Tuple{AbstractArray}' href='#NeuroJ.s_dft-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_dft`** &mdash; *Method*.



```julia
s_dft(signal; fs)
```

Return FFT and DFT sample frequencies for a DFT for the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate

**Returns**

Named tuple containing:

  * `s_fft::Vector{ComplexF64}`
  * `s_sf::Vector{Float64}`

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
s_spectrum(signal; pad)
```

Calculates FFT, amplitudes, powers and phases of the `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `pad::Int64=0`: pad the `signal` with `pad` zeros

**Returns**

Named tuple containing:

  * `s_fft::Vector{ComplexF64}`
  * `s_amplitudes::Vector{Float64}`
  * `s_powers::Vector{Float64}`
  * `s_phases::Vector{Float64}`

<a id='NeuroJ.s_total_power-Tuple{AbstractArray}' href='#NeuroJ.s_total_power-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_total_power`** &mdash; *Method*.



```julia
s_total_power(signal; fs)
```

Calculates `signal` total power.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate
  * `mt::Bool=false`: if true use multi-tapered periodogram

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
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds

**Returns**

  * `sbp::Float64`: signal band power

<a id='NeuroJ.s_taper-Tuple{AbstractArray}' href='#NeuroJ.s_taper-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_taper`** &mdash; *Method*.



```julia
s_taper(signal; taper)
```

Taper the `signal` with `taper`.

**Arguments**

  * `signal::AbstractArray`
  * `taper::Union{Vector{<:Real}, Vector{ComplexF64}}`

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
  * `offset::Real=0`: constant for :constant detrending
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
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`

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
      * `:iirnotch`
      * `:remez`
      * `:mavg`: moving average (with threshold)
      * `:mmed`: moving median (with threshold)
      * `:poly`: polynomial of `order` order
  * `ftype::Symbol`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64=8`: filter order, number of taps for :remez filter, k-value for :mavg and :mmed (window length = 2 × k + 1)
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
  * `bw::Real=-1`: bandwidth for :iirnotch and :remez filters
  * `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass), for causal filter use :onepass
  * `t::Real`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
  * `window::Union{Vector{<:Real}, Nothing} - window, required for FIR filter, weighting window for :mavg and :mmed

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
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `psd_pow::Vector{Float64}`
  * `psd_frq::Vector{Float64}`

<a id='NeuroJ.s_psd-Tuple{Matrix{Float64}}' href='#NeuroJ.s_psd-Tuple{Matrix{Float64}}'>#</a>
**`NeuroJ.s_psd`** &mdash; *Method*.



```julia
s_psd(signal; fs, norm=false)
```

Calculate power spectrum density of the `signal`.

**Arguments**

  * `signal::Matrix{Float64}`
  * `fs::Int64`: sampling rate
  * `norm::Bool`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

named tuple containing:

  * `psd_pow::Matrix{Float64}`
  * `psd_frq::Matrix{Float64}`

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
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

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

<a id='NeuroJ.s_negentropy-Tuple{AbstractArray}' href='#NeuroJ.s_negentropy-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_negentropy`** &mdash; *Method*.



```julia
s_negentropy(signal)
```

Calculate negentropy of `signal`.

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

<a id='NeuroJ.s2_tcoherence-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s2_tcoherence-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s2_tcoherence`** &mdash; *Method*.



```julia
s2_tcoherence(signal1, signal2)
```

Calculate coherence (mean over time), IC (imaginary coherence) and MSC (magnitude-squared coherence) between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

Named tuple containing:

  * `c::Vector{Float64}`: coherence
  * `msc::Vector{Float64}`: magnitude-squares coherence
  * `ic::Vector{Float64}`: imaginary part of coherence

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
  * `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
  * `pc_m::PCA{Float64}`: PC mean

<a id='NeuroJ.s_pca_reconstruct-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_pca_reconstruct-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_pca_reconstruct`** &mdash; *Method*.



```julia
s_pca_reconstruct(signal, pc, pcm)
```

Reconstructs `signal` using PCA components.

**Arguments**

  * `signal::Array{Float64, 3}`
  * `pc::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
  * `pc_m::PCA{Float64}:`: IC(1)..IC(n) × epoch

**Returns**

  * `s_reconstructed::Array{Float64, 3}`

<a id='NeuroJ.s_fconv-Tuple{AbstractArray}' href='#NeuroJ.s_fconv-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_fconv`** &mdash; *Method*.



```julia
s_fconv(signal; kernel, norm)
```

Perform convolution in the frequency domain between `signal` and `kernel`.

**Arguments**

  * `signal::AbstractArray`
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`
  * `norm::Bool=false`: normalize kernel

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
  * `norm::Bool=true`: normalize powers to dB
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `demean::Bool=true`: demean signal prior to analysis

**Returns**

Named tuple containing:

  * `s_pow::Matrix{Float64}`: powers
  * `s_frq::Vector{Float64}`: frequencies
  * `s_t::Vector{Float64}`: time

<a id='NeuroJ.s_detect_epoch_flat-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_flat-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_flat`** &mdash; *Method*.



```julia
s_detect_epoch_flat(signal)
```

Detect bad `signal` epochs based on: flat channel(s)

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.s_detect_epoch_rmse-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_rmse-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_rmse`** &mdash; *Method*.



```julia
s_detect_epoch_rmse(signal)
```

Detect bad `signal` epochs based on: RMSE vs average channel > 95%CI.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.s_detect_epoch_rmsd-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_rmsd-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_rmsd`** &mdash; *Method*.



```julia
detect_epoch_rmsd(signal)
```

Detect bad `signal` epochs based on: RMSD vs average channel > 95%CI.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.s_detect_epoch_euclid-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_euclid-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_euclid`** &mdash; *Method*.



```julia
s_detect_epoch_euclid(signal)
```

Detect bad `signal` epochs based on: Euclidean distance vs median channel > 95% CI.

**Arguments**

  * `signal::Array{Float64, 3}`

**Returns**

  * `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch

<a id='NeuroJ.s_detect_epoch_p2p-Tuple{Array{Float64, 3}}' href='#NeuroJ.s_detect_epoch_p2p-Tuple{Array{Float64, 3}}'>#</a>
**`NeuroJ.s_detect_epoch_p2p`** &mdash; *Method*.



```julia
s_detect_epoch_p2p(signal)
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

<a id='NeuroJ.s_findpeaks-Tuple{AbstractArray}' href='#NeuroJ.s_findpeaks-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_findpeaks`** &mdash; *Method*.



```julia
s_findpeaks(signal; d)
```

Find peaks in `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `d::Int64=32`: distance between peeks in samples

**Returns**

  * `p_idx::Vector{Int64}`

<a id='NeuroJ.s_wdenoise-Tuple{AbstractArray}' href='#NeuroJ.s_wdenoise-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_wdenoise`** &mdash; *Method*.



```julia
s_wdenoise(signal; wt)
```

Perform wavelet denoising.

**Arguments**

  * `signal::AbstractArray`
  * `wt::Symbol=:db4`: wavelet type: :db2, :db4, :db8, :db10, :haar, :coif2, :coif4, :coif8

**Returns**

  * `signal_denoised::Vector{Float64}`

<a id='NeuroJ.effsize-Tuple{Vector{<:Real}, Vector{<:Real}}' href='#NeuroJ.effsize-Tuple{Vector{<:Real}, Vector{<:Real}}'>#</a>
**`NeuroJ.effsize`** &mdash; *Method*.



```julia
effsize(x1, x2)
```

Calculate Cohen's d and Hedges g effect sizes.

**Arguments**

  * `x1::Vector{Float64}`
  * `x2::Vector{Float64}`

**Returns**

Named tuple containing:

  * `d::Float64`: Cohen's d
  * `g::Float64`: Hedges g

<a id='NeuroJ.s_ispc-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s_ispc-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s_ispc`** &mdash; *Method*.



```julia
s_ispc(signal1, signal2)
```

Calculate ISPC (Inter-Site-Phase Clustering) between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

Named tuple containing:

  * `ispc::Float64`: ISPC value
  * `ispc_angle::Float64`: ISPC angle
  * `signal_diff::Vector{Float64}`: signal difference (signal2 - signal1)
  * `phase_diff::Vector{Float64}`: phase difference (signal2 - signal1)
  * `s1_phase::Vector{Float64}`: signal 1 phase
  * `s2_phase::Vector{Float64}`: signal 2 phase

<a id='NeuroJ.s_itpc-Tuple{AbstractArray}' href='#NeuroJ.s_itpc-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_itpc`** &mdash; *Method*.



```julia
s_itpc(signal; t)
```

Calculate ITPC (Inter-Trial-Phase Clustering) over epochs/trials at time `t` of `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `t::Int64`: time point (sample number) at which ITPC is calculated
  * `w::Union{Vector{<:Real}, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc::Float64`: ITPC value
  * `itpcz::Float64`: Rayleigh's ITPC Z value
  * `itpc_angle::Float64`: ITPC angle
  * `itpc_phases::Vector{Float64}`: phases at time `t` averaged across trials/epochs

<a id='NeuroJ.s_pli-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s_pli-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s_pli`** &mdash; *Method*.



```julia
s_pli(signal1, signal2)
```

Calculate PLI (Phase-Lag Index) between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

Named tuple containing:

  * `pli::Float64`: PLI value
  * `signal_diff::Vector{Float64}`: signal difference (signal2 - signal1)
  * `phase_diff::Vector{Float64}`: phase difference (signal2 - signal1)
  * `s1_phase::Vector{Float64}`: signal 1 phase
  * `s2_phase::Vector{Float64}`: signal 2 phase

<a id='NeuroJ.s_ged-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s_ged-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s_ged`** &mdash; *Method*.



```julia
s_ged(signal1, signal2)
```

Perform generalized eigendecomposition between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`: signal to be analyzed
  * `signal2::AbstractArray`: original signal

**Returns**

Named tuple containing:

  * `sged::AbstractArray`
  * `ress::AbstractArray`
  * `ress_normalized::AbstractArray`: RESS normalized to -1..1

<a id='NeuroJ.s_frqinst-Tuple{AbstractArray}' href='#NeuroJ.s_frqinst-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_frqinst`** &mdash; *Method*.



```julia
s_frqinst(signal; fs)
```

Calculate instantaneous frequency `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`

**Returns**

  * `frqinst::Vector{Float64}`

<a id='NeuroJ.s_hspectrum-Tuple{AbstractArray}' href='#NeuroJ.s_hspectrum-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_hspectrum`** &mdash; *Method*.



```julia
s_hspectrum(signal; pad=0)
```

Calculate amplitudes, powers and phases of the `signal` using Hilbert transform.

**Arguments**

  * `signal::AbstractArray`
  * `pad::Int64`: pad the `signal` with `pad` zeros

**Returns**

Named tuple containing:

  * `h::Vector(ComplexF64}`: Hilbert components
  * `h_amplitudes::Vector{Float64}`
  * `h_powers::Vector{Float64}`
  * `h_phases::Vector{Float64}`

<a id='NeuroJ.t2f-Tuple{Real}' href='#NeuroJ.t2f-Tuple{Real}'>#</a>
**`NeuroJ.t2f`** &mdash; *Method*.



```julia
t2f(t)
```

Convert cycle length in ms `t` to frequency.

**Arguments**

  * `t::Real`: cycle length in ms

**Returns**

  * `f::Float64`: frequency in Hz

<a id='NeuroJ.f2t-Tuple{Real}' href='#NeuroJ.f2t-Tuple{Real}'>#</a>
**`NeuroJ.f2t`** &mdash; *Method*.



```julia
f2t(f)
```

Convert frequency `f` to cycle length in ms.

**Arguments**

  * `f::Real`: frequency in Hz

**Returns**

  * `f::Float64`: cycle length in ms

<a id='NeuroJ.s_wspectrogram-Tuple{AbstractArray}' href='#NeuroJ.s_wspectrogram-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_wspectrogram`** &mdash; *Method*.



```julia
s_wspectrogram(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc, demean)
```

Calculate spectrogram of the `signal` using wavelet convolution.

**Arguments**

  * `signal::AbstractArray`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `fs::Int64`: sampling rate
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet
  * `demean::Bool=true`: demean signal prior to analysis

**Returns**

Named tuple containing:

  * `w_conv::Matrix(ComplexF64}`: convoluted signal
  * `w_powers::Matrix{Float64}`
  * `w_phases::Matrix{Float64}`
  * `frq_list::Vector{Float64}`

<a id='NeuroJ.s_fftdenoise-Tuple{AbstractArray}' href='#NeuroJ.s_fftdenoise-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_fftdenoise`** &mdash; *Method*.



s_fftdenoise(signal::AbstractArray; pad::Int64=0, threshold::Int64=100) 

Perform FFT denoising.

**Arguments**

  * `signal::AbstractArray`
  * `pad::Int64=0`: pad the `signal` with `pad` zeros
  * `threshold::Int64=100`: PSD threshold for keeping frequency components

**Returns**

  * `signal_denoised::Vector{Float64}`

<a id='NeuroJ.s_gfilter-Tuple{Vector{Float64}}' href='#NeuroJ.s_gfilter-Tuple{Vector{Float64}}'>#</a>
**`NeuroJ.s_gfilter`** &mdash; *Method*.



```julia
s_gfilter(signal, fs, f, gw)
```

Filter `signal` using Gaussian in the frequency domain.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate
  * `f::Real`: filter frequency
  * `gw::Real=5`: Gaussian width in Hz

**Returns**

Named tuple containing:

  * `s_f::Vector{Float64}`

<a id='NeuroJ.s_ghspectrogram-Tuple{AbstractArray}' href='#NeuroJ.s_ghspectrogram-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_ghspectrogram`** &mdash; *Method*.



```julia
s_ghspectrogram(signal; pad, norm, frq_lim, frq_n, frq, fs, gw, demean)
```

Calculate spectrogram of the `signal` using Gaussian and Hilbert transform.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `gw::Real=5`: Gaussian width in Hz
  * `demean::Bool`=true: demean signal prior to analysis

**Returns**

Named tuple containing:

  * `h_powers::Matrix{Float64}`
  * `frq_list::Vector{Float64}`

<a id='NeuroJ.s_tkeo-Tuple{AbstractArray}' href='#NeuroJ.s_tkeo-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_tkeo`** &mdash; *Method*.



```julia
s_tkeo(signal)
```

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1)x(t+1)

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_new::Vector{Float64}`

<a id='NeuroJ.s_wspectrum-Tuple{AbstractArray}' href='#NeuroJ.s_wspectrum-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_wspectrum`** &mdash; *Method*.



```julia
s_wspectrum(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc)
```

Calculate power spectrum of the `signal` using wavelet convolution.

**Arguments**

  * `signal::AbstractArray`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `fs::Int64`: sampling rate
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet

**Returns**

Named tuple containing:

  * `w_powers::Matrix{Float64}`
  * `frq_list::Vector{Float64}`

<a id='NeuroJ.a2_cmp-Tuple{Array{<:Real, 3}, Array{<:Real, 3}}' href='#NeuroJ.a2_cmp-Tuple{Array{<:Real, 3}, Array{<:Real, 3}}'>#</a>
**`NeuroJ.a2_cmp`** &mdash; *Method*.



```julia
a2_cmp(a1, a2; p, perm_n)
```

Compare two 3-dimensional arrays `a1` and `a2` (e.g. two spectrograms), using permutation based statistic.

**Arguments**

  * `a1::Array{<:Real, 3}`: first array
  * `a2::Array{<:Real, 3}`: second array
  * `p::Float64=0.05`: p-value
  * `perm_n::Int64=1000`: number of permutations

**Returns**

Named tuple containing:

  * `zmap::Array{Float64, 3}`: array of Z-values
  * `zmap_b::Array{Float64, 3}`: binarized mask of statistically significant positions

<a id='NeuroJ.s_fcoherence-Tuple{AbstractArray}' href='#NeuroJ.s_fcoherence-Tuple{AbstractArray}'>#</a>
**`NeuroJ.s_fcoherence`** &mdash; *Method*.



```julia
s_fcoherence(signal; fs, frq)
```

Calculate coherence (mean over all frequencies) and MSC (magnitude-squared coherence) between channels of `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`
  * `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

**Returns**

  * `c::Array{Float64, 3}`: coherence
  * `msc::Array{Float64, 3}`: MSC
  * `f::Vector{Float64}`: frequencies

<a id='NeuroJ.s2_fcoherence-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.s2_fcoherence-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.s2_fcoherence`** &mdash; *Method*.



```julia
s2_fcoherence(signal1, signal2; fs, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)
```

Calculate coherence (mean over all frequencies) and MSC (magnitude-squared coherence) between channels of `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`
  * `fs::Int64`
  * `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

**Returns**

  * `c::Array{Float64, 3}`: coherence
  * `msc::Array{Float64, 3}`: MSC
  * `f::Vector{Float64}`: frequencies

<a id='NeuroJ.a2_l1-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.a2_l1-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.a2_l1`** &mdash; *Method*.



```julia
a2_l1(a1, a2)
```

Compare two 3-dimensional arrays `a1` and `a2` (e.g. two spectrograms), using L1 (Manhattan) distance.

**Arguments**

  * `a1::AbstractArray`: first array
  * `a2::AbstractArray`: second array

**Returns**

  * `l1::Float64`

<a id='NeuroJ.a2_l2-Tuple{AbstractArray, AbstractArray}' href='#NeuroJ.a2_l2-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroJ.a2_l2`** &mdash; *Method*.



```julia
a2_l2(a1, a2)
```

Compare two 3-dimensional arrays `a1` and `a2` (e.g. two spectrograms), using L2 (Euclidean) distance.

**Arguments**

  * `a1::AbstractArray`: first array
  * `a2::AbstractArray`: second array

**Returns**

  * `l2::Float64`

<a id='NeuroJ.s_cums-Tuple{Vector{<:Real}}' href='#NeuroJ.s_cums-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.s_cums`** &mdash; *Method*.



```julia
s_cums(signal)
```

Calculate cumulative sum of the `signal`.

**Arguments**

  * `signal::Vector{<:Real}`

**Returns**

  * `signal_cs::Vector{Float64}`

<a id='NeuroJ.s_cums-Tuple{Array{<:Real, 3}}' href='#NeuroJ.s_cums-Tuple{Array{<:Real, 3}}'>#</a>
**`NeuroJ.s_cums`** &mdash; *Method*.



```julia
s_cums(signal)
```

Calculate cumulative sum of the `signal`.

**Arguments**

  * `signal::Array{<:Real, 3}`

**Returns**

  * `signal_cs::Array{Float64, 3}`

<a id='NeuroJ.s_gfp-Tuple{Vector{<:Real}}' href='#NeuroJ.s_gfp-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.s_gfp`** &mdash; *Method*.



```julia
s_gfp(signal)
```

Calculate GFP (Global Field Power) of the `signal`.

**Arguments**

  * `signal::Vector{<:Real}`

**Returns**

  * `gfp::Float64`

<a id='NeuroJ.s_gfp_norm-Tuple{Vector{<:Real}}' href='#NeuroJ.s_gfp_norm-Tuple{Vector{<:Real}}'>#</a>
**`NeuroJ.s_gfp_norm`** &mdash; *Method*.



```julia
s_gfp_norm(signal)
```

Calculate `signal` values normalized for GFP (Global Field Power) of that signal.

**Arguments**

  * `signal::Vector{<:Real}`

**Returns**

  * `gfp_norm::Float64`

<a id='NeuroJ.s2_diss-Tuple{Vector{<:Real}, Vector{<:Real}}' href='#NeuroJ.s2_diss-Tuple{Vector{<:Real}, Vector{<:Real}}'>#</a>
**`NeuroJ.s2_diss`** &mdash; *Method*.



```julia
s2_diss(signal1, signal2)
```

Calculate DISS (global dissimilarity) and spatial correlation between `signal1` and `signal2`.

**Arguments**

  * `signal1::Vector{<:Real}`
  * `signal2::Vector{<:Real}`

**Returns**

Named tuple containing:

  * `diss::Float64`: global dissimilarity
  * `c::Float64`: spatial correlation


<a id='NSTIM'></a>

<a id='NSTIM-1'></a>

## NSTIM

<a id='NeuroJ.tes_dose-Tuple{Real, Real, Int64}' href='#NeuroJ.tes_dose-Tuple{Real, Real, Int64}'>#</a>
**`NeuroJ.tes_dose`** &mdash; *Method*.



```julia
tes_dose(current, pad_area, duration)
```

Converts `current`, `pad_area` and stimulation `duration` into `charge`, `current_density` and `charge_ density`.

**Arguments**

  * `current::Real`: stimulation current [mA]
  * `pad_area::Real`: electrode pad area [cm^2]
  * `duration::Int64`: stimulation duration [s]

**Returns**

  * `charge::Float64`: charge [C]
  * `current_density::Float64`: current density [A/m^2]
  * `charge_density::Float64`: delibvered charge density [kC/m^2]

**Source**

Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.

