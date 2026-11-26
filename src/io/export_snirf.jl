export export_snirf

# NeuroAnalyzer channel types → SNIRF data type labels (processed data, dataType 99999)
const _SNIRF_LABELS = Dict(
    "nirs_od" => "dOD", "nirs_dmean" => "dMean", "nirs_dvar" => "dVar",
    "nirs_dskew" => "dSkew", "nirs_mua" => "mua", "nirs_musp" => "musp",
    "nirs_hbo" => "HbO", "nirs_hbr" => "HbR", "nirs_hbt" => "HbT",
    "nirs_h2o" => "H2O", "nirs_lipid" => "Lipid", "nirs_bfi" => "BFi",
    "nirs_hrf_dod" => "HRF dOD", "nirs_hrf_dmean" => "HRF dMean",
    "nirs_hrf_dvar" => "HRF dVar", "nirs_hrf_dskew" => "HRF dSkew",
    "nirs_hrf_hbo" => "HRF HbO", "nirs_hrf_hbr" => "HRF HbR",
    "nirs_hrf_hbt" => "HRF HbT", "nirs_hrf_bfi" => "HRF BFi",
)
# data types that do not depend on the wavelength (chromophore concentrations, flow)
const _SNIRF_NO_WL = [
    "nirs_hbo", "nirs_hbr", "nirs_hbt", "nirs_h2o", "nirs_lipid", "nirs_bfi",
    "nirs_hrf_hbo", "nirs_hrf_hbr", "nirs_hrf_hbt", "nirs_hrf_bfi",
]

"""
    export_snirf(obj; <keyword arguments>)

Export a NIRS object to a SNIRF file (Shared Near Infrared Spectroscopy Format, HDF5).

The optical channels (intensity: data type 1, continuous wave amplitude; OD and concentration channels: data type 99999, processed, with the SNIRF data type labels `dOD`, `HbO`, `HbR`, `HbT`, etc.) are written to `/nirs/data1`, the optode labels and 3-D positions (in the length unit of the recording, mm if unknown) to `/nirs/probe`, the markers to one stimulus group per marker value (`/nirs/stim1`, …: onset, duration, value 1) and the auxiliary channels to `/nirs/aux1`, …

The object must be continuous (one epoch).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: file name (`.snirf`)
- `overwrite::Bool=false`: if `true`, overwrite an existing file

# Returns

- `Nothing`

# References

1. Tucker S, Dubb J, Kura S, et al. Introduction to the shared near infrared spectroscopy format. Neurophotonics. 2023;10(1):013507.
"""
function export_snirf(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    overwrite::Bool = false,
)::Nothing
    # validate
    _check_export_datatype(obj, :snirf)
    nepochs(obj) == 1 || throw(ArgumentError("OBJ must be continuous (one epoch)."))
    lowercase(splitext(file_name)[2]) == ".snirf" ||
        throw(ArgumentError("File name extension must be .snirf"))
    (isfile(file_name) && !overwrite) && throw(
        ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."),
    )

    rec = obj.header.recording
    chp = rec[:optode_pairs]
    n_opt = size(chp, 1)
    n_opt > 0 || throw(ArgumentError("OBJ does not contain NIRS optical channels."))
    all(>(0), chp) || throw(ArgumentError("All optical channels must have optode pairs."))
    ct = rec[:channel_type]
    aux = findall(==("nirs_aux"), ct)
    other = setdiff((n_opt + 1):nchannels(obj), aux)
    isempty(other) ||
        _warn("Channels $(labels(obj)[other]) are not NIRS channels and are not exported.")

    # optode positions
    function _pos(l)
        i = findfirst(==(l), obj.locs[!, :label])
        isnothing(i) && return [0.0, 0.0, 0.0]
        return [obj.locs[i, :loc_x], obj.locs[i, :loc_y], obj.locs[i, :loc_z]]
    end
    src_pos = reduce(hcat, [_pos(l) for l in rec[:src_labels]])
    det_pos = reduce(hcat, [_pos(l) for l in rec[:det_labels]])
    length_unit = string(get(rec, :length_unit, ""))
    length_unit == "" && (length_unit = "mm")

    rm(file_name; force = true)
    JLD2.jldopen(file_name, "w") do f
        f["formatVersion"] = "1.1"
        g = "nirs"
        # metadata
        f["$g/metaDataTags/SubjectID"] = string(obj.header.subject[:id])
        f["$g/metaDataTags/MeasurementDate"] =
            rec[:recording_date] == "" ? "unknown" : string(rec[:recording_date])
        f["$g/metaDataTags/MeasurementTime"] =
            rec[:recording_time] == "" ? "unknown" : string(rec[:recording_time])
        f["$g/metaDataTags/LengthUnit"] = length_unit
        f["$g/metaDataTags/TimeUnit"] = "s"
        f["$g/metaDataTags/FrequencyUnit"] = "Hz"

        # data: (channels, samples) in Julia = (time, channels) in HDF5 (row-major)
        f["$g/data1/dataTimeSeries"] = obj.data[1:n_opt, :, 1]
        f["$g/data1/time"] = collect(obj.time_pts)
        for j = 1:n_opt
            m = "$g/data1/measurementList$j"
            f["$m/sourceIndex"] = chp[j, 1]
            f["$m/detectorIndex"] = chp[j, 2]
            f["$m/wavelengthIndex"] = rec[:wavelength_index][j]
            if ct[j] == "nirs_int"
                f["$m/dataType"] = 1
            else
                f["$m/dataType"] = 99999
                f["$m/dataTypeLabel"] = get(_SNIRF_LABELS, ct[j], replace(ct[j], "nirs_" => ""))
            end
            f["$m/dataTypeIndex"] = 1
            rec[:unit][j] != "" && (f["$m/dataUnit"] = rec[:unit][j])
        end

        # probe
        f["$g/probe/wavelengths"] = Float64.(rec[:wavelengths])
        f["$g/probe/sourcePos3D"] = src_pos
        f["$g/probe/detectorPos3D"] = det_pos
        f["$g/probe/sourceLabels"] = rec[:src_labels]
        f["$g/probe/detectorLabels"] = rec[:det_labels]

        # stimuli: one group per marker value (onset, duration, value)
        if _has_markers(obj)
            for (s_idx, v) in enumerate(unique(string.(obj.markers[!, :value])))
                m = obj.markers[string.(obj.markers[!, :value]) .== v, :]
                f["$g/stim$s_idx/name"] = v
                f["$g/stim$s_idx/data"] =
                    permutedims(hcat(m[!, :start], m[!, :length], ones(DataFrames.nrow(m))))
            end
        end

        # auxiliary channels
        for (a_idx, a) in enumerate(aux)
            f["$g/aux$a_idx/name"] = labels(obj)[a]
            f["$g/aux$a_idx/dataTimeSeries"] = obj.data[a, :, 1]
            f["$g/aux$a_idx/time"] = collect(obj.time_pts)
        end
    end

    return nothing
end
