export import_nirs

"""
    import_nirs(file_name)

Load a Homer3-compatible NIRS file and return a `NeuroAnalyzer.NEURO` object.

The NIRS format is a MATLAB `.mat` file containing raw intensity data (`d`), a time vector (`t`), a probe structure (`SD`), stimulus onset matrix (`s`), and optional auxiliary channels (`aux`).

# Arguments

- `file_name::String`: path to the `.nirs` file

# Returns

- `NeuroAnalyzer.NEURO`

# Throws
- `ArgumentError` if the file does not exist, is not a NIRS file, or is missing required fields

# References

1. https://github.com/BUNPC/Homer3/wiki/HOMER3-file-formats
"""
function import_nirs(file_name::String)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) == ".nirs" ||
        throw(ArgumentError("$file_name is not a NIRS file."))

    nirs = try
        matread(file_name)
    catch e
        throw(ArgumentError("File $file_name cannot be loaded: $e"))
    end

    "d" in keys(nirs) || throw(ArgumentError("$file_name is missing required 'd' field."))
    "t" in keys(nirs) || throw(ArgumentError("$file_name is missing required 't' field."))

    file_type = "NIRS"

    # ------------------------------------------------------------------ #
    # time axis                                                          #
    # ------------------------------------------------------------------ #
    time_pts = Float64.(nirs["t"][:])
    time_pts .-= time_pts[1]
    epoch_time = time_pts
    sampling_rate = round(Int64, 1 / (time_pts[2] - time_pts[1]))

    # ------------------------------------------------------------------ #
    # probe geometry                                                     #
    # ------------------------------------------------------------------ #
    probes = nirs["SD"]
    meas = Int.(probes["MeasList"]) # ch_n × 4
    ch_n = size(meas, 1)

    # columns: 1=source index, 2=detector index, 3=unused, 4=wavelength index
    opt_pairs = meas[:, 1:2]
    wavelengths = Float64.(probes["Lambda"][:])
    wavelength_index = meas[:, 4]

    # optode labels: sources first, then detectors
    n_src = Int(probes["nSrcs"])
    n_det = Int(probes["nDets"])
    src_labels = ["S$i" for i in 1:n_src]
    det_labels = ["D$i" for i in 1:n_det]
    opt_labels = vcat(src_labels, det_labels)

    # channel labels: "S<src>_D<det> <wavelength_nm>"
    clabels = [
        replace(
            "S$(meas[i,1])_D$(meas[i,2]) $(wavelengths[wavelength_index[i]])",
            ".0" => "")
        for i in 1:ch_n]

    # ------------------------------------------------------------------ #
    # signal data (intensity, raw)                                       #
    # ------------------------------------------------------------------ #
    data     = Matrix(nirs["d"]')          # (ch_n × n_samples)
    ch_type  = repeat(["nirs_int"], ch_n)
    data_unit = repeat(["V"], ch_n)

    # ------------------------------------------------------------------ #
    # stimuli → markers                                                   #
    # ------------------------------------------------------------------ #
    stim_raw = nirs["s"]
    markers = if all(iszero, stim_raw)
        DataFrame(
            :id => String[],
            :start => Float64[],
            :length => Float64[],
            :value => String[],
            :channel => Int64[]
        )
    else
        stim = ndims(stim_raw) == 1 ? reshape(stim_raw, :, 1) : Matrix(stim_raw)
        s_n = size(stim, 2)
        df = DataFrame(
            :id => String[],
            :start => Float64[],
            :length => Float64[],
            :value => String[],
            :channel => Int64[]
        )
        for col in 1:s_n
            for samp in findall(!iszero, stim[:, col])
                push!(df, ("", time_pts[samp], 0.0, "stim", 0))
            end
        end
        # assign sequential IDs per unique value
        for (i, v) in enumerate(unique(df[!, :value]))
            df[df[!, :value] .== v, :id] .= string(i)
        end
        df
    end

    # ------------------------------------------------------------------ #
    # auxiliary channels                                                  #
    # ------------------------------------------------------------------ #
    aux = Matrix(nirs["aux"]')
    if length(aux) > 0
        data = vcat(data, aux)
        for i in axes(aux, 1)
            push!(clabels, "AUX$i")
            push!(ch_type, "nirs_aux")
            push!(data_unit, "")
        end
    end

    data = reshape(data, size(data, 1), size(data, 2), 1)

    # ------------------------------------------------------------------ #
    # pptode locations                                                   #
    # SrcPos/DetPos are (n × 3) in the file; transpose to (3 × n)        #
    # ------------------------------------------------------------------ #
    src_pos = Matrix(probes["SrcPos"]') # 3 × n_src
    det_pos = Matrix(probes["DetPos"]') # 3 × n_det
    pos3d = hcat(src_pos, det_pos) # 3 × (n_src + n_det)
    x = Float64.(pos3d[1, :])
    y = Float64.(pos3d[2, :])
    z = Float64.(pos3d[3, :])

    if !all(iszero, z)
        x, y, z = _locs_norm(x, y, z)
    else
        x_n, y_n = _locs_norm(x, y)
        x = x_n; y = y_n
    end

    n_opt = length(opt_labels)
    locs = DataFrame(
        :label => opt_labels,
        :loc_radius => zeros(n_opt),
        :loc_theta => zeros(n_opt),
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => zeros(n_opt),
        :loc_theta_sph => zeros(n_opt),
        :loc_phi_sph => zeros(n_opt)
    )
    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                               #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(
        id = "", first_name = "",
        middle_name = "",
        last_name = "",
        head_circumference = -1,
        handedness = "",
        weight = -1,
        height = -1)
    r = _create_recording_nirs(
        data_type = "nirs",
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = file_type,
        recording = "", recording_date = "", recording_time = "",
        recording_notes = "",
        wavelengths = wavelengths,
        wavelength_index = wavelength_index,
        optode_pairs = opt_pairs,
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        clabels = clabels,
        units = data_unit,
        src_labels = src_labels,
        det_labels = det_labels,
        opt_labels = opt_labels,
        sampling_rate = sampling_rate,
        bad_channels = zeros(Bool, size(data, 1)))
    e   = _create_experiment(name = "", notes = "", design = "")
    hdr = _create_header(subject = s, recording = r, experiment = e)

    obj = NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, epoch_time, data)

    _info("Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end]; digits=2)) s)")

    return obj

end
