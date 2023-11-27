export import_npy

"""
    import_npy(file_name; sampling_rate)

Load NPY file (exported from MNE) and return `NeuroAnalyzer.NEURO` object. Data type and channel types are set as  is EEG.

# Arguments

- `file_name::String`: name of the file to load
- `sampling_rate::Int64`: NPY file contains only signal data, therefore its sampling rate must be provided upon importing

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_npy(file_name::String; sampling_rate::Int64)

    @assert sampling_rate > 1 "Sampling rate must be ≥ 1."
    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert splitext(file_name)[2] == ".npy" "This is not NPY file."

    file_type = "NPY"

    data = npzread(file_name)
    ch_n = size(data, 1)
    clabels = String[]
    for idx in 1:ch_n
        push!(clabels, "ch_$idx")
    end
    ch_type = repeat(["eeg"], ch_n)
    units = repeat(["μV"], ch_n)
    channel_order = collect(1:ch_n)

    x = zeros(ch_n)
    y = zeros(ch_n)
    z = zeros(ch_n)
    theta = zeros(ch_n)
    radius = zeros(ch_n)
    phi_sph = zeros(ch_n)
    radius_sph = zeros(ch_n)
    theta_sph = zeros(ch_n)
    radius_sph == zeros(ch_n) && (radius_sph = radius)
    locs = DataFrame(:labels=>clabels,
                     :loc_theta=>theta,
                     :loc_radius=>radius,
                     :loc_x=>x,
                     :loc_y=>y,
                     :loc_z=>z,
                     :loc_radius_sph=>radius_sph,
                     :loc_theta_sph=>theta_sph,
                     :loc_phi_sph=>phi_sph)

    markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])

    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    epoch_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    data_type = "eeg"

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name="",
                        head_circumference=-1,
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_eeg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type=file_type,
                              recording="",
                              recording_date="",
                              recording_time="",
                              recording_notes="",
                              channel_type=ch_type,
                              reference="",
                              clabels=clabels,
                              transducers=repeat([""], ch_n),
                              units=repeat(["μV"], ch_n),
                              prefiltering=repeat([""], ch_n),
                              sampling_rate=sampling_rate,
                              gain=repeat([1.0], ch_n))
    e = _create_experiment(name="",
                           notes="",
                           design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()
    history = String[]

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data[channel_order, :, :], components, markers, locs, history)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(obj.time_pts[end]) s)")

    return obj

end
