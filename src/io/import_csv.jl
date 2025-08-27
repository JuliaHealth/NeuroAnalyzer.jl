export import_csv

"""
    import_csv(file_name; <keyword arguments>)

Load CSV file (e.g. exported from EEGLAB) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

CSV first row or first column must contain channel names.
Shape of data array will be detected automatically. Sampling rate will be detected.
If file is gzip-ed, it will be uncompressed automatically while reading.
"""
function import_csv(file_name::String; detect_type::Bool=true)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."

    file_type = "CSV"

    df = CSV.read(file_name, stringtype=String, DataFrame)

    if df[:, 1] isa Vector{Float64}
        # time by channels
        time_pts = df[!, 1]
        data = Array(df[:, 2:end])'
        ch_n = ncol(df) - 1
        clabels_tmp = String.(names(df)[2:end])
    else
        # channels by time
        time_pts = parse.(Float64, names(df)[2:end])
        data = Array(df[!, 2:end])
        ch_n = nrow(df)
        clabels_tmp = String.(df[:, 1])
    end
    data = Array(reshape(data, size(data, 1), size(data, 2), 1))

    clabels = _clean_labels(clabels_tmp)
    if detect_type
        ch_type = _set_channel_types(clabels, "eeg")
    else
        ch_type = repeat(["eeg"], ch_n)
    end
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :value=>String[],
                        :channel=>Int64[])

    sampling_rate = round(Int64, 1 / time_pts[2] * 1000)
    gain = ones(ch_n)
    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :value=>String[],
                        :channel=>Int64[])

    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    ep_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)

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
                              channel_order=_sort_channels(ch_type),
                              reference=_detect_montage(clabels, ch_type, data_type),
                              clabels=clabels,
                              transducers=repeat([""], ch_n),
                              units=units,
                              prefiltering=repeat([""], ch_n),
                              line_frequency=50,
                              sampling_rate=sampling_rate,
                              gain=gain,
                              bad_channels=zeros(Bool, size(data, 1), 1))
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data, components, markers, locs, history)
    _initialize_locs!(obj)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=3)) s)")

    return obj

end
