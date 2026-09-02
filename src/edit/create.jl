export create_object
export create_time
export create_time!
export create_data
export create_data!

"""
    create_object(; <keyword arguments>)

Create an empty `NeuroAnalyzer.NEURO` object of the specified data type.

All data, time, and header fields are initialized to empty/zero values. Use `create_data!` to populate with signal data and `create_time!` to set the sampling rate and build time vectors.

# Arguments

- `data_type::String`: data type of the new object (must be a recognized type)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function create_object(; data_type::String)::NeuroAnalyzer.NEURO
    # validate
    _check_var(data_type, data_types, "data_type")

    markers = DataFrame(
        :id => String[],
        :start => Float64[],
        :length => Float64[],
        :value => String[],
        :channel => Int64[],
    )

    time_pts = Float64[]
    ep_time  = Float64[]
    data     = Array{Float64, 3}(undef, 0, 0, 0)

    s = _create_subject(;
        id                 = "",
        first_name         = "",
        middle_name        = "",
        last_name          = "",
        handedness         = "",
        head_circumference = -1,
        weight             = -1,
        height             = -1,
    )

    # common fields shared by most recording types
    _common = (
        file_name       = "",
        file_size_mb    = 0,
        file_type       = "",
        recording       = "",
        recording_date  = "",
        recording_time  = "",
        recording_notes = "",
        channel_type    = String[],
        channel_order   = Int64[],
        clabels         = String[],
        units           = String[],
        sampling_rate   = 0,
        bad_channels    = Bool[],
    )

    if data_type == "eeg"
        r = _create_recording_eeg(;
            data_type = "eeg",
            _common...,
            reference = "",
            transducers = String[],
            prefiltering = String[],
            line_frequency = 50,
            gain = Float64[],
        )
    elseif data_type == "seeg"
        r = _create_recording_seeg(;
            data_type = "seeg",
            _common...,
            reference = "",
            transducers = String[],
            prefiltering = String[],
            line_frequency = 50,
            gain = Float64[],
        )
    elseif data_type == "ecog"
        r = _create_recording_ecog(;
            data_type = "ecog",
            _common...,
            reference = "",
            transducers = String[],
            prefiltering = String[],
            line_frequency = 50,
            gain = Float64[],
        )
    elseif data_type == "meg"
        r = _create_recording_meg(;
            data_type = "meg",
            _common...,
            reference      = "",
            prefiltering   = String[],
            line_frequency = 50,
            magnetometers  = Int[],
            gradiometers   = Int[],
            coil_type      = String[],
            ssp_labels     = String[],
            ssp_channels   = Bool[],
            ssp_data       = Matrix{Float64}(undef, 0, 0),
        )
    elseif data_type == "nirs"
        r = _create_recording_nirs(;
            data_type = "nirs",
            _common...,
            wavelengths      = Float64[],
            wavelength_index = Int64[],
            optode_pairs     = Matrix{Int64}(undef, 0, 0),
            src_labels       = String[],
            det_labels       = String[],
            opt_labels       = String[],
        )
    elseif data_type == "sensors"
        r = _create_recording_sensors(;
            data_type = "sensors",
            _common...,
            prefiltering = String[],
        )
    elseif data_type == "mep"
        r = _create_recording_mep(;
            data_type = "mep",
            _common...,
            stimulation_intensity = Int64[],
            coil_type             = String[],
            stimulation_sample    = Int64[],
            markers_pos           = Int64[],
            markers_neg           = Int64[],
        )
    elseif data_type == "eda"
        r = _create_recording_eda(;
            data_type = "eda",
            _common...,
            prefiltering = String[],
        )
    elseif data_type == "tpt"
        r = _create_recording_tpt(;
            data_type = "tpt",
            _common...,
            prefiltering = String[],
        )
    end

    e = _create_experiment(; name = "", notes = "", design = "")
    hdr = _create_header(; subject = s, recording = r, experiment = e)
    locs = _initialize_locs()

    return NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, ep_time, data)
end

"""
    create_time(obj; <keyword arguments>)

Return a copy of `obj` with time vectors built from a specified sampling rate.

`obj` must already contain data and must not already have time points assigned.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `fs::Int64`: sampling rate in Hz; must be > 0

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function create_time(obj::NeuroAnalyzer.NEURO; fs::Int64)::NeuroAnalyzer.NEURO
    # validate
    length(obj.data) > 0 || throw(ArgumentError("OBJ does not contain data."))
    length(obj.time_pts) == 0 || throw(ArgumentError("OBJ already has time points."))
    fs > 0 || throw(ArgumentError("fs must be > 0."))

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.header.recording[:sampling_rate] = fs
    obj_tmp.time_pts, obj_tmp.epoch_time = _get_t(obj_tmp)
    push!(obj_tmp.history, "create_time(obj, fs=$fs)")

    return obj_tmp
end

"""
    create_time!(obj; <keyword arguments>)

Build time vectors from a specified sampling rate in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `fs::Int64`: sampling rate in Hz; must be > 0

# Returns

- `Nothing`
"""
function create_time!(obj::NeuroAnalyzer.NEURO; fs::Int64)::Nothing
    obj_tmp = create_time(obj; fs = fs)
    obj.header = obj_tmp.header
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end

"""
    create_data(obj; <keyword arguments>)

Return a copy of `obj` with data, channel metadata and time vectors populated.

Auto-generates channel labels of the form `"ch-1"`, `"ch-2"`, … and sets channel type and units uniformly for all new channels.

`obj` must be empty (no data, no time points).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Array{Float64, 3}`: signal data, shape (ch_n, epoch_len, n_epochs)
- `fs::Int64`: sampling rate in Hz; must be > 0
- `type::String`: channel type applied to all channels (must be a recognized type)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function create_data(
    obj::NeuroAnalyzer.NEURO;
    data::Array{Float64, 3},
    fs::Int64,
    type::String,
)::NeuroAnalyzer.NEURO
    # validate
    length(obj.data) == 0 || throw(ArgumentError("OBJ already contains data."))
    length(obj.time_pts) == 0 || throw(ArgumentError("OBJ already has time points."))
    _check_var(type, channel_types, "type")
    fs > 0 || throw(ArgumentError("fs must be > 0."))

    # number of channels
    ch_n = size(data, 1)

    # channel labels
    clabels = ["ch-$i" for i in 1:ch_n]

    # create new dataset
    obj_tmp                                  = deepcopy(obj)
    obj_tmp.data                             = data
    obj_tmp.header.recording[:label]         = clabels
    obj_tmp.header.recording[:channel_type]  = fill(type, ch_n)
    obj_tmp.header.recording[:unit]          = fill(_ch_units(type), ch_n)
    obj_tmp.header.recording[:channel_order] = collect(1:ch_n)
    obj_tmp.header.recording[:bad_channel]   = zeros(Bool, ch_n)
    obj_tmp.header.recording[:sampling_rate] = fs
    obj_tmp.time_pts, obj_tmp.epoch_time     = _get_t(obj_tmp)
    push!(obj_tmp.history, "create_data(obj; data, fs=$fs, type=$type)")

    return obj_tmp
end

"""
    create_data!(obj; <keyword arguments>)

Populate `obj` with data, channel metadata and time vectors in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Array{Float64, 3}`: signal data, shape (ch_n, epoch_len, n_epochs)
- `fs::Int64`: sampling rate in Hz; must be > 0
- `type::String`: channel type applied to all channels (must be a recognized type)

# Returns

- `Nothing`
"""
function create_data!(
    obj::NeuroAnalyzer.NEURO;
    data::Array{Float64, 3},
    fs::Int64,
    type::String,
)::Nothing
    obj_tmp = create_data(obj; data = data, fs = fs, type = type)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
