export create
export create_time
export create_time!

"""
    create(; data_type)

Create an empty `NeuroAnalyzer.NEURO` object.

# Arguments

- `data_type::String`: data type of the new object ("eeg", "meg", "nirs", "ecog")

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function create(; data_type::String)

    _check_var(data_type, data_types, "data_type")

    markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])

    time_pts = Float64[]
    ep_time = Float64[]

    data = Array{Float64, 3}(undef, 0, 0, 0)

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name="",
                        handedness="",
                        head_circumference=-1,
                        weight=-1,
                        height=-1)
    r = _create_recording_eeg(data_type=data_type,
                              file_name="",
                              file_size_mb=0,
                              file_type="",
                              recording="",
                              recording_date="",
                              recording_time="",
                              recording_notes="",
                              channel_type=String[],
                              reference="",
                              clabels=String[],
                              transducers=String[],
                              units=String[],
                              prefiltering=String[],
                              sampling_rate=0,
                              gain=Float64[])
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = DataFrame(:channel=>Int64,
                     :labels=>String[],
                     :loc_radius=>Float64[],
                     :loc_theta=>Float64[],
                     :loc_x=>Float64[],
                     :loc_y=>Float64[],
                     :loc_z=>Float64[],
                     :loc_radius_sph=>Float64[],
                     :loc_theta_sph=>Float64[],
                     :loc_phi_sph=>Float64[])

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data, components, markers, locs, history)

    return obj
    
end

"""
    create_time(obj)

Create time points vector for `NeuroAnalyzer.NEURO` object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `fs::Int64

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function create_time(obj::NeuroAnalyzer.NEURO; fs::Int64)

    @assert length(obj.data) > 0 "OBJ does not contain data."
    @assert length(obj.time_pts) == 0 "OBJ already has time points."

    obj_new = deepcopy(obj)
    obj_new.header.recording[:sampling_rate] = fs
    obj_new.time_pts, obj_new.epoch_time = NeuroAnalyzer._get_t(obj_new)
    reset_components!(obj_new)
    push!(obj_new.history, "create_time(OBJ, fs=$fs)")

    return obj_new

end

"""
    create_time!(obj)

Create time points vector for `NeuroAnalyzer.NEURO` object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `fs::Int64
"""
function create_time!(obj::NeuroAnalyzer.NEURO; fs::Int64)

    obj_new = create_time(obj, fs=fs)
    obj.header = obj_new.header
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing

end
