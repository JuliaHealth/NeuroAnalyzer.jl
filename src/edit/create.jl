export create

"""
    create(; data_type)

Create an empty `NeuroAnalyzer.NEURO` object.

# Arguments

- `data_type::String`: data type of the new object ("eeg", "meg", "nirs")

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function create(; data_type::String)

    NeuroAnalyzer._check_var(data_type, data_types, "data_type")

    markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])

    time_pts = Float64[]
    ep_time = Float64[]

    s = NeuroAnalyzer._create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name="",
                        handedness="",
                        weight=-1,
                        height=-1)
    r = NeuroAnalyzer._create_recording_eeg(data_type=data_type,
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
    e = NeuroAnalyzer._create_experiment(name="", notes="", design="")

    hdr = NeuroAnalyzer._create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = DataFrame(:channel=>Int64,
                     :labels=>String[],
                     :loc_theta=>Float64[],
                     :loc_radius=>Float64[],
                     :loc_x=>Float64[],
                     :loc_y=>Float64[],
                     :loc_z=>Float64[],
                     :loc_radius_sph=>Float64[],
                     :loc_theta_sph=>Float64[],
                     :loc_phi_sph=>Float64[])

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, Float64[], components, markers, locs, history)

    return obj
    
end

