# ---------------------------------------------------------------------------
# internal helper: build the fields shared by every recording type
# ---------------------------------------------------------------------------
function _common_recording_fields(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    clabels::Vector{String},
    units::Vector{String},
    sampling_rate::Int64,
    bad_channels::Vector{Bool},
)::Dict
    return Dict(
        :data_type => data_type,
        :file_name => file_name,
        :file_size_mb => file_size_mb,
        :file_type => file_type,
        :recording => recording,
        :recording_date => recording_date,
        :recording_time => recording_time,
        :recording_notes => recording_notes,
        :channel_type => channel_type,
        :channel_order => channel_order,
        :label => clabels,
        :unit => units,
        :sampling_rate => sampling_rate,
        :bad_channel => bad_channels,
    )
end

# ---------------------------------------------------------------------------
# subject
# ---------------------------------------------------------------------------

function _create_subject(;
    id::String,
    first_name::String,
    middle_name::String,
    last_name::String,
    handedness::String,
    head_circumference::Real,
    weight::Real,
    height::Real,
)::Dict
    return Dict(
        :id => id,
        :first_name => first_name,
        :middle_name => middle_name,
        :last_name => last_name,
        :handedness => handedness,
        :head_circumference => head_circumference,
        :weight => weight,
        :height => height,
    )
end

# ---------------------------------------------------------------------------
# EEG / sEEG / ECoG  (identical schema; data_type distinguishes them)
# ---------------------------------------------------------------------------

function _create_recording_eeg_like(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    reference::String,
    clabels::Vector{String},
    transducers::Vector{String},
    units::Vector{String},
    prefiltering::Vector{String},
    line_frequency::Real,
    sampling_rate::Int64,
    gain::Vector{Float64},
    bad_channels::Vector{Bool},
)::Dict
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
    )
    merge!(
        d,
        Dict(
            :reference => reference,
            :transducers => transducers,
            :prefiltering => prefiltering,
            :line_frequency => line_frequency,
            :gain => gain,
            :epoch_id => "",
        ),
    )
    return d
end

# public thin wrappers - callers use the modality-specific name; all three delegate to the shared implementation above.
_create_recording_eeg(; kwargs...) = _create_recording_eeg_like(; kwargs...)
_create_recording_seeg(; kwargs...) = _create_recording_eeg_like(; kwargs...)
_create_recording_ecog(; kwargs...) = _create_recording_eeg_like(; kwargs...)

# ---------------------------------------------------------------------------
# MEG
# ---------------------------------------------------------------------------

function _create_recording_meg(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    reference::String,
    clabels::Vector{String},
    units::Vector{String},
    prefiltering::Vector{String},
    line_frequency::Real,
    sampling_rate::Int64,
    magnetometers::Vector{Int64},
    gradiometers::Vector{Int64},
    coil_type::Vector{String},
    bad_channels::Vector{Bool},
    ssp_labels::Vector{String},
    ssp_channels::Vector{Bool},
    ssp_data::Matrix{Float64},
    maxshield::Bool = false,
    coil_loc::Matrix{Float64} = Matrix{Float64}(undef, 0, 12),
    dev_head_t::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
    coil_id::Vector{Int64} = Int64[],
    head_origin::Vector{Float64} = Float64[],
    hpi::Dict{Symbol, Any} = Dict{Symbol, Any}(),
    first_time::Real = 0.0,
)::Dict
    # sensor geometry (FIFF/MNE `loc` convention), one row per channel:
    # r0 (1:3, m), ex (4:6), ey (7:9), ez (10:12; coil normal)
    # MEG channels: device frame; EEG channels: head frame (4:6 = reference electrode)
    (isempty(coil_loc) || size(coil_loc) == (length(clabels), 12)) || throw(
        ArgumentError("coil_loc must be an n_channels × 12 matrix."),
    )
    (isempty(dev_head_t) || size(dev_head_t) == (4, 4)) ||
        throw(ArgumentError("dev_head_t must be a 4 × 4 matrix."))
    # FIFF coil type code per channel (selects coil integration points)
    (isempty(coil_id) || length(coil_id) == length(clabels)) ||
        throw(ArgumentError("coil_id must have one entry per channel."))
    # head origin (sphere fitted to digitization, head frame, m)
    (isempty(head_origin) || length(head_origin) == 3) ||
        throw(ArgumentError("head_origin must have 3 elements."))
    # hpi: :freqs (Hz), :dig (n_coils × 3, head frame, m), :coil_no
    # first_time: time (s) of the first sample from the start of the acquisition (FIFF
    # first_samp / sampling rate); MaxFilter / MNE-Python head positions are relative to it
    # normalize time separator (MEG files sometimes use '.' instead of ':')
    recording_time = replace(recording_time, '.' => ':')
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
    )
    merge!(
        d,
        Dict(
            :reference => reference,
            :prefiltering => prefiltering,
            :line_frequency => line_frequency,
            :magnetometers => magnetometers,
            :gradiometers => gradiometers,
            :coil_type => coil_type,
            :ssp_labels => ssp_labels,
            :ssp_channels => ssp_channels,
            :ssp_data => ssp_data,
            :maxshield => maxshield,
            :coil_loc => coil_loc,
            :dev_head_t => dev_head_t,
            :coil_id => coil_id,
            :head_origin => head_origin,
            :hpi => hpi,
            :first_time => Float64(first_time),
            :sss => false,
            :tsss => false,
            :epoch_id => "",
        ),
    )
    return d
end

# ---------------------------------------------------------------------------
# NIRS
# ---------------------------------------------------------------------------

function _create_recording_nirs(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    wavelengths::Vector{Float64},
    wavelength_index::Vector{Int64},
    optode_pairs::Matrix{Int64},
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    clabels::Vector{String},
    units::Vector{String},
    src_labels::Vector{String},
    det_labels::Vector{String},
    opt_labels::Vector{String},
    sampling_rate::Int64,
    bad_channels::Vector{Bool},
)::Dict
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
    )
    merge!(
        d,
        Dict(
            :wavelengths => wavelengths,
            :wavelength_index => wavelength_index,
            :optode_pairs => optode_pairs,
            :src_labels => src_labels,
            :det_labels => det_labels,
            :optode_labels => opt_labels,
            :epoch_id => "",
        ),
    )
    return d
end

# ---------------------------------------------------------------------------
# sensors / EDA / TPT  (identical schema; data_type distinguishes them)
# ---------------------------------------------------------------------------

function _create_recording_prefiltered(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    clabels::Vector{String},
    units::Vector{String},
    prefiltering::Vector{String},
    sampling_rate::Int64,
    bad_channels::Vector{Bool},
)::Dict
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
    )
    merge!(
        d, Dict(
        :prefiltering => prefiltering,
        :epoch_id => "",
    ),
    )
    return d
end

_create_recording_sensors(; kwargs...) = _create_recording_prefiltered(; kwargs...)
_create_recording_eda(; kwargs...) = _create_recording_prefiltered(; kwargs...)
_create_recording_tpt(; kwargs...) = _create_recording_prefiltered(; kwargs...)

# ---------------------------------------------------------------------------
# MEP
# ---------------------------------------------------------------------------

function _create_recording_mep(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    clabels::Vector{String},
    units::Vector{String},
    sampling_rate::Int64,
    stimulation_intensity::Vector{Int64},
    coil_type::Vector{String},
    stimulation_sample::Vector{Int64},
    markers_pos::Vector{Int64},
    markers_neg::Vector{Int64},
    bad_channels::Vector{Bool},
)::Dict
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
    )
    merge!(
        d,
        Dict(
            :stimulation_intensity => stimulation_intensity,
            :coil_type => coil_type,
            :stimulation_sample => stimulation_sample,
            :markers_pos => markers_pos,
            :markers_neg => markers_neg,
            # note: MEP does not use :epoch_id (intentional - MEP data is not epoched)
        ),
    )
    return d
end

# ---------------------------------------------------------------------------
# header
# ---------------------------------------------------------------------------

function _create_experiment(; name::String, notes::String, design::String)::Dict
    return Dict(:name => name, :notes => notes, :design => design)
end

function _create_header(;
    subject::Dict,
    recording::Dict,
    experiment::Dict,
)::NeuroAnalyzer.HEADER
    return NeuroAnalyzer.HEADER(subject, recording, experiment)
end
