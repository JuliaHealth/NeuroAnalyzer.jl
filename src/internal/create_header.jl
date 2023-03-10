function _create_subject(;id::String, first_name::String, middle_name::String, last_name::String, handedness::String, weight::Real, height::Real)

    return Dict(:subject_id=>id,
                :subject_first_name=>first_name,
                :subject_middle_name=>middle_name,
                :subject_last_name=>last_name,
                :subject_handedness=>handedness,
                :subject_weight=>weight,
                :subject_height=>height)
end

function _create_recording_eeg(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_n::Int64, channel_type::Vector{String}, reference::String, duration_samples::Int64, duration_seconds::Float64, epoch_n::Int64, epoch_duration_samples::Int64, epoch_duration_seconds::Float64, clabels::Vector{String}, transducers::Vector{String},units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64, gain::Vector{Float64})

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_n=>channel_n,
                :channel_type=>channel_type,
                :reference=>reference,
                :duration_samples=>duration_samples,
                :duration_seconds=>duration_seconds,
                :epoch_n=>epoch_n,
                :epoch_duration_samples=>epoch_duration_samples,
                :epoch_duration_seconds=>epoch_duration_seconds,
                :labels=>clabels,
                :transducers=>transducers,
                :units=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate,
                :gain=>gain)
end

function _create_recording_meg(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_n::Int64, channel_type::Vector{String}, reference::String, duration_samples::Int64, duration_seconds::Float64, epoch_n::Int64, epoch_duration_samples::Int64, epoch_duration_seconds::Float64, clabels::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64, magnetometers::Vector{Int64}, gradiometers::Vector{Int64}, gradiometers_planar::Vector{Int64}, gradiometers_axial::Vector{Int64}, coils::Vector{Int64})

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>replace(recording_time, '.'=>':'),
                :recording_notes=>recording_notes,
                :channel_n=>channel_n,
                :channel_type=>channel_type,
                :reference=>reference,
                :duration_samples=>duration_samples,
                :duration_seconds=>duration_seconds,
                :epoch_n=>epoch_n,
                :epoch_duration_samples=>epoch_duration_samples,
                :epoch_duration_seconds=>epoch_duration_seconds,
                :labels=>clabels,
                :units=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate,
                :magnetometers=>magnetometers,
                :gradiometers=>gradiometers,
                :gradiometers_planar=>gradiometers_planar,
                :gradiometers_axial=>gradiometers_axial,
                :coils=>coils)
end

function _create_experiment(;experiment_name::String, experiment_notes::String, experiment_design::String)

    return Dict(:experiment_name=>experiment_name,
                :experiment_notes=>experiment_notes,
                :experiment_design=>experiment_design)
end

function _create_header(subject::Dict, recording::Dict, experiment::Dict; component_names::Vector{Symbol}, has_markers::Bool, has_locs::Bool, history::Vector{String})

    return NeuroAnalyzer.HEADER(subject,
                                recording,
                                experiment,
                                component_names,
                                has_markers,
                                has_locs,
                                history)
end
