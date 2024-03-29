function _create_subject(;id::String, first_name::String, middle_name::String, last_name::String, handedness::String, head_circumference::Real, weight::Real, height::Real)

    return Dict(:id=>id,
                :first_name=>first_name,
                :middle_name=>middle_name,
                :last_name=>last_name,
                :handedness=>handedness,
                :head_circumference=>head_circumference,
                :weight=>weight,
                :height=>height)
end

function _create_recording_eeg(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, reference::String, clabels::Vector{String}, transducers::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64, gain::Vector{Float64})

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :reference=>reference,
                :labels=>clabels,
                :transducers=>transducers,
                :units=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate,
                :gain=>gain,
                :epoch_id=>"")
end

function _create_recording_seeg(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, reference::String, clabels::Vector{String}, transducers::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64, gain::Vector{Float64})

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :reference=>reference,
                :labels=>clabels,
                :transducers=>transducers,
                :units=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate,
                :gain=>gain,
                :epoch_id=>"")
end

function _create_recording_ecog(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, reference::String, clabels::Vector{String}, transducers::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64, gain::Vector{Float64})

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :reference=>reference,
                :labels=>clabels,
                :transducers=>transducers,
                :units=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate,
                :gain=>gain,
                :epoch_id=>"")
end

function _create_recording_meg(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, reference::String, clabels::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64, magnetometers::Vector{Int64}, gradiometers::Vector{Int64}, gradiometers_planar::Vector{Int64}, gradiometers_axial::Vector{Int64}, coils::Vector{Int64})

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>replace(recording_time, '.'=>':'),
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :reference=>reference,
                :labels=>clabels,
                :units=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate,
                :magnetometers=>magnetometers,
                :gradiometers=>gradiometers,
                :gradiometers_planar=>gradiometers_planar,
                :gradiometers_axial=>gradiometers_axial,
                :coils=>coils,
                :epoch_id=>"")
end

function _create_recording_nirs(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, wavelengths::Vector{Float64}, wavelength_index::Vector{Int64}, optode_pairs::Matrix{Int64}, ch_type::Vector{String}, clabels::Vector{String}, units::Vector{String}, src_labels::Vector{String}, det_labels::Vector{String}, opt_labels::Vector{String}, sampling_rate::Int64)

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :wavelengths=>wavelengths,
                :wavelength_index=>wavelength_index,
                :optode_pairs=>optode_pairs,
                :channel_type=>ch_type,
                :labels=>clabels,
                :units=>units,
                :src_labels=>src_labels,
                :det_labels=>det_labels,
                :optode_labels=>opt_labels,
                :sampling_rate=>sampling_rate,
                :epoch_id=>"")
end

function _create_recording_sensors(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, clabels::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64)

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :labels=>clabels,
                :units=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate)
end

function _create_experiment(;name::String, notes::String, design::String)

    return Dict(:name=>name,
                :notes=>notes,
                :design=>design)
end

function _create_header(subject::Dict, recording::Dict, experiment::Dict)

    return NeuroAnalyzer.HEADER(subject,
                                recording,
                                experiment)
end

function _create_recording_mep(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, clabels::Vector{String}, units::Vector{String}, sampling_rate::Int64, stimulation_intensity::Vector{Int64}, coil_type::Vector{String}, stimulation_sample::Vector{Int64}, markers_pos::Vector{Int64}, markers_neg::Vector{Int64})

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :labels=>clabels,
                :units=>units,
                :sampling_rate=>sampling_rate,
                :stimulation_intensity=>stimulation_intensity,
                :coil_type=>coil_type,
                :stimulation_sample=>stimulation_sample,
                :markers_pos=>markers_pos,
                :markers_neg=>markers_neg)
end

function _create_recording_eda(;data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, clabels::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64)

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :labels=>clabels,
                :units=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate,
                :epoch_id=>"")
end
