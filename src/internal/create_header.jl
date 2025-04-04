function _create_subject(; id::String, first_name::String, middle_name::String, last_name::String, handedness::String, head_circumference::Real, weight::Real, height::Real)::Dict

    return Dict(:id=>id,
                :first_name=>first_name,
                :middle_name=>middle_name,
                :last_name=>last_name,
                :handedness=>handedness,
                :head_circumference=>head_circumference,
                :weight=>weight,
                :height=>height)

end

function _create_recording_eeg(; data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, channel_order::Vector{Int64}, reference::String, clabels::Vector{String}, transducers::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, line_frequency::Real, sampling_rate::Int64, gain::Vector{Float64}, bad_channels::Matrix{Bool})::Dict

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :channel_order=>channel_order,
                :reference=>reference,
                :label=>clabels,
                :transducers=>transducers,
                :unit=>units,
                :prefiltering=>prefiltering,
                :line_frequency=>line_frequency,
                :sampling_rate=>sampling_rate,
                :gain=>gain,
                :bad_channel=>bad_channels,
                :epoch_id=>"")

end

function _create_recording_seeg(; data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, channel_order::Vector{Int64}, reference::String, clabels::Vector{String}, transducers::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, line_frequency::Real, sampling_rate::Int64, gain::Vector{Float64}, bad_channels::Matrix{Bool})::Dict

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :channel_order=>channel_order,
                :reference=>reference,
                :label=>clabels,
                :transducers=>transducers,
                :unit=>units,
                :prefiltering=>prefiltering,
                :line_frequency=>line_frequency,
                :sampling_rate=>sampling_rate,
                :gain=>gain,
                :bad_channel=>bad_channels,
                :epoch_id=>"")

end

function _create_recording_ecog(; data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, channel_order::Vector{Int64}, reference::String, clabels::Vector{String}, transducers::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, line_frequency::Real, sampling_rate::Int64, gain::Vector{Float64}, bad_channels::Matrix{Bool})::Dict

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :channel_order=>channel_order,
                :reference=>reference,
                :label=>clabels,
                :transducers=>transducers,
                :unit=>units,
                :prefiltering=>prefiltering,
                :line_frequency=>line_frequency,
                :sampling_rate=>sampling_rate,
                :gain=>gain,
                :bad_channel=>bad_channels,
                :epoch_id=>"")

end

function _create_recording_meg(; data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, channel_order::Vector{Int64}, reference::String, clabels::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, line_frequency::Real, sampling_rate::Int64, magnetometers::Vector{Int64}, gradiometers::Vector{Int64}, coil_type::Vector{String}, bad_channels::Matrix{Bool}, ssp_labels::Vector{String}, ssp_channels::Vector{Bool}, ssp_data::Matrix{Float64})::Dict

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>replace(recording_time, '.'=>':'),
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :channel_order=>channel_order,
                :reference=>reference,
                :label=>clabels,
                :unit=>units,
                :prefiltering=>prefiltering,
                :line_frequency=>line_frequency,
                :sampling_rate=>sampling_rate,
                :magnetometers=>magnetometers,
                :gradiometers=>gradiometers,
                :coil_type=>coil_type,
                :bad_channel=>bad_channels,
                :ssp_labels=>ssp_labels,
                :ssp_channels=>ssp_channels,
                :ssp_data=>ssp_data,
                :epoch_id=>"")

end

function _create_recording_nirs(; data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, wavelengths::Vector{Float64}, wavelength_index::Vector{Int64}, optode_pairs::Matrix{Int64}, channel_type::Vector{String}, channel_order::Vector{Int64}, clabels::Vector{String}, units::Vector{String}, src_labels::Vector{String}, det_labels::Vector{String}, opt_labels::Vector{String}, sampling_rate::Int64, bad_channels::Matrix{Bool})::Dict

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
                :channel_type=>channel_type,
                :channel_order=>channel_order,
                :label=>clabels,
                :unit=>units,
                :src_labels=>src_labels,
                :det_labels=>det_labels,
                :optode_labels=>opt_labels,
                :sampling_rate=>sampling_rate,
                :bad_channel=>bad_channels,
                :epoch_id=>"")

end

function _create_recording_sensors(; data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, channel_order::Vector{Int64}, clabels::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64, bad_channels::Matrix{Bool})::Dict

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :channel_order=>channel_order,
                :label=>clabels,
                :unit=>units,
                :prefiltering=>prefiltering,
                :bad_channel=>bad_channels,
                :sampling_rate=>sampling_rate,
                :epoch_id=>"")

end

function _create_experiment(; name::String, notes::String, design::String)::Dict

    return Dict(:name=>name,
                :notes=>notes,
                :design=>design)

end

function _create_header(subject::Dict, recording::Dict, experiment::Dict)::NeuroAnalyzer.HEADER

    return NeuroAnalyzer.HEADER(subject,
                                recording,
                                experiment)

end

function _create_recording_mep(; data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, channel_order::Vector{Int64}, clabels::Vector{String}, units::Vector{String}, sampling_rate::Int64, stimulation_intensity::Vector{Int64}, coil_type::Vector{String}, stimulation_sample::Vector{Int64}, markers_pos::Vector{Int64}, markers_neg::Vector{Int64}, bad_channels::Matrix{Bool})::Dict

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :channel_order=>channel_order,
                :label=>clabels,
                :unit=>units,
                :sampling_rate=>sampling_rate,
                :stimulation_intensity=>stimulation_intensity,
                :coil_type=>coil_type,
                :stimulation_sample=>stimulation_sample,
                :markers_pos=>markers_pos,
                :markers_neg=>markers_neg,
                :bad_channel=>bad_channels)

end

function _create_recording_eda(; data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, channel_order::Vector{Int64}, clabels::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64, bad_channels::Matrix{Bool})::Dict

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :channel_order=>channel_order,
                :label=>clabels,
                :unit=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate,
                :epoch_id=>"",
                :bad_channel=>bad_channels)

end

function _create_recording_tpt(; data_type::String, file_name::String, file_size_mb::Real, file_type::String, recording::String, recording_date::String, recording_time::String, recording_notes::String, channel_type::Vector{String}, channel_order::Vector{Int64}, clabels::Vector{String}, units::Vector{String}, prefiltering::Vector{String}, sampling_rate::Int64, bad_channels::Matrix{Bool})::Dict

    return Dict(:data_type=>data_type,
                :file_name=>file_name,
                :file_size_mb=>file_size_mb,
                :file_type=>file_type,
                :recording=>recording,
                :recording_date=>recording_date,
                :recording_time=>recording_time,
                :recording_notes=>recording_notes,
                :channel_type=>channel_type,
                :channel_order=>channel_order,
                :label=>clabels,
                :unit=>units,
                :prefiltering=>prefiltering,
                :sampling_rate=>sampling_rate,
                :epoch_id=>"",
                :bad_channel=>bad_channels)

end
