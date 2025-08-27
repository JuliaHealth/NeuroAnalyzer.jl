export import_duomag

"""
    import_duomag(file_name)

Load DuoMAG TMS MEP recording file (.ascii or .m) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_duomag(file_name::String)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert (splitext(file_name)[2] == ".ascii" || splitext(file_name)[2] == ".m") "This is not DuoMAG file."

    f = nothing

    if splitext(file_name)[2] == ".ascii"
        try
            f = open(file_name, "r")
        catch
            error("File $file_name cannot be opened!")
        end
        _ = readline(f)
        # sw_name = split(strip(readline(f)), '=')[2]
        _ = readline(f)
        sw_version = split(strip(readline(f)), '=')[2]
        subject = split(strip(readline(f)), '=')[2]
        subject_id = split(strip(readline(f)), '=')[2]
        record_id = split(strip(readline(f)), '=')[2]
        record_created = split(strip(readline(f)), '=')[2]
        _ = readline(f)
        # method = split(strip(readline(f)), '=')[2]
        _ = readline(f)
        _ = readline(f)
        # marker_latency_unit = split(strip(readline(f)), '=')[2]
        _ = readline(f)
        sampling_interval = parse(Int, (split(strip(readline(f)), '=')[2]))
        sampling_interval_unit = split(strip(readline(f)), '=')[2]
        sampling_interval_unit = replace(sampling_interval_unit, "\xb5"=>"μ")
        sensitivity = parse(Float64, replace(split(strip(readline(f)), '=')[2], ',' => '.'))
        sensitivity_unit = split(strip(readline(f)), '=')[2]
        sensitivity_unit = replace(sensitivity_unit, "\xb5"=>"μ")
        signal_count = parse(Int, split(strip(readline(f)), '=')[2])
        _ = parse.(Int, split(strip(readline(f)), ' ')[2:end])
        _ = parse.(Int, split(strip(readline(f)), ' ')[2:end])
        samples_count = parse.(Int, split(strip(readline(f)), ' ')[2:end])
        stim_sample = parse.(Int, split(strip(readline(f)), ' ')[2:end])
        stim_intens = parse.(Int, split(strip(readline(f)), ' ')[2:end])
        coil_type = split(strip(readline(f)), ' ')[2:end]
        _ = readline(f)
        markers_pos = replace.(split(strip(readline(f)), ' ')[2:end], ',' => '.')
        markers_pos = replace.(markers_pos, "N/A" => "0")
        markers_pos = parse.(Float64, markers_pos)
        markers_neg = replace.(split(strip(readline(f)), ' ')[2:end], ',' => '.')
        markers_neg = replace.(markers_neg, "N/A" => "0")
        markers_neg = parse.(Float64, markers_neg)
        _ = readline(f)

        # data matrix: signals × samples
        mep_signal = zeros((samples_count[1], signal_count))
        for idx in axes(mep_signal)[1]
            mep_signal[idx, :] = parse.(Float64, replace.(split(strip(readline(f)), ' '), ',' => '.'))
        end

        close(f)
    elseif splitext(file_name)[2] == ".m"
        try
            f = open(file_name, "r")
        catch
            error("File $file_name cannot be opened!")
        end

        _ = readline(f)
        record_info = readline(f)
        _ = readline(f)
        signals_info = readline(f)
        _ = readline(f)
        stimulations_info = readline(f)
        _ = readline(f)
        signals_data = readline(f)
        _ = readline(f)
        markers = readline(f)
        close(f)

        record_info = replace(record_info, "RecordInfo = struct("=>"")
        record_info = replace(record_info, ");"=>"")
        record_info = replace(record_info, "\xb5"=>"μ")
        record_info = replace(record_info, "'"=>"")
        record_info = replace(record_info, "{"=>"")
        record_info = replace(record_info, "}"=>"")
        record_info = split(record_info, ',')
        record_id = record_info[2]
        subject = record_info[4]
        subject_id = record_info[6]
        record_created = record_info[8]
        sw_version = record_info[10]
        signal_count = parse(Int64, record_info[16])
        samples_count = parse(Int64, record_info[18])
        sampling_interval = parse(Int64, record_info[20])
        sampling_interval_unit = record_info[22]

        signals_info = replace(signals_info, "SignalsInfo = struct("=>"")
        signals_info = replace(signals_info, ");"=>"")
        signals_info = replace(signals_info, "\xb5"=>"μ")
        signals_info = replace(signals_info, "'"=>"")
        signals_info = replace(signals_info, "{"=>"")
        signals_info = replace(signals_info, "}"=>"")
        signals_info = replace(signals_info, ";"=>" ")
        signals_info = split(signals_info, ',')
        popfirst!(signals_info)
        popfirst!(signals_info)
        sensitivity = parse(Float64, split(signals_info[2], ' ')[1])
        sensitivity_unit = split(signals_info[2], ' ')[2]

        stimulations_info = replace(stimulations_info, "StimulationInfo = struct("=>"")
        stimulations_info = replace(stimulations_info, ");"=>"")
        stimulations_info = replace(stimulations_info, "\xb5"=>"μ")
        stimulations_info = replace(stimulations_info, "'"=>"")
        stimulations_info = replace(stimulations_info, "{"=>"")
        stimulations_info = replace(stimulations_info, "}"=>"")
        stimulations_info = replace(stimulations_info, ";"=>" ")
        stimulations_info = split(stimulations_info, ',')
        popfirst!(stimulations_info)
        popfirst!(stimulations_info)
        stim_intens = Int64[]
        coil_type = String[]
        for idx in 2:2:(2 * signal_count)
            push!(stim_intens, parse(Int64, split(stimulations_info[idx], ' ')[1]))
            push!(coil_type, split(stimulations_info[idx], ' ')[2])
        end

        signals_data = replace(signals_data, "SignalsData = ["=>"")
        signals_data = replace(signals_data, "];"=>"")
        signals_data = split(signals_data, ';')
        mep_signal = zeros(length(signals_data), length(split(signals_data[1], ' ')))
        for idx in eachindex(signals_data)
            mep_signal[idx, :] = parse.(Float64, split(signals_data[idx], ' '))
        end

        markers = replace(markers, "Markers = struct("=>"")
        markers = replace(markers, ");"=>"")
        markers = replace(markers, "\xb5"=>"μ")
        markers = replace(markers, "'"=>"")
        markers = replace(markers, "{"=>"")
        markers = replace(markers, "}"=>"")
        markers = replace(markers, ";"=>" ")
        markers = split(markers, ',')
        popfirst!(markers)
        popfirst!(markers)
        replace.(markers, '{'=>"")
        replace.(markers, '}'=>"")

        stim_sample = zeros(Int64, signal_count)
        markers_neg = zeros(Int64, signal_count)
        markers_pos = zeros(Int64, signal_count)
        for idx in 2:2:length(markers)
            tmp = split(markers[idx], ' ')
            stim_number = parse(Int64, tmp[3])
            tmp[1] == "Stim" && (stim_sample[stim_number] = parse(Int64, tmp[2]))
            tmp[1] == "A-" && (markers_pos[stim_number] = parse(Int64, tmp[2]))
            tmp[1] == "A+" && (markers_neg[stim_number] = parse(Int64, tmp[2]))
        end
    end

    data = zeros(size(mep_signal, 2), size(mep_signal, 1), 1)
    for idx in axes(data, 1)
        data[idx, :, 1] = @views -mep_signal[:, idx]
        # reduce stimulation amplitude
        data[idx, stim_sample[1]-10:stim_sample[1]+10, 1] .*= 0.05
        # remove DC offset
        data[idx, :] = remove_dc(data[idx, :], stim_sample[idx]-10)
    end
    data .*= sensitivity

    ch_n = signal_count
    clabels = repeat(["MEP"], ch_n)
    for idx in 1:ch_n
        clabels[idx] *= string(idx)
    end

    # convert to mV
    sensitivity_unit == "mV" && (data ./= 10^3)
    sensitivity_unit == "V" && (data .*= 10^3)

    sampling_interval_unit == "μs" && (sampling_interval *= 10^-6)
    sampling_interval_unit == "ms" && (sampling_interval *= 10^-3)
    sampling_rate = round(Int64, 1 / sampling_interval)

    time_pts = collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1]
    time_pts = round.(time_pts .- time_pts[stim_sample[1]], digits=4)
    ep_time = time_pts

    if splitext(file_name)[2] == ".ascii"
        for idx in eachindex(markers_pos)
            markers_pos[idx] = vsearch(markers_pos[idx] / 1000, time_pts)
            markers_pos[idx] == stim_sample[1] && (markers_pos[idx] = 0)
            markers_neg[idx] = vsearch(markers_neg[idx] / 1000, time_pts)
            markers_neg[idx] == stim_sample[1] && (markers_neg[idx] = 0)
        end
        markers_pos = Int64.(markers_pos)
        markers_neg = Int64.(markers_neg)
    end

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    s = _create_subject(id=string(subject_id),
                        first_name="",
                        middle_name="",
                        last_name=string(subject),
                        head_circumference=-1,
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_mep(data_type="mep",
                              file_name=file_name,
                              file_size_mb=0,
                              file_type="DuoMAG",
                              recording="",
                              recording_date=string(split(record_created, ' ')[1]),
                              recording_time=replace(string(split(record_created, ' ')[2]), '.'=>':'),
                              recording_notes="",
                              channel_type=repeat(["mep"], ch_n),
                              channel_order=_sort_channels(repeat(["mep"], ch_n)),
                              clabels=clabels,
                              units=repeat(["μV"], ch_n),
                              sampling_rate=sampling_rate,
                              stimulation_intensity=stim_intens,
                              coil_type=string.(coil_type),
                              stimulation_sample=stim_sample,
                              markers_pos=markers_pos,
                              markers_neg=markers_neg,
                              bad_channels=zeros(Bool, size(data, 1), 1))
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :value=>String[],
                        :channel=>Int64[])

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data, components, markers, locs, history)
    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=3)) s)")

    return obj

end
