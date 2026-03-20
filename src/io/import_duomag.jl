export import_duomag

"""
    import_duomag(file_name)

Load a DuoMAG TMS MEP recording file (`.ascii` or `.m`) and return a
`NeuroAnalyzer.NEURO` object.

Both file formats carry metadata (subject, recording info, stimulation parameters), per-channel MEP signals, and positive/negative peak markers. Signal data are baseline-corrected, stimulation artifacts suppressed, and
amplitudes scaled to μV on import.

# Arguments

- `file_name::String`: path to the `.ascii` or `.m` file

# Returns

- `NeuroAnalyzer.NEURO`

# Throws
- `ArgumentError` if the file does not exist or has an unsupported extension
"""
function import_duomag(file_name::String)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    # extract extension once - used in every branch below.
    ext = splitext(file_name)[2]
    ext in (".ascii", ".m") ||
        throw(ArgumentError("$file_name is not a DuoMAG file (expected .ascii or .m)."))

    # ------------------------------------------------------------------ #
    # parse format-specific header and signal data                       #
    # both branches use open(...) do to guarantee close on exception     #
    # ------------------------------------------------------------------ #

    if ext == ".ascii"
        sw_version, subject, subject_id, record_id, record_created,
        sampling_interval, sampling_interval_unit,
        sensitivity, sensitivity_unit,
        signal_count, samples_count,
        stim_sample, stim_intens, coil_type,
        markers_pos, markers_neg,
        mep_signal = open(file_name, "r") do f

            readline(f) # software name (unused)
            readline(f) # blank
            sw_version = split(strip(readline(f)), '=')[2]
            subject = split(strip(readline(f)), '=')[2]
            subject_id = split(strip(readline(f)), '=')[2]
            record_id = split(strip(readline(f)), '=')[2]
            record_created = split(strip(readline(f)), '=')[2]
            readline(f) # method (unused)
            readline(f)
            readline(f) # marker latency unit (unused)
            sampling_interval = parse(Int, split(strip(readline(f)), '=')[2])
            sampling_interval_unit = replace(split(strip(readline(f)), '=')[2], "\xb5" => "μ")
            sensitivity = parse(Float64,
                replace(split(strip(readline(f)), '=')[2], ',' => '.'))
            sensitivity_unit = replace(split(strip(readline(f)), '=')[2], "\xb5" => "μ")
            signal_count = parse(Int, split(strip(readline(f)), '=')[2])
            readline(f) # reserved field
            readline(f) # reserved field
            samples_count = parse.(Int, split(strip(readline(f)), ' ')[2:end])
            stim_sample = parse.(Int, split(strip(readline(f)), ' ')[2:end])
            stim_intens = parse.(Int, split(strip(readline(f)), ' ')[2:end])
            coil_type = split(strip(readline(f)), ' ')[2:end]
            readline(f)

            # positive and negative peak marker positions (ms, may be "N/A")
            parse_marker_line(l) = begin
                raw = replace.(split(strip(l), ' ')[2:end], ',' => '.')
                parse.(Float64, replace.(raw, "N/A" => "0"))
            end
            markers_pos = parse_marker_line(readline(f))
            markers_neg = parse_marker_line(readline(f))
            readline(f)

            # signal matrix: rows = samples, columns = signals
            mep_signal = zeros(samples_count[1], signal_count)
            for idx in axes(mep_signal, 1)
                mep_signal[idx, :] = parse.(Float64,
                    replace.(split(strip(readline(f)), ' '), ',' => '.'))
            end

            sw_version, subject, subject_id, record_id, record_created,
            sampling_interval, sampling_interval_unit,
            sensitivity, sensitivity_unit,
            signal_count, samples_count,
            stim_sample, stim_intens, coil_type,
            markers_pos, markers_neg, mep_signal
        end
        # file closed here

    elseif ext == ".m"
        sw_version, subject, subject_id, record_id, record_created,
        sampling_interval, sampling_interval_unit,
        sensitivity, sensitivity_unit,
        signal_count, samples_count,
        stim_sample, stim_intens, coil_type,
        markers_pos, markers_neg,
        mep_signal = open(file_name, "r") do f

            readline(f)
            record_info = readline(f)
            readline(f)
            signals_info = readline(f)
            readline(f)
            stimulations_info = readline(f)
            readline(f)
            signals_data = readline(f)
            readline(f)
            markers_raw = readline(f)

            # -------------------------------------------------------- #
            # helper: strip MATLAB struct wrapper and split on comma   #
            # -------------------------------------------------------- #
            parse_struct(s, prefix) = begin
                s = replace(s, prefix => "", ");" => "", "\xb5" => "μ",
                              "'" => "", "{" => "", "}" => "", ";" => " ")
                split(s, ',')
            end

            ri = parse_struct(record_info, "RecordInfo = struct(")
            record_id = ri[2]
            subject = ri[4]
            subject_id = ri[6]
            record_created = ri[8]
            sw_version = ri[10]
            signal_count = parse(Int64, ri[16])
            samples_count = [parse(Int64, ri[18])]
            sampling_interval = parse(Int64, ri[20])
            sampling_interval_unit = ri[22]

            si = parse_struct(signals_info, "SignalsInfo = struct(")
            popfirst!(si)
            popfirst!(si)
            sensitivity = parse(Float64, split(si[2], ' ')[1])
            sensitivity_unit = split(si[2], ' ')[2]

            sti = parse_struct(stimulations_info, "StimulationInfo = struct(")
            popfirst!(sti); popfirst!(sti)
            stim_intens = Int64[]
            coil_type = String[]
            for idx in 2:2:(2 * signal_count)
                push!(stim_intens, parse(Int64, split(sti[idx], ' ')[1]))
                push!(coil_type, split(sti[idx], ' ')[2])
            end

            sd = replace(signals_data, "SignalsData = [" => "", "];" => "")
            sd_rows = split(sd, ';')
            mep_signal = zeros(length(sd_rows), length(split(sd_rows[1], ' ')))
            for idx in eachindex(sd_rows)
                mep_signal[idx, :] = parse.(Float64, split(sd_rows[idx], ' '))
            end

            mk = parse_struct(markers_raw, "Markers = struct(")
            popfirst!(mk)
            popfirst!(mk)

            stim_sample = zeros(Int64, signal_count)
            markers_neg = zeros(Int64, signal_count) # negative (A-) peak
            markers_pos = zeros(Int64, signal_count) # positive (A+) peak

            for idx in 2:2:length(mk)
                tmp = split(mk[idx], ' ')
                length(tmp) < 3 && continue
                stim_number = parse(Int64, tmp[3])
                tmp[1] == "Stim" && (stim_sample[stim_number] = parse(Int64, tmp[2]))
                tmp[1] == "A+" && (markers_pos[stim_number] = parse(Int64, tmp[2]))
                tmp[1] == "A-" && (markers_neg[stim_number] = parse(Int64, tmp[2]))
            end

            sw_version, subject, subject_id, record_id, record_created,
            sampling_interval, sampling_interval_unit,
            sensitivity, sensitivity_unit,
            signal_count, samples_count,
            stim_sample, stim_intens, coil_type,
            markers_pos, markers_neg, mep_signal
        end
        # file closed here

    end

    # ------------------------------------------------------------------ #
    # signal processing: transpose, invert polarity, suppress artifact,  #
    # remove DC, apply sensitivity gain                                  #
    # ------------------------------------------------------------------ #
    data = zeros(size(mep_signal, 2), size(mep_signal, 1), 1)
    for idx in axes(data, 1)
        data[idx, :, 1] = @views -mep_signal[:, idx]
        # suppress stimulation artefact (±10 samples around stim sample)
        data[idx, (stim_sample[idx] - 10):(stim_sample[idx] + 10), 1] .*= 0.05
        # remove DC offset relative to pre-stim baseline
        data[idx, :] = remove_dc(data[idx, :], stim_sample[idx] - 10)
    end
    data .*= sensitivity

    # ------------------------------------------------------------------ #
    # unit conversion to μV                                              #
    # ------------------------------------------------------------------ #
    sensitivity_unit == "mV" && (data .*= 1e3)
    sensitivity_unit == "V"  && (data .*= 1e6)

    # sampling interval to sampling rate
    sampling_interval_unit == "μs" && (sampling_interval *= 1e-6)
    sampling_interval_unit == "ms" && (sampling_interval *= 1e-3)
    sampling_rate = round(Int64, 1 / sampling_interval)

    # ------------------------------------------------------------------ #
    # channel labels                                                      #
    # ------------------------------------------------------------------ #
    ch_n    = signal_count
    clabels = ["MEP$idx" for idx in 1:ch_n]

    # ------------------------------------------------------------------ #
    # time axis - zero-aligned to the stimulation sample                  #
    # ------------------------------------------------------------------ #
    n_samples  = size(data, 2) * size(data, 3)
    time_pts   = round.(
        range(0; step = 1/sampling_rate, length = n_samples) .- (stim_sample[1] / sampling_rate);
        digits = 4)
    epoch_time = time_pts

    # convert .ascii marker positions (ms) to sample indices
    if ext == ".ascii"
        for idx in eachindex(markers_pos)
            markers_pos[idx] = vsearch(markers_pos[idx] / 1000, time_pts)
            markers_pos[idx] == stim_sample[1] && (markers_pos[idx] = 0)
            markers_neg[idx] = vsearch(markers_neg[idx] / 1000, time_pts)
            markers_neg[idx] == stim_sample[1] && (markers_neg[idx] = 0)
        end
        markers_pos = Int64.(markers_pos)
        markers_neg = Int64.(markers_neg)
    end

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                               #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(
        id = string(subject_id),
        first_name = "",
        middle_name = "",
        last_name = string(subject),
        head_circumference = -1,
        handedness = "",
        weight = -1,
        height = -1)
    r = _create_recording_mep(
        data_type = "mep",
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = "DuoMAG",
        recording = "",
        recording_date = string(split(record_created, ' ')[1]),
        recording_time = replace(string(split(record_created, ' ')[2]), '.' => ':'),
        recording_notes = "",
        channel_type = repeat(["mep"], ch_n),
        channel_order = _sort_channels(repeat(["mep"], ch_n)),
        clabels = clabels,
        units = repeat(["μV"], ch_n),
        sampling_rate = sampling_rate,
        stimulation_intensity = stim_intens,
        coil_type = string.(coil_type),
        stimulation_sample = stim_sample,
        markers_pos = markers_pos,
        markers_neg = markers_neg,
        bad_channels = zeros(Bool, size(data, 1)))
    e   = _create_experiment(name = "", notes = "", design = "")
    hdr = _create_header(subject = s, recording = r, experiment = e)

    markers_df = DataFrame(
        :id => String[],
        :start => Float64[],
        :length => Float64[],
        :value => String[],
        :channel => Int64[])

    locs = _initialize_locs()
    obj  = NeuroAnalyzer.NEURO(hdr, String[], markers_df, locs, time_pts, epoch_time, data)

    _info("Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end]; digits=2)) s)")

    return obj

end