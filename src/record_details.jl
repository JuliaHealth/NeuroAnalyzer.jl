export sr
export channel_n
export epoch_n
export signal_len
export epoch_len
export history
export labels
export info

"""
    sr(record)

Return sampling rate.

# Arguments

- `record::NeuroAnalyzer.RECORD`

# Returns

- `sr::Int64`
"""
function sr(record::NeuroAnalyzer.RECORD)
    return record.header.recording[:sampling_rate]
end

"""
    channel_n(record; type=:signal)

Return number of channels of `type`.

# Arguments

- `record::NeuroAnalyzer.RECORD`
- `type::Vector{Symbol}=:all`: channel type: `:all`, `:eeg`, `:meg`, `:ecg`, `:eog`, `:emg`, `:ref`, `:mrk`

# Returns

- `channel_n::Int64`
"""
function channel_n(record::NeuroAnalyzer.RECORD; type::Symbol=:all)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")
    length(record.header.recording[:channel_type]) == 0 && throw(ArgumentError("RECORD has no defined channel types."))
    channel_n = 0
    for idx in 1:record.header.recording[:channel_n]
        record.header.recording[:channel_type][idx] == string(type) && (channel_n += 1)
    end
    type === :all && (channel_n = size(record.data, 1))

    return channel_n
end

"""
    epoch_n(eeg)

Return number of epochs.

# Arguments

- `record::NeuroAnalyzer.RECORD`

# Returns

- `epoch_n::Int64`
"""
function epoch_n(record::NeuroAnalyzer.RECORD)
    ndims(record.data) < 3 && throw(ArgumentError("Record data is either a vector or a matrix."))
    return record.header.recording[:epoch_n]
end

"""
    signal_len(eeg)

Return signal length.

# Arguments

- `record::NeuroAnalyzer.RECORD`

# Returns

- `signal_len::Int64`
"""
function signal_len(record::NeuroAnalyzer.RECORD)
    return record.header.recording[:duration_samples]
end

"""
    epoch_len(eeg)

Return epoch length.

# Arguments

- `record::NeuroAnalyzer.RECORD`

# Returns

- `epoch_len::Int64`
"""
function epoch_len(record::NeuroAnalyzer.RECORD)
    return record.header.recording[:epoch_duration_samples]
end

"""
    signal_channels(eeg)

Return all signal (e.g. EEG or MEG) channels; signal is determined by `:data_type` variable in `record.header.recording`).

# Arguments

- `record::NeuroAnalyzer.RECORD`:

# Returns
 
- `channels::Vector{Int64}`
"""
function signal_channels(record::NeuroAnalyzer.RECORD)
    return get_channel_bytype(record, type=Symbol(record.header.recording[:data_type]))
end

"""
    get_channel_bytype(record; type=:eeg)

Return channel number(s) for channel of `type` type.

# Arguments

- `record::NeuroAnalyzer.RECORD`
- `type::Vector{Symbol}=:all`: channel type

# Returns

- `ch_n::Vector{Int64}`
"""
function get_channel_bytype(record::NeuroAnalyzer.RECORD; type::Symbol=:all)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")
    type === :all && return collect(1:channel_n(record))
    ch_idx = Vector{Int64}()
    for idx in 1:channel_n(record)
        record.header.recording[:channel_type][idx] == string(type) && (push!(ch_idx, idx))
    end

    return ch_idx
end

"""
    history(eeg)

Show processing history.

# Arguments

- `record::NeuroAnalyzer.RECORD`

# Returns

- `history::Vector{String}`
"""
function history(record::NeuroAnalyzer.RECORD)
    return record.header.history
end

"""
    labels(record)

Return channel labels.

# Arguments

- `record::NeuroAnalyzer.RECORD`

# Returns

- `labels::Vector{String}`
"""
function labels(record::NeuroAnalyzer.RECORD)
    length(record.header.recording[:labels]) == 0 && throw(ArgumentError("RECORD has no labels."))
    return record.header.recording[:labels]
end

"""
    info(eeg)

Show info.

# Arguments

- `record::NeuroAnalyzer.RECORD`
"""
function info(record::NeuroAnalyzer.RECORD)

    println("            Signal type: $(uppercase(record.header.recording[:data_type]))")
    println("            File format: $(record.header.recording[:file_type])")
    println("            Source file: $(record.header.recording[:file_name])")
    println("                Subject: $(record.header.subject[:first_name] * " " * record.header.subject[:last_name])")
    println("              Recording: $(record.header.recording[:recording])")
    println("        Recording notes: $(record.header.recording[:recording_notes])")
    println("         Recording date: $(record.header.recording[:recording_date])")
    println("         Recording time: $(record.header.recording[:recording_time])")
    println("     Sampling rate (Hz): $(sr(record))")
    println("Signal length [samples]: $(signal_len(record))")
    println("Signal length [seconds]: $(round(record.header.recording[:duration_seconds], digits=2))")
    println("     Number of channels: $(channel_n(record))")
    println("       Number of epochs: $(epoch_n(record))")
    println(" Epoch length [samples]: $(epoch_len(record))")
    println(" Epoch length [seconds]: $(round(record.header.recording[:epoch_duration_seconds], digits=2))")
    println("         File size [MB]: $(record.header.recording[:file_size_mb])")
    println("       Memory size [MB]: $(round(Base.summarysize(record) / 1024^2, digits=2))")
    if record.header.recording[:reference] == ""
        println("         Reference type: unknown")
    else
        println("         Reference type: $(record.header[:reference])")
    end
    if length(labels(record)) == 0
        println("                 Labels: no")
    else
        println("                 Labels: yes")
    end
    if record.header.markers == false
        println("                Markers: no")
    else
        println("                Markers: yes")
    end
    if record.header.locations == false
        println("      Channel locations: no")
    else
        println("      Channel locations: yes")
    end
    if record.header.components != []
        print("             Components: ")
        c = list_components(record)
        if length(c) == 1
            println(c[1])
        else
            for idx in 1:(length(c) - 1)
                print(c[idx], ", ")
            end
            println(c[end])
        end
    else
        println("             Components: no")
    end
    println("Channels:")
    for idx in eachindex(record.header.recording[:labels])
        println("\tchannel: $idx\tlabel: $(rpad(record.header.recording[:labels][idx], 16, " "))\ttype: $(uppercase(record.header.recording[:channel_type][idx]))")
    end
end