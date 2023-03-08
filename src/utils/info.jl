export sr
export channel_n
export epoch_n
export signal_len
export epoch_len
export signal_channels
export get_channel_bytype
export history
export labels
export info

function _info(s::String)
    verbose == true && @info s
end

function _warn(s::String)
    verbose == true && @warn s
end

"""
    sr(obj)

Return sampling rate.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `sr::Int64`
"""
function sr(obj::NeuroAnalyzer.NEURO)
    return obj.header.recording[:sampling_rate]
end

"""
    channel_n(obj; type=:signal)

Return number of channels of `type`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Vector{Symbol}=:all`: channel type: `:all`, `:eeg`, `:meg`, `:ecg`, `:eog`, `:emg`, `:ref`, `:mrk`

# Returns

- `channel_n::Int64`
"""
function channel_n(obj::NeuroAnalyzer.NEURO; type::Symbol=:all)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")
    length(obj.header.recording[:channel_type]) == 0 && throw(ArgumentError("RECORD has no defined channel types."))
    channel_n = 0
    for idx in 1:obj.header.recording[:channel_n]
        obj.header.recording[:channel_type][idx] == string(type) && (channel_n += 1)
    end
    type === :all && (channel_n = size(obj.data, 1))

    return channel_n
end

"""
    epoch_n(eeg)

Return number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `epoch_n::Int64`
"""
function epoch_n(obj::NeuroAnalyzer.NEURO)
    ndims(obj.data) < 3 && throw(ArgumentError("Record data is either a vector or a matrix."))
    return obj.header.recording[:epoch_n]
end

"""
    signal_len(eeg)

Return signal length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `signal_len::Int64`
"""
function signal_len(obj::NeuroAnalyzer.NEURO)
    return obj.header.recording[:duration_samples]
end

"""
    epoch_len(eeg)

Return epoch length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `epoch_len::Int64`
"""
function epoch_len(obj::NeuroAnalyzer.NEURO)

    ndims(obj.data) < 3 && throw(ArgumentError("Record data is either a vector or a matrix."))

    return obj.header.recording[:epoch_duration_samples]
end

"""
    signal_channels(eeg)

Return all signal (e.g. EEG or MEG) channels; signal is determined by `:data_type` variable in `obj.header.recording`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`:

# Returns
 
- `channels::Vector{Int64}`
"""
function signal_channels(obj::NeuroAnalyzer.NEURO)
    return get_channel_bytype(obj, type=Symbol(obj.header.recording[:data_type]))
end

"""
    get_channel_bytype(obj; type=:eeg)

Return channel number(s) for channel of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Vector{Symbol}=:all`: channel type

# Returns

- `ch_n::Vector{Int64}`
"""
function get_channel_bytype(obj::NeuroAnalyzer.NEURO; type::Symbol=:all)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")
    type === :all && return collect(1:channel_n(obj))
    ch_idx = Vector{Int64}()
    for idx in 1:channel_n(obj)
        obj.header.recording[:channel_type][idx] == string(type) && (push!(ch_idx, idx))
    end

    return ch_idx
end

"""
    history(eeg)

Show processing history.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `history::Vector{String}`
"""
function history(obj::NeuroAnalyzer.NEURO)
    return obj.header.history
end

"""
    labels(obj)

Return channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `labels::Vector{String}`
"""
function labels(obj::NeuroAnalyzer.NEURO)
    length(obj.header.recording[:labels]) == 0 && throw(ArgumentError("RECORD has no labels."))

    return obj.header.recording[:labels]
end

"""
    info(eeg)

Show info.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function info(obj::NeuroAnalyzer.NEURO)

    println("            Signal type: $(uppercase(obj.header.recording[:data_type]))")
    println("            File format: $(obj.header.recording[:file_type])")
    println("            Source file: $(obj.header.recording[:file_name])")
    println("                Subject: $(obj.header.subject[:first_name] * " " * obj.header.subject[:last_name])")
    println("              Recording: $(obj.header.recording[:recording])")
    println("        Recording notes: $(obj.header.recording[:recording_notes])")
    println("         Recording date: $(obj.header.recording[:recording_date])")
    println("         Recording time: $(obj.header.recording[:recording_time])")
    println("     Sampling rate (Hz): $(sr(obj))")
    println("Signal length [samples]: $(signal_len(obj))")
    println("Signal length [seconds]: $(round(obj.header.recording[:duration_seconds], digits=2))")
    println("     Number of channels: $(channel_n(obj))")
    println("       Number of epochs: $(epoch_n(obj))")
    println(" Epoch length [samples]: $(epoch_len(obj))")
    println(" Epoch length [seconds]: $(round(obj.header.recording[:epoch_duration_seconds], digits=2))")
    println("         File size [MB]: $(obj.header.recording[:file_size_mb])")
    println("       Memory size [MB]: $(round(Base.summarysize(obj) / 1024^2, digits=2))")
    if obj.header.recording[:reference] == ""
        println("         Reference type: unknown")
    else
        println("         Reference type: $(obj.header[:reference])")
    end
    if length(labels(obj)) == 0
        println("                 Labels: no")
    else
        println("                 Labels: yes")
    end
    if obj.header.markers == false
        println("                Markers: no")
    else
        println("                Markers: yes")
    end
    if obj.header.locations == false
        println("      Channel locations: no")
    else
        println("      Channel locations: yes")
    end
    if obj.header.components != []
        print("             Components: ")
        c = list_components(obj)
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
    for idx in eachindex(obj.header.recording[:labels])
        println("\tchannel: $idx\tlabel: $(rpad(obj.header.recording[:labels][idx], 16, " "))\ttype: $(uppercase(obj.header.recording[:channel_type][idx]))")
    end
end