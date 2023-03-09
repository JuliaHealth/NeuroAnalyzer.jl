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
export channel_cluster
export band_frq

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

"""
    channels_cluster(obj, cluster)

Return channels belonging to a `cluster` of channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `cluster::Symbol`: available clusters are:
    - `:f1`: left frontal (F1, F3, F5, F7, AF3, AF7)
    - `:f2`: right frontal (F2, F4, F6, F8, AF4, AF8)
    - `:t1`: left temporal (C3, C5, T7, FC3, FC5, FT7)
    - `:t2`: right temporal (C4, C6, T8, FC4, FC6, FT8)
    - `:c1`: anterior central (Cz, C1, C2, FC1, FC2, FCz)
    - `:c2`: posterior central (Pz, P1, P2, CP1, CP2, CPz)
    - `:p1`: left parietal (P3, P5, P7, CP3, CP5, TP7)
    - `:p2`: right parietal (P4, P6, P8, CP4, CP6, TP8)
    - `:o`: occipital (Oz, O1, O2, POz, PO3, PO4)

# Returns

- `channels::Vector{Int64}`: list of channel numbers belonging to a given cluster of channels
"""
function channel_cluster(obj::NeuroAnalyzer.NEURO; cluster::Symbol)

    length(labels(obj)) == 0 && throw(ArgumentError("OBJ does not contain channel labels."))

    _check_var(cluster, [:f1, :f2, :t1, :t2, :c1, :c2, :p1, :p2, :o], "cluster")
    labels = lowercase.(labels(obj))
    channels = Int64[]

    cluster === :f1 && (cluster = ["fp1", "f1", "f3", "f5", "f7", "f9", "af3", "af7"])
    cluster === :f2 && (cluster = ["fp2", "f2", "f4", "f6", "f8", "f10", "af4", "af8"])
    cluster === :t1 && (cluster = ["c3", "c5", "t7", "t9", "fc3", "fc5", "ft7", "ft9"])
    cluster === :t2 && (cluster = ["c4", "c6", "t8", "t10", "fc4", "fc6", "ft8", "ft10"])
    cluster === :c1 && (cluster = ["cz", "c1", "c2", "fc1", "fc2", "fcz"])
    cluster === :c2 && (cluster = ["pz", "p1", "p2", "cp1", "cp2", "cpz"])
    cluster === :p1 && (cluster = ["p3", "p5", "p7", "p9", "cp3", "cp5", "tp7", "tp9"])
    cluster === :p2 && (cluster = ["p4", "p6", "p8", "p10", "cp4", "cp6", "tp8", "tp10"])
    cluster === :o && (cluster = ["o1", "o2", "poz", "po3", "po4", "po7", "po8", "po9", "po10"])

    for idx in cluster
        idx in labels && push!(channels, get_channel(obj, channel=idx))
    end

    return channels
end

"""
    band_frq(obj, band)

Return frequency limits for a `band`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `band::Symbol`: band range name:
    - `:list`
    - `:total`
    - `:delta`
    - `:theta`
    - `:alpha`
    - `:beta`
    - `:beta_high`
    - `:gamma`
    - `:gamma_1`
    - `:gamma_2`
    - `:gamma_lower`
    - `:gamma_higher`.
# Returns

- `band_frequency::Tuple{Real, Real}`
"""
function band_frq(obj::NeuroAnalyzer.NEURO; band::Symbol)

    bands = [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]
    _check_var(band, [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "band")
    if band === :list
        print("Available band names: ")
        [print(":$(bands[x]), ") for x in 2:(length(bands) - 1)]
        println(":$(bands[end])")
        return nothing
    end

    band === :total && (band_frq = (0.1, round(sr(obj) / 2, digits=1)))
    band === :delta && (band_frq = (0.1, 4.0))
    band === :theta && (band_frq = (4.0, 8.0))
    band === :alpha && (band_frq = (8.0, 13.0))
    band === :beta && (band_frq = (14.0, 30.0))
    band === :beta_high && (band_frq = (25.0, 30.0))
    band === :gamma && (band_frq = (30.0, 150.0))
    band === :gamma_1 && (band_frq = (31.0, 40.0))
    band === :gamma_2 && (band_frq = (41.0, 50.0))
    band === :gamma_lower && (band_frq = (30.0, 80.0))
    band === :gamma_higher && (band_frq = (80.0, 150.0))
    
    if band_frq[1] > sr(obj) / 2
        _info("Nyquist frequency based on sampling rate ($(sr(obj) / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(sr(obj) / 2 - 0.2), $(sr(obj) / 2 - 0.1))")
        band_frq = (sr(obj) / 2 - 0.2, sr(obj) / 2 - 0.1)
    end
    if band_frq[2] > sr(obj) / 2
        _info("Nyquist frequency based on sampling rate ($(sr(obj) / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(band_frq[1]), $(sr(obj) / 2 - 0.1))")
        band_frq = (band_frq[1], sr(obj) / 2 - 0.1)
    end

    return band_frq
end

"""
    band_frq(fs, band)

Return frequency limits of a `band`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `band::Symbol`: band range name:
    - `:list`
    - `:total`
    - `:delta`
    - `:theta`
    - `:alpha`
    - `:beta`
    - `:beta_high`
    - `:gamma`
    - `:gamma_1`
    - `:gamma_2`
    - `:gamma_lower`
    - `:gamma_higher`.

# Returns

- `band_frq::Tuple{Real, Real}`
"""
function band_frq(fs::Int64; band::Symbol)

    bands = [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]
    _check_var(band, [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "band")
    if band === :list
        print("Available band names: ")
        [print(":$(bands[x]), ") for x in 2:(length(bands) - 1)]
        println(":$(bands[end])")
        return nothing
    end

    band === :total && (band_frq = (0.1, round(fs / 2, digits=1)))
    band === :delta && (band_frq = (0.1, 4.0))
    band === :theta && (band_frq = (4.0, 8.0))
    band === :alpha && (band_frq = (8.0, 13.0))
    band === :beta && (band_frq = (14.0, 30.0))
    band === :beta_high && (band_frq = (25.0, 30.0))
    band === :gamma && (band_frq = (30.0, 150.0))
    band === :gamma_1 && (band_frq = (31.0, 40.0))
    band === :gamma_2 && (band_frq = (41.0, 50.0))
    band === :gamma_lower && (band_frq = (30.0, 80.0))
    band === :gamma_higher && (band_frq = (80.0, 150.0))
    
    if band_frq[1] > fs / 2
        _info("Nyquist frequency based on EEG sampling rate ($(fs / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(fs / 2 - 0.2), $(fs / 2 - 0.1))")
        band_frq = (fs / 2 - 0.2, fs / 2 - 0.1)
    end
    if band_frq[2] > fs / 2
        _info("Nyquist frequency based on EEG sampling rate ($(fs / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(band_frq[1]), $(fs / 2 - 0.1))")
        band_frq = (band_frq[1], fs / 2 - 0.1)
    end

    return band_frq
end
