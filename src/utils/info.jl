export sr
export channel_n
export epoch_n
export signal_len
export epoch_len
export history
export labels
export info
export channel_cluster
export band_frq
export describe
export size

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
    channel_n(obj; type)

Return number of channels of `type`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Vector{Symbol}=:all`: channel type (stored in the global channel_types constant variable)

# Returns

- `channel_n::Int64`
"""
function channel_n(obj::NeuroAnalyzer.NEURO; type::Symbol=:all)

    _check_var(type, channel_types, "type")

    if type === :all
        if ndims(obj.data) == 1
            ch_n = 1
        else
            ch_n = size(obj.data, 1)
        end
    else
        length(obj.header.recording[:channel_type]) == 0 && throw(ArgumentError("OBJ has no defined channel types."))
        ch_n = 0
        for idx in 1:channel_n(obj)
            obj.header.recording[:channel_type][idx] == string(type) && (ch_n += 1)
        end
    end

    return ch_n

end

"""
    epoch_n(obj)

Return number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `epoch_n::Int64`
"""
function epoch_n(obj::NeuroAnalyzer.NEURO)

    if ndims(obj.data) < 3
        throw(ArgumentError("Record data is either a vector or a matrix."))
    else
        ep_n = size(obj.data, 3)
    end

    return ep_n

end

"""
    signal_len(obj)

Return signal length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `s_len::Int64`
"""
function signal_len(obj::NeuroAnalyzer.NEURO)

    if ndims(obj.data) == 1
        s_len = length(obj.data)
    elseif ndims(obj.data) == 2
        s_len = size(obj.data, 2)
    else
        s_len = size(obj.data, 2) * size(obj.data, 3)
    end

    return s_len

end

"""
    epoch_len(obj)

Return epoch length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `ep_len::Int64`
"""
function epoch_len(obj::NeuroAnalyzer.NEURO)

    if ndims(obj.data) < 3
        throw(ArgumentError("Record data is either a vector or a matrix."))
    else
        ep_len = size(obj.data, 2)
    end

    return ep_len

end

"""
    history(obj)

Show processing history.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `history::Vector{String}`
"""
function history(obj::NeuroAnalyzer.NEURO)

    return obj.history

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

    if length(obj.header.recording[:labels]) == 0
        throw(ArgumentError("RECORD has no labels."))
    else
        return obj.header.recording[:labels]
    end

end

"""
    info(obj)

Show info.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function info(obj::NeuroAnalyzer.NEURO)

    println("              Data type: $(uppercase(obj.header.recording[:data_type]))")
    println("            File format: $(obj.header.recording[:file_type])")
    println("            Source file: $(obj.header.recording[:file_name])")
    println("         File size [MB]: $(obj.header.recording[:file_size_mb])")
    println("       Memory size [MB]: $(round(Base.summarysize(obj) / 1024^2, digits=2))")
    if length(obj.header.subject[:id]) > 0
        println("                Subject: $(obj.header.subject[:id] * ": " * obj.header.subject[:first_name] * " " * obj.header.subject[:last_name])")
    else
        println("                Subject: $(obj.header.subject[:first_name] * " " * obj.header.subject[:last_name])")
    end
    println("              Recording: $(obj.header.recording[:recording])")
    println("        Recording notes: $(obj.header.recording[:recording_notes])")
    println("         Recording date: $(obj.header.recording[:recording_date])")
    println("         Recording time: $(obj.header.recording[:recording_time])")
    println("     Sampling rate (Hz): $(sr(obj))")
    println("Signal length [samples]: $(signal_len(obj))")
    println("Signal length [seconds]: $(round(signal_len(obj) / sr(obj), digits=2))")
    println("     Number of channels: $(channel_n(obj))")
    println("       Number of epochs: $(epoch_n(obj))")
    println(" Epoch length [samples]: $(epoch_len(obj))")
    println(" Epoch length [seconds]: $(round(epoch_len(obj) / sr(obj), digits=2))")
    if obj.header.recording[:data_type] == "eeg"
        if obj.header.recording[:reference] == ""
            println("         Reference type: unknown")
        else
            println("         Reference type: $(obj.header.recording[:reference])")
        end
    end
    if obj.header.recording[:data_type] == "nirs"
        println("        Wavelength [nm]: $(obj.header.recording[:wavelengths])")
    end
    if length(labels(obj)) == 0
        println("                 Labels: no")
    else
        println("                 Labels: yes")
    end
    if _has_markers(obj) == false
        println("                Markers: no")
    else
        println("                Markers: yes")
    end
    if _has_locs(obj) == false
        println("      Channel locations: no")
    else
        println("      Channel locations: yes")
    end
    if length(keys(obj.components)) > 0
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
    if obj.header.recording[:data_type] != "nirs"
        println(rpad(" ch", 8) * 
                rpad("label", 16) * 
                rpad("type", 12) * 
                rpad("unit", 8))
        for idx in eachindex(obj.header.recording[:labels])
            println(rpad(" $idx", 8) * 
                    rpad("$(obj.header.recording[:labels][idx])", 16) * 
                    rpad("$(uppercase(obj.header.recording[:channel_type][idx]))", 12) * 
                    rpad("$(obj.header.recording[:units][idx])", 8))
        end
    else
        if obj.header.recording[:channel_type][idx] !== "nirs_aux"
            println(rpad(" ch", 8) * 
                    rpad("label", 16) * 
                    rpad("type", 12) * 
                    rpad("unit", 8) * 
                    rpad("wavelength", 8))
            for idx in eachindex(obj.header.recording[:labels])
                println(rpad(" $idx", 8) * 
                        rpad("$(obj.header.recording[:labels][idx])", 16) * 
                        rpad("$(uppercase(obj.header.recording[:channel_type][idx]))", 12) * 
                        rpad("$(obj.header.recording[:units][idx])", 8) * 
                        rpad("$(obj.header.recording[:wavelength_index][idx])", 8))
            end
        else
            println(rpad(" ch", 8) * 
                    rpad("label", 16) * 
                    rpad("type", 12))
            for idx in eachindex(obj.header.recording[:labels])
                println(rpad(" $idx", 8) * 
                        rpad("$(obj.header.recording[:labels][idx])", 16) * 
                        rpad("$(uppercase(obj.header.recording[:channel_type][idx]))", 12))
            end
        end
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

- `ch::Vector{Int64}`: list of channel numbers belonging to a given cluster of channels
"""
function channel_cluster(obj::NeuroAnalyzer.NEURO; cluster::Symbol)

    length(labels(obj)) == 0 && throw(ArgumentError("OBJ does not contain channel labels."))

    _check_var(cluster, [:f1, :f2, :t1, :t2, :c1, :c2, :p1, :p2, :o], "cluster")
    clabels = lowercase.(labels(obj))
    ch = Int64[]

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
        idx in clabels && push!(ch, get_channel(obj, ch=idx))
    end

    return ch
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
    _check_var(band, bands, "band")
    if band === :list
        print("Available band names: ")
        [print(":$(bands[x]), ") for x in 2:(length(bands) - 1)]
        println(":$(bands[end])")
        return nothing
    end

    band === :total && (bf = (0.1, round(sr(obj) / 2, digits=1)))
    band === :delta && (bf = (0.1, 4.0))
    band === :theta && (bf = (4.0, 8.0))
    band === :alpha && (bf = (8.0, 13.0))
    band === :beta && (bf = (14.0, 30.0))
    band === :beta_high && (bf = (25.0, 30.0))
    band === :gamma && (bf = (30.0, 150.0))
    band === :gamma_1 && (bf = (31.0, 40.0))
    band === :gamma_2 && (bf = (41.0, 50.0))
    band === :gamma_lower && (bf = (30.0, 80.0))
    band === :gamma_higher && (bf = (80.0, 150.0))
    
    if bf[1] > sr(obj) / 2
        _info("Nyquist frequency based on sampling rate ($(sr(obj) / 2)) is lower than $band range: $bf, band frequency truncated to: ($(sr(obj) / 2 - 0.2), $(sr(obj) / 2 - 0.1))")
        bf = (sr(obj) / 2 - 0.2, sr(obj) / 2 - 0.1)
    end
    if bf[2] > sr(obj) / 2
        _info("Nyquist frequency based on sampling rate ($(sr(obj) / 2)) is lower than $band range: $bf, band frequency truncated to: ($(bf[1]), $(sr(obj) / 2 - 0.1))")
        bf = (bf[1], sr(obj) / 2 - 0.1)
    end

    return bf
end

"""
    band_frq(fs, band)

Return frequency limits of a `band`.

# Arguments

- `fs::Int64`: sampling rate
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
    _check_var(band, bands, "band")
    if band === :list
        print("Available band names: ")
        [print(":$(bands[x]), ") for x in 2:(length(bands) - 1)]
        println(":$(bands[end])")
        return nothing
    end

    band === :total && (bf = (0.1, round(fs / 2, digits=1)))
    band === :delta && (bf = (0.1, 4.0))
    band === :theta && (bf = (4.0, 8.0))
    band === :alpha && (bf = (8.0, 13.0))
    band === :beta && (bf = (14.0, 30.0))
    band === :beta_high && (bf = (25.0, 30.0))
    band === :gamma && (bf = (30.0, 150.0))
    band === :gamma_1 && (bf = (31.0, 40.0))
    band === :gamma_2 && (bf = (41.0, 50.0))
    band === :gamma_lower && (bf = (30.0, 80.0))
    band === :gamma_higher && (bf = (80.0, 150.0))
    
    if bf[1] > fs / 2
        _info("Nyquist frequency based on sampling rate ($(fs / 2)) is lower than $band range: $bf, band frequency truncated to: ($(fs / 2 - 0.2), $(fs / 2 - 0.1))")
        bf = (fs / 2 - 0.2, fs / 2 - 0.1)
    end
    if bf[2] > fs / 2
        _info("Nyquist frequency based on sampling rate ($(fs / 2)) is lower than $band range: $bf, band frequency truncated to: ($(bf[1]), $(fs / 2 - 0.1))")
        bf = (bf[1], fs / 2 - 0.1)
    end

    return bf
end

"""
    describe(obj)

Return basic descriptive statistics of `obj.data`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function describe(obj::NeuroAnalyzer.NEURO)
    println("< " * uppercase(obj.header.recording[:data_type]) * ", $(channel_n(obj)) × $(epoch_len(obj)) × $(epoch_n(obj)) ($(signal_len(obj) / sr(obj)) s) >")
    println(rpad("ch", 4) * 
            rpad("label", 16) * 
            rpad("type", 12) * 
            rpad("unit", 8) * 
            rpad("range", 10) * 
            rpad("mean", 10) * 
            rpad("sd", 10) * 
            rpad("min", 10) * 
            rpad("Q1", 10) * 
            rpad("median", 10) * 
            rpad("Q3", 10) * 
            rpad("max", 10))
    for idx in 1:channel_n(obj)
        println(rpad(string(idx), 4) * 
                rpad(labels(obj)[idx], 16) * 
                rpad(uppercase(obj.header.recording[:channel_type][idx]), 12) * 
                rpad(obj.header.recording[:units][idx], 8) * 
                rpad(round(rng(obj.data[idx, :, :]), digits=2), 10) * 
                rpad(round(mean(obj.data[idx, :, :]), digits=2), 10) * 
                rpad(round(std(obj.data[idx, :, :]), digits=2), 10) * 
                rpad(round(minimum(obj.data[idx, :, :]), digits=2), 10) * 
                rpad(round(quantile(obj.data[idx, :, :][:], 0.5), digits=2), 10) * 
                rpad(round(median(obj.data[idx, :, :]), digits=2), 10) * 
                rpad(round(quantile(obj.data[idx, :, :][:], 0.95), digits=2), 10) * 
                rpad(round(maximum(obj.data[idx, :, :]), digits=2), 10))
    end
end

"""
    size(obj)

Return size of the `obj` data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `size::Tuple{Int64, Int64, Int64}`
"""
function Base.size(obj::NeuroAnalyzer.NEURO)
    
    return size(obj.data)

end