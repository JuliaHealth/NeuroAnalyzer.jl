export sr
export nchannels
export nepochs
export signal_len
export epoch_len
export signal_duration
export epoch_duration
export history
export labels
export optode_labels
export source_labels
export detector_labels
export chtypes
export header
export info
export channel_info
export channel_pick
export channel_cluster
export band_frq
export describe
export size
export datatype
export channel_order

"""
    sr(obj)

Return sampling rate.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `fs::Int64`
"""
function sr(obj::NeuroAnalyzer.NEURO)::Int64

    fs = obj.header.recording[:sampling_rate]

    return fs

end

"""
    nchannels(obj; <keyword arguments>)

Return number `type` channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::String="all"`: channel type (stored in the global `channel_types` constant variable)

# Returns

- `ch_n::Int64`
"""
function nchannels(obj::NeuroAnalyzer.NEURO; type::String="all")::Int64

    _check_var(type, channel_types, "type")
    @assert length(obj.header.recording[:channel_type]) != 0 "OBJ has no defined channel types."

    if type == "all"
        if ndims(obj.data) == 1
            ch_n = 1
        else
            ch_n = size(obj.data, 1)
        end
    else
        ch_n = 0
        for idx in 1:nchannels(obj)
            obj.header.recording[:channel_type][idx] == type && (ch_n += 1)
        end
    end

    return ch_n

end

"""
    nepochs(obj)

Return number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `ep_n::Int64`
"""
function nepochs(obj::NeuroAnalyzer.NEURO)::Int64

    @assert ndims(obj.data) == 3 "Record data is either a vector or a matrix."

    ep_n = size(obj.data, 3)

    return ep_n

end

"""
    signal_len(obj)

Return signal length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `sl::Int64`
"""
function signal_len(obj::NeuroAnalyzer.NEURO)::Int64

    if ndims(obj.data) == 1
        sl = length(obj.data)
    elseif ndims(obj.data) == 2
        sl = size(obj.data, 2)
    else
        sl = size(obj.data, 2) * size(obj.data, 3)
    end

    return sl

end

"""
    epoch_len(obj)

Return epoch length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `len::Int64`
"""
function epoch_len(obj::NeuroAnalyzer.NEURO)::Int64

    @assert ndims(obj.data) == 3 "Record data is either a vector or a matrix."

    len = size(obj.data, 2)

    return len

end

"""
    signal_duration(obj)

Return signal duration.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `sd::Float64`
"""
function signal_duration(obj::NeuroAnalyzer.NEURO)::Float64

    sd = obj.time_pts[end]

    return sd

end

"""
    epoch_duration(obj)

Return epoch length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `ed::Float64`
"""
function epoch_duration(obj::NeuroAnalyzer.NEURO)::Float64

    ed = obj.epoch_time[end]

    return ed

end

"""
    history(obj)

Show processing history.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `h::Vector{String}`
"""
function history(obj::NeuroAnalyzer.NEURO)::Vector{String}

    h = obj.history

    return h

end

"""
    labels(obj)

Return channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `l::Vector{String}`
"""
function labels(obj::NeuroAnalyzer.NEURO)::Vector{String}

    @assert length(obj.header.recording[:label]) > 0 "OBJ has no channel labels."

    l = obj.header.recording[:label]

    return l

end

"""
    optode_labels(obj)

Return optode labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `l::Vector{String}`
"""
function optode_labels(obj::NeuroAnalyzer.NEURO)::Vector{String}

    _check_datatype(obj, "nirs")
    @assert length(obj.header.recording[:optode_labels]) > 0 "OBJ has no optode labels."

    l = obj.header.recording[:optode_labels]

    return l

end

"""
    source_labels(obj)

Return NIRS source labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `l::Vector{String}`
"""
function source_labels(obj::NeuroAnalyzer.NEURO)::Vector{String}

    _check_datatype(obj, "nirs")
    @assert length(obj.header.recording[:src_labels]) > 0 "OBJ has no source labels."

    l = obj.header.recording[:src_labels]

    return l

end

"""
    detector_labels(obj)

Return NIRS detector labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `l::Vector{String}`
"""
function detector_labels(obj::NeuroAnalyzer.NEURO)::Vector{String}

    _check_datatype(obj, "nirs")
    @assert length(obj.header.recording[:det_labels]) > 0 "OBJ has no detector labels."

    l = obj.header.recording[:det_labels]

    return l

end

"""
    chtypes(obj)

Return channel types.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `cht::Vector{String}`
"""
function chtypes(obj::NeuroAnalyzer.NEURO)::Vector{String}

    @assert length(obj.header.recording[:channel_type]) > 0 "OBJ has no channel chtypes."

    cht = obj.header.recording[:channel_type]

    return cht

end

"""
    header(obj; <keyword arguments>)

Show object header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Nothing
"""
function header(obj::NeuroAnalyzer.NEURO)::Nothing

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
    println("Signal length [seconds]: $(round(signal_len(obj) / sr(obj), digits=3))")
    println("     Number of channels: $(nchannels(obj))")
    if !(datatype(obj) in ["mep", "sensors", "eda"])
        println("              Epochs ID: $(obj.header.recording[:epoch_id])")
        println("       Number of epochs: $(nepochs(obj))")
        println(" Epoch length [samples]: $(epoch_len(obj))")
        println(" Epoch length [seconds]: $(round(epoch_len(obj) / sr(obj), digits=3))")
    end
    if datatype(obj) == "eeg"
        if obj.header.recording[:reference] == ""
            println("         Reference type: unknown")
        else
            println("         Reference type: $(obj.header.recording[:reference])")
        end
    end
    if datatype(obj) == "meg"
        ssp_labels = obj.header.recording[:ssp_labels]
        if length(ssp_labels) > 0
            if length(ssp_labels) == 1
                println("       SSP projection: $(obj.header.recording[:ssp_labels][1])")
            else
                print("        SSP projections: ")
                for idx in 1:(length(ssp_labels) - 1)
                    print("$(obj.header.recording[:ssp_labels][idx]), ")
                end
                println("$(obj.header.recording[:ssp_labels][end])")
            end
        end
    end
    if datatype(obj) in ["eeg", "meg", "ecog", "seeg"]
        println("         Line frequency: $(obj.header.recording[:line_frequency]) Hz")
    end
    if datatype(obj) == "nirs"
        println("        Wavelength [nm]: $(obj.header.recording[:wavelengths])")
    end
    if !(datatype(obj) in ["mep", "sensors", "eda"])
        if _has_markers(obj)
            println("                Markers: yes")
        else
            println("                Markers: no")
        end
        if DataFrames.nrow(obj.locs) > 0
            println("      Channel locations: yes")
        else
            println("      Channel locations: no")
        end
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
    if datatype(obj) in ["eeg", "ecog", "seeg", "ieeg", "erp"]
        nch = 0
        for idx in ["eeg", "ecog", "seeg", "ieeg", "erp"]
            nch += count(x -> isequal(x, idx), obj.header.recording[:channel_type])
        end
        println(" Number of EEG channels: $nch")
    elseif datatype(obj) in ["meg", "erf"]
        nch = count(x -> isequal(x, "mag"), obj.header.recording[:channel_type])
        println(" Number of MAG channels: $nch")
        nch = count(x -> isequal(x, "grad"), obj.header.recording[:channel_type])
        println("Number of GRAD channels: $nch")
        nch = count(x -> isequal(x, "eeg"), obj.header.recording[:channel_type])
        println(" Number of EEG channels: $nch")
    elseif datatype(obj) == "nirs"
        nch = 0
        for idx in ["nirs", "nirs_int", "nirs_od", "nirs_dmean", "nirs_dvar", "nirs_dskew", "nirs_mua", "nirs_musp", "nirs_hbo", "nirs_hbr", "nirs_hbt", "nirs_h2o", "nirs_lipid", "nirs_bfi", "nirs_hrf_dod", "nirs_hrf_dmean", "nirs_hrf_dvar", "nirs_hrf_dskew", "nirs_hrf_hbo", "nirs_hrf_hbr", "nirs_hrf_hbt", "nirs_hrf_bfi", "nirs_aux"]
            nch += count(x -> isequal(x, idx), obj.header.recording[:channel_type])
        end
        println("Number of NIRS channels: $nch")
    end

    return nothing

end

"""
    info(obj; <keyword arguments>)

Show info.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `df::Bool=false`: if true, return object data as a DataFrame containing time points and channels

# Returns

- `Union{Nothing, DataFrame}`
"""
function info(obj::NeuroAnalyzer.NEURO; df::Bool=false)::Union{Nothing, DataFrame}

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
    println("Signal length [seconds]: $(round(signal_len(obj) / sr(obj), digits=3))")
    println("     Number of channels: $(nchannels(obj))")
    if !(datatype(obj) in ["mep", "sensors", "eda"])
        println("              Epochs ID: $(obj.header.recording[:epoch_id])")
        println("       Number of epochs: $(nepochs(obj))")
        println(" Epoch length [samples]: $(epoch_len(obj))")
        println(" Epoch length [seconds]: $(round(epoch_len(obj) / sr(obj), digits=3))")
    end
    if datatype(obj) == "eeg"
        if obj.header.recording[:reference] == ""
            println("         Reference type: unknown")
        else
            println("         Reference type: $(obj.header.recording[:reference])")
        end
    end
    if datatype(obj) == "meg"
        ssp_labels = obj.header.recording[:ssp_labels]
        if length(ssp_labels) > 0
            if length(ssp_labels) == 1
                println("       SSP projection: $(obj.header.recording[:ssp_labels][1])")
            else
                print("        SSP projections: ")
                for idx in 1:(length(ssp_labels) - 1)
                    print("$(obj.header.recording[:ssp_labels][idx]), ")
                end
                println("$(obj.header.recording[:ssp_labels][end])")
            end
        end
    end
    if datatype(obj) in ["eeg", "meg", "ecog", "seeg"]
        println("         Line frequency: $(obj.header.recording[:line_frequency]) Hz")
    end
    if datatype(obj) == "nirs"
        println("        Wavelength [nm]: $(obj.header.recording[:wavelengths])")
    end
    if !(datatype(obj) in ["mep", "sensors", "eda"])
        if _has_markers(obj)
            println("                Markers: yes")
        else
            println("                Markers: no")
        end
        if DataFrames.nrow(obj.locs) > 0
            println("      Channel locations: yes")
        else
            println("      Channel locations: no")
        end
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
    if datatype(obj) in ["eeg", "ecog", "seeg", "ieeg", "erp"]
        nch = 0
        for idx in ["eeg", "ecog", "seeg", "ieeg", "erp"]
            nch += count(x -> isequal(x, idx), obj.header.recording[:channel_type])
        end
        println(" Number of EEG channels: $nch")
    elseif datatype(obj) in ["meg", "erf"]
        nch = count(x -> isequal(x, "mag"), obj.header.recording[:channel_type])
        println(" Number of MAG channels: $nch")
        nch = count(x -> isequal(x, "grad"), obj.header.recording[:channel_type])
        println("Number of GRAD channels: $nch")
        nch = count(x -> isequal(x, "eeg"), obj.header.recording[:channel_type])
        println(" Number of EEG channels: $nch")
    elseif datatype(obj) == "nirs"
        nch = 0
        for idx in ["nirs", "nirs_int", "nirs_od", "nirs_dmean", "nirs_dvar", "nirs_dskew", "nirs_mua", "nirs_musp", "nirs_hbo", "nirs_hbr", "nirs_hbt", "nirs_h2o", "nirs_lipid", "nirs_bfi", "nirs_hrf_dod", "nirs_hrf_dmean", "nirs_hrf_dvar", "nirs_hrf_dskew", "nirs_hrf_hbo", "nirs_hrf_hbr", "nirs_hrf_hbt", "nirs_hrf_bfi", "nirs_aux"]
            nch += count(x -> isequal(x, idx), obj.header.recording[:channel_type])
        end
        println("Number of NIRS channels: $nch")
    end
    println()
    println("Channels:")
    if obj.header.recording[:data_type] != "nirs"
        println(rpad(" ch", 8) *
                rpad("label", 16) *
                rpad("type", 12) *
                rpad("unit", 8) *
                rpad("bad", 8))
        println(" " * repeat("-", 6) * " " *
                repeat("-", 15) * " " *
                repeat("-", 11) * " " *
                repeat("-", 7) * " " *
                repeat("-", 7))
        for idx in eachindex(obj.header.recording[:label])
            println(rpad(" $idx", 8) *
                    rpad("$(obj.header.recording[:label][idx])", 16) *
                    rpad("$(uppercase(obj.header.recording[:channel_type][idx]))", 12) *
                    rpad("$(obj.header.recording[:unit][idx])", 8) *
                    rpad("$(obj.header.recording[:bad_channel][idx])", 8))
        end
    else
        println(rpad(" ch", 8) *
                rpad("label", 16) *
                rpad("type", 12) *
                rpad("unit", 8) *
                rpad("wavelength", 8))
        println(" " * repeat("-", 6) * " " *
                repeat("-", 15) * " " *
                repeat("-", 11) * " " *
                repeat("-", 7) * " " *
                repeat("-", 12))
        for idx in eachindex(obj.header.recording[:label])
            if !(obj.header.recording[:channel_type][idx] in ["nirs_aux", "nirs_hbo", "nirs_hbr", "nirs_hbt"])
                println(rpad(" $idx", 8) *
                        rpad("$(obj.header.recording[:label][idx])", 16) *
                        rpad("$(uppercase(obj.header.recording[:channel_type][idx]))", 12) *
                        rpad("$(obj.header.recording[:unit][idx])", 8) *
                        rpad("$(obj.header.recording[:wavelengths][obj.header.recording[:wavelength_index][idx]])", 8))
            else
                println(rpad(" $idx", 8) *
                        rpad("$(obj.header.recording[:label][idx])", 16) *
                        rpad("$(uppercase(obj.header.recording[:channel_type][idx]))", 12))
            end
        end
    end

    if df
        df = DataFrame(hcat(obj.time_pts, reshape(obj.data, nchannels(obj), :, 1)[:, :]'), :auto)
        DataFrames.rename!(df, vcat(:time, Symbol.(labels(obj))))
        return df
    else
        return nothing
    end

end

"""
    channel_info(obj; <keyword arguments>)

Show channel info.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`

# Returns

Nothing
"""
function channel_info(obj::NeuroAnalyzer.NEURO; ch::Int64)::Nothing

    _check_channels(obj, labels(obj)[ch])

    if obj.header.recording[:data_type] != "nirs"
        println(" ch: $(rpad(string(ch), 4)) label: $(rpad(obj.header.recording[:label][ch], 8)) type: $(rpad(uppercase(obj.header.recording[:channel_type][ch]), 8)) unit: $(rpad(obj.header.recording[:unit][ch], 8)) bad: $(obj.header.recording[:bad_channel][ch])")
    else
        if obj.header.recording[:channel_type][ch] !== "nirs_aux"
            println(" ch: $(rpad(string(ch), 4)) label: $(rpad(obj.header.recording[:label][ch], 8)) type: $(rpad(uppercase(obj.header.recording[:channel_type][ch]), 8)) unit: $(rpad(obj.header.recording[:unit][ch], 8)) wavelength: $(rpad((obj.header.recording[:wavelength_index][ch]), 8))")
        else
            println(" ch: $(rpad(string(ch), 4)) label: $(rpad(obj.header.recording[:label][ch], 8)) type: $(rpad(uppercase(obj.header.recording[:channel_type][ch]), 8))")
        end
    end

    return nothing

end

"""
    channel_pick(obj; <keyword arguments>)

Return set of channel indices corresponding to a set of electrodes ("pick", e.g. left or frontal electrodes).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `p::Vector{Symbol}`: pick of electrodes; picks may be combined, e.g. `[:left, :frontal]`
    - `:central` (or `:c`)
    - `:left` (or `:l`)
    - `:right` (or `:r`)
    - `:frontal` (or `:f`)
    - `:temporal` (or `:t`)
    - `:parietal` (or `:p`)
    - `:occipital` (or `:o`)

# Returns

- `channels::Vector{Int64}`: channel numbers
"""
function channel_pick(obj::NeuroAnalyzer.NEURO; p::Union{Symbol, Vector{Symbol}})::Vector{Int64}

    _check_datatype(obj, "eeg")

    @assert length(labels(obj)) != 0 "OBJ does not contain channel labels."

    if p isa Vector{Symbol}
        for idx in p
            _check_var(idx, [:central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "p")
        end

        # convert picks to channel labels
        c = Vector{Char}()
        for idx in p
            (idx === :central || idx === :c) && push!(c, 'z')
            (idx === :frontal || idx === :f) && push!(c, 'F')
            (idx === :temporal || idx === :t) && push!(c, 'T')
            (idx === :parietal || idx === :p) && push!(c, 'P')
            (idx === :occipital || idx === :o) && push!(c, 'O')
        end

        # check which channels are in the picks list
        clabels = get_channel(obj, type="eeg")
        channels = Vector{Int64}()
        for idx1 in eachindex(clabels)
            for idx2 in eachindex(c)
                in(c[idx2], clabels[idx1]) && push!(channels, idx1)
            end
        end

        # check for both :l and :r
        for idx1 in eachindex(p)
            if (p[idx1] === :left || p[idx1] === :l)
                for idx2 in eachindex(p)
                    if (p[idx2] === :right || p[idx2] === :r)
                        return channels
                    end
                end
            end
            if (p[idx1] === :right || p[idx1] === :r)
                for idx2 in eachindex(p)
                    if (p[idx2] === :left || p[idx2] === :l)
                        return channels
                    end
                end
            end
        end

        clabels = get_channel(obj, type="eeg")
        clabels = clabels[channels]
        pat = nothing
        for idx in p
            # for :right remove lefts
            (idx === :right || idx === :r) && (pat = r"[z13579]$")
            # for :left remove rights
            (idx === :left || idx === :l) && (pat = r"[z02468]$")
        end
        if typeof(pat) == Regex
            for idx in length(channels):-1:1
                !isnothing(match(pat, clabels[idx])) && deleteat!(channels, idx)
            end
        end
        return channels
    else
        _check_var(p, [:central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "p")

        c = Vector{Char}()
        (p === :central || p === :c) && (c = ['z'])
        (p === :left || p === :l) && (c = ['1', '3', '5', '7', '9'])
        (p === :right || p === :r) && (c = ['2', '4', '6', '8'])
        (p === :frontal || p === :f) && (c = ['F'])
        (p === :temporal || p === :t) && (c = ['T'])
        (p === :parietal || p === :p) && (c = ['P'])
        (p === :occipital || p === :o) && (c = ['O'])

        clabels = get_channel(obj, type="eeg")
        channels = Vector{Int64}()
        for idx1 in eachindex(c)
            for idx2 in eachindex(clabels)
                in(c[idx1], clabels[idx2]) && push!(channels, idx2)
            end
        end

        return channels
    end
end

"""
    channels_cluster(obj; <keyword arguments>)

Return channels belonging to a cluster of channels.

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

- `ch::Vector{Int64}`: channel numbers
"""
function channel_cluster(obj::NeuroAnalyzer.NEURO; cluster::Symbol)::Vector{Int64}

    @assert length(labels(obj)) != 0 "OBJ does not contain channel labels."

    _check_var(cluster, [:f1, :f2, :t1, :t2, :c1, :c2, :p1, :p2, :o], "cluster")
    clabels = labels(obj)
    ch = String[]

    cluster === :f1 && (cluster = ["Fp1", "F1", "F3", "F5", "F7", "F9", "AF3", "AF7"])
    cluster === :f2 && (cluster = ["Fp2", "F2", "F4", "F6", "F8", "F10", "AF4", "AF8"])
    cluster === :t1 && (cluster = ["C3", "C5", "T7", "T9", "FC3", "FC5", "FT7", "FT9"])
    cluster === :t2 && (cluster = ["C4", "C6", "T8", "T10", "FC4", "FC6", "FT8", "FT10"])
    cluster === :c1 && (cluster = ["Cz", "C1", "C2", "FC1", "FC2", "FCz"])
    cluster === :c2 && (cluster = ["Pz", "P1", "P2", "CP1", "CP2", "CPz"])
    cluster === :P1 && (cluster = ["P3", "P5", "P7", "P9", "CP3", "CP5", "TP7", "TP9"])
    cluster === :P2 && (cluster = ["P4", "P6", "P8", "P10", "CP4", "CP6", "TP8", "TP10"])
    cluster === :o && (cluster = ["O1", "O2", "POz", "PO3", "PO4", "PO7", "PO8", "PO9", "PO10"])

    for idx in cluster
        idx in clabels && push!(ch, idx)
    end

    return get_channel(obj, ch=ch)

end

"""
    band_frq(obj, band)

Return band frequency limits.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `band::Symbol`: band range name:
    - `:list`
    - `:total`
    - `:delta`: 0.1 - 4.0 Hz
    - `:theta`: 4.0 - 8.0 Hz
    - `:alpha`: 8.0 - 13.0 Hz
    - `:alpha_lower`: 8.0 - 10.5 Hz
    - `:alpha_higher`: 10.5 - 13.0 Hz
    - `:beta`: 14.0 - 30.0 Hz
    - `:beta_lower`: 14.0 - 25.0 Hz
    - `:beta_higher`: 25.0 - 30.0 Hz
    - `:gamma`: 30.0 - 150.0 Hz
    - `:gamma_1`: 30.0 - 40.0 Hz
    - `:gamma_2`: 40.0 - 50.0 Hz
    - `:gamma_lower`: 30.0 - 80.0 Hz
    - `:gamma_higher`: 80.0 - 150.0 Hz
# Returns

- `bf::Tuple{Float64, Float64}`
"""
function band_frq(obj::NeuroAnalyzer.NEURO; band::Symbol)::Tuple{Float64, Float64}

    bands = [:list, :total, :delta, :theta, :alpha, :alpha_lower, :alpha_higher, :beta, :beta_lower, :beta_higher, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]
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
    band === :alpha_lower && (bf = (8.0, 10.5))
    band === :alpha_higher && (bf = (10.5, 13.0))
    band === :beta && (bf = (14.0, 30.0))
    band === :beta_lower && (bf = (14.0, 25.0))
    band === :beta_higher && (bf = (25.0, 30.0))
    band === :gamma && (bf = (30.0, 150.0))
    band === :gamma_1 && (bf = (30.0, 40.0))
    band === :gamma_2 && (bf = (40.0, 50.0))
    band === :gamma_lower && (bf = (30.0, 80.0))
    band === :gamma_higher && (bf = (80.0, 150.0))

    if bf[1] > sr(obj) / 2
        _warn("Nyquist frequency based on sampling rate ($(sr(obj) / 2)) is lower than $band range: $bf, band frequency truncated to: ($(sr(obj) / 2 - 0.2), $(sr(obj) / 2 - 0.1)).")
        bf = (sr(obj) / 2 - 0.2, sr(obj) / 2 - 0.1)
    end
    if bf[2] > sr(obj) / 2
        _warn("Nyquist frequency based on sampling rate ($(sr(obj) / 2)) is lower than $band range: $bf, band frequency truncated to: ($(bf[1]), $(sr(obj) / 2 - 0.1)).")
        bf = (bf[1], sr(obj) / 2 - 0.1)
    end

    return bf

end

"""
    band_frq(fs, band)

Return band frequency limits.

# Arguments

- `fs::Int64`: sampling rate
- `band::Symbol`: band range name:
    - `:list`
    - `:total`
    - `:delta`: 0.1 - 4.0 Hz
    - `:theta`: 4.0 - 8.0 Hz
    - `:alpha`: 8.0 - 13.0 Hz
    - `:alpha_lower`: 8.0 - 10.5 Hz
    - `:alpha_higher`: 10.5 - 13.0 Hz
    - `:beta`: 14.0 - 30.0 Hz
    - `:beta_lower`: 14.0 - 25.0 Hz
    - `:beta_higher`: 25.0 - 30.0 Hz
    - `:gamma`: 30.0 - 150.0 Hz
    - `:gamma_1`: 30.0 - 40.0 Hz
    - `:gamma_2`: 40.0 - 50.0 Hz
    - `:gamma_lower`: 30.0 - 80.0 Hz
    - `:gamma_higher`: 80.0 - 150.0 Hz

# Returns

- `bf::Tuple{Float64, Float64}`
"""
function band_frq(fs::Int64; band::Symbol)::Tuple{Float64, Float64}

    bands = [:list, :total, :delta, :theta, :alpha, :alpha_lower, :alpha_higher, :beta, :beta_lower, :beta_higher, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]
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
    band === :alpha_lower && (bf = (8.0, 10.5))
    band === :alpha_higher && (bf = (10.5, 13.0))
    band === :beta && (bf = (14.0, 30.0))
    band === :beta_lower && (bf = (14.0, 25.0))
    band === :beta_higher && (bf = (25.0, 30.0))
    band === :gamma && (bf = (30.0, 150.0))
    band === :gamma_1 && (bf = (30.0, 40.0))
    band === :gamma_2 && (bf = (40.0, 50.0))
    band === :gamma_lower && (bf = (30.0, 80.0))
    band === :gamma_higher && (bf = (80.0, 150.0))

    if bf[1] > fs / 2
        _warn("Nyquist frequency based on sampling rate ($(fs / 2)) is lower than $band range: $bf, band frequency truncated to: ($(fs / 2 - 0.2), $(fs / 2 - 0.1)).")
        bf = (fs / 2 - 0.2, fs / 2 - 0.1)
    end
    if bf[2] > fs / 2
        _warn("Nyquist frequency based on sampling rate ($(fs / 2)) is lower than $band range: $bf, band frequency truncated to: ($(bf[1]), $(fs / 2 - 0.1)).")
        bf = (bf[1], fs / 2 - 0.1)
    end

    return bf

end

"""
    describe(obj; <keyword arguments>)

Return basic descriptive statistics of the object data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `df::Bool=false`: if true, return statistics as a DataFrame

# Returns

- `Union{Nothing, DataFrame}`
"""
function describe(obj::NeuroAnalyzer.NEURO; df::Bool=false)::Union{Nothing, DataFrame}
    d = zeros(8, nchannels(obj))
    @inbounds for idx in 1:nchannels(obj)
        d[1, idx] = round(rng(obj.data[idx, :, :]), digits=2)
        d[2, idx] = round(mean(obj.data[idx, :, :]), digits=2)
        d[3, idx] = round(std(obj.data[idx, :, :]), digits=2)
        d[4, idx] = round(minimum(obj.data[idx, :, :]), digits=2)
        d[5, idx] = round(quantile(obj.data[idx, :, :][:], 0.5), digits=2)
        d[6, idx] = round(median(obj.data[idx, :, :]), digits=2)
        d[7, idx] = round(quantile(obj.data[idx, :, :][:], 0.95), digits=2)
        d[8, idx] = round(maximum(obj.data[idx, :, :]), digits=2)
    end
    if df
        df = DataFrame(:ch=>collect(1:nchannels(obj)),
                       :label=>labels(obj),
                       :type=>uppercase.(obj.header.recording[:channel_type]),
                       :unit=>obj.header.recording[:unit],
                       :range=>d[1, :],
                       :mean=>d[2, :],
                       :sd=>d[3, :],
                       :min=>d[4, :],
                       :Q1=>d[5, :],
                       :median=>d[6, :],
                       :Q3=>d[7, :],
                       :max=>d[8, :])
        return df
    else        
        println("< " * uppercase(obj.header.recording[:data_type]) * ", $(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)) ($(signal_len(obj) / sr(obj)) s) >")
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
        for idx in 1:nchannels(obj)
            println(rpad(string(idx), 4) *
                    rpad(labels(obj)[idx], 16) *
                    rpad(uppercase(obj.header.recording[:channel_type][idx]), 12) *
                    rpad(obj.header.recording[:unit][idx], 8) *
                    rpad(string(d[1, idx]), 10) *
                    rpad(string(d[2, idx]), 10) *
                    rpad(string(d[3, idx]), 10) *
                    rpad(string(d[4, idx]), 10) *
                    rpad(string(d[5, idx]), 10) *
                    rpad(string(d[6, idx]), 10) *
                    rpad(string(d[7, idx]), 10) *
                    rpad(string(d[8, idx]), 10))
        end
        return nothing
    end
end

"""
    size(obj)

Return size of the object data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `s::Tuple{Int64, Int64, Int64}`
"""
function Base.size(obj::NeuroAnalyzer.NEURO)::Tuple{Int64, Int64, Int64}

    s = size(obj.data)

    return s

end

"""
    size(obj, n)

Return size of the object data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `d::Int64`: compute size along dimension `d`

# Returns

- `s::Int64`
"""
function Base.size(obj::NeuroAnalyzer.NEURO, d::Int64)::Int64

    @assert d in 1:ndims(obj.data) "d must be in [1, $(ndims(obj.data))]."
    s = size(obj.data, d)

    return s

end

"""
    datatype(obj)

Return data type of the object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `dt::String`
"""
function datatype(obj::NeuroAnalyzer.NEURO)::String

    dt = obj.header.recording[:data_type]

    return dt

end

"""
    channel_order(obj)

Return channel order of the object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `co::Vector{Int64}`
"""
function channel_order(obj::NeuroAnalyzer.NEURO)::Vector{Int64}

    co = obj.header.recording[:channel_order]

    return co

end
