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

# internal helper: shared band-frequency table lookup
# returns the raw (bf_low, bf_high) tuple, or nothing for :list
function _band_table(band::Symbol, nqf::Float64)::Union{Tuple{Float64, Float64}, Nothing}
    band === :list         && return nothing
    band === :total        && return (0.1, round(nqf, digits=1))
    band === :delta        && return (0.1,  4.0)
    band === :theta        && return (4.0,  8.0)
    band === :alpha        && return (8.0, 13.0)
    band === :alpha_lower  && return (8.0, 10.5)
    band === :alpha_higher && return (10.5,13.0)
    band === :beta         && return (14.0,30.0)
    band === :beta_lower   && return (14.0,25.0)
    band === :beta_higher  && return (25.0,30.0)
    band === :gamma        && return (30.0,150.0)
    band === :gamma_1      && return (30.0, 40.0)
    band === :gamma_2      && return (40.0, 50.0)
    band === :gamma_lower  && return (30.0, 80.0)
    band === :gamma_higher && return (80.0,150.0)
end


"""
    sr(obj)

Return the sampling rate of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Int64`: sampling rate in Hz
"""
function sr(obj::NeuroAnalyzer.NEURO)::Int64

    return obj.header.recording[:sampling_rate]

end

"""
    nchannels(obj; <keyword arguments>)

Return the number of channels of a given type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `type::String="all"`: channel type string (must be one of the values in the global `channel_types` constant); use `"all"` to count every channel.

# Returns

- `Int64`: number of channels of the requested type

# Throws

- `ArgumentError`: if `type` is not in `channel_types`, or if no channel types are defined in the object.
"""
function nchannels(obj::NeuroAnalyzer.NEURO; type::String = "all")::Int64

    _check_var(type, channel_types, "type")
    !(length(obj.header.recording[:channel_type]) != 0) && throw(ArgumentError("OBJ has no defined channel types."))

    if type == "all"
        return size(obj.data, 1)
    else
        return count(==(type), obj.header.recording[:channel_type])
    end

end

"""
    nepochs(obj)

Return the number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Int64`: number of epochs (third data dimension)
"""
function nepochs(obj::NeuroAnalyzer.NEURO)::Int64

    return size(obj.data, 3)

end

"""
    signal_len(obj)

Return the total signal length in samples (across all epochs).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Int64`: total number of samples; for a 3-D array this equals `epoch_len × nepochs`
"""
function signal_len(obj::NeuroAnalyzer.NEURO)::Int64

    return size(obj.data, 2) * size(obj.data, 3)

end

"""
    epoch_len(obj)

Return the epoch length in samples.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Int64`: number of samples per epoch (second data dimension)
"""
function epoch_len(obj::NeuroAnalyzer.NEURO)::Int64

    return size(obj.data, 2)

end

"""
    signal_duration(obj)

Return the total signal duration in seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Float64`: duration in seconds
"""
function signal_duration(obj::NeuroAnalyzer.NEURO)::Float64

    return signal_len(obj) / sr(obj) 

end

"""
    epoch_duration(obj)

Return the epoch duration in seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Float64`: duration of one epoch in seconds
"""
function epoch_duration(obj::NeuroAnalyzer.NEURO)::Float64

    return size(obj, 2) / sr(obj) 

end

"""
    history(obj)

Return the processing history log.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Vector{String}`: list of processing steps applied to the object
"""
function history(obj::NeuroAnalyzer.NEURO)::Vector{String}

    return obj.history

end

"""
    labels(obj)

Return channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Vector{String}`: channel label strings

# Throws

- `ArgumentError`: if the object has no channel labels
"""
function labels(obj::NeuroAnalyzer.NEURO)::Vector{String}

    !(length(obj.header.recording[:label]) > 0) && throw(ArgumentError("OBJ has no channel labels."))
    return obj.header.recording[:label]

end

"""
    optode_labels(obj)

Return optode labels (NIRS objects only).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"nirs"`

# Returns

- `Vector{String}`: optode label strings

# Throws

- `ArgumentError`: if the object is not NIRS, or has no optode labels
"""
function optode_labels(obj::NeuroAnalyzer.NEURO)::Vector{String}

    _check_datatype(obj, "nirs")
    !(length(obj.header.recording[:optode_labels]) > 0) && throw(ArgumentError("OBJ has no optode labels."))
    return obj.header.recording[:optode_labels]

end

"""
    source_labels(obj)

Return source labels (NIRS objects only).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"nirs"`

# Returns

- `Vector{String}`: source label strings

# Throws

- `ArgumentError`: if the object is not NIRS, or has no source labels
"""
function source_labels(obj::NeuroAnalyzer.NEURO)::Vector{String}

    _check_datatype(obj, "nirs")
    !(length(obj.header.recording[:src_labels]) > 0) && throw(ArgumentError("OBJ has no source labels."))
    return obj.header.recording[:src_labels]

end

"""
    detector_labels(obj)

Return detector labels (NIRS objects only).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"nirs"`

# Returns

- `Vector{String}`: detector label strings

# Throws

- `ArgumentError`: if the object is not NIRS, or has no source labels
"""
function detector_labels(obj::NeuroAnalyzer.NEURO)::Vector{String}

    _check_datatype(obj, "nirs")
    !(length(obj.header.recording[:det_labels]) > 0) && throw(ArgumentError("OBJ has no detector labels."))
    return obj.header.recording[:det_labels]

end

"""
    chtypes(obj)

Return channel type strings.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Vector{String}`: channel type strings (one per channel)

# Throws

- `ArgumentError`: if no channel types are defined in the object
"""
function chtypes(obj::NeuroAnalyzer.NEURO)::Vector{String}

    !(length(obj.header.recording[:channel_type]) > 0) && throw(ArgumentError("OBJ has no channel types."))
    return obj.header.recording[:channel_type]

end

# internal helper: print the common header/info block
# both `header()` and `info()` display identical metadata; extracted here to avoid duplication
function _print_header(obj::NeuroAnalyzer.NEURO)
    println("              Data type: $(uppercase(obj.header.recording[:data_type]))")
    println("            File format: $(obj.header.recording[:file_type])")
    println("            Source file: $(obj.header.recording[:file_name])")
    println("         File size [MB]: $(obj.header.recording[:file_size_mb])")
    println("       Memory size [MB]: $(round(Base.summarysize(obj) / 1024^2, digits=2))")

    subj = obj.header.subject
    subj_str = if length(subj[:id]) > 0
        "$(subj[:id]): $(subj[:first_name]) $(subj[:last_name])"
    else
        "$(subj[:first_name]) $(subj[:last_name])"
    end
    println("                Subject: $subj_str")

    rec = obj.header.recording
    println("              Recording: $(rec[:recording])")
    println("        Recording notes: $(rec[:recording_notes])")
    println("         Recording date: $(rec[:recording_date])")
    println("         Recording time: $(rec[:recording_time])")
    println("     Sampling rate (Hz): $(sr(obj))")
    println("Signal length [samples]: $(signal_len(obj))")
    println("Signal length [seconds]: $(round(signal_len(obj) / sr(obj), digits=4))")
    println("     Number of channels: $(nchannels(obj))")

    if !(datatype(obj) in ["mep", "sensors", "eda"])
        println("              Epochs ID: $(rec[:epoch_id])")
        println("       Number of epochs: $(nepochs(obj))")
        println(" Epoch length [samples]: $(epoch_len(obj))")
        println(" Epoch length [seconds]: $(round(epoch_len(obj) / sr(obj), digits=4))")
    end

    if datatype(obj) == "eeg"
        ref = rec[:reference] == "" ? "unknown" : rec[:reference]
        println("         Reference type: $ref")
    end

    if datatype(obj) == "meg"
        ssp_labels = rec[:ssp_labels]
        if length(ssp_labels) == 1
            println("       SSP projection: $(ssp_labels[1])")
        elseif length(ssp_labels) > 1
            print("        SSP projections: ")
            for idx in 1:(length(ssp_labels) - 1)
                print("$(ssp_labels[idx]), ")
            end
            println("$(ssp_labels[end])")
        end
    end

    if datatype(obj) in ["eeg", "meg", "ecog", "seeg"]
        println("         Line frequency: $(rec[:line_frequency]) Hz")
    end
    if datatype(obj) == "nirs"
        println("        Wavelength [nm]: $(rec[:wavelengths])")
    end

    if !(datatype(obj) in ["mep", "sensors", "eda"])
        println("                Markers: $(_has_markers(obj) ? "yes" : "no")")
        println("      Channel locations: $(DataFrames.nrow(obj.locs) > 0 ? "yes" : "no")")
    end

    # channel-type counts, grouped by modality
    dt = datatype(obj)
    if dt in ["eeg", "ecog", "seeg", "ieeg", "erp"]
        nch = count(x -> x in ["eeg", "ecog", "seeg", "ieeg", "erp"], rec[:channel_type])
        println(" Number of EEG channels: $nch")
    elseif dt in ["meg", "erf"]
        println(" Number of MAG channels: $(count(==("mag"),  rec[:channel_type]))")
        println("Number of GRAD channels: $(count(==("grad"), rec[:channel_type]))")
        println(" Number of EEG channels: $(count(==("eeg"),  rec[:channel_type]))")
    elseif dt == "nirs"
        nirs_types = [
            "nirs", "nirs_int", "nirs_od", "nirs_dmean", "nirs_dvar", "nirs_dskew",
            "nirs_mua", "nirs_musp", "nirs_hbo", "nirs_hbr", "nirs_hbt", "nirs_h2o",
            "nirs_lipid", "nirs_bfi", "nirs_hrf_dod", "nirs_hrf_dmean", "nirs_hrf_dvar",
            "nirs_hrf_dskew", "nirs_hrf_hbo", "nirs_hrf_hbr", "nirs_hrf_hbt",
            "nirs_hrf_bfi", "nirs_aux",
        ]
        nch = count(x -> x in nirs_types, rec[:channel_type])
        println("Number of NIRS channels: $nch")
    end
end

"""
    header(obj; <keyword arguments>)

Print object header metadata to stdout.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Nothing`

# See also
[`info`](@ref)
"""
function header(obj::NeuroAnalyzer.NEURO)::Nothing

    _print_header(obj)
    return nothing

end

"""
    info(obj; <keyword arguments>)

Print object metadata and channel table. Optionally return data as a DataFrame.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `df::Bool=false`: if `true`, return the signal data reshaped into a `DataFrame` with columns `:time` and one column per channel label

# Returns

- `Nothing` if `df=false`, or a `DataFrame` if `df=true`.

# See also

[`header`](@ref), [`describe`](@ref)
"""
function info(obj::NeuroAnalyzer.NEURO; df::Bool = false)::Union{Nothing, DataFrame}

    _print_header(obj)
    println()
    println("Channels:")

    rec = obj.header.recording
    if rec[:data_type] != "nirs"
        # standard (non-NIRS) channel table
        println(rpad(" ch", 8) * rpad("label", 16) * rpad("type", 12) *
                rpad("unit", 8) * rpad("bad", 8))
        println(" " * repeat("-", 6) * " " * repeat("-", 15) * " " *
                repeat("-", 11) * " " * repeat("-", 7) * " " * repeat("-", 7))
        for idx in eachindex(rec[:label])
            println(rpad(" $idx", 8) *
                    rpad(rec[:label][idx], 16) *
                    rpad(uppercase(rec[:channel_type][idx]), 12) *
                    rpad(rec[:unit][idx], 8) *
                    rpad(string(rec[:bad_channel][idx]), 8))
        end
    else
        # NIRS channel table (includes wavelength column for non-derived channels)
        println(rpad(" ch", 8) * rpad("label", 16) * rpad("type", 12) *
                rpad("unit", 8) * rpad("wavelength", 12))
        println(" " * repeat("-", 6) * " " * repeat("-", 15) * " " *
                repeat("-", 11) * " " * repeat("-", 7) * " " * repeat("-", 12))
        derived = ["nirs_aux", "nirs_hbo", "nirs_hbr", "nirs_hbt"]
        for idx in eachindex(rec[:label])
            if !(rec[:channel_type][idx] in derived)
                wl = rec[:wavelengths][rec[:wavelength_index][idx]]
                println(rpad(" $idx", 8) *
                        rpad(rec[:label][idx], 16) *
                        rpad(uppercase(rec[:channel_type][idx]), 12) *
                        rpad(rec[:unit][idx], 8) *
                        rpad(string(wl), 12))
            else
                println(rpad(" $idx", 8) *
                        rpad(rec[:label][idx], 16) *
                        rpad(uppercase(rec[:channel_type][idx]), 12))
            end
        end
    end

    if df
        result = DataFrame(
            hcat(obj.time_pts, reshape(obj.data, nchannels(obj), :, 1)[:, :]'),
            :auto,
        )
        DataFrames.rename!(result, vcat(:time, Symbol.(labels(obj))))
        return result
    else
        return nothing
    end

end


"""
    channel_info(obj; <keyword arguments>)

Return or print information for a single channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `pr::Bool=true`: if `true`, print to stdout and return `nothing`; if `false`, return the info string

# Returns

- `Nothing` when `pr=true`, or `String` when `pr=false`.
"""
function channel_info(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    pr::Bool = true
)::Union{Nothing, String}

    ch = get_channel(obj, ch=ch)
    !(length(ch) == 1) && throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    rec = obj.header.recording
    local chi::String

    if rec[:data_type] != "nirs"

        chi = " ch: $(rpad(string(ch), 4))" *
              " label: $(rpad(rec[:label][ch], 8))" *
              " type: $(rpad(uppercase(rec[:channel_type][ch]), 8))" *
              " unit: $(rpad(rec[:unit][ch], 8))" *
              " bad: $(rec[:bad_channel][ch])"

    elseif rec[:channel_type][ch] != "nirs_aux"

        chi = " ch: $(rpad(string(ch), 4))" *
              " label: $(rpad(rec[:label][ch], 8))" *
              " type: $(rpad(uppercase(rec[:channel_type][ch]), 8))" *
              " unit: $(rpad(rec[:unit][ch], 8))" *
              " wavelength: $(rec[:wavelength_index][ch])"

    else

        chi = " ch: $(rpad(string(ch), 4))" *
              " label: $(rpad(rec[:label][ch], 8))" *
              " type: $(rpad(uppercase(rec[:channel_type][ch]), 8))"

    end

    if pr
        println(chi)
        return nothing
    else
        return chi
    end

end

"""
    channel_pick(obj; <keyword arguments>)

Return set of channel indices corresponding to a set of electrodes ("pick", e.g. left or frontal electrodes).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`
- `pick::Vector{Symbol}`: pick of electrodes; picks may be combined, e.g. `[:left, :frontal]`
    - `:central` / `:c`
    - `:left` / `:l`
    - `:right` / `:r`
    - `:frontal` / `:f`
    - `:temporal` / `:t`
    - `:parietal` / `:p`
    - `:occipital` / `:o`

# Returns

- `Vector{String}`: channel names matching the pick

# Throws

- `ArgumentError`: if `pick` contains an unrecognised symbol or the object has no channel labels
"""
function channel_pick(
    obj::NeuroAnalyzer.NEURO;
    pick::Union{Symbol, Vector{Symbol}}
)::Vector{String}

    _check_datatype(obj, "eeg")
    !(length(labels(obj)) != 0) && throw(ArgumentError("OBJ does not contain channel labels."))

    valid = [:central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o]

    if pick isa Vector{Symbol}

        for idx in pick
            _check_var(idx, valid, "pick")
        end

        # map non-lateralised picks to their standard 10-20 letter prefix
        c = Vector{Char}()
        for idx in pick
            (idx === :central  || idx === :c) && push!(c, 'z')
            (idx === :frontal  || idx === :f) && push!(c, 'F')
            (idx === :temporal || idx === :t) && push!(c, 'T')
            (idx === :parietal || idx === :p) && push!(c, 'P')
            (idx === :occipital|| idx === :o) && push!(c, 'O')
        end

        clabels = get_channel(obj, type="eeg")
        ch = Vector{Int64}()
        for idx1 in eachindex(clabels), idx2 in eachindex(c)
            in(c[idx2], clabels[idx1]) && push!(ch, idx1)
        end

        # when both :left and :right are requested simultaneously, return the full bilateral set without further laterality filtering
        has_left  = any(p -> p === :left  || p === :l, pick)
        has_right = any(p -> p === :right || p === :r, pick)
        has_left && has_right && return labels(obj)[ch]

        # single-laterality filtering: remove contralateral electrode numbers
        clabels = get_channel(obj, type="eeg")[ch]
        pat = nothing
        for idx in pick
            # remove left-side numbers
            (idx === :right || idx === :r) && (pat = r"[z13579]$")
            # remove right-side numbers
            (idx === :left  || idx === :l) && (pat = r"[z02468]$")
        end
        if pat isa Regex
            for idx in length(ch):-1:1
                !isnothing(match(pat, clabels[idx])) && deleteat!(ch, idx)
            end
        end

        return labels(obj)[ch]

    else

        _check_var(pick, valid, "pick")

        # map single pick to the relevant 10-20 digit/letter characters
        c = Vector{Char}()
        (pick === :central  || pick === :c) && (c = ['z'])
        (pick === :left     || pick === :l) && (c = ['1', '3', '5', '7', '9'])
        (pick === :right    || pick === :r) && (c = ['2', '4', '6', '8'])
        (pick === :frontal  || pick === :f) && (c = ['F'])
        (pick === :temporal || pick === :t) && (c = ['T'])
        (pick === :parietal || pick === :p) && (c = ['P'])
        (pick === :occipital|| pick === :o) && (c = ['O'])

        clabels = get_channel(obj, type="eeg")
        ch = Vector{Int64}()
        for idx1 in eachindex(c), idx2 in eachindex(clabels)
            in(c[idx1], clabels[idx2]) && push!(ch, idx2)
        end

        return labels(obj)[ch]

    end

end


"""
    channels_cluster(obj; <keyword arguments>)

Return channel names belonging to a predefined spatial cluster.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `cluster::Symbol`: cluster name:
    - `:f1`: left frontal  (Fp1, F1, F3, F5, F7, F9, AF3, AF7)
    - `:f2`: right frontal (Fp2, F2, F4, F6, F8, F10, AF4, AF8)
    - `:t1`: left temporal  (C3, C5, T7, T9, FC3, FC5, FT7, FT9)
    - `:t2`: right temporal (C4, C6, T8, T10, FC4, FC6, FT8, FT10)
    - `:c1`: anterior central (Cz, C1, C2, FC1, FC2, FCz)
    - `:c2`: posterior central (Pz, P1, P2, CP1, CP2, CPz)
    - `:p1`: left parietal  (P3, P5, P7, P9, CP3, CP5, TP7, TP9)
    - `:p2`: right parietal (P4, P6, P8, P10, CP4, CP6, TP8, TP10)
    - `:o`:  occipital (O1, O2, POz, PO3, PO4, PO7, PO8, PO9, PO10)

# Returns

- `Vector{String}`: channel names present in the object that belong to the cluster

# Throws

- `ArgumentError`: if `cluster` is not a recognised symbol or the object has no labels.
"""
function channel_cluster(obj::NeuroAnalyzer.NEURO; cluster::Symbol)::Vector{String}

    !(length(labels(obj)) != 0) && throw(ArgumentError("OBJ does not contain channel labels."))
    _check_var(cluster, [:f1, :f2, :t1, :t2, :c1, :c2, :p1, :p2, :o], "cluster")

    clabels = labels(obj)

    # Map cluster symbol to the standard electrode names it encompasses
    cluster_map = Dict(
        :f1 => ["Fp1", "F1",  "F3",  "F5",  "F7",  "F9",  "AF3", "AF7"],
        :f2 => ["Fp2", "F2",  "F4",  "F6",  "F8",  "F10", "AF4", "AF8"],
        :t1 => ["C3",  "C5",  "T7",  "T9",  "FC3", "FC5", "FT7", "FT9"],
        :t2 => ["C4",  "C6",  "T8",  "T10", "FC4", "FC6", "FT8", "FT10"],
        :c1 => ["Cz",  "C1",  "C2",  "FC1", "FC2", "FCz"],
        :c2 => ["Pz",  "P1",  "P2",  "CP1", "CP2", "CPz"],
        :p1 => ["P3",  "P5",  "P7",  "P9",  "CP3", "CP5", "TP7", "TP9"],
        :p2 => ["P4",  "P6",  "P8",  "P10", "CP4", "CP6", "TP8", "TP10"],
        :o  => ["O1",  "O2",  "POz", "PO3", "PO4", "PO7", "PO8", "PO9", "PO10"],
    )

    return Base.filter(l -> l in clabels, cluster_map[cluster])

end

# internal helper: print the :list of available bands
function _band_list(bands::Vector{Symbol})
    print("Available band names: ")
    for x in 2:(length(bands) - 1)
        print(":$(bands[x]), ")
    end
    println(":$(bands[end])")
end

# internal helper: clamp bf to the Nyquist frequency with a warning
function _clamp_band(bf::Tuple{Float64, Float64}, nqf::Float64, band::Symbol, label::String)
    bf_low, bf_high = bf
    if bf_low > nqf
        _warn("Nyquist frequency ($nqf Hz) is lower than $band range: $bf. " *
              "Band truncated to: ($(nqf - 0.2), $(nqf - 0.1)).")
        return (nqf - 0.2, nqf - 0.1)
    end
    if bf_high > nqf
        _warn("Nyquist frequency ($nqf Hz) is lower than $band range: $bf. " *
              "Band truncated to: ($bf_low, $(nqf - 0.1)).")
        return (bf_low, nqf - 0.1)
    end
    return bf
end

"""
    band_frq(obj, band)

Return the frequency limits for a named EEG band.

When `band = :list`, the available band names are printed to stdout and the function returns `(0.0, 0.0)` as a sentinel (rather than `nothing`, which would violate the declared return type).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `band::Symbol`: band name (see list below or pass `:list` to print all names):
    - `:list`
    - `:total`
    - `:delta`: 0.1-4.0 Hz
    - `:theta`: 4.0-8.0 Hz
    - `:alpha`: 8.0-3.0 Hz
    - `:alpha_lower`: 8.0-10.5 Hz
    - `:alpha_higher`: 10.5-13.0 Hz
    - `:beta`: 14.0-30.0 Hz
    - `:beta_lower`: 14.0-25.0 Hz
    - `:beta_higher`: 25.0-30.0 Hz
    - `:gamma`: 30.0-150.0 Hz
    - `:gamma_1`: 30.0-40.0 Hz
    - `:gamma_2`: 40.0-50.0 Hz
    - `:gamma_lower`: 30.0-80.0 Hz
    - `:gamma_higher`: 80.0-150.0 Hz

# Returns

- `Tuple{Float64, Float64}`: `(low_Hz, high_Hz)` limits, clamped to the Nyquist frequency if necessary

# See also

[`band_frq(::Int64)`](@ref)
"""
function band_frq(obj::NeuroAnalyzer.NEURO; band::Symbol)::Tuple{Float64, Float64}

    bands = [:list, :total, :delta, :theta, :alpha, :alpha_lower, :alpha_higher,
             :beta, :beta_lower, :beta_higher, :gamma, :gamma_1, :gamma_2,
             :gamma_lower, :gamma_higher]
    _check_var(band, bands, "band")

    nqf = sr(obj) / 2.0
    bf = _band_table(band, nqf)

    if isnothing(bf)
        _band_list(bands)
        # sentinel: caller should check for :list before using the result
        return (0.0, 0.0)
    end

    return _clamp_band(bf, nqf, band, "obj")

end

"""
    band_frq(fs, band)

Return the frequency limits for a named EEG band.

When `band = :list`, the available band names are printed to stdout and the function returns `(0.0, 0.0)` as a sentinel.

# Arguments

- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `band::Symbol`: band name (see list below or pass `:list` to print all names):
    - `:list`
    - `:total`
    - `:delta`: 0.1-4.0 Hz
    - `:theta`: 4.0-8.0 Hz
    - `:alpha`: 8.0-3.0 Hz
    - `:alpha_lower`: 8.0-10.5 Hz
    - `:alpha_higher`: 10.5-13.0 Hz
    - `:beta`: 14.0-30.0 Hz
    - `:beta_lower`: 14.0-25.0 Hz
    - `:beta_higher`: 25.0-30.0 Hz
    - `:gamma`: 30.0-150.0 Hz
    - `:gamma_1`: 30.0-40.0 Hz
    - `:gamma_2`: 40.0-50.0 Hz
    - `:gamma_lower`: 30.0-80.0 Hz
    - `:gamma_higher`: 80.0-150.0 Hz

# Returns

- `Tuple{Float64, Float64}`: `(low_Hz, high_Hz)` limits, clamped to the Nyquist frequency if necessary

# See also

[`band_frq(::NeuroAnalyzer.NEURO)`](@ref)
"""
function band_frq(fs::Int64; band::Symbol)::Tuple{Float64, Float64}

    bands = [:list, :total, :delta, :theta, :alpha, :alpha_lower, :alpha_higher,
             :beta, :beta_lower, :beta_higher, :gamma, :gamma_1, :gamma_2,
             :gamma_lower, :gamma_higher]
    _check_var(band, bands, "band")

    nqf = fs / 2.0
    bf = _band_table(band, nqf)

    if isnothing(bf)
        _band_list(bands)
        # sentinel: caller should check for :list before using the result
        return (0.0, 0.0)
    end

    return _clamp_band(bf, nqf, band, "fs=$fs")
end

"""
    describe(obj; <keyword arguments>)

Print or return descriptive statistics for each channel.

Statistics reported: range, mean, SD, minimum, Q1 (25th percentile), median, Q3 (75th percentile), maximum.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `df::Bool=false`: if `true`, return statistics as a `DataFrame`; otherwise print a formatted table and return `nothing`

# Returns

- `Nothing` when `df=false`, or a `DataFrame` when `df=true`

# See also

[`info`](@ref)
"""
function describe(obj::NeuroAnalyzer.NEURO; df::Bool = false)::Union{Nothing, DataFrame}

    d = zeros(8, nchannels(obj))
    @inbounds for idx in 1:nchannels(obj)
        d[1, idx] = round(rng(obj.data[idx, :, :]), digits = 2)
        d[2, idx] = round(mean(obj.data[idx, :, :]), digits = 2)
        d[3, idx] = round(std(obj.data[idx, :, :]), digits = 2)
        d[4, idx] = round(minimum(obj.data[idx, :, :]), digits = 2)
        d[5, idx] = round(quantile(obj.data[idx, :, :][:], 0.5), digits = 2)
        d[6, idx] = round(median(obj.data[idx, :, :]), digits = 2)
        d[7, idx] = round(quantile(obj.data[idx, :, :][:], 0.95), digits = 2)
        d[8, idx] = round(maximum(obj.data[idx, :, :]), digits = 2)
    end

    if df

        df = DataFrame(
            :ch => collect(1:nchannels(obj)),
            :label => labels(obj),
            :type => uppercase.(obj.header.recording[:channel_type]),
            :unit => obj.header.recording[:unit],
            :range => d[1, :],
            :mean => d[2, :],
            :sd => d[3, :],
            :min => d[4, :],
            :Q1 => d[5, :],
            :median => d[6, :],
            :Q3 => d[7, :],
            :max => d[8, :],
        )

        return df

    else

        println("< $(uppercase(obj.header.recording[:data_type])), " *
                "$(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)) " *
                "($(signal_len(obj) / sr(obj)) s) >")
        hdr = rpad("ch", 4) * rpad("label", 16) * rpad("type", 12) * rpad("unit", 8) *
              rpad("range", 10) * rpad("mean", 10) * rpad("sd", 10) *
              rpad("min", 10) * rpad("Q1", 10) * rpad("median", 10) *
              rpad("Q3", 10) * rpad("max", 10)
        println(hdr)
        for idx in 1:nchannels(obj)
            println(
                rpad(string(idx), 4) *
                rpad(labels(obj)[idx], 16) *
                rpad(uppercase(obj.header.recording[:channel_type][idx]), 12) *
                rpad(obj.header.recording[:unit][idx], 8) *
                Base.join(rpad.(string.(d[:, idx]), 10)),
            )
        end

        return nothing

    end
end

"""
    size(obj)

Return the size of the object data array.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Tuple{Int64, Int64, Int64}`: `(channels, samples, epochs)`.
"""
function Base.size(obj::NeuroAnalyzer.NEURO)::Tuple{Int64, Int64, Int64}

    return size(obj.data)

end

"""
    size(obj, d)

Return the size of the object data array along dimension `d`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `d::Int64`: dimension index; must be in `[1, ndims(obj.data)]`

# Returns

- `Int64`: size along dimension `d`

# Throws

- `ArgumentError`: If `d` is out of range.
"""
function Base.size(obj::NeuroAnalyzer.NEURO, d::Int64)::Int64

    !(d in 1:ndims(obj.data)) && throw(ArgumentError("d must be in [1, $(ndims(obj.data))]."))
    return size(obj.data, d)

end

"""
    datatype(obj)

Return the data type string of the object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `String`: data type (e.g. `"eeg"`, `"meg"`, `"nirs"`).
"""
function datatype(obj::NeuroAnalyzer.NEURO)::String

    return obj.header.recording[:data_type]

end

"""
    channel_order(obj)

Return the channel order indices stored in the object header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Vector{Int64}`: channel order indices
"""
function channel_order(obj::NeuroAnalyzer.NEURO)::Vector{Int64}

    return obj.header.recording[:channel_order]

end
