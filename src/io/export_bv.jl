export export_bv

# NeuroAnalyzer marker types → BrainVision marker types (read back by `_marker_type()`)
const _BV_MARKER_TYPE = Dict(
    "stim" => "Stimulus", "response" => "Response", "annotation" => "Comment",
    "boundary" => "New Segment", "BAD" => "Bad Interval", "epoch" => "Epoch",
    "other" => "Other",
)

# BrainVision text: commas in labels and descriptions are written as "\1"
_bv_escape(s::AbstractString)::String = replace(string(s), "," => "\\1", "\n" => " ")

# BrainVision unit: micro sign (U+00B5), as written by BrainVision Recorder
_bv_unit(u::AbstractString)::String = replace(string(u), "μ" => "µ")

# number for the header: integers without a decimal point, other values with up to 6 digits
function _bv_num(x::Real)::String
    isinteger(x) && return string(Int64(x))
    return string(round(Float64(x); sigdigits = 6))
end

"""
    export_bv(obj; <keyword arguments>)

Export a recording to the BrainVision Core Data Format (BrainVision Recorder/Analyzer, EEGLAB, MNE-Python, FieldTrip; one of the formats of BIDS EEG and iEEG): a text header (`.vhdr`), the signal (`.eeg`) and the markers (`.vmrk`).

- signal: 32-bit floating point values (`IEEE_FLOAT_32`), multiplexed, in the units of the channels (resolution 1)
- channels: labels, units (e.g. `µV`); the channel types are not stored by the format (`import_bv()` infers them from the labels)
- filters: the high-pass, low-pass and notch frequencies of the channels are written in the channel table of the `[Comment]` section, as BrainVision Recorder does (read by `import_bv()` and MNE-Python)
- markers: type (`stim` → `Stimulus`, `response` → `Response`, `annotation` → `Comment`, `boundary` → `New Segment`, `BAD` → `Bad Interval`, `epoch` → `Epoch`, `other` → `Other`), value (description), start and length (in samples), channel (one channel; markers of several channels refer to all channels); extra marker columns are not stored
- epoched objects: the epochs are written one after another, with a `New Segment` marker at the start of each epoch
- channel locations are not stored (use `export_locs()`)

# Arguments

- `obj::NeuroAnalyzer.NEURO`: EEG, ECoG, sEEG, iEEG, CSD or ERP object
- `file_name::String`: header file name (`.vhdr`); the signal and marker files get the same name with the extensions `.eeg` and `.vmrk`
- `overwrite::Bool=false`: if `true`, overwrite existing files

# Returns

- `Vector{String}`: names of the written files (header, signal, markers)
"""
function export_bv(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    overwrite::Bool = false,
)::Vector{String}
    # validate
    _check_export_datatype(obj, :bv)
    lowercase(splitext(file_name)[2]) == ".vhdr" ||
        throw(ArgumentError("File name extension must be .vhdr"))
    base = splitext(file_name)[1]
    files = [base * ".vhdr", base * ".eeg", base * ".vmrk"]
    for f in files
        (isfile(f) && !overwrite) && throw(
            ArgumentError("File $f cannot be saved, to overwrite use overwrite=true."),
        )
    end
    rec = obj.header.recording
    ch_n = nchannels(obj)
    fs = sr(obj)
    ep_n = nepochs(obj)
    ep_len = epoch_len(obj)
    clabels = labels(obj)
    units = string.(rec[:unit])
    nyq = fs / 2

    # ------------------------------------------------------------------ #
    # signal: channels × samples (all epochs), multiplexed float32       #
    # ------------------------------------------------------------------ #
    x = reshape(obj.data, ch_n, :)
    any(v -> isfinite(v) && abs(v) > floatmax(Float32), x) &&
        throw(ArgumentError("Values out of the Float32 range cannot be exported."))
    open(files[2], "w") do io
        write(io, htol.(Float32.(x)))
    end

    # ------------------------------------------------------------------ #
    # header                                                             #
    # ------------------------------------------------------------------ #
    hp = [rec[:highpass][i] > 0 ? _bv_num(rec[:highpass][i]) : "DC" for i = 1:ch_n]
    lp = [rec[:lowpass][i] < nyq ? _bv_num(rec[:lowpass][i]) : "Off" for i = 1:ch_n]
    nt = [isempty(rec[:notch][i]) ? "Off" : join(_bv_num.(rec[:notch][i]), "/") for i = 1:ch_n]
    w_l = max(6, maximum(length.(clabels)) + 2)
    open(files[1], "w") do io
        println(io, "Brain Vision Data Exchange Header File Version 1.0")
        println(io, "; Data exported by NeuroAnalyzer $(NeuroAnalyzer.VER)")
        println(io)
        println(io, "[Common Infos]")
        println(io, "Codepage=UTF-8")
        println(io, "DataFile=$(basename(files[2]))")
        println(io, "MarkerFile=$(basename(files[3]))")
        println(io, "DataFormat=BINARY")
        println(io, "; Data orientation: MULTIPLEXED=ch1,pt1, ch2,pt1 ...")
        println(io, "DataOrientation=MULTIPLEXED")
        println(io, "NumberOfChannels=$ch_n")
        println(io, "; Sampling interval in microseconds")
        println(io, "SamplingInterval=$(_bv_num(1e6 / fs))")
        if ep_n > 1
            println(io, "SegmentationType=MARKERBASED")
            println(io, "SegmentDataPoints=$ep_len")
        end
        println(io)
        println(io, "[Binary Infos]")
        println(io, "BinaryFormat=IEEE_FLOAT_32")
        println(io)
        println(io, "[Channel Infos]")
        println(io, "; Each entry: Ch<Channel number>=<Name>,<Reference channel name>,")
        println(io, "; <Resolution in \"Unit\">,<Unit>, Future extensions..")
        for i = 1:ch_n
            println(io, "Ch$i=$(_bv_escape(clabels[i])),,1,$(_bv_unit(units[i]))")
        end
        println(io)
        println(io, "[Comment]")
        println(io)
        println(io, "A m p l i f i e r  S e t u p")
        println(io, "============================")
        println(io, "Number of channels: $ch_n")
        println(io, "Sampling Rate [Hz]: $fs")
        println(io, "Sampling Interval [µS]: $(_bv_num(1e6 / fs))")
        println(io)
        println(io, "Channels")
        println(io, "--------")
        println(io, rpad("#", 6) * rpad("Name", w_l) * rpad("Phys. Chn.", 14) *
                    rpad("Resolution / Unit", 20) * rpad("Low Cutoff [Hz]", 17) *
                    rpad("High Cutoff [Hz]", 18) * "Notch [Hz]")
        for i = 1:ch_n
            println(io, rpad(string(i), 6) * rpad(clabels[i], w_l) * rpad(string(i), 14) *
                        # (a unit for every channel: readers split the rows at spaces)
                        rpad("1 " * (isempty(units[i]) ? "n/a" : _bv_unit(units[i])), 20) *
                        rpad(hp[i], 17) *
                        rpad(lp[i], 18) * nt[i])
        end
    end

    # ------------------------------------------------------------------ #
    # markers                                                            #
    # ------------------------------------------------------------------ #
    mk = DataFrame(:type => String[], :value => String[], :pos => Int64[],
                   :len => Int64[], :ch => Int64[])
    # epochs: a New Segment marker at the start of each epoch
    if ep_n > 1
        for e = 1:ep_n
            push!(mk, ("New Segment", "", (e - 1) * ep_len + 1, 1, 0))
        end
    end
    multi = 0
    for m in eachrow(obj.markers)
        pos = round(Int64, m.start * fs) + 1
        len = max(0, round(Int64, m.length * fs))
        chs = m.channel
        c = if length(chs) == 1
            something(findfirst(==(chs[1]), clabels), 0)
        else
            length(chs) > 1 && (multi += 1)
            0
        end
        push!(mk, (get(_BV_MARKER_TYPE, m.type, "Other"), string(m.value), pos, len, c))
    end
    sort!(mk, :pos)
    multi > 0 && _info("$multi marker(s) of several channels are written for all channels.")
    isempty(_marker_extra_cols(obj.markers)) ||
        _info("Extra marker columns ($(join(_marker_extra_cols(obj.markers), ", "))) are not exported.")
    open(files[3], "w") do io
        println(io, "Brain Vision Data Exchange Marker File, Version 1.0")
        println(io, "; Data exported by NeuroAnalyzer $(NeuroAnalyzer.VER)")
        println(io)
        println(io, "[Common Infos]")
        println(io, "Codepage=UTF-8")
        println(io, "DataFile=$(basename(files[2]))")
        println(io)
        println(io, "[Marker Infos]")
        println(io, "; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,")
        println(io, "; <Size in data points>, <Channel number (0 = marker is related to all channels)>")
        for (i, m) in enumerate(eachrow(mk))
            println(io, "Mk$i=$(_bv_escape(m.type)),$(_bv_escape(m.value)),$(m.pos),$(m.len),$(m.ch)")
        end
    end

    return files
end
