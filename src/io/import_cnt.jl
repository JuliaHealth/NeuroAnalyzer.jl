export import_cnt

"""
    import_cnt(file_name; <keyword arguments>)

Load a Neuroscan continuous signal (CNT) file and return a
`NeuroAnalyzer.NEURO` object.

CNT files store a fixed 900-byte global header, a per-channel header block (75 bytes × number of channels), multiplexed signal data, and an event table. Three event-table variants are supported (TEEG types 1, 2, and 3).

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: infer channel type from channel label
- `data_format::Symbol=:i32`: sample width; `:i16` for 16-bit or `:i32` for 32-bit; Neuroscan does not encode this in the file header - you must supply it based on your recording setup

# Returns

- `NeuroAnalyzer.NEURO`

# Notes

Based on `loadcnt.m` by Sean Fitzgibbon and Arnaud Delorme (https://cnl.salk.edu/~arno/cntload/index.html).

# Throws

- `ArgumentError` if the file does not exist or is not a `.cnt` file
"""
function import_cnt(
    file_name::String;
    data_format::Symbol = :i32,
    detect_type::Bool = true
)::NeuroAnalyzer.NEURO

    _check_var(data_format, [:i32, :i16], "data_format")
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) == ".cnt" ||
        throw(ArgumentError("$file_name is not a CNT file."))

    # ------------------------------------------------------------------ #
    # parse fixed-length global header (256 bytes)                       #
    # ------------------------------------------------------------------ #
    # all header reads share a single open/close via the `do` block
    imported_object = open(file_name, "r") do fid

        file_type = "CNT"

        # ------------------------------------------------------------ #
        # global 900-byte header                                        #
        # ------------------------------------------------------------ #
        rev = _vint2str(_fread(fid, 12, :c))
        rev != "Version 3.0" && _warn(
            "CNT files with version other than 3.0 may not import correctly; " *
            "please send this file to adam.wysokinski@neuroanalyzer.org")

        nextfile          = _fread(fid, 1, :l)
        prevfile          = _fread(fid, 1, :ul) # high-order bits of event table offset
        hdr_type          = _fread(fid, 1, :c) # renamed: avoids collision with event `evt_type`
        id                = _vint2str(_fread(fid, 20, :c))
        oper              = _vint2str(_fread(fid, 20, :c))
        doctor            = _vint2str(_fread(fid, 20, :c))
        referral          = _vint2str(_fread(fid, 20, :c))
        hospital          = _vint2str(_fread(fid, 20, :c))
        patient           = _vint2str(_fread(fid, 20, :c))
        age               = _fread(fid, 1, :i16)
        sex               = _fread(fid, 1, :c)
        hand              = _fread(fid, 1, :c)
        med               = _vint2str(_fread(fid, 20, :c))
        category          = _vint2str(_fread(fid, 20, :c))
        state             = _vint2str(_fread(fid, 20, :c))
        label             = _vint2str(_fread(fid, 20, :c))
        recording_date    = _vint2str(_fread(fid, 10, :c))
        recording_time    = _vint2str(_fread(fid, 12, :c))
        mean_age          = _fread(fid, 1, :f32)
        stdev             = _fread(fid, 1, :f32)
        n                 = _fread(fid, 1, :i16)
        compfile          = _fread(fid, 38, :c)
        spectwincomp      = _fread(fid, 1, :f32)
        meanaccuracy      = _fread(fid, 1, :f32)
        meanlatency       = _fread(fid, 1, :f32)
        sortfile          = _fread(fid, 46, :c)
        numevents         = _fread(fid, 1, :i)
        compoper          = _fread(fid, 1, :c)
        avgmode           = _fread(fid, 1, :c)
        review            = _fread(fid, 1, :c)
        nsweeps           = _fread(fid, 1, :ui16)
        compsweeps        = _fread(fid, 1, :ui16)
        acceptcnt         = _fread(fid, 1, :ui16)
        rejectcnt         = _fread(fid, 1, :ui16)
        pnts              = _fread(fid, 1, :ui16)
        ch_n              = _fread(fid, 1, :ui16)
        avgupdate         = _fread(fid, 1, :ui16)
        domain            = _fread(fid, 1, :c)
        variance          = _fread(fid, 1, :c)
        sampling_rate     = _fread(fid, 1, :ui16)
        scale             = _fread(fid, 1, :f64)
        veogcorrect       = _fread(fid, 1, :c)
        heogcorrect       = _fread(fid, 1, :c)
        aux1correct       = _fread(fid, 1, :c)
        aux2correct       = _fread(fid, 1, :c)
        veogtrig          = _fread(fid, 1, :f32)
        heogtrig          = _fread(fid, 1, :f32)
        aux1trig          = _fread(fid, 1, :f32)
        aux2trig          = _fread(fid, 1, :f32)
        heogchnl          = _fread(fid, 1, :i16)
        veogchnl          = _fread(fid, 1, :i16)
        aux1chnl          = _fread(fid, 1, :i16)
        aux2chnl          = _fread(fid, 1, :i16)
        veogdir           = _fread(fid, 1, :c)
        heogdir           = _fread(fid, 1, :c)
        aux1dir           = _fread(fid, 1, :c)
        aux2dir           = _fread(fid, 1, :c)
        veog_n            = _fread(fid, 1, :i16)
        heog_n            = _fread(fid, 1, :i16)
        aux1_n            = _fread(fid, 1, :i16)
        aux2_n            = _fread(fid, 1, :i16)
        veogmaxcnt        = _fread(fid, 1, :i16)
        heogmaxcnt        = _fread(fid, 1, :i16)
        aux1maxcnt        = _fread(fid, 1, :i16)
        aux2maxcnt        = _fread(fid, 1, :i16)
        veogmethod        = _fread(fid, 1, :c)
        heogmethod        = _fread(fid, 1, :c)
        aux1method        = _fread(fid, 1, :c)
        aux2method        = _fread(fid, 1, :c)
        ampsensitivity    = _fread(fid, 1, :f32)
        lowpass           = _fread(fid, 1, :c)
        highpass          = _fread(fid, 1, :c)
        notch             = _fread(fid, 1, :c)
        autoclipadd       = _fread(fid, 1, :c)
        baseline          = _fread(fid, 1, :c)
        offstart          = _fread(fid, 1, :f32)
        offstop           = _fread(fid, 1, :f32)
        reject            = _fread(fid, 1, :c)
        rejstart          = _fread(fid, 1, :f32)
        rejstop           = _fread(fid, 1, :f32)
        rejmin            = _fread(fid, 1, :f32)
        rejmax            = _fread(fid, 1, :f32)
        trigtype          = _fread(fid, 1, :c)
        trigval           = _fread(fid, 1, :f32)
        trigchnl          = _fread(fid, 1, :c)
        trigmask          = _fread(fid, 1, :i16)
        trigisi           = _fread(fid, 1, :f32)
        trigmin           = _fread(fid, 1, :f32)
        trigmax           = _fread(fid, 1, :f32)
        trigdir           = _fread(fid, 1, :c)
        autoscale         = _fread(fid, 1, :c)
        n2                = _fread(fid, 1, :i16)
        dir               = _fread(fid, 1, :c)
        dispmin           = _fread(fid, 1, :f32)
        dispmax           = _fread(fid, 1, :f32)
        xmin              = _fread(fid, 1, :f32)
        xmax              = _fread(fid, 1, :f32)
        automin           = _fread(fid, 1, :f32)
        automax           = _fread(fid, 1, :f32)
        zmin              = _fread(fid, 1, :f32)
        zmax              = _fread(fid, 1, :f32)
        lowcut            = _fread(fid, 1, :f32)
        highcut           = _fread(fid, 1, :f32)
        common            = _fread(fid, 1, :c)
        savemode          = _fread(fid, 1, :c)
        manmode           = _fread(fid, 1, :c)
        ref               = _fread(fid, 10, :c)
        rectify           = _fread(fid, 1, :c)
        displayxmin       = _fread(fid, 1, :f32)
        displayxmax       = _fread(fid, 1, :f32)
        phase             = _fread(fid, 1, :c)
        screen            = _fread(fid, 16, :c)
        calmode           = _fread(fid, 1, :i16)
        calmethod         = _fread(fid, 1, :i16)
        calupdate         = _fread(fid, 1, :i16)
        calbaseline       = _fread(fid, 1, :i16)
        calsweeps         = _fread(fid, 1, :i16)
        calattenuator     = _fread(fid, 1, :f32)
        calpulsevolt      = _fread(fid, 1, :f32)
        calpulsestart     = _fread(fid, 1, :f32)
        calpulsestop      = _fread(fid, 1, :f32)
        calfreq           = _fread(fid, 1, :f32)
        taskfile          = _fread(fid, 34, :c)
        seqfile           = _fread(fid, 34, :c)
        spectmethod       = _fread(fid, 1, :c)
        spectscaling      = _fread(fid, 1, :c)
        spectwindow       = _fread(fid, 1, :c)
        spectwinlength    = _fread(fid, 1, :f32)
        spectorder        = _fread(fid, 1, :c)
        notchfilter       = _fread(fid, 1, :c)
        headgain          = _fread(fid, 1, :i16)
        additionalfiles   = _fread(fid, 1, :i)
        unused            = _fread(fid, 5, :c)
        fspstopmethod     = _fread(fid, 1, :i16)
        fspstopmode       = _fread(fid, 1, :i16)
        fspfvalue         = _fread(fid, 1, :f32)
        fsppoint          = _fread(fid, 1, :i16)
        fspblocksize      = _fread(fid, 1, :i16)
        fspp1             = _fread(fid, 1, :ui16)
        fspp2             = _fread(fid, 1, :ui16)
        fspalpha          = _fread(fid, 1, :f32)
        fspnoise          = _fread(fid, 1, :f32)
        fspv1             = _fread(fid, 1, :i16)
        montage           = _fread(fid, 40, :c)
        eventfile         = _fread(fid, 40, :c)
        fratio            = _fread(fid, 1, :f32)
        minor_rev         = _fread(fid, 1, :c)
        eegupdate         = _fread(fid, 1, :i16)
        compressed        = _fread(fid, 1, :c)
        xscale            = _fread(fid, 1, :f32)
        yscale            = _fread(fid, 1, :f32)
        xsize             = _fread(fid, 1, :f32)
        ysize             = _fread(fid, 1, :f32)
        acmode            = _fread(fid, 1, :c)
        commonchnl        = _fread(fid, 1, :c)
        xtics             = _fread(fid, 1, :c)
        xrange            = _fread(fid, 1, :c)
        ytics             = _fread(fid, 1, :c)
        yrange            = _fread(fid, 1, :c)
        xscalevalue       = _fread(fid, 1, :f32)
        xscaleinterval    = _fread(fid, 1, :f32)
        yscalevalue       = _fread(fid, 1, :f32)
        yscaleinterval    = _fread(fid, 1, :f32)
        scaletoolx1       = _fread(fid, 1, :f32)
        scaletooly1       = _fread(fid, 1, :f32)
        scaletoolx2       = _fread(fid, 1, :f32)
        scaletooly2       = _fread(fid, 1, :f32)
        port              = _fread(fid, 1, :i16)
        numsamples        = _fread(fid, 1, :ui32)
        filterflag        = _fread(fid, 1, :c)
        lowcutoff         = _fread(fid, 1, :f32)
        lowpoles          = _fread(fid, 1, :i16)
        highcutoff        = _fread(fid, 1, :f32)
        highpoles         = _fread(fid, 1, :i16)
        filtertype        = _fread(fid, 1, :c)
        filterdomain      = _fread(fid, 1, :c)
        snrflag           = _fread(fid, 1, :c)
        coherenceflag     = _fread(fid, 1, :c)
        continuoustype    = _fread(fid, 1, :c)
        eventtablepos     = _fread(fid, 1, :ui32) # low-order bits of event table offset
        continuousseconds = _fread(fid, 1, :f32)
        channeloffset     = _fread(fid, 1, :i32)
        autocorrectflag   = _fread(fid, 1, :c)
        dcthreshold       = _fread(fid, 1, :ui8)

        # ------------------------------------------------------------ #
        # oer-channel header (75 bytes × ch_n)                         #
        # ------------------------------------------------------------ #
        clabels = String[]
        creference = Int64[]
        cskip = Int64[]
        creject = Int64[]
        cdisplay = Int64[]
        cbad = Int64[]
        cn = Int64[]
        cavg_reference = Int64[]
        cclipadd = Int64[]
        cx_coord = Float64[]
        cy_coord = Float64[]
        cveog_wt = Float64[]
        cveog_std = Float64[]
        csnr = Float64[]
        cheog_wt = Float64[]
        cheog_std = Float64[]
        cbaseline = Int64[]
        cfiltered = Int64[]
        cfsp = Int64[]
        caux1_wt = Float64[]
        caux1_std = Float64[]
        csenstivity = Float64[]
        cgain = Int64[]
        chipass = Int64[]
        clopass = Int64[]
        cpage = Int64[]
        csize = Int64[]
        cimpedance = Int64[]
        cphysicalchnl = Int64[]
        crectify = Int64[]
        ccalib = Float64[]

        for _ in 1:ch_n
            push!(clabels,        _vint2str(_fread(fid, 10, :c)))
            push!(creference,     _fread(fid, 1, :c))
            push!(cskip,          _fread(fid, 1, :c))
            push!(creject,        _fread(fid, 1, :c))
            push!(cdisplay,       _fread(fid, 1, :c))
            push!(cbad,           _fread(fid, 1, :c))
            push!(cn,             _fread(fid, 1, :ui16))
            push!(cavg_reference, _fread(fid, 1, :c))
            push!(cclipadd,       _fread(fid, 1, :c))
            push!(cx_coord,       _fread(fid, 1, :f32))
            push!(cy_coord,       _fread(fid, 1, :f32))
            push!(cveog_wt,       _fread(fid, 1, :f32))
            push!(cveog_std,      _fread(fid, 1, :f32))
            push!(csnr,           _fread(fid, 1, :f32))
            push!(cheog_wt,       _fread(fid, 1, :f32))
            push!(cheog_std,      _fread(fid, 1, :f32))
            push!(cbaseline,      _fread(fid, 1, :i16))
            push!(cfiltered,      _fread(fid, 1, :c))
            push!(cfsp,           _fread(fid, 1, :c))
            push!(caux1_wt,       _fread(fid, 1, :f32))
            push!(caux1_std,      _fread(fid, 1, :f32))
            push!(csenstivity,    _fread(fid, 1, :f32))
            push!(cgain,          _fread(fid, 1, :c))
            push!(chipass,        _fread(fid, 1, :c))
            push!(clopass,        _fread(fid, 1, :c))
            push!(cpage,          _fread(fid, 1, :ui8))
            push!(csize,          _fread(fid, 1, :ui8))
            push!(cimpedance,     _fread(fid, 1, :ui8))
            push!(cphysicalchnl,  _fread(fid, 1, :ui8))
            push!(crectify,       _fread(fid, 1, :c))
            push!(ccalib,         _fread(fid, 1, :f32))
        end

        # ------------------------------------------------------------ #
        # signal data                                                  #
        # ------------------------------------------------------------ #
        begdata = position(fid)
        # byte offset where signal data ends
        enddata = eventtablepos

        bytes_per_samp = data_format === :i32 ? 4 : 2
        nums = floor(Int64, (enddata - begdata) / ch_n / bytes_per_samp)

        seek(fid, begdata)
        # CNT stores samples in multiplexed order: ch1_s1, ch2_s1, ..., chN_s1, ch1_s2, ...
        # read into a flat buffer then reshape - avoids millions of individual _fread calls
        raw = [_fread(fid, 1, data_format) for _ in 1:(nums * ch_n)]
        data = zeros(ch_n, nums)
        @inbounds for s in 1:nums, ch in 1:ch_n
            data[ch, s] = raw[(s - 1) * ch_n + ch]
        end

        # scale to μV: value_μV = (raw - baseline) × sensitivity × (calibration / 204.8)
        @inbounds for idx in 1:ch_n
            mf = csenstivity[idx] * (ccalib[idx] / 204.8)
            data[idx, :] = @views (data[idx, :] .- cbaseline[idx]) .* mf
        end
        data = reshape(data, ch_n, nums, 1)

        # ------------------------------------------------------------ #
        # event table                                                  #
        # event table offset = (prevfile << 32) | eventtablepos        #
        # ------------------------------------------------------------ #
        # `prevfile` holds the high-order 32 bits; `eventtablepos` the low
        et_offset = (Int64(prevfile) * Int64(2)^32) + Int64(eventtablepos)
        seek(fid, et_offset)

        teeg = _fread(fid, 1, :ui8)
        et_fsize = _fread(fid, 1, :ui32)
        _et_reserved  = _fread(fid, 1, :ui32)

        sizeEvent1 = 8
        sizeEvent2 = 19
        sizeEvent3 = 19
        evt_stimtype = Int64[]
        evt_offset = Int64[]
        evt_type = Int64[]
        evt_code = Int64[]
        nevents = 0

        if teeg == 2
            nevents = Int(et_fsize ÷ sizeEvent2)
            for _ in 1:nevents
                push!(evt_stimtype, _fread(fid, 1, :ui16))
                _fread(fid, 1, :c); _fread(fid, 1, :ui8)
                push!(evt_offset, _fread(fid, 1, :i32))
                push!(evt_type, _fread(fid, 1, :i16))
                push!(evt_code, _fread(fid, 1, :i16))
                _fread(fid, 1, :f32)
                _fread(fid, 1, :c); _fread(fid, 1, :c); _fread(fid, 1, :c)
            end

        elseif teeg == 3
            # type 3: offset field encodes the global sample frame directly
            nevents = Int(et_fsize ÷ sizeEvent3)
            for _ in 1:nevents
                push!(evt_stimtype, _fread(fid, 1, :ui16))
                _fread(fid, 1, :c); _fread(fid, 1, :ui8)
                os = _fread(fid, 1, :ui32)
                push!(evt_offset, os * bytes_per_samp * ch_n)
                push!(evt_type, _fread(fid, 1, :i16))
                push!(evt_code, _fread(fid, 1, :i16))
                _fread(fid, 1, :f32)
                _fread(fid, 1, :c); _fread(fid, 1, :c); _fread(fid, 1, :c)
            end

        elseif teeg == 1
            nevents = Int(et_fsize ÷ sizeEvent1)
            for _ in 1:nevents
                push!(evt_stimtype, _fread(fid, 1, :ui16))
                _fread(fid, 1, :c); _fread(fid, 1, :ui8)
                push!(evt_offset,   _fread(fid, 1, :i32))
                # event1 has no type/code fields; vectors stay empty
            end
        else
            _warn("Unknown event table type (TEEG=$teeg); events will not be imported.")
        end

        # convert byte offsets → sample indices
        # initial_offset = 900-byte global header + 75-byte-per-channel header
        if nevents > 0
            initial_offset = 900 + (ch_n * 75)
            @inbounds for idx in 1:nevents
                evt_offset[idx] = if teeg == 3
                    evt_offset[idx] ÷ (bytes_per_samp * ch_n)
                else
                    (evt_offset[idx] - initial_offset) ÷ (bytes_per_samp * ch_n)
                end
            end
        end

        (; patient, recording_date, recording_time, ch_n, sampling_rate,
           clabels, cbaseline, csenstivity, ccalib, cx_coord, cy_coord,
           begdata, eventtablepos, prevfile, data_format, data,
           nevents, evt_stimtype, evt_offset, evt_type, evt_code)
    end
    # file closed here in all cases, including exceptions

    # unpack named tuple
    (; patient, recording, recording_date, recording_time, file_type, ch_n, clabels,
       transducers, units, prefiltering, ch_type, annotation_channels, sampling_rate,
       gain, annotations, data) = imported_object
    # reuse binding

    # ------------------------------------------------------------------ #
    # channel types and units                                            #
    # ------------------------------------------------------------------ #
    ch_type = detect_type ? _set_channel_types(clabels, "eeg") : repeat(["eeg"], ch_n)
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    # ------------------------------------------------------------------ #
    # build markers DataFrame from parsed event table                    #
    # ------------------------------------------------------------------ #
    markers = if nevents == 0
        DataFrame(
            :id => String[],
            :start => Float64[],
            :length => Float64[],
            :value => String[],
            :channel => Int64[])
    else
        DataFrame(
            :id => [isempty(evt_stimtype) ? "mrk" : string(evt_stimtype[i]) for i in 1:nevents],
            :start => evt_offset ./ sampling_rate,
            :length => fill(0.0, nevents),
            :value => [isempty(evt_code) ? "" : string(evt_code[i]) for i in 1:nevents],
            :channel => fill(0, nevents),
        )
    end

    # ------------------------------------------------------------------ #
    # time axes                                                          #
    # ------------------------------------------------------------------ #
    nums = size(data, 2)
    # FIX: original had identical expressions for time_pts and epoch_time;
    # they are only equal because the data has a single epoch. Unified.
    time_pts   = round.(range(0; step = 1/sampling_rate, length = nums); digits = 4)
    epoch_time = time_pts   # single epoch → identical

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                              #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(
        id = "", first_name = "", middle_name = "",
        last_name          = patient,
        head_circumference = -1,
        handedness         = "",
        weight             = -1,
        height             = -1)
    r = _create_recording_eeg(
        data_type = "eeg",
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = "CNT",
        recording = "",
        recording_date = recording_date,
        recording_time = replace(recording_time, '.' => ':'),
        recording_notes = "",
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = _detect_montage(clabels, ch_type, "eeg"),
        clabels = clabels,
        transducers = repeat([""], ch_n),
        units = units,
        prefiltering = repeat([""], ch_n),
        line_frequency = 50, # TODO: make this a keyword argument
        sampling_rate = sampling_rate,
        gain = ones(ch_n),
        bad_channels = zeros(Bool, size(data, 1)))
    e   = _create_experiment(name = "", notes = "", design = "")
    hdr = _create_header(subject = s, recording = r, experiment = e)

    locs = _initialize_locs()
    obj  = NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, epoch_time, data)
    _initialize_locs!(obj)

    _info("Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end]; digits=2)) s)")

    return obj

end