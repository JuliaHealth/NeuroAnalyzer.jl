export import_cnt

"""
    import_cnt(file_name; <keyword arguments>)

Load Neuroscan continuous signal file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label
- `data_format::Symbol=:i32`: Neuroscan stores data either in 16-bit (`:i16`) or 32-bit (`:i32`) representation, but does not say so in the file format

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

Based on loadcnt.m by Sean Fitzgibbon and Arnaud Delorme (https://cnl.salk.edu/~arno/cntload/index.html)
"""
function import_cnt(file_name::String; data_format::Symbol=:i32, detect_type::Bool=true)

    _check_var(data_format, [:i32, :i16], "data_format")

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) == ".cnt" "This is not CNT file."

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        @error "File $file_name cannot be loaded."
    end

    file_type = "CNT"

    rev = _vint2str(_fread(fid, 12, :c))
    rev != "Version 3.0" && _warn("CNT files with version number other than 3.0 may not be imported properly; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
    nextfile = _fread(fid, 1, :l)
    prevfile = _fread(fid, 1, :ul)
    type = _fread(fid, 1, :c)
    id = _vint2str(_fread(fid, 20, :c))
    oper = _vint2str(_fread(fid, 20, :c))
    doctor = _vint2str(_fread(fid, 20, :c))
    referral = _vint2str(_fread(fid, 20, :c))
    hospital = _vint2str(_fread(fid, 20, :c))
    patient = _vint2str(_fread(fid, 20, :c))
    age = _fread(fid, 1, :i16)
    sex = _fread(fid, 1, :c)
    hand = _fread(fid, 1, :c)
    med = _vint2str(_fread(fid, 20 , :c))
    category = _vint2str(_fread(fid, 20 , :c))
    state = _vint2str(_fread(fid, 20 , :c))
    label = _vint2str(_fread(fid, 20 , :c))
    recording_date = _vint2str(_fread(fid, 10 , :c))
    recording_time = _vint2str(_fread(fid, 12 , :c))
    mean_age = _fread(fid, 1, :f32)
    stdev = _fread(fid, 1, :f32)
    n = _fread(fid, 1, :i16)
    compfile = _fread(fid, 38 ,:c)
    spectwincomp = _fread(fid, 1, :f32)
    meanaccuracy = _fread(fid, 1, :f32)
    meanlatency = _fread(fid, 1, :f32)
    sortfile = _fread(fid, 46 ,:c)
    numevents = _fread(fid, 1, :i)
    compoper = _fread(fid, 1, :c)
    avgmode = _fread(fid, 1, :c)
    review = _fread(fid, 1, :c)
    nsweeps = _fread(fid, 1, :ui16)
    compsweeps = _fread(fid, 1, :ui16)
    acceptcnt = _fread(fid, 1, :ui16)
    rejectcnt = _fread(fid, 1, :ui16)
    pnts = _fread(fid, 1, :ui16)
    ch_n = _fread(fid, 1, :ui16)
    avgupdate = _fread(fid, 1, :ui16)
    domain = _fread(fid, 1, :c)
    variance = _fread(fid, 1, :c)
    sampling_rate = _fread(fid, 1, :ui16)
    scale = _fread(fid, 1, :f64)
    veogcorrect = _fread(fid, 1, :c)
    heogcorrect = _fread(fid, 1, :c)
    aux1correct = _fread(fid, 1, :c)
    aux2correct = _fread(fid, 1, :c)
    veogtrig = _fread(fid, 1, :f32)
    heogtrig = _fread(fid, 1, :f32)
    aux1trig = _fread(fid, 1, :f32)
    aux2trig = _fread(fid, 1, :f32)
    heogchnl = _fread(fid, 1, :i16)
    veogchnl = _fread(fid, 1, :i16)
    aux1chnl = _fread(fid, 1, :i16)
    aux2chnl = _fread(fid, 1, :i16)
    veogdir = _fread(fid, 1, :c)
    heogdir = _fread(fid, 1, :c)
    aux1dir = _fread(fid, 1, :c)
    aux2dir = _fread(fid, 1, :c)
    veog_n = _fread(fid, 1, :i16)
    heog_n = _fread(fid, 1, :i16)
    aux1_n = _fread(fid, 1, :i16)
    aux2_n = _fread(fid, 1, :i16)
    veogmaxcnt = _fread(fid, 1, :i16)
    heogmaxcnt = _fread(fid, 1, :i16)
    aux1maxcnt = _fread(fid, 1, :i16)
    aux2maxcnt = _fread(fid, 1, :i16)
    veogmethod = _fread(fid, 1, :c)
    heogmethod = _fread(fid, 1, :c)
    aux1method = _fread(fid, 1, :c)
    aux2method = _fread(fid, 1, :c)
    ampsensitivity = _fread(fid, 1, :f32)
    lowpass = _fread(fid, 1, :c)
    highpass = _fread(fid, 1, :c)
    notch = _fread(fid, 1, :c)
    autoclipadd = _fread(fid, 1, :c)
    baseline = _fread(fid, 1, :c)
    offstart = _fread(fid, 1, :f32)
    offstop = _fread(fid, 1, :f32)
    reject = _fread(fid, 1, :c)
    rejstart = _fread(fid, 1, :f32)
    rejstop = _fread(fid, 1, :f32)
    rejmin = _fread(fid, 1, :f32)
    rejmax = _fread(fid, 1, :f32)
    trigtype = _fread(fid, 1, :c)
    trigval = _fread(fid, 1, :f32)
    trigchnl = _fread(fid, 1, :c)
    trigmask = _fread(fid, 1, :i16)
    trigisi = _fread(fid, 1, :f32)
    trigmin = _fread(fid, 1, :f32)
    trigmax = _fread(fid, 1, :f32)
    trigdir = _fread(fid, 1, :c)
    autoscale = _fread(fid, 1, :c)
    n2 = _fread(fid, 1, :i16)
    dir = _fread(fid, 1, :c)
    dispmin = _fread(fid, 1, :f32)
    dispmax = _fread(fid, 1, :f32)
    xmin = _fread(fid, 1, :f32)
    xmax = _fread(fid, 1, :f32)
    automin = _fread(fid, 1, :f32)
    automax = _fread(fid, 1, :f32)
    zmin = _fread(fid, 1, :f32)
    zmax = _fread(fid, 1, :f32)
    lowcut = _fread(fid, 1, :f32)
    highcut = _fread(fid, 1, :f32)
    common = _fread(fid, 1, :c)
    savemode = _fread(fid, 1, :c)
    manmode = _fread(fid, 1, :c)
    ref = _fread(fid, 10 ,:c)
    rectify = _fread(fid, 1, :c)
    displayxmin = _fread(fid, 1, :f32)
    displayxmax = _fread(fid, 1, :f32)
    phase = _fread(fid, 1, :c)
    screen = _fread(fid, 16 ,:c)
    calmode = _fread(fid, 1, :i16)
    calmethod = _fread(fid, 1, :i16)
    calupdate = _fread(fid, 1, :i16)
    calbaseline = _fread(fid, 1, :i16)
    calsweeps = _fread(fid, 1, :i16)
    calattenuator = _fread(fid, 1, :f32)
    calpulsevolt = _fread(fid, 1, :f32)
    calpulsestart = _fread(fid, 1, :f32)
    calpulsestop = _fread(fid, 1, :f32)
    calfreq = _fread(fid, 1, :f32)
    taskfile = _fread(fid, 34 ,:c)
    seqfile = _fread(fid, 34 ,:c)
    spectmethod = _fread(fid, 1, :c)
    spectscaling = _fread(fid, 1, :c)
    spectwindow = _fread(fid, 1, :c)
    spectwinlength = _fread(fid, 1, :f32)
    spectorder = _fread(fid, 1, :c)
    notchfilter = _fread(fid, 1, :c)
    headgain = _fread(fid, 1, :i16)
    additionalfiles = _fread(fid, 1, :i)
    unused = _fread(fid, 5, :c)
    fspstopmethod = _fread(fid, 1, :i16)
    fspstopmode = _fread(fid, 1, :i16)
    fspfvalue = _fread(fid, 1, :f32)
    fsppoint = _fread(fid, 1, :i16)
    fspblocksize = _fread(fid, 1, :i16)
    fspp1 = _fread(fid, 1, :ui16)
    fspp2 = _fread(fid, 1, :ui16)
    fspalpha = _fread(fid, 1, :f32)
    fspnoise = _fread(fid, 1, :f32)
    fspv1 = _fread(fid, 1, :i16)
    montage = _fread(fid, 40 ,:c)
    eventfile = _fread(fid, 40 ,:c)
    fratio = _fread(fid, 1, :f32)
    minor_rev = _fread(fid, 1, :c)
    eegupdate = _fread(fid, 1, :i16)
    compressed = _fread(fid, 1, :c)
    xscale = _fread(fid, 1, :f32)
    yscale = _fread(fid, 1, :f32)
    xsize = _fread(fid, 1, :f32)
    ysize = _fread(fid, 1, :f32)
    acmode = _fread(fid, 1, :c)
    commonchnl = _fread(fid, 1, :c)
    xtics = _fread(fid, 1, :c)
    xrange = _fread(fid, 1, :c)
    ytics = _fread(fid, 1, :c)
    yrange = _fread(fid, 1, :c)
    xscalevalue = _fread(fid, 1, :f32)
    xscaleinterval = _fread(fid, 1, :f32)
    yscalevalue = _fread(fid, 1, :f32)
    yscaleinterval = _fread(fid, 1, :f32)
    scaletoolx1 = _fread(fid, 1, :f32)
    scaletooly1 = _fread(fid, 1, :f32)
    scaletoolx2 = _fread(fid, 1, :f32)
    scaletooly2 = _fread(fid, 1, :f32)
    port = _fread(fid, 1, :i16)
    numsamples = _fread(fid, 1, :ui32)
    filterflag = _fread(fid, 1, :c)
    lowcutoff = _fread(fid, 1, :f32)
    lowpoles = _fread(fid, 1, :i16)
    highcutoff = _fread(fid, 1, :f32)
    highpoles = _fread(fid, 1, :i16)
    filtertype = _fread(fid, 1, :c)
    filterdomain = _fread(fid, 1, :c)
    snrflag = _fread(fid, 1, :c)
    coherenceflag = _fread(fid, 1, :c)
    continuoustype = _fread(fid, 1, :c)
    eventtablepos = _fread(fid, 1, :ui32)
    continuousseconds = _fread(fid, 1, :f32)
    channeloffset = _fread(fid, 1, :i32)
    autocorrectflag = _fread(fid, 1, :c)
    dcthreshold = _fread(fid, 1, :ui8)

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

    for idx in 1:ch_n
        push!(clabels, _vint2str(_fread(fid, 10, :c)))
        push!(creference, _fread(fid, 1, :c))
        push!(cskip, _fread(fid, 1, :c))
        push!(creject, _fread(fid, 1, :c))
        push!(cdisplay, _fread(fid, 1, :c))
        push!(cbad, _fread(fid, 1, :c))
        push!(cn, _fread(fid, 1, :ui16))
        push!(cavg_reference, _fread(fid, 1, :c))
        push!(cclipadd, _fread(fid, 1, :c))
        push!(cx_coord, _fread(fid, 1, :f32))
        push!(cy_coord, _fread(fid, 1, :f32))
        push!(cveog_wt, _fread(fid, 1, :f32))
        push!(cveog_std, _fread(fid, 1, :f32))
        push!(csnr, _fread(fid, 1, :f32))
        push!(cheog_wt, _fread(fid, 1, :f32))
        push!(cheog_std, _fread(fid, 1, :f32))
        push!(cbaseline, _fread(fid, 1, :i16))
        push!(cfiltered, _fread(fid, 1, :c))
        push!(cfsp, _fread(fid, 1, :c))
        push!(caux1_wt, _fread(fid, 1, :f32))
        push!(caux1_std, _fread(fid, 1, :f32))
        push!(csenstivity, _fread(fid, 1, :f32))
        push!(cgain, _fread(fid, 1, :c))
        push!(chipass, _fread(fid, 1, :c))
        push!(clopass, _fread(fid, 1, :c))
        push!(cpage, _fread(fid, 1, :ui8))
        push!(csize, _fread(fid, 1, :ui8))
        push!(cimpedance, _fread(fid, 1, :ui8))
        push!(cphysicalchnl, _fread(fid, 1, :ui8))
        push!(crectify, _fread(fid, 1, :c))
        push!(ccalib, _fread(fid, 1, :f32))
    end

    begdata = position(fid)
    enddata = eventtablepos

    seek(fid, begdata)

    if data_format === :i16
        nums = floor(Int64, (enddata - begdata) / ch_n / 2)
    else
        nums = floor(Int64, (enddata - begdata) / ch_n / 4)
    end

    data = zeros(ch_n, nums)
    for idx1 in 1:nums, idx2 in 1:ch_n
        data[idx2, idx1] = _fread(fid, 1, data_format)
    end

    # scaling to μV
    for idx in 1:ch_n
        bas = cbaseline[idx]
        sen = csenstivity[idx]
        cal = ccalib[idx]
        mf = sen * (cal / 204.8)
        data[idx, :] = @views (data[idx, :] .- bas) .* mf
    end
    data = reshape(data, size(data, 1), size(data, 2), 1)

    # events table

    # prevfile contains high order bits of event table offset
    # eventtablepos contains the low order bits
    et_offset = (prevfile * (2^32)) + eventtablepos
    seek(fid, et_offset)

    teeg = _fread(fid, 1, :ui8)
    fsize = _fread(fid, 1, :ui32)
    offset = _fread(fid, 1, :ui32)

    sizeEvent1 = 8  # 8  bytes for Event1
    sizeEvent2 = 19 # 19 bytes for Event2
    sizeEvent3 = 19 # 19 bytes for Event3

    if teeg == 2
        nevents = fsize ÷ sizeEvent2
        stimtype = Int64[]
        offset = Int64[]
        type = Int64[]
        code = Int64[]
        if nevents > 0
            for i in 1:nevents
                push!(stimtype, _fread(fid, 1, :ui16))
               _fread(fid, 1, :c)
               _fread(fid, 1, :ui8)
                push!(offset, _fread(fid, 1, :i32))
                push!(type, _fread(fid, 1, :i16))
                push!(code, _fread(fid, 1, :i16))
               _fread(fid, 1, :f32)
               _fread(fid, 1, :c)
               _fread(fid, 1, :c)
               _fread(fid, 1, :c)
            end
        end
    elseif teeg == 3 # type 3 is similar to type 2 except the offset field encodes the global sample frame
        nevents = fsize ÷ sizeEvent3
        stimtype = Int64[]
        offset = Int64[]
        type = Int64[]
        code = Int64[]
        if nevents > 0
            bytes_per_samp = data_format === :i32 ? 4 : 2
            for i=1:nevents
                push!(stimtype, _fread(fid, 1, :ui16))
               _fread(fid, 1, :c)
               _fread(fid, 1, :ui8)
                os = _fread(fid, 1, :ui32)
                push!(offset, os * bytes_per_samp * ch_n)
                push!(type, _fread(fid, 1, :i16))
                push!(code, _fread(fid, 1, :i16))
               _fread(fid, 1, :f32)
               _fread(fid, 1, :c)
               _fread(fid, 1, :c)
               _fread(fid, 1, :c)
            end
        end
    elseif teeg == 1
        nevents = fsize ÷ sizeEvent1
        stimtype = Int64[]
        offset = Int64[]
        type = Int64[]
        code = Int64[]
        if nevents > 0
            for i in 1:nevents
                push!(stimtype, _fread(fid, 1, :ui16))
               _fread(fid, 1, :c)
               _fread(fid, 1, :ui8)
                push!(offset, _fread(fid, 1, :i32))
            end
        end
    end

    if nevents > 0
        # initial offset : header + electrodes desc
        initial_offset = 900 + (ch_n * 75)
        if data_format === :i16
            if teeg == 3
                for idx in 1:nevents
                    offset[idx] = offset[idx] ÷ (2 * ch_n)
                end
            else
                for idx in 1:nevents
                    offset[idx] = (offset[idx] - initial_offset) ÷ (2 * ch_n)
                end
            end
        else
            if teeg == 3
                for idx in 1:nevents
                    offset[idx] = offset[idx] ÷ (4 * ch_n)
                end
            else
                for idx in 1:nevents
                    offset[idx] = (offset[idx] - initial_offset) ÷ (4 * ch_n)
                end
            end
        end
    end

    close(fid)

    if detect_type
        ch_type = _set_channel_types(clabels, "eeg")
    else
        ch_type = repeat(["eeg"], ch_n)
    end
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :description=>String[],
                        :channel=>Int64[])

    time_pts = round.(collect(0:1/sampling_rate:nums / sampling_rate)[1:end-1], digits=3)
    ep_time = round.((collect(0:1/sampling_rate:nums / sampling_rate))[1:end-1], digits=3)

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    data_type = "eeg"

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name=patient,
                        head_circumference=-1,
                        handedness="",
                        weight=-1,
                        height=-1)

    r = _create_recording_eeg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type=file_type,
                              recording="",
                              recording_date=recording_date,
                              recording_time=recording_time,
                              recording_notes="",
                              channel_type=ch_type,
                              channel_order=_sort_channels(ch_type),
                              reference=_detect_montage(clabels, ch_type, data_type),
                              clabels=clabels,
                              transducers=repeat([""], ch_n),
                              units=units,
                              prefiltering=repeat([""], ch_n),
                              sampling_rate=sampling_rate,
                              gain=ones(ch_n),
                              bad_channels=zeros(Bool, size(data, 1), 1))
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data, components, markers, locs, history)
    _initialize_locs!(obj)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj

end
