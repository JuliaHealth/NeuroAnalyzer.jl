module EDFPlus
using DSP
using Dates
using IterTools
using Core.Intrinsics

if VERSION < v"0.7.0"
    @warn("This version of EDFPlus requires Julia 0.7 or above.")
end

export ChannelParam, BEDFPlus, Annotation, DataFormat, FileStatus, version,
       loadfile, writefile!, closefile!, samplerate, addannotation!,
       epoch_iterator, annotation_epoch_iterator, readBiosemiStatus,
       digitalchanneldata, physicalchanneldata,
       channeltimesegment, multichanneltimesegment,
       highpassfilter, lowpassfilter, notchfilter, trim

#
# Except for some of the data structures and constants, most of this code is a
# complete rewrite of the C version. Other C ports will need different API calls.
#==============================================================================
* The original EDFlib C code was
* Copyright (c) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017 Teunis van Beelen
* All rights reserved.
*
* email: teuniz@gmail.com
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY Teunis van Beelen ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Teunis van Beelen BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES
* LOSS OF USE, DATA, OR PROFITS OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* For more info about the EDF and EDF+ format, visit: http://edfplus.info/specs/
* For more info about the BDF and BDF+ format, visit: http://www.teuniz.net/edfbrowser/bdfplus%20format%20description.html
=================================================================================#
#
# See also: https://www.edfplus.info/specs/edffaq.html and
#           https://www.biosemi.com/faq/file_format.htm
#

const EDFPLUS_VERSION = 0.59
const MAX_CHANNELS =          512
const MAX_ANNOTATION_LENGTH = 512

"""
    DataFormat
enum for types this package handles. Current format for a potential translation is also /same/.
"""
@enum DataFormat bdf bdfplus edf edfplus same

"""
    FileStatus
enum for type or state of file: type of data detected, whether any errors
"""
@enum FileStatus EDF EDFPLUS BDF BDFPLUS READ_ERROR FORMAT_ERROR CLOSED

"""
    type Int24

# 24-bit integer routines for BDF format signal data.
BDF and BDF+ files use 24 bits per data signal point.
The module caches these after reading as Int32 to fit LLVM CPU registers.
"""
primitive type Int24 24 end
Int24(x::Int) = Core.Intrinsics.trunc_int(Int24, x)
Int(x::Int24) = Core.Intrinsics.zext_int(Int, x)
function writei24(stream::IO, x)
    b1::UInt8 = x & 0xff
    b2::UInt8 = (x >> 8) & 0xff
    b3::UInt8 = (x >> 16) & 0xff
    write(stream, [b1, b2, b3])
end

""" static function to state version of module """
version() = EDFPLUS_VERSION

"""
    mutable struct ChannelParam
Parameters for each channel in the EEG record.
"""
mutable struct ChannelParam      # this structure contains all the relevant EDF-signal parameters of one signal
  label::String                  # label (name) of the signal, eg "C4" if in 10-20 labeling terms
  transducer::String             # signal transducer type
  physdimension::String          # physical dimension (uV, bpm, mA, etc.)
  physmax::Float64               # physical maximum, usually the maximum input of the ADC
  physmin::Float64               # physical minimum, usually the minimum input of the ADC
  digmax::Int                    # digital maximum, usually the maximum output of the ADC, cannot be higher than 32767 for EDF or 8388607 for BDF
  digmin::Int                    # digital minimum, usually the minimum output of the ADC, cannot be lower than -32768 for EDF or -8388608 for BDF
  smp_per_record::Int            # number of samples of this signal in a datarecord
  prefilter::String              # channel prefiltering settings if any
  reserved::String               # header reserved ascii text, 32 bytes
  offset::Float64                # offset of center of physical data value from center of digital values
  bufoffset::Int                 # bytes from start of record to start of this channel (zero for first channel)
  bitvalue::Float64              # physical data value of one unit change in digital value
  annotation::Bool               # true if is an annotation not a binary mapped signal data channel
  ChannelParam() = new("","","",0.0,0.0,0,0,0,"","",0.0,0,0.0,false)
end

"""
    mutable struct Annotation

These are text strings denoting a time, optionally duration, and a list of notes
about the signal at that particular time in the recording. The first onset time
of the annotation channel gives a fractional second offset adjustment of the
start time of that record, which is specified in whole seconds in the header.
"""
mutable struct Annotation
    onset::Float64
    duration::String
    annotation::Array{String,1}
    Annotation() = new(0.0,"",[])
    Annotation(o,d,arr) = new(o, d, typeof(arr) == String ? [arr] : arr)
end
# max size of annotationtext
const MAX_ANNOTATION_TEXT_LENGTH = 40
# minimum length per record of annotation channnel (for writing a new file)
const MIN_ANNOTATION_CHANNEL_LENGTH = 120

"""
    mutable struct BEDFPlus
Data struct for EDF, EDF+, BDF, and BDF+ EEG type signal files.
"""
mutable struct BEDFPlus                   # signal file data for EDF, BDF, EDF+, and BDF+ files
    ios::IOStream                         # file handle for the file containing the data
    path::String                          # file pathname
    writemode::Bool                       # true if is intended for writing to file
    version::String                       # version of the file format
    edf::Bool                             # EDF?
    edfplus::Bool                         # EDF+?
    bdf::Bool                             # BDF?
    bdfplus::Bool                         # BDF+?
    discontinuous::Bool                   # discontinuous (EDF+D?)
    filetype::FileStatus                  # @enum FileStatus as above
    channelcount::Int                     # total number of EDF signal bands in the file INCLUDING annotation channels
    file_duration::Float64                # duration of the file in seconds expressed as 64-bit floating point
    startdate_day::Int                    # startdate of study, day of month of startdate of study
    startdate_month::Int                  # startdate month
    startdate_year::Int                   # startdate year
    starttime_subsecond::Float64          # starttime offset in seconds, should be < 1 sec in size. Only used by EDFplus and BDFplus
    starttime_second::Int                 # this is in integer seconds, the field above makes it more precise
    starttime_minute::Int                 # startdate and time, minutes
    starttime_hour::Int                   # 0 to 23, midnight is 00:00:00
    # next 11 fields are for EDF+ and BDF+ files only
    patient::String                       # contains patientfield of header, is always empty when filetype is EDFPLUS or BDFPLUS
    recording::String                     # contains recordingfield of header, is always empty when filetype is EDFPLUS or BDFPLUS
    patientcode::String                   # empty when filetype is EDF or BDF
    gender::String                        # empty when filetype is EDF or BDF
    birthdate::String                     # empty when filetype is EDF or BDF
    patientname::String                   # empty when filetype is EDF or BDF
    patient_additional::String            # empty when filetype is EDF or BDF
    admincode::String                     # empty when filetype is EDF or BDF
    technician::String                    # empty when filetype is EDF or BDF
    equipment::String                     # empty when filetype is EDF or BDF
    recording_additional::String          # empty when filetype is EDF or BDF
    datarecord_duration::Float64          # duration of one datarecord in units of seconds
    datarecords::Int64                    # number of datarecords in the file
    startdatestring::String               # date recording started in dd-uuu-yyyy format
    reserved::String                      # reserved, 32 byte string
    headersize::Int                       # size of header in bytes
    recordsize::Int                       # size of one data record in bytes, these follow header
    annotationchannel::Int                # position in record of annotation channel
    mapped_signals::Array{Int,1}          # positions in record of channels carrying data
    signalparam::Array{ChannelParam,1}    # Array of structs which contain the per-signal parameters
    annotations::Array{Array{Annotation,1},1} # Array of lists of annotations
    EDFsignals::Array{Int16,2}    # 2D array, each row a record, columns are channels including annotations
    BDFsignals::Array{Int32,2}    # Note that either EDFsignals or BDFsignals is used
    BEDFPlus() = new(IOStream("nothing"),"",false,"",false,false,false,false,false,READ_ERROR,0,0.0,0,0,0,0.0,0,0,0,
                        "","","","","","","","","","","",0.0,0,"","",0,0,0,
                        Array{Int,1}(undef,0),Array{ChannelParam,1}(undef,0),
                        Array{Array{Annotation,1},1}(undef,0),Array{Int16,2}(undef,0,0),Array{Int32,2}(undef,0,0))
end

"""
    loadfile(path::String, read_annotations=true)

Load a BDF+ or EDF+ type file.
Takes a pathname. Will ignore annotations if second argument is set false.
Returns a BEDFPlus structure including header and data.
"""
function loadfile(path::String, read_annotations=true)
    edfh = BEDFPlus()
    fh = open(path, "r")
    edfh.path = path
    edfh.ios = fh
    checkfile!(edfh)
    if edfh.filetype == FORMAT_ERROR
        throw("Bad EDF/BDF file format at file $path")
    end
    edfh.writemode = false
    if edfh.edf
        edfh.filetype = EDF
    end
    if edfh.edfplus
        edfh.filetype = EDFPLUS
    end
    if edfh.bdf
        edfh.filetype = BDF
    end
    if edfh.bdfplus
        edfh.filetype = BDFPLUS
    end
    edfh.file_duration = edfh.datarecord_duration * edfh.datarecords
    if edfh.edfplus == false && edfh.bdfplus == false
        edfh.patientcode = ""
        edfh.gender = ""
        edfh.birthdate = ""
        edfh.patientname = ""
        edfh.patient_additional = ""
        edfh.admincode = ""
        edfh.technician = ""
        edfh.equipment = ""
        edfh.recording_additional = ""
    else
        # EDF+ and BDF+ use different ID information so blank other fields
        edfh.patient = ""
        edfh.recording = ""
        if read_annotations
            readannotations!(edfh)
        end
    end
    readdata!(edfh)
    edfh.path = path
    return edfh
end

"""
    writefile!(edfh, newpath; acquire=dummyacquire, sigformat=same)

Write to data in the edfh struct to the file indicated by newpath
Returns the file handle of the file written, opened for reading

NOTE: The header needs to be completely specified at function start except for
the final number of records, which will be updated after all data records
are written. For a system that is recording the data as it is written, the
acquire(edfh) function should write the data according the the header parameters.
NB: iff the function converts from BDF to EDF or EDF to BDF, the edfh struct is changed.
"""
function writefile!(edfh, newpath; acquire=dummyacquire, sigformat=same)
    if sigformat == same || sigformat == bdfplus || sigformat == edfplus
        if sigformat == bdfplus && edfh.edfplus
            translate16to24bits!(edfh)
            edfh.edfplus = false
            edfh.bdfplus = true
            edfh.edf = false
            edfh.bdf = false
            edfh.filetype = BDFPLUS
            boff = 0
            for chan in edfh.signalparam
                chan.digmin = -8388608
                chan.digmax = 8388607
                chan.bufoffset = boff
                boff += chan.smp_per_record * 3
            end
        elseif sigformat == edfplus && edfh.bdfplus
            translate24to16bits!(edfh)
            edfh.edfplus = true
            edfh.bdfplus = false
            edfh.edf = false
            edfh.bdf = false
            edfh.filetype = EDFPLUS
            boff = 0
            for chan in edfh.signalparam
                chan.digmin = -32768
                chan.digmax = 32767
                chan.bufoffset = boff
                boff += chan.smp_per_record * 2
            end
        end
        fh = open(newpath,"w+")
        written = writeheader(edfh, fh)
        acquirewritten = acquire(edfh)
        # if we did an acquire, the acquire function wrote the channel data
        if acquirewritten > 0
            written += acquirewritten
            seekstart(fh)
            seek(fh, 236)
            writeleftjust(fh, edfh.datarecords, 8)
        else # otherwise write what is there in memory for writing
            written += (edfh.edfplus || edfh.edf) ? writeEDFrecords!(edfh, fh) :
                                                    writeBDFrecords!(edfh, fh)
        end
        # close handle, reopen as a read handle, load/check file, return new handle
        close(fh)
        newedfh = loadfile(newpath, true)
        println("$newpath written successfully, $written bytes.")
        return newedfh
    elseif sigformat == edf
        @warn("Converting file to $newpath as an EDF compatible EDF+ file.")
        writefile!(edfh, newpath, acquire=acquire, sigformat=edfplus)
        return edfh
    elseif sigformat == bdf
        @warn("Converting file to $newpath as a BDF compatible BDF+ file.")
        writefile!(edfh, newpath, acquire=acquire, sigformat=bdfplus)
        return edfh
    else
        throw("Unknown signal file write format request: $sigformat")
    end
end

"""
    epoch_iterator(edfh, epochsecs; channels, startsec, endsec, physical)

Make an iterator for EEG epochs of a given duration between start and stop times.
# Required arguments
- edfh BEDFPlus struct
- epochsecs second duration of each epoch

# Optional arguments
- channels List of channel numbers for data, defaults to all signal channels
- startsec Starting position from 0 at start of file, defaults to file start
- endsec Ending position in seconds from start of _file_, defaults to file end
- physical Whether to return data as translated to the physical units, defaults to true
"""
function epoch_iterator(edfh, epochsecs; channels=edfh.mapped_signals,
                              startsec=0, endsec=edfh.file_duration, physical=true)
    epochs = collect(startsec:epochsecs:endsec)[1:end-1]
    epochwidth = epochs[2] - epochs[1]
    return imap(x -> multichanneltimesegment(edfh,channels,x,x+epochwidth, physical), epochs)
end

"""
    annotation_epoch_iterator(edfh, epochsecs; startsec, endsec)

Return an iterator for a group of annotations for a given epoch as in epoch_iterator
"""
function annotation_epoch_iterator(edfh, epochsecs; startsec=0, endsec=edfh.file_duration)
    epochs = collect(startsec:epochsecs:endsec)
    achan = edfh.annotationchannel
    if length(epochs) < 2
        markers = [signalat(edfh,startsec, achan), signalat(edfh,endsec, achan)]
    else
        markers = map(t->signalat(edfh,t, achan), epochs)
    end
    epochwidth = markers[2][1] - markers[1][1]
    return imap(x -> edfh.annotations[x[1]:x[1]+epochwidth], markers[1:end-1])
end

"""
    dummyacquire(edfh)

Dummy function for call in writefile! for optional acquire function
If using package for data acquisition will need to custom write the acquire function
for your calls to writefile!
"""
dummyacquire(edfh) = 0

"""
    channeltimesegment(edfh, channel, startsec, endsec, physical)

get the channel's data between the time points
"""
function channeltimesegment(edfh, channel, startsec, endsec, physical)
    sigdata = signaldata(edfh)
    if startsec >= endsec || startsec > edfh.file_duration
        @warn("bad parameters for channeltimesegment")
        return sigdata[1,end:1]  # empty but type correct
    elseif endsec > edfh.file_duration
        @warn("bad end parameter for channeltimesegment")
        endsec = edfh.file_duration
        startsec = edfh.file_duration - edfh.datarecord_duration
    end
    row1, col1 = signalat(edfh, startsec, channel)
    row2, col2 = signalat(edfh, endsec, channel)
    multiplier = edfh.signalparam[channel].bitvalue
    if row1 == row2
        return physical ? sigdata[row1,col1:col2] .* multiplier : sigdata[row1,col1:col2]
    end
    startpos, endpos = signalindices(edfh, channel)
    row1data = sigdata[row1, col1:endpos]
    row2data = sigdata[row2, startpos:col2]
    if row2 - row1 > 1
        otherdata = sigdata[row1+1:row2-1, startpos:endpos]
        row2data = vcat(vec(otherdata), row2data)
    end
    if physical
        return collect(Base.Iterators.flatten((vcat(row1data, row2data) .* multiplier)))
    else
        return collect(Base.Iterators.flatten(vcat(row1data, row2data)))
    end
end

"""
    multichanneltimesegment(edfh, chanlist, startsec, endsec, physical)

Get an multichannel array of lists of datapoints over time segment
NB: best if all datapoint signal rates are the same
"""
function multichanneltimesegment(edfh, chanlist, startsec, endsec, physical)
    mdata::Array{Array{Float64,1}} = []
    for chan in chanlist
        push!(mdata, channeltimesegment(edfh, chan, startsec, endsec, physical))
    end
    return mdata
end

"""
    signalindices(edfh, channelnumber)

Get a pair of indices of a channel's bytes within each of the data records
"""
function signalindices(edfh, channelnumber)
    startcol = Int(floor(edfh.signalparam[channelnumber].bufoffset / bytesperdatapoint(edfh))) + 1
    endcol = startcol + edfh.signalparam[channelnumber].smp_per_record - 1
    return (startcol, endcol)
end

"""
    digitalchanneldata(edfh, channelnumber)

Get a single digital channel of data in its entirety.
# Arguments:
- edfh          the BEDFPlus struct
- channelnumber the channel number in the records
"""
function digitalchanneldata(edfh, channelnumber)
    span = signalindices(edfh, channelnumber)
    data = signaldata(edfh)[:, span[1]:span[2]]
    return vec(data)
end

"""
    physicalchanneldata(edfh, channelnumber)

Get a single data channel in its entirely, in the physical units stated in the header
# Arguments
- edfh          the BEDFPlus struct
- channelnumber the channel number in the records-- a channel in the mapped_signals list
"""
function physicalchanneldata(edfh, channel)
    if !(channel in edfh.mapped_signals)
        throw("physicalchanneldata($channel) is not a signal channel")
    end
    digdata = digitalchanneldata(edfh, channel)
    if length(digdata) < 1
        return digdata
    end
    return digdata * edfh.signalparam[channel].bitvalue
end

"""
    samplerate(edfh, channel)

Get sample (sampling) rate (fs) on the channel in sec^-1 units
"""
samplerate(edfh, channel) = edfh.signalparam[channel].smp_per_record / edfh.datarecord_duration


"""
    notchfilter(signals, fs, notchfreq=60, q = 35)

Notch filter signals in array signals, return filtered signals
"""
function notchfilter(signals, fs, notchfreq=60, q = 35)
    wdo = 2.0notchfreq/fs
    return filtfilt(iirnotch(wdo, wdo/q), signals)
end

"""
    highpassfilter(signals, fs, cutoff=1.0, order=4)

Apply high pass filter to signals, return filtered data
"""
function highpassfilter(signals, fs, cutoff=1.0, order=4)
    wdo = 2.0cutoff/fs
    filth = digitalfilter(Highpass(wdo), Butterworth(order))
    return filtfilt(filth, signals)
end

"""
    lowpassfilter(signals, fs, cutoff=25.0, order=4)

Apply low pass filter to signals, return filtered data
"""
function lowpassfilter(signals, fs, cutoff=25.0, order=4)
    wdo = 2.0cutoff/fs
    filtl = digitalfilter(Lowpass(wdo), Butterworth(order))
    return filtfilt(filtl, signals)
end

"""
    closefile!(edfh)

Close the file opened by loadfile and loaded to the BEDFPlus struct
Releases memory from read data in edfh
"""
function closefile!(edfh)
    edfh.EDFsignals = Array{Int16,2}(undef, 0, 0)
    edfh.BDFsignals = Array{Int32,2}(undef, 0, 0)
    close(edfh.ios)
    edfh.filetype = CLOSED
    return 0
end

"""
    readdata!(edfh)

Helper function for loadfile, reads signal data into the BEDFPlus struct
"""
function readdata!(edfh)
    signalpoints = 0
    for chan in 1:edfh.channelcount
        signalpoints += edfh.signalparam[chan].smp_per_record
    end
    if edfh.bdf || edfh.bdfplus
        edfh.BDFsignals = zeros(Int32, (edfh.datarecords, signalpoints))
    else
        edfh.EDFsignals = zeros(Int16, (edfh.datarecords, signalpoints))
    end
    seek(edfh.ios, edfh.headersize)
    for i in 1:edfh.datarecords
        columnstart = 1
        for j in 1:edfh.channelcount
            cbuf = read(edfh.ios, (edfh.signalparam[j].smp_per_record * bytesperdatapoint(edfh)))
            if edfh.bdf || edfh.bdfplus
                for k in 1:3:length(cbuf)-1
                    edfh.BDFsignals[i,columnstart] = Int(reinterpret(Int24, cbuf[k:k+2])[1])
                    columnstart += 1
                end
            else
                intarray = reinterpret(Int16, cbuf)
                edfh.EDFsignals[i, columnstart:columnstart+length(intarray)-1] = intarray
                columnstart += length(intarray)
            end
        end
    end
    return 0
end

"""
    signaldata(edfh)

Return which BEDFPlus variable holds the signal data
"""
signaldata(edfh) = (edfh.bdf || edfh.bdfplus) ? edfh.BDFsignals : edfh.EDFsignals

"""
    recordslice(edfh, startpos, endpos)

Get a slice of the data in the recording from one data entry position to another
"""
recordslice(edfh, startpos, endpos) = signaldata(edfh)[startpos:endpos, :]

"""
    bytesperdatapoint(edfh)

Return how many bytes used per data point entry: 2 for EDF (16-bit), 3 for BDF (24-bit) data.
"""
bytesperdatapoint(edfh) = (edfh.bdfplus || edfh.bdf ) ? 3 : 2

"""
    datapointinterval(edfh, channel)

Time interval in fractions of a second between individual signal data points
"""
function datapointinterval(edfh, channel=edfh.mapped_signals[1])
    return edfh.datarecord_duration / edfh.signalparam[channel].smp_per_record
end


"""
    recordindexat(edfh, secondsafterstart)

Index of the record point at or closest just before a given time from recording start
Translates a values in seconds to a position in the signal data matrix,
returns that record's position
"""
function recordindexat(edfh, secondsafterstart)
    if edfh.discontinuous && edfh.edfplus
        # for EDF+D need to go on annotations about times
        for i in 2:edfh.datarecords
            firstannot = edfh.annotations[i][1]
            if secondsafterstart < firstannot.onset
                return i - 1
            end
        end
        return edfh.datarecords  # last one is at end
    end
    # for continuous files, we do not need to check any annotation channel
    pos = Int(floor(secondsafterstart/Float64(edfh.datarecord_duration))) + 1
    return pos > edfh.datarecords ? edfh.datarecords : pos
end


"""
    signalat(edfh, secondsafter, channel)

Get the position in the signal data of the data point at or closest after a
given time from recording start. Translates a value in seconds to a position
in the signal channel matrix, returns that signal data point's 2D position as list
"""
function signalat(edfh, secondsafter, channel=edfh.mapped_signals[1])
    row = recordindexat(edfh, secondsafter)
    seconddiff = round(secondsafter - edfh.datarecord_duration * (row-1), digits=4)
    seconddif = seconddiff < 0.0 ? 0.0 : seconddiff
    startpos, endpos = signalindices(edfh, channel)
    col = startpos + Int(floor(seconddiff / datapointinterval(edfh, channel)))
    return (row, col > endpos ? endpos : col)
end

"""
    epochmarkers(edfh, secs)

Get a set of (start, stop) positional markers for epochs (sequential windows)
given an epoch duration in seconds
"""
epochmarkers(edfh, secs) = map(t->signalat(edfh,t), 0:secs:edfh.file_duration)

""" replace underslashes in string with spaces """
dash2space(x) = replace(x, "_" => " ")

"""
    checkfile!(edfh)

Helper function for loadfile and writefile!
"""
function checkfile!(edfh)
    function throwifhasforbiddenchars(bytes)
        if something(findfirst(c -> (Int(c) < 32) || (Int(c) > 126), bytes), 0) > 0
            throw("Control type forbidden chars in header")
        end
    end
    seekstart(edfh.ios)
    try
        hdrbuf = read!(edfh.ios, Array{UInt8}(undef, 256))      # check header
        if hdrbuf[1:8] == b"\xffBIOSEMI"                        # version bdf
            edfh.bdf = true
            edfh.edf = false
            edfh.version = "BIOSEMI"
            edfh.filetype = BDF
        elseif hdrbuf[1:8] == b"0       "                       # version edf
            edfh.bdf = false
            edfh.edf = true
            edfh.version = ""
            edfh.filetype = EDF
        else
            throw("identification code error")
        end
        throwifhasforbiddenchars(hdrbuf[9:88])
        edfh.patient = String(trim(hdrbuf[9:88]))      # patient
        throwifhasforbiddenchars(hdrbuf[89:168])
        edfh.recording = String(trim(hdrbuf[89:168]))  # recording
        throwifhasforbiddenchars(hdrbuf[169:176])
        datestring = String(trim(hdrbuf[169:176]))     # start date
        date = Date(datestring, "dd.mm.yy")
        if Dates.year(date) < 84
            date += Dates.Year(2000)
        else
            date += Dates.Year(1900)
        end
        edfh.startdate_day = Dates.day(date)
        edfh.startdate_month = Dates.month(date)
        edfh.startdate_year = Dates.year(date)
        throwifhasforbiddenchars(hdrbuf[177:184])
        timestring = String(hdrbuf[177:184])           # start time
        mat = match(r"(\d\d).(\d\d).(\d\d)", timestring)
        starttime_hour, starttime_minute, starttime_second = mat.captures
        edfh.starttime_hour = parse(Int, trim(starttime_hour))
        edfh.starttime_minute = parse(Int, trim(starttime_minute))
        edfh.starttime_second = parse(Int, trim(starttime_second))
        throwifhasforbiddenchars(hdrbuf[185:192])
        headersize = parse(Int32, trim(hdrbuf[185:192]))        # edf header size, changes with channels
        edfh.headersize = headersize
        throwifhasforbiddenchars(hdrbuf[193:236])
        subtype = String(trim(hdrbuf[193:197]))        # subtype or version of data format
        if edfh.edf
            if subtype == "EDF+C"
                edfh.filetype = EDFPLUS
                edfh.edfplus = true
                edfh.edf = false
                edfh.bdf = false
                edfh.bdfplus = false
                edfh.discontinuous = false
            elseif subtype == "EDF+D"
                edfh.filetype = EDFPLUS
                edfh.edfplus = true
                edfh.edf = false
                edfh.bdf = false
                edfh.bdfplus = false
                edfh.discontinuous = true
            else
                edfh.edfplus = false
            end
        elseif edfh.bdf
            if subtype == "BDF+C"
                edfh.filetype = BDFPLUS
                edfh.bdfplus = true
                edfh.edfplus = false
                edfh.edf = false
                edfh.bdf = false
                edfh.discontinuous = false
            elseif subtype == "BDF+D"
                edfh.filetype = BDFPLUS
                edfh.bdfplus = true
                edfh.edfplus = false
                edfh.edf = false
                edfh.bdf = false
                edfh.discontinuous = true
            else
                edfh.bdfplus = false
            end
        end
        throwifhasforbiddenchars(hdrbuf[237:244])
        edfh.datarecords = parse(Int32, trim(hdrbuf[237:244]))  # number of data records
        if edfh.datarecords < 1
            println("Record count was unknown or invalid: $(trim(hdrbuf[237:244]))")
        end
        throwifhasforbiddenchars(hdrbuf[245:252])
        edfh.datarecord_duration = parse(Float32, trim(hdrbuf[245:252])) # datarecord duration in seconds
        throwifhasforbiddenchars(hdrbuf[253:256])
        edfh.channelcount = parse(Int16, trim(hdrbuf[253:256]))  # number of data signals or records (channels) in file
        if edfh.channelcount < 0
            throw("bad channel count")
        end
        calcheadersize = (edfh.channelcount + 1) * 256
        if calcheadersize != headersize
            throw("Bad header size: in file as $headersize, calculates to be $calcheadersize")
        end

    catch y
        @warn("$y\n")
        edfh.filetype = FORMAT_ERROR
        return edfh
    end
    # process the signal characteristics in the header after reading header into hdrbuf
    seekstart(edfh.ios)
    edfh.recordsize = 0
    multiplier = bytesperdatapoint(edfh)
    try
        hdrbuf = read!(edfh.ios, Array{UInt8}(undef, (edfh.channelcount + 1) * 256))
        for i in 1:edfh.channelcount  # loop over channel signal parameters
            pblock = ChannelParam()
            # channel label gets special handling since it might indicate an annotations channel
            pos = 257 + (i-1) * 16
            channellabel = String(hdrbuf[pos:pos+15])
            throwifhasforbiddenchars(channellabel)
            pblock.label = channellabel                          # channel label in ASCII, eg "Fp1"
            if (edfh.edfplus && channellabel == "EDF Annotations ") ||
               (edfh.bdfplus && channellabel == "BDF Annotations ")
                edfh.annotationchannel = i
                pblock.annotation = true
            else
                push!(edfh.mapped_signals, i)
            end
            pos = 257 + edfh.channelcount * 16 + (i-1) * 80
            transducertype = String(hdrbuf[pos:pos+79])
            throwifhasforbiddenchars(transducertype)
            pblock.transducer = transducertype                    # transducer type eg "active electrode"
            if pblock.annotation && something(findfirst(c->!isspace(c), transducertype), 0) > 0
                throw("Transducer field should be blank in annotation channels")
            end
            pos = 257 + edfh.channelcount * 96 + (i-1) * 8
            pblock.physdimension = String(trim(hdrbuf[pos:pos+7]))   # physical dimensions eg. "uV"
            pos = 257 + edfh.channelcount * 104 + (i-1) * 8
            pblock.physmin = parse(Float32, trim(hdrbuf[pos:pos+7])) # physical minimum in above dimensions
            pos = 257 + edfh.channelcount * 112 + (i-1) * 8
            pblock.physmax = parse(Float32, trim(hdrbuf[pos:pos+7]))  # physical maximum in above dimensions
            pos = 257 + edfh.channelcount * 120 + (i-1) * 8
            pblock.digmin = parse(Float32, trim(hdrbuf[pos:pos+7]))   # digital minimum in above dimensions
            pos = 257 + edfh.channelcount * 128 + (i-1) * 8
            pblock.digmax = parse(Float32, trim(hdrbuf[pos:pos+7]))   # digital maximum in above dimensions
            if edfh.edfplus && pblock.annotation && (pblock.digmin != -32768 || pblock.digmax != 32767)
                throw("edfplus annotation data entry should have the digital min parameters set to extremes")
            elseif edfh.bdfplus && pblock.annotation && (pblock.digmin != -8388608 || pblock.digmax != 8388607)
                throw("bdf annotation data entry should have the digital max parameters set to extremes")
            elseif edfh.edf && (pblock.digmin < -32768 || pblock.digmin > 32767 ||
                   pblock.digmax < -32768 || pblock.digmax > 32767)
                throw("edf digital parameter out of range")
            elseif edfh.bdf && (pblock.digmin < -8388608 || pblock.digmin > 8388607 ||
                   pblock.digmax < -8388608 || pblock.digmax > 8388607)
                throw("bdf digital parameter out of range")
            end
            pos = 257 + edfh.channelcount * 136 + (i-1) * 80
            pfchars = String(hdrbuf[pos:pos+79])
            throwifhasforbiddenchars(pfchars)
            pblock.prefilter = pfchars                            # prefilter field eg "HP:DC"
            if pblock.annotation && something(findfirst(c->!isspace(c), pfchars), 0) > 0
                throw("Prefilter field should be blank in annotation channels")
            end
            pos = 257 + edfh.channelcount * 216 + (i-1) * 8
            pblock.smp_per_record = parse(Int, trim(hdrbuf[pos:pos+7])) # number of samples of this channel per data record
            edfh.recordsize += pblock.smp_per_record * multiplier
            pos = 257 + edfh.channelcount * 224 + (i-1) * 32
            reserved = String(hdrbuf[pos:pos+31])        # reserved text field
            throwifhasforbiddenchars(reserved)
            pblock.reserved = reserved
            if pblock.annotation == false && !edfh.edfplus && !edfh.bdfplus &&
                                             edfh.datarecord_duration < 0.0000001
                    throw("signal data may be mislabeled")
            end
            push!(edfh.signalparam, pblock)
        end
        if edfh.bdf && edfh.recordsize > 1578640
            throw("record size too large for a bdf file")
        elseif edfh.recordsize > 10485760
            throw("Record size too large for an edf file")
        end

    catch y
        @warn("checkfile!\n$y\n")
        edfh.filetype = FORMAT_ERROR
        return edfh
    end

    #=
    from https://www.edfplus.info/specs/edfplus.html#header, December 2017:
    The 'local patient identification' field must start with the subfields
         (subfields do not contain, but are separated by, spaces):
    - the code by which the patient is known in the hospital administration.
    - sex (English, so F or M).
    - birthdate in dd-MMM-yyyy format using the English 3-character abbreviations
      of the month in capitals. 02-AUG-1951 is OK, while 2-AUG-1951 is not.
    - the patients name.
     Any space inside the hospital code or the name of the patient must be replaced
     by a different character, for instance an underscore. For instance, the
     'local patient identification' field could start with:
     MCH-0234567 F 02-MAY-1951 Haagse_Harry.
    Subfields whose contents are unknown, not applicable or must be made
    anonymous are replaced by a single character 'X'.
    Additional subfields may follow the ones described here.
    =#
    try
        if edfh.edfplus || edfh.bdfplus
            patienttxt = String(edfh.patient)
            subfield = split(patienttxt)
            if length(subfield) < 4
                throw("Plus patient identification lacking enough fields")
            end
            edfh.patientcode = subfield[1][1] == 'X' ? "" : subfield[1]
            edfh.patientcode = dash2space(edfh.patientcode)
            if subfield[2] != "M" && subfield[2] != "F" && subfield[2] != "X"
                throw("patient identification second field must be X, F or M")
            elseif subfield[2] == "M"
                edfh.gender = "Male"
            elseif subfield[2] == "F"
                edfh.gender = "Female"
            else
                edfh.gender = ""
            end
            if subfield[3] == "X"
                edfh.birthdate = ""
            elseif !occursin(subfield[3][4:6], "JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC") ||
                   !(Date(subfield[3], "dd-uuu-yyyy") isa Date)
                @warn(subfield[3][4:6])
                throw("Bad birthdate field in patient identification")
            else
                edfh.birthdate = subfield[3]
            end
            edfh.patientname = dash2space(subfield[4])
            if length(subfield) > 4
                edfh.patient_additional = join(subfield[5:end], " ")
            end
    #=
     The 'local recording identification' field must start with the subfields
     (subfields do not contain, but are separated by, spaces):
     - The text 'Startdate'.
     - The startdate itself in dd-MMM-yyyy format using the English 3-character
       abbreviations of the month in capitals.
     - The hospital administration code of the investigation, i.e. EEG number or PSG number.
     - A code specifying the responsible investigator or technician.
     - A code specifying the used equipment.
     Any space inside any of these codes must be replaced by a different character,
     for instance an underscore. The 'local recording identification' field could
     contain: Startdate 02-MAR-2002 PSG-1234/2002 NN Telemetry03.
     Subfields whose contents are unknown, not applicable or must be made anonymous
     are replaced by a single character 'X'. So, if everything is unknown then the
     'local recording identification' field would start with:
     Startdate X X X X. Additional subfields may follow the ones described here.
    =#
            subfield = split(trim(edfh.recording))
            if length(subfield) < 5
                throw("Not enough fields in plus recording data")
            elseif subfield[2] == "X"
                edfh.startdatestring = ""
            else
                edfh.startdatestring = subfield[2]
                dor = Date(subfield[2], "dd-uuu-yyyy")
                if Dates.year(dor) < 1970
                    throw("bad startdate year in recording data: got $(Dates.year(dor)) from $dor")
                else
                    edfh.startdate_year = Dates.year(dor)
                end
            end
            if subfield[3] == ""
                edfh.admincode = ""
            else
                edfh.admincode = dash2space(subfield[3])
            end
            if subfield[4] == ""
                edfh.technician = ""
            else
                edfh.technician = dash2space(subfield[4])
            end
            if subfield[5] == ""
                edfh.equipment = ""
            else
                edfh.equipment = dash2space(subfield[5])
            end
            if length(subfield) > 5
                edfh.recording_additional = dash2space(join(subfield[6:end], " "))
            end
        end

        edfh.headersize = (edfh.channelcount + 1) * 256
        filsiz = filesize(edfh.path)   # get actual file size
        filszbyhdr = edfh.recordsize * edfh.datarecords + edfh.headersize
        if filsiz != filszbyhdr
            throw("file size is not compatible with header information: header says $filszbyhdr, filesystem $filsiz")
        end
    catch y
        @warn("$y\n")
        edfh.filetype = FORMAT_ERROR
        return 0
    end

    n = 0
    for i in 1:edfh.channelcount
        edfh.signalparam[i].bufoffset = n
        n += edfh.signalparam[i].smp_per_record * bytesperdatapoint(edfh)
        edfh.signalparam[i].bitvalue = (edfh.signalparam[i].physmax - edfh.signalparam[i].physmin) /
                                      (edfh.signalparam[i].digmax - edfh.signalparam[i].digmin)
        edfh.signalparam[i].offset = edfh.signalparam[i].physmax / edfh.signalparam[i].bitvalue -
                                    edfh.signalparam[i].digmax
    end
    return edfh
end

"""
    readannotations!(edfh)

Helper function for loadfile
"""
function readannotations!(edfh)
    samplesize = bytesperdatapoint(edfh)
    chan = edfh.annotationchannel
    max_tal_ln = edfh.signalparam[chan].smp_per_record * samplesize
    seek(edfh.ios, (edfh.channelcount + 1) * 256)
    edfh.annotations = Array{Array{Annotation,1},1}(undef, edfh.datarecords)
    fill!(edfh.annotations,[])
    added = 0
    for i in 1:edfh.datarecords
        try
            cnvbuf = read(edfh.ios, (edfh.recordsize))
            startpos = edfh.signalparam[chan].bufoffset + 1
            endpos = startpos + edfh.signalparam[chan].smp_per_record * samplesize -1
            annotbuf = String(cnvbuf[startpos:endpos])
            for (j, tal) in enumerate(split(annotbuf, '\x00'))
                if tal == "" || something(findfirst(isequal('\x14'), tal), 0) < 1
                    break # padding zeroes at end or bad txt
                end
                (times, annottxt) = split(tal, '\x14', limit=2)
                if something(findfirst(isequal('\x15'), times), 0) > 0
                    (ons, duration) = split(times, '\x15')
                    onset = parse(Float64, ons)
                else
                    onset = parse(Float64, times)
                    duration = ""
                end
                if j == 1 && i == 1 # first TAL of first record
                    edfh.starttime_subsecond = onset
                end
                if length(annottxt) > MAX_ANNOTATION_LENGTH
                    @warn("Annotation at $i $j $k too long, truncated")
                    annottxt = annottxt[1:MAX_ANNOTATION_LENGTH]
                end
                newannot = Annotation(onset, duration, split(annottxt, '\x14'))
                edfh.annotations[i] = vcat(edfh.annotations[i], newannot)
                added += 1
            end
        catch y
            @warn("Error reading annotation in record $i: $y\n")
            return -1
        end
    end
    return added
end

"""
    translate24to16bits!(edfh)

Translate data in 24-bit BDF to 16-bit EDF format
Helper function for writefile!
"""
function translate24to16bits!(edfh)
    data = edfh.BDFsignals
    cvrtfactor = min(abs(32767/maximum(data)), abs(-32768/minimum(data)))
    if cvrtfactor < 1.0
        edfh.EDFsignals = map(x->Int16(floor(x * cvrtfactor)), data)
        for chan in edfh.mapped_signals
            edfh.signalparam[chan].physmin /= cvrtfactor
            edfh.signalparam[chan].physmax /= cvrtfactor
        end
    else
        edfh.EDFsignals = map(x->Int16(x), data)
    end

    achan = edfh.annotationchannel
    if achan == 0
        @warn("No annotation channel in source file")
        return -1
    end
    startcol = Int(edfh.signalparam[achan].bufoffset / 3) + 1
    endcol = startcol + edfh.signalparam[achan].smp_per_record - 1
    for rec in 1:edfh.datarecords
        arr::Array{UInt8} = []
        oby = reinterpret(UInt8, edfh.BDFsignals[rec, startcol:endcol])
        for (i, cha) in enumerate(oby)
            if i % 4 != 0
                push!(arr, cha)
            end
        end
        newspace = endcol - startcol + 1
        if length(arr) > 2newspace
            arr = arr[1:2newspace]
            arr[end] = '\x00'
        else
            while length(arr) < 2newspace
                push!(arr, '\x00')
            end
        end
        edfh.EDFsignals[rec, startcol:endcol] .= reinterpret(Int16, arr)
    end
    return 0
end

"""
    translate16to24bits!(edfh)

Translate 16 bit data to 32-bit width, for change to 24-bit data for writefile!
"""
function translate16to24bits!(edfh)
    edfh.BDFsignals = map(x->Int32(x), edfh.EDFsignals)
    chan = edfh.annotationchannel
    if chan == 0
        @warn("No annotation channel in source file")
        return -1
    end
    startcol = Int(edfh.signalparam[chan].bufoffset / 2) + 1
    endcol = startcol + edfh.signalparam[chan].smp_per_record - 1
    for rec in 1:edfh.datarecords
        arr::Array{Int32,1} = []
        oby = vec(reinterpret(UInt8, edfh.EDFsignals[rec, startcol:endcol]))
        for k in 1:3:length(oby) - 1
            push!(arr, reinterpret(Int32,[oby[k], oby[k+1], oby[k+2], UInt8(0)])[1])
        end
        while length(arr) < endcol - startcol + 1
            push!(arr,0)
        end
        edfh.BDFsignals[rec, startcol:endcol] .= arr
    end
    return 0
end

"""
    writeEDFsignalchannel(edfh, fh, record, channel)

Helper function for writefile!
write a record's worth of a signal channel at given record and channel number
"""
function writeEDFsignalchannel(edfh, fh, record, channel)
    (startpos, endpos) = signalindices(edfh, channel)
    signals = edfh.EDFsignals[record,startpos:endpos]
    return write(fh, signals)
end

"""
    writeBDFsignalchannel(edfh,fh, record, channel)

Helper function for writefile!
write a BDF record's worth of a signal channel at given record and channel number
"""
function writeBDFsignalchannel(edfh,fh, record, channel)
    (startpos, endpos) = signalindices(edfh, channel)
    signals = vec(edfh.BDFsignals[record,startpos:endpos])
    written = 0
    for sig in signals
        written += writei24(fh, sig)
    end
    return written
end

"""
    writeEDFrecords!(edfh, fh)

Helper function for writefile!
Write a record's worth of all channels of a given record
"""
function writeEDFrecords!(edfh, fh)
    if isempty(edfh.EDFsignals)
        # write data as EDF -- if was BDF adjust width if needed for 24 to 16 bits
        if (edfh.bdf || edfh.bdfplus) && !isempty(edfh.BDFsignals)
            translate24to16bits!(edfh)
        else
            return 0
        end
    end
    written = 0
    for i in 1:edfh.datarecords
        for j in 1:edfh.channelcount
            written += writeEDFsignalchannel(edfh, fh, i, j)
        end
    end
    return written
end

"""
    writeBDFrecords!(edfh, fh)

Write an BEDFPlus format file
Helper file for writefile!
"""
function writeBDFrecords!(edfh, fh)
    if isempty(edfh.BDFsignals)
        # write data as BDF -- if was EDF adjust width if needed for 16 to 24 bits
        if (edfh.ef || edfh.edfplus) && !isempty(edfh.EDFsignals)
            translate16to24bits!(edfh)
        else
            return 0
        end
    end
    written = 0
    for i in 1:edfh.datarecords
        for j in eachindex(edfh.signalparam)
            written += writeBDFsignalchannel(edfh, fh, i, j)
        end
    end
    return written
end

"""
    writeheader(edfh, fh)

Helper function for writefile!
"""
function writeheader(edfh::BEDFPlus, fh::IOStream)
    written = 0
    channelcount = edfh.channelcount
    if channelcount < 0
        throw("Channel count is negative")
    elseif channelcount > MAX_CHANNELS
         throw("Channel count $channelcount is too large")
    end
    for i in 1:channelcount
        if edfh.signalparam[i].smp_per_record < 1
            throw("negative samples per record")
        elseif edfh.signalparam[i].digmax <= edfh.signalparam[i].digmin
            throw("digmax must be > digmin")
        elseif edfh.signalparam[i].physmax <= edfh.signalparam[i].physmin
            throw("physical max must be > physical min")
        end
    end
    seekstart(fh)
    if edfh.edf || edfh.edfplus
        written += write(fh, "0       ")
    else
        written += write(fh, b"\xffBIOSEMI")
    end
    pidbytes = edfh.patientcode == "" ? "X " : dash2space(edfh.patientcode) * " "
    if edfh.gender == ""
        pidbytes *= "X "
    elseif edfh.gender[1] == 'M'
        pidbytes *= "M "
    elseif edfh.gender[1] == 'F'
        pidbytes *= "F "
    else
        pidbytes *= "X "
    end
    if edfh.birthdate != ""
        pidbytes *= edfh.birthdate * " "
    else
        pidbytes *= "X "
    end
    if edfh.patientname == ""
        pidbytes *= "X "
    else
        pidbytes *= dash2space(edfh.patientname) * " "
    end
    if edfh.patient_additional != ""
        pidbytes *= dash2space(edfh.patient_additional)
    end
    if length(pidbytes) > 80
        pidbytes *= pidbytes[1:80]
    else
        while length(pidbytes) < 80
            pidbytes *= " "
        end
    end
    written += write(fh, pidbytes)

    if edfh.startdate_year != 0
        date = DateTime(edfh.startdate_year, edfh.startdate_month, edfh.startdate_day,
                        edfh.starttime_hour, edfh.starttime_minute, edfh.starttime_second)
    else
        date = now()
    end

    ridbytes = "Startdate " * uppercase(Dates.format(date, "dd-uuu-yyyy")) * " "
    if edfh.admincode == ""
        ridbytes *= "X "
    else
        ridbytes *= dash2space(edfh.admincode) * " "
    end
    if edfh.technician == ""
        ridbytes *= "X "
    else
        ridbytes *= dash2space(edfh.technician) * " "
    end
    if edfh.equipment == ""
        ridbytes *= "X "
    else
        ridbytes *= dash2space(edfh.equipment) * " "
    end
    if edfh.recording_additional != ""
        ridbytes *= dash2space(edfh.recording_additional)
    end
    written += writeleftjust(fh, ridbytes, 80)
    startdate = Dates.format(date, "dd.mm.yy")
    written += write(fh, startdate)
    starttime = Dates.format(date, "HH.MM.SS")
    written += write(fh, starttime)
    hdsize = Int((channelcount + 1) * 256)           # bytes in header
    written += writeleftjust(fh, hdsize, 8)

    if edfh.edfplus
        if edfh.discontinuous
            written += writeleftjust(fh, "EDF+D", 44)
        else
            written += writeleftjust(fh, "EDF+C", 44)
        end
    elseif edfh.bdfplus
        written += writeleftjust(fh, "BDF+C", 44)
    else
        written += writeleftjust(fh, "     ", 44)
    end
    # header initialized to -1 in duration in case data is not finalized yet
    # This must be updated when final channel data is written.
    if edfh.datarecords > 0
        written += writeleftjust(fh, edfh.datarecords, 8)
    else
        written += writeleftjust(fh, "-1      ", 8)
    end
    if floor(edfh.datarecord_duration) == edfh.datarecord_duration
        written += writeleftjust(fh, Int(edfh.datarecord_duration), 8)
    else
        durstr = trimrightzeros("$(edfh.datarecord_duration)")
        written += writeleftjust(fh, durstr, 8)
    end
    written += writeleftjust(fh, channelcount, 4)
    for i in 1:channelcount
        if edfh.signalparam[i].annotation
            if edfh.edfplus
                written += write(fh, "EDF Annotations ")
            elseif edfh.bdfplus
                written += write(fh, "BDF Annotations ")
            end
        else
            written += writeleftjust(fh, edfh.signalparam[i].label, 16)
        end
    end
    for i in 1:channelcount
        written += writeleftjust(fh, edfh.signalparam[i].transducer, 80)
    end
    for i in 1:channelcount
        written += writeleftjust(fh, edfh.signalparam[i].physdimension, 8)
    end
    for i in 1:channelcount
        written += writeleftjust(fh, edfh.signalparam[i].physmin, 8)
    end
    for i in 1:channelcount
        written += writeleftjust(fh, edfh.signalparam[i].physmax, 8)
    end
    for i in 1:channelcount
        written += writeleftjust(fh, edfh.signalparam[i].digmin, 8)
    end
    for i in 1:channelcount
        written += writeleftjust(fh, edfh.signalparam[i].digmax, 8)
    end
    for i in 1:channelcount
        written += writeleftjust(fh, edfh.signalparam[i].prefilter, 80)
    end
    for i in 1:channelcount
        written += writeleftjust(fh, edfh.signalparam[i].smp_per_record, 8)
    end
    for i in 1:channelcount
        written += writeleftjust(fh, edfh.reserved, 32)
    end
    return written
end

"""
   annotationtoTAL(ann)

Create a TAL (timestamped annotation list) text entry out of an annotation
"""
function annotationtoTAL(ann)
    txt = trimrightzeros((ann.onset >= 0 ? "+" : "") * "$(ann.onset)")
    if ann.duration != ""
        txt *= "\x15" * ann.duration
    end
    if length(ann.annotation) > 1
        anntxt = join(ann.annotation, "\x14")
    else
        anntxt = first(ann.annotation)
    end
    txt *= "\x14" * anntxt * "\x14\x00"
    return txt
end

"""
    addannotation!(edfh, onset, duration, description)

Add an annotation at the given onset timepoint IF there is room
Note the description arg is a text string, not an array argument here
"""
function addannotation!(edfh, onset, duration, description)
    if isempty(edfh.annotations) && edfh.annotationchannel == 0
        throw("No annotation channels in file")
    end
    if typeof(duration) != String
        duration = "$duration"
    end
    if typeof(description) <: Array
        if length(description) > 1
            anntxt = join(description, "\x14")
        else
            anntxt = description[1]
        end
    else
       anntxt = description
    end
    if length(anntxt) > MAX_ANNOTATION_LENGTH
        anntxt = anntxt[1:MAX_ANNOTATION_LENGTH]
    end
    anntxt = latintoascii(replace(anntxt, r"[\x00-\x13]" => "."))
    newannot = Annotation(onset, duration, anntxt)
    neartimeindex = recordindexat(edfh, onset)
    toadd = annotationtoTAL(newannot)
    additionalbytes = length(toadd)
    iwidth = bytesperdatapoint(edfh)
    (startpos,endpos) = signalindices(edfh, edfh.annotationchannel)
    chanlen = iwidth*(endpos-startpos+1)
    if additionalbytes > chanlen - 6
        @warn("TAL is too large for adding to a channel of length $chanlen bytes")
    else
        for recordnum in neartimeindex:edfh.datarecords-1
            ints = vec(signaldata(edfh)[recordnum, startpos:endpos])
            ctxt = reinterpret(UInt8, ints)
            addindex = something(findlast(c -> c != 0, ctxt), 0) + 1
            if addindex > 0 && addindex + additionalbytes < length(ctxt)
                ctxt[addindex+1:addindex+additionalbytes] .= Array{UInt8,1}(toadd)
                if iwidth == 2
                    signaldata(edfh)[recordnum, startpos:endpos] = reinterpret(Int16, ctxt)
                else
                    for k in 1:3:length(ctxt)-1
                        signaldata[recordnum, startpos] = reinterpret(Int24, ctxt[k:k+2])
                        startpos += 1
                    end
                end
                edfh.annotations[recordnum] = vcat(edfh.annotations[recordnum], newannot)
                break
            end
        end
    end
    return 0
end

"""
    trim(str)
trim whitespace fore and aft, as in trim in java etc
"""
trim(str::String) = strip(str)
trim(ch::Char) = strip(string(ch))
trim(bytes) = strip(String(bytes))

""" trimrightzeros compact number string by trimming nonsignificant decimal places/point when not zero """
trimrightzeros(fstr) =  (x = parse(Float64, fstr); return x == floor(x) ? string(Int64(x)) : string(x))

"""
    writeleftjust(fh, x, len, fillchar=' ')

Write a stringified object to a file in the leftmost portion of chars written,
filling with fillchar to len length as needed, truncate if too long for field
"""
function writeleftjust(fh, x, len, fillchar=' ')
    str = "$x"
    if length(str) > len
        str = str[1:len]
    else
        str = rpad(str, len, fillchar)
    end
    return write(fh, str)
end

""" map table for translation of latin extended ascii to plain ascii chars """
const latin_dict = Dict(
'¡'=> '!', '¢'=> 'c', '£'=> 'L', '¤'=> 'o', '¥'=> 'Y',
'¦'=> '|', '§'=> 'S', '¨'=> '`', '©'=> 'c', 'ª'=> 'a',
'«'=> '<', '¬'=> '-', '®'=> 'R', '¯'=> '-',
'°'=> 'o', '±'=> '+', '²'=> '2', '³'=> '3', '´'=> '`',
'µ'=> 'u', '¶'=> 'P', '·'=> '.', '¸'=> ',', '¹'=> '1',
'º'=> 'o', '»'=> '>', '¼'=> '/', '½'=> '/', '¾'=> '/',
'¿'=> '?', 'À'=> 'A', 'Á'=> 'A', 'Â'=> 'A', 'Ã'=> 'A',
'Ä'=> 'A', 'Å'=> 'A', 'Æ'=> 'A', 'Ç'=> 'C', 'È'=> 'E',
'É'=> 'E', 'Ê'=> 'E', 'Ë'=> 'E', 'Ì'=> 'I', 'Í'=> 'I',
'Î'=> 'I', 'Ï'=> 'I', 'Ð'=> 'D', 'Ñ'=> 'N', 'Ò'=> 'O',
'Ó'=> 'O', 'Ô'=> 'O', 'Õ'=> 'O', 'Ö'=> 'O', '×'=> '*',
'Ø'=> 'O', 'Ù'=> 'U', 'Ú'=> 'U', 'Û'=> 'U', 'Ü'=> 'U',
'Ý'=> 'Y', 'Þ'=> 'p', 'ß'=> 'b', 'à'=> 'a', 'á'=> 'a',
'â'=> 'a', 'ã'=> 'a', 'ä'=> 'a', 'å'=> 'a', 'æ'=> 'a',
'ç'=> 'c', 'è'=> 'e', 'é'=> 'e', 'ê'=> 'e', 'ë'=> 'e',
'ì'=> 'i', 'í'=> 'i', 'î'=> 'i', 'ï'=> 'i', 'ð'=> 'd',
'ñ'=> 'n', 'ò'=> 'o', 'ó'=> 'o', 'ô'=> 'o', 'õ'=> 'o',
'ö'=> 'o', '÷'=> '/', 'ø'=> 'o', 'ù'=> 'u', 'ú'=> 'u',
'û'=> 'u', 'ü'=> 'u', 'ý'=> 'y', 'þ'=> 'p', 'ÿ'=> 'y')

"""
    latintoascii(str)

Helper function for writefile! related functions
"""
function latintoascii(str)
    len = length(str)
    if len >= 1
        for i in 1:len
            if haskey(latin_dict, str[i])
                arr = split(str, "")
                arr[i] = string(latin_dict[str[i]])
                str = join(arr, "")
            end
        end
    end
    return str
end

"""
   signalchannel(edfh)

Returns channel number of the Signal channel or nothing if none labeled Signal
"""
signalchannel(edfh) = findfirst(x->x.label == "Signal", edfh.signalparam)

"""
    readBiosemiStatus(edfh)

Export BDF Status channel data.
# Returns a Dict structure containing trigger data in the format:
-   Code => vector of triggerbit (Int),
-   Index => vector of index data (Int)
-   Onset => vector of onset times (Float)
-   Duration => vector of durations (Float)
"""
function readBiosemiStatus(edfh, channel=signalchannel(edfh))
    dict = Dict{String, Vector{Any}}()
    timeperdatum =  edfh.datarecord_duration / edfh.signalparam[channel].smp_per_record
    data = digitalchanneldata(edfh, channel)
    dif = findall(x -> x != 0, diff(data))
    starts = vcat(1, dif)
    dict["Code"] = data[starts]
    dict["Index"] = starts
    dict["Onset"] = starts .* timeperdatum
    dict["Duration"] = (starts .- vcat(dif, length(data))) .* timeperdatum
    return dict
end

end # module
