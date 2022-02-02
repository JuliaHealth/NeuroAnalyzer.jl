# EDFPlus.jl

<img src="https://github.com/wherrera10/EDFPlus.jl/blob/master/docs/src/eeg.png">

Julia for handling BDF+ and EDF+ EEG and similar signal data files.

Heavily influenced by the C EEG library edflib.

## Data Structs and Functions

### Structs:	
    
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

Parameters for each channel in the EEG record.
<br /><br /><br />

    mutable struct Annotation
        onset::Float64
        duration::String
        annotation::Array{String,1}
        Annotation() = new(0.0,"",[])
        Annotation(o,d,arr) = new(o, d, typeof(arr) == String ? [arr] : arr)
    end

These are text strings within the file denoting a time, optionally duration, and a list of notes
about the signal at that particular time in the recording. The first onset time
of the annotation channel gives a fractional second offset adjustment of the
start time of that record, which is specified in whole seconds in the header.
<br /><br /><br />

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

Data struct for EDF, EDF+, BDF, and BDF+ EEG type signal files.
<br /><br /><br />

    DataFormat

enum for types this package handles. Current format for a potential translation is also /same/.
<br /><br /><br />

    FileStatus

enum for type or state of file: type of data detected, whether any errors
<br /><br /><br />

### Functions:

    loadfile(path::String, read_annotations=true)

Load a BDF+ or EDF+ type file.
Takes a pathname. Will ignore annotations if second argument is set false.
Returns a BEDFPlus structure including header and data.
<br /><br /><br />


    writefile!(edfh, newpath; acquire=dummyacquire, sigformat=same)

Write to data in the edfh struct to the file indicated by newpath
Returns the file handle of the file written, opened for reading
NOTE: The header needs to be completely specified at function start except for
the final number of records, which will be updated after all data records
are written. For a system that is recording the data as it is written, the
acquire(edfh) function should write the data according the the header parameters.
Note that if the function converts from BDF to EDF or EDF to BDF, the edfh struct is changed.
<br /><br /><br />


    epoch_iterator(edfh, epochsecs; channels, startsec, endsec, physical)

Make an iterator for EEG epochs of a given duration between start and stop times.
Required arguments:

- edfh: BEDFPlus struct
- epochsecs: second duration of each epoch

Optional arguments:

- channels: List of channel numbers for data, defaults to all signal channels
- startsec: Starting position from 0 at start of file, defaults to file start
- endsec: Ending position in seconds from start of _file_, defaults to file end
- physical: Whether to return data as translated to the physical units, defaults to true
<br /><br /><br />


    annotation_epoch_iterator(edfh, epochsecs; startsec, endsec)

Return an iterator for a group of annotations for a given epoch as in epoch_iterator
<br /><br /><br />


    dummyacquire(edfh)

Dummy function for call in writefile! for optional acquire function
If using package for data acquisition will need to custom write the acquire function
for your calls to writefile!
<br /><br /><br />

    channeltimesegment(edfh, channel, startsec, endsec, physical)

Get the channel's data between the time points
<br /><br /><br />


    multichanneltimesegment(edfh, chanlist, startsec, endsec, physical)

Get an multichannel array of lists of datapoints over time segment
Works best if all datapoint signal rates are the same
<br /><br /><br />


    signalindices(edfh, channelnumber)

Get a pair of indices of a channel's bytes within each of the data records
<br /><br /><br />



    digitalchanneldata(edfh, channelnumber)

Get a single digital channel of data in its entirety.
Arguments:

- edfh:          the BEDFPlus struct
- channelnumber: the channel number in the records
<br /><br /><br />




    physicalchanneldata(edfh, channelnumber)
    
Get a single data channel in its entirely, in the physical units stated in the header
Arguments:

- edfh:          the BEDFPlus struct
- channelnumber: the channel number in the records-- a channel in the mapped_signals list
<br /><br /><br />




    samplerate(edfh, channel)

Get sample (sampling) rate (fs) on the channel in sec^-1 units
<br /><br /><br />


    notchfilter(signals, fs, notchfreq=60, q = 35)

Notch filter signals in array signals, return filtered signals
<br /><br /><br />



    highpassfilter(signals, fs, cutoff=1.0, order=4)

Apply high pass filter to signals, return filtered data
<br /><br /><br />




    lowpassfilter(signals, fs, cutoff=25.0, order=4)

Apply low pass filter to signals, return filtered data
<br /><br /><br />




    closefile!(edfh)

Close the file opened by loadfile and loaded to the BEDFPlus struct
May therefore let GC release memory from read data in edfh
<br /><br /><br />


    readdata!(edfh)

Helper function for loadfile, reads signal data into the BEDFPlus struct
<br /><br /><br />


    signaldata(edfh)

Return which BEDFPlus variable holds the signal data
<br /><br /><br />


    recordslice(edfh, startpos, endpos)

Get a slice of the data in the recording from one data entry position to another
<br /><br /><br />


    bytesperdatapoint(edfh)

Return how many bytes used per data point entry: 2 for EDF (16-bit), 3 for BDF (24-bit) data.
<br /><br /><br />


    datapointinterval(edfh, channel)

Time interval in fractions of a second between individual signal data points
<br /><br /><br />


    recordindexat(edfh, secondsafterstart)

Index of the record point at or closest just before a given time from recording start
Translates a values in seconds to a position in the signal data matrix,
returns that record's position
<br /><br /><br />



    signalat(edfh, secondsafter, channel)

Get the position in the signal data of the data point at or closest after a
given time from recording start. Translates a value in seconds to a position
in the signal channel matrix, returns that signal data point's 2D position as list
<br /><br /><br />



    epochmarkers(edfh, secs)

Get a set of (start, stop) positional markers for epochs (sequential windows)
given an epoch duration in seconds
<br /><br /><br />



    checkfile!(edfh)

Check an input file to be valid EDF/BDF/+ format.
<br /><br /><br />



    readannotations!(edfh)

Read the annotations of an input file into the EDFPlus struct.
<br /><br /><br />



    translate24to16bits!(edfh)

Translate data in 24-bit BDF to 16-bit EDF format
<br /><br /><br />


    translate16to24bits!(edfh)

Translate 16 bit data to 32-bit width, for change to 24-bit data for writefile!
<br /><br /><br />



    addannotation!(edfh, onset, duration, description)

Add an annotation at the given onset timepoint IF there is room
Note the description arg is expected to be a String or similar
<br /><br /><br />


    latintoascii(str)

Change non-ascii characters to similar ascii chars
<br /><br /><br />


    readBiosemiStatus(edfh)

Export BDF Status channel data.

Returns a Dict structure containing trigger data in the format:

-   Code: => vector of triggerbit (Int),
-   Index: => vector of index data (Int)
-   Onset: => vector of onset times (Float)
-   Duration: => vector of durations (Float)
<br /><br /><br />




Installation:

To install from a Julia REPL command line session:

    using Pkg
    Pkg.add("EDFPlus"))

Note that the test files include a 23 mb test file. You may need to allow extra time for that to download when installing.
