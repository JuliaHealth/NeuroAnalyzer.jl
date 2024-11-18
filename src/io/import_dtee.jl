export import_dtee

"""
    import_dtee(file_name; <keyword arguments>)

Load Elmiko dtEE EEG file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""

function import_dtee(file_name::String; detect_type::Bool=true)::NeuroAnalyzer.NEURO

    _wip()

    # using XML

    file_name = "20210805-092420-{db06ce21-b871-45b8-bd00-e2b2f9dfaf91}.dtEE"

    @assert isdir(file_name) "file_name must be a recording folder name."
    cd(file_name)
    files = readdir()
    eef = ""
    for idx in files
        splitext(idx)[2] == ".eef" && (eef = idx)
    end
    @assert eef != "" ".eef file could not be find in the recording folder."
    @assert isfile(eef) ".eef file cannot be loaded."
    eef = XML.read(eef, Node)
    id = children(children(children(children(children(eef)[end])[1])[1])[1])[1]["id"]
    recording_date = children(children(children(children(children(eef)[end])[1])[1])[1])[1]["creationTime"]
    age = children(children(children(children(children(eef)[end])[1])[1])[2])[1]["age"]
    if occursin(r"[0-9]+", age)
        patient_age = parse(Int64, match(r"[0-9]+", age).match)
    else
        patient_age = -1
    end
    duration = children(children(children(children(children(eef)[end])[1])[1])[4])[2]["value"]
    duration = 3600 * parse.(Float64, match(r"PT([0-9]+)H([0-9]+)M([0-9]+\.[0-9]+)S", duration).captures[1]) + 60 * parse.(Float64, match(r"PT([0-9]+)H([0-9]+)M([0-9]+\.[0-9]+)S", duration).captures[2]) + parse.(Float64, match(r"PT([0-9]+)H([0-9]+)M([0-9]+\.[0-9]+)S", duration).captures[3])
    recording = children(children(children(children(eef)[end])[1])[2])[1][1]["name"]

    eeg_file = "exam-\$ver$id.1"
    @assert isfile(eeg_file) "File $eeg_file cannot be loaded."

    file_type = ""

    fid = nothing
    try
        fid = open(eeg_file, "r")
    catch
        error("File $eeg_file cannot be loaded.")
    end

    header = zeros(UInt8, 72)
    readbytes!(fid, header, 72)
    ch_n = ntoh.(reinterpret(Int8, header))[end]
    buf = zeros(UInt8, 6)
    readbytes!(fid, buf, 6)

    labels = String[]
    for idx in 1:ch_n
        buf = zeros(UInt8, 26)
        readbytes!(fid, buf, 26)
        strip(String(Char.(buf[1:4])))
        push!(labels, replace(strip(String(Char.(buf[1:4]))), '\0'=>""))
    end

end