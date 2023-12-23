export import_dat

"""
    import_dat(file_name)

Load Neuroscan DAT file.

# data_format = :i16
# data_format = :i32
# Arguments

- `file_name::String`: name of the file to load

# Returns

- `dat::DataFrame`
"""
function import_dat(file_name)

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert splitext(file_name)[2] == ".dat" "This is not DAT file."

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        @error "File $file_name cannot be loaded."
    end
    
    buffer = readlines(fid)[21:end]

    m = zeros(Int64, length(buffer), 5)
    for idx in 1:length(buffer)
        m[idx, 1] = parse(Int64, strip(buffer[idx][1:5]))
        m[idx, 2] = parse(Int64, strip(buffer[idx][6:12]))
        m[idx, 3] = parse(Int64, strip(buffer[idx][13:19]))
        m[idx, 4] = parse(Int64, strip(buffer[idx][20:25]))
        m[idx, 5] = parse(Int64, strip(buffer[idx][26:end]))
    end

    return DataFrame(:event=>1:size(m, 1), :trial=>m[:, 1], :response=>m[:, 2], :type=>m[:, 3], :correct=>m[:, 4], :latency=>m[:, 5])

end
