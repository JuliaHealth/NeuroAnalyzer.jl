export import_dat

"""
    import_dat(file_name)

Load Neuroscan DAT file.

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `dat::DataFrame`
"""
function import_dat(file_name)

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) == ".dat" "This is not DAT file."

    dat = CSV.read(file_name, stringtype=String, delim=' ', ignorerepeated=true, skipto=21, header=0, DataFrame)
    DataFrames.rename!(dat, [:event, :trial, :response, :type, :correct])
    dat[!, :event] .-= 1

    return dat

end
