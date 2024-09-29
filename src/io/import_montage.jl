export import_montage

"""
    import_montage(file_name)

Load montage from a text file. Example montage files are located in the `montages/` folder. The structure of the file is:

- first line: name of the montage, e.g. `longitudinal-BIP`
- next lines: channel pairs or individual channels, e.g. `Fz-Cz` or `Fp1`

Each channel/channel pair must be in a separate line

# Arguments

- `file_name::String`: name of the file to load

# Returns

Named tuple containing:
- `ref_list::Vector{String}`: list of channel pairs
- `ref_name::String`: name of the montage
"""
function import_montage(file_name::String)::NamedTuple{ref_list::Vector{String}, ref_name::String}

    @assert isfile(file_name) "File $file_name cannot be loaded."

    f = open(file_name, "r")
    montage_file = readlines(f)
    close(f)

    @assert length(montage_file) >= 2 "File does not contain proper montage structure."
    ref_name = montage_file[1]
    ref_list = montage_file[2:end]

    return (ref_list=ref_list, ref_name=ref_name)

end
