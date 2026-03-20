export import_montage

"""
    import_montage(file_name)

Load a montage definition from a plain-text file.

The file format is:

- **Line 1**: montage name (e.g. `longitudinal-BIP`)
- **Subsequent lines**: one channel or channel pair per line (e.g. `Fp1`, `Fz-Cz`)

Blank lines and lines beginning with `#` are ignored, so comments and trailing newlines do not corrupt the channel list.

Example montage files are located in the `montages/` folder.

Each channel/channel pair must be in a separate line

# Arguments

- `file_name::String`: path to the montage text file

# Returns

Named tuple:

- `ref_list::Vector{String}`: list of channel / channel-pair strings
- `ref_name::String`: name of the montage

# Throws

- `ArgumentError` if the file does not exist or does not contain a valid montage structure (at least a name line and one channel entry)
"""
function import_montage(
        file_name::String
    )::@NamedTuple{ref_list::Vector{String}, ref_name::String}

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    montage_file = open(readlines, file_name)

    Base.filter!(l -> !isempty(strip(l)) && !startswith(l, '#'), montage_file)

    length(montage_file) >= 2 ||
        throw(ArgumentError(
            "$file_name does not contain a valid montage structure " *
            "(expected: name on line 1, at least one channel entry)."))

    return (ref_list = montage_file[2:end], ref_name = montage_file[1])

end
