export import_dat

"""
    import_dat(file_name)

Load a Neuroscan DAT behavioural data file and return it as a `DataFrame`.

The DAT format has a 20-line text header followed by space-separated data rows with five columns: event index, trial number, response, type, and correctness flag. Event indices are stored as 1-based in the file and are converted to 0-based on import to match Neuroscan's convention.

# Arguments

- `file_name::String`: path to the `.dat` file

# Returns
- `DataFrame`: table with columns `:event`, `:trial`, `:response`, `:type`, `:correct`

# Throws
- `ArgumentError` if the file does not exist, is not a `.dat` file, or does not contain exactly 5 data columns
"""
function import_dat(file_name)::DataFrame

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) == ".dat" ||
        throw(ArgumentError("$file_name is not a DAT file."))

    # read space-delimited data starting after the 20-line text header
    # `ignorerepeated = true` collapses multiple consecutive spaces
    dat = CSV.read(
        file_name;
        stringtype = String,
        delim = ' ',
        ignorerepeated = true,
        skipto = 21, # 20-line header; data starts at line 21
        header = 0, # no column-header row in the data block
        DataFrame,
    )

    expected_cols = [:event, :trial, :response, :type, :correct]
    DataFrames.ncol(dat) == length(expected_cols) ||
        throw(ArgumentError(
            "$file_name has $(DataFrames.ncol(dat)) data columns; " *
            "expected $(length(expected_cols)) " *
            "($(join(expected_cols, ", ")))."))
    DataFrames.rename!(dat, expected_cols)

    # convert event index from 1-based (file convention) to 0-based
    dat[!, :event] = Int64.(dat[!, :event]) .- 1

    return dat

end
