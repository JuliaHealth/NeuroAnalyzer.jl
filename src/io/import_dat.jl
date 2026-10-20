export import_dat

"""
    import_dat(file_name)

Load a Neuroscan (Stim2 / Gentask) DAT behavioural data file and return it as a `DataFrame`.

The DAT format has a text header (`key.....= value` lines), a row of column names (e.g. `Trial  Resp  Type  Correct  Latency`), a row of dashes, and space-separated data rows (one per trial). The column names are taken from the file and converted to lowercase; `Resp` is renamed to `response`. Latencies are in ms, as stored.

# Arguments

- `file_name::String`: path to the `.dat` file

# Returns

- `DataFrame`: one row per trial, e.g. with columns `:trial`, `:response`, `:type`, `:correct`, `:latency`
"""
function import_dat(file_name)::DataFrame
    # validate
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) == ".dat" ||
        throw(ArgumentError("$file_name is not a DAT file."))

    lines = rstrip.(readlines(file_name), '\r')
    # the column names are in the line above the row of dashes
    sep = findfirst(l -> occursin(r"^\s*-+(\s+-+)*\s*$", l), lines)
    (isnothing(sep) || sep == 1) &&
        throw(ArgumentError("$file_name: column header (names and a row of dashes) not found."))
    cols = Symbol.(replace.(lowercase.(split(strip(lines[sep - 1]))), "resp" => "response"))

    rows = [split(strip(l)) for l in lines[(sep + 1):end] if !isempty(strip(l))]
    all(r -> length(r) == length(cols), rows) || throw(
        ArgumentError("$file_name: data rows do not have $(length(cols)) columns ($(Base.join(cols, ", "))).")
    )
    dat = DataFrame()
    for (idx, c) in enumerate(cols)
        v = [r[idx] for r in rows]
        # integer columns stay integer, decimal columns become Float64, anything else stays text
        dat[!, c] = if all(x -> !isnothing(tryparse(Int64, x)), v)
            parse.(Int64, v)
        elseif all(x -> !isnothing(tryparse(Float64, x)), v)
            parse.(Float64, v)
        else
            String.(v)
        end
    end

    return dat
end
