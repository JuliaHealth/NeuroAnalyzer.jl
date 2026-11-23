export import_markers

"""
    import_markers(file_name; <keyword arguments>)

Import markers from a CSV or TSV file and return a markers `DataFrame`, which can be added to an object with `add_markers()`.

Two layouts are recognized (the delimiter, comma, tab or semicolon, is detected from the file; `n/a` and empty cells are missing values):

- markers exported by `export_markers()`: columns `type`, `start`, `length`, `value`, `channel` (in any order; `id` instead of `type` in files of NeuroAnalyzer before v0.26.12), other columns are kept as extra marker columns; `channel` holds channel labels separated by `;` (empty: all channels)
- BIDS `events.tsv`: `onset` and `duration` in seconds (or `sample`, 0-based, if there is no `onset`); a `trial_type` that names a marker type (e.g. `stimulus`, `response`, `BAD_…`) gives the marker type and the value is the event code (`value`); a condition (e.g. `face`) gives the type `stim` and the value is the condition; without `trial_type` the type is `stim` and the value is the event code; all other columns are kept as extra marker columns (`value` as `bids_value`)

# Arguments

- `file_name::String`: path to the CSV/TSV file
- `fs::Union{Nothing, Real}=nothing`: sampling rate in Hz; required for BIDS files with `sample` but no `onset`

# Returns

- `DataFrame`: markers (`type`, `start`, `length`, `value`, `channel` and extra columns), sorted by `start`

# Notes

- The columns `type`, `value`, `channel` and `trial_type` are read as text (e.g. the code `001` keeps its zeros); other columns are parsed (numbers, text).
- The marker channels are checked against the channels of the object by `add_markers()`; channel numbers (`0`: all channels) of old files are also converted by `add_markers()`.
"""
function import_markers(
    file_name::String;
    fs::Union{Nothing, Real} = nothing,
)::DataFrame
    # validate
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    isnothing(fs) || fs > 0 || throw(ArgumentError("fs must be > 0."))

    # codes, labels and conditions are read as text (e.g. "001" or channel labels "1", "2")
    txt = (:type, :id, :value, :channel, :trial_type)
    df = CSV.read(
        file_name,
        DataFrame;
        stringtype = String,
        missingstring = ["n/a", ""],
        types = (_, nm) -> Symbol(nm) in txt ? String : nothing,
    )
    cols = names(df)

    # BIDS events.tsv
    if "onset" in cols || ("sample" in cols && !("start" in cols))
        if !("onset" in cols)
            isnothing(fs) && throw(
                ArgumentError("fs must be given for events with sample numbers but no onset."),
            )
        end
        "onset" in cols && any(ismissing, df[!, :onset]) &&
            throw(ArgumentError("Column onset contains missing values."))
        return sort!(_bids_events(df, isnothing(fs) ? 1.0 : fs), :start)
    end

    # markers exported by export_markers()
    old_id = "id" in cols && !("type" in cols)
    old_id && DataFrames.rename!(df, :id => :type)
    cols = names(df)
    all(in(cols), _MARKER_FIXED) || throw(
        ArgumentError(
            "$file_name must contain the columns $(_MARKER_FIXED) (markers exported by export_markers()) or onset/sample (BIDS events.tsv); got $cols.",
        ),
    )
    any(ismissing, df[!, :start]) && throw(ArgumentError("Column start contains missing values."))
    # fixed columns first, extra columns after them
    df = DataFrames.select(df, _MARKER_FIXED, :)

    value = String[ismissing(v) ? "" : string(v) for v in df[!, :value]]
    df[!, :type] = if old_id
        String[_marker_type_old(ismissing(t) ? "" : string(t), v) for (t, v) in zip(df[!, :type], value)]
    else
        String[ismissing(t) ? "" : string(t) for t in df[!, :type]]
    end
    df[!, :value] = value
    df[!, :start] = Float64.(df[!, :start])
    df[!, :length] = Float64[ismissing(l) ? 0.0 : l for l in df[!, :length]]
    # channel labels separated by ";"; files of NeuroAnalyzer before v0.26.12 store channel
    # numbers (0: all channels), converted to labels by add_markers()
    if old_id
        df[!, :channel] = Int64[ismissing(c) ? 0 : parse(Int64, c) for c in df[!, :channel]]
    else
        df[!, :channel] = Vector{String}[
            ismissing(c) ? String[] :
            String[strip(x) for x in split(string(c), ';') if !isempty(strip(x))] for
            c in df[!, :channel]
        ]
    end

    return sort!(df, :start)
end
