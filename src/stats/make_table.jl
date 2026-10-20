export make_table

"""
    make_table(; <keyword arguments>)

Display data as a formatted table using `PrettyTables.jl`.

The header row is prepended to the data and printed as a single table body with a horizontal separator after the first row.

# Arguments

- `header::Matrix{String}`: single-row header matrix, e.g. `["Group" "A" "B"]`; must have the same number of columns as `data`
- `data::Matrix{Any}`: data matrix, e.g. `["var1" 1.0 2.0; "var2" 3.0 4.0]`; integer values are shown without decimals
- `digits::Int64=3`: number of decimal places shown for non-integer numbers

# Returns

- `Nothing`

# Notes

Values are rounded for display only; `data` is not modified.
"""
function make_table(; header::Matrix{String}, data::Matrix{Any}, digits::Int64 = 3)::Nothing
    # validate
    size(header, 1) == 1 || throw(ArgumentError("header must be a single-row matrix."))
    size(header, 2) == size(data, 2) ||
        throw(ArgumentError("header and data must have the same number of columns."))

    digits >= 0 || throw(ArgumentError("digits must be ≥ 0."))

    # integers are shown as they are (converted on a copy, the caller's data are kept)
    cells = Any[x isa Integer ? string(x) : x for x in data]

    # PrettyTables ≥ 3: the header is shown as column labels, the first (label) column left-aligned
    pretty_table(
        cells;
        column_labels = vec(header),
        alignment = vcat(:l, fill(:r, size(cells, 2) - 1)),
        # applies to numeric cells only
        formatters = [fmt__printf("%1.$(digits)f")],
    )

    return nothing
end
