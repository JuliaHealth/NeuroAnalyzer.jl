export make_table

"""
    make_table(; <keyword arguments>)

Display data as a formatted table using `PrettyTables.jl`.

The header row is prepended to the data and printed as a single table body with a horizontal separator after the first row.

# Arguments

- `header::Matrix{String}`: single-row header matrix, e.g. `["Group" "A" "B"]`; must have the same number of columns as `data`
- `data::Matrix{Any}`: data matrix, e.g. `["var1" 1.0 2.0; "var2" 3.0 4.0]`; integer values are converted to strings before display

# Returns

- `Nothing`

# Throws

- `ArgumentError`: if `size(header, 1) ≠ 1` or `size(header, 2) ≠ size(data, 2)`

# Notes

- The `formatters = ft_printf("%1.3f", 2:3)` format is hardcoded to columns 2 and 3; callers with a different column count should adjust the source accordingly.
- Integer elements in `data` are converted to `String` in-place to ensure uniform display.
"""
function make_table(; header::Matrix{String}, data::Matrix{Any})::Nothing

    !(size(header, 1) == 1) && throw(ArgumentError("header must be a single-row matrix."))
    !(size(header, 2) == size(data, 2)) && throw(ArgumentError("header and data must have the same number of columns."))

    # convert any Integer cells to String to avoid type-display inconsistencies;
    # iterate with CartesianIndices to cover all dimensions safely
    for idx in CartesianIndices(data)
        data[idx] isa Integer && (data[idx] = string(data[idx]))
    end

    pretty_table(
        cat(header, data, dims=1);
        body_hlines = [1],
        cell_alignment = Dict((1, 1) => :l),
        formatters = ft_printf("%1.3f", 2:3),
        show_header = false
    )

    return nothing

end
