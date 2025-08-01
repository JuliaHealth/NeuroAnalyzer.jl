export make_table

"""
    make_table(; <keyword arguments>)

Display data as a table.

# Arguments

- `header::Matrix{String}`: table header, e.g. `header = ["Group" "A" "B"]`
- `data::Matrix{Any}`: table data, e.g. `data = ["var1" 1.0 2.0; "var2" 3.0 4.0]`

# Returns

Nothing
"""
function make_table(; header::Matrix{String}, data::Matrix{Any})::Nothing

    @assert size(header, 2) == size(data, 2) "Number of columns must be equal in the header and in the data."

    for idx in eachindex(data)
        data[idx] isa Int && (data[idx] = string(data[idx]))
    end

    pretty_table(
                 cat(header, data, dims=1);
                 body_hlines        = [1],
                 cell_alignment     = Dict((1, 1) => :l),
                 formatters         = ft_printf("%1.3f", 2:3),
                 show_header        = false
                 )

    return nothing

end
