export view_header

"""
    header(obj)

Print all keys and values stored in the object header.

Iterates over the `:subject`, `:recording`, and `:experiment` sub-dictionaries of `obj.header` and prints each key–value pair in the format `header.<section>[:<key>]: <value>`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
    
# Returns

- `Nothing`

# See also

[`header`](@ref), [`info`](@ref)
"""
function view_header(obj::NeuroAnalyzer.NEURO)::Nothing

    # print the top-level header field names (struct fields, not dict keys)
    fields = join(fieldnames(typeof(obj.header)), ", ")
    println("Header fields: $fields")

    for section in (:subject, :recording, :experiment)
        for (key, value) in getfield(obj.header, section)
            println("header.$section[:$key]: $value")
        end
    end

    return nothing

end
