export view_header

"""
    view_header(obj)

Print all keys and values stored in the object header.

Iterates over the `:subject` and `:recording` sub-dictionaries of `obj.header` and prints each key–value pair in the format `header.<section>[:<key>]: <value>`, followed by the header notes (`header.notes`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
   
# Returns

- `Nothing`
"""
function view_header(obj::NeuroAnalyzer.NEURO)::Nothing
    # print the top-level header field names (struct fields, not dict keys)
    fields = Base.join(fieldnames(typeof(obj.header)), ", ")
    println("Header fields: $fields")

    for section in (:subject, :recording)
        for (key, value) in getfield(obj.header, section)
            println("header.$section[:$key]: $value")
        end
    end
    println("header.notes: $(obj.header.notes)")

    return nothing
end
