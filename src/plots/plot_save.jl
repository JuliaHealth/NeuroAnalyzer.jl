export plot_save

"""
    plot_save(p; <keyword arguments>)

Saves plot as file (PNG/PDF). File format is determined using `file_name` extension.

# Arguments

- `p::Plots.Plot{Plots.GRBackend}`
- `file_name::String`

# Returns

Nothing
"""
function plot_save(p::Plots.Plot{Plots.GRBackend}; file_name::String)::Nothing

    ext = splitext(file_name)[2]
    _check_var(ext, [".png", ".pdf"], "File format")

    (isfile(file_name) && verbose) && _warn("File $file_name will be overwritten.")
    try
        savefig(p, file_name)
    catch err
        if isa(err, SystemError)
            @error "File $file_name cannot be written."
        end
    end

    return nothing

end
