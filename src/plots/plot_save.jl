export plot_save

"""
    plot_save(fig; <keyword arguments>)

Saves plot as file (PNG/PDF). File format is determined using `file_name` extension.

# Arguments

- `fig::Union{Plots.Plot{Plots.GRBackend}, Makie.Figure}`
- `file_name::String`

# Returns

- `Nothing`
"""
function plot_save(fig::Union{Plots.Plot{Plots.GRBackend}, Makie.Figure}; file_name::String)::Nothing

    ext = splitext(file_name)[2]
    _check_var(ext, [".png", ".pdf"], "File format")

    (isfile(file_name) && verbose) && _warn("File $file_name will be overwritten.")
    try
        if fig isa Plots.Plot{Plots.GRBackend}
            Plots.savefig(fig, file_name)
        else
            GLMakie.save(file_name, fig)
        end
    catch err
        if isa(err, SystemError)
            @error "File $file_name cannot be written."
        end
    end

    return nothing

end
