export plot_save

"""
    plot_save(fig; <keyword arguments>)

Saves plot as file (PNG/PDF).

# Arguments

- `fig::Union{Plots.Plot{Plots.GRBackend}, Makie.Figure}`: figure or plot object to save
- `file_name::String`: output file path and name; file format is determined by the extension:
    - `.png`: Portable Network Graphics (lossless)
    - `.pdf`: Portable Document Format (vector graphics)

# Returns

- `Nothing`
"""
function plot_save(
    fig::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure};
    file_name::String,
)::Nothing
    # validate
    ext = splitext(file_name)[2]
    ext in (".png", ".pdf") ||
        throw(ArgumentError("Unsupported file format \"$ext\"; must be .png or .pdf."))
    isfile(file_name) && _warn("File $file_name will be overwritten.")

    try
        if fig isa Plots.Plot{Plots.GRBackend}
            Plots.savefig(fig, file_name)
        else
            GLMakie.save(file_name, fig)
        end
    catch err
        @error "File $file_name cannot be written." exception=err
    end

    return nothing
end
