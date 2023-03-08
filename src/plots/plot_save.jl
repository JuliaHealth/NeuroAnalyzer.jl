export plot_save

"""
    plot_save(p; file_name::String)

Saves plot as file (PDF/PNG/TIFF). File format is determined using `file_name` extension.

# Arguments

- `p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}`
- `file_name::String`
"""
function plot_save(p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}; file_name::String)

    ext = splitext(file_name)[2]
    _check_var(ext, [".png", ".pdf", ".jpg", ".tiff"], "File format")
    (isfile(file_name) && verbose == true) && _info("File $file_name will be overwritten.")
    if typeof(p) == Plots.Plot{Plots.GRBackend}
        savefig(p, file_name)
    else
        save(file_name, p)
    end

    nothing
end
