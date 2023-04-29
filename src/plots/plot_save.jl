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
    
    if p isa Plots.Plot{Plots.GRBackend}
        savefig(p, file_name)
    else
        FileIO.save(file_name, p)
    end

    return nothing

end
