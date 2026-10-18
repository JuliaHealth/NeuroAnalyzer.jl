export plot_save

"""
    plot_save(fig; <keyword arguments>)

Saves plot as file (PNG/PDF).

# Arguments

- `fig::Union{Plots.Plot{Plots.GRBackend}, Makie.Figure}`: figure or plot object to save
- `file_name::String`: output file path and name; file format is determined by the extension:
    - `.png`: Portable Network Graphics (lossless)
    - `.pdf`: Portable Document Format (vector graphics for `Plots` plots; GLMakie figures are rendered to a raster image embedded in the PDF)

# Returns

- `Nothing`
"""
function plot_save(
    fig::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure};
    file_name::String,
)::Nothing
    # validate
    ext = lowercase(splitext(file_name)[2])
    ext in (".png", ".pdf") ||
        throw(ArgumentError("Unsupported file format \"$ext\"; must be .png or .pdf."))
    isfile(file_name) && _warn("File $file_name will be overwritten.")

    if fig isa Plots.Plot{Plots.GRBackend}
        Plots.savefig(fig, file_name)
    elseif ext == ".png"
        # explicit format: FileIO detects it from the (case-sensitive) extension
        GLMakie.save(FileIO.File{FileIO.format"PNG"}(file_name), fig)
    else
        # GLMakie cannot write PDF: render to PNG and embed the image in a PDF page
        img = plot2canvas(fig)
        pdf = Cairo.CairoPDFSurface(file_name, img.width, img.height)
        cr = Cairo.CairoContext(pdf)
        Cairo.set_source_surface(cr, img, 0, 0)
        Cairo.paint(cr)
        Cairo.finish(pdf)
    end

    return nothing
end
