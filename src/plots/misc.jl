export plot_compose
export plot_empty
export add_locs

"""
    plot_compose(p; <keyword arguments>)

Compose a complex plot of various plots contained in vector `p` using layout `layout`. Layout scheme is:
- `(2, 2)`: 2 × 2 plots, regular layout
- `grid(4, 1, heights=[0.6, 0.1, 0.1, 0.1]`: 4 × 1 plots, irregular layout
- `@layout [a{0.2w} b{0.8w};_ c{0.6}]`: complex layout using Plots.jl `@layout` macro

# Arguments

- `p::Vector{Plots.Plot{Plots.GRBackend}}`: vector of plots
- `layout::Union(Matrix{Any}, Tuple{Int64, Int64}, Plots.GridLayout}`: layout
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for `p` vector plots

# Returns

- `pc::Plots.Plot{Plots.GRBackend}`
"""
function plot_compose(p::Vector{Plots.Plot{Plots.GRBackend}}; layout::Union{Matrix{Any}, Tuple{Int64, Int64}, Plots.GridLayout}, mono::Bool=false, kwargs...)

    pal = mono ? :grays : :darktest
    if typeof(layout) == Tuple{Int64, Int64} && length(p) < layout[1] * layout[2]
        for _ in 1:(layout[1] * layout[2]) - length(p)
            push!(p, plot_empty())
        end
    end

    pc = plot_empty()
    pc = Plots.plot!(p..., layout=layout, palette=pal; kwargs...)
    Plots.plot(pc)

    return pc

end

"""
    plot_empty()

Return an empty plot, useful for filling matrices of plots. 

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_empty()

    return Plots.plot(grid=false,
                      framestyle=:none,
                      border=:none,
                      margins=0Plots.px)
    
end

"""
    add_locs(p1, p2; view, file_name)

Add locations to a plot. Locations are placed in the top right corner. If `file_name` is provided, the plot is saved as PNG file.

# Arguments

- `p1::Plots.Plot{Plots.GRBackend}`: signal plot
- `p2::Plots.Plot{Plots.GRBackend}`: locations plot
- `view::Bool=true`: view the output image
- `file_name::String=""`: output image filename

# Returns

- `c::Cairo.CairoSurfaceBase{UInt32}`
"""
function add_locs(p1::Plots.Plot{Plots.GRBackend}, p2::Plots.Plot{Plots.GRBackend}; view::Bool=true, file_name::String="")

    if file_name != ""
        ext = lowercase(splitext(file_name)[2])
        @assert ext == ".png" "Filename extension must be .png"

        (isfile(file_name) && verbose) && _warn("File $file_name will be overwritten.")
    end

    p1_size = p1.attr[:size]
    p2_size = p2.attr[:size]
    c = CairoRGBSurface(p1_size[1], p1_size[2])
    cr = CairoContext(c)
    show(io, MIME("image/png"), p1)
    img = read_from_png(io)
    Cairo.set_source_surface(cr, img, 0, 0)
    Cairo.paint(cr)
    Cairo.scale(cr, 0.5, 0.5)
    show(io, MIME("image/png"), p2)
    img = read_from_png(io)
    Cairo.set_source_surface(cr, img, (2 * p1_size[1]) - p2_size[1], 0)
    Cairo.paint(cr)

    if file_name != ""
        Cairo.write_to_png(c, file_name)
    else
        view && iview_plot(c)
    end

    return c

end