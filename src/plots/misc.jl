export plot_compose
export plot_empty
export add_plot_locs

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

- `pc::GLMakie.Figure`
"""
function plot_compose(p::Vector{GLMakie.Figure}; layout::Union{Tuple{Int64, Int64}}, mono::Bool=false, kwargs...)::GLMakie.Figure

    @assert layout[1] * layout[2] >= length(p) "Layout size ($(layout[1]) × $(layout[2])) must be ≥ the number of plots ($(length(p)))."

    s =(0, 0)
    for idx in eachindex(p)
        size(p[idx].scene) > s && (s = size(p[idx].scene))
    end
    for idx in eachindex(p)
        size(p[idx].scene) != s && _warn("For best results all plots should have the size of $(s[1])×$(s[2]).")
    end

    layout[1] == layout[2] && (s = (s[1] * layout[1] * 0.75, s[2] * layout[2] * 0.75))
    layout[1] > layout[2] && (s = (s[1] * layout[1] * 0.5, s[2] * layout[2] * 1.5))
    layout[1] < layout[2] && (s = (s[1] * layout[1] * 1.25, s[2] * layout[2] * 0.5))

    if length(p) < layout[1] * layout[2]
        for _ in 1:(layout[1] * layout[2]) - length(p)
            push!(p, plot_empty())
        end
    end

    pc = GLMakie.Figure(size=s)
    gl = pc[1, 1] = GridLayout(layout[1], layout[2])
    p_idx = 1
    for idx1 in 1:layout[1]
        for idx2 in 1:layout[2]
            fname = tempname()*".png"
            GLMakie.save(fname, p[p_idx])
            pp = FileIO.load(fname)
            rm(fname)
            ax = GLMakie.Axis(pc[idx1, idx2],
                              topspinevisible=false,
                              bottomspinevisible=false,
                              leftspinevisible=false,
                              rightspinevisible=false,
                              aspect=DataAspect())
            GLMakie.image!(ax, rotr90(pp))
            hidedecorations!(ax)
            p_idx += 1
        end
    end

    return pc

end

"""
    plot_empty()

Return an empty plot, useful for filling matrices of plots.

# Returns

- `p::GLMakie.Figure`
"""
function plot_empty()::GLMakie.Figure

    p = GLMakie.Figure()

    return p

end

"""
    add_plot_locs(p1, p2; <keyword arguments>)

Add locations to a plot. Locations are placed in the top right corner. If `file_name` is provided, the plot is saved as PNG file.

# Arguments

- `p1::Plots.Plot{Plots.GRBackend}`: signal plot
- `p2::Plots.Plot{Plots.GRBackend}`: locations plot
- `view::Bool=true`: view the output image
- `file_name::String=""`: output image filename

# Returns

- `c::Cairo.CairoSurfaceBase{UInt32}`
"""
function add_plot_locs(p1::Plots.Plot{Plots.GRBackend}, p2::Plots.Plot{Plots.GRBackend}; view::Bool=true, file_name::String="")::Cairo.CairoSurfaceBase{UInt32}

    p1_size = p1.attr[:size]
    p2_size = p2.attr[:size]
    c = CairoRGBSurface(p1_size[1], p1_size[2])
    cr = CairoContext(c)
    withenv("GKSwstype" => "100") do
        png(p1, io)
    end
    img = read_from_png(io)
    Cairo.set_source_surface(cr, img, 0, 0)
    Cairo.paint(cr)
    Cairo.scale(cr, 0.5, 0.5)
    withenv("GKSwstype" => "100") do
        png(p2, io)
    end
    img = read_from_png(io)
    Cairo.set_source_surface(cr, img, (2 * p1_size[1]) - p2_size[1], 0)
    Cairo.paint(cr)

    if file_name != ""
        ext = lowercase(splitext(file_name)[2])
        @assert ext == ".png" "Filename extension must be .png"
        (isfile(file_name) && verbose) && _warn("File $file_name will be overwritten.")
        Cairo.write_to_png(c, file_name)
    else
        view && iview(c)
    end

    return c

end
