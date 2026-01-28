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
    add_plot_locs(p, pl; <keyword arguments>)

Add locations to a plot. Locations are placed in the top right corner.

# Arguments

- `p1::GLMakie.Figure`: primary plot
- `p2::GLMakie.Figure`: locations plot

# Returns

- `p::GLMakie.Figure`
"""
function add_plot_locs(p::GLMakie.Figure, pl::GLMakie.Figure; view::Bool=true, file_name::String="")::GLMakie.Figure

        io = IOBuffer()
        show(io, MIME"image/png"(), pl)
        pp = FileIO.load(io)
        transparent_pp = map(c -> RGBA(color(c), 1.0), pp)
        transparent_pp[transparent_pp .== RGBA(1.0, 1.0, 1.0, 1.0)] .= RGBA(1.0, 1.0, 1.0, 0.0)
        transparent_pp[transparent_pp .== RGBA(0.999, 0.999, 0.999, 1.0)] .= RGBA(0.999, 0.999, 0.999, 0.0)
        transparent_pp[transparent_pp .== RGBA(0.998, 0.998, 0.998, 1.0)] .= RGBA(0.998, 0.998, 0.998, 0.0)
        # top right corner
        ax = contents(p[1, 1])[1]
        ax = ax.targetlimits[].origin .+ ax.targetlimits[].widths
        pos_x = ax[1]
        pos_y = ax[2]
        GLMakie.scatter!(p[1, 1],
                         pos_x,
                         pos_y,
                         marker_offset=size(transparent_pp) ./ -2,
                         marker=transparent_pp,
                         markersize=size(transparent_pp),
                         markerspace=:pixel)

        return p

end
