export plot_compose
export plot_empty

"""
    plot_compose(p; <keyword arguments>)

Compose a complex plot of various plots contained in vector `p` using layout `layout`. Layout scheme is:
- `(2, 2)`: 2 × 2 plots, regular layout
- `grid(4, 1, heights=[0.6, 0.1, 0.1, 0.1]`: 4 × 1 plots, irregular layout
- `@layout [a{0.2w} b{0.8w};_ c{0.6}]`: complex layout using Plots.jl `@layout` macro

# Arguments

- `p::Vector{Plots.Plot{Plots.GRBackend}}`: vector of plots
- `layout::Union(Matrix{Any}, Tuple{Int64, Int64}, Plots.GridLayout}`: layout
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for `p` vector plots

# Returns

- `pc::Plots.Plot{Plots.GRBackend}`
"""
function plot_compose(p::Vector{Plots.Plot{Plots.GRBackend}}; layout::Union{Matrix{Any}, Tuple{Int64, Int64}, Plots.GridLayout}, mono::Bool=false, kwargs...)

    pal = mono ? :grays : :darktest
    if typeof(layout) == Tuple{Int64, Int64} && length(p) < layout[1] * layout[2]
        for _ in 1:(layout[1] * layout[2]) - length(p)
            push!(p, Plots.plot(border=:none, title=""))
        end
    end

    pc = Plots.plot(grid=false,
                    framestyle=:none,
                    border=:none,
                    margins=0Plots.px)
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

    return Plots.plot(grid=false, border=:none, title="")
    
end
