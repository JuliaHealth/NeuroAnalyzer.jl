export plot_efield2d

"""
    plot_efield2d(; <keyword arguments>)

Plot 2-dimensional electric field.

# Arguments

- `q::Vector{Int64}`: charges values, qy::Vector{Int64}::String`: anode location
- `qq::Vector{Vector{Float64}}`: charges positions
- `norm_e::Matrix{Float64}`: normalized electric field
- `ex::Matrix{Float64}`: electric field X axis vector
- `ey::Matrix{Float64}`: electric field Y axis vector
- `d::Int64=2`: density of field vectors

# Returns

- `p::GLMakie.Figure`
"""
function plot_efield2d(q::Vector{Int64}, qq::Vector{Vector{Float64}}, norm_e::Matrix{Float64}, ex::Matrix{Float64}, ey::Matrix{Float64}, d::Int64=2)::GLMakie.Figure

    m, n = size(norm_e)
    x = round.(collect(range(-1, 1, m)), digits=3)
    y = round.(collect(range(-1, 1, n)), digits=3)
    nq = length(q)

    c = floor.(Int64, NeuroAnalyzer.normalize_n(norm_e, 15)) .+ 1
    cl = range(colorant"lightgreen", stop=colorant"yellow", length=16)

    function E(q, rx, ry, x, y)
        d = sqrt((x-rx)^2 + (y-ry)^2)^3
        return (q * (x - rx) / d, q * (y - ry) / d)
    end
    function charges(; nq = 2)
        qs = []
        for i in 1:nq
            q = i % 2 * 2 - 1
            push!(qs, (q, cos(2π * i / nq), sin(2π * i / nq)))
        end
        qs
    end
    function fieldE(x,y)
        Ex, Ey = 0, 0
        for q in qs
            ex, ey = E(q..., x, y)
            Ex += ex
            Ey += ey
        end
        Point(Ex, Ey)
    end

    # prepare plot
    plot_size = (800, 800)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      aspect=DataAspect(),
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0),
                      xzoomlock=true,
                      yzoomlock=true,
                      xpanlock=true,
                      ypanlock=true,
                      xrectzoom=false,
                      yrectzoom=false)
    GLMakie.xlims!(ax, (-1.2, 1.2))
    GLMakie.ylims!(ax, (-1.2, 1.2))
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    for idx1 in 1:d:m
        for idx2 in 1:d:n
            streamplot!(ax,
                        fieldE,
                        -2..2, -2..2;
                        arrow_size = 6,
                        linewidth = 0.5,
                        colorrange = (-3,3),
                        colormap=:bluesreds)
#            p = Plots.plot!([x[idx1], x[idx1] + ex[idx2, idx1] / 30],
#                            [y[idx2], y[idx2] + ey[idx2, idx1] / 30],
#                            lc=cl[c[idx2, idx1]],
#                            lw=0.75,
#                            la=0.75,
#                            arrow=arrow(:closed))
        end
    end

    qs = sqrt.(abs.(q))
    for idx in 1:nq
        p = GLMakie.scatter!(ax,
                             qq[idx][2],
                             qq[idx][1],
                             ms=qs[idx],
                             color=q[idx] > 0 ? :red : :blue,
                             strokewidth=0)
    end

    return p

end
