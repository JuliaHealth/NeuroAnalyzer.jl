export plot_efield2d

"""
    plot_efield2d(; <keyword arguments>)

Plot 2-dimensional electric field.

# Arguments

  - `q::Vector{Int64}`: charges values
  - `qx::Vector{Float64}`: charges x-axis positions
  - `qy::Vector{Float64}`: charges y-axis positions
  - `d::Int64=2`: density of field vectors

# Returns

  - `p::GLMakie.Figure`
"""
function plot_efield2d(q::Vector{Int64}, qx::Vector{Float64}, qy::Vector{Float64}, d::Int64 = 2)::GLMakie.Figure

    # m, n = size(norm_e)
    # x = round.(collect(range(-1, 1, m)), digits=3)
    # y = round.(collect(range(-1, 1, n)), digits=3)
    # nq = length(q)

    # c = floor.(Int64, NeuroAnalyzer.normalize_n(norm_e, 15)) .+ 1
    # cl = range(colorant"lightgreen", stop=colorant"yellow", length=16)

    qs = []
    for idx in eachindex(q)
        push!(qs, (q[idx], qx[idx], qy[idx]))
    end

    function E(q, rx, ry, x, y)
        d = sqrt((x-rx)^2 + (y-ry)^2)^3
        return (q * (x - rx) / d, q * (y - ry) / d)
    end

    function fieldE(x, y)
        Ex, Ey = 0, 0
        for q in qs
            ex, ey = E(q..., x, y)
            Ex += ex
            Ey += ey
        end
        return Point(Ex, Ey)
    end

    # prepare plot
    GLMakie.activate!(title = "plot_efield()")
    plot_size = (800, 800)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1];
        aspect = DataAspect(),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # draw field lines
    streamplot!(ax, fieldE, -5..5, -5..5; arrow_size = 10, linewidth = 1, colorrange = (-3, 3), colormap = :bluesreds)

    # draw charges
    for idx in eachindex(qs)
        GLMakie.scatter!(
            ax, Point(qs[idx][2:3]); markersize = abs(qs[idx][1]) * 20, color = qs[idx][1] > 0 ? :red : :blue
        )
    end

    return p

end
