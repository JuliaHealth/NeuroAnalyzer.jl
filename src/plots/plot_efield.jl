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

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_efield2d(q::Vector{Int64}, qq::Vector{Vector{Float64}}, norm_e::Matrix{Float64}, ex::Matrix{Float64}, ey::Matrix{Float64}, d::Int64=2)::Plots.Plot{Plots.GRBackend}

    m, n = size(norm_e)
    x = round.(collect(range(-1, 1, m)), digits=3)
    y = round.(collect(range(-1, 1, n)), digits=3)
    nq = length(q)

    c = floor.(Int64, NeuroAnalyzer.normalize_n(norm_e, 15)) .+ 1
    cl = range(colorant"lightgreen", stop=colorant"yellow", length=16)

    p = Plots.plot(legend=false,
                   xlims=(-1.2, 1.2),
                   ylims=(-1.2, 1.2),
                   ratio=1,
                   tick_direction=:out,
                   size=(800, 800))

    for idx1 in 1:d:m
        for idx2 in 1:d:n
            p = Plots.plot!([x[idx1], x[idx1] + ex[idx2, idx1] / 30],
                            [y[idx2], y[idx2] + ey[idx2, idx1] / 30],
                            lc=cl[c[idx2, idx1]],
                            lw=0.75,
                            la=0.75,
                            arrow=arrow(:closed))
        end
    end

    qs = sqrt.(abs.(q))
    for idx in 1:nq
        p = Plots.scatter!((qq[idx][2], qq[idx][1]),
                           ms=qs[idx],
                           color=q[idx] > 0 ? :red : :blue,
                           markerstrokewidth=0,
                           markerstrokealpha=0)
    end

    Plots.plot!(p)

    return p

end
