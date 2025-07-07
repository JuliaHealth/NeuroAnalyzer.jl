export plot_efield2d

function plot_efield2d(q::Vector{Int64}, qq::Vector{Vector{Float64}}, norm_e::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64}, ex::Matrix{Float64}, ey::Matrix{Float64}, d::Int64=2)

    m, n = size(norm_e)
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

    for idx in 1:nq
        p = Plots.scatter!((qq[idx][2], qq[idx][1]),
                           ms=NeuroAnalyzer.normalize_n(abs.(q), 7)[idx]+3,
                           color=q[idx] > 0 ? :red : :blue,
                           markerstrokewidth=0,
                           markerstrokealpha=0)
    end

    Plots.plot!(p)

    return p

end
