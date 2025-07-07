export plot_efield2d

function plot_efield2d(norm_e::matrix{Float64}, ex::matrix{Float64}, ey::matrix{Float64})

    m, n = size(norm_e)

    ex ./= round(m * 0.25, digits=2)
    ey ./= round(m * 0.25, digits=2)

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
            p = Plots.plot!([x[idx1], x[idx1] + ex[idx2, idx1]],
                            [y[idx2], y[idx2] + ey[idx2, idx1]],
                            lc=cl[c[idx2, idx1]],
                            lw=0.75,
                            la=0.75,
                            arrow=arrow(:closed))
        end
    end

    for idx in 1:nq
        p = Plots.scatter!((qq[idx][2], qq[idx][1]),
                           ms=5 * abs(q[idx]),
                           color=q[idx] > 0 ? :red : :blue,
                           markerstrokewidth=0,
                           markerstrokealpha=0)
    end

    Plots.plot!(p)

    return p

end
