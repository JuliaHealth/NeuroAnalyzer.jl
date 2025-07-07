export efield2d

function efield2d(; q::Vector{Int64}, qx::Vector{Int64}, qy::Vector{Int64})

    _wip()
    m = 100
    n = 100

    @assert length(qx) == length(q) "Length of qx and number of charges must be equal."
    @assert length(qx) == length(q) "Length of qy and number of charges must be equal."

    # number of charges
    nq = length(q)
    for idx in eachindex(q)
        NeuroAnalyzer._in(qx[idx], (1, m))
        NeuroAnalyzer._in(qy[idx], (1, n))
    end

    x = round.(collect(range(-1, 1, m)), digits=3)
    y = round.(collect(range(-1, 1, n)), digits=3)

    #strength
    ex = zeros(n, m)
    ey = zeros(n, m)

    qq = Vector{Vector{Float64}}()

    @inbounds for idx in 1:nq
        push!(qq, [x[qy[idx]], y[qx[idx]]])
        for idx1 in 1:n
            for idx2 in 1:m
                denom = ((idx1 - qy[idx])^2 + (idx2 - qx[idx])^2)^1.5
                if denom != 0 
                    ex[idx1, idx2] += q[idx] * (idx2 - qx[idx]) / denom
                    ey[idx1, idx2] += q[idx] * (idx1 - qy[idx]) / denom
                end
            end
        end
    end

    # normalized values for arrows to be of equal length
    # norm of E field matrix
    norm_e = hypot.(ex, ey)
    ex = ex ./ norm_e
    ey = ey ./ norm_e

    return (qq=qq, norm_e=norm_e, x=x, y=y, ex=ex, ey=ey)

end
