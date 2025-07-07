export efield2d

function efield2d()

    _wip()

    d = 1
    m = 100
    n = 100
    x = collect(range(-1, 1, m))
    y = collect(range(-1, 1, n))

    #strength
    ex = zeros(n, m)
    ey = zeros(n, m)

    # number of charges
    nq = 50

    q = rand(0.1:0.1:1, nq)
    q = ones(nq)
    qn = rand(Bool, nq)
    q[qn] = -q[qn]

    qq = Vector{Vector{Float64}}()

    @inbounds for idx in 1:nq
        qy = rand(1:n)
        qx = rand(1:m)
        push!(qq, [x[qy], y[qx]])
        for idx1 in 1:n
            for idx2 in 1:m
                denom = ((idx1 - qy)^2 + (idx2 - qx)^2)^1.5
                if denom != 0 
                    ex[idx1, idx2] += q[idx] * (idx2 - qx) / denom
                    ey[idx1, idx2] += q[idx] * (idx1 - qy) / denom
                end
            end
        end
    end

    # normalized values for arrows to be of equal length
    #  norm of E field matrix
    norm_e = hypot.(ex, ey)
    ex = ex ./ norm_e
    ey = ey ./ norm_e

    return (norm_e=norm_e, ex=ex, ey=ey)

end
