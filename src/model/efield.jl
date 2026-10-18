export efield2d

"""
    efield2d(q, qx, qy)

Create model of 2-dimensional electric field of point charges on a 100 × 100 grid spanning `[-1, 1]` in both directions.

The field is computed in grid units (Coulomb's law with k = 1). At the position of a charge the field is undefined and set to 0.

# Arguments

- `q::Vector{Int64}`: charges values
- `qx::Vector{Int64}`: charges X positions (grid column, 1 to 100)
- `qy::Vector{Int64}`: charges Y positions (grid row, 1 to 100)

# Returns

Named tuple:

- `qq::Vector{Vector{Float64}}`: charges positions `[x, y]`
- `norm_e::Matrix{Float64}`: magnitude of the electric field (rows: y, columns: x)
- `ex::Matrix{Float64}`: X component of the unit field direction
- `ey::Matrix{Float64}`: Y component of the unit field direction
"""
function efield2d(
    q::Vector{Int64},
    qx::Vector{Int64},
    qy::Vector{Int64},
)::@NamedTuple{
    qq::Vector{Vector{Float64}},
    norm_e::Matrix{Float64},
    ex::Matrix{Float64},
    ey::Matrix{Float64},
}
    # validate
    length(qx) == length(q) ||
        throw(ArgumentError("Length of qx and number of charges must be equal."))
    length(qy) == length(q) ||
        throw(ArgumentError("Length of qy and number of charges must be equal."))

    m = 100
    n = 100

    # number of charges
    nq = length(q)
    for idx in eachindex(q)
        _in(qx[idx], (1, m), "qx")
        _in(qy[idx], (1, n), "qy")
    end

    x = collect(range(-1, 1, m))
    y = collect(range(-1, 1, n))

    #strength
    ex = zeros(n, m)
    ey = zeros(n, m)

    qq = Vector{Vector{Float64}}()

    @inbounds for idx = 1:nq
        # column qx ↔ x, row qy ↔ y
        push!(qq, [x[qx[idx]], y[qy[idx]]])
        for idx1 = 1:n
            for idx2 = 1:m
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
    # (0 where the field vanishes, instead of 0 / 0 = NaN)
    ex = [e == 0 ? 0.0 : ex[i] / e for (i, e) in pairs(norm_e)]
    ey = [e == 0 ? 0.0 : ey[i] / e for (i, e) in pairs(norm_e)]

    return (; qq, norm_e, ex, ey)
end
