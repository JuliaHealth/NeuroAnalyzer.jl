using LinearAlgebra
using Statistics
using DataFrames
using Optim

using Random
Random.seed!(1)

# TEST
# synthetic 64-"channel" cap on a 90 mm sphere, centered at (0,0,10) mm
n = 64
pts = [begin v = randn(3); v ./= norm(v); v end for _ in 1:n]
pts = filter(v -> v[3] > -0.3, pts)
while length(pts) < n
    v = randn(3); v ./= norm(v); v[3] > -0.3 && push!(pts, v)
end
Rmm = 90.0
ctr = [0.0, 0.0, 10.0]
locs = DataFrame(loc_x=[Rmm*v[1]+ctr[1] for v in pts],
                 loc_y=[Rmm*v[2]+ctr[2] for v in pts],
                 loc_z=[Rmm*v[3]+ctr[3] for v in pts])
U = reduce(vcat, [v' for v in pts])
 
# test 1: dipole at center, radial moment -> V = 3 q cosθ / (4πσR²)
R = Rmm * 1e-3; σ = 0.33
L0 = _leadfield([0.0, 0.0, 0.0], U, R, σ)
q = [0.0, 0.0, 1.0]
V = L0 * q
@show maximum(abs.(V .- 3 .* U[:, 3] ./ (4π * σ * R^2)))
 
# test 2: recover an off-center dipole (fixed mode), with noise
x_true = [0.02, -0.015, 0.03]                    # meters
qt = [1.0, 0.5, -0.8]
H = _avgref(n)
tt = 0:0.01:0.1
wave = sin.(2π .* 10 .* tt)'
S = (H * _leadfield(x_true, U, R, σ) * qt) * wave
S .+= 0.01 * maximum(abs.(S)) .* randn(size(S))
res = dipole_fit(S, locs)
@show res.pos .- (x_true' .* 1e3 .+ ctr')        # error in mm
@show res.gof res.sphere_center res.sphere_radius

# ---------------------------------------------------------------------------
# Draft: single equivalent current dipole fit, homogeneous spherical head.
# Coordinates: NeuroAnalyzer uses RAS (x=right, y=anterior, z=superior).
# Forward model: Brody series for the surface potential of a dipole in a
# homogeneous conducting sphere (Neumann BC), electrodes projected onto sphere.
# ---------------------------------------------------------------------------

"""
    _fit_sphere(P)

Least-squares sphere fit to n×3 points. Returns (center, radius).
"""
function _fit_sphere(P::AbstractMatrix{<:Real})
    A = hcat(2 .* P, ones(size(P, 1)))
    b = vec(sum(abs2, P; dims=2))
    s = A \ b
    c = s[1:3]
    return c, sqrt(s[4] + sum(abs2, c))
end

"""
    _leadfield(r0, U, R, σ; nmax=400, tol=1e-10)

Lead field (n_channels × 3) for a dipole at `r0` (3-vector, relative to sphere center, same length unit as `R`). `U` is n×3 matrix of unit vectors toward the electrodes. V = L * q.
"""
function _leadfield(
    r0::AbstractVector{<:Real},
    U::AbstractMatrix{<:Real},
    R::Real,
    σ::Real;
    nmax::Int=400,
    tol::Real=1e-10
)::Matrix{Float64}
    nch = size(U, 1)
    ρ = norm(r0)
    L = zeros(nch, 3)
    r̂0 = ρ > 1e-12 ? r0 ./ ρ : zeros(3)
    for i in 1:nch
        u = U[i, :]
        # n = 1 term (independent of r̂0)
        g = 3 .* u
        if ρ > 1e-12
            c = clamp(dot(u, r̂0), -1.0, 1.0)
            # tangential part
            t = u .- c .* r̂0
            # P_0, P_1
            Pm, P = 1.0, c
            # P_0', P_1'
            dPm, dP = 0.0, 1.0
            for n in 1:(nmax - 1)
                Pn = ((2n + 1) * c * P - n * Pm) / (n + 1)
                dPn = dPm + (2n + 1) * P
                Pm, P = P, Pn
                dPm, dP = dP, dPn
                m = n + 1
                w = (2m + 1) / m * (ρ / R)^(m - 1)
                g .+= w .* (m * P .* r̂0 .+ dP .* t)
                w * m^2 < tol && break
            end
        end
        L[i, :] = g
    end
    return L ./ (4π * σ * R^2)
end

# average-reference operator
_avgref(n::Int) = Matrix{Float64}(I, n, n) .- 1 / n

# relative residual variance for dipole position x, with optimal moments
function _cost(x, U, R, σ, H, D, ρmax)
    ρ = norm(x)
    ρ > ρmax * R && return 1.0 + (ρ / R - ρmax) * 10      # soft wall
    L = H * _leadfield(x, U, R, σ)
    Q = L \ D
    return sum(abs2, D .- L * Q) / sum(abs2, D)
end

"""
    dipole_fit(s, locs; mode=:fixed, σ=0.33, loc_scale=1e-3, ρmax=0.95,
               n_grid=7, n_starts=5)

Fit a single current dipole to EEG data.

# Arguments

- `s::AbstractMatrix{<:Real}`: channels × time points (signal units, e.g. µV)
- `locs::DataFrame`: columns `loc_x`, `loc_y`, `loc_z`, one row per channel (same order as `s`)
- `mode::Symbol=:fixed`:
    - `:fixed`: one position for all time points, moment varies with time;
    - `:moving`: independent fit per time point
- `σ::Real=0.33`: conductivity [S/m]
- `loc_scale::Real=1e-3`: factor converting `locs` units to meters (1e-3 = locs in mm)
- `ρmax::Real=0.95`: max dipole eccentricity (fraction of sphere radius)
- `n_grid::Int=7`: grid points per axis for the initial search
- `n_starts::Int=5`: best grid points refined by Nelder-Mead

# Returns

Named tuple:

- `pos` (n×3 for `:moving`, else 1×3, in `locs` units, original coordinates)
- `q` (3 × time points)
- `gof` (goodness of fit, 0–1, per fit)
- `sphere_center`
- `sphere_radius`.
"""
function dipole_fit(
    s::AbstractMatrix{<:Real},
    locs::DataFrame;
    mode::Symbol=:fixed,
    σ::Real=0.33,
    loc_scale::Real=1e-3,
    ρmax::Real=0.95,
    n_grid::Int=7,
    n_starts::Int=5
)
    # validate
    @assert mode in (:fixed, :moving) "mode must be :fixed or :moving."
    nch, nt = size(s)
    @assert nrow(locs) == nch "number of rows in locs must equal number of channels."
    @assert nch >= 8 "too few channels for a stable fit."
    @assert 0 < ρmax < 1 "ρmax must be in (0, 1)."

    P = Matrix{Float64}(locs[:, [:loc_x, :loc_y, :loc_z]])
    c, R_loc = _fit_sphere(P)
    U = P .- c'
    # project electrodes to sphere surface
    U ./= sqrt.(sum(abs2, U; dims=2))
    R = R_loc * loc_scale

    H = _avgref(nch)
    # average-referenced data
    D = H * Float64.(s)

    # candidate grid inside the sphere (in meters)
    ax = range(-ρmax, ρmax, length=n_grid) .* R
    grid = [[x, y, z] for x in ax for y in ax for z in ax if norm([x, y, z]) <= ρmax * R]

    function fit_window(Dw, seeds)
        costs = [_cost(g, U, R, σ, H, Dw, ρmax) for g in grid]
        starts = vcat(seeds, grid[sortperm(costs)[1:min(n_starts, length(grid))]])
        best_x, best_f = starts[1], Inf
        for x0 in starts
            res = optimize(x -> _cost(x, U, R, σ, H, Dw, ρmax), x0, NelderMead(),
                           Optim.Options(iterations=500, g_abstol=1e-8))
            if Optim.minimum(res) < best_f
                best_f, best_x = Optim.minimum(res), Optim.minimizer(res)
            end
        end
        L = H * _leadfield(best_x, U, R, σ)
        return best_x, L \ Dw, 1 - best_f
    end

    if mode == :fixed
        x, Q, g = fit_window(D, Vector{Float64}[])
        return (pos = reshape(x ./ loc_scale .+ c, 1, 3), q = Q, gof = [g],
                sphere_center = c, sphere_radius = R_loc)
    else
        pos = zeros(nt, 3); Q = zeros(3, nt); gof = zeros(nt)
        prev = Vector{Float64}[]
        for t in 1:nt
            x, q, g = fit_window(D[:, t:t], prev)
            pos[t, :] = x ./ loc_scale .+ c
            Q[:, t] = vec(q); gof[t] = g
            prev = [x]
        end
        return (pos = pos, q = Q, gof = gof, sphere_center = c, sphere_radius = R_loc)
    end
end
