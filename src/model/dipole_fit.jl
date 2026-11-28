export dipole_fit

# ------------------------------------------------------------------------------
# Draft: single equivalent current dipole fit, homogeneous spherical head.
# Coordinates: NeuroAnalyzer uses RAS (x=right, y=anterior, z=superior).
# Forward model: Brody series for the surface potential of a dipole in a
# homogeneous conducting sphere (Neumann BC), electrodes projected onto sphere.
# ------------------------------------------------------------------------------

"""
    _fit_sphere(P)

Least-squares sphere fit to n×3 points. Returns (center, radius).
"""
function _fit_sphere(P::AbstractMatrix{<:Real})
    A = hcat(2 .* P, ones(size(P, 1)))
    b = vec(sum(abs2, P; dims = 2))
    s = A \ b
    c = s[1:3]
    return c, sqrt(s[4] + sum(abs2, c))
end

"""
    _leadfield(r0, U, R, σ; <keyword arguments>)

Lead field (n_channels × 3) for a dipole in a sphere; the surface potentials are `V = L * q` for dipole moment `q`.

The default is a homogeneous sphere of conductivity `σ`. For concentric multi-shell spheres pass the shell factors `g` (`g[n]` scales the Legendre term of degree `n`, see `_shell_factors()`; `σ` is then the conductivity of the innermost shell).

# Arguments

- `r0::AbstractVector{<:Real}`: dipole position (3-vector, relative to the sphere center, same length unit as `R`)
- `U::AbstractMatrix{<:Real}`: n×3 matrix of unit vectors toward the electrodes
- `R::Real`: sphere radius
- `σ::Real`: conductivity [S/m]
- `nmax::Int=400`: maximum degree of the Legendre series
- `tol::Real=1e-10`: convergence tolerance of the Legendre series
- `g::Union{Nothing, AbstractVector{<:Real}}=nothing`: shell factors for multi-shell spheres; `length(g)` must be ≥ `nmax`

# Returns

- `Matrix{Float64}`: lead field (n_channels × 3)
"""
function _leadfield(
    r0::AbstractVector{<:Real},
    U::AbstractMatrix{<:Real},
    R::Real,
    σ::Real;
    nmax::Int = 400,
    tol::Real = 1e-10,
    g::Union{Nothing, AbstractVector{<:Real}} = nothing,
)::Matrix{Float64}
    isnothing(g) ||
        length(g) >= nmax ||
        throw(ArgumentError("Length of g must be ≥ nmax ($nmax)."))
    nch = size(U, 1)
    ρ = norm(r0)
    L = zeros(nch, 3)
    g1 = isnothing(g) ? 1.0 : Float64(g[1])
    if ρ <= 1e-12
        # only the n = 1 term remains
        @inbounds for i = 1:nch, j = 1:3
            L[i, j] = 3 * g1 * U[i, j]
        end
        return L ./ (4π * σ * R^2)
    end
    r̂1, r̂2, r̂3 = r0[1] / ρ, r0[2] / ρ, r0[3] / ρ
    @inbounds for i = 1:nch
        u1, u2, u3 = U[i, 1], U[i, 2], U[i, 3]
        c = clamp(u1 * r̂1 + u2 * r̂2 + u3 * r̂3, -1.0, 1.0)
        # Σ w m Pₙ (radial part) and Σ w Pₙ' (tangential part) for n ≥ 2
        a = 0.0
        b = 0.0
        # P_0, P_1
        Pm, P = 1.0, c
        # P_0', P_1'
        dPm, dP = 0.0, 1.0
        for n = 1:(nmax - 1)
            Pn = ((2n + 1) * c * P - n * Pm) / (n + 1)
            dPn = dPm + (2n + 1) * P
            Pm, P = P, Pn
            dPm, dP = dP, dPn
            m = n + 1
            w = (2m + 1) / m * (ρ / R)^(m - 1)
            isnothing(g) || (w *= g[m])
            a += w * m * P
            b += w * dP
            abs(w) * m^2 < tol && break
        end
        # n = 1 term (3g₁u, independent of r̂0) + a r̂0 + b t, with tangential part t = u - c r̂0
        ar = a - b * c
        L[i, 1] = (3 * g1 + b) * u1 + ar * r̂1
        L[i, 2] = (3 * g1 + b) * u2 + ar * r̂2
        L[i, 3] = (3 * g1 + b) * u3 + ar * r̂3
    end
    return L ./ (4π * σ * R^2)
end

"""
    _shell_factors(radii, σ, nmax)

Shell factors of a concentric multi-shell sphere: `g[n]` is the ratio of the outer-surface potential of the Legendre term of degree `n` of a source in the innermost shell to that of a homogeneous sphere (conductivity `σ[1]`, outer radius `radii[end]`). `g ≡ 1` for equal conductivities.

# Arguments

- `radii::AbstractVector{<:Real}`: outer radii of the shells, innermost first
- `σ::AbstractVector{<:Real}`: conductivities of the shells, innermost first
- `nmax::Int`: highest degree

# Returns

- `Vector{Float64}`
"""
function _shell_factors(
    radii::AbstractVector{<:Real},
    σ::AbstractVector{<:Real},
    nmax::Int,
)::Vector{Float64}
    N = length(radii)
    length(σ) == N || throw(ArgumentError("radii and σ must have the same length."))
    issorted(radii) && all(>(0), radii) ||
        throw(ArgumentError("radii must be positive and increasing."))
    all(>(0), σ) || throw(ArgumentError("σ must be > 0."))
    # high precision: rⁿ and r⁻⁽ⁿ⁺¹⁾ span many orders of magnitude for large n
    r = big.(radii) ./ big(radii[end])
    s = big.(σ)
    g = zeros(nmax)
    for n = 1:nmax
        # potential of degree n in shell k: Aₖ rⁿ + Bₖ r⁻⁽ⁿ⁺¹⁾; in shell 1 B₁ = 1/σ₁ (primary source)
        # unknowns: A₁, A₂, B₂, …, A_N, B_N
        M = zeros(BigFloat, 2N - 1, 2N - 1)
        rhs = zeros(BigFloat, 2N - 1)
        iA(k) = k == 1 ? 1 : 2k - 2
        iB(k) = 2k - 1
        f(x) = x^n
        h(x) = x^(-(n + 1))
        df(x) = n * x^(n - 1)
        dh(x) = -(n + 1) * x^(-(n + 2))
        row = 0
        for k = 1:(N - 1)
            x = r[k]
            # continuity of the potential
            row += 1
            M[row, iA(k)] += f(x)
            k == 1 ? (rhs[row] -= h(x) / s[1]) : (M[row, iB(k)] += h(x))
            M[row, iA(k + 1)] -= f(x)
            M[row, iB(k + 1)] -= h(x)
            # continuity of the normal current
            row += 1
            M[row, iA(k)] += s[k] * df(x)
            k == 1 ? (rhs[row] -= dh(x)) : (M[row, iB(k)] += s[k] * dh(x))
            M[row, iA(k + 1)] -= s[k + 1] * df(x)
            M[row, iB(k + 1)] -= s[k + 1] * dh(x)
        end
        # no current leaves the outer surface (r = 1)
        row += 1
        M[row, iA(N)] = df(big(1))
        N == 1 ? (rhs[row] = -dh(big(1)) / s[1]) : (M[row, iB(N)] = dh(big(1)))
        x = M \ rhs
        V = x[iA(N)] + (N == 1 ? 1 / s[1] : x[iB(N)])
        # homogeneous sphere: Vₙ = (2n + 1) / (n σ₁)
        g[n] = Float64(V * n * s[1] / (2n + 1))
    end
    return g
end

# average-reference operator
_avgref(n::Int) = Matrix{Float64}(I, n, n) .- 1 / n

# relative residual variance for dipole position x, with optimal moments
function _dipole_cost(x, U, R, σ, H, D, ρmax)
    ρ = norm(x)
    ρ > ρmax * R && return 1.0 + (ρ / R - ρmax) * 10      # soft wall
    L = H * _leadfield(x, U, R, σ)
    Q = L \ D
    return sum(abs2, D .- L * Q) / sum(abs2, D)
end

"""
    dipole_fit(s, locs; <keyword arguments>)

Fit a single current dipole to EEG data.

# Arguments

- `s::AbstractMatrix{<:Real}`: channels × time points (signal units, e.g. μV)
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

- `NeuroAnalyzer.DIPOLE`:
    - `pos` (n×3 for `:moving`, else 1×3, in `locs` units, original coordinates)
    - `q` (3 × time points)
    - `gof` (goodness of fit, 0–1, per fit)
    - `sphere_center`
    - `sphere_radius`
"""
function dipole_fit(
    s::AbstractMatrix{<:Real},
    locs::DataFrame;
    mode::Symbol = :fixed,
    σ::Real = 0.33,
    loc_scale::Real = 1e-3,
    ρmax::Real = 0.95,
    n_grid::Int = 7,
    n_starts::Int = 5,
)::NeuroAnalyzer.DIPOLE
    _wip()

    # validate
    _check_var(mode, [:fixed, :moving], "mode")
    nch, nt = size(s)
    DataFrames.nrow(locs) == nch ||
        throw(ArgumentError("Number of rows in locs must equal number of channels."))
    nch >= 8 || throw(ArgumentError("Too few channels for a stable fit (≥ 8 needed)."))
    0 < ρmax < 1 || throw(ArgumentError("ρmax must be in (0, 1)."))
    σ > 0 || throw(ArgumentError("σ must be > 0."))
    loc_scale > 0 || throw(ArgumentError("loc_scale must be > 0."))
    n_grid >= 2 || throw(ArgumentError("n_grid must be ≥ 2."))
    n_starts >= 1 || throw(ArgumentError("n_starts must be ≥ 1."))

    P = Matrix{Float64}(locs[:, [:loc_x, :loc_y, :loc_z]])
    c, R_loc = _fit_sphere(P)
    U = P .- c'
    # project electrodes to sphere surface
    U ./= sqrt.(sum(abs2, U; dims = 2))
    R = R_loc * loc_scale

    H = _avgref(nch)
    # average-referenced data
    D = H * Float64.(s)

    # candidate grid inside the sphere (in meters)
    ax = range(-ρmax, ρmax; length = n_grid) .* R
    grid = [[x, y, z] for x in ax for y in ax for z in ax if norm([x, y, z]) <= ρmax * R]

    function fit_window(Dw, seeds)
        costs = [_dipole_cost(g, U, R, σ, H, Dw, ρmax) for g in grid]
        starts = vcat(seeds, grid[sortperm(costs)[1:min(n_starts, length(grid))]])
        best_x, best_f = starts[1], Inf
        for x0 in starts
            res = optimize(x -> _dipole_cost(x, U, R, σ, H, Dw, ρmax), x0, NelderMead(),
                Optim.Options(; iterations = 500, g_abstol = 1e-8))
            if Optim.minimum(res) < best_f
                best_f, best_x = Optim.minimum(res), Optim.minimizer(res)
            end
        end
        L = H * _leadfield(best_x, U, R, σ)
        return best_x, L \ Dw, 1 - best_f
    end

    if mode == :fixed
        x, Q, g = fit_window(D, Vector{Float64}[])
        return NeuroAnalyzer.DIPOLE(reshape(x ./ loc_scale .+ c, 1, 3), Q, [g], c, R_loc)
    else
        pos = zeros(nt, 3)
        Q = zeros(3, nt)
        gof = zeros(nt)
        prev = Vector{Float64}[]
        for t = 1:nt
            x, q, g = fit_window(D[:, t:t], prev)
            pos[t, :] = x ./ loc_scale .+ c
            Q[:, t] = vec(q)
            gof[t] = g
            prev = [x]
        end
        return NeuroAnalyzer.DIPOLE(pos, Q, gof, c, R_loc)
    end
end
