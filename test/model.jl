using Test
using DataFrames
using LinearAlgebra
using Random
using GLMakie

@info "Test: efield2d()"
r = efield2d([1, -1], [25, 75], [50, 50])
@test size(r.norm_e) == (100, 100)
# charge positions [x, y]: columns 25 and 75 of row 50
@test r.qq[1] ≈ [range(-1, 1, 100)[25], range(-1, 1, 100)[50]]
@test r.qq[2] ≈ [range(-1, 1, 100)[75], range(-1, 1, 100)[50]]
# unit field direction, no NaN at the charge positions
@test !any(isnan, r.ex) && !any(isnan, r.ey)
@test hypot(r.ex[10, 10], r.ey[10, 10]) ≈ 1
@test_throws ArgumentError efield2d([1, -1], [25, 75], [50])
@test_throws ArgumentError efield2d([1], [101], [50])

@info "Test: tes_model()"
@test tes_model(anode = "F3", cathode = "Fp2") isa GLMakie.Figure
@test_throws ArgumentError tes_model(anode = "XX", cathode = "Fp2")

@info "Test: dipole_fit()"
# synthetic 64-channel cap on a 90 mm sphere centered at (0, 0, 10) mm
rng = Random.Xoshiro(1)
n = 64
pts = Vector{Float64}[]
while length(pts) < n
    local v = randn(rng, 3)
    v ./= norm(v)
    v[3] > -0.3 && push!(pts, v)
end
Rmm = 90.0
ctr = [0.0, 0.0, 10.0]
locs = DataFrame(
    loc_x = [Rmm * v[1] + ctr[1] for v in pts],
    loc_y = [Rmm * v[2] + ctr[2] for v in pts],
    loc_z = [Rmm * v[3] + ctr[3] for v in pts],
)
U = reduce(vcat, [v' for v in pts])
R = Rmm * 1e-3
σ = 0.33
# dipole at the center, radial moment: V = 3 q cosθ / (4πσR²)
V = NeuroAnalyzer._leadfield([0.0, 0.0, 0.0], U, R, σ) * [0.0, 0.0, 1.0]
@test maximum(abs.(V .- 3 .* U[:, 3] ./ (4π * σ * R^2))) < 1e-8
# recover an off-center dipole (fixed mode) from noisy data
x_true = [0.02, -0.015, 0.03] # m
qt = [1.0, 0.5, -0.8]
H = NeuroAnalyzer._avgref(n)
wave = sin.(2π .* 10 .* (0:0.01:0.1))'
S = (H * NeuroAnalyzer._leadfield(x_true, U, R, σ) * qt) * wave
S .+= 0.01 * maximum(abs.(S)) .* randn(rng, size(S))
d = dipole_fit(S, locs)
@test d isa NeuroAnalyzer.DIPOLE
@test size(d.pos) == (1, 3) && size(d.q) == (3, length(wave))
# position error < 2 mm, goodness of fit close to 1
@test norm(vec(d.pos) .- (x_true .* 1e3 .+ ctr)) < 2
@test d.gof[1] > 0.99
@test d.sphere_radius ≈ Rmm atol = 1e-6
@test_throws ArgumentError dipole_fit(S[1:4, :], locs[1:4, :])

@info "Test: recorder / tester helpers"
@test NeuroAnalyzer._parse_eda_sample("gsr:1.25") == 1.25
@test isnothing(NeuroAnalyzer._parse_eda_sample("noise"))
smp = NeuroAnalyzer._parse_tpt_sample("tpt: 1 -2 3 0.5 -1.25 9.81")
@test (smp.x, smp.y, smp.z, smp.accx, smp.accy, smp.accz) == (1.0, -2.0, 3.0, 0.5, -1.25, 9.81)
@test isnothing(NeuroAnalyzer._parse_tpt_sample("tpt: 1 2"))
# taps: the last press has no release (still pressed at the end) and is dropped
n_taps, t_ms, d_ms = NeuroAnalyzer._ftt_taps([0.1, 0.5, 0.9], [0.2, 0.65])
@test n_taps == 2 && t_ms ≈ [100.0, 500.0] && d_ms ≈ [100.0, 150.0]
res = NeuroAnalyzer._ftt_result([[0.1], [0.2, 0.4]], [[0.2], [0.3, 0.5]], [Float64[], [0.1]], [Float64[], [0.2]], 2)
@test res.taps == [1, 2] && res.taps_int == [0, 1]
@test res.tap_t[2] ≈ [200.0, 400.0] && res.tap_d_int[2] ≈ [100.0]


@info "Test: dipole fit helpers"
# least-squares sphere through the electrodes: the cap's center and radius
c_fit, r_fit = NeuroAnalyzer._fit_sphere(Matrix(locs[!, [:loc_x, :loc_y, :loc_z]]))
@test c_fit ≈ ctr atol = 1e-8
@test r_fit ≈ Rmm
# relative residual of the fit: ≈ 0 at the true position (noise only), larger elsewhere, soft wall outside
D = H * S
c_true = NeuroAnalyzer._dipole_cost(x_true, U, R, σ, H, D, 0.95)
@test c_true < 1e-3
@test NeuroAnalyzer._dipole_cost([-0.03, 0.02, 0.0], U, R, σ, H, D, 0.95) > 10 * c_true
@test NeuroAnalyzer._dipole_cost([0.0, 0.0, 0.099], U, R, σ, H, D, 0.95) > 1.0

true
