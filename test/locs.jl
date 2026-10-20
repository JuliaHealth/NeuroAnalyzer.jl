using NeuroAnalyzer
using Test
using DataFrames

# location values are stored exactly (no rounding); reference values have 2 decimals
_approx₂(a::AbstractString, b::AbstractString) = a == b
_approx₂(a, b) = isapprox(a, b; atol = 0.01)
≈₂(a, b) = all(map(_approx₂, Tuple(a), Tuple(b)))


@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
locs = import_locs(joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"))

@info "Test: add_locs()"
eeg_tmp = add_locs(eeg; locs = locs)
@test nrow(eeg_tmp.locs) == 23
add_locs!(eeg; locs = locs)
@test nrow(eeg.locs) == 23

@info "Test: locs_details()"
@test locs_details(eeg, ch = "Fp1", out = false) ≈₂ (
    label = "Fp1",
    theta_pl = 108.0,
    radius_pl = 1.0,
    x = -0.31,
    y = 0.95,
    z = -0.03,
    theta_sph = 108.0,
    radius_sph = 1.0,
    phi_sph = -2.0,
)

@info "Test: cart2pol()"
@test cart2pol(1.0, 1.0) ≈₂ (1.41, 45.0)

@info "Test: pol2cart()"
@test pol2cart(1.41, 45.0) ≈₂ (1, 1)

@info "Test: sph2cart()"
@test sph2cart(1.73, 45.0, 35.27) ≈₂ (1.0, 1.0, 1.0)

@info "Test: cart2sph()"
@test cart2sph(1.0, 1.0, 1.0) ≈₂ (1.73, 45.0, 35.26)

@info "Test: sph2pol()"
@test sph2pol(1.73, 45.0, 35.26) ≈₂ (1.41, 45.0)

@info "Test: locs_sph2cart()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_tmp = locs_sph2cart(locs_tmp)
@test locs_tmp[1, :loc_x] ≈₂ 1.0

@info "Test: locs_sph2cart!()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_sph2cart!(locs_tmp)
@test locs_tmp[1, :loc_x] ≈₂ 1.0

@info "Test: locs_cart2sph()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_x] = 1
locs_tmp[1, :loc_y] = 0
locs_tmp = locs_cart2sph(locs_tmp)
@test locs_tmp[1, :loc_theta_sph] == 0.0

@info "Test: locs_cart2sph!()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_x] = 1
locs_tmp[1, :loc_y] = 0
locs_cart2sph!(locs_tmp)
@test locs_tmp[1, :loc_theta_sph] == 0.0

@info "Test: locs_cart2pol()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_x] = 1
locs_tmp[1, :loc_y] = 0
locs_tmp = locs_cart2pol(locs_tmp)
@test locs_tmp[1, :loc_theta] == 0.0

@info "Test: locs_cart2pol!()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_x] = 1
locs_tmp[1, :loc_y] = 0
locs_cart2pol!(locs_tmp)
@test locs_tmp[1, :loc_theta] == 0.0

@info "Test: locs_sph2pol()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_tmp[1, :loc_phi_sph] = 0
locs_tmp = locs_sph2pol(locs_tmp)
@test locs_tmp[1, :loc_radius] == 1.0
@test locs_tmp[1, :loc_theta] == 0.0

@info "Test: locs_sph2pol!()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_tmp[1, :loc_phi_sph] = 0
locs_sph2pol!(locs_tmp)
@test locs_tmp[1, :loc_radius] == 1.0
@test locs_tmp[1, :loc_theta] == 0.0

@info "Test: edit_locs()"
eeg_tmp = deepcopy(eeg)
@test locs_details(eeg_tmp, ch = "Fp1", out = false) ≈₂ (
    label = "Fp1",
    theta_pl = 108.0,
    radius_pl = 1.0,
    x = -0.31,
    y = 0.95,
    z = -0.03,
    theta_sph = 108.0,
    radius_sph = 1.0,
    phi_sph = -2.0,
)
eeg_tmp = edit_locs(eeg_tmp; ch = "Fp1", x = 0.5)
@test locs_details(eeg_tmp, ch = "Fp1", out = false) ≈₂ (
    label = "Fp1",
    theta_pl = 108.0,
    radius_pl = 1.0,
    x = 0.5,
    y = 0.95,
    z = -0.03,
    theta_sph = 108.0,
    radius_sph = 1.0,
    phi_sph = -2.0,
)
edit_locs!(eeg_tmp; ch = "Fp1", x = 0.5)
@test locs_details(eeg_tmp, ch = "Fp1", out = false) ≈₂ (
    label = "Fp1",
    theta_pl = 108.0,
    radius_pl = 1.0,
    x = 0.5,
    y = 0.95,
    z = -0.03,
    theta_sph = 108.0,
    radius_sph = 1.0,
    phi_sph = -2.0,
)

@info "Test: locs_flipx()"
@test locs[1, :loc_x] ≈₂ -0.31
locs2 = locs_flipx(locs)
@test locs2[1, :loc_x] ≈₂ 0.31

@info "Test: locs_flipy()"
@test locs[1, :loc_y] ≈₂ 0.95
locs2 = locs_flipy(locs)
@test locs2[1, :loc_y] ≈₂ -0.95

@info "Test: locs_flipz()"
@test locs[1, :loc_z] ≈₂ -0.03
locs2 = locs_flipz(locs)
@test locs2[1, :loc_z] ≈₂ 0.03

@info "Test: locs_scale()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_radius_sph] == 1.0
locs2 = locs_scale(locs; r = 1.2)
@test locs2[1, :loc_radius] == 1.2
@test locs2[1, :loc_radius_sph] == 1.2

@info "Test: locs_normalize()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_radius_sph] == 1.0
locs2 = locs_normalize(locs)
@test locs2[1, :loc_radius] == 1.0
@test locs2[1, :loc_radius_sph] == 1.0

@info "Test: locs_swapxy()"
@test locs[1, :loc_theta] == 108.0
@test locs2[1, :loc_theta_sph] == 108.0
locs2 = locs_swapxy(locs)
@test locs2[1, :loc_theta] == 198.0
@test locs2[1, :loc_theta_sph] == -162.0

@info "Test: locs_rotx()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_theta] == 108.0
@test locs[1, :loc_theta_sph] == 108.0
@test locs[1, :loc_phi_sph] == -2.0
locs2 = locs_rotx(locs; a = 20)
@test locs2[1, :loc_radius] == 1.0
@test locs2[1, :loc_theta] == 108.0
@test locs2[1, :loc_theta_sph] ≈₂ 108.84
@test locs2[1, :loc_phi_sph] ≈₂ 16.99

@info "Test: locs_roty()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_theta] == 108.0
@test locs[1, :loc_theta_sph] == 108.0
@test locs[1, :loc_phi_sph] == -2.0
locs2 = locs_roty(locs; a = 20)
@test locs2[1, :loc_radius] == 1.0
@test locs2[1, :loc_theta] == 108.0
@test locs2[1, :loc_theta_sph] ≈₂ 107.63
@test locs2[1, :loc_phi_sph] ≈₂ 4.18

@info "Test: locs_rotz()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_theta] == 108.0
@test locs[1, :loc_theta_sph] == 108.0
@test locs[1, :loc_phi_sph] == -2.0
locs2 = locs_rotz(locs; a = 20)
@test locs2[1, :loc_radius] == 1.0
@test locs2[1, :loc_theta] == 128.0
@test locs2[1, :loc_theta_sph] == 128.0
@test locs2[1, :loc_phi_sph] == -2.0

@info "Test: locs_generate()"
locs_tmp = locs_generate(locs_tmp)
eeg_tmp = locs_generate(eeg_tmp)
@test locs_tmp isa DataFrame
@test eeg_tmp isa NeuroAnalyzer.NEURO

@info "Test: locs_origin()"
@test locs_origin(locs, y = 0.2) isa DataFrame


##
## values (geometry) and in-place variants
##

xyz(l, i) = [l[i, :loc_x], l[i, :loc_y], l[i, :loc_z]]
coords(l) = Matrix(l[!, [:loc_x, :loc_y, :loc_z]])

@info "Test: coordinate conversions values"
# polar/spherical angles in degrees: theta from +x (right) towards +y (nose), phi = elevation
@test all(pol2cart(1.0, 90.0) .≈ (0.0, 1.0) .|| abs.(pol2cart(1.0, 90.0) .- (0.0, 1.0)) .< 1e-12)
@test collect(pol2cart(2.0, 180.0)) ≈ [-2.0, 0.0] atol = 1e-12
@test collect(sph2cart(1.0, 0.0, 90.0)) ≈ [0.0, 0.0, 1.0] atol = 1e-12
@test collect(sph2cart(2.0, 90.0, 0.0)) ≈ [0.0, 2.0, 0.0] atol = 1e-12
r_tmp, t_tmp, p_tmp = cart2sph(0.0, 0.0, 2.0)
@test r_tmp ≈ 2.0 && p_tmp ≈ 90.0
@test collect(cart2sph(sph2cart(1.5, -120.0, 30.0)...)) ≈ [1.5, -120.0, 30.0]
@test collect(cart2pol(pol2cart(1.5, 135.0)...)) ≈ [1.5, 135.0]
# sph2pol: azimuthal projection, radius = horizontal distance
@test collect(sph2pol(1.0, 30.0, 0.0)) ≈ [1.0, 30.0]
@test sph2pol(1.0, 30.0, 90.0)[1] ≈ 0.0 atol = 1e-12
# pol2sph: the polar layout on the equator
@test pol2sph(0.8, 45.0) == (0.8, 45.0, 0.0)

@info "Test: locs_pol2cart() / locs_pol2sph()"
lp = DataFrame(:label => ["a", "b", "c"], :loc_radius => [1.0, 1.0, 0.5], :loc_theta => [0.0, 90.0, 180.0],
               :loc_x => zeros(3), :loc_y => zeros(3), :loc_z => zeros(3), :loc_radius_sph => zeros(3),
               :loc_theta_sph => zeros(3), :loc_phi_sph => zeros(3))
l2 = locs_pol2cart(lp)
@test l2[!, :loc_x] ≈ [1.0, 0.0, -0.5] atol = 1e-12
@test l2[!, :loc_y] ≈ [0.0, 1.0, 0.0] atol = 1e-12
@test lp[!, :loc_x] == zeros(3)
l3 = deepcopy(lp)
locs_pol2cart!(l3)
@test l3 == l2
l2 = locs_pol2sph(lp)
@test l2[!, :loc_radius_sph] == [1.0, 1.0, 0.5]
@test l2[!, :loc_theta_sph] == [0.0, 90.0, 180.0]
@test l2[!, :loc_phi_sph] == zeros(3)
l3 = deepcopy(lp)
locs_pol2sph!(l3)
@test l3 == l2

@info "Test: in-place variants = non-mutating variants"
for (f, f!, kw) in ((locs_flipx, locs_flipx!, ()), (locs_flipy, locs_flipy!, ()), (locs_flipz, locs_flipz!, ()),
                    (locs_rotx, locs_rotx!, (a = 30,)), (locs_roty, locs_roty!, (a = 30,)),
                    (locs_rotz, locs_rotz!, (a = 30,)), (locs_scale, locs_scale!, (r = 2,)),
                    (locs_normalize, locs_normalize!, ()), (locs_swapxy, locs_swapxy!, ()),
                    (locs_center, locs_center!, ()), (locs_origin, locs_origin!, (x = 0.1, z = -0.2)),
                    (locs_cart2sph, locs_cart2sph!, ()), (locs_cart2pol, locs_cart2pol!, ()),
                    (locs_sph2cart, locs_sph2cart!, ()), (locs_sph2pol, locs_sph2pol!, ()),
                    (locs_generate, locs_generate!, ()))
    l1 = f(locs; kw...)
    l2 = deepcopy(locs)
    f!(l2; kw...)
    @test Matrix(l2[!, 2:end]) ≈ Matrix(l1[!, 2:end])
end

@info "Test: locs_flip*() values"
l2 = locs_flipx(locs)
@test l2[!, :loc_x] ≈ -locs[!, :loc_x]
@test l2[!, :loc_y] == locs[!, :loc_y]
# flipping twice restores the locations
@test coords(locs_flipy(locs_flipy(locs))) ≈ coords(locs)
@test mod.(locs_flipx(locs_flipx(locs))[!, :loc_theta], 360) ≈ mod.(locs[!, :loc_theta], 360)
@test locs_flipz(locs)[!, :loc_phi_sph] ≈ -locs[!, :loc_phi_sph]
# Fp1 (left, θ = 108°) → right (θ = 72°)
@test locs_flipx(locs)[1, :loc_theta] ≈ 72.0

@info "Test: locs_rot*() values"
# rotations keep distances from the origin; 360° restores the locations
for f in (locs_rotx, locs_roty, locs_rotz)
    l2 = f(locs; a = 37)
    @test sqrt.(sum(abs2, coords(l2); dims = 2)) ≈ sqrt.(sum(abs2, coords(locs); dims = 2))
    @test coords(f(locs; a = 360)) ≈ coords(locs) atol = 1e-10
end
# 90° about z: (x, y) → (-y, x)
l2 = locs_rotz(locs; a = 90)
@test l2[!, :loc_x] ≈ -locs[!, :loc_y] atol = 1e-10
@test l2[!, :loc_y] ≈ locs[!, :loc_x] atol = 1e-10
# 90° about x: (y, z) → (-z, y)
l2 = locs_rotx(locs; a = 90)
@test l2[!, :loc_y] ≈ -locs[!, :loc_z] atol = 1e-10
@test l2[!, :loc_z] ≈ locs[!, :loc_y] atol = 1e-10
# 90° about y: (z, x) → (-x, z)
l2 = locs_roty(locs; a = 90)
@test l2[!, :loc_z] ≈ -locs[!, :loc_x] atol = 1e-10
@test l2[!, :loc_x] ≈ locs[!, :loc_z] atol = 1e-10

@info "Test: locs_scale() / locs_normalize() values"
l2 = locs_scale(locs; r = 2)
@test coords(l2) ≈ 2 .* coords(locs)
@test l2[!, :loc_radius] ≈ 2 .* locs[!, :loc_radius]
l2 = locs_scale(locs; r = 2, cart = false, spherical = false)
@test coords(l2) == coords(locs)
l2 = locs_normalize(locs_scale(locs; r = 3))
@test maximum(l2[!, :loc_radius]) ≈ 1.0
@test maximum(l2[!, :loc_radius_sph]) ≈ 1.0

@info "Test: locs_swapxy() values"
# (x, y) → (-y, x): a 90° rotation about z, the same as locs_rotz(a = 90)
l2 = locs_swapxy(locs)
@test l2[!, :loc_x] ≈ -locs[!, :loc_y]
@test l2[!, :loc_y] ≈ locs[!, :loc_x]
@test coords(l2) ≈ coords(locs_rotz(locs; a = 90)) atol = 1e-10
@test mod.(l2[!, :loc_theta], 360) ≈ mod.(locs[!, :loc_theta] .+ 90, 360)
@test coords(locs_swapxy(locs_swapxy(locs_swapxy(locs_swapxy(locs))))) ≈ coords(locs)
# swap + flip x = an exact exchange of the axes
l2 = locs_flipx(locs_swapxy(locs))
@test l2[!, :loc_x] ≈ locs[!, :loc_y] && l2[!, :loc_y] ≈ locs[!, :loc_x]

@info "Test: locs_center() values"
l2 = locs_center(locs_origin(locs; x = 0.3, y = -0.2))
cz = findfirst(==("Cz"), l2[!, :label])
@test l2[cz, :loc_x] ≈ 0.0 atol = 1e-12
@test l2[cz, :loc_y] ≈ 0.0 atol = 1e-12
@test_throws ArgumentError locs_center(locs[[1, 2], :])

@info "Test: locs_origin() values"
l2 = locs_origin(locs; x = 0.1, y = 0.2, z = -0.3)
@test coords(l2) ≈ coords(locs) .+ [0.1 0.2 -0.3]
# spherical / polar coordinates are recomputed from the moved Cartesian ones
@test l2[!, :loc_radius_sph] ≈ sqrt.(sum(abs2, coords(l2); dims = 2))[:]
eeg_tmp = deepcopy(eeg)
locs_origin!(eeg_tmp; z = 1)
@test eeg_tmp.locs[!, :loc_z] ≈ eeg.locs[!, :loc_z] .+ 1

@info "Test: locs_generate() values"
lg = DataFrame(:label => ["Cz", "Fpz", "Oz", "T7", "T8", "Fp1-F3", "xyz"], :loc_radius => zeros(7),
               :loc_theta => zeros(7), :loc_x => zeros(7), :loc_y => zeros(7), :loc_z => zeros(7),
               :loc_radius_sph => zeros(7), :loc_theta_sph => zeros(7), :loc_phi_sph => zeros(7))
lg = locs_generate(lg)
# 10/10 positions on the unit sphere: Cz on top, Fpz in front, Oz at the back, T7 left, T8 right
@test xyz(lg, 1) ≈ [0.0, 0.0, 1.0] atol = 1e-3
@test lg[2, :loc_y] > 0.9 && abs(lg[2, :loc_x]) < 1e-3
@test lg[3, :loc_y] < -0.9 && abs(lg[3, :loc_x]) < 1e-3
@test lg[4, :loc_x] < -0.9
@test lg[5, :loc_x] > 0.9
# bipolar label: the first electrode (Fp1) is used
@test lg[6, :loc_x] < 0 && lg[6, :loc_y] > 0.9
# unknown label: not moved
@test xyz(lg, 7) == [0.0, 0.0, 0.0]
eeg_tmp = locs_generate(eeg)
@test eeg_tmp.locs[1, :label] == "Fp1"
locs_generate!(eeg_tmp)
@test eeg_tmp.locs == locs_generate(eeg).locs

@info "Test: add_locs() values"
eeg_tmp = add_locs(eeg; locs = locs)
fp1 = findfirst(==("Fp1"), eeg_tmp.locs[!, :label])
@test xyz(eeg_tmp.locs, fp1) == xyz(locs, findfirst(==("Fp1"), locs[!, :label]))

@info "Test: edit_locs() values"
eeg_tmp = edit_locs(eeg; ch = "Fp1", theta = 100.0, name = "Fp1x", type = "eog")
@test labels(eeg_tmp)[1] == "Fp1x"
@test eeg_tmp.header.recording[:channel_type][1] == "eog"
@test labels(eeg)[1] == "Fp1"
@test_throws ArgumentError edit_locs(eeg; ch = "no_such_channel", x = 1.0)

@info "Test: locs_details() print"
out = mktemp() do path, io
    redirect_stdout(io) do
        locs_details(eeg; ch = "Fp1")
    end
    flush(io)
    read(path, String)
end
@test occursin("Fp1", out)

true
