using NeuroAnalyzer
using Test
using DataFrames

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
locs = import_locs(joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))

@info "Test: add_locs()"
eeg_tmp = add_locs(eeg, locs=locs)
@test nrow(eeg_tmp.locs) == 23
add_locs!(eeg, locs=locs)
@test nrow(eeg.locs) == 23

@info "Test: locs_details()"
@test locs_details(eeg, ch="Fp1", out=false) == (label = "Fp1", theta_pl = 108.0, radius_pl = 1.0, x = -0.31, y = 0.95, z = -0.03, theta_sph = 108.02, radius_sph = 1.0, phi_sph = -1.72)

@info "Test: cart2pol()"
@test cart2pol(1.0, 1.0) == (1.41, 45.0)

@info "Test: pol2cart()"
@test pol2cart(1.41, 45.0) == (1, 1)

@info "Test: sph2cart()"
@test sph2cart(1.73, 45.0, 35.27) == (1.0, 1.0, 1.0)

@info "Test: cart2sph()"
@test cart2sph(1.0, 1.0, 1.0) == (1.73, 45.0, 35.26)

@info "Test: sph2pol()"
@test sph2pol(1.73, 45.0, 35.26) == (1.41, 45.0)

@info "Test: locs_sph2cart()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_tmp = locs_sph2cart(locs_tmp)
@test locs_tmp[1, :loc_x] == 1.0

@info "Test: locs_sph2cart!()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_sph2cart!(locs_tmp)
@test locs_tmp[1, :loc_x] == 1.0

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
@test locs_details(eeg_tmp, ch="Fp1", out=false) == (label = "Fp1", theta_pl = 108.0, radius_pl = 1.0, x = -0.31, y = 0.95, z = -0.03, theta_sph = 108.02, radius_sph = 1.0, phi_sph = -1.72)
eeg_tmp = edit_locs(eeg_tmp, ch="Fp1", x=0.5)
@test locs_details(eeg_tmp, ch="Fp1", out=false) == (label = "Fp1", theta_pl = 108.0, radius_pl = 1.0, x = 0.5, y = 0.95, z = -0.03, theta_sph = 108.02, radius_sph = 1.0, phi_sph = -1.72)
edit_locs!(eeg_tmp, ch="Fp1", x=0.5)
@test locs_details(eeg_tmp, ch="Fp1", out=false) == (label = "Fp1", theta_pl = 108.0, radius_pl = 1.0, x = 0.5, y = 0.95, z = -0.03, theta_sph = 108.02, radius_sph = 1.0, phi_sph = -1.72)

@info "Test: locs_flipx()"
@test locs[1, :loc_x] == -0.31
locs2 = locs_flipx(locs)
@test locs2[1, :loc_x] == 0.31

@info "Test: locs_flipy()"
@test locs[1, :loc_y] == 0.95
locs2 = locs_flipy(locs)
@test locs2[1, :loc_y] == -0.95

@info "Test: locs_flipz()"
@test locs[1, :loc_z] == -0.03
locs2 = locs_flipz(locs)
@test locs2[1, :loc_z] == 0.03

@info "Test: locs_scale()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_radius_sph] == 1.0
locs2 = locs_scale(locs, r=1.2)
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
@test locs2[1, :loc_theta_sph] == 108.02
locs2 = locs_swapxy(locs)
@test locs2[1, :loc_theta] == 198.0
@test locs2[1, :loc_theta_sph] == -161.93

@info "Test: locs_rotx()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_theta] == 108.0
@test locs[1, :loc_theta_sph] == 108.02
@test locs[1, :loc_phi_sph] == -1.72
locs2 = locs_rotx(locs, a=20)
@test locs2[1, :loc_radius] == 1.0
@test locs2[1, :loc_theta] == 108.0
@test locs2[1, :loc_theta_sph] == 108.95
@test locs2[1, :loc_phi_sph] == 17.27

@info "Test: locs_roty()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_theta] == 108.0
@test locs[1, :loc_theta_sph] == 108.02
@test locs[1, :loc_phi_sph] == -1.72
locs2 = locs_roty(locs, a=20)
@test locs2[1, :loc_radius] == 1.0
@test locs2[1, :loc_theta] == 108.0
@test locs2[1, :loc_theta_sph] == 107.61
@test locs2[1, :loc_phi_sph] == 4.47

@info "Test: locs_rotz()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_theta] == 108.0
@test locs[1, :loc_theta_sph] == 108.02
@test locs[1, :loc_phi_sph] == -1.72
locs2 = locs_rotz(locs, a=20)
@test locs2[1, :loc_radius] == 1.0
@test locs2[1, :loc_theta] == 128.0
@test locs2[1, :loc_theta_sph] == 132.1
@test locs2[1, :loc_phi_sph] == -1.87

@info "Test: locs_generate()"
locs_tmp = locs_generate(locs_tmp)
eeg_tmp = locs_generate(eeg_tmp)
@test locs_tmp isa DataFrame
@test eeg_tmp isa NeuroAnalyzer.NEURO

@info "Test: locs_origin()"
@test locs_origin(locs, y=0.2) isa DataFrame

true