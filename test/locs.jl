using NeuroAnalyzer
using Test
using DataFrames

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
locs = import_locs(joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))

@info "test 1/21: add_locs()"
eeg_tmp = add_locs(eeg, locs=locs)
@test NeuroAnalyzer._has_locs(eeg_tmp) == true
@test nrow(eeg_tmp.locs) == 19
@test NeuroAnalyzer._has_locs(eeg) == false
add_locs!(eeg, locs=locs)
@test NeuroAnalyzer._has_locs(eeg) == true
@test nrow(eeg.locs) == 19

@info "test 2/21: locs_details()"
@test locs_details(eeg, ch=1, out=false) == (ch = 1, label = "Fp1", theta_pl = 108.0, radius_pl = 1.0, x = -0.309, y = 0.95, z = -0.035, theta_sph = 108.0, radius_sph = 1.0, phi_sph = -2.0)

@info "test 3/21: cart2pol()"
@test cart2pol(1.0, 1.0) == (1.414, 45.0)

@info "test 4/21: pol2cart()"
@test pol2cart(1.414, 45.0) == (1, 1)

@info "test 5/21: sph2cart()"
@test sph2cart(1.732, 45.0, 35.266) == (1.0, 1.0, 1.0)

@info "test 6/21: cart2sph()"
@test cart2sph(1.0, 1.0, 1.0) == (1.732, 45.0, 54.734)

@info "test 7/21: sph2pol()"
@test sph2pol(1.732, 45.0, 35.266) == (1.414, 45.0)

@info "test 8/21: locs_sph2cart()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_tmp = locs_sph2cart(locs_tmp)
@test locs_tmp[1, :loc_x] == 0.999

@info "test 9/21: locs_sph2cart!()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_sph2cart!(locs_tmp)
@test locs_tmp[1, :loc_x] == 0.999

@info "test 8/21: locs_cart2sph()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_x] = 1
locs_tmp[1, :loc_y] = 0
locs_tmp = locs_cart2sph(locs_tmp)
@test locs_tmp[1, :loc_theta_sph] == 0.0

@info "test 9/21: locs_cart2sph!()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_x] = 1
locs_tmp[1, :loc_y] = 0
locs_cart2sph!(locs_tmp)
@test locs_tmp[1, :loc_theta_sph] == 0.0

@info "test 10/21: locs_cart2pol()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_x] = 1
locs_tmp[1, :loc_y] = 0
locs_tmp = locs_cart2pol(locs_tmp)
@test locs_tmp[1, :loc_theta] == 0.0

@info "test 11/21: locs_cart2pol!()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_x] = 1
locs_tmp[1, :loc_y] = 0
locs_cart2pol!(locs_tmp)
@test locs_tmp[1, :loc_theta] == 0.0

@info "test 12/21: locs_sph2pol()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_tmp[1, :loc_phi_sph] = 0
locs_tmp = locs_sph2pol(locs_tmp)
@test locs_tmp[1, :loc_radius] == 1.0
@test locs_tmp[1, :loc_theta] == 0.0

@info "test 13/21: locs_sph2pol!()"
locs_tmp = deepcopy(locs)
locs_tmp[1, :loc_theta_sph] = 0
locs_tmp[1, :loc_phi_sph] = 0
locs_sph2pol!(locs_tmp)
@test locs_tmp[1, :loc_radius] == 1.0
@test locs_tmp[1, :loc_theta] == 0.0

@info "test 14/21: edit_electrode()"
eeg_tmp = deepcopy(eeg)
@test locs_details(eeg_tmp, ch=1, out=false) == (ch = 1, label = "Fp1", theta_pl = 108.0, radius_pl = 1.0, x = -0.309, y = 0.95, z = -0.035, theta_sph = 108.0, radius_sph = 1.0, phi_sph = -2.0)
eeg_tmp = edit_locs(eeg_tmp, ch=1, x=0.5)
@test locs_details(eeg_tmp, ch=1, out=false) == (ch = 1, label = "Fp1", theta_pl = 108.0, radius_pl = 1.0, x = 0.5, y = 0.95, z = -0.035, theta_sph = 108.0, radius_sph = 1.0, phi_sph = -2.0)
edit_locs!(eeg_tmp, ch=1, x=0.5)
@test locs_details(eeg_tmp, ch=1, out=false) == (ch = 1, label = "Fp1", theta_pl = 108.0, radius_pl = 1.0, x = 0.5, y = 0.95, z = -0.035, theta_sph = 108.0, radius_sph = 1.0, phi_sph = -2.0)

@info "test 15/21: locs_flipx()"
@test locs[1, :loc_x] == -0.309
locs2 = locs_flipx(locs)
@test locs2[1, :loc_x] == 0.309

@info "test 16/21: locs_flipy()"
@test locs[1, :loc_y] == 0.95
locs2 = locs_flipy(locs)
@test locs2[1, :loc_y] == -0.95

@info "test 17/21: locs_flipz()"
@test locs[1, :loc_z] == -0.035
locs2 = locs_flipz(locs)
@test locs2[1, :loc_z] == 0.035

@info "test 18/21: locs_scale()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_radius_sph] == 1.0
locs2 = locs_scale(locs, r=1.2)
@test locs2[1, :loc_radius] == 1.2
@test locs2[1, :loc_radius_sph] == 1.2

@info "test 19/21: locs_normalize()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_radius_sph] == 1.0
locs2 = locs_normalize(locs)
@test locs2[1, :loc_radius] == 1.0
@test locs2[1, :loc_radius_sph] == 1.0

@info "test 20/21: locs_swapxy()"
@test locs[1, :loc_theta] == 108.0
@test locs2[1, :loc_theta_sph] == 108.0
locs2 = locs_swapxy(locs)
@test locs2[1, :loc_theta] == 198.0
@test locs2[1, :loc_theta_sph] == 198.0

@info "test 21/23: locs_rotx()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_theta] == 108.0
@test locs[1, :loc_theta_sph] == 108.0
@test locs[1, :loc_phi_sph] == -2.0
locs2 = locs_rotx(locs, a=20)
@test locs2[1, :loc_radius] == 0.956
@test locs2[1, :loc_theta] == 108.858
@test locs2[1, :loc_theta_sph] == 108.858
@test locs2[1, :loc_phi_sph] == 73.02

@info "test 22/23: locs_roty()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_theta] == 108.0
@test locs[1, :loc_theta_sph] == 108.0
@test locs[1, :loc_phi_sph] == -2.0
locs2 = locs_roty(locs, a=20)
@test locs2[1, :loc_radius] == 0.997
@test locs2[1, :loc_theta] == 107.654
@test locs2[1, :loc_theta_sph] == 107.654
@test locs2[1, :loc_phi_sph] == 85.825

@info "test 23/23: locs_rotz()"
@test locs[1, :loc_radius] == 1.0
@test locs[1, :loc_theta] == 108.0
@test locs[1, :loc_theta_sph] == 108.0
@test locs[1, :loc_phi_sph] == -2.0
locs2 = locs_rotz(locs, a=20)
@test locs2[1, :loc_radius] == 0.999
@test locs2[1, :loc_theta] == 128.018
@test locs2[1, :loc_theta_sph] == 128.018
@test locs2[1, :loc_phi_sph] == 92.006

true