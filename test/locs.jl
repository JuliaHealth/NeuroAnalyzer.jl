using NeuroAnalyzer
using Test

eeg2 = edit_electrode(eeg, channel=1, x=2)
@test eeg2.locs[!, :loc_x][1] == 2.0
_, _, x, _, _, _, _, _ = electrode_loc(eeg2, channel=1, output=false)
@test x == 2.0

ch1 = electrode_loc(eeg, channel=1, output=false)
@test ch1[1] == 108.0

locs = locs_import_ced("../locs/standard-10-20-cap19-elmiko.ced")
locs2 = locs_flipx(locs)
@test locs2[1, 3] == 72.0
locs2 = locs_flipy(locs)
@test locs2[1, 3] == 252.0
locs2 = locs_flipz(locs)
@test locs2[1, 3] == 108.0
locs2 = locs_swapxy(locs)
@test locs2[1, 3] == 198.0
locs2 = locs_sph2cart(locs)
@test locs2[1, 5] == -0.309
locs2 = locs_cart2sph(locs)
@test locs2[1, 3] == 108.0
locs2 = locs_maximize(locs)
@test locs2[1, :loc_radius] == 1.0

true