using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets
using DataFrames

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a0 = zeros(2, 3, 2)

@info "Test: channel_type() / set_channel_type()"
@test channel_type(e10, ch = "F3") == "eeg"
e10_tmp = set_channel_type(e10; ch = "F3", type = "mrk")
@test channel_type(e10_tmp, ch = "F3") == "mrk"
set_channel_type!(e10_tmp; ch = "F3", type = "eeg")
@test channel_type(e10_tmp, ch = "F3") == "eeg"

@info "Test: get_channel()"
# labels of the selected type, in the order of the channels in the object
@test get_channel(e10, type = "eeg") ==
      labels(e10)[e10.header.recording[:channel_type] .== "eeg"]
@test get_channel(e10, type = ["ecg", "eeg"]) ==
      labels(e10)[in.(e10.header.recording[:channel_type], Ref(["ecg", "eeg"]))]

@info "Test: rename_channel()"
e10_tmp = rename_channel(e10; ch = "Fp1", name = "FP1")
@test labels(e10_tmp)[1] == "FP1"
rename_channel!(e10_tmp; ch = "FP1", name = "Fp1")
@test labels(e10_tmp)[1] == "Fp1"

@info "Test: replace_channel()"
e10_tmp = replace_channel(e10; ch = "F3", s = ones(1, epoch_len(e10), nepochs(e10)));
@test e10_tmp.data[3, :, :] == ones(epoch_len(e10), nepochs(e10))
replace_channel!(e10_tmp; ch = "F3", s = zeros(1, epoch_len(e10), nepochs(e10)));
@test e10_tmp.data[3, :, :] == zeros(epoch_len(e10), nepochs(e10))

@info "Test: add_labels()"
l = string.(1:24)
e10_tmp = add_label(e10; clabels = l)
@test labels(e10_tmp) == l
add_label!(e10_tmp; clabels = l)
@test labels(e10_tmp) == l

@info "Test: delete_channel()"
e10_tmp = delete_channel(e10; ch = "F3")
@test nchannels(e10_tmp) == nchannels(e10) - 1
@test length(e10_tmp.header.recording[:lowpass]) == nchannels(e10_tmp)
@test_throws ArgumentError delete_channel(e10; ch = String[])

@info "Test: delete_channel!()"
e10_copy = deepcopy(e10)
delete_channel!(e10_copy; ch = "F4")
@test nchannels(e10_copy) == nchannels(e10) - 1

@info "Test: keep_channel()"
e10_tmp = keep_channel(e10; ch = ["eeg"])
@test nchannels(e10_tmp) == 19
e10_tmp = keep_channel(e10; ch = labels(e10))
@test nchannels(e10_tmp) == nchannels(e10)

@info "Test: keep_channel!()"
e10_copy = deepcopy(e10)
kept_channels = ["Fp1", "Fp2", "F3", "F4"]
keep_channel!(e10_copy; ch = kept_channels)
@test nchannels(e10_copy) == length(kept_channels)

@info "Test: create_data()"
e = create_object(; data_type = "eda")
create_data!(e; data = rand(10, 100, 1), fs = 10, type = "eda")
@test size(e.data) == (10, 100, 1)

@info "Test: delete_epoch()"
e10_tmp = delete_epoch(e10; ep = 1)
@test nepochs(e10_tmp) == 9
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 23040 # 25600 - 2560
e10_tmp = deepcopy(e10)
delete_epoch!(e10_tmp; ep = 1)
@test nepochs(e10_tmp) == 9
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 23040 # 25600 - 2560

@info "Test: keep_epoch()"
e10_tmp = keep_epoch(e10; ep = 1:2)
@test nepochs(e10_tmp) == 2
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 5120 # 2 × 2560
e10_tmp = deepcopy(e10)
keep_epoch!(e10_tmp; ep = 1:2)
@test nepochs(e10_tmp) == 2
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 5120 # 2 × 2560

@info "Test: channel_reject()"
bc = channel_reject(e10; ch = ["Fp1", "Fp2"], method = [:flat, :amp])
@test bc isa Vector{Bool}
@test length(bc) == 2

@info "Test: epoch_reject()"
be = epoch_reject(e10; ch = ["Fp1", "Fp2"], method = [:flat, :amp])
@test be isa Vector{Int64}

@info "Test: epoch(), subepoch()"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
@test epoch_len(e10) == 10 * sr(eeg)
@test nepochs(e10) == 120
e2 = subepoch(e10; ep_start = 2, ep_end = 8)
@test nepochs(e2) == 120
e10o = epoch(eeg; ep_len = 10, overlap = 5)
@test epoch_len(e10o) == 10 * sr(eeg)
@test nepochs(e10o) == 240
@test e10o.data[:, :, 2] == eeg.data[:, (5 * sr(eeg) + 1):(15 * sr(eeg)), 1]
@test e10o.data[:, 1:(5 * sr(eeg)), 2] == e10o.data[:, (5 * sr(eeg) + 1):end, 1]
@test epoch(eeg; ep_len = 10, overlap = 0).data == e10.data
@test_throws ArgumentError epoch(eeg; ep_len = 10, overlap = 10)
@test_throws ArgumentError epoch(eeg; ep_len = 10, overlap = -1)
@test_throws ArgumentError epoch(eeg; overlap = 5)

@info "Test: epoch_time()"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
@test e10.epoch_time[1] == 0.0
e10_tmp = epoch_ts(e10; ts = -1.0)
@test e10_tmp.epoch_time[1] == -1.0
epoch_ts!(e10_tmp; ts = -2.0)
@test e10_tmp.epoch_time[1] == -3.0

@info "Test: extract_channel()"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
s = extract_channel(e10; ch = "F3")
@test size(s) == (1, 2560, 120)

@info "Test: extract_epoch()"
e10_tmp = extract_epoch(e10; ep = 1)
@test size(e10_tmp.data) == (24, 2560, 1)
@test length(e10_tmp.time_pts) == 2560
@test length(e10_tmp.epoch_time) == 2560

@info "Test: extract_data()"
d = extract_data(e10; ch = "all")
@test size(d) == (24, 2560, 120)

@info "Test: trim()"
s = collect(1:100)
@test NeuroAnalyzer._trim(s, seg = (1, 10)) == 11:100
m = rand(10, 100)
@test size(NeuroAnalyzer._trim(m, seg = (1, 10))) == (10, 90)
a = rand(10, 100, 10)
@test size(NeuroAnalyzer._trim(a, seg = (1, 10))) == (10, 90, 10)
e_tmp = NeuroAnalyzer.trim(eeg; seg = (0, 11))
@test size(e_tmp) == (24, 305663, 1)
e_tmp = NeuroAnalyzer.trim(eeg; seg = (0, 21), keep = true)
@test size(e_tmp) == (24, 5377, 1)

@info "Test: crop()"
e_tmp = NeuroAnalyzer.crop(eeg; seg = (0, 21))
@test size(e_tmp) == (24, 5377, 1)

@info "Test: delete_marker()"
eeg_mrk = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edfplus.edf"))
@test nrow(eeg_mrk.markers) == 45
delete_marker!(eeg_mrk; n = 40)
@test nrow(eeg_mrk.markers) == 44

@info "Test: add_marker()"
add_marker!(eeg_mrk; type = "test", start = 1988, len = 1, value = "test")
@test nrow(eeg_mrk.markers) == 45

@info "Test: edit_marker()"
@test eeg_mrk.markers[45, :type] == "test"
edit_marker!(eeg_mrk; n = 45, type = "TEST", start = 1989, len = 1, value = "test")
@test eeg_mrk.markers[45, :type] == "TEST"

@info "Test: channel2marker()"
@test nrow(eeg_mrk.markers) == 45
unique(eeg_mrk.data[28, :, :])
# the marker level is the highest value of the Mark2 channel
idx = getindex.(findall(eeg_mrk.data[28, :, :] .== maximum(eeg_mrk.data[28, :, :])), 1)
eeg_mrk.data[28, idx, :] .= 1.0
idx = getindex.(findall(eeg_mrk.data[28, :, :] .!= 1.0), 1)
eeg_mrk.data[28, idx, :] .= 0.0
channel2marker!(eeg_mrk; ch = "Mark2", type = "stim")
@test nrow(eeg_mrk.markers) == 1094

@info "Test: epoch()"
eeg_mrk2 = epoch(eeg_mrk; marker_value = "Mark2", offset = 0.2, ep_len = 1.2)
@test size(eeg_mrk2) == (29, 240, 1049)
# each epoch contains its marker at the offset (other markers inside the epoch are kept)
m_tmp = eeg_mrk2.markers[eeg_mrk2.markers.value .== "Mark2", :]
@test any(isapprox.(m_tmp.start, 0.2; atol = 1e-9)) &&
      any(isapprox.(m_tmp.start, 1.2 + 0.2; atol = 1e-9))
@test all(0 .<= eeg_mrk2.markers.start .< 1.2 * 1049)
e10 = epoch(eeg; ep_len = 10)
@test size(e10) == (24, 2560, 120)
e2 = subepoch(e10; ep_start = 2.0, ep_end = 3.996)
@test size(e2) == (24, 512, 120)

@info "Test: join()"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edfplus.edf"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)
e10_tmp = NeuroAnalyzer.join(e10, e10)
@test size(e10_tmp) == (29, 4000, 10)

@info "Test: create_object()"
@test create_object(data_type = "eeg") isa NeuroAnalyzer.NEURO
@test create_object(data_type = "ecog") isa NeuroAnalyzer.NEURO
@test create_object(data_type = "seeg") isa NeuroAnalyzer.NEURO
@test create_object(data_type = "meg") isa NeuroAnalyzer.NEURO
@test create_object(data_type = "nirs") isa NeuroAnalyzer.NEURO
@test create_object(data_type = "mep") isa NeuroAnalyzer.NEURO
@test create_object(data_type = "sensors") isa NeuroAnalyzer.NEURO
@test create_object(data_type = "eda") isa NeuroAnalyzer.NEURO

@info "Test: add_channel()"
e = create_object(; data_type = "ecog")
e.header.recording[:sampling_rate] = 10
add_channel!(
    e;
    data = rand(2, 10, 1),
    label = ["ecog1", "ecog2"],
    type = repeat(["ecog"], 2),
    unit = repeat(["μV"], 2),
)
@test size(e.data) == (2, 10, 1)
# new channels: no filter (Nyquist frequency, 0 Hz)
@test e.header.recording[:lowpass] == [5.0, 5.0] && e.header.recording[:highpass] == [0.0, 0.0]
@test e.header.recording[:notch] == [Float64[], Float64[]]
@test e.header.recording[:notch][1] !== e.header.recording[:notch][2]

@info "Test: create_time()"
e = create_object(; data_type = "ecog")
e.header.recording[:sampling_rate] = 10
add_channel!(
    e;
    data = rand(2, 10, 1),
    label = ["ecog1", "ecog2"],
    type = repeat(["ecog"], 2),
    unit = repeat(["μV"], 2),
)
create_time!(e; fs = 5)
@test e.header.recording[:sampling_rate] == 5
@test length(e.time_pts) == 10

@info "Test: delete_optode()"
n = import_snirf(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf"))
@test size(n.header.recording[:optode_pairs]) == (72, 2)
@test nrow(n.locs) == 36
@test size(n) == (72, 500, 1)
delete_optode!(n; opt = 1)
@test size(n.header.recording[:optode_pairs]) == (70, 2)
@test nrow(n.locs) == 35
@test size(n) == (70, 500, 1)

@info "Test: delete_channel() / keep_channel() / delete_pair() (NIRS)"
n2 = import_nirs(joinpath(testfiles_path, "nirs", "fnirs-test-nirs.nirs"))
intensity2od!(n2)
od2conc!(n2)
has(o, l) = l in labels(o)
nlab(o) = length(o.header.recording[:wavelength_index])
# a single channel is deleted with the other wavelengths of its optode pair
d2 = delete_channel(n2; ch = "S1_D2 690")
@test !has(d2, "S1_D2 690") && !has(d2, "S1_D2 830") && has(d2, "S1_D2 690 OD")
@test nchannels(d2) == nchannels(n2) - 2
@test nlab(d2) == nlab(n2) - 2 && size(d2.header.recording[:optode_pairs], 1) == nlab(d2)
@test od2conc(d2) isa NeuroAnalyzer.NEURO
# keeping one wavelength keeps the other wavelengths of the pair
@test sort(labels(keep_channel(n2; ch = "S1_D1 690"))) == ["S1_D1 690", "S1_D1 830"]
# keeping only the concentration channels
hb = keep_channel(n2; ch = get_channel(n2; type = ["nirs_hbo", "nirs_hbr"]))
@test labels(hb) == ["S1_D1 HbO", "S1_D1 HbR", "S1_D2 HbO", "S1_D2 HbR"]
@test nlab(hb) == 4 && all(iszero, hb.header.recording[:wavelength_index])
# all channels of a source-detector pair; optodes are kept
dp = delete_pair(n2; pair = "S1_D2")
@test !any(startswith("S1_D2"), labels(dp))
@test nchannels(dp) == nchannels(n2) - 7
@test dp.header.recording[:optode_labels] == n2.header.recording[:optode_labels]
@test nrow(dp.locs) == nrow(n2.locs)
@test_throws ArgumentError delete_pair(n2; pair = "S1-D2")
@test_throws ArgumentError delete_pair(n2; pair = "S9_D1")
dp2 = deepcopy(n2)
delete_pair!(dp2; pair = ["S1_D2"])
@test labels(dp2) == labels(dp)
# deleting an optode removes its concentration channels too
d_opt = delete_optode(n2; opt = findfirst(==("D2"), n2.header.recording[:optode_labels]))
@test !any(startswith("S1_D2"), labels(d_opt))
@test nlab(d_opt) == size(d_opt.header.recording[:optode_pairs], 1)

##
## values: synthetic signals and marker bookkeeping
##

eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)

@info "Test: channel_reject() / epoch_reject() (injected faults)"
# synthetic recording: a common signal plus small independent noise in each channel
e_rej = keep_channel(crop(eeg; seg = (0, 60)); ch = "eeg")
load_locs!(
    e_rej;
    file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"),
)
e_rej = epoch(e_rej; ep_len = 2)
rng_ = NeuroAnalyzer.Random.MersenneTwister(7)
ch_n, ep_l, ep_n = size(e_rej.data)
t_ = collect(0:(ep_l * ep_n - 1)) ./ sr(e_rej)
common =
    20 .* sin.(2π .* 10 .* t_) .+ 10 .* sin.(2π .* 3 .* t_) .+ 5 .* randn(rng_, length(t_))
for idx = 1:ch_n
    e_rej.data[idx, :, :] = reshape(common .+ 5 .* randn(rng_, length(t_)), ep_l, ep_n)
end
l_ = labels(e_rej)
ix_(l) = findfirst(==(l), l_)
@test !any(channel_reject(e_rej; ch = "eeg"))
e_bad = deepcopy(e_rej)
e_bad.data[ix_("F3"), :, :] .= 0.0                                      # flat
e_bad.data[ix_("C4"), :, :] .+= 40 .* randn(rng_, ep_l, ep_n)            # noise
for ep_idx = 1:ep_n
    e_bad.data[ix_("P3"), 100:110, ep_idx] .+= 300                      # pops
end
@test findall(channel_reject(e_bad; ch = "eeg", method = :flat)) == [ix_("F3")]
@test findall(channel_reject(e_bad; ch = "eeg", method = :amp, amp_t = 200)) == [ix_("P3")]
# the flat channel is also far from the median signal
for m in (:rmse, :rmsd, :euclid)
    @test findall(channel_reject(e_bad; ch = "eeg", method = m)) ==
          sort([ix_("F3"), ix_("C4"), ix_("P3")])
end
for m in (:p2p, :tkeo)
    @test findall(channel_reject(e_bad; ch = "eeg", method = m)) ==
          sort([ix_("C4"), ix_("P3")])
end
@test findall(channel_reject(e_bad; ch = "eeg", method = :var)) ==
      sort([ix_("F3"), ix_("C4"), ix_("P3")])
# the flat channel is not peaked
bc = channel_reject(e_bad; ch = "eeg", method = :kurt)
@test bc[ix_("P3")] && !bc[ix_("F3")]
@test findall(channel_reject(e_bad; ch = "eeg", method = :z)) == [ix_("P3")]
@test ix_("F3") in findall(channel_reject(e_bad; ch = "eeg", method = :ransac))
@test ix_("C4") in findall(channel_reject(e_bad; ch = "eeg", method = :ransac))
e_bad = deepcopy(e_rej)
e_bad.data[ix_("C3"), :, 2] .= 0.0                                      # flat epoch
e_bad.data[ix_("O1"), 256:end, 3] .+= 200                               # jump
e_bad.data[ix_("T3"), 200:400, 4] .+= 30 .* randn(rng_, 201)            # muscle
@test epoch_reject(e_bad; ch = "eeg", method = :flat) == [2]
@test epoch_reject(e_bad; ch = "eeg", method = :amp, amp_t = 150) == [3]
@test epoch_reject(e_bad; ch = "eeg", method = :var) == [2, 3, 4]
@test epoch_reject(e_bad; ch = "eeg", method = :p2p) == [3, 4]
@test 4 in epoch_reject(e_bad; ch = "eeg", method = :tkeo)
@test epoch_reject(e_bad; ch = "eeg", method = :z) == [3, 4]
@test epoch_reject(e_bad; ch = "eeg", method = :var, nbad = 2) == Int64[]
e_rej = deepcopy(e_bad)
epoch_reject!(e_rej; ch = "eeg", method = :flat)
@test findall(e_rej.epochs.reject) == [2]
epoch_reject!(e_rej; ch = "eeg", method = :amp, amp_t = 150)
@test findall(e_rej.epochs.reject) == [2, 3]
@test !any(e_bad.epochs.reject)
@test_throws ArgumentError epoch_reject(e_bad; ch = "eeg", nbad = 0)
@test_throws ArgumentError epoch_reject(e_bad; ch = "eeg", nbad = 20)
@test_throws ArgumentError epoch_reject(e_bad; ch = "eeg", method = :var, w = 1)
@test_throws ArgumentError channel_reject(e_bad; ch = "eeg", method = :var, z = 0)

@info "Test: channel_reject() values"
e_tmp = deepcopy(e10)
e_tmp.data[2, 100, 3] = 1e4
@test channel_reject(
    e_tmp;
    ch = ["Fp1", "Fp2", "F3"],
    method = :amp,
    amp_t = 9000,
    tc = 0.1,
) == [false, true, false]
@test !any(channel_reject(e_tmp; ch = ["Fp1", "Fp2", "F3"], method = :amp, amp_t = 9000))
e_tmp.data[5, :, :] .= 1.0
@test channel_reject(e_tmp; ch = ["Fp1", "C3"], method = :flat)[2]
channel_reject!(e_tmp; ch = ["Fp1", "Fp2", "F3"], method = :amp, amp_t = 9000, tc = 0.1)
@test e_tmp.header.recording[:bad_channel][1:3] == [false, true, false]
@test_throws ArgumentError channel_reject(e10; ch = "Fp1", method = :rmse)
@test_throws ArgumentError epoch_reject(
    keep_epoch(e10; ep = 1:2);
    ch = "Fp1",
    method = :var,
)
@test_throws ArgumentError channel_reject(e10; ch = "Fp1", method = :x)

@info "Test: epoch_reject() values"
e_tmp = deepcopy(e10)
be = epoch_reject(e_tmp; ch = ["Fp1", "Fp2"], method = :amp, amp_t = 1e6)
@test isempty(be)
e_tmp.data[1, 10, 4] = 2e6
@test epoch_reject(e_tmp; ch = ["Fp1", "Fp2"], method = :amp, amp_t = 1e6) == [4]

@info "Test: vch()"
v = vch(e10; f = "Fp1 - Fp2")
@test size(v) == (1, 2560, 10)
@test v[1, :, :] == e10.data[1, :, :] .- e10.data[2, :, :]
# the formula is element-wise: mean() of a single sample is the sample itself
v = vch(e10; f = "mean(Fp1) + 0 * Cz")
@test v[1, :, :] ≈ e10.data[1, :, :]
v = vch(e10; f = "2 * F3")
@test v[1, :, :] == 2 .* e10.data[3, :, :]
v = vch(e10; f = "Cz - (C3 + C4 + Fz + Pz) / 4")
@test v[1, :, 3] ≈
      e10.data[get_channel(e10; ch = "Cz")[1], :, 3] .-
      sum(e10.data[get_channel(e10; ch = ["C3", "C4", "Fz", "Pz"]), :, 3]; dims = 1)[:] ./ 4
@test_throws ArgumentError vch(e10; f = "Fp1 +* Fp2")
# unknown channel label
@test_throws ArgumentError vch(e10; f = "Fp1 - Xyz")

@info "Test: markers"
mrk = DataFrame(
    :type => ["a", "b"],
    :start => [1.0, 12.5],
    :length => [0.0, 0.5],
    :value => ["x", "y"],
    :channel => [0, 0],
)
e_tmp = add_markers(eeg; markers = mrk)
# channel numbers are converted to labels (0: all channels = empty vector)
@test e_tmp.markers[!, 1:4] == mrk[!, 1:4]
@test e_tmp.markers[!, :channel] == [String[], String[]]
@test isempty(eeg.markers)
add_markers!(e_tmp; markers = mrk[1:1, :])
@test nrow(e_tmp.markers) == 1
@test_throws ArgumentError add_markers(eeg; markers = DataFrame(:a => [1]))
@test NeuroAnalyzer._check_marker_cols(mrk) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_marker_cols(mrk[!, 1:3])
e_tmp = add_markers(eeg; markers = mrk)
out = mktemp() do path, io
    redirect_stdout(io) do
        return view_marker(e_tmp)
    end
    flush(io)
    return read(path, String)
end
@test occursin("'x'", out) && occursin("12.5", out)
@test_throws ArgumentError view_marker(eeg)
e_tmp = add_marker(e_tmp; type = "c", start = 5.0, value = "z")
# markers are kept sorted by start
@test e_tmp.markers[!, :value] == ["x", "z", "y"]
@test_throws ArgumentError add_marker(e_tmp; type = "c", start = -1.0, value = "z")
@test_throws ArgumentError add_marker(e_tmp; type = "c", start = 1e6, value = "z")
e_tmp2 = delete_marker(e_tmp; n = 2)
@test e_tmp2.markers[!, :value] == ["x", "y"]
@test_throws ArgumentError delete_marker(e_tmp; n = 9)
e_tmp2 = edit_marker(e_tmp; n = 1, type = "a", start = 20.0, value = "x2")
@test e_tmp2.markers[!, :value] == ["z", "y", "x2"]
# marker channels: labels (empty vector = all channels)
@test e_tmp.markers[!, :channel] == [String[], String[], String[]]
e_tmp2 = add_marker(e_tmp; type = "BAD", start = 0.5, len = 1.0, value = "", ch = "Fp1")
@test e_tmp2.markers[1, :channel] == ["Fp1"]
e_tmp2 = add_marker(e_tmp; type = "BAD", start = 0.5, len = 1.0, value = "", ch = ["Fp1", "F7"])
@test e_tmp2.markers[1, :channel] == ["Fp1", "F7"]
@test_throws ArgumentError add_marker(e_tmp; type = "c", start = 2.0, value = "z", ch = "XX")
e_tmp3 = edit_marker(e_tmp2; n = 1, type = "BAD", start = 0.5, len = 1.0, value = "", ch = "O1")
@test e_tmp3.markers[1, :channel] == ["O1"]
out = mktemp() do path, io
    redirect_stdout(io) do
        return view_marker(e_tmp2)
    end
    flush(io)
    return read(path, String)
end
@test occursin("Fp1;F7", out) && occursin("all", out)
# deleting channels: removed from the markers; a marker of deleted channels only is removed
e_tmp3 = add_marker(e_tmp2; type = "BAD", start = 3.0, len = 1.0, value = "", ch = "F7")
e_del = delete_channel(e_tmp3; ch = "F7")
@test nrow(e_del.markers) == nrow(e_tmp3.markers) - 1
@test e_del.markers[1, :channel] == ["Fp1"]
@test count(isempty, e_del.markers[!, :channel]) == 3
# markers follow renamed channels; rows do not share their vectors
e_ren = rename_channel(e_tmp3; ch = "Fp1", name = "FP1")
@test e_ren.markers[1, :channel] == ["FP1", "F7"]
@test e_tmp3.markers[1, :channel] == ["Fp1", "F7"]
e_ren = add_label(e_tmp3; clabels = uppercase.(labels(e_tmp3)))
@test e_ren.markers[1, :channel] == ["FP1", "F7"]
# referencing appends a suffix to the labels
e_ref = reference_avg(e_tmp3)
@test e_ref.markers[1, :channel] == ["Fp1-avg", "F7-avg"]
# montage: markers move to the montage channels using the channel
e_mon = reference_custom(e_tmp3; ref_list = ["Fp1-F7", "F7-T3", "O1-O2"], ref_name = "bip")
@test e_mon.markers[1, :channel] == ["Fp1-F7", "F7-T3"]
@test e_mon.markers[3, :channel] == ["Fp1-F7", "F7-T3"]
@test count(isempty, e_mon.markers[!, :channel]) == 3
e_mon = reference_custom(e_tmp3; ref_list = ["O1-O2"], ref_name = "bip")
@test nrow(e_mon.markers) == 3
# add_markers(): labels separated by ";" (as written by export_markers()) and channel numbers
mrk2 = DataFrame(:type => ["a", "b", "c"], :start => [1.0, 2.0, 3.0], :length => zeros(3),
                 :value => ["x", "y", "z"], :channel => ["", "Fp1;F7", "O1"])
@test add_markers(eeg; markers = mrk2).markers[!, :channel] == [String[], ["Fp1", "F7"], ["O1"]]
mrk2[!, :channel] = [0, 1, 2]
@test add_markers(eeg; markers = mrk2).markers[!, :channel] == [String[], [labels(eeg)[1]], [labels(eeg)[2]]]
mrk2[!, :channel] = [0, 1, 999]
@test_throws ArgumentError add_markers(eeg; markers = mrk2)
# export_markers() → CSV → add_markers() round trip
f_tmp = tempname() * ".csv"
export_markers(e_tmp3; file_name = f_tmp)
m_csv = NeuroAnalyzer.CSV.read(
    f_tmp,
    DataFrame;
    types = Dict(:type => String, :value => String, :channel => String),
)
@test add_markers(eeg; markers = m_csv).markers[!, :channel] == e_tmp3.markers[!, :channel]
# import_markers(): exact round trip (also extra columns)
e_rt = deepcopy(e_tmp3)
e_rt.markers[!, :rt] = [missing, 0.5, missing, 0.25, missing]
export_markers(e_rt; file_name = f_tmp, overwrite = true)
m_imp = import_markers(f_tmp)
@test isequal(m_imp, e_rt.markers)
@test isequal(add_markers(eeg; markers = m_imp).markers, e_rt.markers)
rm(f_tmp)
# BIDS events.tsv: onset/duration in seconds, the other columns are kept
f_tsv = tempname() * "_events.tsv"
write(f_tsv, "onset\tduration\ttrial_type\tvalue\tresponse_time\n" *
             "12.5\t0\tface\t71\t0.512\n2.0\tn/a\tresponse\t1\tn/a\n")
m_bids = import_markers(f_tsv)
@test m_bids[!, :start] == [2.0, 12.5] && m_bids[!, :length] == [0.0, 0.0]
@test m_bids[!, :type] == ["response", "stim"]
@test m_bids[!, :value] == ["1", "face"]
@test m_bids[!, :bids_value] == ["1", "71"] && isequal(m_bids[!, :response_time], [missing, 0.512])
@test nrow(add_markers(eeg; markers = m_bids).markers) == 2
# events with sample numbers only (0-based) need the sampling rate
# codes keep their leading zeros
write(f_tsv, "onset\tvalue\n1.0\t001\n")
@test import_markers(f_tsv)[!, :value] == ["001"]
write(f_tsv, "sample\tvalue\n256\tA\n")
@test_throws ArgumentError import_markers(f_tsv)
@test import_markers(f_tsv; fs = 256)[!, :start] == [1.0]
rm(f_tsv)
# CSV of NeuroAnalyzer before v0.26.12: id instead of type, channel numbers
f_old = tempname() * ".csv"
write(f_old, "id,start,length,value,channel\nmrk,1.0,0.0,s1,0\nNA,2.0,0.0,DELETED,2\n")
m_old = import_markers(f_old)
@test m_old[!, :type] == ["stim", "boundary"]
@test add_markers(eeg; markers = m_old).markers[!, :channel] == [String[], [labels(eeg)[2]]]
write(f_old, "a,b\n1,2\n")
@test_throws ArgumentError import_markers(f_old)
rm(f_old)
@test_throws ArgumentError import_markers("no_such_file.csv")
# objects saved with channel numbers are converted on loading
e_old = deepcopy(e_tmp)
e_old.markers[!, :channel] = [0, 2, 0]
NeuroAnalyzer._markers_upgrade!(e_old)
@test e_old.markers[!, :channel] == [String[], [labels(eeg)[2]], String[]]
# marker types: names used by formats and programs → standard types
@test NeuroAnalyzer._marker_type("Stimulus") == "stim"
@test NeuroAnalyzer._marker_type("Response") == "response"
@test NeuroAnalyzer._marker_type("Comment") == "annotation"
@test NeuroAnalyzer._marker_type("New Segment") == "boundary"
@test NeuroAnalyzer._marker_type("BAD boundary") == "boundary"
@test NeuroAnalyzer._marker_type("BAD_ACQ_SKIP") == "BAD"
@test NeuroAnalyzer._marker_type("Badania") == "other"
@test NeuroAnalyzer._marker_type("SyncStatus"; default = "stim") == "stim"
@test NeuroAnalyzer._gdf_marker_type(0x0101) == "BAD"
@test NeuroAnalyzer._gdf_marker_type(0x7ffe) == "boundary"
@test NeuroAnalyzer._gdf_marker_type(0x0301) == "stim"
@test NeuroAnalyzer._gdf_marker_type(0x8301) == "stim"
@test NeuroAnalyzer._gdf_marker_type(0x0410) == "annotation"
# objects and tables with the `id` column (before v0.26.12) are converted
e_old = deepcopy(e_tmp)
DataFrames.rename!(e_old.markers, :type => :id)
e_old.markers[!, :id] = ["NA", "mrk", "3"]
e_old.markers[1, :value] = "DELETED"
e_old.markers[!, :channel] = [0, 0, 0]
NeuroAnalyzer._markers_upgrade!(e_old)
@test names(e_old.markers) == ["type", "start", "length", "value", "channel"]
@test e_old.markers[!, :type] == ["boundary", "stim", "stim"]
mrk_old = DataFrame(:id => ["Stimulus", "x"], :start => [1.0, 2.0], :length => zeros(2),
                    :value => ["S  1", "y"], :channel => [0, 0])
@test add_markers(eeg; markers = mrk_old).markers[!, :type] == ["stim", "other"]
# trim(): the cut is marked by a boundary marker
@test NeuroAnalyzer.trim(e_tmp; seg = (4, 6)).markers[!, :type] == ["a", "boundary", "b"]
# epoch(): markers can be selected by type
e_ep = add_marker(eeg; type = "stim", start = 10.0, value = "A")
e_ep = add_marker(e_ep; type = "stim", start = 30.0, value = "A")
e_ep = add_marker(e_ep; type = "annotation", start = 50.0, value = "A")
@test nepochs(epoch(e_ep; marker_value = "A", ep_len = 2)) == 3
@test nepochs(epoch(e_ep; marker_value = "A", marker_type = "stim", ep_len = 2)) == 2
@test_throws ArgumentError epoch(e_ep; marker_value = "A", marker_type = "response", ep_len = 2)
# epoching by several markers at once: condition labels in the epoch metadata
e_mc = add_marker(eeg; type = "stim", start = 10.0, value = "s70")
add_marker!(e_mc; type = "stim", start = 20.0, value = "s71")
add_marker!(e_mc; type = "stim", start = 30.0, value = "s130")
add_marker!(e_mc; type = "stim", start = 40.0, value = "s71")
e_mc.markers[!, :trial_type] = ["face", "face", "house", "face"]
ep_mc = epoch(e_mc; marker_value = ["s70", "s71"], offset = 0.5, ep_len = 2)
@test nepochs(ep_mc) == 3 && ep_mc.epochs.condition == ["s70", "s71", "s71"]
@test ep_mc.epochs.onset ≈ [9.5, 19.5, 39.5]
cond_map = Dict("s70" => "face", "s71" => "face", "s130" => "house")
ep_mc = epoch(e_mc; marker_value = r"^s", offset = 0.5, ep_len = 2, condition = cond_map)
@test ep_mc.epochs.condition == ["face", "face", "house", "face"]
@test epoch_conditions(ep_mc) == ["face", "house"]
@test get_epochs(ep_mc; condition = "face") == [1, 2, 4]
@test get_epochs(ep_mc; condition = r"ou") == [3]
@test get_epochs(ep_mc; condition = ["house", "car"]) == [3]
@test get_epochs(ep_mc) == 1:4
# condition from a marker column (e.g. trial_type of BIDS events)
@test epoch(e_mc; marker_value = r"^s", offset = 0.5, ep_len = 2,
            condition = "trial_type").epochs.condition == ep_mc.epochs.condition
@test_throws ArgumentError epoch(e_mc; marker_value = r"^s", ep_len = 2, condition = "xx")
@test_throws ArgumentError epoch(e_mc; marker_value = String[], ep_len = 2)
@test_throws ArgumentError epoch(e_mc; marker_value = ["s1", "s2"], ep_len = 2)
@test isempty(epoch_conditions(epoch(eeg; ep_len = 10)))
# rejected epochs are excluded by get_epochs() and average_epochs()
ep_mc.epochs.reject[2] = true
@test get_epochs(ep_mc; condition = "face") == [1, 4]
@test get_epochs(ep_mc; condition = "face", rejected = true) == [1, 2, 4]
a_mc = average_epochs(ep_mc; condition = "face")
@test nepochs(a_mc) == 3 && a_mc.epochs.condition == ["face", "face", "face"]
@test a_mc.data[:, :, 1] ≈ (a_mc.data[:, :, 2] .+ a_mc.data[:, :, 3]) ./ 2
@test nepochs(average_epochs(ep_mc)) == 4
@test nepochs(average_epochs(ep_mc; rejected = true)) == 5
@test_throws ArgumentError average_epochs(ep_mc; condition = "car")
# statistics: conditions compared within one object (paired tests need equal numbers of epochs)
md = mdiff(ep_mc, ep_mc; ch1 = "Fp1", ch2 = "Fp1",
           ep1 = get_epochs(ep_mc; condition = "face")[1:1],
           ep2 = get_epochs(ep_mc; condition = "house"), n = 10)
@test all(0 .<= md.p .<= 1)
# extra marker columns are kept by all marker functions (missing in markers added by NeuroAnalyzer)
mrk_x = DataFrame(:type => ["stim", "stim"], :start => [10.0, 30.0], :length => zeros(2),
                  :value => ["A", "B"], :channel => [String[], String[]],
                  :trial_type => ["face", "house"], :response_time => [0.512, 0.634])
e_x = add_markers(eeg; markers = mrk_x)
@test names(e_x.markers)[6:7] == ["trial_type", "response_time"]
add_marker!(e_x; type = "BAD", start = 20.0, len = 1.0, value = "", ch = "Fp1")
@test nrow(e_x.markers) == 3 && ismissing(e_x.markers[2, :trial_type])
@test e_x.markers[[1, 3], :response_time] == [0.512, 0.634]
e_x2 = edit_marker(e_x; n = 1, type = "stim", start = 11.0, value = "A2")
@test e_x2.markers[1, :trial_type] == "face"
@test delete_marker(e_x; n = 2).markers[!, :trial_type] == ["face", "house"]
e_x2 = NeuroAnalyzer.trim(e_x; seg = (15, 16))
@test "boundary" in e_x2.markers[!, :type] && "response_time" in names(e_x2.markers)
@test count(ismissing, e_x2.markers[!, :trial_type]) == 2
@test nrow(delete_channel(e_x; ch = "Fp1").markers) == 2
@test "trial_type" in names(rename_channel(e_x; ch = "Fp1", name = "FP1").markers)
@test epoch(e_x; marker_value = "B", offset = 1, ep_len = 2).markers[1, :response_time] == 0.634
e_j = NeuroAnalyzer.join(e_x, eeg)
@test "trial_type" in names(e_j.markers) && nrow(e_j.markers) == 3
out = mktemp() do path, io
    redirect_stdout(io) do
        return view_marker(e_x)
    end
    flush(io)
    return read(path, String)
end
@test occursin("trial_type=face", out) && occursin("response_time=0.634", out)
@test_throws ArgumentError add_markers(eeg; markers = mrk_x[!, [:start, :type, :length, :value, :channel]])
# channel → markers: rising and falling edges of a pulse channel
e_tmp = deepcopy(eeg)
e_tmp.data[24, :, 1] .= 0.0
e_tmp.data[24, 257:512, 1] .= 1.0
e_tmp.data[24, 2561:2600, 1] .= 1.0
e_tmp = channel2marker(e_tmp; ch = labels(eeg)[24], value = "pulse")
@test e_tmp.markers[!, :start] ≈ [1.0, 10.0]
@test e_tmp.markers[!, :length] ≈ [1.0, 40 / 256]
@test e_tmp.markers[!, :value] == ["pulse", "pulse"]
@test_throws ArgumentError channel2marker(e_tmp; ch = labels(eeg)[24], v = 7)
# trigger channel with several event codes: an event may end by a change to another code
e_tmp.data[24, :, 1] .= 0.0
e_tmp.data[24, 257:300, 1] .= 1.0
e_tmp.data[24, 301:350, 1] .= 32.0
e_tmp.data[24, 1025:1100, 1] .= 1.0
e_tmp.markers = e_tmp.markers[1:0, :]
e_tmp2 = channel2marker(e_tmp; ch = labels(eeg)[24], v = 1, value = "code 1")
@test e_tmp2.markers[!, :start] ≈ [1.0, 4.0]
@test e_tmp2.markers[!, :length] ≈ [44 / 256, 76 / 256]
e_tmp2 = channel2marker(e_tmp; ch = labels(eeg)[24], v = 32, value = "code 32")
@test e_tmp2.markers[!, :start] ≈ [300 / 256]

@info "Test: reflect() / chop()"
e_tmp = add_markers(
    e10;
    markers = DataFrame(
        :type => ["a"],
        :start => [12.0],
        :length => [0.0],
        :value => ["x"],
        :channel => [0],
    ),
)
er = NeuroAnalyzer.reflect(e_tmp; n = 256)
@test epoch_len(er) == 2560 + 2 * 256
# reflected edges: mirrored first / last n samples
@test er.data[1, 1:256, 1] == reverse(e10.data[1, 1:256, 1])
@test er.data[1, 257:2816, 1] == e10.data[1, :, 1]
# the marker at 12 s (epoch 2, 2 s into the epoch) moves with its epoch
@test er.markers[1, :start] ≈ (2560 + 512) / 256 + 1 + 2.0
ec = NeuroAnalyzer.chop(er; n = 256)
@test ec.data == e10.data
@test ec.markers[!, :start] ≈ [12.0]
er2 = deepcopy(e_tmp)
NeuroAnalyzer.reflect!(er2; n = 256)
@test er2.data == er.data
NeuroAnalyzer.chop!(er2; n = 256)
@test er2.data == e10.data
@test_throws ArgumentError NeuroAnalyzer.chop(e10; n = 2000)

@info "Test: join() values"
e1 = keep_epoch(e10; ep = 1:2)
j = NeuroAnalyzer.join(e1, e1)
@test j.data[:, 1:2560, 1] == e1.data[:, :, 1]
@test j.data[:, 2561:5120, 1] == e1.data[:, :, 1]
j2 = deepcopy(e1)
NeuroAnalyzer.join!(j2, e1)
@test j2.data == j.data
@test_throws ArgumentError NeuroAnalyzer.join(e1, keep_epoch(e10; ep = 1:3))
@test_throws ArgumentError NeuroAnalyzer.join(
    e1,
    rename_channel(e1; ch = "Fp1", name = "X"),
)

@info "Test: subepoch() / extract_epoch() / trim() / crop() in-place variants"
e_tmp = deepcopy(e10)
subepoch!(e_tmp; ep_start = 2.0, ep_end = 3.996)
@test e_tmp.data == subepoch(e10; ep_start = 2.0, ep_end = 3.996).data
@test e_tmp.data[:, 1, 1] == e10.data[:, 513, 1]
@test_throws ArgumentError subepoch(e10; ep_start = 3.0, ep_end = 2.0)
e_tmp = deepcopy(e10)
extract_epoch!(e_tmp; ep = 3)
@test e_tmp.data[:, :, 1] == e10.data[:, :, 3]
@test_throws ArgumentError extract_epoch(e10; ep = 11)
e_tmp = deepcopy(eeg)
NeuroAnalyzer.trim!(e_tmp; seg = (0, 11))
@test e_tmp.data == NeuroAnalyzer.trim(eeg; seg = (0, 11)).data
@test e_tmp.data[:, 1, 1] == eeg.data[:, 11 * 256 + 2, 1]
e_tmp = deepcopy(eeg)
NeuroAnalyzer.crop!(e_tmp; seg = (1, 2))
@test e_tmp.data[:, :, 1] == eeg.data[:, 257:513, 1]
@test_throws ArgumentError NeuroAnalyzer.crop(e10; seg = (1, 2))
@test_throws ArgumentError NeuroAnalyzer.trim!(e10; seg = (1, 2))
# trim removes the markers in the segment and shifts the later ones
e_tmp = add_markers(
    eeg;
    markers = DataFrame(
        :type => ["a", "b"],
        :start => [5.0, 20.0],
        :length => [0.0, 0.0],
        :value => ["x", "y"],
        :channel => [0, 0],
    ),
)
e_tmp = NeuroAnalyzer.trim(e_tmp; seg = (4, 6))
# a "DELETED" boundary marker is placed at the cut; [4 s, 6 s] = 513 samples are removed
@test e_tmp.markers[!, :value] == ["DELETED", "y"]
@test e_tmp.markers[1, :start] == 4.0
@test e_tmp.markers[2, :start] ≈ 20.0 - 513 / 256

@info "Test: channel helpers"
@test NeuroAnalyzer._single_channel(e10, "Fp2") == 2
@test_throws ArgumentError NeuroAnalyzer._single_channel(e10, r"Fp" |> string)
@test !NeuroAnalyzer._has_ssp(e10)
lr = DataFrame(:label => ["a", "b", "a"])
NeuroAnalyzer._rename_locs!(lr, "a", "c")
@test lr[!, :label] == ["c", "b", "c"]
rec = Dict{Symbol, Any}(:magnetometers => [1, 4], :gradiometers => [2, 3])
NeuroAnalyzer._update_meg_sensor_lists!(rec, 4, "grad")
@test rec[:magnetometers] == [1] && rec[:gradiometers] == [2, 3, 4]
NeuroAnalyzer._update_meg_sensor_lists!(rec, 2, "eeg")
@test rec[:gradiometers] == [3, 4]
e_tmp = NeuroAnalyzer._delete_channel(e10, [2, 3])
@test labels(e_tmp) == labels(e10)[[1; 4:24]]
@test e_tmp.data == e10.data[[1; 4:24], :, :]
@test_throws ArgumentError NeuroAnalyzer._delete_channel(e10, collect(1:24))
@test_throws ArgumentError NeuroAnalyzer._delete_channel(e10, [30])

@info "Test: delete/keep channel values"
e_tmp = delete_channel(e10; ch = ["F3", "Fp1"])
@test labels(e_tmp)[1:2] == ["Fp2", "F4"]
@test e_tmp.data[1, :, :] == e10.data[2, :, :]
e_tmp = keep_channel(e10; ch = ["F3", "Fp1"])
@test labels(e_tmp) == ["Fp1", "F3"]
@test e_tmp.data == e10.data[[1, 3], :, :]

@info "Test: delete_epoch() values"
e_tmp = delete_epoch(e10; ep = [2, 5])
@test e_tmp.data == e10.data[:, :, [1, 3, 4, 6, 7, 8, 9, 10]]
@test_throws ArgumentError delete_epoch(e10; ep = 1:10)
@test_throws ArgumentError keep_epoch(e10; ep = 11)

@info "Test: extract_channel() / extract_data() values"
@test extract_channel(e10; ch = "F3") == e10.data[3:3, :, :]
@test extract_data(e10; ch = ["Fp1", "F3"], ep = 2) == e10.data[[1, 3], :, 2:2]

@info "Test: epoch_ts() values"
e_tmp = epoch_ts(e10; ts = -0.5)
@test e_tmp.epoch_time ≈ e10.epoch_time .- 0.5

@info "Test: remove_pops()"
# synthetic: 10 Hz sine + noise (σ = 2), a positive pop (200, τ = 40 samples) and a negative pop (-150, τ = 20 samples)
fs_p = 256
n_p = 10 * fs_p
t_p = (0:(n_p - 1)) ./ fs_p
clean_p = 10 .* sin.(2π * 10 .* t_p) .+ 2 .* randn(NeuroAnalyzer.Random.Xoshiro(1), n_p)
pops_p = zeros(n_p)
pops_p[1200:end] .+= 200 .* exp.(-(0:(n_p - 1200)) ./ 40)
pops_p[1900:end] .-= 150 .* exp.(-(0:(n_p - 1900)) ./ 20)
x = clean_p .+ pops_p
x0 = copy(x)
p = NeuroAnalyzer._remove_pops(x; r = 128)
@test x == x0
@test p.pop_loc == [1200, 1900]
@test isapprox(p.amp, [200, -150]; rtol = 0.05)
@test isapprox(p.tau, [40, 20]; rtol = 0.1)
# the pops are removed, the underlying signal is preserved (residual below the noise level)
@test maximum(abs.(p.s .- clean_p)) < 10
@test sqrt(sum(abs2, p.s .- clean_p) / n_p) < 1
# detection only
p = NeuroAnalyzer._remove_pops(x; r = 128, repair = false)
@test p.s == x && p.pop_loc == [1200, 1900]
# no pops in a clean signal; a single-sample spike is not a pop
p = NeuroAnalyzer._remove_pops(clean_p; r = 128)
@test isempty(p.pop_loc) && p.s == clean_p
x_sp = copy(clean_p)
x_sp[1500] += 300
@test isempty(NeuroAnalyzer._remove_pops(x_sp; r = 128).pop_loc)
# exponential fit
@test all(
    isapprox.(
        values(NeuroAnalyzer._pop_fit(50 .* exp.(-(0:99) ./ 12))),
        (50, 12);
        rtol = 0.02,
    ),
)
# a step that does not recover is not a pop: not repaired (a warning is shown)
x_st = copy(clean_p)
x_st[1200:end] .+= 300
p = NeuroAnalyzer._remove_pops(x_st; r = 128)
@test isempty(p.pop_loc) && p.s == x_st
# acceptance: a small pop (relative to the signal SD) or a shape that does not follow the model is rejected
@test isempty(NeuroAnalyzer._remove_pops(x; r = 128, a_min = 100).pop_loc)
x_bl = copy(clean_p)
x_bl[1200:1327] .+= 300 .* sin.(π .* (0:127) ./ 127) .^ 8   # a fast bump that returns (not a step with exponential recovery)
@test isempty(NeuroAnalyzer._remove_pops(x_bl; r = 128).pop_loc)
# a pop at the very beginning of the signal (e.g. a filter transient) is not analyzed
x_b = copy(clean_p)
x_b[2:end] .+= 300 .* exp.(-(0:(n_p - 2)) ./ 40)
@test isempty(NeuroAnalyzer._remove_pops(x_b; r = 128).pop_loc)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x; r = 128, a_min = -1)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x; r = 128, r2_min = 2)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x; r = 4)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x; r = 128, z = 0)
@test_throws ArgumentError NeuroAnalyzer._remove_pops(x[1:100]; r = 128)
# NEURO object: channel 1 with the pops, channel 2 clean
o_p = create_object(; data_type = "eeg")
o_p.header.recording[:sampling_rate] = fs_p
add_channel!(
    o_p;
    data = reshape(vcat(x', clean_p'), 2, n_p, 1),
    label = ["P1", "P2"],
    type = ["eeg", "eeg"],
    unit = ["μV", "μV"],
)
create_time!(o_p; fs = fs_p)
o_r, pops = remove_pops(o_p; ch = "all", window = 5)
@test pops[!, :channel] == ["P1", "P1"]
@test pops[!, :sample] == [1200, 1900]
@test pops[!, :time] ≈ [1199, 1899] ./ fs_p
@test isapprox(pops[!, :tau], [40, 20] ./ fs_p; rtol = 0.1)
@test o_r.data[2, :, 1] == clean_p
@test sqrt(sum(abs2, o_r.data[1, :, 1] .- clean_p) / n_p) < 1
@test o_p.data[1, :, 1] == x
# detection only: same pops (amplitudes of later pops differ: measured from the unrepaired level)
@test remove_pops(o_p; ch = "all", repair = false)[!, [:channel, :sample]] ==
      pops[!, [:channel, :sample]]
@test_throws ArgumentError remove_pops(o_p; ch = "all", window = 0.5)
@test_throws ArgumentError remove_pops(o_p; ch = "all", window = 20)
@test_throws ArgumentError remove_pops(epoch(o_p; ep_len = 5); ch = "all")
# in-place
rp_out = remove_pops(o_p; ch = "all")
o_tmp = deepcopy(o_p)
@test remove_pops!(o_tmp; ch = "all") == rp_out[2]
@test o_tmp.data == rp_out[1].data

true
