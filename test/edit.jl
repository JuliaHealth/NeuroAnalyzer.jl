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
@test get_channel(e10, type = "eeg") == [
    "C3",
    "C4",
    "Cz",
    "F3",
    "F4",
    "F7",
    "F8",
    "Fp1",
    "Fp2",
    "Fz",
    "O1",
    "O2",
    "P3",
    "P4",
    "Pz",
    "T3",
    "T4",
    "T5",
    "T6",
]

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
bc = channel_reject(e10; ch = ["Fp1", "Fp2"])
@test sum(bc) == 2

@info "Test: epoch_reject()"
be = epoch_reject(e10; ch = ["Fp1", "Fp2"])
@test be == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

@info "Test: epoch(), subepoch()"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
@test epoch_len(e10) == 10 * sr(eeg)
@test nepochs(e10) == 120
e2 = subepoch(e10; ep_start = 2, ep_end = 8)
@test nepochs(e2) == 120

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
add_marker!(eeg_mrk; id = "test", start = 1988, len = 1, value = "test", ch = 0)
@test nrow(eeg_mrk.markers) == 45

@info "Test: edit_marker()"
@test eeg_mrk.markers[45, :id] == "test"
edit_marker!(eeg_mrk; n = 45, id = "TEST", start = 1989, len = 1, value = "test", ch = 0)
@test eeg_mrk.markers[45, :id] == "TEST"

@info "Test: channel2marker()"
@test nrow(eeg_mrk.markers) == 45
unique(eeg_mrk.data[28, :, :])
# the marker level is the highest value of the Mark2 channel
idx = getindex.(findall(eeg_mrk.data[28, :, :] .== maximum(eeg_mrk.data[28, :, :])), 1)
eeg_mrk.data[28, idx, :] .= 1.0
idx = getindex.(findall(eeg_mrk.data[28, :, :] .!= 1.0), 1)
eeg_mrk.data[28, idx, :] .= 0.0
channel2marker!(eeg_mrk; ch = "Mark2", id = "mrk")
@test nrow(eeg_mrk.markers) == 1094

@info "Test: epoch()"
eeg_mrk2 = epoch(eeg_mrk; marker = "Mark2", offset = 0.2, ep_len = 1.2)
@test size(eeg_mrk2) == (29, 240, 1049)
# each epoch contains its marker at the offset (other markers inside the epoch are kept)
m_tmp = eeg_mrk2.markers[eeg_mrk2.markers.value .== "Mark2", :]
@test any(isapprox.(m_tmp.start, 0.2; atol = 1e-9)) && any(isapprox.(m_tmp.start, 1.2 + 0.2; atol = 1e-9))
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


##
## values: synthetic signals and marker bookkeeping
##

eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)

@info "Test: artifact detectors (synthetic signals)"
rng_ = NeuroAnalyzer.Random.MersenneTwister(7)
s = randn(rng_, 8, 2000) .* 10
# amplitude: channel 3 exceeds ±400
s_a = copy(s)
s_a[3, 1000] = 500
@test NeuroAnalyzer.detect_amp(s_a) == [false, false, true, false, false, false, false, false]
@test NeuroAnalyzer.detect_amp(s_a; amp_t = 600) == falses(8)
# distance from the cross-channel median: channel 5 is 20× larger
s_d = copy(s)
s_d[5, :] .*= 20
for f in (NeuroAnalyzer.detect_rmse, NeuroAnalyzer.detect_rmsd, NeuroAnalyzer.detect_euclid)
    @test f(s_d)[5]
    @test count(f(s_d)) <= 2
end
# peak-to-peak: a step of 2000 μV in channel 2; DC offsets do not matter
s_p = copy(s)
s_p[2, 1001:end] .+= 2000
@test NeuroAnalyzer.detect_p2p(s_p; w = 10)[2]
@test !any(NeuroAnalyzer.detect_p2p(s_p; w = 10)[[1, 3, 4, 5, 6, 7, 8]])
@test !any(NeuroAnalyzer.detect_p2p(s .+ 1000; w = 10))
@test_throws ArgumentError NeuroAnalyzer.detect_p2p(s; w = 5000)
# TKEO: bursts of high-frequency activity
t_ = collect(0:(1 / 256):(2000 / 256))[1:2000]
@test NeuroAnalyzer.detect_tkeo(s, t_) isa Vector{Bool}
@test length(NeuroAnalyzer.detect_tkeo(s, t_)) == 8
# RANSAC: channels 1-4 share a signal with their nearest neighbors, channel 4's neighbor is noise
base = randn(rng_, 2000)
s_r = vcat(base', base' .+ 0.05 .* randn(rng_, 1, 2000), randn(rng_, 1, 2000), base' .+ 0.05 .* randn(rng_, 1, 2000))
bad = NeuroAnalyzer.detect_ransac(s_r; loc_x = [0.0, 0.1, 1.0, 0.2], loc_y = [0.0, 0.0, 1.0, 0.0], w = 100)
@test bad[3]
@test !bad[1] && !bad[2]

@info "Test: channel_reject() values"
e_tmp = deepcopy(e10)
e_tmp.data[2, 100, 3] = 1e4
@test channel_reject(e_tmp; ch = ["Fp1", "Fp2", "F3"], method = :amp, amp_t = 9000) == [false, true, false]
e_tmp.data[5, :, :] .= 1.0
@test channel_reject(e_tmp; ch = ["Fp1", "C3"], method = :flat)[2]
channel_reject!(e_tmp; ch = ["Fp1", "Fp2", "F3"], method = :amp, amp_t = 9000)
@test e_tmp.header.recording[:bad_channel][1:3] == [false, true, false]
@test_throws ArgumentError channel_reject(e10; ch = "Fp1", method = :rmse)
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
v = vch(e10; f = "mean(Fp1) + 0 * Cz")
@test v[1, 5, 2] ≈ sum(e10.data[1, :, 2]) / 2560 || v[1, 5, 2] ≈ e10.data[1, 5, 2]
v = vch(e10; f = "2 * F3")
@test v[1, :, :] == 2 .* e10.data[3, :, :]
@test_throws ArgumentError vch(e10; f = "Fp1 +* Fp2")

@info "Test: markers"
mrk = DataFrame(:id => ["a", "b"], :start => [1.0, 12.5], :length => [0.0, 0.5], :value => ["x", "y"], :channel => [0, 0])
e_tmp = add_markers(eeg; markers = mrk)
@test e_tmp.markers == mrk
@test isempty(eeg.markers)
add_markers!(e_tmp; markers = mrk[1:1, :])
@test nrow(e_tmp.markers) == 1
@test_throws ArgumentError add_markers(eeg; markers = DataFrame(:a => [1]))
@test NeuroAnalyzer._check_marker_cols(mrk) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_marker_cols(mrk[!, 1:3])
e_tmp = add_markers(eeg; markers = mrk)
out = mktemp() do path, io
    redirect_stdout(io) do
        view_marker(e_tmp)
    end
    flush(io)
    read(path, String)
end
@test occursin("'x'", out) && occursin("12.5", out)
@test_throws ArgumentError view_marker(eeg)
e_tmp = add_marker(e_tmp; id = "c", start = 5.0, value = "z")
# markers are kept sorted by start
@test e_tmp.markers[!, :value] == ["x", "z", "y"]
@test_throws ArgumentError add_marker(e_tmp; id = "c", start = -1.0, value = "z")
@test_throws ArgumentError add_marker(e_tmp; id = "c", start = 1e6, value = "z")
e_tmp2 = delete_marker(e_tmp; n = 2)
@test e_tmp2.markers[!, :value] == ["x", "y"]
@test_throws ArgumentError delete_marker(e_tmp; n = 9)
e_tmp2 = edit_marker(e_tmp; n = 1, id = "a", start = 20.0, value = "x2")
@test e_tmp2.markers[!, :value] == ["z", "y", "x2"]
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
e_tmp = add_markers(e10; markers = DataFrame(:id => ["a"], :start => [12.0], :length => [0.0], :value => ["x"], :channel => [0]))
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
@test_throws ArgumentError NeuroAnalyzer.join(e1, rename_channel(e1; ch = "Fp1", name = "X"))

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
e_tmp = add_markers(eeg; markers = DataFrame(:id => ["a", "b"], :start => [5.0, 20.0], :length => [0.0, 0.0], :value => ["x", "y"], :channel => [0, 0]))
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

true
