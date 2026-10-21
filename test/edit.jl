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

@info "Test: channel_reject() / epoch_reject() (injected faults)"
# synthetic recording: a common signal plus small independent noise in each channel
e_rej = keep_channel(crop(eeg; seg = (0, 60)); ch = "eeg")
load_locs!(e_rej; file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"))
e_rej = epoch(e_rej; ep_len = 2)
rng_ = NeuroAnalyzer.Random.MersenneTwister(7)
ch_n, ep_l, ep_n = size(e_rej.data)
t_ = collect(0:(ep_l * ep_n - 1)) ./ sr(e_rej)
common = 20 .* sin.(2π .* 10 .* t_) .+ 10 .* sin.(2π .* 3 .* t_) .+ 5 .* randn(rng_, length(t_))
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
    @test findall(channel_reject(e_bad; ch = "eeg", method = m)) == sort([ix_("F3"), ix_("C4"), ix_("P3")])
end
for m in (:p2p, :tkeo)
    @test findall(channel_reject(e_bad; ch = "eeg", method = m)) == sort([ix_("C4"), ix_("P3")])
end
@test findall(channel_reject(e_bad; ch = "eeg", method = :var)) == sort([ix_("F3"), ix_("C4"), ix_("P3")])
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
@test_throws ArgumentError epoch_reject(e_bad; ch = "eeg", nbad = 0)
@test_throws ArgumentError epoch_reject(e_bad; ch = "eeg", nbad = 20)
@test_throws ArgumentError epoch_reject(e_bad; ch = "eeg", method = :var, w = 1)
@test_throws ArgumentError channel_reject(e_bad; ch = "eeg", method = :var, z = 0)

@info "Test: channel_reject() values"
e_tmp = deepcopy(e10)
e_tmp.data[2, 100, 3] = 1e4
@test channel_reject(e_tmp; ch = ["Fp1", "Fp2", "F3"], method = :amp, amp_t = 9000, tc = 0.1) == [false, true, false]
@test !any(channel_reject(e_tmp; ch = ["Fp1", "Fp2", "F3"], method = :amp, amp_t = 9000))
e_tmp.data[5, :, :] .= 1.0
@test channel_reject(e_tmp; ch = ["Fp1", "C3"], method = :flat)[2]
channel_reject!(e_tmp; ch = ["Fp1", "Fp2", "F3"], method = :amp, amp_t = 9000, tc = 0.1)
@test e_tmp.header.recording[:bad_channel][1:3] == [false, true, false]
@test_throws ArgumentError channel_reject(e10; ch = "Fp1", method = :rmse)
@test_throws ArgumentError epoch_reject(keep_epoch(e10; ep = 1:2); ch = "Fp1", method = :var)
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
@test v[1, :, 3] ≈ e10.data[get_channel(e10; ch = "Cz")[1], :, 3] .-
                   sum(e10.data[get_channel(e10; ch = ["C3", "C4", "Fz", "Pz"]), :, 3]; dims = 1)[:] ./ 4
@test_throws ArgumentError vch(e10; f = "Fp1 +* Fp2")
# unknown channel label
@test_throws ArgumentError vch(e10; f = "Fp1 - Xyz")

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
