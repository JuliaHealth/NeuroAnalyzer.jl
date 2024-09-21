using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets
using DataFrames

ntests = 29

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a0 = zeros(2, 3, 2)

@info "Test 1/$ntests: channel_type() / set_channel_type()"
@test channel_type(e10, ch="F3") == "eeg"
e10_tmp = set_channel_type(e10, ch="F3", type="mrk")
@test channel_type(e10_tmp, ch="F3") == "mrk"
set_channel_type!(e10_tmp, ch="F3", type="eeg")
@test channel_type(e10_tmp, ch="F3") == "eeg"

@info "Test 2/$ntests: get_channel()"
@test get_channel(e10, type="eeg") == ["Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2", "F7", "F8", "T3", "T4", "T5", "T6", "Fz", "Cz", "Pz"]

@info "Test 3/$ntests: rename_channel()"
e10_tmp = rename_channel(e10, ch="Fp1", name="FP1")
@test labels(e10_tmp)[1] == "FP1"
rename_channel!(e10_tmp, ch="FP1", name="Fp1")
@test labels(e10_tmp)[1] == "Fp1"

@info "Test 4/$ntests: replace_channel()"
e10_tmp = replace_channel(e10, ch="F3", s=ones(1, epoch_len(e10), nepochs(e10)));
@test e10_tmp.data[3, :, :] == ones(epoch_len(e10), nepochs(e10))
replace_channel!(e10_tmp, ch="F3", s=zeros(1, epoch_len(e10), nepochs(e10)));
@test e10_tmp.data[3, :, :] == zeros(epoch_len(e10), nepochs(e10))

@info "Test 5/$ntests: add_labels()"
l = string.(1:24)
e10_tmp = add_label(e10, clabels=l)
@test labels(e10_tmp) == l
add_label!(e10_tmp, clabels=l)
@test labels(e10_tmp) == l

@info "Test 6/$ntests: delete_channel()"
e10_tmp = delete_channel(e10, ch="F3")
@test nchannels(e10_tmp) == 23

@info "Test 7/$ntests: delete_channel!()"
e10_copy = deepcopy(e10)
delete_channel!(e10_copy, ch="F4")
@test nchannels(e10_tmp) == nchannels(e10) - 1

@info "Test 8/$ntests: keep_channel()"
e10_tmp = keep_channel(e10, ch=["eeg"])
@test nchannels(e10_tmp) == 19

@info "Test 9/$ntests: keep_channel!()"
e10_copy = deepcopy(e10)
kept_channels = ["Fp1", "Fp2", "F3", "F4"]
keep_channel!(e10_copy, ch=kept_channels)
@test nchannels(e10_copy) == length(kept_channels)

@info "Test 10/$ntests: create_data()"
e = create_object(data_type="eda")
create_data!(e, data=rand(10, 100, 1), fs=10, type="eda")
@test size(e.data) == (10, 100, 1)

@info "Test 11/$ntests: delete_epoch()"
e10_tmp = delete_epoch(e10, ep=1)
@test nepochs(e10_tmp) == 9
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 23040 # 25600 - 2560
e10_tmp = deepcopy(e10)
delete_epoch!(e10_tmp, ep=1)
@test nepochs(e10_tmp) == 9
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 23040 # 25600 - 2560

@info "Test 12/$ntests: keep_epoch()"
e10_tmp = keep_epoch(e10, ep=1:2)
@test nepochs(e10_tmp) == 2
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 5120 # 2 × 2560
e10_tmp = deepcopy(e10)
keep_epoch!(e10_tmp, ep=1:2)
@test nepochs(e10_tmp) == 2
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 5120 # 2 × 2560

@info "Test 13/$ntests: detect_bad()"
bm, be = detect_bad(e10, ch=["Fp1", "Fp2"])
@test sum(bm) == 20
@test be == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

@info "Test 14/$ntests: epoch(), subepoch()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
@test epoch_len(e10) == 10 * sr(eeg)
e10 = epoch(eeg, ep_n=10)
@test nepochs(e10) == 10
e2 = subepoch(e10, ep_start=2, ep_end=8)
@test nepochs(e2) == 10

@info "Test 15/$ntests: epoch_time()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
@test e10.epoch_time[1] == 0.0
e10_tmp = epoch_ts(e10, ts=-1.0)
@test e10_tmp.epoch_time[1] == -1.0
epoch_ts!(e10_tmp, ts=-2.0)
@test e10_tmp.epoch_time[1] == -3.0

@info "Test 16/$ntests: extract_channel()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
s = extract_channel(e10, ch="F3")
@test size(s) == (1, 2560, 120)

@info "Test 17/$ntests: extract_epoch()"
e10_tmp = extract_epoch(e10, ep=1)
@test size(e10_tmp.data) == (24, 2560, 1)
@test length(e10_tmp.time_pts) == 2560
@test length(e10_tmp.epoch_time) == 2560

@info "Test 18/$ntests: extract_data()"
d = extract_data(e10, ch="all")
@test size(d) == (24, 2560, 120)

@info "Test 19/$ntests: extract_time()"
tpts = extract_time(e10)
@test length(tpts) == 307200

@info "Test 20/$ntests: extract_eptime()"
et = extract_eptime(e10)
@test length(et) == 2560

@info "Test 21/$ntests: trim()"
s = collect(1:100)
@test NeuroAnalyzer.trim(s, seg=(1, 10)) == 11:100
m = rand(10, 100)
@test size(NeuroAnalyzer.trim(m, seg=(1, 10))) == (10, 90)
a = rand(10, 100, 10)
@test size(NeuroAnalyzer.trim(a, seg=(1, 10))) == (10, 90, 10)
e10_tmp = NeuroAnalyzer.trim(e10, seg=(0, 11), remove_epochs=true)
@test size(e10_tmp) == (24, 2560, 118)
e10_tmp = NeuroAnalyzer.trim(e10, seg=(0, 21), remove_epochs=false)
@test size(e10_tmp) == (24, 2560, 117)

@info "Test 22/$ntests: delete_marker()"
eeg_mrk = import_edf(joinpath(testfiles_path, "eeg-test-edfplus.edf"))
@test nrow(eeg_mrk.markers) == 45
delete_marker!(eeg_mrk, n=40)
@test nrow(eeg_mrk.markers) == 44

@info "Test 23/$ntests: add_marker()"
add_marker!(eeg_mrk, id="test", start=1988, len=1, desc="test", ch=0)
@test nrow(eeg_mrk.markers) == 45

@info "Test 24/$ntests: edit_marker()"
@test eeg_mrk.markers[45, :id] == "test"
edit_marker!(eeg_mrk, n=45, id="TEST", start=1989, len=1, desc="test", ch=0)
@test eeg_mrk.markers[45, :id] == "TEST"

@info "Test 25/$ntests: channel2marker()"
@test nrow(eeg_mrk.markers) == 45
unique(eeg_mrk.data[28, :, :])
idx = getindex.(findall(eeg_mrk.data[28, :, :] .== -11502.913868619822), 1)
eeg_mrk.data[28, idx, :] .= 1.0
idx = getindex.(findall(eeg_mrk.data[28, :, :] .!= 1.0), 1)
eeg_mrk.data[28, idx, :] .= 0.0
channel2marker!(eeg_mrk, ch="Mark2", id="mrk")
@test nrow(eeg_mrk.markers) == 1094

@info "Test 26/$ntests: epoch()"
eeg_mrk2 = epoch(eeg_mrk, marker="Mark2", offset=0.2, ep_len=1.2)
@test size(eeg_mrk2) == (29, 240, 1049)
e10 = epoch(eeg, ep_len=10)
@test size(e10) == (24, 2560, 120)
e2 = subepoch(e10, ep_start=2.0, ep_end=3.996)
@test size(e2) == (24, 512, 120)

@info "Test 27/$ntests: join()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edfplus.edf"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
e10_tmp = NeuroAnalyzer.join(e10, e10)
@test size(e10_tmp) == (29, 4000, 10)

@info "Test 28/$ntests: create_object()"
@test create_object(data_type="eeg") isa NeuroAnalyzer.NEURO
@test create_object(data_type="ecog") isa NeuroAnalyzer.NEURO
@test create_object(data_type="seeg") isa NeuroAnalyzer.NEURO
@test create_object(data_type="meg") isa NeuroAnalyzer.NEURO
@test create_object(data_type="nirs") isa NeuroAnalyzer.NEURO
@test create_object(data_type="mep") isa NeuroAnalyzer.NEURO
@test create_object(data_type="sensors") isa NeuroAnalyzer.NEURO
@test create_object(data_type="eda") isa NeuroAnalyzer.NEURO

@info "Test 29/$ntests: add_channel()"
e = create_object(data_type="ecog")
add_channel!(e, data=rand(2, 10, 1), label=["ecog1", "ecog2"], type=repeat(["ecog"], 2), unit=repeat(["µV"], 2))
@test size(e.data) == (2, 10, 1)

@info "Test 30/$ntests: create_time()"
create_time!(e, fs=1)
@test e.header.recording[:sampling_rate] == 1
@test length(e.time_pts) == 10

@info "Test 31/$ntests: delete_optode()"
n = import_snirf(joinpath(testfiles_path, "fnirs-test-snirf.snirf"))
@test size(n.header.recording[:optode_pairs]) == (72, 2)
@test nrow(n.locs) == 36
@test size(n) == (72, 500, 1)
delete_optode!(n, opt=1)
@test size(n.header.recording[:optode_pairs]) == (70, 2)
@test nrow(n.locs) == 35
@test size(n) == (70, 500, 1)

true
