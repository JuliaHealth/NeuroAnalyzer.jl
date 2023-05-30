using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets
using DataFrames

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

@info "test 1/26: channel_type()"
@test e10.header.recording[:channel_type][1] == "eeg"
e10_tmp = channel_type(e10, ch=1, type="???")
@test e10_tmp.header.recording[:channel_type][1] == "???"
channel_type!(e10_tmp, ch=1, type="eeg")
@test e10_tmp.header.recording[:channel_type][1] == "eeg"

@info "test 2/26: get_channel()"
@test get_channel(e10, ch=1) == "Fp1"
@test get_channel(e10, ch="Fp1") == 1

@info "test 3/26: rename_channel()"
e10_tmp = rename_channel(e10, ch=1, name="FP1")
@test get_channel(e10_tmp, ch=1) == "FP1"
rename_channel!(e10_tmp, ch=1, name="Fp1")
@test get_channel(e10_tmp, ch=1) == "Fp1"

@info "test 4/26: replace_channel()"
e10_tmp = replace_channel(e10, ch=1, s=ones(1, epoch_len(e10), epoch_n(e10)));
@test e10_tmp.data[1, :, :] == ones(epoch_len(e10), epoch_n(e10))
replace_channel!(e10_tmp, ch=1, s=zeros(1, epoch_len(e10), epoch_n(e10)));
@test e10_tmp.data[1, :, :] == zeros(epoch_len(e10), epoch_n(e10))

@info "test 5/26: add_labels()"
l = string.(1:24)
e10_tmp = add_labels(e10, clabels=l)
@test labels(e10_tmp) == l
add_labels!(e10_tmp, clabels=l)
@test labels(e10_tmp) == l

@info "test 6/26: add_labels()"
e10_tmp = delete_channel(e10, ch=1)
@test channel_n(e10_tmp) == 23
delete_channel!(e10_tmp, ch=1)
@test channel_n(e10_tmp) == 22

@info "test 7/26: keep_channel()"
e10_tmp = keep_channel(e10, ch=10:24)
@test channel_n(e10_tmp) == 15
keep_channel!(e10_tmp, ch=5:15)
@test channel_n(e10_tmp) == 11

@info "test 8/26: keep_channel_type()"
e10_tmp = keep_channel_type(e10, type=:eog)
@test channel_n(e10_tmp) == 2
e10_tmp = deepcopy(e10)
keep_channel_type!(e10_tmp, type=:eog)
@test channel_n(e10_tmp) == 2

@info "test 9/26: delete_epoch()"
e10_tmp = delete_epoch(e10, ep=1)
@test epoch_n(e10_tmp) == 9
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 23040 # 25600 - 2560
e10_tmp = deepcopy(e10)
delete_epoch!(e10_tmp, ep=1)
@test epoch_n(e10_tmp) == 9
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 23040 # 25600 - 2560

@info "test 10/26: keep_epoch()"
e10_tmp = keep_epoch(e10, ep=1:2)
@test epoch_n(e10_tmp) == 2
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 5120 # 2 × 2560
e10_tmp = deepcopy(e10)
keep_epoch!(e10_tmp, ep=1:2)
@test epoch_n(e10_tmp) == 2
@test length(e10.time_pts) == 25600
@test length(e10_tmp.time_pts) == 5120 # 2 × 2560

@info "test 11/26: detect_bad()"
bm, be = detect_bad(e10)
@test bm == Bool[1 1 1 1 1 1 1 1 1 1; 
                 1 1 1 1 1 1 1 1 1 1; 
                 0 0 0 0 1 1 1 1 1 1; 
                 0 1 0 0 1 1 1 1 1 1; 
                 1 1 1 1 1 1 1 1 1 1; 
                 1 1 1 1 1 1 1 1 1 1; 
                 1 1 1 1 1 1 1 1 1 1; 
                 1 1 1 1 1 1 0 0 0 0; 
                 1 1 1 1 1 1 1 1 1 1; 
                 1 1 1 1 1 1 1 1 1 1; 
                 1 1 1 1 1 1 1 1 1 1; 
                 1 1 1 1 1 1 1 1 1 1; 
                 1 0 1 1 1 1 1 1 1 1; 
                 1 1 1 1 1 0 0 0 0 0; 
                 1 1 1 1 1 1 1 0 1 1; 
                 0 0 0 1 0 0 1 1 1 1; 
                 1 1 1 1 0 0 0 1 1 1; 
                 0 1 1 1 1 1 1 1 1 1; 
                 1 1 1 1 1 1 1 1 1 1]
@test be == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

@info "test 12/26: epoch()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
@test epoch_len(e10) == 10*sr(eeg)
e10 = epoch(eeg, ep_n=10)
@test epoch_n(e10) == 10

@info "test 13/26: epoch_time()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
@test e10.epoch_time[1] == 0.0
e10_tmp = epoch_ts(e10, ts=-1.0)
@test e10_tmp.epoch_time[1] == -1.0
epoch_ts!(e10_tmp, ts=-2.0)
@test e10_tmp.epoch_time[1] == -3.0

@info "test 14/26: extract_channel()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10 = epoch(eeg, ep_len=10)
s = extract_channel(e10, ch=1)
@test size(s) == (1, 2560, 121)

@info "test 15/26: extract_epoch()"
e10_tmp = extract_epoch(e10, ep=1)
@test size(e10_tmp.data) == (24, 2560, 1)
@test length(e10_tmp.time_pts) == 2560
@test length(e10_tmp.epoch_time) == 2560

@info "test 16/26: extract_data()"
d = extract_data(e10)
@test size(d) == (19, 2560, 121)

@info "test 17/26: extract_time()"
tpts = extract_time(e10)
@test length(tpts) == 309760

@info "test 18/26: extract_eptime()"
et = extract_eptime(e10)
@test length(et) == 2560

@info "test 19/26: trim()"
s = collect(1:100)
@test trim(s, seg=(1, 10)) == 11:100
m = rand(10, 100)
@test size(trim(m, seg=(1, 10))) == (10, 90)
a = rand(10, 100, 10)
@test size(trim(a, seg=(1, 10))) == (10, 90, 10)
e10_tmp = trim(e10, seg=(0, 11), remove_epochs=true)
@test size(e10_tmp) == (24, 2560, 119)
e10_tmp = trim(e10, seg=(0, 21), remove_epochs=false)
@test size(e10_tmp) == (24, 2560, 118)

@info "test 20/26: delete_marker()"
eeg_mrk = import_edf(joinpath(testfiles_path, "eeg-test-edfplus.edf"))
@test nrow(eeg_mrk.markers) == 45 
delete_marker!(eeg_mrk, n=40)
@test nrow(eeg_mrk.markers) == 44

@info "test 21/26: add_marker()"
add_marker!(eeg_mrk, id="test", start=1988, len=1, desc="test", ch=0)
@test nrow(eeg_mrk.markers) == 45 

@info "test 22/26: edit_marker()"
@test eeg_mrk.markers[45, :id] == "test"
edit_marker!(eeg_mrk, n=45, id="TEST", start=1989, len=1, desc="test", ch=0)
@test eeg_mrk.markers[45, :id] == "TEST"

@info "test 23/26: channel2marker()"
@test nrow(eeg_mrk.markers) == 45
unique(eeg_mrk.data[28, :, :])
idx = getindex.(findall(eeg_mrk.data[28, :, :] .== -11502.913868619822), 1)
eeg_mrk.data[28, idx, :] .= 1.0
idx = getindex.(findall(eeg_mrk.data[28, :, :] .!= 1.0), 1)
eeg_mrk.data[28, idx, :] .= 0.0
channel2marker!(eeg_mrk, ch=28, id="mrk")
@test nrow(eeg_mrk.markers) == 227

@info "test 24/26: epoch()"
eeg_mrk2 = epoch(eeg_mrk, marker="Mark2", offset=0.2, ep_len=1.2)
@test size(eeg_mrk2) == (29, 240, 182)

@info "test 25/26: join()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edfplus.edf"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
e10_tmp = NeuroAnalyzer.join(e10, e10)
@test size(e10_tmp) == (29, 4000, 10)

@info "test 26/26: create()"
@test typeof(create(data_type="eeg")) == NeuroAnalyzer.NEURO
@test typeof(create(data_type="meg")) == NeuroAnalyzer.NEURO
@test typeof(create(data_type="nirs")) == NeuroAnalyzer.NEURO

true