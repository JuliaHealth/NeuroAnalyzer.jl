using NeuroAnalyzer
using Test

eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
epoch!(eeg, ep_len=5)

eeg1 = deepcopy(eeg)
eeg2 = deepcopy(eeg)
eeg2.data .*= 0.75

@info "Test 1/8: create_study()"
s = create_study([eeg1, eeg2], [:a, :b])
@test s isa NeuroAnalyzer.STUDY

@info "Test 2/8: obj_n()"
@test obj_n(s) == 2

@info "Test 3/8: nchannels()"
@test nchannels(s) == 24

@info "Test 4/8: nepochs()"
@test nepochs(s) == 241

@info "Test 5/8: epoch_len()"
@test epoch_len(s) == 1280

@info "Test 6/8: sr()"
@test sr(s) == 256

@info "Test 7/8: save_study()"
isfile("test.hdf") && rm("test.hdf")
NeuroAnalyzer.save_study(s, file_name="test.hdf")
@test isfile("test.hdf")

@info "Test 8/8: load_study()"
s = NeuroAnalyzer.load_study("test.hdf")
@test s isa NeuroAnalyzer.STUDY
isfile("test.hdf") && rm("test.hdf")

true