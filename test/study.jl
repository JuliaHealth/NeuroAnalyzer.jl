using NeuroAnalyzer
using Test

ntests = 8

eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
epoch!(eeg, ep_len=5)

eeg1 = deepcopy(eeg)
eeg2 = deepcopy(eeg)
eeg2.data .*= 0.75

@info "Test 1/$ntests: create_study()"
s = create_study([eeg1, eeg2], [:a, :b])
@test s isa NeuroAnalyzer.STUDY

@info "Test 2/$ntests: obj_n()"
@test obj_n(s) == 2

@info "Test 3/$ntests: nchannels()"
@test nchannels(s) == 24

@info "Test 4/$ntests: nepochs()"
@test nepochs(s) == 241

@info "Test 5/$ntests: epoch_len()"
@test epoch_len(s) == 1280

@info "Test 6/$ntests: sr()"
@test sr(s) == 256

@info "Test 7/$ntests: save_study()"
isfile("test.hdf") && rm("test.hdf")
NeuroAnalyzer.save_study(s, file_name="test.hdf")
@test isfile("test.hdf")

@info "Test 8/$ntests: load_study()"
s = NeuroAnalyzer.load_study("test.hdf")
@test s isa NeuroAnalyzer.STUDY
isfile("test.hdf") && rm("test.hdf")

true