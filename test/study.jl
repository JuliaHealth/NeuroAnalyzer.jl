using NeuroAnalyzer
using Test

eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
epoch!(eeg, ep_len=5)

eeg1 = deepcopy(eeg)
eeg2 = deepcopy(eeg)
eeg2.data .*= 0.75

@info "Test: create_study()"
s = create_study([eeg1, eeg2], [:a, :b])
@test s isa NeuroAnalyzer.STUDY

@info "Test: obj_n()"
@test obj_n(s) == 2

@info "Test: nchannels()"
@test nchannels(s) == 24

@info "Test: nepochs()"
@test nepochs(s) == 241

@info "Test: epoch_len()"
@test epoch_len(s) == 1280

@info "Test: sr()"
@test sr(s) == 256

@info "Test: save_study()"
isfile("test.hdf") && rm("test.hdf")
NeuroAnalyzer.save_study(s, file_name="test.hdf")
@test isfile("test.hdf")

@info "Test: load_study()"
s = NeuroAnalyzer.load_study("test.hdf")
@test s isa NeuroAnalyzer.STUDY
isfile("test.hdf") && rm("test.hdf")

true