using NeuroAnalyzer
using Test

eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
epoch!(eeg, ep_len=5)

eeg1 = deepcopy(eeg)
eeg2 = deepcopy(eeg)
eeg2.data .*= 0.75

@info "test 1/6: create_study()"
s = create_study([eeg1, eeg2], [:a, :b])
@test s isa NeuroAnalyzer.STUDY

@info "test 2/6: obj_n()"
@test obj_n(s) == 2

@info "test 3/6: nchannels()"
@test nchannels(s) == 24

@info "test 4/6: nepochs()"
@test nepochs(s) == 241

@info "test 5/6: epoch_len()"
@test epoch_len(s) == 1280

@info "test 6/6: sr()"
@test sr(s) == 256

true