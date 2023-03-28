using NeuroAnalyzer
using Test

eeg = import_edf("files/eeg-test-edf.edf")
epoch!(eeg, ep_len=5*256)

eeg1 = deepcopy(eeg)
eeg2 = deepcopy(eeg)
eeg2.data .*= 0.75

@info "test 1/6: create_study()"
s = create_study([eeg1, eeg2], [:a, :b])
@test typeof(s) == NeuroAnalyzer.STUDY

@info "test 2/6: obj_n()"
@test obj_n(s) == 2

@info "test 3/6: channel_n()"
@test channel_n(s) == 24

@info "test 4/6: epoch_n()"
@test epoch_n(s) == 242

@info "test 5/6: epoch_len()"
@test epoch_len(s) == 1280

@info "test 6/6: sr()"
@test sr(s) == 256

true