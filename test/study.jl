using NeuroAnalyzer
using Test

eeg = import_edf("eeg-test-edf.edf")
eeg = epoch(eeg, ep_len=5*256)

eeg1 = deepcopy(eeg)
eeg2 = deepcopy(eeg)
eeg2.data .*= 0.75

s = study([eeg1, eeg2], [:a, :b])

@test typeof(s) == NeuroAnalyzer.STUDY
@test obj_n(s) == 2
@test channel_n(s) == 24
@test epoch_n(s) == 242
@test epoch_len(s) == 1280
@test sr(s) == 256

true