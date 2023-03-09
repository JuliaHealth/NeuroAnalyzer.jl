using NeuroAnalyzer
using Test

eeg = import_edf("eeg-test-edf.edf")
eeg = epoch(eeg, epoch_len=5*256)

eeg1 = copy(eeg)
eeg2 = copy(eeg)
eeg2.signals .*= 0.75

s = study_create([eeg1, eeg2], [:a, :b])

@test typeof(s) == NeuroAnalyzer.STUDY
@test study_n(s) == 2
@test study_channel_n(s) == 24
@test study_epoch_n(s) == 242
@test study_epoch_len(s) == 1280
@test study_sr(s) == 256

true