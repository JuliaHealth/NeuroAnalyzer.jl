using NeuroAnalyzer
using Test

eeg = eeg_import_edf("eeg-test-edf.edf")
eeg = eeg_epoch(eeg, epoch_len=5*256)

eeg1 = eeg_copy(eeg)
eeg2 = eeg_copy(eeg)
eeg2.eeg_signals .*= 0.75

s = eeg_study_create([eeg1, eeg2], [:a, :b])

@test typeof(s) == NeuroAnalyzer.STUDY
@test eeg_study_n(s) == 2
@test eeg_study_epoch_n(s) == 242
@test eeg_study_epoch_len(s) == 1280
@test eeg_study_channel_n(s) == 24

true