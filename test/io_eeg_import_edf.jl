using NeuroJ
using Test

edf = eeg_import_edf("eeg-test-edf.edf")
@test edf.eeg_header[:version] == 0
@test edf.eeg_header[:channels_no] == 19
@test edf.eeg_header[:channel_locations] == false

true