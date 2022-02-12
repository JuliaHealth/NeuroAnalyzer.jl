using NeuroJ
using Test

edf = eeg_import_edf("eeg-test-edf.edf")
@test edf.eeg_header[:version] == 0
@test edf.eeg_header[:channel_n] == 19
@test edf.eeg_header[:channel_locations] == false
@test edf.eeg_header[:channel_locations] == false

edf = eeg_load_electrode_positions(edf, "standard-10-20-cap19.ced")
@test typeof(edf) == NeuroJ.EEG
@test edf.eeg_header[:channel_locations] == true

isfile("test.hdf5") && rm("test.hdf5")

@test eeg_save(edf, "test.hdf5"; overwrite=true) == true
@test isfile("test.hdf5") == true

@test eeg_load("test.hdf5") != false
edf_new = eeg_load("test.hdf5")
@test typeof(edf_new) == NeuroJ.EEG

isfile("test.hdf5") && rm("test.hdf5")

true