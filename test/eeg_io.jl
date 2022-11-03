using NeuroAnalyzer
using Test
using DataFrames

bdf = eeg_import_bdf("eeg-test-bdf.bdf")
@test typeof(bdf) == NeuroAnalyzer.EEG
bdf = eeg_import_bdf("eeg-test-bdfplus.bdf")
@test typeof(bdf) == NeuroAnalyzer.EEG
edf = eeg_import_edf("eeg-test-edf.edf")
@test typeof(edf) == NeuroAnalyzer.EEG
edf = eeg_import_edf("eeg-test-edfplus.edf")
@test typeof(edf) == NeuroAnalyzer.EEG
dt = eeg_import_digitrack("eeg-test-digitrack.txt")
@test typeof(dt) == NeuroAnalyzer.EEG
bv = eeg_import_bv("eeg-test-bv.vhdr")
@test typeof(bv) == NeuroAnalyzer.EEG

edf = eeg_import("eeg-test-edf.edf")
ecg = eeg_extract_channel(edf, channel=24)
eeg_delete_channel!(edf, channel=24)
eog2 = eeg_extract_channel(edf, channel=23)
eeg_delete_channel!(edf, channel=23)
eog1 = eeg_extract_channel(edf, channel=22)
eeg_delete_channel!(edf, channel=22)
a2 = eeg_extract_channel(edf, channel=18)
eeg_delete_channel!(edf, channel=18)
a1 = eeg_extract_channel(edf, channel=17)
eeg_delete_channel!(edf, channel=17)

@test edf.eeg_header[:eeg_filetype] == "EDF"
@test edf.eeg_header[:channel_n] == 19
@test edf.eeg_header[:channel_locations] == false
@test edf.eeg_header[:channel_locations] == false

s = eeg_import_ced("test.ced")
@test typeof(s) == DataFrame
s = eeg_import_locs("test.locs")
@test typeof(s) == DataFrame
s = eeg_import_elc("test.elc")
@test typeof(s) == DataFrame
s = eeg_import_tsv("test.tsv")
@test typeof(s) == DataFrame
s = eeg_import_sfp("test.sfp")
@test typeof(s) == DataFrame
s = eeg_import_csd("test.csd")
@test typeof(s) == DataFrame

edf = eeg_load_electrodes(edf, file_name="standard-10-20-cap19-elmiko.ced")
@test typeof(edf) == NeuroAnalyzer.EEG
@test edf.eeg_header[:channel_locations] == true

isfile("test.hdf5") && rm("test.hdf5")
eeg_save(edf, file_name="test.hdf5")
@test isfile("test.hdf5") == true

edf_new = eeg_load("test.hdf5")
@test typeof(edf_new) == NeuroAnalyzer.EEG
isfile("test.hdf5") && rm("test.hdf5")

isfile("edf.csv") && rm("edf.csv")
eeg_export_csv(edf, file_name="edf.csv", header=false)
@test isfile("edf.csv") == true
isfile("edf.csv") && rm("edf.csv")

isfile("test_out.ced") && rm("test_out.ced")
eeg_save_electrodes(edf, file_name="test_out.ced")
@test isfile("test_out.ced") == true
isfile("test_out.ced") && rm("test_out.ced")

isfile("test_out.locs") && rm("test_out.locs")
eeg_save_electrodes(edf, file_name="test_out.locs")
@test isfile("test_out.locs") == true
isfile("test_out.locs") && rm("test_out.locs")

locs = eeg_import_ced("standard-10-20-cap19-elmiko.ced")
edf2 = eeg_add_electrodes(edf, locs=locs)
@test typeof(edf2) == NeuroAnalyzer.EEG

true