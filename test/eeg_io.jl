using NeuroAnalyzer
using Test
using DataFrames

eeg = eeg_import_bdf("eeg-test-bdf.bdf")
@test typeof(eeg) == NeuroAnalyzer.EEG
eeg = eeg_import_bdf("eeg-test-bdfplus.bdf")
@test typeof(eeg) == NeuroAnalyzer.EEG
eeg = eeg_import_edf("eeg-test-edf.edf")
@test typeof(eeg) == NeuroAnalyzer.EEG
eeg = eeg_import_edf("eeg-test-edfplus.edf")
@test typeof(eeg) == NeuroAnalyzer.EEG
eeg = eeg_import_digitrack("eeg-test-digitrack.txt")
@test typeof(eeg) == NeuroAnalyzer.EEG
eeg = eeg_import_bv("eeg-test-bv.vhdr")
@test typeof(eeg) == NeuroAnalyzer.EEG
eeg = eeg_import_csv("eeg-test_txch.csv.gz")
@test typeof(eeg) == NeuroAnalyzer.EEG
eeg = eeg_import_csv("eeg-test_chxt.csv.gz")
@test typeof(eeg) == NeuroAnalyzer.EEG
eeg = eeg_import_set("eeg-test.set")
@test typeof(eeg) == NeuroAnalyzer.EEG

eeg = eeg_import("eeg-test-edf.edf")
ecg = eeg_extract_channel(eeg, channel=24)
eeg_delete_channel!(eeg, channel=24)
eog2 = eeg_extract_channel(eeg, channel=23)
eeg_delete_channel!(eeg, channel=23)
eog1 = eeg_extract_channel(eeg, channel=22)
eeg_delete_channel!(eeg, channel=22)
a2 = eeg_extract_channel(eeg, channel=21)
eeg_delete_channel!(eeg, channel=21)
a1 = eeg_extract_channel(eeg, channel=20)
eeg_delete_channel!(eeg, channel=20)

@test eeg.eeg_header[:eeg_filetype] == "EDF"
@test eeg.eeg_header[:channel_n] == 19
@test eeg.eeg_header[:channel_locations] == false
@test eeg.eeg_header[:channel_locations] == false

s = locs_import_ced("test.ced")
@test typeof(s) == DataFrame
s = locs_import_locs("test.locs")
@test typeof(s) == DataFrame
s = locs_import_elc("test.elc")
@test typeof(s) == DataFrame
s = locs_import_tsv("test.tsv")
@test typeof(s) == DataFrame
s = locs_import_sfp("test.sfp")
@test typeof(s) == DataFrame
s = locs_import_csd("test.csd")
@test typeof(s) == DataFrame
s = locs_import_geo("test.geo")
@test typeof(s) == DataFrame

eeg = eeg_load_electrodes(eeg, file_name="standard-10-20-cap19-elmiko.ced")
@test typeof(eeg) == NeuroAnalyzer.EEG
@test eeg.eeg_header[:channel_locations] == true

isfile("test.hdf5") && rm("test.hdf5")
eeg_save(eeg, file_name="test.hdf5")
@test isfile("test.hdf5") == true

eeg_new = eeg_load("test.hdf5")
@test typeof(eeg_new) == NeuroAnalyzer.EEG
isfile("test.hdf5") && rm("test.hdf5")

isfile("eeg.csv") && rm("eeg.csv")
eeg_export_csv(eeg, file_name="eeg.csv", header=false)
@test isfile("eeg.csv") == true
isfile("eeg.csv") && rm("eeg.csv")

isfile("test_out.ced") && rm("test_out.ced")
eeg_save_electrodes(eeg, file_name="test_out.ced")
@test isfile("test_out.ced") == true
isfile("test_out.ced") && rm("test_out.ced")

isfile("test_out.locs") && rm("test_out.locs")
eeg_save_electrodes(eeg, file_name="test_out.locs")
@test isfile("test_out.locs") == true
isfile("test_out.locs") && rm("test_out.locs")

locs = locs_import_ced("standard-10-20-cap19-elmiko.ced")
eeg2 = eeg_add_electrodes(eeg, locs=locs)
@test typeof(eeg2) == NeuroAnalyzer.EEG

true