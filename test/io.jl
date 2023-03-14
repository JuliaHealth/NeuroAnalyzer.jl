using NeuroAnalyzer
using Test
using DataFrames

@info "test 1/22: import_bdf()"
eeg = import_bdf("files/eeg-test-bdf.bdf")
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "BDF"
@test channel_n(eeg) == 16

@info "test 2/22: import_bdf()"
eeg = import_bdf("files/eeg-test-bdfplus.bdf")
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "BDF+"
@test channel_n(eeg) == 11

@info "test 3/22: import_edf()"
eeg = import_edf("files/eeg-test-edf.edf")
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "EDF"
@test channel_n(eeg) == 24

@info "test 4/22: import_edf()"
eeg = import_edf("files/eeg-test-edfplus.edf")
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "EDF+"
@test channel_n(eeg) == 29

@info "test 5/22: import_digitrack()"
eeg = import_digitrack("files/eeg-test-digitrack.txt")
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "Digitrack"
@test channel_n(eeg) == 24

@info "test 6/22: import_bv()"
eeg = import_bv("files/eeg-test-bv.vhdr")
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "BrainVision"
@test channel_n(eeg) == 2

@info "test 7/22: import_csv()"
eeg = import_csv("files/eeg-test_txch.csv.gz")
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "CSV"
@test channel_n(eeg) == 24

@info "test 8/22: import_csv()"
eeg = import_csv("files/eeg-test_chxt.csv.gz")
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "CSV"
@test channel_n(eeg) == 24

@info "test 9/22: import_set()"
eeg = import_set("files/eeg-test-eeglab.set")
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "SET"
@test channel_n(eeg) == 24

@info "test 10/22: import_locs_ced()"
l = import_locs_ced("files/locs.ced")
@test typeof(l) == DataFrame

@info "test 11/22: import_locs_locs()"
l = import_locs_locs("files/locs.locs")
@test typeof(l) == DataFrame

@info "test 12/22: import_locs_elc()"
l = import_locs_elc("files/locs.elc")
@test typeof(l) == DataFrame

@info "test 13/22: import_locs_tsv()"
l = import_locs_tsv("files/locs.tsv")
@test typeof(l) == DataFrame

@info "test 14/22: import_locs_sfp()"
l = import_locs_sfp("files/locs.sfp")
@test typeof(l) == DataFrame

@info "test 15/22: import_locs_csd()"
l = import_locs_csd("files/locs.csd")
@test typeof(l) == DataFrame

@info "test 16/22: import_locs_geo()"
l = import_locs_geo("files/locs.geo")
@test typeof(l) == DataFrame

@info "test 17/22: import_locs_mat()"
l = import_locs_mat("files/locs.mat")
@test typeof(l) == DataFrame

@info "test 18/22: load_locs()"
eeg = load_locs(eeg, file_name="files/standard-10-20-cap19-elmiko.ced")
@test eeg.header.has_locs == true

@info "test 19/22: save()"
isfile("test.hdf5") && rm("test.hdf5")
save(eeg, file_name="test.hdf5")
@test isfile("test.hdf5") == true

@info "test 20/22: load()"
new = load("test.hdf5")
@test typeof(new) == NeuroAnalyzer.NEURO
isfile("test.hdf5") && rm("test.hdf5")

@info "test 21/22: export_csv()"
isfile("eeg.csv") && rm("eeg.csv")
export_csv(eeg, file_name="eeg.csv", header=false)
@test isfile("eeg.csv") == true
isfile("eeg.csv") && rm("eeg.csv")

@info "test 22/22: export_locs()"
isfile("test_out.ced") && rm("test_out.ced")
export_locs(eeg, file_name="test_out.ced")
@test isfile("test_out.ced") == true
isfile("test_out.ced") && rm("test_out.ced")
isfile("test_out.locs") && rm("test_out.locs")
export_locs(eeg, file_name="test_out.locs")
@test isfile("test_out.locs") == true
isfile("test_out.locs") && rm("test_out.locs")

true