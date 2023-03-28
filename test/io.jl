using NeuroAnalyzer
using Test
using DataFrames

@info "test 1/25: import_bdf()"
eeg = import_bdf(joinpath(testfiles_path, "eeg-test-bdf.bdf"))
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "BDF"
@test channel_n(eeg) == 16

@info "test 2/25: import_bdf()"
eeg = import_bdf(joinpath(testfiles_path, "eeg-test-bdfplus.bdf"))
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "BDF+"
@test channel_n(eeg) == 11

@info "test 3/25: import_edf()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "EDF"
@test channel_n(eeg) == 24

@info "test 4/25: import_edf()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edfplus.edf"))
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "EDF+"
@test channel_n(eeg) == 29

@info "test 5/25: import_digitrack()"
eeg = import_digitrack(joinpath(testfiles_path, "eeg-test-digitrack.txt"))
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "Digitrack"
@test channel_n(eeg) == 24

@info "test 6/25: import_bv()"
eeg = import_bv(joinpath(testfiles_path, "eeg-test-bv.vhdr"))
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "BrainVision"
@test channel_n(eeg) == 2

@info "test 7/25: import_csv()"
eeg = import_csv(joinpath(testfiles_path, "eeg-test_txch.csv.gz"))
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "CSV"
@test channel_n(eeg) == 24

@info "test 8/25: import_csv()"
eeg = import_csv(joinpath(testfiles_path, "eeg-test_chxt.csv.gz"))
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "CSV"
@test channel_n(eeg) == 24

@info "test 9/25: import_set()"
eeg = import_set(joinpath(testfiles_path, "eeg-test-eeglab.set"))
@test typeof(eeg) == NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "SET"
@test channel_n(eeg) == 24

@info "test 10/25: import_locs_ced()"
l = import_locs_ced(joinpath(testfiles_path, "locs.ced"))
@test typeof(l) == DataFrame

@info "test 11/25: import_locs_locs()"
l = import_locs_locs(joinpath(testfiles_path, "locs.locs"))
@test typeof(l) == DataFrame

@info "test 12/25: import_locs_elc()"
l = import_locs_elc(joinpath(testfiles_path, "locs.elc"))
@test typeof(l) == DataFrame

@info "test 13/25: import_locs_tsv()"
l = import_locs_tsv(joinpath(testfiles_path, "locs.tsv"))
@test typeof(l) == DataFrame

@info "test 14/25: import_locs_sfp()"
l = import_locs_sfp(joinpath(testfiles_path, "locs.sfp"))
@test typeof(l) == DataFrame

@info "test 15/25: import_locs_csd()"
l = import_locs_csd(joinpath(testfiles_path, "locs.csd"))
@test typeof(l) == DataFrame

@info "test 16/25: import_locs_geo()"
l = import_locs_geo(joinpath(testfiles_path, "locs.geo"))
@test typeof(l) == DataFrame

@info "test 17/25: import_locs_mat()"
l = import_locs_mat(joinpath(testfiles_path, "locs.mat"))
@test typeof(l) == DataFrame

@info "test 18/25: load_locs()"
eeg = load_locs(eeg, file_name=joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))
@test NeuroAnalyzer._has_locs(eeg) == true

@info "test 19/25: save()"
isfile("test.hdf5") && rm("test.hdf5")
save(eeg, file_name="test.hdf5")
@test isfile("test.hdf5") == true

@info "test 20/25: load()"
new = load("test.hdf5")
@test typeof(new) == NeuroAnalyzer.NEURO
isfile("test.hdf5") && rm("test.hdf5")

@info "test 21/25: export_csv()"
isfile("eeg.csv") && rm("eeg.csv")
export_csv(eeg, file_name="eeg.csv", header=false)
@test isfile("eeg.csv") == true
isfile("eeg.csv") && rm("eeg.csv")

@info "test 22/25: export_locs()"
isfile("test_out.ced") && rm("test_out.ced")
export_locs(eeg, file_name="test_out.ced")
@test isfile("test_out.ced") == true
isfile("test_out.ced") && rm("test_out.ced")
isfile("test_out.locs") && rm("test_out.locs")
export_locs(eeg, file_name="test_out.locs")
@test isfile("test_out.locs") == true
isfile("test_out.locs") && rm("test_out.locs")

@info "test 23/25: import_snirf()"
n = import_snirf(joinpath(testfiles_path, "fnirs-test-snirf.snirf"))
@test typeof(n) == NeuroAnalyzer.NEURO
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "SNIRF"

@info "test 24/25: import_nirs()"
n = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
@test typeof(n) == NeuroAnalyzer.NEURO
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "NIRS"

@info "test 25/25: import_nirx()"
n = import_nirx(joinpath(testfiles_path, "nirx", "NIRS-2020-08-18_001.hdr"))
@test typeof(n) == NeuroAnalyzer.NEURO
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "NIRX"

true