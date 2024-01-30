using NeuroAnalyzer
using Test
using DataFrames

@info "test 1/37: import_bdf()"
eeg = import_bdf(joinpath(testfiles_path, "eeg-test-bdf.bdf"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "BDF"
@test nchannels(eeg) == 17

@info "test 2/37: import_bdf()"
eeg = import_bdf(joinpath(testfiles_path, "eeg-test-bdfplus.bdf"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "BDF+"
@test nchannels(eeg) == 11

@info "test 3/37: import_edf()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "EDF"
@test nchannels(eeg) == 24

@info "test 4/37: import_edf()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edfplus.edf"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "EDF+"
@test nchannels(eeg) == 29

@info "test 5/37: import_digitrack()"
eeg = import_digitrack(joinpath(testfiles_path, "eeg-test-digitrack.txt"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "Digitrack"
@test nchannels(eeg) == 24

@info "test 6/37: import_bv()"
eeg = import_bv(joinpath(testfiles_path, "eeg-test-bv.vhdr"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "BrainVision"
@test nchannels(eeg) == 2

@info "test 7/37: import_csv()"
eeg = import_csv(joinpath(testfiles_path, "eeg-test_txch.csv.gz"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "CSV"
@test nchannels(eeg) == 24

@info "test 8/37: import_csv()"
eeg = import_csv(joinpath(testfiles_path, "eeg-test_chxt.csv.gz"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "CSV"
@test nchannels(eeg) == 24

@info "test 9/37: import_set()"
eeg = import_set(joinpath(testfiles_path, "eeg-test-eeglab.set"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "SET"
@test nchannels(eeg) == 24

@info "test 10/37: import_locs_ced()"
l = import_locs_ced(joinpath(testfiles_path, "locs.ced"))
@test l isa DataFrame

@info "test 11/37: import_locs_locs()"
l = import_locs_locs(joinpath(testfiles_path, "locs.locs"))
@test l isa DataFrame

@info "test 12/37: import_locs_elc()"
l = import_locs_elc(joinpath(testfiles_path, "locs.elc"))
@test l isa DataFrame

@info "test 13/37: import_locs_tsv()"
l = import_locs_tsv(joinpath(testfiles_path, "locs.tsv"))
@test l isa DataFrame

@info "test 14/37: import_locs_sfp()"
l = import_locs_sfp(joinpath(testfiles_path, "locs.sfp"))
@test l isa DataFrame

@info "test 15/37: import_locs_csd()"
l = import_locs_csd(joinpath(testfiles_path, "locs.csd"))
@test l isa DataFrame

@info "test 16/37: import_locs_geo()"
l = import_locs_geo(joinpath(testfiles_path, "locs.geo"))
@test l isa DataFrame

@info "test 17/37: import_locs_mat()"
l = import_locs_mat(joinpath(testfiles_path, "locs.mat"))
@test l isa DataFrame

@info "test 18/37: load_locs()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
eeg = load_locs(eeg, file_name=joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))
@test NeuroAnalyzer._has_locs(eeg) == true

@info "test 19/37: save()"
isfile("test.hdf") && rm("test.hdf")
NeuroAnalyzer.save(eeg, file_name="test.hdf")
@test isfile("test.hdf") == true

@info "test 20/37: load()"
new = NeuroAnalyzer.load("test.hdf")
@test new isa NeuroAnalyzer.NEURO
isfile("test.hdf") && rm("test.hdf")

@info "test 21/37: export_csv()"
isfile("eeg.csv") && rm("eeg.csv")
export_csv(eeg, file_name="eeg.csv", header=false)
@test isfile("eeg.csv") == true
isfile("eeg.csv") && rm("eeg.csv")

@info "test 22/37: export_locs()"
isfile("test_out.ced") && rm("test_out.ced")
export_locs(eeg, file_name="test_out.ced")
@test isfile("test_out.ced") == true
isfile("test_out.ced") && rm("test_out.ced")
isfile("test_out.locs") && rm("test_out.locs")
export_locs(eeg, file_name="test_out.locs")
@test isfile("test_out.locs") == true
isfile("test_out.locs") && rm("test_out.locs")

@info "test 23/37: import_snirf()"
n = import_snirf(joinpath(testfiles_path, "fnirs-test-snirf.snirf"))
@test n isa NeuroAnalyzer.NEURO
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "SNIRF"

@info "test 24/37: import_nirs()"
n = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
@test n isa NeuroAnalyzer.NEURO
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "NIRS"

@info "test 25/37: import_nirx()"
n = import_nirx(joinpath(testfiles_path, "nirx", "NIRS-2020-08-18_001.hdr"))
@test n isa NeuroAnalyzer.NEURO
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "NIRX"

@info "test 26/37: export_markers()"
eeg = import_bdf(joinpath(testfiles_path, "eeg-test-bdfplus.bdf"))
isfile("markers.csv") && rm("markers.csv")
export_markers(eeg, file_name="markers.csv")
@test isfile("markers.csv") == true
isfile("markers.csv") && rm("markers.csv")

@info "test 27/37: import_gdf()"
eeg = import_gdf(joinpath(testfiles_path, "eeg-test-gdf_1.25.gdf"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "GDF"
@test nchannels(eeg) == 16
eeg = import_gdf(joinpath(testfiles_path, "eeg-test-gdf_2.20.gdf"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "GDF"
@test nchannels(eeg) == 65

@info "test 28/37: import_montage()"
ref_list, ref_name = import_montage(joinpath(NeuroAnalyzer.PATH, "montages", "bip_long.mnt"))
@test ref_list == ["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]
@test ref_name == "longitudinal-BIP"

@info "test 29/37: import_npy()"
n = import_npy(joinpath(testfiles_path, "eeg-test-npy.npy"), sampling_rate=256)
@test n isa NeuroAnalyzer.NEURO

@info "test 30/37: import_locs_txt()"
l = import_locs_txt(joinpath(testfiles_path, "locs.txt"))
@test l isa DataFrame

@info "test 31/37: import_duomag()"
mep = import_duomag(joinpath(testfiles_path, "mep-duomag.m"))
@test mep isa NeuroAnalyzer.NEURO
mep = import_duomag(joinpath(testfiles_path, "mep-duomag.ascii"))
@test mep isa NeuroAnalyzer.NEURO

@info "test 32/37: import_thymatron()"
mep = import_thymatron(joinpath(testfiles_path, "thymatron/1.png"))
@test mep isa NeuroAnalyzer.NEURO

@info "test 33/37: import_cnt()"
eeg = import_cnt(joinpath(testfiles_path, "eeg-test-cnt.cnt"))
@test eeg isa NeuroAnalyzer.NEURO

@info "test 34/37: import_locs_dat()"
l = import_locs_dat(joinpath(testfiles_path, "locs.dat"))
@test l isa DataFrame

@info "test 35/37: import_locs_asc()"
l = import_locs_asc(joinpath(testfiles_path, "locs.asc"))
@test l isa DataFrame

@info "test 36/37: import_ncs()"
eeg = import_ncs(joinpath(testfiles_path, "eeg-test-ncs.ncs"))
@test eeg isa NeuroAnalyzer.NEURO

@info "test 37/37: import_xdf()"
eeg = import_xdf(joinpath(testfiles_path, "eeg-test-xdf.xdf"))
@test eeg isa NeuroAnalyzer.NEURO

true