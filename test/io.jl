using NeuroAnalyzer
using Test
using DataFrames

ntests = 38

@info "Test 1/$ntests: import_bdf()"
eeg = import_bdf(joinpath(testfiles_path, "eeg-test-bdf.bdf"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "BDF"
@test nchannels(eeg) == 17

@info "Test 2/$ntests: import_bdf()"
eeg = import_bdf(joinpath(testfiles_path, "eeg-test-bdfplus.bdf"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "BDF+"
@test nchannels(eeg) == 11

@info "Test 3/$ntests: import_edf()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "EDF"
@test nchannels(eeg) == 24

@info "Test 4/$ntests: import_edf()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edfplus.edf"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "EDF+"
@test nchannels(eeg) == 29

@info "Test 5/$ntests: import_digitrack()"
eeg = import_digitrack(joinpath(testfiles_path, "eeg-test-digitrack.txt"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "Digitrack"
@test nchannels(eeg) == 24

@info "Test 6/$ntests: import_bv()"
eeg = import_bv(joinpath(testfiles_path, "eeg-test-bv.vhdr"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "BrainVision"
@test nchannels(eeg) == 2

@info "Test 7/$ntests: import_csv()"
eeg = import_csv(joinpath(testfiles_path, "eeg-test_txch.csv.gz"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "CSV"
@test nchannels(eeg) == 24

@info "Test 8/$ntests: import_csv()"
eeg = import_csv(joinpath(testfiles_path, "eeg-test_chxt.csv.gz"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "CSV"
@test nchannels(eeg) == 24

@info "Test 9/$ntests: import_set()"
eeg = import_set(joinpath(testfiles_path, "eeg-test-eeglab.set"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "SET"
@test nchannels(eeg) == 24

@info "Test 10/$ntests: import_locs_ced()"
l = import_locs_ced(joinpath(testfiles_path, "locs.ced"));
@test l isa DataFrame

@info "Test 11/$ntests: import_locs_locs()"
l = import_locs_locs(joinpath(testfiles_path, "locs.locs"));
@test l isa DataFrame

@info "Test 12/$ntests: import_locs_elc()"
l = import_locs_elc(joinpath(testfiles_path, "locs.elc"));
@test l isa DataFrame

@info "Test 13/$ntests: import_locs_tsv()"
l = import_locs_tsv(joinpath(testfiles_path, "locs.tsv"));
@test l isa DataFrame

@info "Test 14/$ntests: import_locs_sfp()"
l = import_locs_sfp(joinpath(testfiles_path, "locs.sfp"));
@test l isa DataFrame

@info "Test 15/$ntests: import_locs_csd()"
l = import_locs_csd(joinpath(testfiles_path, "locs.csd"));
@test l isa DataFrame

@info "Test 16/$ntests: import_locs_geo()"
l = import_locs_geo(joinpath(testfiles_path, "locs.geo"));
@test l isa DataFrame

@info "Test 17/$ntests: import_locs_mat()"
l = import_locs_mat(joinpath(testfiles_path, "locs.mat"));
@test l isa DataFrame

@info "Test 18/$ntests: load_locs()"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"));
eeg = load_locs(eeg, file_name=joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"));
@test NeuroAnalyzer._has_locs(eeg)

@info "Test 19/$ntests: save()"
isfile("test.hdf") && rm("test.hdf")
NeuroAnalyzer.save(eeg, file_name="test.hdf")
@test isfile("test.hdf")

@info "Test 20/$ntests: load()"
new = NeuroAnalyzer.load("test.hdf")
@test new isa NeuroAnalyzer.NEURO
isfile("test.hdf") && rm("test.hdf")

@info "Test 21/$ntests: export_csv()"
isfile("eeg.csv") && rm("eeg.csv")
export_csv(eeg, file_name="eeg.csv", header=false)
@test isfile("eeg.csv")
isfile("eeg.csv") && rm("eeg.csv")

@info "Test 22/$ntests: export_locs()"
isfile("test_out.ced") && rm("test_out.ced")
export_locs(eeg, file_name="test_out.ced")
@test isfile("test_out.ced")
isfile("test_out.ced") && rm("test_out.ced")
isfile("test_out.locs") && rm("test_out.locs")
export_locs(eeg, file_name="test_out.locs")
@test isfile("test_out.locs")
isfile("test_out.locs") && rm("test_out.locs")

@info "Test 23/$ntests: import_snirf()"
n = import_snirf(joinpath(testfiles_path, "fnirs-test-snirf.snirf"));
@test n isa NeuroAnalyzer.NEURO
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "SNIRF"

@info "Test 24/$ntests: import_nirs()"
n = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"));
@test n isa NeuroAnalyzer.NEURO
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "NIRS"

@info "Test 25/$ntests: import_nirx()"
n = import_nirx(joinpath(testfiles_path, "nirx", "NIRS-2020-08-18_001.hdr"));
@test n isa NeuroAnalyzer.NEURO
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "NIRX"

@info "Test 26/$ntests: export_markers()"
eeg = import_bdf(joinpath(testfiles_path, "eeg-test-bdfplus.bdf"));
isfile("markers.csv") && rm("markers.csv")
export_markers(eeg, file_name="markers.csv")
@test isfile("markers.csv")
isfile("markers.csv") && rm("markers.csv")

@info "Test 27/$ntests: import_gdf()"
eeg = import_gdf(joinpath(testfiles_path, "eeg-test-gdf_1.25.gdf"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "GDF"
@test nchannels(eeg) == 16
eeg = import_gdf(joinpath(testfiles_path, "eeg-test-gdf_2.20.gdf"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:file_type] == "GDF"
@test nchannels(eeg) == 65

@info "Test 28/$ntests: import_montage()"
ref_list, ref_name = import_montage(joinpath(NeuroAnalyzer.PATH, "montages", "bip_long.mnt"));
@test ref_list == ["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]
@test ref_name == "longitudinal-BIP"

@info "Test 29/$ntests: import_npy()"
n = import_npy(joinpath(testfiles_path, "eeg-test-npy.npy"), sampling_rate=256)
@test n isa NeuroAnalyzer.NEURO

@info "Test 30/$ntests: import_locs_txt()"
l = import_locs_txt(joinpath(testfiles_path, "locs.txt"));
@test l isa DataFrame

@info "Test 31/$ntests: import_duomag()"
mep = import_duomag(joinpath(testfiles_path, "mep-duomag.m"));
@test mep isa NeuroAnalyzer.NEURO
@test mep.header.recording[:data_type] == "mep"
@test mep.header.recording[:file_type] == "DuoMAG"
mep = import_duomag(joinpath(testfiles_path, "mep-duomag.ascii"));
@test mep isa NeuroAnalyzer.NEURO
@test mep.header.recording[:data_type] == "mep"
@test mep.header.recording[:file_type] == "DuoMAG"

@info "Test 32/$ntests: import_thymatron()"
eeg = import_thymatron(joinpath(testfiles_path, "thymatron/1.png"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "Thymatron"

@info "Test 33/$ntests: import_cnt()"
eeg = import_cnt(joinpath(testfiles_path, "eeg-test-cnt.cnt"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "CNT"

@info "Test 34/$ntests: import_locs_dat()"
l = import_locs_dat(joinpath(testfiles_path, "locs.dat"));
@test l isa DataFrame

@info "Test 35/$ntests: import_locs_asc()"
l = import_locs_asc(joinpath(testfiles_path, "locs.asc"));
@test l isa DataFrame

@info "Test 36/$ntests: import_ncs()"
eeg = import_ncs(joinpath(testfiles_path, "eeg-test-ncs.ncs"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "ieeg"
@test eeg.header.recording[:file_type] == "CSC"

@info "Test 37/$ntests: import_xdf()"
eeg = import_xdf(joinpath(testfiles_path, "eeg-test-xdf.xdf"));
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "XDF"

@info "Test 38/$ntests: import_fiff()"
meg = import_fiff(joinpath(testfiles_path, "meg-test-fiff.fif"));
@test meg isa NeuroAnalyzer.NEURO
@test meg.header.recording[:data_type] == "meg"
@test meg.header.recording[:file_type] == "FIFF"

true