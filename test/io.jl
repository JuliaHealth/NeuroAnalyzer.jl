using NeuroAnalyzer
using Test
using DataFrames

# reference values: MNE-Python 1.13.2 (mne.io.read_raw_*), converted to NeuroAnalyzer units
# (μV for EEG, fT for magnetometers, fT/cm for gradiometers); samples are compared at
# positions 1, 2, 101, n ÷ 2 + 1 and n
tmp = mktempdir()
pos(n) = [1, 2, 101, n ÷ 2 + 1, n]

@info "Test: detect_format()"
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf")) === :edf
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test-bdf.bdf")) === :bdf
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test-gdf_2.20.gdf")) === :gdf
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test-bv.vhdr")) === :bv
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test-eeglab.set")) === :set
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test-npy.npy")) === :npy
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test-xdf.xdf")) === :xdf
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test-ncs.ncs")) === :ncs
@test detect_format(joinpath(testfiles_path, "meg", "meg-test-fiff.fif")) === :fiff
@test detect_format(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf")) === :snirf
@test detect_format(joinpath(testfiles_path, "nirs", "fnirs-test-nirs.nirs")) === :nirs
@test detect_format(joinpath(testfiles_path, "nirs", "nirx", "NIRS-2020-08-18_001.hdr")) === :nirx
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test_txch.csv.gz")) === :csv_gz
@test detect_format(joinpath(testfiles_path, "eeg", "eeg-test-fieldtrip.mat")) === :ft
@test detect_format(joinpath(testfiles_path, "mep", "mep-duomag.ascii")) === :duomag
# content wins over the extension (with a warning); HDF5 flavor detected without an extension
cp(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"), joinpath(tmp, "edf.bdf"))
@test (@test_logs (:warn,) detect_format(joinpath(tmp, "edf.bdf"))) === :edf
cp(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf"), joinpath(tmp, "snirf_noext"))
@test detect_format(joinpath(tmp, "snirf_noext")) === :snirf
touch(joinpath(tmp, "empty.edf"))
@test_throws ArgumentError detect_format(joinpath(tmp, "empty.edf"))
write(joinpath(tmp, "unknown.xyz"), UInt8[0x00, 0x01, 0x02])
@test_throws ArgumentError detect_format(joinpath(tmp, "unknown.xyz"))
# a BrainVision data file instead of its header: the error suggests the header
err = try
    detect_format(joinpath(testfiles_path, "eeg", "eeg-test-bv.eeg"))
catch e
    e
end
@test err isa ArgumentError && occursin("eeg-test-bv.vhdr", err.msg)
@test NeuroAnalyzer._sniff_format(joinpath(testfiles_path, "eeg", "eeg-test-fieldtrip.mat")) === :mat73
@test NeuroAnalyzer._sniff_format(joinpath(testfiles_path, "nirs", "fnirs-test-nirs.nirs")) === :mat5
@test NeuroAnalyzer._sniff_format(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf")) === :hdf5
@test NeuroAnalyzer._sniff_format(joinpath(testfiles_path, "locs", "locs.ced")) === :text
@test NeuroAnalyzer._sniff_format(joinpath(tmp, "empty.edf")) === :empty
@test NeuroAnalyzer._ext_format("a.csv.gz") === :csv_gz
@test NeuroAnalyzer._ext_format("a.tar.gz") === nothing
@test NeuroAnalyzer._ext_format("a.FIF") === :fiff
@test NeuroAnalyzer._ext_format("a.xyz") === nothing
@test NeuroAnalyzer._hasprefix(UInt8[0x47, 0x44, 0x46, 0x20, 0x01], "GDF ")
@test !NeuroAnalyzer._hasprefix(UInt8[0x47], "GDF ")
@test NeuroAnalyzer._hdf5_flavor(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf")) === :snirf
@test NeuroAnalyzer._hdf5_flavor(joinpath(testfiles_path, "eeg", "eeg.hdf5")) === nothing
@test NeuroAnalyzer._hdf5_flavor(joinpath(testfiles_path, "locs", "locs.ced")) === nothing
@test NeuroAnalyzer._sibling_hint(joinpath(testfiles_path, "eeg", "eeg-test-bv.vmrk")) ==
      " Did you mean $(joinpath(testfiles_path, "eeg", "eeg-test-bv.vhdr"))?"
@test NeuroAnalyzer._sibling_hint(joinpath(testfiles_path, "locs", "locs.ced")) == ""

@info "Test: import_edf()"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
@test eeg isa NeuroAnalyzer.NEURO
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "EDF"
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (256, 24, 308480)
@test labels(eeg)[1:3] == ["Fp1", "Fp2", "F3"]
@test eeg.header.recording[:unit][1] == "μV"
@test eeg.data[1, pos(308480), 1] ≈ [741.6464417486837, 923.6939362172883, 1149.145857938506, 558.3434473182268, 0.00467856870363903] rtol = 1e-9
@test eeg.data[24, pos(308480), 1] ≈ [-1115.0586446936752, -8.066392614633507, 433.33152254520473, 1103.230430762188, 0.00467856870363903] rtol = 1e-9
@test eeg.time_pts[1:2] == [0.0, 1 / 256]
@test_throws ArgumentError import_edf(joinpath(testfiles_path, "no_such_file.edf"))
@test_throws ArgumentError import_edf(joinpath(testfiles_path, "eeg", "eeg-test-bdf.bdf"))

@info "Test: import_edf() - EDF+"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edfplus.edf"))
@test eeg.header.recording[:file_type] == "EDF+"
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (200, 29, 400200)
@test eeg.data[1, pos(400200), 1] ≈ [49.99982543678964, 49.99982543678964, -0.00017187762250614466, 16.40607724116898, 49.99982543678964] rtol = 1e-9
# "E" is an EOG-type signal channel (in μV), not a marker channel
@test eeg.header.recording[:channel_type][20] == "eog"
@test eeg.data[20, pos(400200), 1] ≈ [49.99982543678964, 49.99982543678964, -0.00017187762250614466, -11.718921248187852, 49.99982543678964] rtol = 1e-9
# annotations (MNE: 45, the same onsets; trailing spaces are stripped)
@test nrow(eeg.markers) == 45
@test eeg.markers[1:3, :start] == [0.0, 1.14, 4.22]
@test eeg.markers[1:3, :value] == ["REC START BioCal CAL", "A1+A2 OFF", "IMP CHECK ON"]

@info "Test: import_alice4()"
# an Alice 4 export: "Alice 4" in the recording field and an unknown (-1) number of data records
raw = read(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
raw[89:168] = codeunits(rpad("Startdate X X X Alice 4", 80))
raw[237:244] = codeunits(rpad("-1", 8))
write(joinpath(tmp, "alice4.edf"), raw)
a4 = import_alice4(joinpath(tmp, "alice4.edf"))
@test a4.data == import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf")).data
@test_throws ArgumentError import_alice4(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
@test_throws ArgumentError import_alice4(joinpath(tmp, "no_such_file.edf"))

@info "Test: import_edf_annotations()"
# annotation-only EDF+: one "EDF Annotations" signal, data record duration 0
function edf_annotations_file(f, recs)
    spr = 30
    hdr = rpad("0", 8) * rpad("X X X X", 80) * rpad("Startdate X X X X", 80) * "01.01.20" * "00.00.00" *
          rpad(string(512), 8) * rpad("EDF+C", 44) * rpad(string(length(recs)), 8) * rpad("0", 8) * rpad("1", 4) *
          rpad("EDF Annotations", 16) * rpad("", 80) * rpad("", 8) * rpad("-1", 8) * rpad("1", 8) *
          rpad("-32768", 8) * rpad("32767", 8) * rpad("", 80) * rpad(string(spr), 8) * rpad("", 32)
    open(f, "w") do io
        write(io, hdr)
        for r in recs
            write(io, r * repeat("\0", 2 * spr - ncodeunits(r)))
        end
    end
end
edf_annotations_file(joinpath(tmp, "ann.edf"), ["+0\x14\x14\0", "+1.5\x150.5\x14stim\x14\0", "+3\x14\x14resp\x14\0"])
ann = import_edf_annotations(joinpath(tmp, "ann.edf"))
@test ann[!, :start] == [1.5, 3.0]
@test ann[!, :length] == [0.5, 0.0]
@test ann[!, :value] == ["stim", "resp"]
# a regular EDF with signals is rejected
@test_throws ArgumentError import_edf_annotations(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
@test_throws ArgumentError import_edf_annotations(joinpath(testfiles_path, "locs", "locs.ced"))

@info "Test: import_bdf()"
eeg = import_bdf(joinpath(testfiles_path, "eeg", "eeg-test-bdf.bdf"))
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "BDF"
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (2048, 17, 122880)
@test eeg.data[1, pos(122880), 1] ≈ [-636.9844129672297, -449.9844018211604, -485.98440396692774, -489.9844042053463, -545.9844075432067] rtol = 1e-9
@test eeg.data[16, pos(122880), 1] ≈ [-261.9843906154865, -260.9843905558819, -245.98438966181217, -240.9843893637889, -201.9843870392076] rtol = 1e-9
# Status: trigger bits 0-15 only
@test eeg.header.recording[:channel_type][17] == "mrk"
@test all(0 .<= eeg.data[17, :, 1] .< 2^16)

@info "Test: import_bdf() - BDF+"
eeg = import_bdf(joinpath(testfiles_path, "eeg", "eeg-test-bdfplus.bdf"))
@test eeg.header.recording[:file_type] == "BDF+"
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (200, 11, 120000)
@test eeg.data[1, pos(120000), 1] ≈ [99.99997019767581, 99.99997019767581, 99.99997019767581, 99.99997019767581, -99.9998509883792] rtol = 1e-9
@test eeg.data[4, pos(120000), 1] ≈ [5.960464829968259e-05, 5.960464829968259e-05, -0.0005364418349813624, 5.960464829968259e-05, 5.960464829968259e-05] rtol = 1e-6
@test eeg.markers[!, :start] == [0.0, 600.0]
@test eeg.markers[!, :value] == ["RECORD START", "REC STOP"]

@info "Test: import_gdf()"
eeg = import_gdf(joinpath(testfiles_path, "eeg", "eeg-test-gdf_1.25.gdf"))
@test eeg.header.recording[:file_type] == "GDF"
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (128, 16, 34888)
@test eeg.data[1, pos(34888), 1] ≈ [-0.002365354838579789, -1.4558388368327606, -196819.04730317753, 10.052126447075345, 2.0511066003997165] rtol = 1e-9
@test eeg.data[16, pos(34888), 1] ≈ [-0.01143420344710269, -2.3150784839536866, -245748.8208853573, -2.714158265248663, -2.646438793674498] rtol = 1e-9
# 963 events; MNE: onsets 40.054688, 44.507812, 44.507812 s, codes 32769 (0x8001), 33025 (0x8101), 33286 (0x8206)
@test nrow(eeg.markers) == 963
@test eeg.markers[1:3, :start] ≈ [40.0546875, 44.5078125, 44.5078125]
@test eeg.markers[1:3, :value] == ["condition 1 (end)", "artifact:EOG (blinks, fast, large amplitude) (end)", "unknown code (0x8206)"]
eeg = import_gdf(joinpath(testfiles_path, "eeg", "eeg-test-gdf_2.20.gdf"))
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (2048, 65, 22477)
@test eeg.data[2, pos(22477), 1] ≈ [-7210.375, -7190.124999999999, -7135.1875, -6964.8125, -7170.968749999999] rtol = 1e-9
@test eeg.data[65, pos(22477), 1] ≈ [8374.6875, 8395.625, 8438.593749999998, 8613.593749999998, 8330.9375] rtol = 1e-9

@info "Test: import_bv()"
eeg = import_bv(joinpath(testfiles_path, "eeg", "eeg-test-bv.vhdr"))
@test eeg.header.recording[:file_type] == "BrainVision"
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (1000, 2, 10000)
@test labels(eeg) == ["chan1", "chan2"]
@test eeg.data[1, pos(10000), 1] ≈ [2.129114723205566, -6.382275772094726, 9.197168731689453, -1.579151916503906, -3.148046112060547] rtol = 1e-9
@test eeg.data[2, pos(10000), 1] ≈ [-23.854678344726562, 16.196929931640625, -4.036828231811524, 6.430135345458984, -8.600655364990235] rtol = 1e-9
@test nrow(eeg.markers) == 10
@test eeg.markers[1:3, :start] == [0.499, 1.499, 2.499]
@test eeg.markers[1:3, :length] == [0.001, 0.001, 0.001]
@test eeg.markers[1, :value] == "S  1"

@info "Test: import_set()"
eeg = import_set(joinpath(testfiles_path, "eeg", "eeg-test-eeglab.set"))
@test eeg.header.recording[:file_type] == "SET"
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (256, 24, 309760)
@test eeg.data[1, pos(309760), 1] ≈ [222.40752895399393, 182.76960158693817, 105.82538963912403, 78.74246189059271, 0.00467856870363903] rtol = 1e-6
# EEGLAB stores all channels in μV: the ECG channel is in μV (not mV)
@test eeg.header.recording[:unit][24] == "μV"
@test eeg.data[24, pos(309760), 1] ≈ [68.87781933318062, 192.45488700694276, -103.48438971541934, 174.33981612878603, 0.00467856870363903] rtol = 1e-6
eeg = import_set(joinpath(testfiles_path, "eeg", "eeg-test-10events.set"))
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (1000, 2, 10000)
@test eeg.data[1, pos(10000), 1] ≈ [2.1291146278381348, -6.382275581359863, 9.197168350219727, -1.5791518688201904, -3.1480460166931152] rtol = 1e-6
@test nrow(eeg.markers) == 10
@test eeg.markers[1:3, :start] == [0.499, 1.499, 2.499]
@test eeg.markers[1, :value] == "S  1"
eeg = import_set(joinpath(testfiles_path, "eeg", "eeg-test-7epochs.set"))
@test nepochs(eeg) == 7

@info "Test: import_cnt()"
eeg = import_cnt(joinpath(testfiles_path, "eeg", "eeg-test-cnt.cnt"))
@test eeg.header.recording[:file_type] == "CNT"
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (500, 34, 270800)
@test labels(eeg)[1:3] == ["O2", "O1", "Pz"]
@test eeg.data[1, pos(270800), 1] ≈ [-42.86018885235535, -43.15317089523887, -40.299264841887634, -46.16379861829046, -11.251637300738366] rtol = 1e-6
# all channels are scaled to μV: EKG in μV (not mV)
@test eeg.header.recording[:unit][32] == "μV"
@test eeg.data[32, pos(270800), 1] ≈ [-90.15840326523175, -93.32797671195294, -140.2630772467819, -13.405559434831957, -497.9152718267869] rtol = 1e-6
# 359 events; value = trigger code (StimType); MNE places the events 1 sample (2 ms) earlier
@test nrow(eeg.markers) == 359
@test eeg.markers[1:3, :value] == ["4", "4", "1"]
@test eeg.markers[1:3, :start] ≈ [12.996, 14.594, 16.194] .+ 1 / 500

@info "Test: import_fiff()"
meg = import_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fif"))
@test meg.header.recording[:data_type] == "meg"
@test meg.header.recording[:file_type] == "FIFF"
@test (sr(meg), nchannels(meg), signal_len(meg)) == (1000, 327, 32000)
@test labels(meg)[1] == "IASX+"
@test meg.header.recording[:channel_type][[12, 13]] == ["mag", "grad"]
@test meg.header.recording[:unit][[12, 13]] == ["fT", "fT/cm"]
@test meg.data[12, pos(32000), 1] ≈ [37475.75757133061, 47371.71169613672, 37014.06476161138, -2739.7155741580623, -558.0902095507164] rtol = 1e-9
@test meg.data[13, pos(32000), 1] ≈ [-73.67706281630497, 193.15338089679952, 93.58978249638739, -304.6646111052611, 205.10101270484896] rtol = 1e-9
# DACQ events = MNE find_events(STI101): 10 events, 17.944 … 26.985 s, codes 1-10
@test nrow(meg.markers) == 10
@test meg.markers[[1, 2, 10], :start] ≈ [17.944, 18.945, 26.985]
@test meg.markers[!, :value] == string.(1:10)
@test meg.header.recording[:first_time] == 9.0
meg = import_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fiff"))
@test nchannels(meg) == 327
# non-integer sampling rate (600.615 Hz): stored as 601 Hz without resampling (MEG data must not be
# resampled before Maxwell filtering); the markers stay at their samples
meg = import_fiff(joinpath(testfiles_path, "meg", "test_raw.fif"))
@test (sr(meg), nchannels(meg), signal_len(meg)) == (601, 376, 14400)
@test labels(meg)[1] == "MEG 113"
@test meg.data[3, pos(14400), 1] ≈ [101.07422236953386, 63.17138898095866, -12.634277796191732, 12.634277796191732, -214.78272253525944] rtol = 1e-9
@test meg.data[1, pos(14400), 1] ≈ [96.43554807784935, 0.0, -28.930664423354806, 57.86132884670961, -38.574219231139736] rtol = 1e-9 atol = 1e-9
@test meg.data[316, pos(14400), 1] ≈ [11.398925963135298, 9.850158848578872, -10.283813640654671, 0.6195068458225705, 29.05487106907856] rtol = 1e-9
# MNE find_events(STI 014): 32 events at samples 2177, 2545, 2971, 3419, 3852 (from the first sample), codes 2, 3, 1, 4, 2
@test nrow(meg.markers) == 32
@test round.(Int64, meg.markers[1:5, :start] .* sr(meg)) == [2177, 2545, 2971, 3419, 3852]
@test meg.markers[1:5, :value] == ["2", "3", "1", "4", "2"]
@test_throws ArgumentError import_fiff(joinpath(testfiles_path, "locs", "locs.ced"))

@info "Test: load_fiff()"
fiff, fobj, fb = load_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fif"))
@test fiff[:meas_info][:sfreq] == 1000.0
@test fiff[:meas_info][:nchan] == 327
@test size(fiff[:raw_data][:data]) == (327, 32000)
@test size(fb, 2) == 6
fiff, _, _ = load_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fif"); read_data = false)
@test isnothing(fiff[:raw_data])

@info "Test: FIFF import helpers"
fiff, fobj, fb = load_fiff(joinpath(testfiles_path, "meg", "test_raw.fif"))
mi = fiff[:meas_info]
@test NeuroAnalyzer._fiff_get(fobj, "nchan"; block = "meas_info") == 376
@test NeuroAnalyzer._fiff_get(fobj, "sfreq") ≈ 600.614990234375
@test length(NeuroAnalyzer._fiff_getall(fobj, "ch_info")) == 376
@test isnothing(NeuroAnalyzer._fiff_get(fobj, "no_such_tag_name_xyz"))
@test length(NeuroAnalyzer._fiff_instances(fobj, "meas_info")) == 1
@test NeuroAnalyzer._fiff_name_list("A:B{COLON}C") == ["A", "B:C"]
@test NeuroAnalyzer._fiff_name_list(nothing) == String[]
@test NeuroAnalyzer._fiff_vec(nothing) == Any[]
@test NeuroAnalyzer._fiff_vec(5) == [5]
@test NeuroAnalyzer._fiff_vec([1, 2]) == [1, 2]
@test isnothing(NeuroAnalyzer._fiff_nan2nothing(NaN))
@test NeuroAnalyzer._fiff_nan2nothing(1.5) == 1.5
@test NeuroAnalyzer._fiff_labels(mi) == labels(import_fiff(joinpath(testfiles_path, "meg", "test_raw.fif")))
# sphere fit: points on a known sphere
θ = range(0.1, 3.0; length = 20)
φ = range(0, 2π; length = 20)
sp_pts = reduce(vcat, [[0.1 + 0.09 * sin(a) * cos(b) -0.02 + 0.09 * sin(a) * sin(b) 0.04 + 0.09 * cos(a)] for (a, b) in zip(θ, φ)])
c_s, r_s = NeuroAnalyzer._fiff_fit_sphere(sp_pts)
@test c_s ≈ [0.1, -0.02, 0.04] && r_s ≈ 0.09
@test_throws ArgumentError NeuroAnalyzer._fiff_fit_sphere(sp_pts[1:3, :])
# MNE-Python: info["dev_head_t"], fit_sphere_to_headshape(info, dig_kinds = ("eeg", "extra"))
R_t, t_t = NeuroAnalyzer._fiff_dev_head_t(mi)
@test R_t ≈ [0.9914199709892273 -0.039936363697052 -0.12446729093790054;
             0.06066115200519562 0.984011709690094 0.1674562394618988;
             0.11578965932130814 -0.17356973886489868 0.9779911637306213] rtol = 1e-6
@test t_t ≈ [-0.006129309069365263, 6.360560655593872e-05, 0.06474152207374573] rtol = 1e-6
@test NeuroAnalyzer._fiff_head_origin(mi) ≈ [-0.003430428724524391, 0.014118320507103961, 0.05298561185550042] atol = 1e-6
@test all(p -> length(p) == 3, NeuroAnalyzer._fiff_hsp(mi))
@test isnothing(NeuroAnalyzer._fiff_next_file(fobj, joinpath(testfiles_path, "meg", "test_raw.fif")))
@test NeuroAnalyzer._fiff_locs(mi, NeuroAnalyzer._fiff_labels(mi)) == import_locs_fiff(joinpath(testfiles_path, "meg", "test_raw.fif"))
# raw data: uncalibrated buffers = the data read by load_fiff()
buf, fb2, insts = open(NeuroAnalyzer._create_fiff_block, joinpath(testfiles_path, "meg", "test_raw.fif"))
raw_inst = NeuroAnalyzer._fiff_instances(fobj, "raw_data")[1]
rd = NeuroAnalyzer._fiff_raw_data(buf, fb2, insts, raw_inst, 376)
@test rd[1] == fiff[:raw_data][:data]
# HPI coils: MNE-Python get_chpi_info(): 83, 143, 203, 263, 323 Hz
fiff_m, _, _ = load_fiff(joinpath(testfiles_path, "meg", "SSS", "test_move_anon_raw.fif"))
hpi = NeuroAnalyzer._fiff_hpi(fiff_m[:meas_info])
@test hpi[:freqs] == [83.0, 143.0, 203.0, 263.0, 323.0]
@test size(hpi[:dig]) == (5, 3)

@info "Test: import_locs_fiff()"
l = import_locs_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fif"))
@test l isa DataFrame
@test all(in(labels(import_fiff(joinpath(testfiles_path, "meg", "meg-test-fiff.fif")))), l[!, :label])
# centered and scaled: the outermost channel has radius 1.2
@test maximum(sqrt.(l.loc_x .^ 2 .+ l.loc_y .^ 2 .+ l.loc_z .^ 2)) ≈ 1.2

@info "Test: import_digitrack()"
eeg = import_digitrack(joinpath(testfiles_path, "eeg", "eeg-test-digitrack.txt"))
@test eeg.header.recording[:file_type] == "Digitrack"
@test (sr(eeg), nchannels(eeg)) == (256, 24)
# first rows of the file: 230, 4,83e+03, … 193 / 339, …, -116
@test eeg.data[[1, 2, 24], 1, 1] == [230.0, 4830.0, 193.0]
@test eeg.data[24, 2, 1] == -116.0
# Digitrack exports all channels in μV
@test eeg.header.recording[:unit][24] == "μV"

@info "Test: import_csv()"
# both files: the EEGLAB recording, time in ms (256 Hz); Fp1 = 222.4075, 182.7696, …
eeg = import_csv(joinpath(testfiles_path, "eeg", "eeg-test_txch.csv.gz"))
@test eeg.header.recording[:file_type] == "CSV"
@test (sr(eeg), nchannels(eeg)) == (256, 24)
@test eeg.data[1, 1:2, 1] == [222.4075, 182.7696]
eeg2 = import_csv(joinpath(testfiles_path, "eeg", "eeg-test_chxt.csv.gz"))
@test eeg2.data == eeg.data
@test labels(eeg2) == labels(eeg)

@info "Test: import_npy()"
n = import_npy(joinpath(testfiles_path, "eeg", "eeg-test-npy.npy"); sampling_rate = 256)
# numpy.load: shape (16, 34888); a[0, :3], a[-1, -1]
@test size(n.data) == (16, 34888, 1)
@test n.data[1, 1:3, 1] == [-0.0023653548385797894, -1.4558388368327606, -28.302241687116272]
@test n.data[16, end, 1] == -2.646438793674498
@test sr(n) == 256
@test_throws ArgumentError import_npy(joinpath(testfiles_path, "eeg", "eeg-test-npy.npy"); sampling_rate = 0)

@info "Test: import_xdf()"
eeg = import_xdf(joinpath(testfiles_path, "eeg", "eeg-test-xdf.xdf"))
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "XDF"
@test sr(eeg) > 0 && size(eeg.data, 2) == length(eeg.time_pts)
@test all(m -> 0 <= m <= eeg.time_pts[end], eeg.markers[!, :start])

@info "Test: import_ncs()"
eeg = import_ncs(joinpath(testfiles_path, "eeg", "eeg-test-ncs.ncs"))
@test eeg.header.recording[:data_type] == "ieeg"
@test eeg.header.recording[:file_type] == "CSC"
@test nchannels(eeg) == 1
@test eeg.header.recording[:unit] == ["μV"]
@test length(eeg.time_pts) == signal_len(eeg)

@info "Test: import_snirf()"
# compared with the HDF5 contents of the file (nirs/data1/dataTimeSeries, time, stim1)
n = import_snirf(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf"))
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "SNIRF"
@test (sr(n), nchannels(n), signal_len(n)) == (50, 72, 500)
@test n.data[1:2, 1:3, 1] == [4.81625 1.834 1.834; 4.81625 1.74825 1.74825]
@test n.markers[1, :start] == 4.98
@test n.markers[1, :value] == "test"

@info "Test: import_nirs()"
n = import_nirs(joinpath(testfiles_path, "nirs", "fnirs-test-nirs.nirs"))
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "NIRS"
mat = NeuroAnalyzer.matread(joinpath(testfiles_path, "nirs", "fnirs-test-nirs.nirs"))
# 25.04 Hz (t) → resampled to 25 Hz, the same duration
@test 1 / (mat["t"][2] - mat["t"][1]) ≈ 25.040064102564102
@test sr(n) == 25
@test signal_len(n) == round(Int64, size(mat["d"], 1) * 25 / 25.040064102564102)
@test n.data[1:size(mat["d"], 2), 1, 1] ≈ mat["d"][1, :] rtol = 1e-3
# markers: one per non-zero stimulus sample; value = condition (column of s)
s_mat = mat["s"] isa AbstractVector ? reshape(mat["s"], :, 1) : mat["s"]
@test nrow(n.markers) == count(!iszero, s_mat)
@test all(in(string.(1:size(s_mat, 2))), n.markers[!, :value])

@info "Test: import_nirx()"
# MNE read_raw_nirx: the same samples (12.5 Hz, NeuroAnalyzer: exact ×2 upsampling to 25 Hz), channels
# grouped by wavelength (MNE: S1_D2 760, S1_D2 850, S1_D9 760, …)
n = import_nirx(joinpath(testfiles_path, "nirs", "nirx", "NIRS-2020-08-18_001.hdr"))
@test n.header.recording[:data_type] == "nirs"
@test n.header.recording[:file_type] == "NIRX"
@test (sr(n), nchannels(n), signal_len(n)) == (25, 26, 440)
@test n.header.recording[:optode_pairs][1:2, :] == [1 2; 1 9]
@test n.header.recording[:wavelengths] == [760.0, 850.0]
# MNE: onsets 0.0, 7.52, 10.64 s, conditions 4, 2, 1
@test n.markers[!, :start] == [0.0, 7.52, 10.64]
@test n.markers[!, :value] == ["4", "2", "1"]

@info "Test: import_ft()"
# compared with the MATLAB contents (data.trial{1})
mat = NeuroAnalyzer.matread(joinpath(testfiles_path, "eeg", "eeg-test-fieldtrip.mat"))["data"]
eeg = import_ft(joinpath(testfiles_path, "eeg", "eeg-test-fieldtrip.mat"); type = :eeg)
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "FT"
@test (sr(eeg), nchannels(eeg), signal_len(eeg)) == (256, 24, 308480)
@test eeg.data[:, :, 1] == mat["trial"][1]
@test labels(eeg)[1] == "Fp1"
meg = import_ft(joinpath(testfiles_path, "meg", "meg-test-fieldtrip.mat"); type = :meg)
@test meg.header.recording[:data_type] == "meg"
@test meg.header.recording[:file_type] == "FT"
# the same recording as meg-test-fiff.fif: T → fT, T/m → fT/cm
@test meg.data[12, pos(32000), 1] ≈ [37475.75757133061, 47371.71169613672, 37014.06476161138, -2739.7155741580623, -558.0902095507164] rtol = 1e-9
@test meg.data[13, pos(32000), 1] ≈ [-73.67706281630497, 193.15338089679952, 93.58978249638739, -304.6646111052611, 205.10101270484896] rtol = 1e-9
# fNIRS: 25.04 Hz → resampled to 25 Hz (the same duration)
nirs = import_ft(joinpath(testfiles_path, "nirs", "fnirs-test-fieldtrip.mat"); type = :nirs)
@test nirs.header.recording[:data_type] == "nirs"
@test nirs.header.recording[:file_type] == "FT"
@test sr(nirs) == 25
@test signal_len(nirs) == round(Int64, 9015 * 25 / 25.040064102561335)
m = import_ft(joinpath(testfiles_path, "events", "events-test-fieldtrip.mat"); type = :events)
@test m isa DataFrame
# start and length in samples
@test m[1:3, :start] == [0.0, 10528.0, 10742.0]
@test m[1:3, :value] == ["Impedance", "s111", "s130"]
m_new = markers_s2t(m; fs = 256)
@test m_new[2, :start] == NeuroAnalyzer.s2t(10528, 256)
@test m[2, :start] == 10528.0
eeg.markers = m
markers_s2t!(eeg)
@test eeg.markers[2, :start] == m_new[2, :start]
@test_throws ArgumentError import_ft(joinpath(testfiles_path, "eeg", "eeg-test-fieldtrip.mat"); type = :x)

@info "Test: import_duomag()"
mep = import_duomag(joinpath(testfiles_path, "mep", "mep-duomag.m"))
@test mep.header.recording[:data_type] == "mep"
@test mep.header.recording[:file_type] == "DuoMAG"
mep2 = import_duomag(joinpath(testfiles_path, "mep", "mep-duomag.ascii"))
@test mep2.header.recording[:data_type] == "mep"
# the same recording in both formats
@test size(mep2.data) == size(mep.data)
@test sr(mep2) == sr(mep)
@test mep2.data ≈ mep.data rtol = 1e-3

@info "Test: import_thymatron()"
eeg = import_thymatron(joinpath(testfiles_path, "thymatron/1.png"))
@test eeg.header.recording[:data_type] == "eeg"
@test eeg.header.recording[:file_type] == "Thymatron"
@test nchannels(eeg) >= 1 && signal_len(eeg) > 0

@info "Test: import_dat()"
d = import_dat(joinpath(testfiles_path, "eeg", "eeg-test-dat.dat"))
# columns from the file header: Trial, Resp, Type, Correct, Latency (ms); 320 trials (2 … 321)
@test names(d) == ["trial", "response", "type", "correct", "latency"]
@test nrow(d) == 320
@test d[1, :] == (trial = 2, response = 0, type = 4, correct = 1, latency = 1600)
@test d[4, [:trial, :response, :type, :latency]] == (trial = 5, response = 2, type = 2, latency = 282)
@test d[end, :trial] == 321
@test_throws ArgumentError import_dat(joinpath(testfiles_path, "locs", "locs.ced"))

@info "Test: import_nwb()"
# minimal NWB-like file: acquisition/EEG data + timestamps, NWB root datasets
nwb = joinpath(tmp, "rec_eeg.nwb")
NeuroAnalyzer.JLD2.jldopen(nwb, "w") do f
    f["acquisition/EEG/data"] = [1.0 2.0 3.0 4.0; 5.0 6.0 7.0 8.0]
    f["acquisition/EEG/timestamps"] = [0.0, 0.01, 0.02, 0.03]
    f["file_create_date"] = "2020-01-01T10:00:00+00:00"
    f["identifier"] = "sub-01"
    f["session_description"] = "test"
end
write(joinpath(tmp, "rec_events.tsv"), "onset\tduration\ttrial_type\tvalue\n0.01\t0\tstim\ta\n")
n = import_nwb(nwb)
@test n.header.recording[:file_type] == "NWB"
@test sr(n) == 100
@test n.data[:, :, 1] == [1.0 2.0 3.0 4.0; 5.0 6.0 7.0 8.0]
@test n.header.subject[:id] == "sub-01"
@test n.header.recording[:recording_date] == "2020-01-01"
@test n.markers[!, :start] == [0.01]
@test n.markers[!, :value] == ["a"]
@test_throws ArgumentError import_nwb(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))

@info "Test: import_recording()"
eeg = import_recording(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
@test eeg.data == import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf")).data
eeg = import_recording(joinpath(testfiles_path, "eeg", "eeg-test-bv.vhdr"))
@test eeg.header.recording[:file_type] == "BrainVision"
n = import_recording(joinpath(testfiles_path, "eeg", "eeg-test-npy.npy"); sampling_rate = 256)
@test sr(n) == 256
@test_throws ArgumentError import_recording(joinpath(testfiles_path, "eeg", "eeg-test-npy.npy"))
@test import_recording(joinpath(testfiles_path, "eeg", "eeg-test-fieldtrip.mat"); type = :eeg) isa NeuroAnalyzer.NEURO
@test import_recording(joinpath(testfiles_path, "events", "events-test-fieldtrip.mat"); type = :events) isa DataFrame
eeg = import_recording(joinpath(tmp, "edf.bdf"); format = :edf)
@test eeg.header.recording[:file_type] == "EDF"
@test_throws ArgumentError import_recording(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"); format = :xyz)
@test_throws ArgumentError import_recording(joinpath(tmp, "no_such_file.edf"))

@info "Test: import_locs()"
# EEGLAB polar convention (0° = nose, positive to the right) → NeuroAnalyzer (0° = right, 90° = nose):
# θ = 90° - θ_EEGLAB; x+ right, y+ nose, z+ up
l = import_locs_ced(joinpath(testfiles_path, "locs", "locs.ced"))
@test nrow(l) == 19
# Fp1: theta -18, X 0.95 (nose), Y 0.309 (left), Z -0.0349
@test l[1, :label] == "Fp1"
@test l[1, :loc_theta] == 108.0
@test [l[1, :loc_x], l[1, :loc_y], l[1, :loc_z]] ≈ [-0.309, 0.95, -0.0349] rtol = 1e-2
locs2 = import_locs_locs(joinpath(testfiles_path, "locs", "locs.locs"))
@test nrow(locs2) == 81
@test locs2[1, :label] == "Fp1" && locs2[1, :loc_theta] == 108.0
# spherical coordinates from the polar ones (radius 0.511 → elevation -2°): the same position as in .ced
@test locs2[1, :loc_theta_sph] == 108.0
@test locs2[1, :loc_phi_sph] ≈ 90 - 180 * 0.51111
@test [locs2[1, :loc_x], locs2[1, :loc_y], locs2[1, :loc_z]] ≈ [l[1, :loc_x], l[1, :loc_y], l[1, :loc_z]] rtol = 1e-2
l = import_locs_tsv(joinpath(testfiles_path, "locs", "locs.tsv"))
# C3: x -0.5878, y 0, z 0.8090 (left)
@test l[1, :label] == "C3"
@test [l[1, :loc_x], l[1, :loc_y], l[1, :loc_z]] ≈ [-0.5878, 0.0, 0.809] rtol = 1e-3 atol = 1e-6
@test l[1, :loc_theta] == 180.0
l = import_locs_csd(joinpath(testfiles_path, "locs", "locs.csd"))
# Nose: theta 90 (front), phi -33.75, radius 1; Cartesian 0, 0.8315, -0.5556
@test nrow(l) == 465
@test l[1, :label] == "Nose"
@test (l[1, :loc_theta_sph], l[1, :loc_phi_sph], l[1, :loc_radius_sph]) == (90.0, -33.75, 1.0)
@test [l[1, :loc_x], l[1, :loc_y], l[1, :loc_z]] ≈ [0.0, 0.8315, -0.5556] atol = 1e-3
l = import_locs_txt(joinpath(testfiles_path, "locs", "locs.txt"))
@test nrow(l) == 67
@test l[1, :label] == "Fp1" && round(l[1, :loc_theta]) == 108.0
for (f, fn, nr) in (("locs.elc", import_locs_elc, 68), ("locs.sfp", import_locs_sfp, 36), ("locs.geo", import_locs_geo, 80),
                    ("locs.mat", import_locs_mat, 63), ("locs.dat", import_locs_dat, 120), ("locs.asc", import_locs_asc, 128))
    lf = fn(joinpath(testfiles_path, "locs", f))
    @test lf isa DataFrame
    @test nrow(lf) == nr
    @test names(lf) == ["label", "loc_radius", "loc_theta", "loc_x", "loc_y", "loc_z", "loc_radius_sph", "loc_theta_sph", "loc_phi_sph"]
    # polar and Cartesian coordinates agree
    @test lf.loc_radius .* cosd.(lf.loc_theta) ≈ lf.loc_x .* lf.loc_radius ./ max.(hypot.(lf.loc_x, lf.loc_y), eps()) atol = 1e-3
end
# the left frontal electrodes lie in the left frontal quadrant
l = import_locs_geo(joinpath(testfiles_path, "locs", "locs.geo"))
@test l[1, :label] == "Fp1" && l[1, :loc_x] < 0 && l[1, :loc_y] > 0
l = import_locs_sfp(joinpath(testfiles_path, "locs", "locs.sfp"))
@test l[1, :label] == "FidNz"
# generic importer: by extension
@test import_locs(joinpath(testfiles_path, "locs", "locs.ced")) == import_locs_ced(joinpath(testfiles_path, "locs", "locs.ced"))
@test_throws ArgumentError import_locs(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
@test_throws ArgumentError import_locs_ced(joinpath(testfiles_path, "locs", "locs.locs"))

@info "Test: load_locs()"
eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
eeg = load_locs(eeg; file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"))
# every locatable channel of the recording (EEG, EOG, reference, …) gets a row
@test all(in(labels(eeg)), eeg.locs[!, :label])
@test all(in(eeg.locs[!, :label]), get_channel(eeg; type = "eeg"))
eeg_tmp = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
load_locs!(eeg_tmp; file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"))
@test eeg_tmp.locs == eeg.locs

@info "Test: export_locs()"
@test_throws ArgumentError export_locs(eeg; file_name = joinpath(tmp, "locs_out.sfp"))
for ext in ("ced", "locs", "tsv", "csv")
    f = joinpath(tmp, "locs_out.$ext")
    export_locs(eeg; file_name = f)
    @test isfile(f)
    @test_throws ArgumentError export_locs(eeg; file_name = f)
    export_locs(eeg; file_name = f, overwrite = true)
end
# round trips: the coordinates are kept
l = import_locs_ced(joinpath(tmp, "locs_out.ced"))
@test l[!, :label] == eeg.locs[!, :label]
# the same directions (the importer normalizes the radius)
unitrows(m) = m ./ sqrt.(sum(abs2, m; dims = 2))
@test unitrows(Matrix(l[!, [:loc_x, :loc_y, :loc_z]])) ≈ unitrows(Matrix(eeg.locs[!, [:loc_x, :loc_y, :loc_z]])) atol = 1e-2
@test mod.(l[!, :loc_theta], 360) ≈ mod.(eeg.locs[!, :loc_theta], 360) atol = 1e-2
l = import_locs_locs(joinpath(tmp, "locs_out.locs"))
@test mod.(l[!, :loc_theta], 360) ≈ mod.(eeg.locs[!, :loc_theta], 360) atol = 1e-2
l = import_locs_tsv(joinpath(tmp, "locs_out.tsv"))
@test unitrows(Matrix(l[!, [:loc_x, :loc_y, :loc_z]])) ≈ unitrows(Matrix(eeg.locs[!, [:loc_x, :loc_y, :loc_z]])) atol = 1e-2
l = import_locs_csv(joinpath(tmp, "locs_out.csv"))
@test l == eeg.locs
@test_throws ArgumentError import_locs_csv(joinpath(tmp, "no_such_file.csv"))

@info "Test: save_neuro() / load_neuro()"
f = joinpath(tmp, "test.hdf")
save_neuro(eeg; file_name = f)
@test isfile(f)
@test_throws ArgumentError save_neuro(eeg; file_name = f)
save_neuro(eeg; file_name = f, overwrite = true)
new = load_neuro(f)
@test new.data == eeg.data
@test new.header.recording == eeg.header.recording
@test new.locs == eeg.locs
@test new.markers == eeg.markers
@test new.time_pts == eeg.time_pts
@test_throws ArgumentError load_neuro(joinpath(tmp, "no_such_file.hdf"))
# file shipped with the test data: an EEG object saved by an old version of NeuroAnalyzer → clear error
@test_throws ArgumentError load_neuro(joinpath(testfiles_path, "eeg", "eeg.hdf5"))

@info "Test: export_csv()"
f = joinpath(tmp, "eeg.csv")
export_csv(eeg; file_name = f, header = true, epoch_time = true, locs = true)
@test isfile(f)
@test isfile(joinpath(tmp, "eeg_header.txt"))
@test isfile(joinpath(tmp, "eeg_locs.csv"))
@test_throws ArgumentError export_csv(eeg; file_name = f)
# round trip: time × channels, with channel names
eeg_csv = import_csv(f)
@test labels(eeg_csv) == labels(eeg)
@test sr(eeg_csv) == sr(eeg)
@test eeg_csv.data ≈ eeg.data

@info "Test: export_markers()"
eeg = import_bdf(joinpath(testfiles_path, "eeg", "eeg-test-bdfplus.bdf"))
f = joinpath(tmp, "markers.csv")
export_markers(eeg; file_name = f)
@test isfile(f)
@test_throws ArgumentError export_markers(eeg; file_name = f)
m = NeuroAnalyzer.CSV.read(f, DataFrame; stringtype = String)
@test m[!, :start] == eeg.markers[!, :start]
@test m[!, :value] == eeg.markers[!, :value]

@info "Test: import_montage()"
ref_list, ref_name = import_montage(joinpath(NeuroAnalyzer.PATH, "montages", "bip_long.mnt"))
@test ref_list == [
    "Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3",
    "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2",
]
@test ref_name == "longitudinal-BIP"
@test_throws ArgumentError import_montage(joinpath(tmp, "no_such_file.mnt"))

rm(tmp; recursive = true)

true
