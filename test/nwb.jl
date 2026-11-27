using Test
using DataFrames

# NWB test files (nwb/README.md): synthetic (pynwb), DANDI 000576 (trials, coordinates),
# DANDI 001347 (sEEG depth electrodes)
H = NeuroAnalyzer.HDF5
nwb_root = joinpath(testfiles_path, "nwb")
f_syn = joinpath(nwb_root, "synthetic", "synthetic_eeg.nwb")
f_ndx = joinpath(nwb_root, "synthetic", "synthetic_ndx_events.nwb")
f_576 = joinpath(nwb_root, "dandi-000576", "sub-01", "sub-01_ses-20140828T132700_ecephys.nwb")
f_1347 = joinpath(nwb_root, "dandi-001347", "sub-SD186", "sub-SD186_ses-puffControl_ecephys.nwb")

@info "Test: NWB objects"
H.h5open(f_syn, "r") do f
    o = NeuroAnalyzer._nwb_objects(f)
    @test NeuroAnalyzer._nwb_type(f) == ("NWBFile", "core")
    es = NeuroAnalyzer._nwb_objects(f; types = ["ElectricalSeries"])
    @test sort(es.path) ==
          ["/acquisition/EEG", "/acquisition/EOG", "/processing/ecephys/LFP/EEG_downsampled"]
    # the device is linked from both electrode groups: listed once
    @test count(==("Device"), o.type) == 1 && !any(startswith("/specifications"), o.path)
    @test sort(NeuroAnalyzer._nwb_objects(f; types = ["EventsTable"]).path) ==
          ["/events/button_presses", "/events/stimuli"]
end

@info "Test: NWB tables"
H.h5open(f_syn, "r") do f
    e = NeuroAnalyzer._nwb_table(f, f["general/extracellular_ephys/electrodes"])
    @test nrow(e) == 82 && names(e)[1] == "id"
    @test e.label[1:3] == ["Fp1", "Fpz", "Fp2"] && e.x[1] ≈ -28.296923
    @test eltype(e.bad) == Bool && e.label[e.bad] == ["T7"]
    # references: paths of the electrode groups
    @test unique(e.group) == ["/general/extracellular_ephys/EEG", "/general/extracellular_ephys/EOG"]
    t = NeuroAnalyzer._nwb_table(f, f["intervals/trials"])
    @test nrow(t) == 22 && count(==("target"), t.condition) == 5
    @test all(isnan, t.response_time[t.condition .== "standard"])
    @test t.stop_time .- t.start_time ≈ ones(22)
    s = NeuroAnalyzer._nwb_table(f, f["events/stimuli"])
    @test s.timestamp ≈ t.start_time && s.stimulus[t.condition .== "target"] == fill("tone 2000 Hz", 5)
    # ragged column (tags of invalid_times)
    bad_t = NeuroAnalyzer._nwb_table(f, f["intervals/invalid_times"])
    @test bad_t.tags == [["movement"]] && bad_t.start_time == [12.0]
end
# legacy ndx-events: ragged event times, LabeledEvents labels
H.h5open(f_ndx, "r") do f
    a = NeuroAnalyzer._nwb_table(f, f["acquisition/annotated_events"])
    @test a.label == ["seizure onset", "artifact"] && a.event_times == [[8.2], [3.3, 6.6]]
    @test NeuroAnalyzer._nwb_type(f["acquisition/stimulus_codes"]) == ("LabeledEvents", "ndx-events")
    @test NeuroAnalyzer._nwb_attr(f["acquisition/stimulus_codes/data"], "labels") == ["", "face", "house"]
end
# compound references (trials → signal), DANDI 000576
H.h5open(f_576, "r") do f
    t = NeuroAnalyzer._nwb_table(f, f["intervals/trials"])
    @test nrow(t) == 17 && Set(t.condition) == Set(["Neutral", "Aversive"])
    @test t.timeseries[1][1].timeseries == "/processing/ecephys/LFP/ecephys.ieeg"
    @test t.timeseries[1][1].count == 52000
    e = NeuroAnalyzer._nwb_table(f, f["general/extracellular_ephys/electrodes"])
    @test e.label == ["AL2-AL1", "AR2-AR1"] && e.is_inside_soz == [0, 1]
end

@info "Test: NWB series"
H.h5open(f_syn, "r") do f
    s = NeuroAnalyzer._nwb_series(f, f["acquisition/EEG"])
    @test s.type == "ElectricalSeries" && size(s.data) == (80, 7500) && s.rate == 250.0
    @test s.electrodes == 1:80 && s.electrodes_table == "/general/extracellular_ephys/electrodes"
    @test s.unit == "volts" && s.conversion == 1e-7 && s.channel_conversion[1:4] == [2, 2, 2, 1]
    raw = NeuroAnalyzer._nwb_series(f, f["acquisition/EEG"]; scaled = false)
    @test eltype(raw.data) == Int16
    # scaled: value × channel_conversion × conversion (V)
    @test s.data[1, :] ≈ Float64.(raw.data[1, :]) .* 2 .* 1e-7
    @test s.data[4, :] ≈ Float64.(raw.data[4, :]) .* 1e-7
    @test 10e-6 < maximum(abs.(s.data[1, :])) < 500e-6
    eog = NeuroAnalyzer._nwb_series(f, f["acquisition/EOG"])
    @test isnothing(eog.rate) && length(eog.timestamps) == 7500 && size(eog.data) == (2, 7500)
    @test eog.electrodes == [81, 82] && !issorted(diff(eog.timestamps))
    r = NeuroAnalyzer._nwb_series(f, f["acquisition/respiration"])
    @test r.type == "TimeSeries" && r.rate == 25.0 && size(r.data) == (750,) && isnothing(r.electrodes)
end
H.h5open(f_1347, "r") do f
    es = NeuroAnalyzer._nwb_objects(f; types = ["ElectricalSeries"])
    @test es.path == ["/processing/ecephys/LFP/v_trace"]
    s = NeuroAnalyzer._nwb_series(f, f[es.path[1]])
    @test size(s.data) == (88, 240000) && s.rate == 200.0 && length(s.electrodes) == 88
    e = NeuroAnalyzer._nwb_table(f, f[s.electrodes_table])
    @test length(unique(e.group)) == 7 && all(isnan, e.x)
end

@info "Test: nwb_series(), import_nwb(; series)"
l = nwb_series(f_syn)
@test l.path == ["/acquisition/EEG", "/acquisition/EOG", "/processing/ecephys/LFP/EEG_downsampled",
                 "/acquisition/respiration"]
@test l.channels == [80, 2, 80, 1] && l.rate[[1, 3, 4]] == [250.0, 125.0, 25.0] && isnan(l.rate[2])
@test l.unit[1] == "volts" && l.type[4] == "TimeSeries"
# default: the first ElectricalSeries in acquisition
n1 = import_nwb(f_syn)
@test nchannels(n1) == 80 && sr(n1) == 250 && signal_len(n1) == 7500
@test labels(n1)[1:3] == ["Fp1", "Fpz", "Fp2"]
# volts → μV (data × channel_conversion × conversion)
H.h5open(f_syn, "r") do f
    raw = read(f["acquisition/EEG/data"])
    @test n1.data[1, :, 1] ≈ Float64.(raw[1, :]) .* 2 .* 1e-7 .* 1e6
end
# by name and by path
n2 = import_nwb(f_syn; series = "EEG_downsampled")
@test nchannels(n2) == 80 && sr(n2) == 125
@test import_nwb(f_syn; series = "/processing/ecephys/LFP/EEG_downsampled").data == n2.data
# irregular timestamps: rate from the timestamps
n3 = import_nwb(f_syn; series = "EOG")
@test nchannels(n3) == 2 && labels(n3) == ["VEOG", "HEOG"] && sr(n3) == 250
# a plain TimeSeries
n4 = import_nwb(f_syn; series = "respiration")
@test nchannels(n4) == 1 && sr(n4) == 25 && signal_len(n4) == 750
@test_throws ArgumentError import_nwb(f_syn; series = "x")
# signal in processing/ecephys (no ElectricalSeries in acquisition)
l2 = nwb_series(f_1347)
@test l2.path == ["/processing/ecephys/LFP/v_trace", "/acquisition/dc_trace"]
n5 = import_nwb(f_1347)
@test nchannels(n5) == 88 && sr(n5) == 200 && labels(n5)[1] == "LaTH1"
@test_throws ArgumentError nwb_series(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
@test NeuroAnalyzer._nwb_datetime("2016-12-17T12:02:00+01:00") == ("2016-12-17", "11:02:00")
@test NeuroAnalyzer._nwb_datetime("2025-02-28T16:19:55.220891-05:00") == ("2025-02-28", "21:19:55")
@test NeuroAnalyzer._nwb_datetime("2020-01-01") == ("2020-01-01", "")

@info "Test: nwb_series(), import_nwb(; series)"
l = nwb_series(f_syn)
@test l.path == ["/acquisition/EEG", "/acquisition/EOG", "/processing/ecephys/LFP/EEG_downsampled",
                 "/acquisition/respiration"]
@test l.channels == [80, 2, 80, 1] && l.samples == [7500, 7500, 3750, 750]
@test l.rate[[1, 3, 4]] == [250.0, 125.0, 25.0] && isnan(l.rate[2])
@test l.unit[1] == "volts" && l.type == ["ElectricalSeries", "ElectricalSeries", "ElectricalSeries", "TimeSeries"]
# default: the first ElectricalSeries in acquisition
n1 = import_nwb(f_syn)
@test nchannels(n1) == 80 && sr(n1) == 250 && signal_len(n1) == 7500
@test labels(n1)[1:3] == ["Fp1", "Fpz", "Fp2"] && n1.header.subject[:id] == "synthetic01"
@test (n1.header.recording[:recording_date], n1.header.recording[:recording_time]) ==
      ("2020-01-01", "12:00:00")
# volts → μV (data × channel_conversion × conversion)
H.h5open(f_syn, "r") do f
    raw = read(f["acquisition/EEG/data"])
    @test n1.data[1, :, 1] ≈ Float64.(raw[1, :]) .* 2 .* 1e-7 .* 1e6
    @test n1.data[4, :, 1] ≈ Float64.(raw[4, :]) .* 1e-7 .* 1e6
end
# by name and by path
n2 = import_nwb(f_syn; series = "EEG_downsampled")
@test nchannels(n2) == 80 && sr(n2) == 125 && signal_len(n2) == 3750
@test import_nwb(f_syn; series = "/processing/ecephys/LFP/EEG_downsampled").data == n2.data
@test import_nwb(f_syn; series = "processing/ecephys/LFP/EEG_downsampled").data == n2.data
# timestamps: the sampling rate is estimated
n3 = import_nwb(f_syn; series = "EOG")
@test nchannels(n3) == 2 && labels(n3) == ["VEOG", "HEOG"] && sr(n3) == 250
# a plain TimeSeries
n4 = import_nwb(f_syn; series = "respiration")
@test nchannels(n4) == 1 && sr(n4) == 25 && signal_len(n4) == 750 && labels(n4) == ["ch_1"]
@test_throws ArgumentError import_nwb(f_syn; series = "x")
# signal in processing/ecephys, a TimeSeries in acquisition
l2 = nwb_series(f_1347)
@test l2.path == ["/processing/ecephys/LFP/v_trace", "/acquisition/dc_trace"]
n5 = import_nwb(f_1347)
@test nchannels(n5) == 88 && sr(n5) == 200 && labels(n5)[1] == "LaTH1"
@test n5.header.subject[:id] == "SD186"
# timestamps (2 kHz), DANDI 000576
n6 = import_nwb(f_576)
@test labels(n6) == ["AL2-AL1", "AR2-AR1"] && sr(n6) == 2000 && signal_len(n6) == 884000
@test_throws ArgumentError nwb_series(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
@test NeuroAnalyzer._nwb_datetime("2016-12-17T12:02:00+01:00") == ("2016-12-17", "11:02:00")
@test NeuroAnalyzer._nwb_datetime("2025-02-28T16:19:55.220891-05:00") == ("2025-02-28", "21:19:55")
@test NeuroAnalyzer._nwb_datetime("2020-01-01") == ("2020-01-01", "")

@info "Test: import_nwb(): electrodes, resampling, metadata, events"
rec1 = n1.header.recording
@test datatype(n1) == "eeg" && all(==("eeg"), rec1[:channel_type])
@test labels(n1)[rec1[:bad_channel]] == ["T7"]
# locations: MNI coordinates (RAS); filters from "0.1-100 Hz"
@test nrow(n1.locs) == 80 && n1.locs.loc_x[1] ≈ -28.296923 && n1.locs.loc_y[1] ≈ 81.688141
@test locs_coord(n1).unit == "mm" && locs_coord(n1).guessed
@test all(rec1[:highpass] .≈ 0.1) && all(rec1[:lowpass] .≈ 100)
@test occursin("Subject: sex U, age P30Y, species Homo sapiens", n1.header.notes)
@test occursin("Experimenter: NeuroAnalyzer", n1.header.notes)
# NWB schema axes (x: posterior, y: inferior, z: right) → RAS
n1p = import_nwb(f_syn; locs_axes = :pir, events = false)
@test n1p.locs.loc_x ≈ n1.locs.loc_z && n1p.locs.loc_y ≈ -n1.locs.loc_x && n1p.locs.loc_z ≈ -n1.locs.loc_y
@test nrow(n1p.markers) == 0
@test_throws ArgumentError import_nwb(f_syn; locs_axes = :x)
# events: trials (condition, response time), EventsTable (stimuli, button presses),
# invalid times, annotations
mk = n1.markers
@test nrow(mk) == 53 && issorted(mk.start)
@test count(==("BAD"), mk.type) == 1 && mk.value[mk.type .== "BAD"] == ["BAD_movement"]
@test mk.start[mk.type .== "BAD"] == [12.0] && mk.length[mk.type .== "BAD"] == [1.0]
@test count(==("response"), mk.type) == 5 && count(==("annotation"), mk.type) == 3
tr = mk[.!ismissing.(mk.condition), :]
@test nrow(tr) == 22 && tr.value == tr.condition && all(tr.length .≈ 1)
@test count(==("target"), tr.value) == 5 && all(isnan, tr.response_time[tr.value .== "standard"])
@test count(==("tone 2000 Hz"), mk.value) == 5
# irregular timestamps: resampled to 250 Hz (linear interpolation)
@test datatype(n3) == "eeg" && all(==("eog"), n3.header.recording[:channel_type])
H.h5open(f_syn, "r") do f
    ts = read(f["acquisition/EOG/timestamps"])
    v = Float64.(read(f["acquisition/EOG/data"])) .* 1e6
    # samples of the regular time axis between the original samples
    t = (0:(signal_len(n3) - 1)) ./ 250
    k = [clamp(searchsortedlast(ts, x), 1, length(ts) - 1) for x in t]
    w = (t .- ts[k]) ./ (ts[k .+ 1] .- ts[k])
    @test n3.data[1, :, 1] ≈ (1 .- w) .* v[1, k] .+ w .* v[1, k .+ 1]
end
# depth electrodes (group descriptions), ndx-events 0.2 events
n7 = import_nwb(f_ndx)
@test datatype(n7) == "seeg" && labels(n7)[1:2] == ["LHIP1", "LHIP2"] && nrow(n7.locs) == 8
m7 = n7.markers
@test nrow(m7) == 14
@test sort(unique(m7.value)) == ["artifact", "face", "house", "seizure onset", "spikes_detected", "start", "stop"]
@test m7.start[m7.value .== "face"] == [1.0, 3.0, 4.0, 6.0]
@test m7.hemisphere[m7.value .== "artifact"] == ["both", "both"]
@test m7.start[m7.value .== "artifact"] == [3.3, 6.6]
# iEEG (000576): reference, filters, trials with lengths
@test datatype(n6) == "ieeg" && n6.header.recording[:reference] == "Common intracranial reference"
@test all(n6.header.recording[:highpass] .≈ 0.5)
@test nrow(n6.markers) == 17 && n6.markers.length[1:2] == [26.0, 24.0]
@test n6.markers.value[1:2] == ["Neutral", "Aversive"]
@test (n6.header.recording[:recording_date], n6.header.recording[:recording_time]) ==
      ("2014-08-28", "11:27:00")
# sEEG shanks (001347): no coordinates, no events
@test datatype(n5) == "seeg" && nrow(n5.markers) == 0 && locs_coord(n5).unit == ""
@test occursin("Lab: Karl Deisseroth Lab", n5.header.notes)
# helpers
@test NeuroAnalyzer._nwb_filtering("0.1-100 Hz") == "HP: 0.1 Hz LP: 100 Hz"
@test NeuroAnalyzer._nwb_filtering("Passband, 0.5 to 1000 Hz") == "HP: 0.5 Hz LP: 1000 Hz"
@test NeuroAnalyzer._nwb_filtering("unknown") == "unknown"
@test NeuroAnalyzer._nwb_ch_type("LAMY1", "LAMY electrode group for shank LAMY") == "seeg"
@test NeuroAnalyzer._nwb_ch_type("G1", "grid subdural") == "ecog"
@test NeuroAnalyzer._nwb_ch_type("ECG", "EEG 10-10 EEG cap scalp") == "ecg"
@test NeuroAnalyzer._nwb_ch_type("A1", "ieeg iEEG electrodes Intracranial") == "ieeg"
@test NeuroAnalyzer._nwb_ch_type("X1", "") == "eeg"
@test NeuroAnalyzer._nwb_ch_type("A1", "EEG 10-10 EEG cap scalp") == "ref"
@test NeuroAnalyzer._nwb_datatype(["seeg", "ecog"]) == "ieeg"
@test NeuroAnalyzer._nwb_datatype(["ecog", "eog"]) == "ecog"
@test NeuroAnalyzer._nwb_datatype(["eog"]) == "eeg"
@test NeuroAnalyzer._nwb_resample([0.0 1.0 3.0], [0.0, 0.1, 0.4], 10) ≈ [0.0 1.0 1.6666666666666667 2.3333333333333335 3.0]

@info "Test: import_recording() (NWB)"
@test all(f -> detect_format(f) === :nwb, (f_syn, f_ndx, f_576, f_1347))
@test import_recording(f_syn).data == n1.data
@test import_recording(f_syn; series = "EOG").data == n3.data
# NWB 2.x identified by the root neurodata type only (no NWB root datasets)
fr = tempname() * ".h5"
H.h5open(fr, "w") do f
    NeuroAnalyzer._nwb_attrs!(f, "NWBFile")
    f["x"] = [1.0]
end
@test detect_format(fr) === :nwb
@test_throws ArgumentError import_recording(fr)
rm(fr)

@info "Test: export_nwb()"
e_x = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edfplus.edf"))
NeuroAnalyzer.filter!(e_x; ch = "all", fprototype = :butterworth, ftype = :hp, cutoff = 0.5, order = 4)
e_x.header.recording[:bad_channel][3] = true
e_x.header.recording[:notch][1] = [50.0]
add_marker!(e_x; type = "response", start = 1.0, len = 0.5, value = "a, b", ch = labels(e_x)[3])
add_marker!(e_x; type = "BAD", start = 2.0, len = 1.0, value = "noise")
load_locs!(e_x; file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"))
fx = tempname() * ".nwb"
export_nwb(e_x; file_name = fx, subject_id = "S01", sex = "F", age = "P35Y",
           experimenter = ["Doe, Jane"], institution = "Test")
el_x = findall(u -> u in ("μV", "mV", "V"), e_x.header.recording[:unit])
b = import_nwb(fx)
@test labels(b) == labels(e_x)[el_x] && sr(b) == sr(e_x) && signal_len(b) == signal_len(e_x)
@test b.data ≈ e_x.data[el_x, :, :] rtol = 1e-6
@test b.header.recording[:channel_type] == e_x.header.recording[:channel_type][el_x]
@test b.header.recording[:bad_channel] == e_x.header.recording[:bad_channel][el_x]
@test b.header.recording[:highpass] ≈ e_x.header.recording[:highpass][el_x]
@test b.header.recording[:notch][1] == [50.0]
@test nrow(b.locs) == nrow(e_x.locs) && b.locs.loc_x ≈ e_x.locs.loc_x
@test (b.header.recording[:recording_date], b.header.recording[:recording_time]) == ("1988-09-02", "09:11:38")
@test b.header.subject[:id] == "S01" && occursin("sex F, age P35Y", b.header.notes)
mb = sort(b.markers, :start)
me = sort(e_x.markers, :start)
@test mb.start ≈ me.start && mb.length ≈ me.length && mb.value == me.value
@test mb.type == me.type && mb.channel == me.channel
# channels without voltage units: auxiliary TimeSeries
l = nwb_series(fx)
@test l.path == ["/acquisition/ElectricalSeries", "/acquisition/auxiliary_no_unit"]
@test l.channels == [length(el_x), nchannels(e_x) - length(el_x)]
@test import_nwb(fx; series = "auxiliary_no_unit").data[:, :, 1] ≈
      e_x.data[setdiff(1:nchannels(e_x), el_x), :, 1] rtol = 1e-6
# NWB file structure: root, cached schema, links
H.h5open(fx, "r") do f
    @test NeuroAnalyzer._nwb_type(f) == ("NWBFile", "core")
    @test NeuroAnalyzer._nwb_attr(f, "nwb_version") == NeuroAnalyzer._NWB_EXPORT_VERSION
    @test NeuroAnalyzer._nwb_ref_path(f, NeuroAnalyzer._nwb_attr(f, ".specloc")) == "/specifications"
    @test sort(keys(f["specifications"])) == ["core", "hdmf-common", "hdmf-experimental"]
    @test keys(f["specifications/core"]) == [NeuroAnalyzer._NWB_EXPORT_VERSION]
    @test all(in(keys(f)), ["acquisition", "analysis", "processing", "stimulus", "general",
                            "file_create_date", "identifier", "session_start_time"])
    @test read(f["general/subject/sex"]) == "F" && read(f["general/experimenter"]) == ["Doe, Jane"]
    d = f["general/extracellular_ephys/EEG/device"]
    @test NeuroAnalyzer._nwb_type(d)[1] == "Device"
    @test NeuroAnalyzer._nwb_objid(d) == NeuroAnalyzer._nwb_objid(f["general/devices/amplifier"])
    @test NeuroAnalyzer._nwb_type(f["events/markers/timestamp"]) == ("TimestampVectorData", "core")
    @test read(f["acquisition/ElectricalSeries/channel_conversion"])[1] == 1e-6
end
# epochs: intervals/trials
e_ep = epoch(e_x; ep_len = 10)
fe = tempname() * ".nwb"
export_nwb(e_ep; file_name = fe)
H.h5open(fe, "r") do f
    t = NeuroAnalyzer._nwb_table(f, f["intervals/trials"])
    @test nrow(t) == nepochs(e_ep) && t.start_time[2] == 10.0 && t.stop_time[1] == 10.0
end
@test signal_len(import_nwb(fe; events = false)) == epoch_len(e_ep) * nepochs(e_ep)
@test export_recording(e_x; file_name = fe, overwrite = true) == [fe]
# errors
@test_throws ArgumentError export_nwb(e_x; file_name = fx)
@test_throws ArgumentError export_nwb(e_x; file_name = fx * ".h5")
@test_throws ArgumentError export_nwb(e_x; file_name = fx, overwrite = true, sex = "X")
@test_throws ArgumentError export_nwb(e_x; file_name = fx, overwrite = true, age = "35")
@test_throws ArgumentError export_nwb(import_snirf(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf"));
                                      file_name = fx, overwrite = true)
rm(fx)
rm(fe)
@test NeuroAnalyzer._nwb_ch_type("X1", "reference reference electrodes scalp") == "ref"

@info "Test: NWB writing"
fw = tempname() * ".nwb"
H.h5open(fw, "w") do f
    g = H.create_group(f, "groups")
    for nm in ("A", "B")
        NeuroAnalyzer._nwb_attrs!(H.create_group(g, nm), "ElectrodeGroup"; description = "shaft $nm")
    end
    df = DataFrame(:label => ["a1", "a2", "b1"], :x => [1.0, 2.0, 3.0], :bad => [false, true, false],
                   :group => ["/groups/A", "/groups/A", "/groups/B"],
                   :tags => [["x"], String[], ["y", "z"]])
    NeuroAnalyzer._nwb_table!(f, f, "electrodes", df; neurodata_type = "ElectrodesTable",
                              namespace = "core", description = "test electrodes",
                              references = ["group"])
    tr = DataFrame(:start_time => [0.0, 1.0], :electrode => [3, 1])
    NeuroAnalyzer._nwb_table!(f, f, "trials", tr; neurodata_type = "TimeIntervals",
                              namespace = "core", regions = Dict("electrode" => "/electrodes"))
    x = Float32.(reshape(1:4000, 4, 1000))
    ds = NeuroAnalyzer._nwb_dataset!(f, "data", x; compress = 4, shuffle = true)
    NeuroAnalyzer._nwb_attrs!(ds, "VectorData"; unit = "volts")
end
H.h5open(fw, "r") do f
    e = NeuroAnalyzer._nwb_table(f, f["electrodes"])
    @test e.id == [0, 1, 2] && e.label == ["a1", "a2", "b1"] && e.bad == [false, true, false]
    @test NeuroAnalyzer._nwb_isbool(f["electrodes/bad"])
    @test e.group == ["/groups/A", "/groups/A", "/groups/B"]
    @test e.tags == [["x"], String[], ["y", "z"]]
    @test NeuroAnalyzer._nwb_attr(f["electrodes"], "colnames") == ["label", "x", "bad", "group", "tags"]
    @test NeuroAnalyzer._nwb_type(f["electrodes"]) == ("ElectrodesTable", "core")
    @test occursin(r"^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$",
                   NeuroAnalyzer._nwb_attr(f["electrodes"], "object_id"))
    t = NeuroAnalyzer._nwb_table(f, f["trials"])
    @test t.electrode == [3, 1] && DataFrames.colmetadata(t, "electrode", "table") == "/electrodes"
    @test read(f["trials/electrode"]) == [2, 0]
    @test read(f["data"]) == Float32.(reshape(1:4000, 4, 1000))
    @test H.get_create_properties(f["data"]).filters ==
          [H.Filters.Shuffle(), H.Filters.Deflate(level = 4)]
end
rm(fw)
@test_throws ArgumentError NeuroAnalyzer._nwb_dataset!(nothing, "x", [1]; compress = 10)

true
