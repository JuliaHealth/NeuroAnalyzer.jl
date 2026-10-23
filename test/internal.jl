using NeuroAnalyzer
using Test
using DataFrames
using Cairo
using Random

@info "Initializing"

Random.seed!(123)

eeg = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
eeg.data = rand(-100:0.01:100, nchannels(eeg), sr(eeg) * 100, nepochs(eeg))
eeg.time_pts = eeg.time_pts[1:(sr(eeg) * 100)]
eeg.epoch_time = eeg.time_pts[1:(sr(eeg) * 100)]

locs = import_locs(joinpath(testfiles_path, "locs", "locs.ced"))

add_locs!(eeg; locs = locs)

e10 = epoch(eeg; ep_len = 10)
keep_epoch!(e10; ep = 1:10)

v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a = rand(2, 3, 4)
a1 = ones(2, 3, 4)
a2 = zeros(2, 3, 4)

##
## Channel Operations
##

@test NeuroAnalyzer._check_channels(a1, 1) === nothing
@test NeuroAnalyzer._check_channels(a1, 2) === nothing
@test NeuroAnalyzer._check_channels(e10, "F3") === nothing
@test NeuroAnalyzer._check_channels(e10, ["F3", "F4"]) === nothing
@test NeuroAnalyzer._check_channels(1:19, 2) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_channels(a1, 0)
@test_throws ArgumentError NeuroAnalyzer._check_channels(e10, "InvalidChannel")

@test NeuroAnalyzer._ch_units(e10, "Pz") == "μV"
@test_throws ArgumentError NeuroAnalyzer._ch_units(e10, "PZ")

@test NeuroAnalyzer._ch_rename("nirs_hbo") == "NIRS HbO concentration"

@test NeuroAnalyzer._set_channel_types(["fp1", "stim"]) == ["eeg", "mrk"]
@test NeuroAnalyzer._sort_channels(["emg", "eeg"]) == [2, 1]

@test NeuroAnalyzer._select_channels(e10, 1) == 1

l = labels(e10)
@test NeuroAnalyzer._ch_idx(e10, l[2]) == [2]
@test NeuroAnalyzer._ch_idx(e10, "") == Int64[]
@test NeuroAnalyzer._ch_idx(e10, "all") == collect(1:length(l))
@test NeuroAnalyzer._ch_idx(e10, "eeg") ==
      NeuroAnalyzer._ch_idx(e10, get_channel(e10; type = "eeg"))
# regex selects full labels matching anywhere
@test NeuroAnalyzer._ch_idx(e10, r"Fp") == findall(occursin.(r"Fp", l))
@test NeuroAnalyzer._ch_idx(e10, r"^Fp1$") == [findfirst(isequal("Fp1"), l)]
@test NeuroAnalyzer._ch_idx(e10, r"no_such_channel") == Int64[]
@test NeuroAnalyzer._ch_idx(l, r"Fp") == findall(occursin.(r"Fp", l))
@test NeuroAnalyzer._ch_idx(l, "all") == collect(1:length(l))
# keywords are expanded in place, every occurrence
e10_tmp = deepcopy(e10)
e10_tmp.header.recording[:bad_channel][[3, 5]] .= true
@test NeuroAnalyzer._ch_idx(e10_tmp, [l[1], "bad", l[2]]) == [1, 3, 5, 2]
@test NeuroAnalyzer._ch_idx(e10_tmp, ["bad", "bad"]) == [3, 5]
@test NeuroAnalyzer._ch_idx(e10, "bad") == Int64[]
@test_throws ArgumentError NeuroAnalyzer._ch_idx(e10, "no_such_channel")

##
## Epoch Operations
##

@test NeuroAnalyzer._check_epochs(e10, 1:10) === nothing
@test NeuroAnalyzer._check_segment(e10, 1, 10) === nothing
@test NeuroAnalyzer._check_segment(v1, 1, 2) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_epochs(e10, 0:10)
@test_throws ArgumentError NeuroAnalyzer._check_segment(e10, 1, 0)
@test_throws ArgumentError NeuroAnalyzer._check_segment(e10, 1, 0)

@test NeuroAnalyzer._len(e10, 0, 20) == 2560

a = NeuroAnalyzer._make_epochs(rand(10, 1000); ep_len = 100)
@test size(a) == (10, 100, 10)
a = NeuroAnalyzer._make_epochs(rand(10, 1000, 2); ep_len = 100)
@test size(a) == (10, 100, 20)
@test NeuroAnalyzer._get_epoch_markers(e10) == 0.0:10.0:90.0

@test NeuroAnalyzer._s2epoch(e10, 256, 512) == 1
@test NeuroAnalyzer._s2epoch(e10, 3256, 3512) == 2
@test NeuroAnalyzer._epoch2s(e10, 2) == (2561, 5120)

@test NeuroAnalyzer._select_epochs(e10, 1) == 1

@test NeuroAnalyzer._s2epoch(e10, 1, 256) == 1
@test NeuroAnalyzer._epoch2s(e10, 1) == (1, 2560)

##
## Data Operations
##

@test NeuroAnalyzer._check_var(:a, [:a], "a") === nothing
@test NeuroAnalyzer._check_var("a", ["a"], "a") === nothing
@test NeuroAnalyzer._check_markers(["aa", "bb"], "aa") === nothing

@test NeuroAnalyzer._check_svec("[1, 2]") == true
@test NeuroAnalyzer._check_svec("[1, 2]]") == true
@test NeuroAnalyzer._check_svec("1.2") == false
@test NeuroAnalyzer._check_srange("1:2") == true
@test NeuroAnalyzer._check_srange("1:2:3") == false
@test NeuroAnalyzer._check_srange("1.2") == false
@test NeuroAnalyzer._check_stuplei("(1, 2)") == true
@test NeuroAnalyzer._check_stuplei("(1.1, 2)") == false
@test NeuroAnalyzer._check_stuplef("(1, 2)") == true
@test NeuroAnalyzer._check_stuplef("(1.1, 2.9)") == true
@test NeuroAnalyzer._check_sfloat("2.9") == true
@test NeuroAnalyzer._check_sfloat("2") == true
@test NeuroAnalyzer._check_sint("2.9") == false
@test NeuroAnalyzer._check_sint("2") == true

@test NeuroAnalyzer._tuple_max((2, 1)) == (-2, 2)

@test NeuroAnalyzer._in(1, (1, 2.0)) == true
@test NeuroAnalyzer._in(0.9, (1, 2.0)) == false
@test NeuroAnalyzer._bin(1.1, (1, 2.0)) == true
@test NeuroAnalyzer._bin(1, (1, 2.0)) == false
@test NeuroAnalyzer._bin(0.9, (1.0, 2)) == false

@test NeuroAnalyzer._zeros([1.0, 2.0, -1.0, 0.0, -1.0, 2.0, 0.0, 1.0, -1.0]) == 7

s, x, y = NeuroAnalyzer._interpolate2d(rand(10), rand(10), rand(10))
@test size(s) == (100, 100)
@test length(x) == 100
@test length(y) == 100
@test_throws ArgumentError NeuroAnalyzer._interpolate2d(rand(10), rand(10), rand(11))

@test NeuroAnalyzer._tlength((1, 1)) == 1
@test_throws ArgumentError NeuroAnalyzer._tlength((1, 0))
@test NeuroAnalyzer._tlength((1, 10)) == 10
@test NeuroAnalyzer._tlength((1, 10.5)) == 10

@test NeuroAnalyzer._s2i("1, 2, 3") == [1, 2, 3]
@test NeuroAnalyzer._i2s([1, 2, 3]) == "1, 2, 3"
@test NeuroAnalyzer._s2tf("(1,2)") == (1.0, 2.0)
@test NeuroAnalyzer._s2ti("(1, 2)") == (1, 2)
@test NeuroAnalyzer._v2r([1, 2, 3]) == 1:3
@test NeuroAnalyzer._v2r(1:3) == 1:3
@test NeuroAnalyzer._v2r(1) == 1

##
## Create Operations
##

s = NeuroAnalyzer._create_subject(;
    id = "001",
    first_name = "A",
    middle_name = "B",
    last_name = "C",
    head_circumference = 64,
    handedness = "left",
    weight = 90,
    height = 180,
)
@test s isa Dict{Symbol, Any}

r = NeuroAnalyzer._create_recording_eeg(;
    data_type = "a",
    file_name = "a",
    file_size_mb = 1,
    file_type = "a",
    recording = "a",
    recording_date = "a",
    recording_time = "a",
    recording_notes = "a",
    channel_type = ["a"],
    reference = "a",
    clabels = ["a"],
    transducers = ["a"],
    units = ["a"],
    prefiltering = ["a"],
    sampling_rate = 1,
    gain = [0.0],
    channel_order = [1],
    line_frequency = 50,
    bad_channels = [false],
)
@test r isa Dict{Symbol, Any}

r = NeuroAnalyzer._create_recording_meg(;
    data_type = "a",
    file_name = "a",
    file_size_mb = 1,
    file_type = "a",
    recording = "a",
    recording_date = "a",
    recording_time = "a",
    recording_notes = "a",
    channel_type = ["a"],
    reference = "a",
    clabels = ["a"],
    units = ["a"],
    prefiltering = ["a"],
    sampling_rate = 1,
    magnetometers = [0],
    gradiometers = [0],
    coil_type = [""],
    channel_order = [1],
    line_frequency = 50,
    bad_channels = [false],
    ssp_labels = String[],
    ssp_data = Array{Float64}(undef, 0, 0),
    ssp_channels = Bool[],
)
@test r isa Dict{Symbol, Any}

e = NeuroAnalyzer._create_experiment(; name = "a", notes = "a", design = "a")
@test e isa Dict{Symbol, String}

hdr = NeuroAnalyzer._create_header(; subject = s, recording = r, experiment = e)
@test hdr isa NeuroAnalyzer.HEADER

##
## Markers Operations
##

@test NeuroAnalyzer._has_markers(["eeg", "mrk"]) == (true, 2)

df = NeuroAnalyzer._a2df(["1.0\x14stim\x14"])
@test names(df) == ["id", "start", "length", "value", "channel"]

##
## Locs Operations
##

@test NeuroAnalyzer._initialize_locs(e10) isa DataFrame
@test NeuroAnalyzer._initialize_locs() isa DataFrame

add_locs!(e10; locs = locs)
x, y, z, = e10.locs[!, :loc_x], e10.locs[!, :loc_y], e10.locs[!, :loc_z]
xn, yn = NeuroAnalyzer._locs_norm(x, y)
@test xn[1] ≈ -0.31 atol = 0.01
@test yn[1] ≈ 0.95 atol = 0.01
xn, yn, zn = NeuroAnalyzer._locs_norm(x, y, z)
@test xn[1] ≈ -0.31 atol = 0.01
@test yn[1] ≈ 0.95 atol = 0.01
@test zn[1] ≈ -0.03 atol = 0.01
locs[1, :loc_theta] = 108.12
locs[1, :loc_theta] = 108.12

locs = NeuroAnalyzer._locs_round(locs)
@test locs[1, :loc_theta] == 108.12

@test NeuroAnalyzer._angle_quadrant(45) == 1
@test NeuroAnalyzer._angle_quadrant(90) == 1
@test NeuroAnalyzer._angle_quadrant(90+45) == 2
@test NeuroAnalyzer._angle_quadrant(180) == 2
@test NeuroAnalyzer._angle_quadrant(180+45) == 3
@test NeuroAnalyzer._angle_quadrant(270) == 3
@test NeuroAnalyzer._angle_quadrant(270+45) == 4

@test NeuroAnalyzer._locs_remove_nans(DataFrame(:a=>[0.0, 1.0, NaN])) ==
      DataFrame(:a=>[0.0, 1.0, 0.0])

@test NeuroAnalyzer._flipx([1.0, 2.0, -1.0, 0.0, -1.0, 2.0, 0.0, 1.0, -1.0]) == [
    -0.6666666666666667,
    -1.6666666666666667,
    1.3333333333333333,
    0.3333333333333333,
    1.3333333333333333,
    -1.6666666666666667,
    0.3333333333333333,
    -0.6666666666666667,
    1.3333333333333333,
]

@test NeuroAnalyzer._find_bylabel(eeg.locs, "fp1") == 1

##
## Plotting Operations
##

@test NeuroAnalyzer._xlims(1:10) == (1.0, 10.0)
@test NeuroAnalyzer._ticks(1:10) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
@test NeuroAnalyzer._ticks((1, 10)) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
@test NeuroAnalyzer._erpticks(1:10) ==
      [1.0, 0.5, 0.0, 1.25, 2.5, 3.75, 5.0, 6.25, 7.5, 8.75, 10.0]
@test NeuroAnalyzer._erpticks((1, 10)) ==
      [1.0, 0.5, 0.0, 1.25, 2.5, 3.75, 5.0, 6.25, 7.5, 8.75, 10.0]
@test NeuroAnalyzer._set_defaults("a", "b", "c", "d", "e", "f") == ("a", "b", "c")
@test NeuroAnalyzer._set_defaults("default", "default", "default", "a", "b", "c") ==
      ("a", "b", "c")

@test NeuroAnalyzer._def_ylabel("eeg", "μV") == "Amplitude [μV]"

@test NeuroAnalyzer._convert_t(1.0, 2.0) == (1.0, "1.0 s", 2.0, "2.0 s")
@test NeuroAnalyzer._convert_t(0.22, 11.0) == (0.22, "220.0 ms", 11.0, "11.0 s")

##
## Miscellaneous Operations
##

r = NeuroAnalyzer._fir_response(rand(100), range(0; stop = π, length = 1024))
@test r isa Vector{ComplexF64}
@test length(r) == 1024

lmt = NeuroAnalyzer._labeled_matrix2dict(["a"], [[1.0]])
@test lmt == Dict("a"=>[1.0])
@test NeuroAnalyzer._dict2labeled_matrix(lmt) == (["a"], [[1.0]])
@test NeuroAnalyzer._clean_labels(["eeg fp1"]) == ["Fp1"]

df1, df2 = NeuroAnalyzer._split(DataFrame(:a=>1:10))
@test nrow(df1) == 8
@test nrow(df2) == 2

@test NeuroAnalyzer._copy_lt2ut([1 0 0; 1 1 0; 1 1 1]) == ones(3, 3)

t, et = NeuroAnalyzer._get_t(e10)
@test length(t) == 25600
@test length(et) == 2560
@test NeuroAnalyzer._get_t(1, 10, 10) == [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

@test NeuroAnalyzer._wl2ext(760) ≈ [645.5, 1669.0] .* 2.303
@test NeuroAnalyzer._wl2ext(761) ≈ [(645.5 + 650.0) / 2, (1669.0 + 1660.0) / 2] .* 2.303

# descriptions come verbatim from BioSig eventcodes.txt (whitespace runs collapsed)
@test NeuroAnalyzer._gdf_etp([0x01, 0x01]) == "artifact:EOG (blinks, fast, large amplitude)"
# event type is a little-endian UInt16: bytes [low, high]
@test NeuroAnalyzer._gdf_etp([0x04, 0x05]) == "ecg:R-point"
@test NeuroAnalyzer._gdf_etp([0x01, 0x03]) == "class1, Left hand - cue onset (BCI experiment)"
@test NeuroAnalyzer._gdf_etp([0x84, 0x00]) == "condition 131"
# bit 0x8000 marks the end of an event
@test NeuroAnalyzer._gdf_etp([0x03, 0x81]) == "artifact:EMG/Muscle (end)"
# end codes with their own (reserved) name in the table
@test NeuroAnalyzer._gdf_etp([0x02, 0x85]) == "ecg:P-wave-end, MDC_ECG_TIME_END_P"
@test NeuroAnalyzer._gdf_etp([0x20, 0x84]) == "Lights off"
@test NeuroAnalyzer._gdf_etp([0x80, 0x85]) == "end"
# reserved code, and the range comment line is not parsed as a code
@test NeuroAnalyzer._gdf_etp([0x0f, 0x01]) == "calibration"
@test NeuroAnalyzer._gdf_etp([0x1b, 0x03]) == "unknown code (0x031b)"
@test NeuroAnalyzer._gdf_etp([0x99, 0x09]) == "unknown code (0x0999)"
@test length(NeuroAnalyzer._GDF_ETP_TABLE) == 200

@test NeuroAnalyzer._midxy(1, 1, 4, 4) == (2.5, 2.5)

@test length(NeuroAnalyzer._split(1:55, wlen = 32, wstep = 8)) == 4
@test length(NeuroAnalyzer._fsplit(1:55, wlen = 32, wstep = 8)) == 3

@test NeuroAnalyzer._chunks(55, wlen = 32, wstep = 8) == [1 32; 9 40; 17 48; 25 55]
@test NeuroAnalyzer._chunks(1:55, wlen = 32, wstep = 8) == [1 32; 9 40; 17 48; 25 55]
@test NeuroAnalyzer._fchunks(55, wlen = 32, wstep = 8) == [1 32; 9 40; 17 48]
@test NeuroAnalyzer._fchunks(1:55, wlen = 32, wstep = 8) == [1 32; 9 40; 17 48]


##
## Check helpers
##

@test NeuroAnalyzer._chk2d(rand(2, 2)) === nothing
@test NeuroAnalyzer._chk3d(rand(2, 2, 2)) === nothing
@test NeuroAnalyzer._chk4d(rand(2, 2, 2, 2)) === nothing
@test_throws ArgumentError NeuroAnalyzer._chk2d(rand(2))
@test_throws ArgumentError NeuroAnalyzer._chk3d(rand(2, 2))
@test_throws ArgumentError NeuroAnalyzer._chk4d(rand(2, 2, 2))
@test NeuroAnalyzer._check_tuple((1, 2), (0, 10)) === nothing
@test NeuroAnalyzer._check_tuple((0, 10), (0, 10), "t", :bin) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_tuple((2, 1), (0, 10))
@test_throws ArgumentError NeuroAnalyzer._check_tuple((1, 1), (0, 10))
@test_throws ArgumentError NeuroAnalyzer._check_tuple((0, 11), (0, 10))
@test_throws ArgumentError NeuroAnalyzer._check_tuple((0, 5), (0, 10), "t", :bin)
@test_throws ArgumentError NeuroAnalyzer._check_tuple((1, 5), (0, 10), "t", :x)
@test NeuroAnalyzer._in(5, (1, 10), "x") === nothing
@test_throws ArgumentError NeuroAnalyzer._in(0, (1, 10), "x")
@test_throws ArgumentError NeuroAnalyzer._in(11, (1, 10), "x")
@test NeuroAnalyzer._bin(5, (1, 10), "x") === nothing
@test_throws ArgumentError NeuroAnalyzer._bin(1, (1, 10), "x")
@test_throws ArgumentError NeuroAnalyzer._bin(10, (1, 10), "x")
@test NeuroAnalyzer._check_datatype(e10, "eeg") === nothing
@test NeuroAnalyzer._check_datatype(e10, ["eeg", "erp"]) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_datatype(e10, "meg")
@test_throws ArgumentError NeuroAnalyzer._check_datatype(e10, ["meg", "nirs"])
@test NeuroAnalyzer._check_channels(e10, "Fp1", "eeg") === nothing
@test_throws ArgumentError NeuroAnalyzer._check_channels(e10, "Fp1", "mag")
@test NeuroAnalyzer._check_channels(["a", "b"], r"a") === nothing
@test_throws ArgumentError NeuroAnalyzer._check_channels(["a", "b"], r"z")
@test_throws ArgumentError NeuroAnalyzer._check_channels(String[], "a")
@test_throws ArgumentError NeuroAnalyzer._check_segment(v1, 0, 2)
@test_throws ArgumentError NeuroAnalyzer._check_segment(v1, 2, 1)
@test_throws ArgumentError NeuroAnalyzer._check_segment(v1, 1, 6)
@test NeuroAnalyzer._check_segment(e10, (0, 10)) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_segment(e10, (0, 1000))
@test_throws ArgumentError NeuroAnalyzer._check_var(:b, [:a], "a")
@test_throws ArgumentError NeuroAnalyzer._check_var("b", ["a"], "a")
@test_throws ArgumentError NeuroAnalyzer._check_markers(["aa", "bb"], "cc")
@test NeuroAnalyzer._check_sfloat(".5") == true
@test NeuroAnalyzer._check_sfloat(".") == false
@test NeuroAnalyzer._check_sfloat("") == false
@test NeuroAnalyzer._check_sint("") == false
@test NeuroAnalyzer._check_stuplei("(1,)") == false
@test NeuroAnalyzer._check_srange("2:1") == false

##
## Channel helpers
##

@test NeuroAnalyzer._ch_regex(["Fp1", "Fp2", "Cz"], r"Fp") == ["Fp1", "Fp2"]
@test NeuroAnalyzer._ch_units("csd") == "μV/m²"
@test NeuroAnalyzer._ch_units("mag") == "fT"
@test NeuroAnalyzer._ch_units("grad") == "fT/cm"
@test NeuroAnalyzer._ch_units("ecg") == "mV"
@test NeuroAnalyzer._ch_units("unknown_type") == ""
@test NeuroAnalyzer._ch_rename("eeg") == "EEG"
@test NeuroAnalyzer._ch_rename("unknown_type") == "unknown_type"
@test NeuroAnalyzer._def_ylabel("nirs_od", "") == "OD []"
@test NeuroAnalyzer._def_ylabel("nirs_hbo", "μM/mm") == "HbO concentration [μM/mm]"
@test NeuroAnalyzer._set_channel_types(["ecg", "eog", "emg", "a1", "acc x", "status", "760nm", "xyz", "Cz"]) ==
      ["ecg", "eog", "emg", "ref", "accel", "mrk", "nirs_od", "other", "eeg"]
# priorities: eeg (3) < ecg (6) < unknown (8) < mrk (9)
@test NeuroAnalyzer._sort_channels(["mrk", "eeg", "ecg", "unknown"]) == [2, 3, 4, 1]
@test NeuroAnalyzer._ch_idx(["a", "b", "c"], ["c", "a", "c"]) == [3, 1]
@test_throws ArgumentError NeuroAnalyzer._ch_idx(["a", "b"], "z")
@test NeuroAnalyzer._v2r([3, 1, 2]) == 1:3
@test NeuroAnalyzer._v2r([1, 3]) == [1, 3]
@test NeuroAnalyzer._v2r(Int64[]) == Int64[]

# labels
# format prefixes, duplicated prefixes, padding dots, reference suffixes, 10-20 capitalization
@test NeuroAnalyzer._clean_labels(["EEG Fp1", "EDF Cz", "EEG-O1", "EOG EOG", "ECG EKG", "Fc5.", "T9..", "EEG FP2-REF", "C3-LE", "PZ", "afz", "fpz"]) ==
      ["Fp1", "Cz", "O1", "EOG", "ECG", "FC5", "T9", "Fp2", "C3", "Pz", "AFz", "Fpz"]
# numbered channels: spacing only, leading zeros kept
@test NeuroAnalyzer._clean_labels(["MEG0113", "MEG 2443", "MEG  0111", "EEG001", "EEG 060", "EOG061", "STI 014", "STI101", "MEG000"]) ==
      ["MEG 0113", "MEG 2443", "MEG 0111", "EEG 001", "EEG 060", "EOG 061", "STI 014", "STI 101", "MEG 000"]
# other labels unchanged (apart from surrounding whitespace)
@test NeuroAnalyzer._clean_labels([" Fp1 ", "A12", "Pulse", "EDF Annotations", "EEG", "ECG1", "EOG 2"]) ==
      ["Fp1", "A12", "Pulse", "Annotations", "EEG", "ECG1", "EOG 2"]
# labels that would become identical keep their original form
@test NeuroAnalyzer._clean_labels(["EEG Fp1", "Fp1-REF", "Cz"]) == ["EEG Fp1", "Fp1-REF", "Cz"]
@test NeuroAnalyzer._clean_labels(["MEG0113", "MEG 0113"]) == ["MEG0113", "MEG 0113"]
@test NeuroAnalyzer._clean_labels(String[]) == String[]
# layout labels renamed to the recording labels
@test NeuroAnalyzer._label_key("MEG 0113") == NeuroAnalyzer._label_key("meg113") == "MEG113"
@test NeuroAnalyzer._match_labels(["MEG 113", "MEG 2443", "X"], ["MEG 0113", "MEG 2443"]) == ["MEG 0113", "MEG 2443", "X"]
@test NeuroAnalyzer._match_labels(["MEG 113"], ["MEG0113"]) == ["MEG0113"]

##
## Selection, mapping, conversion helpers
##

@test NeuroAnalyzer._select_indices(0, 5) == [1, 2, 3, 4, 5]
@test NeuroAnalyzer._select_indices(0, 5, 3) == [1, 2, 3]
@test NeuroAnalyzer._select_indices(0, 5, 9) == [1, 2, 3, 4, 5]
@test NeuroAnalyzer._select_indices(2, 5) == 2
@test NeuroAnalyzer._select_indices(2:4, 5) == [2, 3, 4]
sel = [3, 1]
@test NeuroAnalyzer._select_indices(sel, 5) == [1, 3]
# the caller's vector is not reordered
@test sel == [3, 1]
@test NeuroAnalyzer._select_epochs(e10, 0) == collect(1:10)
@test NeuroAnalyzer._select_epochs(e10, 0, 2) == [1, 2]
@test NeuroAnalyzer._select_channels(e10, 0, 3) == [1, 2, 3]
@test NeuroAnalyzer._map_channels([3, 5], [1, 3, 5]) == ([2, 3], [3, 5])
@test NeuroAnalyzer._map_channels(5, [1, 3, 5]) == (3, 5)
@test NeuroAnalyzer._map_channels(3:5, [3, 4, 5]) == ([1, 2, 3], [3, 4, 5])

@test NeuroAnalyzer._pl([1, 2]) == "s"
@test NeuroAnalyzer._pl([1]) == ""
@test NeuroAnalyzer._pl(2) == "s"
@test NeuroAnalyzer._pl(1) == ""
@test NeuroAnalyzer._get_range([1.2, 5.6]) == 4.0
@test NeuroAnalyzer._c(3) == [1, 2, 3]
@test NeuroAnalyzer._n2v(3) == [3]
@test NeuroAnalyzer._n2v([3, 4]) == [3, 4]
@test NeuroAnalyzer._n2v(1:3) == [1, 2, 3]
@test NeuroAnalyzer._v2s(["a", "b"]) == "ab"
@test NeuroAnalyzer._v2s([1, 2]) == "1, 2"
@test NeuroAnalyzer._swap(1, 2) == (2, 1)
@test NeuroAnalyzer._veqlen([1, 2], [1, 2, 3, 4]) == ([1, 2, 0, 0], [1, 2, 3, 4])
@test NeuroAnalyzer._veqlen([1, 2, 3], [1]) == ([1, 2, 3], [1, 0, 0])
@test NeuroAnalyzer._veqlen([1, 2], [3, 4]) == ([1, 2], [3, 4])
@test NeuroAnalyzer._vint2str([97, 98, 0, 32]) == "ab"
@test NeuroAnalyzer._fmem() > 0
@test NeuroAnalyzer._s2i("3") == 3
@test NeuroAnalyzer._s2i("[2:4]") == [2, 3, 4]
@test_throws ArgumentError NeuroAnalyzer._s2i("a")
@test NeuroAnalyzer._i2s(2:4) == "2, 3, 4"
# 3-D method: every slice is symmetrized in place
m3 = cat([1.0 0; 2 3], [4.0 0; 5 6]; dims = 3)
@test NeuroAnalyzer._copy_lt2ut(m3) == cat([1.0 2; 2 3], [4.0 5; 5 6]; dims = 3)
@test NeuroAnalyzer._common_ref_name("A1") == "common (A)"
@test NeuroAnalyzer._common_ref_name("m") == "common (M)"
@test NeuroAnalyzer._common_ref_name("avg") == "common"
@test NeuroAnalyzer._detect_montage(["Fp1-A1", "Fp2-A1"], ["eeg", "eeg"], "eeg") == "common (A)"
@test NeuroAnalyzer._detect_montage(["Fp1-F3", "Fp2-F4"], ["eeg", "eeg"], "eeg") == "bipolar"
@test NeuroAnalyzer._detect_montage(["Fp1M1", "Fp2M1"], ["eeg", "eeg"], "eeg") == "common (M)"
@test NeuroAnalyzer._detect_montage(["Fp1", "Fp2"], ["eeg", "eeg"], "eeg") == "physical"
@test NeuroAnalyzer._detect_montage(["Fp1", "Fp2"], ["eeg", "eeg"], "meg") == "physical"

# _fread(): little-endian binary fields
io = IOBuffer()
write(io, htol(Int16(-2)), htol(Float32(1.5)), UInt8[0x61, 0x62, 0x63], htol(UInt32(7)), htol(Float64(-0.25)))
seekstart(io)
@test NeuroAnalyzer._fread(io, 1, :i16) == -2
@test NeuroAnalyzer._fread(io, 1, :f32) == 1.5
@test NeuroAnalyzer._fread(io, 3, :s) == [97, 98, 99]
@test NeuroAnalyzer._fread(io, 1, :ui32) == 7
@test NeuroAnalyzer._fread(io, 1, :f64) == -0.25
@test_throws ArgumentError NeuroAnalyzer._fread(io, 1, :x)

##
## Time and epochs
##

@test NeuroAnalyzer._time_axis(4, 2) == [0.0, 0.5, 1.0, 1.5]
@test NeuroAnalyzer._time_axis(3, 10; t0 = -0.1) ≈ [-0.1, 0.0, 0.1]
# exact: no accumulated step error
@test NeuroAnalyzer._time_axis(1001, 1000)[end] == 1.0
@test NeuroAnalyzer._time_axis(0, 10) == Float64[]
@test_throws ArgumentError NeuroAnalyzer._time_axis(-1, 10)
@test_throws ArgumentError NeuroAnalyzer._time_axis(3, 0)
@test_throws ArgumentError NeuroAnalyzer._s2epoch(e10, 2, 1)
@test NeuroAnalyzer._s2epoch(e10, 0, 2560) == 1
@test NeuroAnalyzer._s2epoch(e10, 2560, 2561) == 1:2
@test NeuroAnalyzer._make_epochs(reshape(1.0:12, 1, 12); ep_len = 5)[1, :, :] == [1.0 6; 2 7; 3 8; 4 9; 5 10]
@test_throws ArgumentError NeuroAnalyzer._make_epochs(rand(2, 10); ep_len = 0)
@test_throws ArgumentError NeuroAnalyzer._make_epochs(rand(2, 10); ep_len = 11)

tps = NeuroAnalyzer._epochs_tps(e10)
@test size(tps) == (2, 10)
@test tps[:, 1] == [0.0, 2559 / 256]
@test tps[1, 2] == 10.0
e10_tmp = deepcopy(e10)
e10_tmp.markers = DataFrame(:id => ["1", "1"], :start => [5.0, 25.0], :length => [0.0, 0.0],
                            :value => ["a", "a"], :channel => [0, 0])
@test NeuroAnalyzer._markers_epochs(e10_tmp) == [1, 3]

# epochs around markers: signal 1:20 at 10 Hz, markers at samples 5, 12 and 18
s_tmp = reshape(collect(1.0:20), 1, 20, 1)
mrk = DataFrame(:id => ["1", "1", "2"], :start => [0.4, 1.1, 1.7], :length => [0.0, 0.0, 0.0],
                :value => ["a", "a", "b"], :channel => [0, 0, 0])
ep_tmp, mrk_tmp = NeuroAnalyzer._make_epochs_bymarkers(
    s_tmp; marker = "a", markers = mrk, marker_start = [5, 12], offset = 2, ep_len = 5, fs = 10,
)
@test ep_tmp[1, :, 1] == 3:7
@test ep_tmp[1, :, 2] == 10:14
# the markers move to the epoched timeline (sample 3 of epoch 1, sample 3 of epoch 2); the marker
# outside every epoch is dropped
@test mrk_tmp[!, :start] ≈ [0.2, 0.7]
@test mrk_tmp[!, :value] == ["a", "a"]
# an epoch reaching beyond the signal is skipped
ep_tmp, _ = NeuroAnalyzer._make_epochs_bymarkers(
    s_tmp; marker = "a", markers = mrk, marker_start = [1, 12], offset = 2, ep_len = 5, fs = 10,
)
@test size(ep_tmp, 3) == 1
@test_throws ArgumentError NeuroAnalyzer._make_epochs_bymarkers(
    s_tmp; marker = "a", markers = mrk, marker_start = [5], offset = -1, ep_len = 5, fs = 10,
)

##
## Markers
##

mrk = DataFrame(:id => ["1", "2", "3"], :start => [0.5, 1.5, 3.0], :length => [0.0, 0.0, 0.0],
                :value => ["a", "b", "c"], :channel => [0, 0, 0])
@test NeuroAnalyzer._delete_markers(mrk, (1, 2))[!, :value] == ["a", "c"]
@test NeuroAnalyzer._shift_markers(mrk, (1, 2))[!, :start] == [0.5, 1.5, 2.0]
# input is not modified
@test mrk[!, :start] == [0.5, 1.5, 3.0]
# remapping by sample number: +10 samples at 10 Hz = +1 s; nothing drops the marker
@test NeuroAnalyzer._remap_markers(mrk, 10, s -> s + 10)[!, :start] ≈ [1.5, 2.5, 4.0]
@test NeuroAnalyzer._remap_markers(mrk, 10, s -> s > 20 ? nothing : s)[!, :value] == ["a", "b"]
@test NeuroAnalyzer._has_markers(["eeg"]) == (false, 0)
@test NeuroAnalyzer._has_markers(e10) == false
@test NeuroAnalyzer._has_markers(e10_tmp) == true
# a single annotation record: start, duration, label
df = NeuroAnalyzer._a2df(["+0\x14\x14\0", "+1.5\x150.5\x14stim\x14"])
@test df[!, :start] == [1.5]
@test df[!, :length] == [0.5]
@test df[!, :value] == ["stim"]
# several records; records without a duration get length 0; ids are per label
df = NeuroAnalyzer._a2df(["+1\x14\x14stim\x14", "+2.5\x150.5\x14resp\x14", "+3\x14\x14stim\x14"])
@test df[!, :start] == [1.0, 2.5, 3.0]
@test df[!, :length] == [0.0, 0.5, 0.0]
@test df[!, :id] == ["1", "2", "1"]
@test nrow(NeuroAnalyzer._a2df(["+0\x14\x14\0"])) == 0

##
## Vectors, signals
##

# the trailing partial window starts where the next step would
@test NeuroAnalyzer._window_indices(10, 4, 3, false) == [1 4; 4 7; 7 10; 10 10]
@test NeuroAnalyzer._window_indices(11, 4, 3, false) == [1 4; 4 7; 7 10; 10 11]
@test NeuroAnalyzer._window_indices(11, 4, 3, true) == [1 4; 4 7; 7 10]
@test size(NeuroAnalyzer._window_indices(3, 4, 1, true)) == (0, 2)
@test_throws ArgumentError NeuroAnalyzer._window_indices(10, 4, 0, true)
@test NeuroAnalyzer._split(collect(1:10); wlen = 4, wstep = 3)[end] == [10.0]
@test NeuroAnalyzer._fsplit(collect(1:10); wlen = 4, wstep = 3)[2] == [4.0, 5, 6, 7]
@test NeuroAnalyzer._flipx([1.0, 3.0]) == [1.0, -1.0]
@test NeuroAnalyzer._zeros([1.0, 1.0, 1.0]) == 0
@test NeuroAnalyzer._group_delay(ones(5)) == 2
@test NeuroAnalyzer._group_delay(ones(1)) == 0
@test NeuroAnalyzer._remove_kernel(collect(1:10), ones(5)) == 3:8
@test NeuroAnalyzer._remove_kernel(collect(1:10), ones(1)) == 1:10
@test NeuroAnalyzer._reflect([1, 2, 3]) == [3, 2, 1, 1, 2, 3, 3, 2, 1]
@test NeuroAnalyzer._chop(NeuroAnalyzer._reflect([1, 2, 3])) == [1, 2, 3]
@test NeuroAnalyzer._reflect([1, 2], [3, 4], [5, 6]) == [4, 3, 1, 2, 6, 5]
@test NeuroAnalyzer._chop(collect(1:10), 2) == 3:8
@test NeuroAnalyzer._chop(collect(1:10), 0) == 1:10
@test_throws ArgumentError NeuroAnalyzer._reflect(Int64[])
@test_throws ArgumentError NeuroAnalyzer._chop(Int64[])
@test_throws ArgumentError NeuroAnalyzer._chop(collect(1:10), -1)
# gradient of a linear ramp is constant
g, gl = NeuroAnalyzer._gradient(2.0 .* (1:5))
@test all(x -> x[1] ≈ 2, g)
@test all(gl .≈ 2)
g, gl = NeuroAnalyzer._gradient(2.0 .* (1:5); rev = true)
@test all(x -> x[1] ≈ -2, g)
g, gl = NeuroAnalyzer._gradient([i + 2.0 * j for i = 1:3, j = 1:3])
@test all(gl .≈ sqrt(5))
# FIR response: H(0) = sum of the coefficients, |H(π)| of a 2-tap average = 0
@test NeuroAnalyzer._fir_response([0.5, 0.5], [0.0, π]) ≈ [1.0 + 0im, 0.0 + 0im] atol = 1e-12
@test NeuroAnalyzer._dict2labeled_matrix(Dict("a" => [1.0]); rev = false) == (["a"], [[1.0]])
@test_throws ArgumentError NeuroAnalyzer._dict2labeled_matrix(Dict())
@test_throws ArgumentError NeuroAnalyzer._dict2labeled_matrix(Dict("a" => [1]))
@test_throws ArgumentError NeuroAnalyzer._labeled_matrix2dict(["a", "b"], [[1.0]])
@test_throws ArgumentError NeuroAnalyzer._wl2ext(600)
@test NeuroAnalyzer._wl2ext(650) ≈ [506.0, 3743.0] .* 2.303
@test NeuroAnalyzer._gdf_etp_desc(0x0101) == "artifact:EOG (blinks, fast, large amplitude)"
@test isnothing(NeuroAnalyzer._gdf_etp_desc(0x0999))
@test_throws ArgumentError NeuroAnalyzer._gdf_etp([0x01])

# mesh scaling factors
msh = NeuroAnalyzer.GeometryBasics.Mesh(
    NeuroAnalyzer.GeometryBasics.Point3f[(1, -2, 0.5), (0, 0.5, -3), (-1.5, 1, 0)],
    [NeuroAnalyzer.GeometryBasics.TriangleFace(1, 2, 3)],
)
@test NeuroAnalyzer._mesh_normalize_xyz(msh) == 3.0
@test NeuroAnalyzer._mesh_normalize_xy(msh) == 2.0f0

# confusion matrix print-out
cm_out = mktemp() do path, io
    redirect_stdout(io) do
        NeuroAnalyzer._cm([1 2; 3 4])
    end
    flush(io)
    read(path, String)
end
@test occursin("prediction", cm_out)
@test occursin("│ 1 │ 2 │", replace(cm_out, r" +" => " "))

# timing helpers
@test NeuroAnalyzer._delay(0.01, 0.0) === nothing
@test (@elapsed NeuroAnalyzer._delay(0.05, 0.0)) >= 0.045
# a lag longer than the delay → no sleep
@test (@elapsed NeuroAnalyzer._delay(0.5, 1.0)) < 0.2
@test_throws ArgumentError NeuroAnalyzer._delay(0.01, -1)
@test NeuroAnalyzer._get_lag(2) >= 0
@test NeuroAnalyzer._mtime("1 + 1"; n = 2) >= 0
@test NeuroAnalyzer._nill() === nothing

##
## Record-based formats (EDF/BDF/GDF helpers)
##

@test NeuroAnalyzer._is_int_rate(256.0)
@test NeuroAnalyzer._is_int_rate(256.01)
@test !NeuroAnalyzer._is_int_rate(256.5)
# Fourier upsampling ×2 keeps the original samples
x = sin.(2π .* (0:63) ./ 16)
@test NeuroAnalyzer._resample_mirror(x, 128)[1:2:end] ≈ x atol = 1e-10
@test NeuroAnalyzer._resample_mirror(x, 64) == x
@test NeuroAnalyzer._resample_mirror([2.0], 3) == [2.0, 2.0, 2.0]
# nearest: sample-and-hold, every sample repeated k times
@test NeuroAnalyzer._resample_nearest([1, 2, 3, 4], 8) == [1, 1, 2, 2, 3, 3, 4, 4]
@test NeuroAnalyzer._resample_nearest([1, 2, 3, 4], 12) == [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
@test NeuroAnalyzer._resample_nearest(1:8, 4) == [1, 3, 5, 7]
# 12.5 Hz → exact ×2 upsampling to 25 Hz; integer rates are not changed
d_tmp = reshape(sin.(2π .* (0:63) ./ 16), 1, 64, 1)
d_new, fs_new = NeuroAnalyzer._integer_rate(d_tmp, 12.5)
@test fs_new == 25
@test size(d_new) == (1, 128, 1)
@test d_new[1, 1:2:end, 1] ≈ d_tmp[1, :, 1] atol = 1e-10
d_new, fs_new = NeuroAnalyzer._integer_rate(d_tmp, 256)
@test fs_new == 256 && d_new == d_tmp
# 99.7337 Hz × 11 = 1097.07 Hz: integer within the relative tolerance → ×11 upsampling
@test NeuroAnalyzer._integer_rate(rand(1, 10, 1), 99.7337)[2] == 1097
# no integer multiple k ≤ 16 → approximate resampling to the nearest integer rate
d_new, fs_new = NeuroAnalyzer._integer_rate(rand(1, 1000, 1), 99.0413)
@test fs_new == 99
@test size(d_new, 2) == round(Int64, 1000 * 99 / 99.0413)
@test_throws ArgumentError NeuroAnalyzer._integer_rate(d_tmp, 0)
@test NeuroAnalyzer._normalize_unit("uV\0") == "μV"
@test NeuroAnalyzer._normalize_unit("µV") == "μV"
@test NeuroAnalyzer._normalize_unit("MV") == "mV"
@test NeuroAnalyzer._normalize_unit(" degC ") == "degC"
u, sc = NeuroAnalyzer._resolve_units(["mV", "", "uV", "mV"], ["eeg", "eeg", "ecg", "mrk"])
@test u == ["μV", "μV", "μV", ""]
@test sc == [1000.0, 1.0, 1.0, 1.0]
@test_throws ArgumentError NeuroAnalyzer._resolve_units(["mV"], ["eeg", "eeg"])
x = [0.0, 50.0, 100.0]
@test NeuroAnalyzer._phys_scale!(x, 0, 100, -1, 1) == 0.02
@test x == [-1.0, 0.0, 1.0]
x = [5.0]
@test NeuroAnalyzer._phys_scale!(x, 3, 3, -1, 1) == 1.0
@test x == [5.0]
@test NeuroAnalyzer._record_count(-1, 1000, 100, "f") == 10
@test NeuroAnalyzer._record_count(5, 1000, 100, "f") == 5
@test NeuroAnalyzer._record_count(12, 1000, 100, "f") == 10
@test_throws ArgumentError NeuroAnalyzer._record_count(5, 1000, 0, "f")
@test NeuroAnalyzer._int24(0xff, 0xff, 0xff) == -1
@test NeuroAnalyzer._int24(0x01, 0x00, 0x00) == 1
@test NeuroAnalyzer._int24(0x00, 0x00, 0x80) == -8388608
@test NeuroAnalyzer._int24(0xff, 0xff, 0x7f) == 8388607
# 2 records, channel 1: 2 samples/record, channel 2: 1 sample/record (EDF, 16-bit)
raw = collect(reinterpret(UInt8, htol.(Int16[1, 2, 3, 4, 5, -6])))
sig, ann = NeuroAnalyzer._read_edf_records(raw, [2, 1], 2, 2)
@test sig == [[1.0, 2.0, 4.0, 5.0], [3.0, -6.0]]
@test isempty(ann)
sig, ann = NeuroAnalyzer._read_edf_records(raw, [2, 1], 2, 2; annotation_channels = [2])
@test sig[2] == Float64[]
@test length(ann) == 2
# BDF: 24-bit samples; the trigger channel keeps only bits 0-15
raw = UInt8[0xff, 0xff, 0xff, 0x01, 0x00, 0x00, 0x05, 0x00, 0xff]
sig, _ = NeuroAnalyzer._read_edf_records(raw, [2, 1], 1, 3; trigger_channels = [2])
@test sig == [[-1.0, 1.0], [5.0]]
@test_throws ArgumentError NeuroAnalyzer._read_edf_records(raw, [1], 1, 4)
d_tmp, rs = NeuroAnalyzer._merge_channels([[1.0, 2, 3, 4], [5.0, 6], Float64[]], [3]; nearest = [2])
@test rs
@test d_tmp[:, :, 1] == [1.0 2 3 4; 5 5 6 6; 0 0 0 0]
d_tmp, rs = NeuroAnalyzer._merge_channels([[1.0, 2], [3.0, 4]])
@test !rs
@test_throws ArgumentError NeuroAnalyzer._merge_channels([[1.0]], [1])

##
## FIFF low-level reader
##

bebytes(T, v) = collect(reinterpret(UInt8, hton.(T.(v))))
@test NeuroAnalyzer._i16i64(bebytes(Int16, [258])) == 258
@test NeuroAnalyzer._i32i64(bebytes(Int32, [-5])) == -5
@test NeuroAnalyzer._i32i32(bebytes(Int32, [7])) === Int32(7)
@test NeuroAnalyzer._ui32i32(bebytes(UInt32, [7])) === Int32(7)
@test NeuroAnalyzer._i8i8(UInt8[0xff]) === Int8(-1)
@test NeuroAnalyzer._f16f64(bebytes(Float16, [0.5])) == 0.5
@test NeuroAnalyzer._f32f64(bebytes(Float32, [1.5])) == 1.5
@test NeuroAnalyzer._f64f64(bebytes(Float64, [-2.25])) == -2.25
@test NeuroAnalyzer._i16f64(bebytes(Int16, [-3])) === -3.0
@test NeuroAnalyzer._i32f64(bebytes(Int32, [9])) === 9.0
@test NeuroAnalyzer._fiff_string(UInt8[0x61, 0x62, 0x00, 0x63]) == "ab"
@test NeuroAnalyzer._fiff_string(UInt8[0x61, 0x62]) == "ab"
# lookups: id ↔ name; unknown codes do not throw
@test NeuroAnalyzer._find_fiff_tag(200) == "nchan"
@test NeuroAnalyzer._find_fiff_tag("nchan") == 200
@test NeuroAnalyzer._find_fiff_tag(99999) == "unknown_99999"
@test_throws ArgumentError NeuroAnalyzer._find_fiff_tag("no_such_tag")
@test NeuroAnalyzer._find_fiff_block(101) == "meas_info"
@test NeuroAnalyzer._find_fiff_block("raw_data") == 102
@test NeuroAnalyzer._find_fiff_dt(4) == "float"
@test NeuroAnalyzer._find_fiff_dt(0x40000000 | 5) == "double"
@test NeuroAnalyzer._find_fiff_unit(107) == "V"
@test NeuroAnalyzer._find_fiff_mul(-6) == "μ"
@test NeuroAnalyzer._find_fiff_chtype(1) == "meg"
@test NeuroAnalyzer._find_fiff_chtype(2) == "eeg"
@test NeuroAnalyzer._find_fiff_gantry_type(2) == "free"
@test NeuroAnalyzer._find_fiff_dacq_system(4) == "triux"
@test NeuroAnalyzer._find_fiff_proj_item(10) == "eeg_avref"
@test NeuroAnalyzer._find_fiff_proj_by(1) == "space"
@test NeuroAnalyzer._find_fiff_coiltype(3012) == "vv_planar_t1"
@test NeuroAnalyzer._find_fiff_aspect(100) == "average"
@test NeuroAnalyzer._find_fiff_sss_job(10) == "sss_job_st"
@test NeuroAnalyzer._fiff_lookup(NeuroAnalyzer.fiff_units, :id, :unit, 112) == "T"
@test isnothing(NeuroAnalyzer._fiff_lookup(NeuroAnalyzer.fiff_units, :id, :unit, -99))
@test NeuroAnalyzer._fiff_name(NeuroAnalyzer.fiff_units, :id, :unit, -99) == "unknown_-99"
@test NeuroAnalyzer._fiff_id(NeuroAnalyzer.fiff_units, :id, :unit, "T") == 112
tbl = NeuroAnalyzer._fiff_table([1 => "a", 2 => "b"], :id, :name)
@test tbl[:id] == [1, 2] && tbl[:name] == ["a", "b"]
@test NeuroAnalyzer._fiff_eltype(4) == Float32
@test NeuroAnalyzer._fiff_eltype(16) == Int16
@test isnothing(NeuroAnalyzer._fiff_eltype(30))
@test NeuroAnalyzer._fiff_elements(3, bebytes(Int32, [1, -2])) == [1, -2]
@test NeuroAnalyzer._fiff_elements(20, bebytes(Float32, [1, 2, 3, 4])) == [1.0 + 2.0im, 3.0 + 4.0im]
@test isnothing(NeuroAnalyzer._fiff_elements(30, UInt8[]))
@test NeuroAnalyzer._fiff_scalar([5]) == 5
@test NeuroAnalyzer._fiff_scalar([1, 2]) == [1, 2]
# payload decoding by data type
@test NeuroAnalyzer._fiff_decode(3, bebytes(Int32, [42])) == 42
@test NeuroAnalyzer._fiff_decode(5, bebytes(Float64, [1.0, 2.0])) == [1.0, 2.0]
@test NeuroAnalyzer._fiff_decode(10, UInt8[0x4d, 0x45, 0x47]) == "MEG"
@test isnothing(NeuroAnalyzer._fiff_decode(0, UInt8[0x00]))
@test isnothing(NeuroAnalyzer._fiff_decode(3, UInt8[]))
# Julian day 2451545 = 2000-01-01
@test NeuroAnalyzer._fiff_decode(6, bebytes(Int32, [2451545])) == NeuroAnalyzer.Date(2000, 1, 1)
# old_pack: offset, scale, Int16 values
@test NeuroAnalyzer._fiff_decode(23, vcat(bebytes(Float32, [1.0, 2.0]), bebytes(Int16, [3, -1]))) == [7.0, -1.0]
# dig_point: kind, ident, r
@test NeuroAnalyzer._fiff_decode(33, vcat(bebytes(Int32, [1, 2]), bebytes(Float32, [0.5, -0.5, 0.25]))) ==
      (1, 2, 0.5, -0.5, 0.25)
# coord_trans: rotation (row-major), translation and their inverses
ct = NeuroAnalyzer._fiff_decode(
    35,
    vcat(bebytes(Int32, [1, 4]), bebytes(Float32, [0, -1, 0, 1, 0, 0, 0, 0, 1, 1, 2, 3, 0, 1, 0, -1, 0, 0, 0, 0, 1, -2, 1, -3])),
)
@test ct[1:2] == (1, 4)
@test ct[3] == [0.0 -1 0; 1 0 0; 0 0 1]
@test ct[4] == [1.0, 2, 3]
@test ct[3] * ct[5] ≈ [1.0 0 0; 0 1 0; 0 0 1]
# ch_info: 96 bytes
chi = NeuroAnalyzer._fiff_decode(
    30,
    vcat(bebytes(Int32, [1, 2, 1]), bebytes(Float32, [1.0, 2.0]), bebytes(Int32, [3012]), bebytes(Float32, collect(1:12)),
         bebytes(Int32, [201, -15]), UInt8.(collect("MEG 0113")), zeros(UInt8, 8)),
)
@test chi.kind == 1 && chi.coil_type == 3012 && chi.cal == 2.0
@test chi.loc == collect(1.0:12)
@test chi.ch_name == "MEG 0113"
@test isnothing(NeuroAnalyzer._fiff_decode(30, zeros(UInt8, 10)))
# dense matrix (C order, dims stored reversed): 2 × 3
dm = vcat(bebytes(Int32, 1:6), bebytes(Int32, [3, 2]), bebytes(Int32, [2]))
@test NeuroAnalyzer._fiff_matrix(NeuroAnalyzer._FIFFT_MATRIX | 3, dm) == [1 2 3; 4 5 6]
@test NeuroAnalyzer._fiff_decode(NeuroAnalyzer._FIFFT_MATRIX | 3, dm) == [1 2 3; 4 5 6]
@test_throws ArgumentError NeuroAnalyzer._fiff_matrix(NeuroAnalyzer._FIFFT_MATRIX | 3, vcat(bebytes(Int32, 1:2), bebytes(Int32, [3, 2]), bebytes(Int32, [2])))
# sparse CCS: [1 0 0; 0 0 2] = values, row indices, column pointers, (nnz, nrow, ncol), ndim
sm = vcat(bebytes(Float64, [1, 2]), bebytes(Int32, [0, 1]), bebytes(Int32, [0, 1, 1, 2]), bebytes(Int32, [2, 2, 3]), bebytes(Int32, [2]))
@test NeuroAnalyzer._fiff_matrix(NeuroAnalyzer._FIFFT_MATRIX | NeuroAnalyzer._FIFFT_SPARSE_CCS | 5, sm) ==
      [1.0 0 0; 0 0 2]
# sparse RCS: same matrix = values, column indices, row pointers
sm = vcat(bebytes(Float64, [1, 2]), bebytes(Int32, [0, 2]), bebytes(Int32, [0, 1, 2]), bebytes(Int32, [2, 2, 3]), bebytes(Int32, [2]))
@test NeuroAnalyzer._fiff_matrix(NeuroAnalyzer._FIFFT_MATRIX | NeuroAnalyzer._FIFFT_SPARSE_RCS | 5, sm) ==
      [1.0 0 0; 0 0 2]
# FIFF tag stream: block_start(meas_info), nchan = 3, block_end; the last tag has next = -1
fiff_tag(kind, dt, data, nxt = 0) = vcat(bebytes(Int32, [kind, dt, length(data), nxt]), data)
fio = IOBuffer(vcat(
    fiff_tag(104, 3, bebytes(Int32, [101])),
    fiff_tag(200, 3, bebytes(Int32, [3])),
    fiff_tag(105, 3, bebytes(Int32, [101]), -1),
))
k, dt_tmp, sz, d_tmp, nxt = NeuroAnalyzer._read_fiff_tag(fio)
@test (k, dt_tmp, sz, nxt) == (104, 3, 4, 0)
@test NeuroAnalyzer._i32i64(d_tmp) == 101
@test_throws EOFError NeuroAnalyzer._read_fiff_tag(IOBuffer(UInt8[0x00, 0x01]))
fd, fb, finst = NeuroAnalyzer._create_fiff_block(fio)
# columns: position, kind, data type, size, block level, block type
@test fb[:, 2] == [104, 200, 105]
@test fb[:, 1] == [0, 20, 40]
@test fb[:, 5] == [2, 2, 2]
@test fb[:, 6] == [101, 101, 101]
@test finst == [1, 1, 1]
@test NeuroAnalyzer._i32i64(fd[2]) == 3
bidx, btypes = NeuroAnalyzer._get_blocks(fb)
@test bidx == [[1, 2, 3]]
@test btypes == [101]
fobj = Any[(1, "nchan", "meas_info", 3), (2, "ch_info", "meas_info", :a), (3, "ch_info", "meas_info", :b),
           (4, "nchan", "raw_data", 9)]
pk = NeuroAnalyzer._pack_fiff_blocks(fobj, "meas_info", ["nchan", "ch_info", "sfreq"])
@test pk[:nchan] == 3
@test pk[:ch_info_1] == :a && pk[:ch_info_2] == :b
@test isnothing(pk[:sfreq])
tree_out = mktemp() do path, io
    redirect_stdout(io) do
        NeuroAnalyzer._fiff_tree(fobj)
        NeuroAnalyzer._view_fiff_block(fb)
    end
    flush(io)
    read(path, String)
end
@test occursin("meas_info - nchan", tree_out)
@test occursin("200 [nchan]", tree_out)

##
## Create header
##

common = (data_type = "x", file_name = "f", file_size_mb = 1, file_type = "t", recording = "r",
          recording_date = "d", recording_time = "12.30.00", recording_notes = "n", channel_type = ["a"],
          channel_order = [1], clabels = ["a"], units = ["u"], sampling_rate = 100, bad_channels = [false])
c_tmp = NeuroAnalyzer._common_recording_fields(; common...)
@test c_tmp[:label] == ["a"] && c_tmp[:sampling_rate] == 100 && c_tmp[:bad_channel] == [false]
for f in (NeuroAnalyzer._create_recording_seeg, NeuroAnalyzer._create_recording_ecog,
          NeuroAnalyzer._create_recording_eeg_like)
    r_tmp = f(; common..., reference = "a", transducers = ["t"], prefiltering = ["p"],
              line_frequency = 50, gain = [1.0])
    @test r_tmp[:reference] == "a" && r_tmp[:gain] == [1.0] && r_tmp[:epoch_id] == ""
end
for f in (NeuroAnalyzer._create_recording_sensors, NeuroAnalyzer._create_recording_eda,
          NeuroAnalyzer._create_recording_tpt, NeuroAnalyzer._create_recording_prefiltered)
    r_tmp = f(; common..., prefiltering = ["p"])
    @test r_tmp[:prefiltering] == ["p"] && r_tmp[:data_type] == "x"
end
r_tmp = NeuroAnalyzer._create_recording_nirs(; common..., wavelengths = [760.0, 850.0], wavelength_index = [1],
                                             optode_pairs = [1 1], src_labels = ["S1"], det_labels = ["D1"],
                                             opt_labels = ["S1", "D1"])
@test r_tmp[:wavelengths] == [760.0, 850.0] && r_tmp[:optode_labels] == ["S1", "D1"]
r_tmp = NeuroAnalyzer._create_recording_mep(; common..., stimulation_intensity = [50], coil_type = ["c"],
                                            stimulation_sample = [10], markers_pos = [1], markers_neg = [2])
@test r_tmp[:stimulation_intensity] == [50] && !haskey(r_tmp, :epoch_id)
# MEG: time separator normalized, geometry validated
r_tmp = NeuroAnalyzer._create_recording_meg(; common..., reference = "", prefiltering = [""], line_frequency = 50,
                                            magnetometers = [1], gradiometers = Int64[], coil_type = [""],
                                            ssp_labels = String[], ssp_channels = Bool[],
                                            ssp_data = zeros(0, 0), first_time = 1.5)
@test r_tmp[:recording_time] == "12:30:00"
@test r_tmp[:first_time] == 1.5 && r_tmp[:sss] == false
@test_throws ArgumentError NeuroAnalyzer._create_recording_meg(; common..., reference = "", prefiltering = [""],
    line_frequency = 50, magnetometers = [1], gradiometers = Int64[], coil_type = [""], ssp_labels = String[],
    ssp_channels = Bool[], ssp_data = zeros(0, 0), coil_loc = zeros(2, 12))
@test_throws ArgumentError NeuroAnalyzer._create_recording_meg(; common..., reference = "", prefiltering = [""],
    line_frequency = 50, magnetometers = [1], gradiometers = Int64[], coil_type = [""], ssp_labels = String[],
    ssp_channels = Bool[], ssp_data = zeros(0, 0), head_origin = [0.0])

##
## Locations
##

@test NeuroAnalyzer._sph_distance_cart(0, 0, 0, 3, 4, 0) == 5.0
# same point → 0; antipodal points on the unit sphere → 2
@test NeuroAnalyzer._sph_distance_sph(1, 30, 20, 1, 30, 20) ≈ 0 atol = 1e-6
@test NeuroAnalyzer._sph_distance_sph(1, 0, 0, 1, 180, 0) ≈ 2
@test NeuroAnalyzer._xyz2r(1, 2, 2) == 3.0
@test NeuroAnalyzer._wrap180([190, -190, 180, 0]) == [-170, 170, -180, 0]
xs, ys, zs = NeuroAnalyzer._sph2cart_raw([1.0], [90.0], [0.0])
@test [xs[1], ys[1], zs[1]] ≈ [0, 1, 0] atol = 1e-12
r_tmp, th_tmp, ph_tmp = NeuroAnalyzer._cart2sph_raw([0.0, 0.0], [1.0, 0.0], [0.0, 0.0])
@test r_tmp == [1.0, 0.0] && th_tmp ≈ [90.0, 0.0] && ph_tmp == [0.0, 0.0]
# round trip spherical → Cartesian → spherical through an identity transformation
lt = DataFrame(:loc_radius_sph => [1.0, 0.8], :loc_theta_sph => [45.0, -120.0], :loc_phi_sph => [10.0, -30.0])
r_tmp, th_tmp, ph_tmp = NeuroAnalyzer._sph_transform(lt, (x, y, z) -> (x, y, z))
@test r_tmp ≈ lt.loc_radius_sph && th_tmp ≈ lt.loc_theta_sph && ph_tmp ≈ lt.loc_phi_sph
# a 90° rotation about z shifts the horizontal angle
r_tmp, th_tmp, ph_tmp = NeuroAnalyzer._sph_transform(lt, (x, y, z) -> (-y, x, z))
@test th_tmp ≈ [135.0, -30.0] && ph_tmp ≈ lt.loc_phi_sph
@test NeuroAnalyzer._has_locs(e10) === nothing
e_tmp = deepcopy(e10)
e_tmp.locs = NeuroAnalyzer._initialize_locs()
@test_throws ArgumentError NeuroAnalyzer._has_locs(e_tmp)
NeuroAnalyzer._initialize_locs!(e_tmp)
@test nrow(e_tmp.locs) == length(get_channel(e10; type = ["meg", "grad", "mag", "eeg", "ecog", "seeg", "ieeg", "nirs_int", "nirs_od", "eog", "ref"]))
@test all(e_tmp.locs[!, :loc_x] .== 0)
@test NeuroAnalyzer._check_ch_locs([1, 2], ["Fp1", "Fp2"], ["Fp1", "Fp2"]) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_ch_locs(3, ["Fp1", "Fp2"], ["Fp1"])
@test_throws ArgumentError NeuroAnalyzer._check_ch_locs(2, ["Fp1", "Fp2"], ["Fp1"])
@test_throws ArgumentError NeuroAnalyzer._check_ch_locs(Int64[], ["Fp1"], ["Fp1"])
@test NeuroAnalyzer._check_locs_cols(e10.locs) === nothing
@test_throws ArgumentError NeuroAnalyzer._check_locs_cols(DataFrame(:label => ["a"]))
lbl1 = labels(e10)[1]
@test NeuroAnalyzer._loc_idx(e10, 1) == NeuroAnalyzer._find_bylabel(e10.locs, lbl1)
@test NeuroAnalyzer._loc_idx(e10, [lbl1]) == [NeuroAnalyzer._find_bylabel(e10.locs, lbl1)]
@test NeuroAnalyzer._idx2lab(e10, 1) == lbl1
@test NeuroAnalyzer._idx2lab(e10, [lbl1]) == [lbl1]
@test NeuroAnalyzer._ch_locs(e10, [1, 2])[!, :label] == labels(e10)[1:2]
@test NeuroAnalyzer._ch_locs(e10, labels(e10)[1:2]) == NeuroAnalyzer._ch_locs(e10, [1, 2])
@test NeuroAnalyzer._find_bylabel(e10.locs, "no_such_label") == Int64[]
@test NeuroAnalyzer._find_bylabel(e10.locs, ["FP2", "no_such_label"]) == [NeuroAnalyzer._find_bylabel(e10.locs, "Fp2")]
@test_throws ArgumentError NeuroAnalyzer._find_bylabel(DataFrame(:a => [1]), "a")
lr = DataFrame(:label => ["a"], :loc_radius => [1.234], :loc_theta => [1.236], :loc_x => [0.111],
               :loc_y => [0.0], :loc_z => [0.0], :loc_radius_sph => [1.0], :loc_theta_sph => [1.005],
               :loc_phi_sph => [0.0])
NeuroAnalyzer._locs_round!(lr)
@test lr[1, :loc_radius] == 1.23 && lr[1, :loc_theta] == 1.24 && lr[1, :loc_x] == 0.11
e_tmp = deepcopy(e10)
e_tmp.locs[1, :loc_x] = 0.123456
NeuroAnalyzer._locs_round!(e_tmp)
@test e_tmp.locs[1, :loc_x] == 0.12
@test NeuroAnalyzer._locs_round(e10)[1, :loc_y] == round(e10.locs[1, :loc_y]; digits = 2)
lnan = DataFrame(:a => [NaN, 1.0], :b => ["x", "y"])
NeuroAnalyzer._locs_remove_nans!(lnan)
@test lnan[!, :a] == [0.0, 1.0]
@test_throws ArgumentError NeuroAnalyzer._locs_remove_nans(DataFrame())
# min-max normalization of all axes jointly: the largest |coordinate| maps to ±1
e_tmp = NeuroAnalyzer._locs_norm(e10)
@test maximum(abs, Matrix(e_tmp.locs[!, [:loc_x, :loc_y, :loc_z]])) ≈ 1
@test e10.locs != e_tmp.locs
lnan = copy(e10.locs)
NeuroAnalyzer._locs_norm!(lnan)
@test lnan == NeuroAnalyzer._locs_norm(e10.locs)
# normalizing again changes nothing
NeuroAnalyzer._locs_norm!(e_tmp)
@test Matrix(e_tmp.locs[!, [:loc_x, :loc_y, :loc_z]]) ≈ Matrix(NeuroAnalyzer._locs_norm(e10).locs[!, [:loc_x, :loc_y, :loc_z]])

##
## Interpolation
##

# Shepard and nearest-neighbor interpolators reproduce the data at the data points
for im in (:sh, :nn)
    itp = NeuroAnalyzer._interpolator(im, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 2.0, 3.0])
    @test itp([1.0], [0.0])[1] ≈ 2.0
end
# natural neighbour interpolators reproduce the data at the data points and
# evaluate several points at once
lx = [0.0, 1.0, 0.0, 1.0, 0.5, 0.2, 0.8]
ly = [0.0, 0.0, 1.0, 1.0, 0.5, 0.7, 0.3]
lv = lx .+ 2 .* ly
for im in (:sib, :sib1, :lap, :far, :hiy, :tri)
    itp = NeuroAnalyzer._interpolator(im, lx, ly, lv)
    @test itp([1.0, 0.5], [0.0, 0.5]) ≈ [1.0, 1.5]
end
# spherical spline (almost) reproduces the data at the data points
itp = NeuroAnalyzer._interpolator(:sph, lx, ly, lv)
@test isapprox(itp(lx, ly), lv; atol = 1e-3)
# the vertex maps to the top of the sphere, r = 1 to the equator
@test NeuroAnalyzer._xy2sph([0.0, 1.0], [0.0, 0.0]) ≈ [0.0 0.0 1.0; 1.0 0.0 0.0] atol = 1e-12
@test_throws ArgumentError NeuroAnalyzer._interpolator(:x, [0.0], [0.0], [1.0])
s, x, y = NeuroAnalyzer._interpolate2d(rand(10), rand(10), rand(10), 50)
@test size(s) == (50, 50) && x[1] == -1.2 && x[end] == 1.2
# grid is not rounded, so dense grids have no duplicate coordinates
s, x, y = NeuroAnalyzer._interpolate2d(rand(10), rand(10), rand(10), 500, :sib)
@test size(s) == (500, 500) && allunique(x) && !any(isnan, s)
@test_throws ArgumentError NeuroAnalyzer._interpolate2d(rand(3), rand(3), rand(3), 0)

##
## Plot helpers
##

@test NeuroAnalyzer._ylims([-3.2, 5.7]) == (-6, 6)
@test NeuroAnalyzer._ylims([0.0, 4.2]) == (-4, 4)
@test NeuroAnalyzer._ylims(zeros(3)) == (-1.0, 1.0)
@test NeuroAnalyzer._ylims([-150.456, 20.0]) == (-151, 151)
@test NeuroAnalyzer._log_flim([0.0, 1.0, 2.0], (0, 2), :lin) == (0, 2)
@test NeuroAnalyzer._log_flim([0.0, 1.0, 2.0], (0, 2), :log) == (1.0, 2)
@test_throws ArgumentError NeuroAnalyzer._log_flim([-1.0, 0.0], (0, 2), :log)
@test NeuroAnalyzer._trace_colors(3; mono = false, black = true) == fill(:black, 3)
@test length(NeuroAnalyzer._trace_colors(4; mono = true, black = false)) == 4
@test NeuroAnalyzer._bernstein(1, 2; steps = 3) ≈ [0.0, 0.5, 0.0]
# quadratic Bézier: the curve passes through the end points, midpoint = (P0 + 2 P1 + P2) / 4
bx, by = NeuroAnalyzer._bernstein_poly([0, 1, 2], [0, 2, 0]; steps = 3)
@test bx ≈ [0.0, 1.0, 2.0] && by ≈ [0.0, 1.0, 0.0]
@test_throws ArgumentError NeuroAnalyzer._bernstein_poly([0, 1], [0])
# Catmull-Rom spline passes through every input point
sp_tmp = NeuroAnalyzer._smooth_path([0.0, 1.0, 2.0], [0.0, 1.0, 0.0]; n = 4)
@test length(sp_tmp) == 9
@test Tuple(sp_tmp[1]) == (0.0f0, 0.0f0) && Tuple(sp_tmp[5]) == (1.0f0, 1.0f0) && Tuple(sp_tmp[9]) == (2.0f0, 0.0f0)
dp = NeuroAnalyzer.DIPOLE([1.0 2.0 3.0], [1.0 2.0; 3.0 4.0; 5.0 6.0], [0.9], zeros(3), 1.0)
@test NeuroAnalyzer._dipole_at(dp, 2) == ([1.0, 2.0, 3.0], [2.0, 4.0, 6.0])
@test_throws ArgumentError NeuroAnalyzer._dipole_at(dp, 3)
# head outline (circle, nose, two ears) and fiducial labels
fig = NeuroAnalyzer.GLMakie.Figure()
ax = NeuroAnalyzer.GLMakie.Axis(fig[1, 1])
NeuroAnalyzer._draw_head_outline!(ax)
@test length(ax.scene.plots) == 4
# :xz is seen from behind: IN, LPA and RPA (no NAS)
NeuroAnalyzer._draw_head_labels!(ax; plane = :xz)
@test length(ax.scene.plots) == 7
@test_throws ArgumentError NeuroAnalyzer._draw_head_labels!(ax; plane = :ab)
@test NeuroAnalyzer._main_axis(fig) === ax
# legend outside the plot, with an inset above it
NeuroAnalyzer.GLMakie.lines!(ax, 1:3, 1:3; label = "a")
NeuroAnalyzer._side_legend!(fig, ax)
@test NeuroAnalyzer._side_layout(fig) isa NeuroAnalyzer.GLMakie.GridLayout
@test count(c -> c isa NeuroAnalyzer.GLMakie.Legend, fig.content) == 1
NeuroAnalyzer._figure_inset!(fig, NeuroAnalyzer.GLMakie.Figure(; size = (50, 40)); width = 25)
@test count(c -> c isa NeuroAnalyzer.GLMakie.Axis, fig.content) == 2
@test NeuroAnalyzer._fit_side_legend!(fig) === nothing
@test_throws ArgumentError NeuroAnalyzer._figure_inset!(fig, NeuroAnalyzer.GLMakie.Figure(); pos = :x)
fig = NeuroAnalyzer.GLMakie.Figure()
NeuroAnalyzer.GLMakie.Axis(fig[1, 1])
NeuroAnalyzer._ch_locs_inset!(fig, e10, [1, 2], [:red, :blue]; pos = :rt)
@test count(c -> c isa NeuroAnalyzer.GLMakie.Axis, fig.content) == 2

##
## Misc
##

@test NeuroAnalyzer._locs_norm_scale([1.0 -3.0; 2.0 0.5]) == 3.0
@test NeuroAnalyzer._locs_norm_scale(zeros(2, 3)) == 1.0
# all-zero locations stay at the origin
@test NeuroAnalyzer._locs_norm([0.0, 0.0], [0.0, 0.0]) == ([0.0, 0.0], [0.0, 0.0])
# scale only: the origin is kept (an asymmetric layout is not shifted)
xn, yn = NeuroAnalyzer._locs_norm([0.0, 1.0, 2.0], [0.0, 0.5, 0.0])
@test xn == [0.0, 0.5, 1.0] && yn == [0.0, 0.25, 0.0]
ax = NeuroAnalyzer.GLMakie.Axis(NeuroAnalyzer.GLMakie.Figure()[1, 1])
@test NeuroAnalyzer._style_axis!(ax) === ax
@test ax.titlesize[] == 18 && ax.xticklabelsize[] == 12

true

@info "Test: _plot_scale()"
rng_ps = NeuroAnalyzer.Random.MersenneTwister(1)
s_ps = 10 .* randn(rng_ps, 4, 10_000)
sc_ps = NeuroAnalyzer._plot_scale(s_ps)
@test isapprox(sc_ps, 60; rtol = 0.05)
# a large artifact and a flat channel do not change the scale
s_ps[1, 100:110] .= 1e4
@test isapprox(NeuroAnalyzer._plot_scale(vcat(s_ps, zeros(1, 10_000))), sc_ps; rtol = 0.01)
@test NeuroAnalyzer._plot_scale(zeros(3, 100)) == 1.0
