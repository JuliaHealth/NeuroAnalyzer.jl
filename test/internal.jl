using NeuroAnalyzer
using Test
using DataFrames
using Cairo

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
locs = import_locs(joinpath(testfiles_path, "locs.ced"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@test NeuroAnalyzer._check_channels(a1, 1) === nothing
@test NeuroAnalyzer._check_channels(a1, 2) === nothing
@test NeuroAnalyzer._check_channels(e10, "F3") === nothing
@test NeuroAnalyzer._check_channels(e10, ["F3", "F4"]) === nothing
@test NeuroAnalyzer._check_channels(1:19, 2) === nothing
@test NeuroAnalyzer._check_epochs(e10, 1:10) === nothing
@test NeuroAnalyzer._check_cidx(a1, 1) === nothing
@test NeuroAnalyzer._check_segment(e10, 1, 256) === nothing
@test NeuroAnalyzer._check_segment(v1, 1, 2) === nothing
@test NeuroAnalyzer._check_var(:a, [:a], "a") === nothing
@test NeuroAnalyzer._check_var("a", ["a"], "a") === nothing
@test NeuroAnalyzer._check_markers(["aa", "bb"], "aa") === nothing
@test NeuroAnalyzer._select_cidx(rand(2, 2), 1) == 1
s = NeuroAnalyzer._create_subject(id="001", first_name="A", middle_name="B", last_name="C", head_circumference=64, handedness="left", weight=90, height=180)
@test s isa Dict{Symbol, Any}
r = NeuroAnalyzer._create_recording_eeg(; data_type="a", file_name="a", file_size_mb=1, file_type="a", recording="a", recording_date="a", recording_time="a", recording_notes="a", channel_type=["a"], reference="a", clabels=["a"], transducers=["a"],units=["a"], prefiltering=["a"], sampling_rate=1, gain=[0.0], channel_order=[1], line_frequency=50, bad_channels=[false;;])
@test r isa Dict{Symbol, Any}
r = NeuroAnalyzer._create_recording_meg(; data_type="a", file_name="a", file_size_mb=1, file_type="a", recording="a", recording_date="a", recording_time="a", recording_notes="a", channel_type=["a"], reference="a", clabels=["a"], units=["a"], prefiltering=["a"], sampling_rate=1, magnetometers=[0], gradiometers=[0], coil_type=[""], channel_order=[1], line_frequency=50, bad_channels=[false;;], ssp_labels=String[], ssp_data=Array{Float64}(undef, 0, 0), ssp_channels=Bool[])
@test r isa Dict{Symbol, Any}
e = NeuroAnalyzer._create_experiment(name="a", notes="a", design="a")
@test e isa Dict{Symbol, String}
hdr = NeuroAnalyzer._create_header(s, r, e)
@test hdr isa NeuroAnalyzer.HEADER
r = NeuroAnalyzer._fir_response(rand(100), range(0, stop=π, length=1024))
@test r isa Vector{ComplexF32}
@test length(r) == 1024
s, x, y = NeuroAnalyzer._interpolate2d(rand(10), rand(10), rand(10))
@test size(s) == (100, 100)
@test length(x) == 100
@test length(y) == 100
@test NeuroAnalyzer._has_markers(["eeg", "mrk"]) == (true, 2)
@test !NeuroAnalyzer._has_markers(e10)
@test NeuroAnalyzer._set_channel_types(["fp1", "stim"]) == ["eeg", "mrk"]
df = NeuroAnalyzer._a2df(["1.0\x14stim\x14"])
@test names(df) == ["id", "start", "length", "value", "channel"]
@test NeuroAnalyzer._sort_channels(["emg", "eeg"]) == [2, 1]
lmt = NeuroAnalyzer._labeled_matrix2dict(["a"], [[1.0]])
@test lmt == Dict("a"=>[1.0])
@test NeuroAnalyzer._dict2labeled_matrix(lmt) == (["a"], [[1.0]])
@test NeuroAnalyzer._clean_labels(["eeg fp1"]) == ["fp1"]
add_component!(e10, c=:x, v=rand(4, 100))
l = NeuroAnalyzer._gen_clabels(e10, :x)
@test l == ["1", "2", "3", "4"]
@test NeuroAnalyzer._gen_clabels(a1) == ["1", "2"]
@test NeuroAnalyzer._len(e10, 0, 20) == 2560
x, y, z, = locs[!, :loc_x], locs[!, :loc_y], locs[!, :loc_z]
@test NeuroAnalyzer._initialize_locs(e10) isa DataFrame
NeuroAnalyzer._initialize_locs!(e10)
@test NeuroAnalyzer._initialize_locs() isa DataFrame
xn, yn = NeuroAnalyzer._locs_norm(x, y)
@test xn[1] ≈ -0.31
@test yn[1] == 0.95
xn, yn, zn = NeuroAnalyzer._locs_norm(x, y, z)
@test xn[1] ≈ -0.31
@test yn[1] == 0.95
@test zn[1] ≈ -0.03
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
@test NeuroAnalyzer._locs_remove_nans(DataFrame(:a=>[0.0, 1.0, NaN])) == DataFrame(:a=>[0.0, 1.0, 0.0])
a = NeuroAnalyzer._make_epochs(rand(10, 1000), ep_len=100)
@test size(a) == (10, 100, 10)
a = NeuroAnalyzer._make_epochs(rand(10, 1000), ep_n=100)
@test size(a) == (10, 10, 100)
a = NeuroAnalyzer._make_epochs(rand(10, 1000, 2), ep_len=100)
@test size(a) == (10, 100, 20)
a = NeuroAnalyzer._make_epochs(rand(10, 1000, 2), ep_n=100)
@test size(a) == (10, 20, 100)
# NeuroAnalyzer._make_epochs_bymarkers(signal::Array{<:Real, 3}; markers::DataFrame, marker_start::Vector{Int64}, epoch_offset::Int64, ep_len::Int64)
# NeuroAnalyzer._delete_markers(markers::DataFrame, segment::Tuple{Int64, Int64})
# NeuroAnalyzer._shift_markers(m::DataFrame, pos::Int64, offset::Int64)
@test NeuroAnalyzer._get_epoch_markers(e10) == 0.0:10.0:90.0
@test NeuroAnalyzer._tuple_max((2, 1)) == (-2, 2)
@test NeuroAnalyzer._s2v(1) == [1]
df1, df2 = NeuroAnalyzer._split(DataFrame(:a=>1:10))
@test nrow(df1) == 8
@test nrow(df2) == 2
@test NeuroAnalyzer._set_defaults("default", "default", "default", "a", "b", "c") == ("a", "b", "c")
@test NeuroAnalyzer._select_channels(e10, 1) == 1
@test NeuroAnalyzer._select_epochs(e10, 1) == 1
@test NeuroAnalyzer._convert_t(1.0, 2.0) == (1.0, "1.0 s", 2.0, "2.0 s")
@test NeuroAnalyzer._s2epoch(e10, 1, 256) == 1
@test NeuroAnalyzer._epoch2s(e10, 1) == (1, 2560)
@test NeuroAnalyzer._copy_lt2ut([1 0 0; 1 1 0; 1 1 1]) == ones(3, 3)
@test NeuroAnalyzer._tlength((1, 10)) == 10
t, et = NeuroAnalyzer._get_t(e10)
@test length(t) == 25600
@test length(et) == 2560
@test NeuroAnalyzer._get_t(1, 10, 10) == [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
@test NeuroAnalyzer._convert_t(0.22, 11.0) == (0.22, "220.0 ms", 11.0, "11.0 s")
@test NeuroAnalyzer._s2epoch(e10, 256, 512) == 1
@test NeuroAnalyzer._s2epoch(e10, 3256, 3512) == 2
@test NeuroAnalyzer._epoch2s(e10, 2) == (2561, 5120)
@test NeuroAnalyzer._ch_units(e10, "Pz") == "μV"
@test NeuroAnalyzer._ch_rename("nirs_hbo") == "NIRS HbO concentration"
@test NeuroAnalyzer._def_ylabel("eeg", "μV") == "Amplitude [μV]"
@test NeuroAnalyzer._wl2ext(760) == [1486.5865, 3843.707]
@test NeuroAnalyzer._gdf_etp([0x01, 0x01]) == "artifact:EOG (blinks)"
@test NeuroAnalyzer._check_svec("[1, 2]")
@test !NeuroAnalyzer._check_svec("[1, 2]]")
@test !NeuroAnalyzer._check_svec("1.2")
@test NeuroAnalyzer._check_srange("1:2")
@test !NeuroAnalyzer._check_srange("1:2:3")
@test !NeuroAnalyzer._check_srange("1.2")
@test NeuroAnalyzer._check_stuplei("(1, 2)")
@test !NeuroAnalyzer._check_stuplei("(1.1, 2)")
@test NeuroAnalyzer._check_stuplef("(1, 2)")
@test NeuroAnalyzer._check_stuplef("(1.1, 2.9)")
@test NeuroAnalyzer._check_sfloat("2.9")
@test NeuroAnalyzer._check_sfloat("2")
@test !NeuroAnalyzer._check_sint("2.9")
@test NeuroAnalyzer._check_sint("2")
@test NeuroAnalyzer._s2i("1, 2, 3") == [1, 2, 3]
@test NeuroAnalyzer._i2s([1, 2, 3]) == "1, 2, 3"
@test NeuroAnalyzer._s2tf("(1,2)") == (1.0, 2.0)
@test NeuroAnalyzer._s2ti("(1, 2)") == (1, 2)
@test NeuroAnalyzer._v2r([1, 2, 3]) == 1:3
@test NeuroAnalyzer._v2r(1:3) == 1:3
@test NeuroAnalyzer._v2r(1) == 1
@test NeuroAnalyzer._find_bylabel(eeg.locs, "fp1") == 1
@test NeuroAnalyzer._xlims(1:10) == (1.0, 10.0)
@test NeuroAnalyzer._ticks(1:10) == [1.0, 1.9, 2.8, 3.7, 4.6, 5.5, 6.4, 7.3, 8.2, 9.1, 10.0]
@test NeuroAnalyzer._ticks((1, 10)) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
@test NeuroAnalyzer._erpticks(1:10) == [1.0, 0.5, 0.0, 1.25, 2.5, 3.75, 5.0, 6.25, 7.5, 8.75, 10.0]
@test NeuroAnalyzer._erpticks((1, 10)) == [1.0, 0.5, 0.0, 1.25, 2.5, 3.75, 5.0, 6.25, 7.5, 8.75, 10.0]
@test NeuroAnalyzer._set_defaults("a", "b", "c", "d", "e", "f") == ("a", "b", "c")
@test NeuroAnalyzer._set_defaults("default", "default", "default", "d", "e", "f") == ("d", "e", "f")
@test NeuroAnalyzer._midxy(1, 1, 4, 4) == (2.5, 2.5)
@test NeuroAnalyzer._in(1, (1, 2.0)) == true
@test NeuroAnalyzer._bin(1, (1, 2.0)) == false
@test NeuroAnalyzer._bin(0.9, (1.0, 2)) == false
@test NeuroAnalyzer._zeros([1.0, 2.0, -1.0, 0.0, -1.0, 2.0, 0.0, 1.0, -1.0]) == 7
@test NeuroAnalyzer._flipx([1.0, 2.0, -1.0, 0.0, -1.0, 2.0, 0.0, 1.0, -1.0]) == [-0.3333333333333334, -1.3333333333333335, 1.6666666666666665, 0.6666666666666666, 1.6666666666666665, -1.3333333333333335, 0.6666666666666666, -0.3333333333333334, 1.6666666666666665]
@test length(NeuroAnalyzer._split(1:55, wlen=32, woverlap=8)) == 4
@test length(NeuroAnalyzer._fsplit(1:55, wlen=32, woverlap=8)) == 3
@test NeuroAnalyzer._chunks(55, wlen=32, woverlap=8) == [1 32; 9 40; 17 48; 25 55]
@test NeuroAnalyzer._chunks(1:55, wlen=32, woverlap=8) == [1 32; 9 40; 17 48; 25 55]
@test NeuroAnalyzer._fchunks(55, wlen=32, woverlap=8) == [1 32; 9 40; 17 48]
@test NeuroAnalyzer._fchunks(1:55, wlen=32, woverlap=8) == [1 32; 9 40; 17 48]

true