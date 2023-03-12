using DataFrames
using NeuroAnalyzer
using Test

@info "Initializing"
eeg = import_edf("files/eeg-test-edf.edf")
locs = locs_import("files/locs.ced")
e10 = epoch(eeg, ep_len=10*sr(eeg))
keep_epoch!(e10, epoch=1:10)
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@test _check_channels(a1, 1) === nothing
@test _check_channels(a1, 2) === nothing
@test _check_channels(e10, 1:24) === nothing
@test _check_channels(e10, 1:19, :eeg) === nothing
@test _check_channels(1:19, 2) === nothing
@test _check_epochs(e10, 1:10) === nothing
@test _check_cidx(a1, 1) === nothing
@test _check_segment(e10, 1, 256) === nothing
@test _check_segment(v1, 1, 2) === nothing
@test _check_var(:a, [:a], "a") === nothing
@test _check_var("a", ["a"], "a") === nothing
@test _check_markers(["aa", "bb"], "aa") === nothing
@test _get_ch_idx(["aa", "bb"], "aa") == 1
@test _get_ch_idx(["aa", "bb"], 1) == 1
@test _select_cidx(rand(2, 2), 1) == 1
s = _create_subject(id="a", first_name="a", middle_name="a", last_name="a", handedness="a", weight=0, height=0)
@test typeof(s) == Dict{Symbol, Any}
r = _create_recording_eeg(;data_type="a", file_name="a", file_size_mb=1, file_type="a", recording="a", recording_date="a", recording_time="a", recording_notes="a", channel_n=1, channel_type=["a"], reference="a", duration_samples=1, duration_seconds=0.0, epoch_n=1, epoch_duration_samples=1, epoch_duration_seconds=0.0, clabels=["a"], transducers=["a"],units=["a"], prefiltering=["a"], sampling_rate=1, gain=[0.0])
@test typeof(r) == Dict{Symbol, Any}
r = _create_recording_meg(;data_type="a", file_name="a", file_size_mb=1, file_type="a", recording="a", recording_date="a", recording_time="a", recording_notes="a", channel_n=1, channel_type=["a"], reference="a", duration_samples=1, duration_seconds=0.0, epoch_n=1, epoch_duration_samples=1, epoch_duration_seconds=0.0, clabels=["a"], units=["a"], prefiltering=["a"], sampling_rate=1, magnetometers=[0], gradiometers=[0], gradiometers_planar=[0], gradiometers_axial=[0], coils=[0])
@test typeof(r) == Dict{Symbol, Any}
e = _create_experiment(experiment_name="a", experiment_notes="a", experiment_design="a")
@test typeof(e) == Dict{Symbol, String}
hdr = _create_header(s, r, e; component_names=Symbol[], has_markers=false, has_locs=false, history=String[])
@test typeof(hdr) == NeuroAnalyzer.HEADER
r = _fir_response(rand(100), range(0, stop=π, length=1024))
@test typeof(r) == Vector{ComplexF32}
@test length(r) == 1024
s, x, y = _interpolate(rand(10), rand(10), rand(10))
@test size(s) == (100, 100)
@test length(x) == 100
@test length(y) == 100
@test _has_markers(["eeg", "mrk"]) == (true, 2)
@test _set_channel_types(["fp1", "stim"]) == ["eeg", "mrk"]
df = _m2df(["\x141.0\x14stim"])
@test names(df) == ["id", "start", "length", "description", "channel"]
@test _sort_channels(["emg", "eeg"]) == [2, 1]
lm = _labeled_matrix2dict(["a"], [[1.0]])
@test lm == Dict("a"=>[1.0])
@test _dict2labeled_matrix(lm) == (["a"], [[1.0]])
@test _clean_labels(["eeg fp1"]) == ["fp1"]
add_component!(e10, c=:x, v=rand(4, 100))
l = _gen_clabels(e10, :x)
@test l == ["1", "2", "3", "4"]
@test _gen_clabels(a1) == ["1", "2"]
@test _channel2channel_name(1:10) == "1:10"
@test _len(e10, 0, 20) == 2560
x, y, z, = locs[!, :loc_x], locs[!, :loc_y], locs[!, :loc_z]
_locnorm(x, y)
_locnorm(x, y, z)
locs[1, :loc_theta] = 108.1239
locs[1, :loc_theta] = 108.1239
locs = _round_locs(locs)
@test locs[1, :loc_theta] == 108.124
@test _angle_quadrant(45) == 1
@test _angle_quadrant(90) == 1
@test _angle_quadrant(90+45) == 2
@test _angle_quadrant(180) == 2
@test _angle_quadrant(180+45) == 3
@test _angle_quadrant(270) == 3
@test _angle_quadrant(270+45) == 4
a = _make_epochs(rand(10, 1000), ep_len=100)
@test size(a) == (10, 100, 10)
a = _make_epochs(rand(10, 1000), ep_n=100)
@test size(a) == (10, 10, 100)
a = _make_epochs(rand(10, 1000, 2), ep_len=100)
@test size(a) == (10, 100, 20)
a = _make_epochs(rand(10, 1000, 2), ep_n=100)
@test size(a) == (10, 20, 100)
# _make_epochs_bymarkers(signal::Array{<:Real, 3}; markers::DataFrame, marker_start::Vector{Int64}, epoch_offset::Int64, ep_len::Int64)
# _map_channels(channel::Union{Int64, Vector{Int64}, AbstractRange}, channels=Vector{Int64})
# _delete_markers(markers::DataFrame, segment::Tuple{Int64, Int64})
# _shift_markers(m::DataFrame, pos::Int64, offset::Int64)
@test _get_epoch_markers(e10) == 0.0:10.0:90.0
@test _tuple_max((2, 1)) == (-2, 2)
@test _s2v(1) == [1]
df = DataFrame(:a=>1:10)
df1, df2 = _split(df)
@test nrow(df1) == 8
@test nrow(df2) == 2
@test _set_defaults("default", "default", "default", "a", "b", "c") == ("a", "b", "c")
@test _select_channels(e10, 1) == 1
@test _select_epochs(e10, 1) == 1
@test _convert_t(1.0, 2.0) == (1.0, "1.0 s", 2.0, "2.0 s")
@test _s2epoch(e10, 1, 256) == 1
@test _epoch2s(e10, 1) == (1, 2560)

# these function are still in work:
## FIFF
# _read_fiff_tag(fid::IOStream)
# _get_fiff_block_type(fid::IOStream, tags::Vector{NTuple{5, Int64}}, id)
# _read_fiff_tag(fid::IOStream, fiff_blocks::Matrix{Int64}, tag_id::Int64)
# _read_fiff_data(fid::IOStream, fiff_blocks::Matrix{Int64}, id::Union{Int64, Nothing})
# _find_fiff_tag(fiff_blocks::Matrix{Int64}, id::Int64)
# _extract_struct(s, id::Int64)
# _fiff_channel_type(channel_types::Vector{Int32})
# _create_fiff_block(file_name="a")
# _view_fiff_block(fiff_block::Matrix{Int64})
## tester
# _get_lag(n::Int64=10)
# _mtime(f::String; n::Int64=10)
# _delay(n::Real, l::Real)
# _check_accuracy()

true