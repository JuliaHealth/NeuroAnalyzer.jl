using Test
using CSV
using DataFrames

# BIDS test datasets:
# - ds002718 (OpenNeuro, CC0): EEG, EEGLAB .set, subject 002; electrodes.tsv, events.tsv,
#   sidecars in the subject folder and the dataset root
# - ds003029 (OpenNeuro, CC0): iEEG, BrainVision, subjects jh101 (run 01 cropped to 90 s) and
#   jh102 (run 01 cropped to 120 s); participants.tsv
# - ds004504 (OpenNeuro, CC0): EEG, EEGLAB .set, subjects 001, 002 (group A), 037, 038 (C),
#   066, 067 (F); participants.tsv and participants.json
# - fnirs_tapping (R. Luke, BIDS-NIRS-Tapping, Zenodo 6575155): NIRS, SNIRF, subjects 01, 02
bids_root = joinpath(testfiles_path, "bids")
d1 = joinpath(bids_root, "ds002718")
d2 = joinpath(bids_root, "ds003029")
d4 = joinpath(bids_root, "ds004504")
f1 = joinpath(d1, "sub-002", "eeg", "sub-002_task-FaceRecognition_eeg.set")

@info "Test: BIDS entities and sidecars"
e = NeuroAnalyzer._bids_entities("sub-01_ses-2_task-x_acq-y_run-3_eeg.vhdr")
@test e == Dict("sub" => "01", "ses" => "2", "task" => "x", "acq" => "y", "run" => "3",
                "suffix" => "eeg")
@test NeuroAnalyzer._bids_entities("task-x_events.json") == Dict("task" => "x",
                                                                  "suffix" => "events")
@test NeuroAnalyzer._bids_root(f1) == abspath(d1)
# inheritance: coordsystem.json in the subject folder, events.json in the dataset root
@test NeuroAnalyzer._bids_sidecar(f1, "coordsystem", ".json") ==
      joinpath(abspath(d1), "sub-002", "sub-002_task-FaceRecognition_coordsystem.json")
@test NeuroAnalyzer._bids_sidecars(f1, "events", ".json") ==
      [joinpath(abspath(d1), "task-FaceRecognition_events.json")]
@test isnothing(NeuroAnalyzer._bids_sidecar(f1, "optodes", ".tsv"))
@test NeuroAnalyzer._bids_json(f1, "eeg")["PowerLineFrequency"] == 50

@info "Test: bids_files()"
bf = bids_files(d1)
@test bf.subject == ["002"] && all(bf.task .== "FaceRecognition")
@test all(bf.modality .== "eeg") && all(bf.session .== "")
@test bf.file[1] == joinpath("sub-002", "eeg", "sub-002_task-FaceRecognition_eeg.set")
@test nrow(bids_files(d1; subject = "002")) == 1 && nrow(bids_files(d1; subject = "003")) == 0
@test nrow(bids_files(d1; task = "x")) == 0
bf2 = bids_files(d2)
@test bf2.subject == ["jh101", "jh102"] && all(bf2.session .== "presurgery")
@test all(bf2.run .== "01") && all(bf2.acq .== "ecog") && all(bf2.modality .== "ieeg")
# participants.tsv; "F" is not read as a Boolean
@test bf2.site == ["JHH", "JHH"] && bf2.outcome == ["F", "NR"] && bf2.engel_score == [4.0, -1.0]
bf4 = bids_files(d4)
@test bf4.subject == ["001", "002", "037", "038", "066", "067"]
@test bf4.Group == ["A", "A", "C", "C", "F", "F"] && bf4.Gender[1] == "F"
@test bf4.Age[1] == 57 && bf4.MMSE[1] == 16 && all(bf4.task .== "eyesclosed")
@test nrow(bids_files(d4; subject = ["001", "066"])) == 2
@test_throws ArgumentError bids_files(joinpath(bids_root, "none"))

@info "Test: import_bids() (EEG, EEGLAB)"
o = import_bids(f1)
rec = o.header.recording
@test o.header.subject[:id] == "002" && startswith(o.header.notes, "Task: FaceRecognition")
@test rec[:line_frequency] == 50 && rec[:reference] == "nose"
# channels.tsv: 70 EEG, HEOG + VEOG, 2 EKG
@test count(==("eeg"), rec[:channel_type]) == 70
@test count(==("eog"), rec[:channel_type]) == 2 && count(==("ecg"), rec[:channel_type]) == 2
# events.tsv
ev = CSV.read(joinpath(d1, "sub-002", "eeg", "sub-002_task-FaceRecognition_events.tsv"),
              DataFrame; missingstring = "n/a")
@test nrow(o.markers) == nrow(ev) && o.markers.start ≈ sort(ev.onset)
@test o.markers.value[1] == "unfamiliar_new" && "face_type" in names(o.markers)
# electrodes.tsv (EEGLAB coordinates, x: nose) → x: right, y: nose; units of the file (cm)
el = CSV.read(joinpath(d1, "sub-002", "eeg", "sub-002_task-FaceRecognition_electrodes.tsv"),
              DataFrame; missingstring = "n/a")
i = findfirst(==("EEG 001"), o.locs.label)
@test o.locs.loc_x[i] ≈ -el.y[1] && o.locs.loc_y[i] ≈ el.x[1] && o.locs.loc_z[i] ≈ el.z[1]
@test o.locs.label == Base.filter(in(o.locs.label), labels(o))
# coordinate system and unit of coordsystem.json
@test locs_coord(o) == (system = "EEGLAB", unit = "cm", guessed = false)
# without sidecars
o0 = import_bids(f1; channels = false, events = false, electrodes = false)
@test count(==("eog"), o0.header.recording[:channel_type]) == 0
@test o0.markers != o.markers

@info "Test: import_bids() (iEEG, BrainVision)"
f2 = joinpath(d2, bf2.file[1])
o2 = import_bids(f2)
rec2 = o2.header.recording
@test datatype(o2) == "ecog" && sr(o2) == 1000 && signal_len(o2) == 90_000
@test o2.header.subject[:id] == "jh101" && startswith(o2.header.notes, "Task: ictal")
@test rec2[:line_frequency] == 60
ct2 = CSV.read(replace(f2, "_ieeg.vhdr" => "_channels.tsv"), DataFrame; missingstring = "n/a")
@test count(rec2[:bad_channel]) == count(==("bad"), ct2.status)
@test count(==("ecog"), rec2[:channel_type]) == count(==("ECOG"), ct2.type)
@test all(rec2[:lowpass] .== 500) && all(rec2[:highpass] .== 0)
ev2 = CSV.read(replace(f2, "_ieeg.vhdr" => "_events.tsv"), DataFrame; missingstring = "n/a")
@test nrow(o2.markers) == nrow(ev2) && o2.markers.start[end] ≈ 82.8
@test_throws ArgumentError import_bids(joinpath(d1, "dataset_description.json"))
@test_throws ArgumentError import_bids(joinpath(d1, "none_eeg.set"))

@info "Test: import_bids() (NIRS, optodes.tsv)"
# BIDS-NIRS dataset from the SNIRF test file
nd = mktempdir()
write(joinpath(nd, "dataset_description.json"), "{\"Name\": \"nirs\", \"BIDSVersion\": \"1.10.0\"}")
mkpath(joinpath(nd, "sub-01", "nirs"))
fn = joinpath(nd, "sub-01", "nirs", "sub-01_task-x_nirs.snirf")
cp(joinpath(testfiles_path, "nirs", "fnirs-test-snirf.snirf"), fn)
n0 = import_recording(fn)
src = n0.header.recording[:src_labels]
det = n0.header.recording[:det_labels]
ol = n0.header.recording[:optode_labels]
no = length(ol)
function write_optodes(names, x, y, z; template = false)
    open(joinpath(nd, "sub-01", "nirs", "sub-01_task-x_optodes.tsv"), "w") do io
        println(io, template ? "name\ttype\tx\ty\tz\ttemplate_x\ttemplate_y\ttemplate_z" :
                               "name\ttype\tx\ty\tz")
        for i in eachindex(names)
            t = i <= length(src) ? "source" : "detector"
            println(io, template ? "$(names[i])\t$t\tn/a\tn/a\tn/a\t$(x[i])\t$(y[i])\t$(z[i])" :
                                   "$(names[i])\t$t\t$(x[i])\t$(y[i])\t$(z[i])")
        end
    end
end
write_coordsystem(sys) = write(joinpath(nd, "sub-01", "nirs", "sub-01_task-x_coordsystem.json"),
                               "{\"NIRSCoordinateSystem\": \"$sys\", \"NIRSCoordinateUnits\": \"cm\"}")
x, y, z = collect(1.0:no), collect(101.0:(100 + no)), fill(5.0, no)
# names of the recording, RAS system
write_optodes(ol, x, y, z)
write_coordsystem("CapTrak")
@test bids_files(nd).modality == ["nirs"]
n1 = import_bids(fn)
@test n1.locs.label == ol && n1.locs.loc_x == x && n1.locs.loc_y == y && n1.locs.loc_z == z
@test n1.header.recording[:length_unit] == "cm"
@test locs_coord(n1) == (system = "CapTrak", unit = "cm", guessed = false)
# ALS system (x: nose, y: left): rotated
write_coordsystem("EEGLAB")
n2 = import_bids(fn)
@test n2.locs.loc_x == -y && n2.locs.loc_y == x
# other names: matched by the order of the sources and detectors; template positions
write_optodes(vcat("S" .* string.(1:length(src)), "D" .* string.(1:length(det))), x, y, z;
              template = true)
write_coordsystem("CapTrak")
n3 = import_bids(fn)
@test n3.locs.label == ol && n3.locs.loc_x == x && n3.locs.loc_z == z
# not applied
@test import_bids(fn; optodes = false).locs == n0.locs

@info "Test: import_bids() (NIRS, fnirs_tapping)"
d3 = joinpath(bids_root, "fnirs_tapping")
bf3 = bids_files(d3)
@test bf3.subject == ["01", "02"] && all(bf3.modality .== "nirs") && all(bf3.task .== "tapping")
@test bf3.age == [34, 32] && bf3.sex == ["M", "F"]
f3 = joinpath(d3, bf3.file[1])
o3 = import_bids(f3)
@test datatype(o3) == "nirs" && o3.header.subject[:id] == "01"
@test startswith(o3.header.notes, "Task: tapping")
# optodes.tsv: CapTrak (RAS), meters
ot3 = CSV.read(joinpath(d3, "sub-01", "nirs", "sub-01_optodes.tsv"), DataFrame)
@test o3.header.recording[:length_unit] == "m"
@test locs_coord(o3) == (system = "CapTrak", unit = "m", guessed = false)
@test o3.locs.label == o3.header.recording[:optode_labels]
i = indexin(o3.locs.label, ot3.name)
@test o3.locs.loc_x ≈ ot3.x[i] && o3.locs.loc_y ≈ ot3.y[i] && o3.locs.loc_z ≈ ot3.z[i]
ev3 = CSV.read(joinpath(d3, "sub-01", "nirs", "sub-01_task-tapping_events.tsv"), DataFrame)
@test nrow(o3.markers) == nrow(ev3) && o3.markers.start ≈ ev3.onset
@test Set(o3.markers.value) ⊇ Set(["Control", "Tapping/Left", "Tapping/Right"])
st3 = study_bids(d3; out_dir = mktempdir(), group = :sex, covariates = [:age])
@test st3.design.group == ["M", "F"] && st3.design.age == [34, 32]

@info "Test: import_bids() (EEG, ds004504)"
f4 = joinpath(d4, bf4.file[1])
o4 = import_bids(f4)
@test o4.header.subject[:id] == "001" && startswith(o4.header.notes, "Task: eyesclosed")
@test o4.header.recording[:line_frequency] == 50 && o4.header.recording[:reference] == "A1 A2"
@test nchannels(o4) == 19 && sr(o4) == 500

@info "Test: study_bids()"
out = mktempdir()
st = study_bids(d4; out_dir = out, task = "eyesclosed", group = :Group,
                covariates = [:Age, :Gender, :MMSE])
@test st.design.subject == bf4.subject && all(st.design.session .== "1")
@test st.design.group == ["A", "A", "C", "C", "F", "F"] && st.design.MMSE == bf4.MMSE
@test propertynames(st.design) ==
      [:subject, :session, :group, :file, :include, :Age, :Gender, :MMSE]
@test startswith(st.name, "A dataset of EEG recordings from")
h = joinpath(out, "sub-001", "eeg", "sub-001_task-eyesclosed_eeg.hdf")
@test isfile(h) && study_files(st)[1] == h
@test labels(load_neuro(h)) == labels(o4)
# existing files are not imported again
mt = mtime(h)
study_bids(d4; out_dir = out, group = :Group)
@test mtime(h) == mt
# all participants.tsv columns as covariates; no group
st4 = study_bids(d4; out_dir = out, subject = ["001", "037"])
@test all(st4.design.group .== "all") && :Group in propertynames(st4.design)
st2 = study_bids(d2; out_dir = mktempdir(), group = :outcome, covariates = [:site, :engel_score])
@test st2.design.group == ["F", "NR"] && all(st2.design.session .== "presurgery")
@test st2.design.site == ["JHH", "JHH"] && st2.design.engel_score == [4.0, -1.0]
@test_throws ArgumentError study_bids(d4; out_dir = out, task = "x")
@test_throws ArgumentError study_bids(d2; out_dir = out, group = :x)
@test_throws ArgumentError study_bids(d2; out_dir = out, covariates = [:x])


@info "Test: export_bids()"
br = mktempdir()
e_b = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edfplus.edf"))
NeuroAnalyzer.filter!(e_b; ch = "all", fprototype = :butterworth, ftype = :hp, cutoff = 0.5, order = 4)
e_b.header.recording[:bad_channel][3] = true
e_b.header.subject[:last_name] = "Secret"
add_marker!(e_b; type = "response", start = 1.0, len = 0.5, value = "a, b", ch = labels(e_b)[3])
load_locs!(e_b; file_name = joinpath(testfiles_path, "locs", "standard-10-20-cap19-elmiko.ced"))
fb = export_bids(e_b; root = br, subject = "01", task = "rest", participant = (; age = 30, sex = "F"),
                 authors = ["A. Author"], license = "CC0")
dd = joinpath(br, "sub-01", "eeg")
@test sort(basename.(fb)) == sort(["sub-01_task-rest_eeg.vhdr", "sub-01_task-rest_eeg.eeg",
                                   "sub-01_task-rest_eeg.vmrk", "sub-01_task-rest_eeg.json",
                                   "sub-01_task-rest_channels.tsv", "sub-01_task-rest_events.tsv",
                                   "sub-01_task-rest_events.json", "sub-01_electrodes.tsv",
                                   "sub-01_coordsystem.json"])
@test isfile(joinpath(br, "dataset_description.json")) && isfile(joinpath(br, "README"))
dsd = NeuroAnalyzer.JSON.parsefile(joinpath(br, "dataset_description.json"))
@test dsd["Authors"] == ["A. Author"] && dsd["License"] == "CC0"
ej = NeuroAnalyzer.JSON.parsefile(joinpath(dd, "sub-01_task-rest_events.json"))
@test haskey(ej, "value") && haskey(ej, "sample") && haskey(ej, "channel") && !haskey(ej, "onset")
js = NeuroAnalyzer.JSON.parsefile(joinpath(dd, "sub-01_task-rest_eeg.json"))
@test js["TaskName"] == "rest" && js["SamplingFrequency"] == 200 && js["RecordingType"] == "continuous"
@test js["EEGChannelCount"] == count(==("eeg"), e_b.header.recording[:channel_type])
ch = CSV.read(joinpath(dd, "sub-01_task-rest_channels.tsv"), DataFrame; missingstring = "n/a")
@test js["MISCChannelCount"] == count(==("MISC"), ch.type)
@test ch.name == labels(e_b) && ch.status[3] == "bad" && ch.units[1] == "uV" && ch.low_cutoff[1] == 0.5
cs = NeuroAnalyzer.JSON.parsefile(joinpath(dd, "sub-01_coordsystem.json"))
@test cs["EEGCoordinateSystem"] == "Other" && cs["EEGCoordinateUnits"] == "n/a"
@test !occursin("Secret", read(joinpath(dd, "sub-01_task-rest_eeg.vhdr"), String))
# round trip
b = import_bids(joinpath(dd, "sub-01_task-rest_eeg.vhdr"))
@test labels(b) == labels(e_b)
@test b.data ≈ e_b.data rtol = 1e-6
@test b.header.recording[:channel_type] == e_b.header.recording[:channel_type]
@test b.header.recording[:bad_channel] == e_b.header.recording[:bad_channel]
@test b.locs.loc_x ≈ e_b.locs.loc_x && b.locs.loc_y ≈ e_b.locs.loc_y
mb = sort(b.markers, :start)
me = sort(e_b.markers, :start)
@test mb.start ≈ me.start && mb.length ≈ me.length && mb.value == me.value
@test mb.type == me.type && mb.channel == me.channel
p = CSV.read(joinpath(br, "participants.tsv"), DataFrame; missingstring = "n/a")
@test p.participant_id == ["sub-01"] && names(p) == ["participant_id", "age", "sex"]
@test_throws ArgumentError export_bids(e_b; root = br, subject = "01", task = "rest")
# EDF, epochs
e_ep = epoch(e_b; ep_len = 10)
export_bids(e_ep; root = br, subject = "01", task = "rest", run = "2", format = :edf)
js = NeuroAnalyzer.JSON.parsefile(joinpath(dd, "sub-01_task-rest_run-2_eeg.json"))
@test js["RecordingType"] == "epoched" && js["EpochLength"] == 10
ev = CSV.read(joinpath(dd, "sub-01_task-rest_run-2_events.tsv"), DataFrame; missingstring = "n/a")
@test count(==("epoch"), ev.trial_type) == nepochs(e_ep)
# iEEG as NWB, with a session; participants: row added
o = import_bids(joinpath(d2, bf2.file[1]))
o = trim(o; seg = (0, 20), keep = true)
fi = export_bids(o; root = br, subject = "jh101", session = "pre", task = "ictal", format = :nwb,
                 sex = "U", age = "P30Y")
@test basename(fi[1]) == "sub-jh101_ses-pre_task-ictal_ieeg.nwb"
# (iEEG: electrodes.tsv is required, n/a positions without locations)
@test "sub-jh101_ses-pre_electrodes.tsv" in basename.(fi)
ii = import_bids(fi[1])
@test datatype(ii) == "ecog" && labels(ii) == labels(o)
@test count(ii.header.recording[:bad_channel]) == count(o.header.recording[:bad_channel])
@test nrow(ii.markers) == nrow(o.markers)
js = NeuroAnalyzer.JSON.parsefile(joinpath(br, "sub-jh101", "ses-pre", "ieeg", "sub-jh101_ses-pre_task-ictal_ieeg.json"))
@test js["ECOGChannelCount"] == count(==("ecog"), o.header.recording[:channel_type])
@test CSV.read(joinpath(br, "participants.tsv"), DataFrame).participant_id == ["sub-01", "sub-jh101"]
@test_throws ArgumentError export_bids(o; root = br, subject = "x", task = "t", format = :fiff)
# NIRS, optodes
n = import_bids(f3)
fn = export_bids(n; root = br, subject = "03", task = "tapping")
@test "sub-03_optodes.tsv" in basename.(fn)
nch = CSV.read(joinpath(br, "sub-03", "nirs", "sub-03_task-tapping_channels.tsv"), DataFrame)
@test names(nch)[1:6] == ["name", "type", "source", "detector", "wavelength_nominal", "units"]
@test !("notch" in names(nch))
nn = import_bids(fn[1])
@test labels(nn) == labels(n)
@test nn.data ≈ n.data rtol = 1e-5
@test nn.locs.loc_x ≈ n.locs.loc_x && nrow(nn.markers) == nrow(n.markers)
fnb = bids_files(br)
@test sort(unique(fnb.modality)) == ["eeg", "ieeg", "nirs"]
@test_throws ArgumentError export_bids(e_b; root = br, subject = "-", task = "rest")
# a study: subjects, groups and covariates
sb = mktempdir()
st = study_bids(d4; out_dir = mktempdir(), subject = ["001", "037"], group = :Group,
                covariates = [:Age])
fs_ = export_bids(st; root = sb, task = "eyesclosed")
@test bids_files(sb).subject == ["001", "037"]
p = CSV.read(joinpath(sb, "participants.tsv"), DataFrame; missingstring = "n/a")
@test p.participant_id == ["sub-001", "sub-037"] && p.group == ["A", "C"] && p.Age == [57, 57]

true
