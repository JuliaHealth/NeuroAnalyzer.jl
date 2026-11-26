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
@test o.header.subject[:id] == "002" && o.header.experiment[:name] == "FaceRecognition"
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
@test o2.header.subject[:id] == "jh101" && o2.header.experiment[:name] == "ictal"
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
@test o3.header.experiment[:name] == "tapping"
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
@test o4.header.subject[:id] == "001" && o4.header.experiment[:name] == "eyesclosed"
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

true
