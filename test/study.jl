using Test
using DataFrames
using Statistics
using GLMakie
using Random

# study files: 4 subjects (groups p and c), subject 01 with 2 sessions; epochs of
# conditions "a" and "b"; channel Fp2 is bad in subject 03
root = mktempdir()
e = epoch(import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf")); ep_len = 2)
keep_channel!(e; ch = ["Fp1", "Fp2", "Cz", "O1"])
keep_epoch!(e; ep = 1:20)
e.epochs.condition = [isodd(i) ? "a" : "b" for i = 1:nepochs(e)]
files = String[]
for (s, ses, gain) in (("01", "1", 1), ("01", "2", 3), ("02", "1", 2), ("03", "1", 4),
                       ("04", "1", 5))
    o = deepcopy(e)
    o.data .*= gain
    s == "03" && (o.header.recording[:bad_channel][2] = true)
    f = "sub-$(s)_ses-$(ses)_eeg.hdf"
    save_neuro(o; file_name = joinpath(root, f))
    push!(files, f)
end

@info "Test: study_create()"
st = study_create(files; group = ["p", "p", "c", "p", "c"], root = root,
                  covariates = (; age = [30, 30, 40, 50, 60]), name = "test")
@test st isa NeuroAnalyzer.STUDY
@test st.design.subject == ["01", "01", "02", "03", "04"]
@test st.design.session == ["1", "2", "1", "1", "1"]
@test propertynames(st.design) == [:subject, :session, :group, :file, :include, :age]
@test all(st.design.include)
@test study_files(st) == joinpath.(root, files)
@test occursin("subjects: 4, sessions: 2", sprint(show, MIME("text/plain"), st))
@test occursin("groups: p (2), c (2)", sprint(show, MIME("text/plain"), st))
# defaults without BIDS names: subject = file name, session "1", group "all"
cp(joinpath(root, files[1]), joinpath(root, "x.hdf"))
s1 = study_create(["x.hdf"]; root = root)
@test (s1.design.subject[1], s1.design.session[1], s1.design.group[1]) == ("x", "1", "all")
@test_throws ArgumentError study_create(["nofile.hdf"]; root = root)
@test_throws ArgumentError study_create([files[1], files[1]]; root = root)
@test_throws ArgumentError study_create(files; root = root, subject = fill("01", 5))
@test_throws ArgumentError study_create(files; root = root, group = ["p"])
@test_throws ArgumentError study_create(files; root = root, covariates = (; group = 1:5))
# from a design table
d = DataFrame(:subject => ["01", "02"], :session => [1, 1], :group => ["p", "c"],
              :file => files[[1, 3]])
s2 = study_create(d; root = root)
@test s2.design.session == ["1", "1"] && all(s2.design.include)
@test_throws ArgumentError study_create(select(d, Not(:file)); root = root)

@info "Test: study_check()"
c = study_check(st)
@test all(c.ok) && c.sr == fill(sr(e), 5) && c.conditions[1] == "a, b"
o = deepcopy(e)
reorder_channel!(o; order = ["Cz", "Fp1", "Fp2", "O1"])
save_neuro(o; file_name = joinpath(root, "sub-05_eeg.hdf"))
o = resample(e; new_sr = 128)
save_neuro(o; file_name = joinpath(root, "sub-06_eeg.hdf"))
s3 = study_create([files[1], "sub-05_eeg.hdf", "sub-06_eeg.hdf"]; root = root)
c = study_check(s3)
@test c.ok == [true, false, false]
@test occursin("reorder_channel()", c.issue[2])
@test occursin("resample()", c.issue[3])
@test_throws ArgumentError study_measure!(s3, :erp; ch = "all")
# continuous recordings of different lengths are compatible, epochs of different
# lengths are not
c1 = import_edf(joinpath(testfiles_path, "eeg", "eeg-test-edf.edf"))
save_neuro(trim(c1; seg = (0, 100), keep = true);
           file_name = joinpath(root, "c1.hdf"))
save_neuro(trim(c1; seg = (0, 200), keep = true);
           file_name = joinpath(root, "c2.hdf"))
save_neuro(epoch(c1; ep_len = 2); file_name = joinpath(root, "c3.hdf"))
save_neuro(epoch(c1; ep_len = 3); file_name = joinpath(root, "c4.hdf"))
@test all(study_check(study_create(["c1.hdf", "c2.hdf"]; root = root)).ok)
c = study_check(study_create(["c3.hdf", "c4.hdf"]; root = root))
@test c.ok == [true, false] && occursin("epoch()", c.issue[2])

@info "Test: study_apply!()"
hp!(obj) = NeuroAnalyzer.filter!(obj; ch = "all", fprototype = :butterworth, ftype = :hp,
                                 cutoff = 1, order = 4)
s4 = deepcopy(st)
r = study_apply!(s4, hp!; out_dir = "derivatives")
@test all(r.ok)
@test s4.design.file == joinpath.("derivatives", files)
@test all(isfile, study_files(s4))
o = load_neuro(study_files(s4)[1])
@test startswith(o.history[end], "filter(obj; ch=\"all\", fprototype=:butterworth")
@test o.data ≈ NeuroAnalyzer.filter(e; ch = "all", fprototype = :butterworth,
                                    ftype = :hp, cutoff = 1, order = 4).data
@test occursin("study_apply!(study, hp!", s4.history[end])
# existing output files
@test_throws ArgumentError study_apply!(deepcopy(st), hp!; out_dir = "derivatives")
# output folder = input folder
@test_throws ArgumentError study_apply!(deepcopy(st), hp!; out_dir = root, overwrite = true)
# a function returning a new object
s5 = deepcopy(st)
r = study_apply!(s5, obj -> keep_channel(obj; ch = ["Fp1", "Cz"]); out_dir = "d2")
@test all(r.ok)
@test nchannels(load_neuro(study_files(s5)[2])) == 2
s5 = deepcopy(st)
# an error in one recording (the first one, gain 1)
r = study_apply!(s5, obj -> obj.data[1, 1, 1] == e.data[1, 1, 1] ? error("first") : obj;
                 out_dir = "d3")
@test r.ok == [false, true, true, true, true]
@test occursin("first", r.error[1])
@test s5.design.file[1] == files[1] && s5.design.file[2] == joinpath("d3", files[2])
# excluded recordings are not processed
s5 = deepcopy(st)
s5.design.include[5] = false
r = study_apply!(s5, hp!; out_dir = "d4")
@test DataFrames.nrow(r) == 4 && s5.design.file[5] == files[5]

@info "Test: study_measure!()"
study_measure!(st, :erp; ch = "all")
m = st.measures[:erp]
@test size(m.data) == (5, 2, 4, epoch_len(e))
@test m.conditions == ["a", "b"] && m.labels == labels(e) && m.n_epochs == fill(10, 5, 2)
@test m.t == e.epoch_time && m.data_type == "eeg" && m.files == files
ea = dropdims(mean(e.data[:, :, 1:2:20]; dims = 3); dims = 3)
@test m.data[1, 1, :, :] ≈ ea
@test m.data[3, 1, :, :] ≈ 2 .* ea
eb = dropdims(mean(e.data[:, :, 2:2:20]; dims = 3); dims = 3)
@test m.data[1, 2, :, :] ≈ eb
# bad channel: NaN
@test all(isnan, m.data[4, :, 2, :]) && !any(isnan, m.data[4, :, [1, 3, 4], :])
# cached: not recomputed
nh = length(st.history)
study_measure!(st, :erp; ch = "all")
@test length(st.history) == nh
study_measure!(st, :erp; ch = "all", force = true)
@test length(st.history) == nh + 1
# changed file: recomputed
touch(joinpath(root, files[5]))
study_measure!(st, :erp; ch = "all")
@test length(st.history) == nh + 2
# rejected epochs are not used
o = load_neuro(joinpath(root, files[3]))
o.epochs.reject[1] = true
save_neuro(o; file_name = joinpath(root, files[3]), overwrite = true)
study_measure!(st, :erp; ch = "all")
@test st.measures[:erp].n_epochs[3, :] == [9, 10]
@test st.measures[:erp].data[3, 1, :, :] ≈
      2 .* dropdims(mean(e.data[:, :, 3:2:20]; dims = 3); dims = 3)
# baseline, condition, channels
s6 = deepcopy(st)
study_measure!(s6, :erp; ch = ["Cz", "Fp1"], condition = "b", bl = (0, 0.5))
m = s6.measures[:erp]
@test m.conditions == ["b"] && m.labels == ["Fp1", "Cz"]
b2 = vsearch(0.5, e.epoch_time)
@test all(abs.(mean(m.data[1, 1, :, 1:b2]; dims = 2)) .< 1e-10)
@test_throws ArgumentError study_measure!(s6, :erp; ch = "all", condition = "x")
@test_throws ArgumentError study_measure!(s6, :psd; ch = "all", db = true)
@test_throws ArgumentError study_measure!(s6, :psd; ch = "all", bl = (0, 0.5))
@test_throws ArgumentError study_measure!(s6, :erp; ch = "all", method = :mt)
@test_throws ArgumentError study_measure!(s6, :xyz; ch = "all")
# spectra and time-frequency power (linear units)
study_measure!(st, :psd; ch = "all", flim = (0, 40))
m = st.measures[:psd]
p = psd(keep_epoch(e; ep = collect(1:2:20)); ch = "all", flim = (0, 40))
@test size(m.data) == (5, 2, 4, length(p.f)) && m.f == p.f && isempty(m.t)
@test m.data[1, 1, :, :] ≈ dropdims(mean(p.p; dims = 3); dims = 3)
@test m.params == (; ch = "all", condition = nothing, flim = (0, 40))
study_measure!(st, :tf; ch = ["Fp1", "Cz"], method = :stft, wlen = 128)
m = st.measures[:tf]
sp = spectrogram(keep_epoch(e; ep = collect(1:2:20)); ch = ["Fp1", "Cz"], db = false,
                 method = :stft, wlen = 128)
@test size(m.data) == (5, 2, 2, length(sp.f), length(sp.t))
@test m.f == sp.f && m.t == sp.t
@test m.data[1, 1, 1, :, :] ≈ dropdims(mean(sp.p[:, :, 1, :]; dims = 3); dims = 3)

@info "Test: study_average()"
a = study_average(st, :erp)
x = st.measures[:erp].data
@test a.levels.group == ["p", "p", "c", "c"] && a.levels.condition == ["a", "b", "a", "b"]
@test a.levels.n == [2, 2, 2, 2]
# group p: subject 01 (sessions 1 and 2 averaged first) and subject 03
s01 = (x[1, 1, :, :] .+ x[2, 1, :, :]) ./ 2
@test a.data[1, [1, 3, 4], :] ≈ (s01[[1, 3, 4], :] .+ x[4, 1, [1, 3, 4], :]) ./ 2
@test a.sd[1, 1, :] ≈ vec(std(hcat(s01[1, :], x[4, 1, 1, :]); dims = 2))
# bad channel Fp2 of subject 03: only subject 01
@test a.data[1, 2, :] ≈ s01[2, :]
@test a.n[1, :] == [2, 1, 2, 2]
@test isnan(a.sd[1, 2, 1])
@test a.levels.n[1] == 2
# all subjects, one condition
a = study_average(st, :erp; by = Symbol[], condition = "a")
@test DataFrames.nrow(a.levels) == 1 && a.levels.n == [4]
@test a.data[1, 1, :] ≈ (s01[1, :] .+ x[3, 1, 1, :] .+ x[4, 1, 1, :] .+ x[5, 1, 1, :]) ./ 4
# sessions as levels: recordings not averaged over sessions
a = study_average(st, :erp; by = [:group, :session], condition = "a")
@test a.levels.session == ["1", "2", "1"] && a.levels.n == [2, 1, 2]
@test a.data[2, :, :] ≈ x[2, 1, :, :]
# power in dB: each subject converted before averaging
a = study_average(st, :psd; by = :group, condition = "b")
xp = st.measures[:psd].data
@test a.db
@test a.data[2, 1, 2:end] ≈
      (10 .* log10.(xp[3, 2, 1, 2:end]) .+ 10 .* log10.(xp[5, 2, 1, 2:end])) ./ 2
a = study_average(st, :psd; db = false, condition = "b")
@test !a.db && a.data[2, 1, :] ≈ (xp[3, 2, 1, :] .+ xp[5, 2, 1, :]) ./ 2
a = study_average(st, :tf)
@test size(a.data) == (4, 2, size(st.measures[:tf].data)[4:5]...)
# excluded recordings are left out
s7 = deepcopy(st)
s7.design.include[4] = false
a = study_average(s7, :erp; condition = "a")
@test a.levels.n == [1, 2] && a.data[1, :, :] ≈ s01
@test_throws ArgumentError study_average(st, :erp; by = :xyz)
@test_throws ArgumentError study_average(st, :erp; condition = "x")
@test_throws ArgumentError study_average(study_create(files; root = root), :erp)
# included recordings without the measure
s7.design.include[4] = true
s7.design.file[4] = "x.hdf"
@test_throws ArgumentError study_average(s7, :erp)

@info "Test: study_erp()"
g = study_erp(st; condition = "a", group = "p")
@test datatype(g) == "erp" && nepochs(g) == 3
@test g.epochs.source == ["grand average", "01", "03"]
@test g.data[:, :, 1] ≈ study_average(st, :erp; condition = "a").data[1, :, :]
@test g.data[:, :, 2] ≈ s01
# bad channel of subject 03: the grand average
@test g.data[2, :, 3] ≈ g.data[2, :, 1]
@test g.epoch_time == e.epoch_time && labels(g) == labels(e)
@test !any(g.header.recording[:bad_channel]) && isempty(g.markers)
@test plot_erp(g; ch = "Fp1", gui = false) isa GLMakie.Figure
@test nepochs(study_erp(st; condition = "b")) == 5
@test_throws ArgumentError study_erp(st)
@test_throws ArgumentError study_erp(st; condition = "a", group = "x")

@info "Test: study_stats()"
x = st.measures[:erp].data
s01a, s01b = (x[1, 1, :, :] .+ x[2, 1, :, :]) ./ 2, (x[1, 2, :, :] .+ x[2, 2, :, :]) ./ 2
# between groups (p: 01, 03; c: 02, 04), one condition, no correction
r = study_stats(st, :erp; condition = "a", correction = :none)
gp = permutedims(cat(s01a, x[4, 1, :, :]; dims = 3), (3, 1, 2))
gc = permutedims(cat(x[3, 1, :, :], x[5, 1, :, :]; dims = 3), (3, 1, 2))
mt = mass_ttest(gp, gc)
@test r.test === :ttest && r.levels == ["p", "c"] && r.n == [2, 2]
@test isequal(r.stat, mt.t) && isequal(r.p, mt.p) && isequal(r.p_adj, mt.p)
@test size(r.stat) == (4, epoch_len(e)) && isempty(r.clusters)
# bad channel Fp2 of subject 03: one subject in group p, no test
@test all(isnan, r.stat[2, :])
r = study_stats(st, :erp; condition = "a", correction = :fdr)
@test isequal(r.p_adj, p_adjust(mt.p)) && r.mask == (.!isnan.(r.p_adj) .& (r.p_adj .< 0.05))
r = study_stats(st, :erp; condition = "a", levels = ["c", "p"], correction = :none)
@test isequal(r.stat, -mt.t)
# within subjects: conditions (paired)
r = study_stats(st, :erp; factor = :condition, correction = :none)
ya = permutedims(cat(s01a, x[3, 1, :, :], x[4, 1, :, :], x[5, 1, :, :]; dims = 3), (3, 1, 2))
yb = permutedims(cat(s01b, x[3, 2, :, :], x[4, 2, :, :], x[5, 2, :, :]; dims = 3), (3, 1, 2))
mt = mass_ttest(ya, yb; paired = true)
@test r.test === :ttest_paired && r.levels == ["a", "b"] && r.n == [4]
@test isequal(r.stat, mt.t)
r = study_stats(st, :erp; factor = :condition, group = "c", correction = :none)
@test r.n == [2]
# cluster-based permutation test
r = study_stats(st, :erp; factor = :condition, nperm = 100, rng = Random.Xoshiro(1),
                ch_adj = trues(4, 4))
@test r.correction === :cluster && isequal(r.stat, mt.t)
@test DataFrames.nrow(r.clusters) == maximum(r.cluster_id)
@test r.mask == (.!isnan.(r.p_adj) .& (r.p_adj .< 0.05))
# spectra in dB
r = study_stats(st, :psd; factor = :condition, correction = :none)
xp = 10 .* log10.(st.measures[:psd].data)
@test r.db && r.stat[1, 3] ≈ mass_ttest(
    [(xp[1, 1, 1, 3] + xp[2, 1, 1, 3]) / 2; xp[3:5, 1, 1, 3]],
    [(xp[1, 2, 1, 3] + xp[2, 2, 1, 3]) / 2; xp[3:5, 2, 1, 3]]; paired = true).t[]
r = study_stats(st, :tf; condition = "a", correction = :holm)
@test size(r.stat) == size(st.measures[:tf].data)[3:end]
# sessions: only subject 01 has two sessions
@test_throws ArgumentError study_stats(st, :erp; factor = :session, condition = "a")
@test_throws ArgumentError study_stats(st, :erp)
@test_throws ArgumentError study_stats(st, :erp; condition = "a", factor = :xyz)
@test_throws ArgumentError study_stats(st, :erp; condition = "a", levels = ["p"])
@test_throws ArgumentError study_stats(st, :erp; factor = :condition, condition = "a")
@test_throws ArgumentError study_stats(st, :erp; condition = "a", group = "p")
@test_throws ArgumentError study_stats(st, :erp; factor = :subject, condition = "a")

@info "Test: save_study() / load_study()"
sf = joinpath(root, "study.hdf")
save_study(st; file_name = sf)
@test_throws ArgumentError save_study(st; file_name = sf)
@test_throws ArgumentError save_study(st; file_name = joinpath(root, "study.jld2"))
s8 = load_study(sf)
@test s8.name == "test" && s8.root == st.root && s8.design == st.design
@test s8.history == st.history
@test sort(collect(keys(s8.measures))) == [:erp, :psd, :tf]
@test isequal(s8.measures[:erp].data, st.measures[:erp].data)
@test s8.measures[:psd].params == st.measures[:psd].params
# measures are still up to date
nh = length(s8.history)
study_measure!(s8, :psd; ch = "all", flim = (0, 40))
@test length(s8.history) == nh
# moved study
s9 = load_study(sf; root = root)
@test s9.root == abspath(root)
@test_throws ArgumentError load_study(joinpath(root, files[1]))
@test_throws ArgumentError load_neuro(sf)

true
