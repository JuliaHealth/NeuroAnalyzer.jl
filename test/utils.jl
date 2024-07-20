using NeuroAnalyzer
using Test

ntests = 79

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
n = import_snirf(joinpath(testfiles_path, "fnirs-test-snirf.snirf"))
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

@info "Test 1/$ntests: apply()"
@test size(apply(e10, f="mean(obj, dims=1)")) == (23, 1, 10)

@info "Test 2/$ntests: l1()"
@test l1(a1, a2) == 12

@info "Test 3/$ntests: l2()"
@test l2(a1, a2) == 3.4641016151377544

@info "Test 4/$ntests: perm_cmp()"
z, b = perm_cmp(a1, a2)
@test round.(z) == -2.0 .* ones(2, 3)
@test size(b) == (2, 3)

@info "Test 5/$ntests: add_component()"
x = collect(1:10)
e10_tmp = add_component(e10, c=:x, v=x)
@test e10_tmp.components[:x] == x

@info "Test 6/$ntests: component_type()"
@test component_type(e10_tmp, c=:x) == Vector{Int64}

@info "Test 7/$ntests: extract_component()"
@test extract_component(e10_tmp, c=:x) == x

@info "Test 8/$ntests: delete_component()"
e10_tmp2 = deepcopy(e10_tmp)
delete_component!(e10_tmp2, c=:x)
@test e10_tmp2.components == Dict()

@info "Test 9/$ntests: rename_component()"
rename_component!(e10_tmp, c_old=:x, c_new=:y)
@test e10_tmp.components[:y] == 1:10

@info "Test 10/$ntests: reset_components()"
reset_components!(e10_tmp)
@test e10_tmp.components == Dict()

@info "Test 11/$ntests: vsearch()"
@test vsearch(2.1, [1, 2, 3, 4]) == 2
@test vsearch([2.1, 2.9], [1, 2, 3, 4]) == [2.0, 3.0]

@info "Test 12/$ntests: trim()"
@test vsplit(1:10, 2) == [[1, 6], [2, 7], [3, 8], [4, 9], [5, 10]]

@info "Test 13/$ntests: fft0()"
x = fft0(v1)
@test length(x) == 5
x = fft0(v1, 10)
@test length(x) == 15
x = rfft0(v1)
@test length(x) == 3
x = rfft0(v1, 10)
@test length(x) == div(length(v1) + 10, 2) + 1

@info "Test 14/$ntests: ifft0()"
x = fft0(v1, 2)
@test round.(real.(ifft0(x, 2))) == [1.0, 2.0, 3.0, 4.0, 5.0]

@info "Test 15/$ntests: fft2()"
x = fft2(v1)
@test length(x) == 8
x = rfft2(v1)
@test length(x) == 5

@info "Test 16/$ntests: nextpow()"
@test nextpow2(5) == 8

@info "Test 17/$ntests: gradient()"
g, g_l = NeuroAnalyzer.gradient(rand(10))
@test size(g) == (10, )
@test size(g_l) == (10, )
g, g_l = NeuroAnalyzer.gradient(rand(10, 10))
@test size(g) == (10, 10)
@test size(g_l) == (10, 10)
g, g_l = NeuroAnalyzer.gradient(rand(10, 10, 10))
@test size(g) == (10, 10, 10)
@test size(g_l) == (10, 10, 10)

@info "Test 18/$ntests: findpeaks()"
x = rand(100)
x[10] *= 1000
@test 10 in findpeaks(x)

@info "Test 19/$ntests: hz2rads()"
@test hz2rads(1) == 2pi

@info "Test 20/$ntests: rads2hz()"
@test rads2hz(2pi) == 1

@info "Test 21/$ntests: t2f()"
@test t2f(1000) == 1.0

@info "Test 22/$ntests: f2t()"
@test f2t(1.0) == 1000.0

@info "Test 23/$ntests: freqs()"
f, nf = NeuroAnalyzer.freqs(0:1/10:10)
@test length(f) == 51
@test nf == 5
f, nf = NeuroAnalyzer.freqs(rand(100), 10)
@test length(f) == 51
@test nf == 5
f, nf = NeuroAnalyzer.freqs(e10)
@test length(f) == 1281
@test nf == 128

@info "Test 24/$ntests: generate_window()"
s = generate_window(:hann, 100)
@test length(s) == 100
s = generate_window(:bh, 100)
@test length(s) == 100
s = generate_window(:bohman, 100)
@test length(s) == 100
s = generate_window(:flat, 100)
@test length(s) == 100
s = generate_window(:bn, 100)
@test length(s) == 100
s = generate_window(:nutall, 100)
@test length(s) == 100
s = generate_window(:triangle, 100)
@test length(s) == 100
s = generate_window(:exp, 100)
@test length(s) == 100

@info "Test 25/$ntests: generate_sine()"
s = generate_sine(10, 1:100)
@test length(s) == 100

@info "Test 26/$ntests: generate_csine()"
s = generate_csine(10, 1:100)
@test length(s) == 100

@info "Test 27/$ntests: generate_sinc()"
s = generate_sinc()
@test length(s) == 401

@info "Test 28/$ntests: generate_morlet()"
s = generate_morlet(100, 10)
@test length(s) == 201

@info "Test 29/$ntests: generate_gaussian()"
s = generate_gaussian(100, 10)
@test length(s) == 201

@info "Test 30/$ntests: generate_noise()"
s = generate_noise(100)
@test length(s) == 100

@info "Test 31/$ntests: generate_morlet_fwhm()"
s = generate_morlet_fwhm(100, 10)
@test length(s) == 201

@info "Test 32/$ntests: sr()"
@test sr(e10) == 256

@info "Test 33/$ntests: nchannels()"
@test nchannels(e10) == 24

@info "Test 34/$ntests: nepochs()"
@test nepochs(e10) == 10

@info "Test 35/$ntests: signal_len()"
@test signal_len(e10) == 25600

@info "Test 36/$ntests: epoch_len()"
@test epoch_len(e10) == 2560

@info "Test 37/$ntests: signal_channels()"
@test signal_channels("eeg", ["eeg", "eeg", "ecg", "mrk"]) == [1, 2]
@test signal_channels(e10) == 1:23

@info "Test 38/$ntests: get_channel_bytype()"
@test get_channel_bytype(["eeg", "ecg", "mrk"], type="eeg") == 1
@test get_channel_bytype(e10, type="eeg") == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 20, 21]

@info "Test 39/$ntests: history()"
@test NeuroAnalyzer.history(e10) isa Vector{String}

@info "Test 40/$ntests: labels()"
@test length(labels(e10)) == 24

@info "Test 41/$ntests: channel_cluster()"
@test channel_cluster(e10, cluster=:f1) == [1, 3, 11]

@info "Test 42/$ntests: band_frq()"
@test band_frq(256, band=:alpha) == (8.0, 13.0)
@test band_frq(e10, band=:alpha) == (8.0, 13.0)

@info "Test 43/$ntests: m_pad0()"
@test m_pad0(m1) == [1 2 3; 4 5 6; 0 0 0]

@info "Test 44/$ntests: m_sortperm()"
@test m_sortperm(m1) == [1 1 1; 2 2 2]

@info "Test 45/$ntests: m_sort()"
@test m_sort(m1, [2, 1]) == [4 5 6; 1 2 3]

@info "Test 46/$ntests: m_norm()"
@test m_norm(m1) == [0.5 1.0 1.5; 2.0 2.5 3.0]

@info "Test 47/$ntests: linspace()"
@test linspace(1 , 10, 10) == 1:10

@info "Test 48/$ntests: logspace()"
@test logspace(1, 2, 2) == [10.0, 100.0]

@info "Test 49/$ntests: cmax()"
@test cmax([1+2im, 10+10im]) == 10+10im

@info "Test 50/$ntests: cmin()"
@test cmin([1+2im, 10+10im]) == 1+2im

@info "Test 51/$ntests: tuple_order()"
@test tuple_order((2, 1)) == (1, 2)

@info "Test 52/$ntests: cums()"
@test cums(a1) == [1.0 2.0 3.0; 1.0 2.0 3.0;;; 1.0 2.0 3.0; 1.0 2.0 3.0]

@info "Test 53/$ntests: f_nearest()"
m = [(1.0, 2.0) (3.0, 4.0); (5.0, 6.0) (7.0, 8.0)]
@test f_nearest(m, (3.2, 3.9)) == (1, 2)

@info "Test 54/$ntests: view_note()"
@test view_note(e10) == ""

@info "Test 55/$ntests: add_note()"
add_note!(e10, note="test")
@test view_note(e10) == "test"

@info "Test 56/$ntests: delete_note()"
delete_note!(e10)
@test view_note(e10) == ""

@info "Test 57/$ntests: pad0()"
@test pad0(v1, 5) == [1, 2, 3, 4, 5, 0, 0, 0, 0, 0]
@test pad0(m1, 5) == [1 2 3 0 0 0 0 0; 4 5 6 0 0 0 0 0]

@info "Test 58/$ntests: pad2()"
@test pad2(v1) == [1, 2, 3, 4, 5, 0, 0, 0]
@test pad2(m1) == [1 2 3 0; 4 5 6 0]

@info "Test 59/$ntests: phases()"
@test length(phases(e10.data[1, :, 1])) == 2560

@info "Test 60/$ntests: pick()"
@test NeuroAnalyzer.channel_pick(e10, p=[:l, :f]) == [1, 3, 11]

@info "Test 61/$ntests: t2s()"
@test t2s(1.0, 256) == 256
@test t2s(e10, t=1.0) == 257

@info "Test 62/$ntests: s2t()"
@test s2t(2560, 256) == 10.0
@test s2t(e10, s=256) == 1.0

@info "Test 63/$ntests: get_channel_bywl()"
@test get_channel_bywl(n, wl=760) == 1:36

@info "Test 64/$ntests: size()"
@test size(e10) == (24, 2560, 10)

@info "Test 65/$ntests: to_df()"
@test to_df(eeg) isa DataFrame

@info "Test 66/$ntests: chtypes()"
@test length(chtypes(e10)) == 24

@info "Test 67/$ntests: optode_labels()"
n = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
@test length(optode_labels(n)) == 3
@test length(source_labels(n)) == 1
@test length(detector_labels(n)) == 2

@info "Test 68/$ntests: delmean()"
@test delmean(v1) == [-2, -1, 0, 1, 2]
@test delmean(a1, dims=3) == zeros(2, 3, 2)

@info "Test 69/$ntests: tavg()"
@test tavg(a1) == ones(2, 3, 1)

@info "Test 70/$ntests: padm()"
@test padm(ones(2), 2) == ones(4)
@test padm(ones(2, 4), 2, mode=:all) == ones(2, 6)
@test padm(ones(2, 4), 2, mode=:row) == ones(2, 6)
@test padm(ones(2, 4, 3), 2, mode=:all) == ones(2, 6, 3)
@test padm(ones(2, 4, 3), 2, mode=:row) == ones(2, 6, 3)

@info "Test 71/$ntests: vec2mat()"
@test size(vec2mat(ones(10), wlen=2, woverlap=2)) == (5, 2)

@info "Test 72/$ntests: arr2mat()"
@test size(arr2mat(rand(1, 10, 10))) == (10, 10)

@info "Test 73/$ntests: minat()"
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
@test minat(v1, v2) == (6, 1)

@info "Test 74/$ntests: maxat()"
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
@test maxat(v1, v2) == (2, 5)

@info "Test 75/$ntests: paired_labels()"
l = ["ch1", "ch2", "ch3"]
@test length(paired_labels(l; unq=true)) == 6
@test length(paired_labels(l; unq=false)) == 9
@test length(paired_labels(l, l)) == 3

@info "Test 76/$ntests: vreduce()"
x = rand(100)
f = linspace(0, 10, 100)
x2, f2 = vreduce(x, f)
@test length(x2) == length(f2) == 21

@info "Test 77/$ntests: areduce()"
x = rand(10, 100)
f = linspace(0, 10, 100)
x2, f2 = areduce(x, f)
@test size(x2, 2) == length(f2) == 21
x = rand(10, 100, 5)
f = linspace(0, 10, 100)
x2, f2 = areduce(x, f)
@test size(x2, 2) == length(f2) == 21

@info "Test 78/$ntests: generate_signal()"
s = generate_signal(100)
@test length(s) == 100

@info "Test 79/$ntests: cwtfrq()"
s = rand(100)
@test length(cwtfrq(s, fs=10)) == 130
@test length(cwtfrq(e10)) == 19

true