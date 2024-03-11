using NeuroAnalyzer
using Test

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

@info "test 1/74: apply()"
@test size(apply(e10, f="mean(obj, dims=1)")) == (23, 1, 10)

@info "test 2/74: l1()"
@test l1(a1, a2) == 12

@info "test 3/74: l2()"
@test l2(a1, a2) == 3.4641016151377544

@info "test 4/74: perm_cmp()"
z, b = perm_cmp(a1, a2)
@test round.(z) == -2.0 .* ones(2, 3)
@test size(b) == (2, 3)

@info "test 5/74: add_component()"
x = collect(1:10)
e10_tmp = add_component(e10, c=:x, v=x)
@test e10_tmp.components[:x] == x

@info "test 6/74: component_type()"
@test component_type(e10_tmp, c=:x) == Vector{Int64}

@info "test 7/74: extract_component()"
@test extract_component(e10_tmp, c=:x) == x

@info "test 8/74: delete_component()"
e10_tmp2 = deepcopy(e10_tmp)
delete_component!(e10_tmp2, c=:x)
@test e10_tmp2.components == Dict()

@info "test 9/74: rename_component()"
rename_component!(e10_tmp, c_old=:x, c_new=:y)
@test e10_tmp.components[:y] == 1:10

@info "test 10/74: reset_components()"
reset_components!(e10_tmp)
@test e10_tmp.components == Dict()

@info "test 11/74: vsearch()"
@test vsearch(2.1, [1, 2, 3, 4]) == 2
@test vsearch([2.1, 2.9], [1, 2, 3, 4]) == [2.0, 3.0]

@info "test 12/74: trim()"
@test vsplit(1:10, 2) == [[1, 6], [2, 7], [3, 8], [4, 9], [5, 10]]

@info "test 13/74: fft0()"
x = fft0(v1)
@test length(x) == 5
x = fft0(v1, 10)
@test length(x) == 15
x = rfft0(v1)
@test length(x) == 3
x = rfft0(v1, 10)
@test length(x) == div(length(v1) + 10, 2) + 1

@info "test 14/74: ifft0()"
x = fft0(v1, 2)
@test round.(real.(ifft0(x, 2))) == [1.0, 2.0, 3.0, 4.0, 5.0]

@info "test 15/74: fft2()"
x = fft2(v1)
@test length(x) == 8
x = rfft2(v1)
@test length(x) == 5

@info "test 16/74: nextpow()"
@test nextpow2(5) == 8

@info "test 17/74: gradient()"
g, g_l = NeuroAnalyzer.gradient(rand(10))
@test size(g) == (10, )
@test size(g_l) == (10, )
g, g_l = NeuroAnalyzer.gradient(rand(10, 10))
@test size(g) == (10, 10)
@test size(g_l) == (10, 10)
g, g_l = NeuroAnalyzer.gradient(rand(10, 10, 10))
@test size(g) == (10, 10, 10)
@test size(g_l) == (10, 10, 10)

@info "test 18/74: findpeaks()"
x = rand(100)
x[10] *= 1000
@test 10 in findpeaks(x)

@info "test 19/74: hz2rads()"
@test hz2rads(1) == 2pi

@info "test 20/74: rads2hz()"
@test rads2hz(2pi) == 1

@info "test 21/74: t2f()"
@test t2f(1000) == 1.0

@info "test 22/74: f2t()"
@test f2t(1.0) == 1000.0

@info "test 23/74: freqs()"
f, nf = freqs(0:1/10:10)
@test length(f) == 51
@test nf == 5
f, nf = freqs(rand(100), 10)
@test length(f) == 51
@test nf == 5
f, nf = freqs(e10)
@test length(f) == 1281
@test nf == 128

@info "test 24/74: generate_window()"
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

@info "test 25/74: generate_sine()"
s = generate_sine(10, 1:100)
@test length(s) == 100

@info "test 26/74: generate_csine()"
s = generate_csine(10, 1:100)
@test length(s) == 100

@info "test 27/74: generate_sinc()"
s = generate_sinc()
@test length(s) == 401

@info "test 28/74: generate_morlet()"
s = generate_morlet(100, 10)
@test length(s) == 201

@info "test 29/74: generate_gaussian()"
s = generate_gaussian(100, 10)
@test length(s) == 201

@info "test 30/74: generate_noise()"
s = generate_noise(100)
@test length(s) == 100

@info "test 31/74: generate_morlet_fwhm()"
s = generate_morlet_fwhm(100, 10)
@test length(s) == 201

@info "test 32/74: sr()"
@test sr(e10) == 256

@info "test 33/74: nchannels()"
@test nchannels(e10) == 24

@info "test 34/74: nepochs()"
@test nepochs(e10) == 10

@info "test 35/74: signal_len()"
@test signal_len(e10) == 25600

@info "test 36/74: epoch_len()"
@test epoch_len(e10) == 2560

@info "test 37/74: signal_channels()"
@test signal_channels("eeg", ["eeg", "eeg", "ecg", "mrk"]) == [1, 2]
@test signal_channels(e10) == 1:23

@info "test 38/74: get_channel_bytype()"
@test get_channel_bytype(["eeg", "ecg", "mrk"], type="eeg") == 1
@test get_channel_bytype(e10, type="eeg") == 1:19

@info "test 39/74: history()"
@test NeuroAnalyzer.history(e10) isa Vector{String}

@info "test 40/74: labels()"
@test length(labels(e10)) == 24

@info "test 41/74: channel_cluster()"
@test channel_cluster(e10, cluster=:f1) == [1, 3, 11]

@info "test 42/74: band_frq()"
@test band_frq(256, band=:alpha) == (8.0, 13.0)
@test band_frq(e10, band=:alpha) == (8.0, 13.0)

@info "test 43/74: m_pad0()"
@test m_pad0(m1) == [1 2 3; 4 5 6; 0 0 0]

@info "test 44/74: m_sortperm()"
@test m_sortperm(m1) == [1 1 1; 2 2 2]

@info "test 45/74: m_sort()"
@test m_sort(m1, [2, 1]) == [4 5 6; 1 2 3]

@info "test 46/74: m_norm()"
@test m_norm(m1) == [0.5 1.0 1.5; 2.0 2.5 3.0]

@info "test 47/74: linspace()"
@test linspace(1 , 10, 10) == 1:10

@info "test 48/74: logspace()"
@test logspace(1, 2, 2) == [10.0, 100.0]

@info "test 49/74: cmax()"
@test cmax([1+2im, 10+10im]) == 10+10im

@info "test 50/74: cmin()"
@test cmin([1+2im, 10+10im]) == 1+2im

@info "test 51/74: tuple_order()"
@test tuple_order((2, 1)) == (1, 2)

@info "test 52/74: cums()"
@test cums(a1) == [1.0 2.0 3.0; 1.0 2.0 3.0;;; 1.0 2.0 3.0; 1.0 2.0 3.0]

@info "test 53/74: f_nearest()"
m = [(1.0, 2.0) (3.0, 4.0); (5.0, 6.0) (7.0, 8.0)]
@test f_nearest(m, (3.2, 3.9)) == (1, 2)

@info "test 54/74: view_note()"
@test view_note(e10) == ""

@info "test 55/74: add_note()"
add_note!(e10, note="test")
@test view_note(e10) == "test"

@info "test 56/74: delete_note()"
delete_note!(e10)
@test view_note(e10) == ""

@info "test 57/74: pad0()"
@test pad0(v1, 5) == [1, 2, 3, 4, 5, 0, 0, 0, 0, 0]
@test pad0(m1, 5) == [1 2 3 0 0 0 0 0; 4 5 6 0 0 0 0 0]

@info "test 58/74: pad2()"
@test pad2(v1) == [1, 2, 3, 4, 5, 0, 0, 0]
@test pad2(m1) == [1 2 3 0; 4 5 6 0]

@info "test 59/74: phases()"
@test length(phases(e10.data[1, :, 1])) == 2560

@info "test 60/74: pick()"
@test NeuroAnalyzer.pick(e10, p=[:l, :f]) == [1, 3, 11]

@info "test 61/74: t2s()"
@test t2s(1.0, 256) == 256
@test t2s(e10, t=1.0) == 257

@info "test 62/74: s2t()"
@test s2t(2560, 256) == 10.0
@test s2t(e10, s=256) == 1.0

@info "test 63/74: get_channel_bywl()"
@test get_channel_bywl(n, wl=760) == 1:36

@info "test 64/74: size()"
@test size(e10) == (24, 2560, 10)

@info "test 65/74: to_df()"
@test to_df(eeg) isa DataFrame

@info "test 66/74: chtypes()"
@test length(chtypes(e10)) == 24

@info "test 67/74: optode_labels()"
n = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
@test length(optode_labels(n)) == 3
@test length(source_labels(n)) == 1
@test length(detector_labels(n)) == 2

@info "test 68/74: delmean()"
@test delmean(v1) == [-2, -1, 0, 1, 2]
@test delmean(a1, dims=3) == zeros(2, 3, 2)

@info "test 69/74: tavg()"
@test tavg(a1) == ones(2, 3, 1)

@info "test 70/74: padm()"
@test padm(ones(2), 2) == ones(4)
@test padm(ones(2, 4), 2, mode=:all) == ones(2, 6)
@test padm(ones(2, 4), 2, mode=:row) == ones(2, 6)
@test padm(ones(2, 4, 3), 2, mode=:all) == ones(2, 6, 3)
@test padm(ones(2, 4, 3), 2, mode=:row) == ones(2, 6, 3)

@info "test 71/74: vec2mat()"
@test size(vec2mat(ones(10), wlen=2, woverlap=2)) == (5, 2)

@info "test 72/74: arr2mat()"
@test size(arr2mat(rand(1, 10, 10))) == (10, 10)

@info "test 73/74: minat()"
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
@test minat(v1, v2) == (6, 1)

@info "test 74/74: maxat()"
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
@test maxat(v1, v2) == (2, 5)

true