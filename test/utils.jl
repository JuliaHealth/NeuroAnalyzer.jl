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

@info "test 1/68: apply()"
@test size(apply(e10, f="mean(obj, dims=1)")) == (23, 1, 10)

@info "test 2/68: l1()"
@test l1(a1, a2) == 12

@info "test 3/68: l2()"
@test l2(a1, a2) == 3.4641016151377544

@info "test 4/68: perm_cmp()"
z, b = perm_cmp(a1, a2)
@test round.(z) == -2.0 .* ones(2, 3)
@test size(b) == (2, 3)

@info "test 5/68: add_component()"
x = collect(1:10)
e10_tmp = add_component(e10, c=:x, v=x)
@test e10_tmp.components[:x] == x

@info "test 6/68: component_type()"
@test component_type(e10_tmp, c=:x) == Vector{Int64}

@info "test 7/68: extract_component()"
@test extract_component(e10_tmp, c=:x) == x

@info "test 8/68: delete_component()"
e10_tmp2 = deepcopy(e10_tmp)
delete_component!(e10_tmp2, c=:x)
@test e10_tmp2.components == Dict()

@info "test 9/68: rename_component()"
rename_component!(e10_tmp, c_old=:x, c_new=:y)
@test e10_tmp.components[:y] == 1:10

@info "test 10/68: reset_components()"
reset_components!(e10_tmp)
@test e10_tmp.components == Dict()

@info "test 11/68: vsearch()"
@test vsearch(2.1, [1, 2, 3, 4]) == 2
@test vsearch([2.1, 2.9], [1, 2, 3, 4]) == [2.0, 3.0]

@info "test 12/68: trim()"
@test vsplit(1:10, 2) == [[1, 6], [2, 7], [3, 8], [4, 9], [5, 10]]

@info "test 13/68: fft0()"
x = fft0(v1, 2)
@test length(x) == 7

@info "test 14/68: ifft0()"
@test round.(real.(ifft0(x, 2))) == [1.0, 2.0, 3.0, 4.0, 5.0]

@info "test 15/68: fft2()"
x = fft2(v1)
@test length(x) == 8

@info "test 16/68: nextpow()"
@test nextpow2(5) == 8

@info "test 17/68: gradient()"
g, g_l = NeuroAnalyzer.gradient(rand(10))
@test size(g) == (10, )
@test size(g_l) == (10, )
g, g_l = NeuroAnalyzer.gradient(rand(10, 10))
@test size(g) == (10, 10)
@test size(g_l) == (10, 10)
g, g_l = NeuroAnalyzer.gradient(rand(10, 10, 10))
@test size(g) == (10, 10, 10)
@test size(g_l) == (10, 10, 10)

@info "test 18/68: findpeaks()"
x = rand(100)
x[10] *= 1000
@test 10 in findpeaks(x)

@info "test 19/68: hz2rads()"
@test hz2rads(1) == 2pi

@info "test 20/68: rads2hz()"
@test rads2hz(2pi) == 1

@info "test 21/68: t2f()"
@test t2f(1000) == 1.0

@info "test 22/68: f2t()"
@test f2t(1.0) == 1000.0

@info "test 23/68: freqs()"
f, nf = freqs(0:1/10:10)
@test length(f) == 51
@test nf == 5
f, nf = freqs(rand(100), 10)
@test length(f) == 51
@test nf == 5
f, nf = freqs(e10)
@test length(f) == 1281
@test nf == 128

@info "test 24/68: generate_window()"
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

@info "test 25/68: generate_sine()"
s = generate_sine(10, 1:100)
@test length(s) == 100

@info "test 26/68: generate_csine()"
s = generate_csine(10, 1:100)
@test length(s) == 100

@info "test 27/68: generate_sinc()"
s = generate_sinc()
@test length(s) == 401

@info "test 28/68: generate_morlet()"
s = generate_morlet(100, 10)
@test length(s) == 201

@info "test 29/68: generate_gaussian()"
s = generate_gaussian(100, 10)
@test length(s) == 201

@info "test 30/68: generate_noise()"
s = generate_noise(100)
@test length(s) == 100

@info "test 31/68: generate_morlet_fwhm()"
s = generate_morlet_fwhm(100, 10)
@test length(s) == 201

@info "test 32/68: sr()"
@test sr(e10) == 256

@info "test 33/68: nchannels()"
@test nchannels(e10) == 24

@info "test 34/68: nepochs()"
@test nepochs(e10) == 10

@info "test 35/68: signal_len()"
@test signal_len(e10) == 25600

@info "test 36/68: epoch_len()"
@test epoch_len(e10) == 2560

@info "test 37/68: signal_channels()"
@test signal_channels("eeg", ["eeg", "eeg", "ecg", "mrk"]) == [1, 2]
@test signal_channels(e10) == 1:23

@info "test 38/68: get_channel_bytype()"
@test get_channel_bytype(["eeg", "ecg", "mrk"], type="eeg") == 1
@test get_channel_bytype(e10, type="eeg") == 1:19

@info "test 39/68: history()"
@test NeuroAnalyzer.history(e10) isa Vector{String}

@info "test 40/68: labels()"
@test length(labels(e10)) == 24

@info "test 41/68: channel_cluster()"
@test channel_cluster(e10, cluster=:f1) == [1, 3, 11]

@info "test 42/68: band_frq()"
@test band_frq(256, band=:alpha) == (8.0, 13.0)
@test band_frq(e10, band=:alpha) == (8.0, 13.0)

@info "test 43/68: m_pad0()"
@test m_pad0(m1) == [1 2 3; 4 5 6; 0 0 0]

@info "test 44/68: m_sortperm()"
@test m_sortperm(m1) == [1 1 1; 2 2 2]

@info "test 45/68: m_sort()"
@test m_sort(m1, [2, 1]) == [4 5 6; 1 2 3]

@info "test 46/68: m_norm()"
@test m_norm(m1) == [0.5 1.0 1.5; 2.0 2.5 3.0]

@info "test 47/68: linspace()"
@test linspace(1 , 10, 10) == 1:10

@info "test 48/68: logspace()"
@test logspace(1, 2, 2) == [10.0, 100.0]

@info "test 49/68: cmax()"
@test cmax([1+2im, 10+10im]) == 10+10im

@info "test 50/68: cmin()"
@test cmin([1+2im, 10+10im]) == 1+2im

@info "test 51/68: tuple_order()"
@test tuple_order((2, 1)) == (1, 2)

@info "test 52/68: cums()"
@test cums(a1) == [1.0 2.0 3.0; 1.0 2.0 3.0;;; 1.0 2.0 3.0; 1.0 2.0 3.0]

@info "test 53/68: f_nearest()"
m = [(1.0, 2.0) (3.0, 4.0); (5.0, 6.0) (7.0, 8.0)]
@test f_nearest(m, (3.2, 3.9)) == (1, 2)

@info "test 54/68: view_note()"
@test view_note(e10) == ""

@info "test 55/68: add_note()"
add_note!(e10, note="test")
@test view_note(e10) == "test"

@info "test 56/68: delete_note()"
delete_note!(e10)
@test view_note(e10) == ""

@info "test 57/68: pad0()"
@test pad0(v1, 5) == [1, 2, 3, 4, 5, 0, 0, 0, 0, 0]
@test pad0(m1, 5) == [1 2 3 0 0 0 0 0; 4 5 6 0 0 0 0 0]

@info "test 58/68: pad2()"
@test pad2(v1) == [1, 2, 3, 4, 5, 0, 0, 0]
@test pad2(m1) == [1 2 3 0; 4 5 6 0]

@info "test 59/68: phases()"
@test length(phases(e10.data[1, :, 1])) == 2560

@info "test 60/68: pick()"
@test NeuroAnalyzer.pick(e10, p=[:l, :f]) == [1, 3, 11]

@info "test 61/68: t2s()"
@test t2s(1.0, 256) == 256
@test t2s(e10, t=1.0) == 257

@info "test 62/68: s2t()"
@test s2t(2560, 256) == 10.0
@test s2t(e10, s=256) == 1.0

@info "test 63/68: get_channel_bywl()"
@test get_channel_bywl(n, wl=760) == 1:36

@info "test 64/68: size()"
@test size(e10) == (24, 2560, 10)

@info "test 65/68: to_df()"
@test to_df(eeg) isa DataFrame

@info "test 66/68: types()"
@test length(types(e10)) == 24

@info "test 67/68: optode_labels()"
n = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
@test length(optode_labels(n)) == 3
@test length(source_labels(n)) == 1
@test length(detector_labels(n)) == 2

@info "test 68/68: datatype()"
@test datatype(e10) == "eeg"

true