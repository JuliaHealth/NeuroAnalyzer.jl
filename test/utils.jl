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

@info "Test: apply()"
@test size(apply(e10, ch="all", f="mean(obj, dims=1)")) == (24, 1, 10)

@info "Test: l1()"
@test l1(a1, a2) == 12

@info "Test: l2()"
@test l2(a1, a2) == 3.4641016151377544

@info "Test: perm_cmp()"
z, b = perm_cmp(a1, a2)
@test round.(z) == -2.0 .* ones(2, 3)
@test size(b) == (2, 3)

@info "Test: add_component()"
x = collect(1:10)
e10_tmp = add_component(e10, c=:x, v=x)
@test e10_tmp.components[:x] == x

@info "Test: component_type()"
@test component_type(e10_tmp, c=:x) == Vector{Int64}

@info "Test: extract_component()"
@test extract_component(e10_tmp, c=:x) == x

@info "Test: delete_component()"
e10_tmp2 = deepcopy(e10_tmp)
delete_component!(e10_tmp2, c=:x)
@test e10_tmp2.components == Dict()

@info "Test: rename_component()"
rename_component!(e10_tmp, c_old=:x, c_new=:y)
@test e10_tmp.components[:y] == 1:10

@info "Test: reset_components()"
reset_components!(e10_tmp)
@test e10_tmp.components == Dict()

@info "Test: vsearch()"
@test vsearch(2.1, [1, 2, 3, 4]) == 2
@test vsearch([2.1, 2.9], [1, 2, 3, 4]) == [2.0, 3.0]
@test vsearch("a", ["a", "b", "c"]) == 1
@test vsearch("d", ["a", "b", "c"]) === nothing

@info "Test: trim()"
@test vsplit(1:10, 2) == [[1, 6], [2, 7], [3, 8], [4, 9], [5, 10]]

@info "Test: fft0()"
x = fft0(v1)
@test length(x) == 5
x = fft0(v1, 10)
@test length(x) == 15
x = rfft0(v1)
@test length(x) == 3
x = rfft0(v1, 10)
@test length(x) == div(length(v1) + 10, 2) + 1

@info "Test: ifft0()"
x = fft0(v1, 2)
@test round.(real.(ifft0(x, 2))) == [1.0, 2.0, 3.0, 4.0, 5.0]

@info "Test: fft2()"
x = fft2(v1)
@test length(x) == 8
x = rfft2(v1)
@test length(x) == 5

@info "Test: nextpow()"
@test nextpow2(5) == 8

@info "Test: gradient()"
g, g_l = NeuroAnalyzer.gradient(rand(10))
@test size(g) == (10,)
@test size(g_l) == (10,)
g, g_l = NeuroAnalyzer.gradient(rand(10, 10))
@test size(g) == (10, 10)
@test size(g_l) == (10, 10)
g, g_l = NeuroAnalyzer.gradient(rand(10, 10, 10))
@test size(g) == (10, 10, 10)
@test size(g_l) == (10, 10, 10)

@info "Test: findpeaks()"
x = rand(100)
x[10] *= 1000
@test 10 in findpeaks(x)

@info "Test: hz2rads()"
@test hz2rads(1) == 2pi

@info "Test: rads2hz()"
@test rads2hz(2pi) == 1

@info "Test: t2f()"
@test t2f(1000) == 1.0

@info "Test: f2t()"
@test f2t(1.0) == 1000.0

@info "Test: freqs()"
f, nf = NeuroAnalyzer.freqs(0:1/10:10)
@test length(f) == 51
@test nf == 5
f, nf = NeuroAnalyzer.freqs(0:1/10:10, nf=true)
@test length(f) == 101
@test nf == 5
f, nf = NeuroAnalyzer.freqs(rand(100), 10)
@test length(f) == 51
@test nf == 5
f, nf = NeuroAnalyzer.freqs(rand(100), 10, nf=true)
@test length(f) == 100
@test nf == 5
f, nf = NeuroAnalyzer.freqs(e10)
@test length(f) == 1281
@test nf == 128
f, nf = NeuroAnalyzer.freqs(e10, nf=true)
@test length(f) == 2560
@test nf == 128

@info "Test: generate_window()"
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

@info "Test: generate_sine()"
s = generate_sine(10, 1:100)
@test length(s) == 100

@info "Test: generate_cosine()"
s = generate_cosine(10, 1:100)
@test length(s) == 100

@info "Test: generate_csine()"
s = generate_csine(10, 1:100)
@test length(s) == 100

@info "Test: generate_sinc()"
s = generate_sinc()
@test length(s) == 401

@info "Test: generate_morlet()"
s = generate_morlet(100, 10)
@test length(s) == 201

@info "Test: generate_gaussian()"
s = generate_gaussian(100, 10)
@test length(s) == 201

@info "Test: generate_noise()"
s = generate_noise(100)
@test length(s) == 100

@info "Test: generate_morlet_fwhm()"
s = generate_morlet_fwhm(100, 10)
@test length(s) == 201

@info "Test: sr()"
@test sr(e10) == 256

@info "Test: nchannels()"
@test nchannels(e10) == 24

@info "Test: nepochs()"
@test nepochs(e10) == 10

@info "Test: signal_len()"
@test signal_len(e10) == 25600

@info "Test: epoch_len()"
@test epoch_len(e10) == 2560

@info "Test: get_channel()"
@test length(get_channel(e10, type=["eeg", "eeg", "ecg", "mrk"])) == 20
@test length(get_channel(e10, ch=r"Fp.*")) == 2

@info "Test: cwtfrq()"
s = rand(100)
@test length(cwtfrq(s, fs=10)) == 130
@test length(cwtfrq(e10)) == 19

@info "Test: history()"
@test NeuroAnalyzer.history(e10) isa Vector{String}

@info "Test: labels()"
@test length(labels(e10)) == 24

@info "Test: channel_cluster()"
@test channel_cluster(e10, cluster=:f1) == [1, 3, 11]

@info "Test: band_frq()"
@test band_frq(256, band=:alpha) == (8.0, 13.0)
@test band_frq(e10, band=:alpha) == (8.0, 13.0)

@info "Test: m_pad0()"
@test m_pad0(m1) == [1 2 3; 4 5 6; 0 0 0]
@test size(m_pad0(m1, 6, 5)) == (6, 5)

@info "Test: m_sortperm()"
@test m_sortperm(m1) == [1 1 1; 2 2 2]

@info "Test: m_sort()"
@test m_sort(m1, [2, 1]) == [4 5 6; 1 2 3]

@info "Test: m_norm()"
@test m_norm(m1) == [0.5 1.0 1.5; 2.0 2.5 3.0]

@info "Test: linspace()"
@test linspace(1, 10, 10) == 1:10

@info "Test: log10space()"
@test log10space(1, 2, 2) == [10.0, 100.0]

@info "Test: log2space()"
@test log2space(1, 2, 2) == [2.0, 4.0]

@info "Test: cmax()"
@test cmax([1 + 2im, 10 + 10im]) == 10 + 10im

@info "Test: cmin()"
@test cmin([1 + 2im, 10 + 10im]) == 1 + 2im

@info "Test: tuple_order()"
@test tuple_order((2, 1)) == (1, 2)

@info "Test: cums()"
@test cums(a1) == [1.0 2.0 3.0; 1.0 2.0 3.0;;; 1.0 2.0 3.0; 1.0 2.0 3.0]

@info "Test: f_nearest()"
m = [(1.0, 2.0) (3.0, 4.0); (5.0, 6.0) (7.0, 8.0)]
@test f_nearest(m, (3.2, 3.9)) == (1, 2)

@info "Test: view_note()"
@test view_note(e10) == ""

@info "Test: add_note()"
add_note!(e10, note="test")
@test view_note(e10) == "test"

@info "Test: delete_note()"
delete_note!(e10)
@test view_note(e10) == ""

@info "Test: pad0()"
@test pad0(v1, 5) == [1, 2, 3, 4, 5, 0, 0, 0, 0, 0]
@test pad0(m1, 5) == [1 2 3 0 0 0 0 0; 4 5 6 0 0 0 0 0]

@info "Test: pad2()"
@test pad2(v1) == [1, 2, 3, 4, 5, 0, 0, 0]
@test pad2(m1) == [1 2 3 0; 4 5 6 0]

@info "Test: phases()"
@test length(phases(e10.data[1, :, 1])) == 2560

@info "Test: pick()"
@test channel_pick(e10, p=[:l, :f]) == [1, 3, 11]

@info "Test: t2s()"
@test t2s(1.0, 256) == 256
@test t2s(e10, t=1.0) == 256

@info "Test: s2t()"
@test s2t(2560, 256) == 9.996
@test s2t(e10, s=256) == 0.996

@info "Test: get_channel(wl)"
@test length(get_channel(n, wl=760)) == 36

@info "Test: size()"
@test size(e10) == (24, 2560, 10)

@info "Test: generate_signal()"
s = generate_signal(100)
@test length(s) == 100

@info "Test: chtypes()"
@test length(chtypes(e10)) == 24

@info "Test: optode_labels()"
n = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
@test length(optode_labels(n)) == 3
@test length(source_labels(n)) == 1
@test length(detector_labels(n)) == 2

@info "Test: delmean()"
@test delmean(v1) == [-2, -1, 0, 1, 2]
@test delmean(a1, dims=3) == zeros(2, 3, 2)

@info "Test: tavg()"
@test tavg(a1) == ones(2, 3, 1)

@info "Test: padm()"
@test padm(ones(2), 2) == ones(4)
@test padm(ones(2, 4), 2, mode=:all) == ones(2, 6)
@test padm(ones(2, 4), 2, mode=:row) == ones(2, 6)
@test padm(ones(2, 4, 3), 2, mode=:all) == ones(2, 6, 3)
@test padm(ones(2, 4, 3), 2, mode=:row) == ones(2, 6, 3)

@info "Test: vec2mat()"
@test size(vec2mat(ones(10), wlen=2, woverlap=2)) == (5, 2)

@info "Test: arr2mat()"
@test size(arr2mat(rand(1, 10, 10))) == (10, 10)

@info "Test: minat()"
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
@test minat(v1, v2) == (6, 1)

@info "Test: maxat()"
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
@test maxat(v1, v2) == (2, 5)

@info "Test: paired_labels()"
l = ["ch1", "ch2", "ch3"]
@test length(paired_labels(l; unq=true)) == 6
@test length(paired_labels(l; unq=false)) == 9
@test length(paired_labels(l, l)) == 3

@info "Test: vreduce()"
x = rand(100)
f = linspace(0, 10, 100)
x2, f2 = vreduce(x, f)
@test length(x2) == length(f2) == 21

@info "Test: areduce()"
x = rand(10, 100)
f = linspace(0, 10, 100)
x2, f2 = areduce(x, f)
@test size(x2, 2) == length(f2) == 21
x = rand(10, 100, 5)
f = linspace(0, 10, 100)
x2, f2 = areduce(x, f)
@test size(x2, 2) == length(f2) == 21

@info "Test: info()"
@test isa(info(e10, df=true), DataFrame)

@info "Test: describe()"
@test isa(NeuroAnalyzer.describe(e10, df=true), DataFrame)

@info "Test: ntapers()"
@test ntapers(e10, df=1) == 9

@info "Test: trtm()"
@test size(trtm(e10, ch="Fp1")) == (10, 2560)

@info "Test: meshgrid()"
@test length(meshgrid(collect(range(-1, 1, 100)), collect(range(-1, 1, 100)))) == 2

@info "Test: aff_mni2tal()"
@test aff_mni2tal([10, 12, 14]) == [8.0, 8.32, 12.48]

@info "Test: mni2tal()"
@test mni2tal([10, 12, 14]) == [9.9, 12.2696, 12.2826]

@info "Test: aff_tal2mni()"
@test aff_tal2mni([9.9, 12.2692, 12.2826]) == [12.15909090909091, 16.071340206185567, 13.544355670103092]

@info "Test: tal2mni()"
@test tal2mni([9.9, 12.2692, 12.2821]) == [10.0, 11.999613921643125, 13.999435493742183]

@info "Test: fir_order()"
@test fir_order(0.1, 39) == 18

true