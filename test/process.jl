using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

@info "Initializing"
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
n = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:10)
load_locs!(e10, file_name=joinpath(testfiles_path, "standard-10-20-cap19-elmiko.ced"))
v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1 2 3; 4 5 6]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "test 1/44: add_signal()"
@test add_signal(v1, v2) == v1 + v2
x = rand(epoch_len(e10))
e10_tmp = add_signal(e10, s=x)
e10_tmp.data[1, :, 1] == e10.data[1, :, 1] + x

@info "test 2/44: average()"
@test average(a1) == ones(1, 3, 2)
@test average(a1, a2) == 0.5 .* ones(2, 1, 2)
e10_tmp = average(e10)
@test size(e10_tmp.data) == (1, 2560, 10)
e10_tmp = average(e10, e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 3/44: cbp()"
@test length(cbp(rand(100), fs=10, frq=4)) == 100
e10_tmp = cbp(e10, frq=4)
@test size(e10_tmp.data) == size(e10.data)

@info "test 4/44: ch_zero()"
e10_tmp = ch_zero(e10)
@test e10_tmp.data[1, 1, 1] == 0
@test e10_tmp.data[1, end, 1] == 0

@info "test 5/44: cw_trans()"
s = rand(100)
ct = cw_trans(s, wt=wavelet(Morlet(2π), β=2))
@test size(ct) == (12, 100)
ct = cw_trans(e10, wt=wavelet(Morlet(2π), β=2));
@test size(ct) == (19, 30, 2560, 10)

@info "test 6/44: icw_trans()"
ct = cw_trans(s, wt=wavelet(Morlet(2π), β=2))
s_new = icw_trans(ct, wt=wavelet(Morlet(2π), β=2), type=:nd)
@test length(s_new) == 100
s_new = icw_trans(ct, wt=wavelet(Morlet(2π), β=2), type=:pd)
@test length(s_new) == 100
s_new = icw_trans(ct, wt=wavelet(Morlet(2π), β=2), type=:df)
@test length(s_new) == 100

@info "test 7/44: denoise_fft()"
s = rand(100)
s2, f = denoise_fft(s)
@test length(s2) == 100
@test length(f) == 100
e10_tmp = denoise_fft(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 8/44: denoise_wavelet())"
s = denoise_wavelet(rand(100), wt=wavelet(WT.haar))
@test length(s) == 100
e10_tmp = denoise_wavelet(e10, wt=wavelet(WT.haar))
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 9/44: denoise_wien()"
s = denoise_wien(a1)
@test size(s) == (2, 3, 2)
e10_tmp = denoise_wien(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 10/44: derivative()"
@test NeuroAnalyzer.derivative(v1) == ones(5)
@test NeuroAnalyzer.derivative(a1) == zeros(2, 3, 2)
e10_tmp = NeuroAnalyzer.derivative(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 11/44: detrend()"
@test detrend(v1) == zeros(5)
e10_tmp = detrend(e10, type=:ls)
@test size(e10_tmp.data) == (24, 2560, 10)
e10_tmp = detrend(e10, type=:linear)
@test size(e10_tmp.data) == (24, 2560, 10)
e10_tmp = detrend(e10, type=:constant)
@test size(e10_tmp.data) == (24, 2560, 10)
e10_tmp = detrend(e10, type=:poly)
@test size(e10_tmp.data) == (24, 2560, 10)
e10_tmp = detrend(e10, type=:loess)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 12/44: dw_trans()"
s = rand(100)
dt = dw_trans(s, wt=wavelet(WT.haar), type=:sdwt)
@test size(dt) == (3, 100)
dt = dw_trans(s, wt=wavelet(WT.haar), type=:acdwt)
@test size(dt) == (3, 100)
dt = dw_trans(e10, wt=wavelet(WT.haar), type=:sdwt)
@test size(dt) == (19, 10, 2560, 10)
dt = dw_trans(e10, wt=wavelet(WT.haar), type=:acdwt)
@test size(dt) == (19, 10, 2560, 10)

@info "test 13/44: idw_trans()"
dt = dw_trans(s, wt=wavelet(WT.haar), type=:sdwt)
s_new = idw_trans(dt, wt=wavelet(WT.haar), type=:sdwt)
@test length(s_new) == 100
s_new = idw_trans(dt, wt=wavelet(WT.haar), type=:acdwt)
@test length(s_new) == 100

@info "test 14/44: dwtsplit()"
s = dwtsplit(e10, ch=1, wt = wavelet(WT.haar), type=:sdwt)
@test size(s) == (10, 2560, 10)

@info "test 15/44: erp()"
e = erp(e10)
@test size(e.data) == (19, 2560, 11)

@info "test 16/44: bpsplit()"
s, bn, bf = bpsplit(e10)
@test length(bn) == 10
@test length(bf) == 10
@test size(s) == (10, 19, 2560, 10)

@info "test 17/44: fconv()"
@test fconv(v1, kernel=v2) == [2.0, 2.0, 3.0, 3.0, 2.0]
@test fconv(a1, kernel=[0.5, 1.0, 0.5]) == ones(2, 3, 2)
e10_tmp = fconv(e10, kernel=[0.0, 0.5, 1.0, 0.5, 0.0])
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 18/44: filter_mavg()"
@test filter_mavg(v1, k=2) == [0.0, 0.0, 3.0, 0.0, 0.0]
e10_tmp = filter_mavg(e10, k=2)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 19/44: filter_mmed()"
@test filter_mmed(v1, k=2) == [0.0, 0.0, 3.0, 0.0, 0.0]
e10_tmp = filter_mmed(e10, k=2)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 20/44: filter_poly()"
s = filter_poly(rand(20))
@test length(s) == 20
e10_tmp = filter_poly(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 21/44: filter_sg()"
s = filter_sg(rand(20))
@test length(s) == 20
e10_tmp = filter_poly(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 22/44: filter()"
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:lp, cutoff=40, order=8)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:hp, cutoff=1, order=12)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:bs, cutoff=(49, 51), order=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:bp, cutoff=(49, 51), order=4)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev1, ftype=:lp, cutoff=40, order=8, rs=35)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev1, ftype=:hp, cutoff=1, order=12, rs=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev1, ftype=:bs, cutoff=(49, 51), order=2, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev1, ftype=:bp, cutoff=(49, 51), order=4, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev2, ftype=:lp, cutoff=40, order=8, rs=35)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev2, ftype=:hp, cutoff=1, order=12, rs=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev2, ftype=:bs, cutoff=(49, 51), order=2, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:chebyshev2, ftype=:bp, cutoff=(49, 51), order=4, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:elliptic, ftype=:lp, cutoff=40, order=8, rs=35)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:elliptic, ftype=:hp, cutoff=1, order=12, rs=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:elliptic, ftype=:bs, cutoff=(49, 51), order=2, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:elliptic, ftype=:bp, cutoff=(49, 51), order=4, rs=50)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:lp, cutoff=40, order=8)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:hp, cutoff=1, order=12)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:bs, cutoff=(49, 51), order=2)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:bp, cutoff=(49, 51), order=4)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:lp, cutoff=40, order=8, bw=10)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:hp, cutoff=1, order=12, bw=0.5)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:bs, cutoff=(49, 51), order=4, bw=0.1)
@test size(eeg_tmp.data) == (24, 2560, 10)
eeg_tmp = NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:bp, cutoff=(49, 51), order=4, bw=0.5)
@test size(eeg_tmp.data) == (24, 2560, 10)

@info "test 23/44: filter_g()"
s = filter_g(rand(20), fs=2, f=4)
@test length(s) == 20
e10_tmp = filter_g(e10, f=20)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 24/44: invert_polarity()"
e10_tmp = invert_polarity(e10)
@test e10_tmp.data == .-(e10.data)

@info "test 25/44: lrinterpolate_channel()"
e10_tmp = deepcopy(e10)
e10_tmp.data[1, :, 1] = zeros(epoch_len(e10))
e10_int = lrinterpolate_channel(e10_tmp, ch=1, ep=1)
@test e10_int.data[1, :, 1] != e10_tmp.data[1, :, 1]

@info "test 26/44: normalize()"
@test normalize(v1, method=:none) == v1
@test normalize_zscore(v1) == [-1.2649110640673518, -0.6324555320336759, 0.0, 0.6324555320336759, 1.2649110640673518]
@test normalize_minmax(v1) == [-1.0, -0.5, 0.0, 0.5, 1.0]
@test normalize_n(v1) == [0.0, 0.25, 0.5, 0.75, 1.0]
@test normalize_log(v1) == [1.0986122886681098, 1.3862943611198906, 1.6094379124341003, 1.791759469228055, 1.9459101490553132]
@test normalize_gauss(v1) == [-0.8047189562170504, -0.3465735902799727, 0.0, 0.3465735902799726, 0.8047189562170504]
@test normalize_log10(v1) == [0.47712125471966244, 0.6020599913279624, 0.6989700043360189, 0.7781512503836436, 0.8450980400142568]
@test normalize_neglog(v1) == [-0.0, -0.6931471805599453, -1.0986122886681098, -1.3862943611198906, -1.6094379124341003]
@test normalize_neglog10(v1) == [-0.0, -0.3010299956639812, -0.47712125471966244, -0.6020599913279624, -0.6989700043360189]
@test normalize_neg(v1) == [-4, -3, -2, -1, 0]
@test normalize_pos(v1) == [2, 3, 4, 5, 6]
@test normalize_perc(v1) == [0.0, 0.25, 0.5, 0.75, 1.0]
@test normalize_invroot(v1) == [0.7071067811865475, 0.5773502691896258, 0.5, 0.4472135954999579, 0.4082482904638631]
e10_tmp = normalize(e10, method=:zscore)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 27/44: lrinterpolate_channel()"
e10_tmp = deepcopy(e10)
e10_tmp.data[1, :, 1] = zeros(epoch_len(e10))
e10_int = plinterpolate_channel(e10_tmp, ch=1, ep=1)
@test e10_int.data[1, :, 1] != e10_tmp.data[1, :, 1]

@info "test 28/44: remove_dc()"
@test remove_dc(v1) == [-2.0, -1.0, 0.0, 1.0, 2.0]
e10_tmp = remove_dc(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 29/44: scale()"
e10_tmp = scale(e10, factor=2.0)
@test e10_tmp.data == e10.data .* 2.0

@info "test 30/44: reference()"
e10_tmp = reference_ch(e10, ch=1)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 31/44: reference_a()"
e10_tmp = reference_a(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 32/44: reference_m()"
edit_channel!(e10, ch=20, field=:labels, value="M1")
edit_channel!(e10, ch=21, field=:labels, value="M2")
e10_tmp = reference_m(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 33/44: reference_car()"
e10_tmp = reference_car(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 34/44: reference_plap()"
e10_tmp = reference_plap(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 35/44: csd()"
g, h = gh(e10.locs)
@test size(g) == (19, 19)
@test size(h) == (19, 19)
e10_tmp = csd(e10)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 36/44: standardize()"
m_s, sc = standardize(a1)
@test length(sc) == 2
m_s, sc = standardize(e10)
@test length(sc) == 10

@info "test 37/44: taper()"
@test taper(v1, t=v1) == [1, 4, 9, 16, 25]
e10_tmp = taper(e10, t=e10.data[1, :, 1])
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 38/44: tconv()"
@test tconv(v1, kernel=[0.2, 0.1, 0.2]) == [0.49999999999999983, 1.0000000000000002, 1.5000000000000002, 2.0, 1.2999999999999998]
@test tconv(a1, kernel=[0.2, 0.1, 0.2]) == [0.30000000000000004 0.5 0.30000000000000004; 0.30000000000000004 0.5 0.30000000000000004;;; 0.30000000000000004 0.5 0.30000000000000004; 0.30000000000000004 0.5 0.30000000000000004]
s = tconv(e10, kernel=[0.2, 0.1, 0.2])
@test size(s) == (19, 2560, 10)

@info "test 39/44: wbp()"
@test length(wbp(e10.data[1, :, 1], fs=10, frq=4)) == 2560
e10_tmp = wbp(e10, frq=4)
@test size(e10_tmp.data) == (24, 2560, 10)

@info "test 40/44: intensity2od()"
n_tmp = intensity2od(n)
@test size(n_tmp.data) == (16, 9015, 1)

@info "test 41/44: od2conc()"
n_tmp2 = od2conc(n_tmp)
@test size(n_tmp2.data) == (22, 9015, 1)

@info "test 42/44: npl()"
e = erp(e10)
npl!(e)
@test size(e.data) == (19, 2560, 11)

@info "test 43/44: remove_pops()"
eeg_tmp, pl, ls, rs = remove_pops(eeg)
@test size(eeg_tmp.data) == size(eeg.data)

@info "test 44/44: remove_powerline()"
e10_tmp = remove_powerline(e10)
@test size(e10_tmp.data) == size(e10.data)


true