using NeuroAnalyzer
using Test
using Wavelets
using ContinuousWavelets

ntests = 54

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

@info "Test 1/$ntests: add_signal()"
@test add_signal(v1, v2) == v1 + v2
x = rand(epoch_len(e10))
e10_tmp = add_signal(e10, ch="Fp1", s=x)
e10_tmp.data[1, :, 1] == e10.data[1, :, 1] + x

@info "Test 2/$ntests: average()"
@test average(a1) == ones(1, 3, 2)
@test average(a1, a2) == 0.5 .* ones(2, 1, 2)
e10_tmp = average(e10, ch="Fp1")
@test size(e10_tmp) == (1, 2560, 10)
e10_tmp = average(e10, e10)
@test size(e10_tmp) == size(e10)

@info "Test 3/$ntests: cbp()"
@test length(cbp(rand(100), fs=10, frq=4)) == 100
e10_tmp = cbp(e10, ch="all", frq=4)
@test size(e10_tmp) == size(e10)

@info "Test 4/$ntests: ch_zero()"
e10_tmp = ch_zero(e10)
@test e10_tmp.data[1, 1, 1] == 0
@test e10_tmp.data[1, end, 1] == 0

@info "Test 5/$ntests: cw_trans()"
s = rand(100)
ct = cw_trans(s, wt=wavelet(Morlet(π), β=32, Q=128))
@test size(ct) == (130, 100)
ct = cw_trans(e10, ch="Fp1", wt=wavelet(Morlet(π), β=32, Q=128));
@test size(ct) == (1, 131, 2560, 10)

@info "Test 6/$ntests: icw_trans()"
ct = cw_trans(s, wt=wavelet(Morlet(π), β=32, Q=128))
s_new = icw_trans(ct, wt=wavelet(Morlet(π), β=32, Q=128), type=:nd)
@test length(s_new) == 100
s_new = icw_trans(ct, wt=wavelet(Morlet(π), β=32, Q=128), type=:pd)
@test length(s_new) == 100
s_new = icw_trans(ct, wt=wavelet(Morlet(π), β=32, Q=128), type=:df)
@test length(s_new) == 100

@info "Test 7/$ntests: denoise_fft()"
s = rand(100)
s2, f = denoise_fft(s)
@test length(s2) == 100
@test length(f) == 100
e10_tmp = denoise_fft(e10, ch="all")
@test size(e10_tmp) == size(e10)

@info "Test 8/$ntests: denoise_dwt())"
s = denoise_dwt(rand(100), wt=wavelet(WT.haar))
@test length(s) == 100
e10_tmp = denoise_dwt(e10, ch="all", wt=wavelet(WT.haar))
@test size(e10_tmp) == size(e10)

@info "Test 9/$ntests: denoise_wien()"
s = denoise_wien(a1)
@test size(s) == (2, 3, 2)
e10_tmp = denoise_wien(e10, ch="all")
@test size(e10_tmp) == size(e10)

@info "Test 10/$ntests: derivative()"
@test NeuroAnalyzer.derivative(v1) == [1, 1, 1, 1, 1]
@test NeuroAnalyzer.derivative(a1) == [0.0 0.0 0.0; 0.0 0.0 0.0;;; 0.0 0.0 0.0; 0.0 0.0 0.0]
e10_tmp = NeuroAnalyzer.derivative(e10, ch="all")
@test size(e10_tmp) == size(e10)

@info "Test 11/$ntests: detrend()"
@test round.(detrend(v1)) == zeros(5)
e10_tmp = detrend(e10, ch="all", type=:ls)
@test size(e10_tmp) == size(e10)
e10_tmp = detrend(e10, ch="all", type=:linear)
@test size(e10_tmp) == size(e10)
e10_tmp = detrend(e10, ch="all", type=:constant)
@test size(e10_tmp) == size(e10)
e10_tmp = detrend(e10, ch="all", type=:poly)
@test size(e10_tmp) == size(e10)
e10_tmp = detrend(e10, ch="all", type=:loess)
@test size(e10_tmp) == size(e10)

@info "Test 12/$ntests: dw_trans()"
s = rand(100)
dt = dw_trans(s, wt=wavelet(WT.haar), type=:sdwt)
@test size(dt) == (3, 100)
dt = dw_trans(s, wt=wavelet(WT.haar), type=:acdwt)
@test size(dt) == (3, 100)
dt = dw_trans(e10, ch="all", wt=wavelet(WT.haar), type=:sdwt)
@test size(dt) == (24, 10, 2560, 10)
dt = dw_trans(e10, ch="all", wt=wavelet(WT.haar), type=:acdwt)
@test size(dt) == (24, 10, 2560, 10)

@info "Test 13/$ntests: idw_trans()"
dt = dw_trans(s, wt=wavelet(WT.haar), type=:sdwt)
s_new = idw_trans(dt, wt=wavelet(WT.haar), type=:sdwt)
@test length(s_new) == 100
s_new = idw_trans(dt, wt=wavelet(WT.haar), type=:acdwt)
@test length(s_new) == 100

@info "Test 14/$ntests: dw_split()"
s = dw_split(e10, ch="Fp1", wt = wavelet(WT.haar), type=:sdwt)
@test size(s) == (10, 2560, 10)

@info "Test 15/$ntests: erp()"
e = erp(e10)
@test size(e.data) == (19, 2560, 11)

@info "Test 16/$ntests: bpsplit()"
s, bn, bf = bpsplit(e10, ch="all")
@test length(bn) == 13
@test length(bf) == 13
@test size(s) == (13, 24, 2560, 10)

@info "Test 17/$ntests: fconv()"
@test fconv(v1, kernel=v2) == [0.8500000000000002 - 9.868649107779169e-17im, 1.5999999999999996 + 0.0im, 2.4999999999999996 + 0.0im, 3.499999999999999 - 3.61217627448181e-17im, 2.6999999999999997 - 3.204937810639273e-17im]
@test fconv(a1, kernel=[0.5, 1.0, 0.5]) == [0.25 + 0.0im 0.75 + 0.0im 1.0 + 0.0im; 0.25 + 0.0im 0.75 + 0.0im 1.0 + 0.0im;;; 0.25 + 0.0im 0.75 + 0.0im 1.0 + 0.0im; 0.25 + 0.0im 0.75 + 0.0im 1.0 + 0.0im]
s_conv = fconv(e10, ch="all", kernel=[0.0, 0.5, 1.0, 0.5, 0.0])
@test size(s_conv) == size(e10)

@info "Test 18/$ntests: filter_mavg()"
@test filter_mavg(vcat(v1, v1), k=2) == [1, 2, 3, 3, 3, 3, 3, 3, 4, 5]
e10_tmp = filter_mavg(e10, ch="all", k=2)
@test size(e10_tmp) == size(e10)

@info "Test 19/$ntests: filter_mmed()"
@test filter_mmed(vcat(v1, v1), k=2) == [1, 2, 3, 3, 3, 3, 3, 3, 4, 5]
e10_tmp = filter_mmed(e10, ch="all", k=2)
@test size(e10_tmp) == size(e10)

@info "Test 20/$ntests: filter_poly()"
s = filter_poly(rand(20))
@test length(s) == 20
e10_tmp = filter_poly(e10, ch="all")
@test size(e10_tmp) == size(e10)

@info "Test 21/$ntests: filter_sg()"
s = filter_sg(rand(20))
@test length(s) == 20
e10_tmp = filter_poly(e10, ch="all")
@test size(e10_tmp) == size(e10)

@info "Test 22/$ntests: filter()"
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:butterworth, ftype=:lp, cutoff=40, order=8)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:butterworth, ftype=:hp, cutoff=1, order=12)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:butterworth, ftype=:bs, cutoff=(49, 51), order=2)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:butterworth, ftype=:bp, cutoff=(49, 51), order=4)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:chebyshev1, ftype=:lp, cutoff=40, order=8, rs=35)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:chebyshev1, ftype=:hp, cutoff=1, order=12, rs=2)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:chebyshev1, ftype=:bs, cutoff=(49, 51), order=2, rs=50)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:chebyshev1, ftype=:bp, cutoff=(49, 51), order=4, rs=50)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:chebyshev2, ftype=:lp, cutoff=40, order=8, rs=35)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:chebyshev2, ftype=:hp, cutoff=1, order=12, rs=2)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:chebyshev2, ftype=:bs, cutoff=(49, 51), order=2, rs=50)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:chebyshev2, ftype=:bp, cutoff=(49, 51), order=4, rs=50)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:elliptic, ftype=:lp, cutoff=40, order=8, rs=35)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:elliptic, ftype=:hp, cutoff=1, order=12, rs=2)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:elliptic, ftype=:bs, cutoff=(49, 51), order=2, rs=50)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:elliptic, ftype=:bp, cutoff=(49, 51), order=4, rs=50)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:fir, ftype=:lp, cutoff=40, order=8)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:fir, ftype=:hp, cutoff=1, order=12)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:fir, ftype=:bs, cutoff=(49, 51), order=2)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:fir, ftype=:bp, cutoff=(49, 51), order=4)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:firls, ftype=:lp, cutoff=40)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:firls, ftype=:hp, cutoff=1)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:firls, ftype=:bs, cutoff=(49, 51))
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:firls, ftype=:bp, cutoff=(49, 51))
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:remez, ftype=:lp, cutoff=40, order=8, bw=10)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:remez, ftype=:hp, cutoff=1, order=12, bw=0.5)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:remez, ftype=:bs, cutoff=(49, 51), order=4, bw=0.1)
@test size(eeg_tmp) == size(e10)
eeg_tmp = NeuroAnalyzer.filter(e10, ch="all", fprototype=:remez, ftype=:bp, cutoff=(49, 51), order=4, bw=0.5)
@test size(eeg_tmp) == size(e10)

@info "Test 23/$ntests: filter_g()"
s = filter_g(rand(20), fs=2, f=4)
@test length(s) == 20
e10_tmp = filter_g(e10, ch="all", f=20)
@test size(e10_tmp) == size(e10)

@info "Test 24/$ntests: invert_polarity()"
e10_tmp = invert_polarity(e10, ch="all")
@test e10_tmp.data == .-(e10.data)

@info "Test 25/$ntests: lrinterpolate_channel()"
e10_tmp = deepcopy(e10)
e10_tmp.data[1, :, 1] = zeros(epoch_len(e10))
e10_int = lrinterpolate_channel(e10_tmp, ch="Fp1", ep=1)
@test e10_int.data[1, :, 1] != e10_tmp.data[1, :, 1]

@info "Test 26/$ntests: normalize()"
@test NeuroAnalyzer.normalize(v1, method=:none) == v1
@test NeuroAnalyzer.normalize(m1, method=:none) == m1
@test NeuroAnalyzer.normalize(a1, method=:none) == a1
@test normalize_zscore(v1) == [-1.2649110640673518, -0.6324555320336759, 0.0, 0.6324555320336759, 1.2649110640673518]
@test normalize_zscore(m1, bych=true) == [-1 0 1; -1 0 1]
@test normalize_zscore(m1, bych=false) == [-1.3363062095621219 -0.8017837257372732 -0.2672612419124244; 0.2672612419124244 0.8017837257372732 1.3363062095621219]
@test normalize_minmax(v1) == [-1.0, -0.5, 0.0, 0.5, 1.0]
@test normalize_minmax(m1, bych=true) == [-1 0 1; -1 0 1]
@test normalize_minmax(m1, bych=false) == [-1.0 -0.6 -0.19999999999999996; 0.19999999999999996 0.6000000000000001 1.0]
@test normalize_n(v1) == [0.0, 0.25, 0.5, 0.75, 1.0]
@test normalize_n(m1, bych=true) == [0.0 0.5 1.0; 0.0 0.5 1.0]
@test normalize_n(m1, bych=false) == [0.0 0.2 0.4; 0.6 0.8 1.0]
@test normalize_log(v1) == [1.0986122886681098, 1.3862943611198906, 1.6094379124341003, 1.791759469228055, 1.9459101490553132]
@test normalize_log(m1, bych=true) == [1.0986122886681098 1.3862943611198906 1.6094379124341003; 2.1972245773362196 2.302585092994046 2.3978952727983707]
@test normalize_log(m1, bych=false) == [1.0986122886681098 1.3862943611198906 1.6094379124341003; 1.791759469228055 1.9459101490553132 2.0794415416798357]
@test normalize_gauss(v1) == [-0.8047189562170504, -0.3465735902799727, 0.0, 0.3465735902799726, 0.8047189562170504]
@test normalize_gauss(m1, bych=true) == [-0.5493061443340549 0.0 0.5493061443340549; -0.5493061443340549 0.0 0.5493061443340549]
@test normalize_gauss(m1, bych=false) == [-0.8958797346140276 -0.45814536593707755 -0.1438410362258905; 0.1438410362258904 0.45814536593707755 0.8958797346140273]
@test normalize_log10(v1) == [0.47712125471966244, 0.6020599913279624, 0.6989700043360189, 0.7781512503836436, 0.8450980400142568]
@test normalize_log10(m1, bych=true) == [0.47712125471966244 0.6020599913279624 0.6989700043360189; 0.9542425094393249 1.0 1.0413926851582251]
@test normalize_log10(m1, bych=false) == [0.47712125471966244 0.6020599913279624 0.6989700043360189; 0.7781512503836436 0.8450980400142568 0.9030899869919435]
@test normalize_neglog(v1) == [-0.0, -0.6931471805599453, -1.0986122886681098, -1.3862943611198906, -1.6094379124341003]
@test normalize_neglog(m1, bych=true) == [-0.0 -0.6931471805599453 -1.0986122886681098; -1.3862943611198906 -1.6094379124341003 -1.791759469228055]
@test normalize_neglog(m1, bych=false) == [-0.0 -0.6931471805599453 -1.0986122886681098; -1.3862943611198906 -1.6094379124341003 -1.791759469228055]
@test normalize_neglog10(v1) == [-0.0, -0.3010299956639812, -0.47712125471966244, -0.6020599913279624, -0.6989700043360189]
@test normalize_neglog10(m1, bych=true) == [-0.0 -0.3010299956639812 -0.47712125471966244; -0.6020599913279624 -0.6989700043360189 -0.7781512503836436]
@test normalize_neglog10(m1, bych=false) == [-0.0 -0.3010299956639812 -0.47712125471966244; -0.6020599913279624 -0.6989700043360189 -0.7781512503836436]
@test normalize_neg(v1) == [-4, -3, -2, -1, 0]
@test normalize_neg(m1, bych=true) == [-2.0 -1.0 0.0; -2.0 -1.0 0.0]
@test normalize_neg(m1, bych=false) == [-5 -4 -3; -2 -1 0]
@test normalize_pos(v1) == [2, 3, 4, 5, 6]
@test normalize_pos(m1, bych=true) == [2.0 3.0 4.0; 8.0 9.0 10.0]
@test normalize_pos(m1, bych=false) == [2 3 4; 5 6 7]
@test normalize_perc(v1) == [0.0, 0.25, 0.5, 0.75, 1.0]
@test normalize_perc(m1, bych=true) == [0.0 0.5 1.0; 0.0 0.5 1.0]
@test normalize_perc(m1, bych=false) == [0.0 0.2 0.4; 0.6 0.8 1.0]
@test normalize_invroot(v1) == [1.0, 0.7071067811865475, 0.5773502691896258, 0.5, 0.4472135954999579]
@test normalize_invroot(m1, bych=true) == [1.0 0.7071067811865475 0.5773502691896258; 0.5 0.4472135954999579 0.4082482904638631]
@test normalize_invroot(m1, bych=false) == [1.0 0.7071067811865475 0.5773502691896258; 0.5 0.4472135954999579 0.4082482904638631]
@test normalize_softmax(v1) == [0.011656230956039607, 0.03168492079612427, 0.0861285444362687, 0.23412165725273662, 0.6364086465588308]
@test normalize_softmax(m1, bych=true) == [0.0042697785452821095 0.011606461431184656 0.03154963320110001; 0.08576079462509834 0.23312200962361299 0.6336913225737218]
@test normalize_softmax(m1, bych=false) == [0.0042697785452821095 0.011606461431184656 0.03154963320110001; 0.08576079462509834 0.23312200962361299 0.6336913225737218]
@test normalize_sigmoid(v1) == [0.7310585786300049, 0.8807970779778823, 0.9525741268224334, 0.9820137900379085, 0.9933071490757153]
@test normalize_sigmoid(m1, bych=true) == [0.7310585786300049 0.8807970779778823 0.9525741268224334; 0.9820137900379085 0.9933071490757153 0.9975273768433653]
@test normalize_sigmoid(m1, bych=false) == [0.7310585786300049 0.8807970779778823 0.9525741268224334; 0.9820137900379085 0.9933071490757153 0.9975273768433653]
@test normalize_mad(v1) == [-0.9098742077378683, -0.45493710386893416, 0.0, 0.45493710386893416, 0.9098742077378683]
@test normalize_mad(m1, bych=true) == [-0.45493710386893416 0.0 0.45493710386893416; -0.45493710386893416 0.0 0.45493710386893416]
@test normalize_mad(m1, bych=false) == [-1.1241495836601363 -0.6744897501960818 -0.22482991673202726; 0.22482991673202726 0.6744897501960818 1.1241495836601363]
@test normalize_rank(v1) == 1.0:5.0
@test normalize_rank(m1, bych=true) == [1.0 2.0 3.0; 1.0 2.0 3.0]
@test normalize_rank(m1, bych=false) == [1.0 2.0 3.0; 4.0 5.0 6.0]
@test normalize_fisher(v1) == [-18.36840028483855, -0.5493061443340549, 0.0, 0.5493061443340549, 18.36840028483855]
@test normalize_fisher(m1, bych=true) == [-18.36840028483855 -0.6931471805599453 -0.20273255405408214; 0.2027325540540821 0.6931471805599454 18.36840028483855]
@test normalize_fisher(m1, bych=false) == [-18.36840028483855 -0.6931471805599453 -0.20273255405408214; 0.2027325540540821 0.6931471805599454 18.36840028483855]
e10_tmp = NeuroAnalyzer.normalize(e10, ch="all", method=:zscore)
@test size(e10_tmp) == size(e10)

@info "Test 27/$ntests: plinterpolate_channel()"
e10_tmp = deepcopy(e10)
e10_tmp.data[1, :, 1] = zeros(epoch_len(e10))
e10_int = plinterpolate_channel(e10_tmp, ch="Fp1", ep=1)
@test e10_int.data[1, :, 1] != e10_tmp.data[1, :, 1]

@info "Test 28/$ntests: remove_dc()"
@test remove_dc(v1) == [-2.0, -1.0, 0.0, 1.0, 2.0]
e10_tmp = remove_dc(e10, ch="all")
@test size(e10_tmp) == size(e10)

@info "Test 29/$ntests: scale()"
e10_tmp = NeuroAnalyzer.scale(e10, ch="all", factor=2.0)
@test e10_tmp.data == e10.data .* 2.0

@info "Test 30/$ntests: reference()"
e10_tmp = reference_ce(e10, ch="Fp1")
@test size(e10_tmp) == size(e10)
e10_tmp = reference_ce(e10, ch=labels(e10)[1:5])
@test size(e10_tmp) == size(e10)
e10_tmp = reference_ce(e10, ch="Fp1", med=true)
@test size(e10_tmp) == size(e10)

@info "Test 31/$ntests: reference_a()"
e10_tmp = reference_a(e10)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_a(e10, med=true)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_a(e10, type=:c)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_a(e10, type=:i)
@test size(e10_tmp) == size(e10)

@info "Test 32/$ntests: reference_m()"
edit_channel!(e10, ch="A1", field=:label, value="M1")
edit_channel!(e10, ch="A2", field=:label, value="M2")
e10_tmp = reference_m(e10)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_m(e10, med=true)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_m(e10, type=:c)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_m(e10, type=:i)
@test size(e10_tmp) == size(e10)

@info "Test 33/$ntests: reference_avg()"
e10_tmp = reference_avg(e10)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_avg(e10, exclude_fpo=true)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_avg(e10, exclude_current=true)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_avg(e10, average=false)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_avg(e10, weighted=true)
@test size(e10_tmp) == size(e10)

@info "Test 34/$ntests: reference_plap()"
e10_tmp = reference_plap(e10)
@test size(e10_tmp) == size(e10)
e10_tmp = reference_plap(e10, weighted=true)
@test size(e10_tmp) == size(e10)

@info "Test 35/$ntests: csd()"
g, h = gh(e10.locs)
@test size(g) == (23, 23)
@test size(h) == (23, 23)
e10_tmp = csd(e10)
@test size(e10_tmp) == size(e10)

@info "Test 36/$ntests: standardize()"
m_s, sc = NeuroAnalyzer.standardize(a1)
@test length(sc) == 2
m_s, sc = NeuroAnalyzer.standardize(e10, ch="all")
@test length(sc) == 10

@info "Test 37/$ntests: taper()"
@test taper(v1, t=v1) == [1, 4, 9, 16, 25]
e10_tmp = taper(e10, ch="all", t=e10.data[1, :, 1])
@test size(e10_tmp) == size(e10)

@info "Test 38/$ntests: tconv()"
@test round.(tconv(v1, kernel=[0.2, 0.1, 0.2]), digits=2) == [0.2, 0.5, 1.0, 1.5, 2.0]
@test round.(tconv(a1, kernel=[0.2, 0.1, 0.2]), digits=2) == [0.2 0.3 0.5; 0.2 0.3 0.5;;; 0.2 0.3 0.5; 0.2 0.3 0.5]
e10_tmp = tconv(e10, ch="all", kernel=[0.2, 0.1, 0.2])
@test size(e10_tmp) == size(e10)

@info "Test 39/$ntests: wbp()"
@test length(wbp(e10.data[1, :, 1], fs=10, frq=4)) == 2560
e10_tmp = wbp(e10, ch="all", frq=4)
@test size(e10_tmp) == size(e10)

@info "Test 40/$ntests: intensity2od()"
n_tmp = intensity2od(n)
@test size(n_tmp.data) == (16, 9015, 1)

@info "Test 41/$ntests: od2conc()"
n_tmp2 = od2conc(n_tmp)
@test size(n_tmp2.data) == (22, 9015, 1)

@info "Test 42/$ntests: npl()"
e = erp(e10)
npl!(e)
@test size(e.data) == (19, 2560, 11)

@info "Test 43/$ntests: remove_pops()"
eeg_tmp, pl, ls, rs = remove_pops(eeg, ch="all")
@test size(eeg_tmp) == size(eeg.data)

@info "Test 44/$ntests: remove_powerline()"
e10_tmp = keep_epoch(e10, ep=1)
remove_powerline!(e10_tmp, pl_frq=50, ch="Fp1")
@test size(e10_tmp) == (24, 2560, 1)

@info "Test 45/$ntests: ica_decompose()"
ic, ic_mw = ica_decompose(rand(10, 1000), n=5)
@test size(ic) == (5, 1000)
@test size(ic_mw) == (10, 5)
ic, ic_mw, ic_var = ica_decompose(eeg, ch="all", n=5, iter=10)
@test size(ic) == (5, 308480)
@test size(ic_mw) == (24, 5)
@test length(ic_var) == 5

@info "Test 46/$ntests: pca_decompose()"
pc, pcv, pcm, pc_model = pca_decompose(rand(4, 4, 2), n=2)
@test size(pc) == (2, 4, 2)
@test size(pcv) == (2, 2)
@test length(pcm) == 4
pc, pcv, pcm, _ = pca_decompose(e10, ch="all", n=4)
@test size(pc) == (2, 2560, 10)
@test size(pcv) == (2, 10)
@test length(pcm) == 24

@info "Test 47/$ntests: pca_reconstruct()"
pc, pcv, pcm, pc_model = pca_decompose(rand(4, 4, 2), n=2)
s = pca_reconstruct(rand(4, 4, 2); pc=pc, pc_model=pc_model)
@test size(s) == (4, 4, 2)
pc, pcv, pcm, pc_model = pca_decompose(e10, ch="all", n=4)
e10_tmp = add_component(e10, c=:pc, v=pc)
add_component!(e10_tmp, c=:pc_model, v=pc_model)
e10_rec = pca_reconstruct(e10_tmp, ch="all");
@test size(e10_rec.data) == size(e10)
e10_rec = pca_reconstruct(e10_tmp, ch="all", pc, pc_model);
@test size(e10_rec.data) == size(e10)

@info "Test 48/$ntests: reference_custom()"
e10_tmp = reference_custom(e10)
@test size(e10_tmp) == (23, 2560, 10)

@info "Test 49/$ntests: ica_reconstruct()"
ic, ic_mw = ica_decompose(rand(10, 1000), n=5)
s = ica_reconstruct(ic=ic, ic_mw=ic_mw, ic_idx=5)
@test size(s) == (10, 1000)
ic, ic_mw = ica_decompose(eeg, ch="all", n=5, iter=10)
eeg_tmp = ica_reconstruct(eeg, ch="all", ic, ic_mw, ic_idx=1)
@test size(eeg_tmp) == size(eeg)
eeg_tmp = deepcopy(eeg)
add_component!(eeg_tmp, c=:ic, v=ic)
add_component!(eeg_tmp, c=:ic_mw, v=ic_mw)
eeg_tmp = ica_reconstruct(eeg_tmp, ch="all", ic_idx=1);
@test size(eeg_tmp) == size(eeg)

@info "Test 50/$ntests: ica_remove()"
ic, ic_mw = ica_decompose(eeg, ch="all", n=5, iter=10)
eeg_tmp = ica_remove(eeg, ch="all", ic, ic_mw, ic_idx=1)
@test size(eeg_tmp) == size(eeg)
eeg_tmp = deepcopy(eeg)
add_component!(eeg_tmp, c=:ic, v=ic)
add_component!(eeg_tmp, c=:ic_mw, v=ic_mw)
eeg_tmp = ica_remove(eeg_tmp, ch="all", ic_idx=1);
@test size(eeg_tmp) == size(eeg)

@info "Test 51/$ntests: normpower()"
@test round.(normpower(1:10)) == [6.0, 12.0, 19.0, 25.0, 31.0, 37.0, 43.0, 50.0, 56.0, 62.0]
@test size(normpower(e10, ch="all")) == size(e10)

@info "Test 52/$ntests: sort_epochs()"
e10_erp = erp(e10)
e = sort_epochs(e10_erp, s=collect((nepochs(e10_erp) - 1):-1:1))
@test size(e.data) == (19, 2560, 11)

@info "Test 53/$ntests: denoise_cwt())"
s = denoise_cwt(rand(100), fs=10, nf=2)
@test length(s) == 100
e10_tmp = denoise_cwt(e10, ch="all", nf=50)
@test size(e10_tmp) == size(e10)

@info "Test 54/$ntests: remove_cwt())"
e10_tmp = remove_cwt(e10, ch="Fp1", ep=1, tseg=(0.2, 0.4), fseg=(10, 12.5))
@test size(e10_tmp) == size(e10)

true