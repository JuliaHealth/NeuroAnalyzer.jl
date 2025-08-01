using DataFrames
using Plots
using StatsKit

v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1.0 2.0 3.0; 4.0 5.0 6.0]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "Test: hildebrand_rule()"
@test hildebrand_rule([1, 2, 3]) == 0.0

@info "Test: jaccard()"
@test jaccard(ones(3), zeros(3)) == 0.0

@info "Test: zscore()"
@test NeuroAnalyzer.zscore([1, 2, 3]) == [-1.0, 0.0, 1.0]
@test NeuroAnalyzer.zscore(1, 2.0, 3) == -0.3333333333333333

@info "Test: k_categories()"
@test round(k_categories(10)[1]) == 3.0

@info "Test: efs()"
@test efs([1, 2, 3], [2, 3, 4]) == (d = -1.0, g = -1.224744871391589, Δ = -1.0)

@info "Test: infcrit()"
x = 1:10
y = 1:10
df = DataFrame(:x=>x, :y=>y)
m = GLM.lm(@formula(y ~ x), df)
R2, R2adj, aic, bic = infcrit(m)
@test R2 == 1.0
@test R2adj == 1.0
@test aic == -651.7952199530463
@test bic == -651.9926348600522

@info "Test: outlier_detect()"
@test !grubbs([1, 2, 3, 4, 5])
@test outlier_detect(ones(10)) == zeros(10)

@info "Test: cmp_test()"
_, _, _, df, _ = cmp_test(ones(5), zeros(5), paired=true, type=:p)
@test df == 4
_, p1, p2 = cmp_test(ones(1000), zeros(1000), paired=false, type=:perm)
@test p1 == 0.0
@test p2 == 0.0

@info "Test: cor_test()"
_, _, _, _, df, _ = cor_test(ones(5), zeros(5))
@test df == 8

@info "Test: binom_prob()"
@test round(binom_prob(0.5, 6, 10), digits=2) == 0.21

@info "Test: cvm()"
@test cvm(ones(10)) == 0.0

@info "Test: cvmd()"
@test cvmd(ones(10)) == 0.0

@info "Test: fano()"
@test fano(ones(10)) == 0.0

@info "Test: meang()"
@test NeuroAnalyzer.meang(ones(5)) == 1.0

@info "Test: meanh()"
@test NeuroAnalyzer.meanh(ones(5)) == 1.0

@info "Test: meanw()"
@test NeuroAnalyzer.meanw(ones(5), [1, 2, 3, 4, 5]) == 3.0

@info "Test: efs_p1g()"
@test efs_p1g(0.5) == 1.5707963267948968

@info "Test: efs_p2g()"
@test efs_p2g(0.5, 0.5) == 0.0

@info "Test: moe()"
@test moe(100) == 0.1
@test moe(rand(100)) == 0.1

@info "Test: rng()"
@test rng(1:5) == 4

@info "Test: mrng()"
@test mrng(1:5) == 2.0

@info "Test: sem()"
@test NeuroAnalyzer.sem(1:5) == 0.7071067811865476

@info "Test: semd()"
@test NeuroAnalyzer.semd(1:5) == 0.886004796826744

@info "Test: pred_int()"
@test pred_int(2) == 15.56

@info "Test: sem_diff()"
@test sem_diff(1:5, 2:6) == 1.0

@info "Test: prank()"
@test round.(NeuroAnalyzer.prank([1,2,3]), digits=1) == [0.0, 0.1, 0.2]

@info "Test: linreg()"
_, _, c, _, _, _, _ = NeuroAnalyzer.linreg(ones(100), zeros(100))
@test c[1] == 0.0
@test isnan(c[2])

@info "Test: dprime()"
@test dprime(0.5, 0.5) == (dp=0.0, rb=-0.0)

@info "Test: norminv()"
@test norminv(0.5) == 0.0

@info "Test: dranks()"
@test dranks(1:4) == [1, 2, 3, 3]

@info "Test: res_norm()"
@test res_norm(ones(2))[2] == [0.5]

@info "Test: mcc()"
@test mcc(tp=90, tn=90, fp=10, fn=10) == 0.8

@info "Test: meancirc()"
@test NeuroAnalyzer.meancirc([10, 350]) == 0.0
@test round(NeuroAnalyzer.meancirc([0.17453292519943295, 6.1086523819801535], rad=true)) == 0.0

@info "Test: summary()"
@test NeuroAnalyzer.summary([1, 2.5, 3, NaN, 4, missing, 5.1]) == (n = 7, ms = 2, m = 3.12, v = 2.397, s = 1.5482247898803325, min = 1.0, q1 = 2.5, me = 3.0, q3 = 4.0, max = 5.1, mo = 1.0)
@test NeuroAnalyzer.summary(rand(10, 3), g=["g1", "g2", "g3"], d=2) isa DataFrame
@test NeuroAnalyzer.summary(rand(10), rand(11), rand(12), g=["g1", "g2", "g3"], d=2) isa DataFrame

@info "Test: cimd()"
@test cimd(collect(1:100)) == (41.0, 60.0)

@info "Test: cir()"
@test cir(r=0.5, n=50) == (0.2574878607682835, 0.6832563020988138)
@test cir(r=-0.5, n=50) == (-0.6832563020988138, -0.2574878607682835)
@test cir([1, 2, 3, 4], [1, 2, 3.1, 4]) == (0.9660280689153002, 0.9999863933445201)

@info "Test: p2z()"
@test p2z(0.05, twotailed=false) == 1.6448536269514717
@test p2z(0.05, twotailed=true) == 1.9599639845400576

@info "Test: slope()"
@test slope((0, 0), (1, 1)) == 1.0

@info "Test: distance()"
@test NeuroAnalyzer.distance((0, 0), (1, 1)) == 1.4142135623730951

@info "Test: friedman()"
m = [1 4 7; 2 5 8; 3 6 9]
@test friedman(m) == (q = 6.0, w = 1.0, p = 0.04978706836786394)

@info "Test: count_thresh()"
m = [1 4 7; 2 5 8; 3 6 9]
@test count_thresh(m, t=4, t_type=:eq) == (x_t = [0 1 0; 0 0 0; 0 0 0], n = 1)
@test count_thresh(m, t=4, t_type=:g) == (x_t = [0 0 1; 0 1 1; 0 1 1], n = 5)
@test count_thresh(m, t=4, t_type=:geq) == (x_t = [0 1 1; 0 1 1; 0 1 1], n = 6)
@test count_thresh(m, t=4, t_type=:l) == (x_t = [1 0 0; 1 0 0; 1 0 0], n = 3)
@test count_thresh(m, t=4, t_type=:leq) == (x_t = [1 1 0; 1 0 0; 1 0 0], n = 4)

@info "Test: crit_z()"
@test crit_z(0.05, twotailed=false) == 1.6448536269514717
@test crit_z(0.05, twotailed=true) == 1.9599639845400576

@info "Test: crit_t()"
@test crit_t(20, 0.05, twotailed=false) == 1.7247182429207868
@test crit_t(20, 0.05, twotailed=true) == 2.0859634472658644

@info "Test: size_c2g()"
@test size_c2g(m1=100, s1=10, m2=120) == (n1 = 6, n2 = 6)
@test size_c2g(m1=100, s1=10, m2=120, r=2) == (n1 = 4, n2 = 8)

@info "Test: size_c1g()"
@test size_c1g(m=100, s=10, xbar=120) == 3

@info "Test: size_p1g()"
@test size_p1g(p1=0.40, p2=0.50) == 257

@info "Test: size_p2g()"
@test size_p2g(p1=0.40, p2=0.50) == (n1 = 519, n2 = 519)

@info "Test: power_c1g()"
@test power_c1g(m=0, s=2, xbar=1, n=42) == 0.8854398137187739

@info "Test: power_c2g()"
@test power_c2g(m1=100, s1=10, n1=40, m2=101, s2=10, n2=40) == 0.06517153743820353

@info "Test: power_p1g()"
@test power_p1g(p1=0.20, p2=0.10, n1=15) == 0.3079297312284095

@info "Test: power_p2g()"
@test power_p2g(p1=0.10, p2=0.20, n1=15, n2=25) == 0.11073439642784569

@info "Test: z2p()"
@test z2p(1.0) == 0.15865525393145702

@info "Test: size_c1diff()"
@test size_c1diff(s1=10, s2=20) == 128

@info "Test: size_p1diff()"
@test size_p1diff(p1=0.12, p2=0.09) == 7352

@info "Test: bootstrap_ci()"
x = rand(10, 100)
s1, s2, s3 = bootstrap_ci(x)
@test length(s1) == 10
@test length(s2) == 10
@test length(s3) == 10

@info "Test: bootstrap_stat()"
x = rand(10, 100)
s = bootstrap_stat(x, f="abs(maximum(obj))")
@test length(s) == 3000

@info "Test: f1()"
@test f1(tp=90, tn=90, fp=10, fn=10) == (f1 = 0.9, p = 0.9, r = 0.9)

@info "Test: mscr()"
@test mscr(tp=90, tn=90, fp=10, fn=10) == (mr = 0.1, acc = 0.9)

@info "Test: cosim()"
x = 1:4
y = 101:104
@test cosim(x, y) == 0.9172693928327048

@info "Test: cip()"
@test cip(0.5, 10) == (0.19010248384771866, 0.8098975161522813)

@info "Test: cl2z()"
@test cl2z(0.95) == 1.9599639845400576

@info "Test: stdp()"
@test stdp([1, 2, 3, 4], [5, 6, 7, 8], type=:cohen) == 1.2909944487358056
@test stdp([1, 2, 3, 4], [5, 6, 7, 8], type=:hedges) == 1.118033988749895
@test stdp(std([1, 2, 3, 4]), std([5, 6, 7, 8]), 4, 4, type=:cohen) == 1.2909944487358056
@test stdp(std([1, 2, 3, 4]), std([5, 6, 7, 8]), 4, 4, type=:hedges) == 1.118033988749895

@info "Test: permute()"
s = NeuroAnalyzer.permute(rand(5), 10)
@test size(s) == (10, 5)
s = NeuroAnalyzer.permute(rand(4, 8), 10)
@test size(s) == (10, 4, 8)
s = NeuroAnalyzer.permute(rand(2, 4, 8), 10)
@test size(s) == (10, 2, 4, 8)

@info "Test: r2f()"
@test rfz(0.201) == 0.20377443815685448

@info "Test: r1r2_zscore()"
@test r1r2_zscore(r1=0.3, r2=0.6, n1=50, n2=50) == -1.8597036746544668

@info "Test: mdiff()"
st, sts, p = mdiff(m1, m2, method=:absdiff)
@test length(st) == 6
@test sts == 3.0
@test p in [0.0, 1.0]
st, sts, p = mdiff(a1, a2, method=:absdiff)
@test size(st) == (2, 6)
@test sts == [1.0, 1.0]
@test p == [0.0, 0.0]
st, sts, p = mdiff(m1, m2, method=:diff2int)
@test length(st) == 6
@test sts == 4.666666666666666
@test p == 1.0 || p == 0.0
st, sts, p = mdiff(a1, a2, method=:diff2int)
@test size(st) == (2, 6)
@test sts == [2.0, 2.0]
@test p == [0.0, 0.0]

@info "Test: msci95()"
@test msci95(v1) == (sm = 3.0, ss = 0.7071067811865476, su = 4.385929291125633, sl = 1.6140707088743669)
@test msci95(v2) == (sm = 4.0, ss = 0.7071067811865476, su = 5.385929291125633, sl = 2.614070708874367)
@test msci95(m1) == (sm = [2.5, 3.5, 4.5], ss = [1.4999999999999998, 1.4999999999999998, 1.4999999999999998], su = [5.4399999999999995, 6.4399999999999995, 7.4399999999999995], sl = [-0.4399999999999995, 0.5600000000000005, 1.5600000000000005])
@test msci95(a1) == (sm = [1.0 1.0 1.0; 1.0 1.0 1.0], ss = [0.0 0.0 0.0; 0.0 0.0 0.0], su = [1.0 1.0 1.0; 1.0 1.0 1.0], sl = [1.0 1.0 1.0; 1.0 1.0 1.0])
@test msci95(v1, v2) == (sm = -1.0, ss = 1.0, su = 0.96, sl = -2.96)
@test msci95(m1, m2) == (sm = [-4.0; 2.0;;], ss = [0.8164965809277261; 0.8164965809277261;;], su = [-2.3996667013816566; 3.6003332986183434;;], sl = [-5.600333298618343; 0.39966670138165683;;])
@test msci95(a1, a2) == (sm = [1.0 1.0; 1.0 1.0], ss = [0.0 0.0; 0.0 0.0], su = [1.0 1.0; 1.0 1.0], sl = [1.0 1.0; 1.0 1.0])

@info "Test: ba()"
x = ones(10)
y = ones(10) .+ 0.5
@test ba(x, y) == (m = -0.5, s_u = 0.0, s_d = -0.0)

@info "Test: logit()"
@test NeuroAnalyzer.logit(0.8) == 1.3862943611198908

@info "Test: ss()"
@test ss(1:5) == 10

@info "Test: varp()"
@test var(0.5, 10) == 0.025

@info "Test: varc()"
@test var([0, 1, 2, 3], [2, 8, 27, 45]) == 15.804878048780488

@info "Test: stdp()"
@test std(0.5, 10) == 0.15811388300841897

@info "Test: stdc()"
@test std([0, 1, 2, 3], [2, 8, 27, 45]) == 3.9755349386944756

@info "Test: meanp()"
@test meanp(0.5, 10) == 5

@info "Test: meanc()"
@test meanc([0, 1, 2, 3], [2, 8, 27, 45]) == 2.402439024390244

@info "Test: sep()"
@test sep(0.5, 10) == 0.15811388300841897

@info "Test: sep_diff()"
@test sep_diff(0.5, 0.6, 10, 15) == 0.2024845673131659

@info "Test: sen()"
@test sen(25) == 5

@info "Test: sen_diff()"
@test sen_diff(10, 15) == 5

@info "Test: mde()"
@test mde(n=20, s=1.0) == 0.392404582723816

@info "Test: normalize()"
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

@info "Test: binom_test()"
@test binom_test(5, 10) == (x0 = 5, x1 = 5, p0 = 0.5, p1 = 0.5, ci0 = (0.19010248384771866, 0.8098975161522813), ci1 = (0.19010248384771866, 0.8098975161522813), p = 1.0)
@test binom_test(0.5, 10) == (x0 = 5, x1 = 5, p0 = 0.5, p1 = 0.5, ci0 = (0.19010248384771866, 0.8098975161522813), ci1 = (0.19010248384771866, 0.8098975161522813), p = 1.0)
@test binom_test([true, true, false, false]) == (x0 = 2, x1 = 2, p0 = 0.5, p1 = 0.5, ci0 = (0.0100090038649856, 0.9899909961350144), ci1 = (0.0100090038649856, 0.9899909961350144), p = 1.0)

@info "Test: op()"
@test op([1, 2, 3], [0.5, 2, 1]) == [0.5 2.0 1.0; 1.0 4.0 2.0; 1.5 6.0 3.0]

@info "Test: angle()"
@test NeuroAnalyzer.angle([1.0+1.0im, 2.0+2.0im]) == [0.7853981633974483, 0.7853981633974483]

@info "Test: na()"
@test NeuroAnalyzer.na([1, NaN, 2.0, missing]) == [1.0, 2.0]

@info "Test: dap()"
@test dap(collect(1:20)) == (zs = 0.0, zk = -1.215271012491787, d = 1.4768836338028128, p = 0.4778579258619382)

@info "Test: cim()"
@test cim([1, 2, 3, 4], d=:t) == (0.4457397432394794, 4.554260256760521)
@test cim([1, 2, 3, 4], d=:z) == (1.2348486881183376, 3.765151311881662)

@info "Test: t2p()"
@test t2p(2, df=2, twotailed=true) == 0.18350341907227394

@info "Test: chi2p()"
@test chi2p(2, df=2) == 0.36787944117144233

@info "Test: f2p()"
@test f2p(2, df1=4, df2=3) == 0.2978022709324748

@info "Test: crit_chi()"
@test crit_chi(10) == 3.94029913611906

@info "Test: cis()"
@test cis([1, 2, 3, 4]) == (0.7313348599303633, 4.813533834942627)

@info "Test: civ()"
@test civ([1, 2, 3, 4]) == (0.5348506773493641, 23.17010798013747)
@test cis([1, 2, 3, 4]) == sqrt.(civ([1, 2, 3, 4]))

@info "Test: meant()"
@test NeuroAnalyzer.meant([-100, 2, 3, 4, 5, 6, 7, 8, 9, 100]) == 5.5

@info "Test: size_p()"
@test size_p(p=0.169, E=0.04) == 338
@test size_p(E=0.04) == 601

@info "Test: size_m()"
@test size_m(sigma=15, E=2) == 217

@info "Test: df()"
@test NeuroAnalyzer.df(rand(10)) == 9

@info "Test: center()"
@test center([1, 2, 3, 4]) == [-1.5, -0.5, 0.5, 1.5]

@info "Test: p2o()"
@test p2o(0.5) == 1.0

@info "Test: o2p()"
@test o2p(1.0) == 0.5

@info "Test: pcacomp()"
m = rand(4, 5)
p = pcacomp(m)
@test size(p.pc) == (4, 3)
df = DataFrame(m, :auto)
p = pcacomp(df, names(df))
@test size(p.pc) == (4, 3)

@info "Test: biplot()"
@test biplot(df, names(df)) isa Plots.Plot{Plots.GRBackend}

@info "Test: screeplot()"
@test screeplot(df, names(df)) isa Plots.Plot{Plots.GRBackend}

@info "Test: sdi()"
@test sdi([1, 2, 3, 4, 5], [1, 2, 3, 6, 7]) == 0.6

@info "Test: npca()"
@test npca(m, type=:var, value=0.9) <= size(m, 2)
@test npca(m, type=:eig, value=1) <= size(m, 2)

@info "Test: arf()"
df = DataFrame(:sex=>["F", "M", "F", "M", "F"], :group=>[1, 2, 1, 2, 2])
@test arf(df, :sex) == [3.0 2.0 5.0; 0.6 0.4 1.0; 60.0 40.0 100.0]
@test arf(df, :group) == [2.0 3.0 5.0; 0.4 0.6 1.0; 40.0 60.0 100.0]

true