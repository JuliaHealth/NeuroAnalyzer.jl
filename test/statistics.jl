using NeuroAnalyzer
using Test
using DataFrames
using GLM

@info "Test 1/63: hildebrand_rule()"
@test NeuroAnalyzer.hildebrand_rule([1, 2, 3]) == 0.0

@info "Test 2/63: jaccard_similarity()"
@test NeuroAnalyzer.jaccard_similarity(ones(3), zeros(3)) == 0.0

@info "Test 3/63: z_score()"
@test NeuroAnalyzer.z_score([1, 2, 3]) == [-1.0, 0.0, 1.0]

@info "Test 4/63: k_categories()"
@test round(NeuroAnalyzer.k_categories(10)[1]) == 3.0

@info "Test 5/63: effsize()"
@test NeuroAnalyzer.effsize([1,2,3], [2,3,4]) == (cohen = 1.0, hedges = 1.0)

@info "Test 6/63: infcrit()"
x = 1:10
y = 1:10
df = DataFrame(:x=>x, :y=>y)
m = GLM.lm(@formula(y ~ x), df)
aic, bic = NeuroAnalyzer.infcrit(m)
@test aic == 2.0
@test bic == 2.302585092994046

@info "Test 7/63: outlier_detect()"
@test !NeuroAnalyzer.grubbs([1, 2, 3, 4, 5])
@test NeuroAnalyzer.outlier_detect(ones(10)) == zeros(10)

@info "Test 8/63: seg_mean()"
@test NeuroAnalyzer.seg_mean(ones(5,5,5)) == ones(5)

@info "Test 9/63: seg_mean()"
@test NeuroAnalyzer.seg_mean(ones(5,5,5), ones(5, 5, 5)) == (seg1=ones(5), seg2=ones(5))

@info "Test 10/63: cmp_test()"
_, _, _, df, _ = NeuroAnalyzer.cmp_test(ones(5), zeros(5), paired=true, type=:p)
@test df == 4
_, p1, p2 = NeuroAnalyzer.cmp_test(ones(1000), zeros(1000), paired=false, type=:perm)
@test p1 == 0.0
@test p2 == 0.0

@info "Test 11/63: cor_test()"
_, _, _, _, df, _ = NeuroAnalyzer.cor_test(ones(5), zeros(5))
@test df == 8

@info "Test 12/63: binom_prob()"
@test round(NeuroAnalyzer.binom_prob(0.5, 6, 10), digits=2) == 0.21

@info "Test 13/63: binom_stat()"
@test NeuroAnalyzer.binom_stat(1.0, 10) == (m = 10.0, s = 0.0)

@info "Test 14/63: cvar_mean()"
@test NeuroAnalyzer.cvar_mean(ones(10)) == 0.0

@info "Test 15/63: cvar_median()"
@test NeuroAnalyzer.cvar_median(ones(10)) == 0.0

@info "Test 16/63: cvar()"
@test NeuroAnalyzer.cvar(0.1, 2.0) == 5.0

@info "Test 17/63: meang()"
@test NeuroAnalyzer.meang(ones(5)) == 1.0

@info "Test 18/63: meanh()"
@test NeuroAnalyzer.meanh(ones(5)) == 1.0

@info "Test 19/63: meanw()"
@test NeuroAnalyzer.meanw(ones(5), [1,2,3,4,5]) == 1.0

@info "Test 20/63: effsize()"
@test NeuroAnalyzer.effsize(0.5, 0.5) == 0.0

@info "Test 21/63: moe()"
@test NeuroAnalyzer.moe(100) == 0.1

@info "Test 22/63: rng()"
@test NeuroAnalyzer.rng(1:5) == 4

@info "Test 23/63: se()"
@test NeuroAnalyzer.se(ones(5)) == 0.0

@info "Test 24/63: pred_int()"
@test NeuroAnalyzer.pred_int(2) == 15.56

@info "Test 25/63: sem_diff()"
@test NeuroAnalyzer.sem_diff(ones(5), zeros(5)) == 0.0

@info "Test 26/63: prank()"
@test NeuroAnalyzer.round.(prank([1,2,3]), digits=1) == [0.0, 0.1, 0.2]

@info "Test 27/63: linreg()"
_, _, c, _, _, _, _ = NeuroAnalyzer.linreg(ones(100), zeros(100))
@test c == [0.0, 0.0]

@info "Test 28/63: dprime()"
@test NeuroAnalyzer.dprime(0.5, 0.5) == (dprime=0.0, rb=-0.0)

@info "Test 29/63: norminv()"
@test NeuroAnalyzer.norminv(0.5) == 0.0

@info "Test 30/63: dranks()"
@test NeuroAnalyzer.dranks(1:4) == [1, 2, 3, 3]

@info "Test 31/63: res_norm()"
@test NeuroAnalyzer.res_norm(ones(2))[2] == [0.5]

@info "Test 32/63: mcc()"
@test NeuroAnalyzer.mcc(tp=90, tn=90, fp=10, fn=10) == 0.8

@info "Test 33/63: meanc()"
@test NeuroAnalyzer.meanc([10, 350]) == 0.0
@test round(NeuroAnalyzer.meanc([0.17453292519943295, 6.1086523819801535], rad=true)) == 0.0

@info "Test 34/63: summary()"
@test NeuroAnalyzer.summary(ones(10)) == (mm = 1.0, s = 0.0, me = 1.0, mo = 1.0)
@test NeuroAnalyzer.summary(ones(10), ones(10)) == (mm1 = 1.0, mm2 = 1.0, s1 = 0.0, s2 = 0.0, me1 = 1.0, me2 = 1.0, mo1 = 1.0, mo2 = 1.0)

@info "Test 35/63: ci_median()"
@test ci_median(collect(1:100)) == (42, 59)

@info "Test 36/63: ci_r()"
@test ci_r(r=0.3, n=50) == (0.07, 0.5)
@test ci_r(ones(10), zeros(10)) == (0.01, 0.56)

@info "Test 37/63: crit_z()"
@test crit_z(0.95) == 1.9599639845400576

@info "Test 38/63: r_test()"
@test r_test(r1=0.3, r2=0.6, n1=50, n2=50) == -1.8566613853904539

@info "Test 39/63: slope()"
@test slope((0, 0), (1, 1)) == 1.0

@info "Test 40/63: distance()"
@test distance((0, 0), (1, 1)) == 1.4142135623730951

@info "Test 41/63: friedman()"
m = [1 4 7; 2 5 8; 3 6 9]
@test friedman(m) == (f = 6.0, k = 1.0, p = 0.04978706836786394)

@info "Test 42/63: count_thresh()"
m = [1 4 7; 2 5 8; 3 6 9]
@test count_thresh(m, t=4, t_type=:eq) == (x_t = [0 1 0; 0 0 0; 0 0 0], n = 1)
@test count_thresh(m, t=4, t_type=:g) == (x_t = [0 0 1; 0 1 1; 0 1 1], n = 5)
@test count_thresh(m, t=4, t_type=:geq) == (x_t = [0 1 1; 0 1 1; 0 1 1], n = 6)
@test count_thresh(m, t=4, t_type=:l) == (x_t = [1 0 0; 1 0 0; 1 0 0], n = 3)
@test count_thresh(m, t=4, t_type=:leq) == (x_t = [1 1 0; 1 0 0; 1 0 0], n = 4)

@info "Test 43/63: crit_t()"
@test crit_t(0.95, 20) == 2.085963447265865

@info "Test 44/63: size_c2g()"
@test size_c2g(m1=100, s1=10, m2=120) == (n1 = 4, n2 = 4)
@test size_c2g(m1=100, s1=10, m2=120, r=2) == (n1 = 3, n2 = 6)

@info "Test 45/63: size_c1g()"
@test size_c1g(m0=100, s0=10, m1=120) == 2

@info "Test 46/63: size_p2g()"
@test size_p2g(p1=0.40, p2=0.50) == (n1 = 387, n2 = 387)

@info "Test 47/63: size_p1g()"
@test size_p1g(p0=0.40, p1=0.50) == 191

@info "Test 48/63: power_c2g()"
@test power_c2g(m1=100, s1=10, n1=40, m2=101, s2=10, n2=40) == 0.06517153743820353

@info "Test 49/63: power_c1g()"
@test power_c1g(m0=100, s0=10, m1=101, n1=40) == 0.09217027287864314

@info "Test 50/63: power_p2g()"
@test power_p2g(p1=0.10, p2=0.20, n1=15, n2=25) == 0.11073439642784569

@info "Test 51/63: power_p1g()"
@test power_p1g(p0=0.10, p1=0.20, n1=15) == 0.3079297312284095

@info "Test 52/63: z2pow()"
@test z2pow(0.44) == 0.6700314463394064

@info "Test 53/63: size_c1diff()"
@test size_c1diff(s0=20, s1=10) == 128

@info "Test 54/63: size_p1diff()"
@test size_p1diff(p0=0.12, p1=0.09) == 7352

@info "Test 55/63: bootstrap_ci()"
x = rand(10, 100)
s1, s2, s3 = bootstrap_ci(x)
@test length(s1) == 10
@test length(s2) == 10
@test length(s3) == 10

@info "Test 56/63: bootstrap_stat()"
x = rand(10, 100)
s = bootstrap_stat(x, f="abs(maximum(OBJ))")
@test length(s) == 3000

@info "Test 57/63: seg_extract()"
x = ones(100, 100)
@test seg_extract(x, (10, 10, 20, 20)) == ones(11, 11)
@test seg_extract(x, (10, 10, 20, 20), v=true) == ones(11 * 11)
@test seg_extract(x, (10, 10, 20, 20), c=true) == ones(496)

@info "Test 58/63: f1()"
@test NeuroAnalyzer.f1(tp=90, tn=90, fp=10, fn=10) == (f1 = 0.9, p = 0.9, r = 0.9)

@info "Test 59/63: mscr()"
@test NeuroAnalyzer.mscr(tp=90, tn=90, fp=10, fn=10) == (mr = 0.1, acc = 0.9)

@info "Test 60/63: var_test()"
f, p = NeuroAnalyzer.vartest(e10)
@test size(f) == (23, 23, 10)
@test size(p) == (23, 23, 10)
f, p = NeuroAnalyzer.vartest(e10, e10)
@test size(f) == (23, 23, 10)
@test size(p) == (23, 23, 10)

@info "Test 61/63: fwhm()"
s = generate_gaussian(256, 10, ncyc=2)
@test fwhm(s) == (247, 257, 267)

@info "Test 62/63: fwhm()"
x = 1:4
y = 101:104
@test cosine_similarity(x, y) == 0.9172693928327048

@info "Test 63/63: ci_prop()"
@test ci_prop(0.2, 10) == (-0.0479180129218251, 0.4479180129218251)

true