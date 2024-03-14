using NeuroAnalyzer
using Test
using DataFrames
using GLM

@info "test 1/60: hildebrand_rule()"
@test NeuroAnalyzer.hildebrand_rule([1, 2, 3]) == 0.0

@info "test 2/60: jaccard_similarity()"
@test NeuroAnalyzer.jaccard_similarity(ones(3), zeros(3)) == 0.0

@info "test 3/60: z_score()"
@test NeuroAnalyzer.z_score([1, 2, 3]) == [-1.0, 0.0, 1.0]

@info "test 4/60: z_score()"
@test round(NeuroAnalyzer.k_categories(10)[1]) == 3.0

@info "test 5/60: effsize()"
@test NeuroAnalyzer.effsize([1,2,3], [2,3,4]) == (cohen = 1.0, hedges = 1.0)

@info "test 6/60: infcrit()"
x = 1:10
y = 1:10
df = DataFrame(:x=>x, :y=>y)
m = GLM.lm(@formula(y ~ x), df)
aic, bic = NeuroAnalyzer.infcrit(m)
@test aic == 2.0
@test bic == 2.302585092994046

@info "test 7/60: outlier_detect()"
@test NeuroAnalyzer.grubbs([1, 2, 3, 4, 5]) == false
@test NeuroAnalyzer.outlier_detect(ones(10)) == zeros(10)

@info "test 8/60: seg_mean()"
@test NeuroAnalyzer.seg_mean(ones(5,5,5)) == ones(5)

@info "test 9/60: seg_mean()"
@test NeuroAnalyzer.seg_mean(ones(5,5,5), ones(5, 5, 5)) == (seg1=ones(5), seg2=ones(5))

@info "test 10/60: cmp_test()"
_, _, _, df, _ = NeuroAnalyzer.cmp_test(ones(5), zeros(5), paired=true, type=:p)
@test df == 4
_, p1, p2 = NeuroAnalyzer.cmp_test(ones(1000), zeros(1000), paired=false, type=:perm)
@test p1 == 0.0
@test p2 == 0.0

@info "test 11/60: cor_test()"
_, _, _, _, df, _ = NeuroAnalyzer.cor_test(ones(5), zeros(5))
@test df == 8

@info "test 12/60: binom_prob()"
@test round(NeuroAnalyzer.binom_prob(0.5, 6, 10), digits=2) == 0.21

@info "test 13/60: binom_stat()"
@test NeuroAnalyzer.binom_stat(1.0, 10) == (m = 10.0, s = 0.0)

@info "test 14/60: cvar_mean()"
@test NeuroAnalyzer.cvar_mean(ones(10)) == 0.0

@info "test 15/60: cvar_median()"
@test NeuroAnalyzer.cvar_median(ones(10)) == 0.0

@info "test 16/60: cvar()"
@test NeuroAnalyzer.cvar(0.1, 2.0) == 5.0

@info "test 17/60: meang()"
@test NeuroAnalyzer.meang(ones(5)) == 1.0

@info "test 18/60: meanh()"
@test NeuroAnalyzer.meanh(ones(5)) == 1.0

@info "test 19/60: meanw()"
@test NeuroAnalyzer.meanw(ones(5), [1,2,3,4,5]) == 1.0

@info "test 20/60: effsize()"
@test NeuroAnalyzer.effsize(0.5, 0.5) == 0.0

@info "test 21/60: moe()"
@test NeuroAnalyzer.moe(100) == 0.1

@info "test 22/60: rng()"
@test NeuroAnalyzer.rng(1:5) == 4

@info "test 23/60: se()"
@test NeuroAnalyzer.se(ones(5)) == 0.0

@info "test 24/60: pred_int()"
@test NeuroAnalyzer.pred_int(2) == 15.56

@info "test 25/60: sem_diff()"
@test NeuroAnalyzer.sem_diff(ones(5), zeros(5)) == 0.0

@info "test 26/60: prank()"
@test NeuroAnalyzer.round.(prank([1,2,3]), digits=1) == [0.0, 0.1, 0.2]

@info "test 27/60: linreg()"
_, _, c, _, _, _, _ = NeuroAnalyzer.linreg(ones(100), zeros(100))
@test c == [0.0, 0.0]

@info "test 28/60: dprime()"
@test NeuroAnalyzer.dprime(0.5, 0.5) == (dprime=0.0, rb=-0.0)

@info "test 29/60: norminv()"
@test NeuroAnalyzer.norminv(0.5) == 0.0

@info "test 30/60: dranks()"
@test NeuroAnalyzer.dranks(1:4) == [1, 2, 3, 3]

@info "test 31/60: res_norm()"
@test NeuroAnalyzer.res_norm(ones(2))[2] == [0.5]

@info "test 32/60: mcc()"
@test NeuroAnalyzer.mcc(tp=90, tn=90, fp=10, fn=10) == 0.8

@info "test 33/60: meanc()"
@test NeuroAnalyzer.meanc([10, 350]) == 0.0
@test round(NeuroAnalyzer.meanc([0.17453292519943295, 6.1086523819801535], rad=true)) == 0.0

@info "test 34/60: summary()"
@test NeuroAnalyzer.summary(ones(10)) == (mm = 1.0, s = 0.0, me = 1.0, mo = 1.0)
@test NeuroAnalyzer.summary(ones(10), ones(10)) == (mm1 = 1.0, mm2 = 1.0, s1 = 0.0, s2 = 0.0, me1 = 1.0, me2 = 1.0, mo1 = 1.0, mo2 = 1.0)

@info "test 35/60: ci_median()"
@test ci_median(collect(1:100)) == (41, 60)

@info "test 36/60: ci_r()"
@test ci_r(r=0.3, n=50) == (0.02, 0.53)
@test ci_r(ones(10), zeros(10)) == (0.01, 0.64)

@info "test 37/60: ci2z()"
@test ci2z(0.95) == 1.9599639845400576

@info "test 38/60: r_test()"
@test r_test(r1=0.3, r2=0.6, n1=50, n2=50) == -1.8566613853904539

@info "test 39/60: slope()"
@test slope((0, 0), (1, 1)) == 1.0

@info "test 40/60: distance()"
@test distance((0, 0), (1, 1)) == 1.4142135623730951

@info "test 41/60: friedman()"
m = [1 4 7; 2 5 8; 3 6 9]
@test friedman(m) == (f = 6.0, k = 1.0, p = 0.04978706836786394)

@info "test 42/60: count_thresh()"
m = [1 4 7; 2 5 8; 3 6 9]
@test count_thresh(m, t=4, t_type=:eq) == (x_t = [0 1 0; 0 0 0; 0 0 0], n = 1)
@test count_thresh(m, t=4, t_type=:g) == (x_t = [0 0 1; 0 1 1; 0 1 1], n = 5)
@test count_thresh(m, t=4, t_type=:geq) == (x_t = [0 1 1; 0 1 1; 0 1 1], n = 6)
@test count_thresh(m, t=4, t_type=:l) == (x_t = [1 0 0; 1 0 0; 1 0 0], n = 3)
@test count_thresh(m, t=4, t_type=:leq) == (x_t = [1 1 0; 1 0 0; 1 0 0], n = 4)

@info "test 43/60: crit_z()"
@test crit_z(1 - 0.05 / 2) == 1.9599639845400576

@info "test 44/60: size_c2g()"
@test size_c2g(m1=100, s1=10, m2=120) == (n1 = 4, n2 = 4)
@test size_c2g(m1=100, s1=10, m2=120, r=2) == (n1 = 3, n2 = 6)

@info "test 45/60: size_c1g()"
@test size_c1g(m0=100, s0=10, m1=120) == 2

@info "test 46/60: size_p2g()"
@test size_p2g(p1=0.40, p2=0.50) == (n1 = 387, n2 = 387)

@info "test 47/60: size_p1g()"
@test size_p1g(p0=0.40, p1=0.50) == 191

@info "test 48/60: power_c2g()"
@test power_c2g(m1=100, s1=10, n1=40, m2=101, s2=10, n2=40) == 0.06517153743820353

@info "test 49/60: power_c1g()"
@test power_c1g(m0=100, s0=10, m1=101, n1=40) == 0.09217027287864314

@info "test 50/60: power_p2g()"
@test power_p2g(p1=0.10, p2=0.20, n1=15, n2=25) == 0.11073439642784569

@info "test 51/60: power_p1g()"
@test power_p1g(p0=0.10, p1=0.20, n1=15) == 0.3079297312284095

@info "test 52/60: z2pow()"
@test z2pow(0.44) == 0.6700314463394064

@info "test 53/60: size_c1diff()"
@test size_c1diff(s0=20, s1=10) == 128

@info "test 54/60: size_p1diff()"
@test size_p1diff(p0=0.12, p1=0.09) == 7352

@info "test 55/60: bootstrap_ci()"
x = rand(10, 100)
s1, s2, s3 = bootstrap_ci(x)
@test length(s1) == 10
@test length(s2) == 10
@test length(s3) == 10

@info "test 56/60: bootstrap_stat()"
x = rand(10, 100)
s = bootstrap_stat(x, f="abs(maximum(OBJ))")
@test length(s) == 3000

@info "test 57/60: seg_extract()"
x = ones(100, 100)
@test seg_extract(x, (10, 10, 20, 20)) == ones(11, 11)
@test seg_extract(x, (10, 10, 20, 20), v=true) == ones(11 * 11)
@test seg_extract(x, (10, 10, 20, 20), c=true) == ones(496)

@info "test 58/60: f1()"
@test NeuroAnalyzer.f1(tp=90, tn=90, fp=10, fn=10) == (f1 = 0.9, p = 0.9, r = 0.9)

@info "test 59/60: mscr()"
@test NeuroAnalyzer.mscr(tp=90, tn=90, fp=10, fn=10) == (mr = 0.1, acc = 0.9)

@info "test 60/60: var_test()"
f, p = NeuroAnalyzer.vartest(e10)
@test size(f) == (23, 23, 10)
@test size(p) == (23, 23, 10)
f, p = NeuroAnalyzer.vartest(e10, e10)
@test size(f) == (23, 23, 10)
@test size(p) == (23, 23, 10)



true