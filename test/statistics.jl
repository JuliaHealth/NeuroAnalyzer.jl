using NeuroAnalyzer
using Test
using DataFrames
using GLM

@info "test 1/38: hildebrand_rule()"
@test NeuroAnalyzer.hildebrand_rule([1, 2, 3]) == 0.0

@info "test 2/38: jaccard_similarity()"
@test NeuroAnalyzer.jaccard_similarity(ones(3), zeros(3)) == 0.0

@info "test 3/38: z_score()"
@test NeuroAnalyzer.z_score([1, 2, 3]) == [-1.0, 0.0, 1.0]

@info "test 4/38: z_score()"
@test round(NeuroAnalyzer.k_categories(10)[1]) == 3.0

@info "test 5/38: effsize()"
@test NeuroAnalyzer.effsize([1,2,3], [2,3,4]) == (cohen = 1.0, hedges = 1.0)

@info "test 6/38: infcrit()"
x = 1:10
y = 1:10
df = DataFrame(:x=>x, :y=>y)
m = GLM.lm(@formula(y ~ x), df)
aic, bic = NeuroAnalyzer.infcrit(m)
@test aic == 2.0
@test bic == 2.302585092994046

@info "test 7/38: outlier_detect()"
@test NeuroAnalyzer.grubbs([1, 2, 3, 4, 5]) == false
@test NeuroAnalyzer.outlier_detect(ones(10)) == zeros(10)

@info "test 8/38: seg_mean()"
@test NeuroAnalyzer.seg_mean(ones(5,5,5)) == ones(5)

@info "test 9/38: seg_mean()"
@test NeuroAnalyzer.seg_mean(ones(5,5,5), ones(5, 5, 5)) == (seg1=ones(5), seg2=ones(5))

@info "test 10/38: cmp_test()"
_, _, _, df, _ = NeuroAnalyzer.cmp_test(ones(5), zeros(5), paired=true, type=:p)
@test df == 4
_, p1, p2 = NeuroAnalyzer.cmp_test(ones(1000), zeros(1000), paired=false, type=:perm)
@test p1 == 0.0
@test p2 == 0.0

@info "test 11/38: cor_test()"
_, _, _, _, df, _ = NeuroAnalyzer.cor_test(ones(5), zeros(5))
@test df == 8

@info "test 12/38: binom_prob()"
@test round(NeuroAnalyzer.binom_prob(0.5, 6, 10), digits=2) == 0.21

@info "test 13/38: binom_stat()"
@test NeuroAnalyzer.binom_stat(1.0, 10) == (m = 10.0, s = 0.0)

@info "test 14/38: cvar_mean()"
@test NeuroAnalyzer.cvar_mean(ones(10)) == 0.0

@info "test 15/38: cvar_median()"
@test NeuroAnalyzer.cvar_median(ones(10)) == 0.0

@info "test 16/38: cvar()"
@test NeuroAnalyzer.cvar(0.1, 2.0) == 5.0

@info "test 17/38: meang()"
@test NeuroAnalyzer.meang(ones(5)) == 1.0

@info "test 18/38: meanh()"
@test NeuroAnalyzer.meanh(ones(5)) == 1.0

@info "test 19/38: meanw()"
@test NeuroAnalyzer.meanw(ones(5), [1,2,3,4,5]) == 1.0

@info "test 20/38: effsize()"
@test NeuroAnalyzer.effsize(0.5, 0.5) == 0.0

@info "test 21/38: moe()"
@test NeuroAnalyzer.moe(100) == 0.1

@info "test 22/38: rng()"
@test NeuroAnalyzer.rng(1:5) == 4

@info "test 23/38: se()"
@test NeuroAnalyzer.se(ones(5)) == 0.0

@info "test 24/38: pred_int()"
@test NeuroAnalyzer.pred_int(2) == 15.56

@info "test 25/38: sem_diff()"
@test NeuroAnalyzer.sem_diff(ones(5), zeros(5)) == 0.0

@info "test 26/38: prank()"
@test NeuroAnalyzer.round.(prank([1,2,3]), digits=1) == [0.0, 0.1, 0.2]

@info "test 27/38: linreg()"
_, _, c, _, _, _, _ = NeuroAnalyzer.linreg(ones(100), zeros(100))
@test c == [0.0, 0.0]

@info "test 28/38: dprime()"
@test NeuroAnalyzer.dprime(0.5, 0.5) == (dprime=0.0, rb=-0.0)

@info "test 29/38: norminv()"
@test NeuroAnalyzer.norminv(0.5) == 0.0

@info "test 30/38: dranks()"
@test NeuroAnalyzer.dranks(1:4) == [1, 2, 3, 3]

@info "test 31/38: res_norm()"
@test NeuroAnalyzer.res_norm(ones(2))[2] == [0.5]

@info "test 32/38: mcc()"
@test NeuroAnalyzer.mcc(90, 90, 10, 10) == 0.8

@info "test 33/38: meanc()"
@test NeuroAnalyzer.meanc([10, 350]) == 0.0
@test round(NeuroAnalyzer.meanc([0.17453292519943295, 6.1086523819801535], rad=true)) == 0.0

@info "test 34/38: summary()"
@test NeuroAnalyzer.summary(ones(10)) == (mm = 1.0, s = 0.0, me = 1.0, mo = 1.0)
@test NeuroAnalyzer.summary(ones(10), ones(10)) == (mm1 = 1.0, mm2 = 1.0, s1 = 0.0, s2 = 0.0, me1 = 1.0, me2 = 1.0, mo1 = 1.0, mo2 = 1.0)

@info "test 35/38: ci_median()"
@test ci_median(collect(1:100)) == (41, 60)

@info "test 36/38: ci_r()"
@test ci_r(r=0.3, n=50) == (0.02, 0.53)
@test ci_r(ones(10), zeros(10)) == (0.01, 0.64)

@info "test 37/38: ci2z()"
@test ci2z(0.95) == 1.9599639845400576

@info "test 38/38: r_test()"
@test r_test(r1=0.3, r2=0.6, n1=50, n2=50) == -1.8566613853904539

true