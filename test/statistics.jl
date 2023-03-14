using NeuroAnalyzer
using Test
using DataFrames
using GLM

@info "test 1/32: hildebrand_rule()"
@test hildebrand_rule([1, 2, 3]) == 0.0

@info "test 2/32: jaccard_similarity()"
@test jaccard_similarity(ones(3), zeros(3)) == 0.0

@info "test 3/32: z_score()"
@test z_score([1, 2, 3]) == [-1.0, 0.0, 1.0]

@info "test 4/32: z_score()"
@test round(k_categories(10)[1]) == 3.0

@info "test 5/32: effsize()"
@test effsize([1,2,3], [2,3,4]) == (cohen = 1.0, hedges = 1.0)

@info "test 6/32: infcrit()"
x = 1:10
y = 1:10
df = DataFrame(:x=>x, :y=>y)
m = lm(@formula(y ~ x), df)
@test infcrit(m) == (aic = 2.0, bic = 2.302585092994046)

@info "test 7/32: outlier_detect()"
@test grubbs([1, 2, 3, 4, 5]) == false
@test outlier_detect(ones(10)) == zeros(10)

@info "test 8/32: seg_mean()"
@test seg_mean(ones(5,5,5)) == ones(5)

@info "test 9/32: seg_mean()"
@test seg_mean(ones(5,5,5), ones(5, 5, 5)) == (seg1=ones(5), seg2=ones(5))

@info "test 10/32: cmp_test()"
_, _, _, df, _ = cmp_test(ones(5), zeros(5), paired=true, type=:p)
@test df == 4
_, p1, p2 = cmp_test(ones(1000), zeros(1000), paired=false, type=:perm)
@test p1 == 0.0
@test p2 == 0.0

@info "test 11/32: cor_test()"
_, _, _, _, df, _ = cor_test(ones(5), zeros(5))
@test df == 8

@info "test 12/32: binom_prob()"
@test round(binom_prob(0.5, 6, 10), digits=2) == 0.21

@info "test 13/32: binom_stat()"
@test binom_stat(1.0, 10) == (m = 10.0, s = 0.0)

@info "test 14/32: cvar_mean()"
@test cvar_mean(ones(10)) == 0.0

@info "test 15/32: cvar_median()"
@test cvar_median(ones(10)) == 0.0

@info "test 16/32: cvar()"
@test cvar(0.1, 2.0) == 5.0

@info "test 17/32: meang()"
@test meang(ones(5)) == 1.0

@info "test 18/32: meanh()"
@test meanh(ones(5)) == 1.0

@info "test 19/32: meanw()"
@test meanw(ones(5), [1,2,3,4,5]) == 1.0

@info "test 20/32: effsize()"
@test effsize(0.5, 0.5) == 0.0

@info "test 21/32: moe()"
@test moe(100) == 0.1

@info "test 22/32: rng()"
@test rng(1:5) == 4

@info "test 23/32: se()"
@test se(ones(5)) == 0.0

@info "test 24/32: pred_int()"
@test pred_int(2) == 15.56

@info "test 25/32: sem_diff()"
@test sem_diff(ones(5), zeros(5)) == 0.0

@info "test 26/32: prank()"
@test round.(prank([1,2,3]), digits=1) == [0.0, 0.1, 0.2]

@info "test 27/32: linreg()"
_, _, c, _, _, _, _ = linreg(ones(100), zeros(100))
@test c == [0.0, 0.0]

@info "test 28/32: dprime()"
@test dprime(0.5, 0.5) == (dprime=0.0, rb=-0.0)

@info "test 29/32: norminv()"
@test norminv(0.5) == 0.0

@info "test 30/32: dranks()"
@test dranks(1:4) == [1, 2, 3, 3]

@info "test 31/32: res_norm()"
@test res_norm(ones(2))[2] == [0.5]

@info "test 32/32: mcc(()"
@test mcc(90, 90, 10, 10) == 0.8

true