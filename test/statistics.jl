using NeuroAnalyzer
using Test
using DataFrames
using GLM

@test hildebrand_rule([1, 2, 3]) == 0.0
@test jaccard_similarity(ones(3), zeros(3)) == 0.0
@test z_score([1, 2, 3]) == [-1.0, 0.0, 1.0]
@test round(k_categories(10)[1]) == 3.0
@test effsize([1,2,3], [2,3,4]) == (cohen = 1.0, hedges = 1.0)

x = rand(10)
y = rand(10)
df = DataFrame(:x=>x, :y=>y)
m = lm(@formula(y ~ x), df)
@test length(infcrit(m)) == 2

@test grubbs([1, 2, 3, 4, 5]) == false
@test outlier_detect(ones(10)) == zeros(10)

@test seg_mean(ones(5,5,5)) == ones(5)
@test seg2_mean(ones(5,5,5), ones(5, 5, 5)) == (seg1=ones(5), seg2=ones(5))
_, _, _, df, _ = s2_cmp(ones(5), zeros(5), paired=true, type=:p)
@test df == 4
_, p1, p2 = s2_cmp(ones(1000), zeros(1000), paired=false, type=:perm)
@test p1 == 0.0
@test p2 == 0.0
_, _, _, _, df, _ = s2_cor(ones(5), zeros(5))
@test df == 8

@test round(binom_prob(0.5, 6, 10), digits=2) == 0.21
@test binom_stat(1.0, 10) == (10.0, 0.0)
@test cvar_mean(ones(10)) == 0.0
@test cvar_median(ones(10)) == 0.0
@test cvar(0.1, 2.0) == 5.0
@test meang(ones(5)) == 1.0
@test meanh(ones(5)) == 1.0
@test effsize(0.5, 0.5) == 0.0
@test meanw(ones(5), [1,2,3,4,5]) == 1.0
@test moe(100) == 0.1
@test rng(1:5) == 4
@test se(ones(5)) == 0.0
@test pred_int(2) == 15.56
@test sem_diff(ones(5), zeros(5)) == 0.0
@test round.(prank([1,2,3]), digits=1) == [0.0, 0.1, 0.2]

_, _, c, _, _, _, _ = linreg(ones(100), zeros(100))
@test c == [0.0, 0.0]

@test dprime(0.5, 0.5) == (dprime=0.0, rb=-0.0)
@test norminv(0.5) == 0.0
@test dranks(1:4) == [1, 2, 3, 3]
@test res_norm(ones(2))[2] == [0.5]

true