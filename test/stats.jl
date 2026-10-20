using DataFrames
using GLMakie
using StatsKit

v = [1, 2, 3, 4, 5]
v1 = [1, 2, 3, 4, 5]
v2 = [6, 5, 4, 3, 2]
m = [1.0 2.0 3.0; 4.0 5.0 6.0]
m1 = [1 2 3; 4 5 6]
m2 = [7 6 5; 4 3 2]
a1 = ones(2, 3, 2)
a2 = zeros(2, 3, 2)

@info "Test: zscore()"
@test zscore([1, 2, 3]) == [-1.0, 0.0, 1.0]
@test zscore(1, 2.0, 3) == -0.3333333333333333

@info "Test: infcrit()"
x = 1:10
y = 1:10
df = DataFrame(:x=>x, :y=>y)
m = GLM.lm(@formula(y ~ x), df)
R2, R2adj, aic, bic = infcrit(m)
@test R2 == 1.0
@test R2adj == 1.0
@test aic == GLM.aicc(m) == -644.2952199530463
@test bic == GLM.bic(m) == -647.3874646740642

@info "Test: outlier_detect()"
@test grubbs([1, 2, 3, 4, 5, 6, 100]) == true
@test outlier_detect(ones(10)) == zeros(10)
@test outlier_detect(ones(10); method = :z) == zeros(10)
@test outlier_detect([zeros(20); 100.0]; method = :z) == [zeros(Bool, 20); true]

@info "Test: cor_test()"
_, _, _, _, df, _ = cor_test(ones(5), zeros(5))
@test df == 3

@info "Test: moe()"
@test moe(100) == 0.1
@test moe(rand(100)) == 0.1

@info "Test: rng()"
@test rng(1:5) == 4

@info "Test: mrng()"
@test mrng(1:5) == 3.0

@info "Test: sem()"
@test sem(1:5) ≈ 0.7071067811865476

@info "Test: semd()"
@test NeuroAnalyzer.semd(1:5) == 0.886004796826744

@info "Test: sem_diff()"
@test sem_diff(1:5, 2:6) == 1.0
# unequal lengths: pooled SD × √(1/n1 + 1/n2)
@test sem_diff(1:5, 1:3) ≈ sqrt((4 * 2.5 + 2 * 1.0) / 6 * (1 / 5 + 1 / 3))

@info "Test: linreg()"
_, _, c, _, _, _, _ = NeuroAnalyzer.linreg(ones(100), zeros(100))
@test c[1] == 0.0
@test isnan(c[2])

@info "Test: norminv()"
@test norminv(0.5) == 0.0

@info "Test: summary_stats() / summary()"
st = summary_stats([1, 2.5, 3, NaN, 4, missing, 5.1])
@test (st.n, st.ms) == (7, 2)
@test st.m ≈ 3.12 && st.v ≈ 2.397 && st.s ≈ sqrt(2.397)
@test (st.mn, st.q1, st.me, st.q3, st.mx, st.mo) == (1.0, 2.5, 3.0, 4.0, 5.1, 1.0)
@test summary_stats(rand(10, 3), g = ["g1", "g2", "g3"], d = 2) isa DataFrame
@test summary(rand(10), rand(11), rand(12), g = ["g1", "g2", "g3"], d = 2) isa DataFrame
# Base.summary(x) and summary(io, x) (used by show()) are unchanged
@test sprint(summary, [1, 2]) == "2-element Vector{Int64}"
@test summary([1, 2]) == "2-element Vector{Int64}"
@test summary(zeros(2, 3)) == "2×3 Matrix{Float64}"

@info "Test: cimd()"
@test cimd(collect(1:100)) == (40.0, 61.0)

@info "Test: cir()"
@test cir(r = 0.5, n = 50) == (0.2574878607682835, 0.6832563020988138)
@test cir(r = -0.5, n = 50) == (-0.6832563020988138, -0.2574878607682835)
@test cir([1, 2, 3, 4], [1, 2, 3.1, 4]) == (0.9660280689153002, 0.9999863933445201)

@info "Test: p2z()"
@test p2z(0.05, twotailed = false) == 1.6448536269514717
@test p2z(0.05, twotailed = true) == 1.9599639845400576

@info "Test: distance()"
@test NeuroAnalyzer.distance((0, 0), (1, 1)) == 1.4142135623730951
@test NeuroAnalyzer.distance((0, 0, 0), (1, 2, 2)) == 3.0

@info "Test: count_thresh()"
m = [1 4 7; 2 5 8; 3 6 9]
@test count_thresh(m, t = 4, t_type = :eq) == (x_t = [0 1 0; 0 0 0; 0 0 0], n = 1)
@test count_thresh(m, t = 4, t_type = :g) == (x_t = [0 0 1; 0 1 1; 0 1 1], n = 5)
@test count_thresh(m, t = 4, t_type = :geq) == (x_t = [0 1 1; 0 1 1; 0 1 1], n = 6)
@test count_thresh(m, t = 4, t_type = :l) == (x_t = [1 0 0; 1 0 0; 1 0 0], n = 3)
@test count_thresh(m, t = 4, t_type = :leq) == (x_t = [1 1 0; 1 0 0; 1 0 0], n = 4)

@info "Test: crit_z()"
@test crit_z(0.05, twotailed = false) == 1.6448536269514717
@test crit_z(0.05, twotailed = true) == 1.9599639845400576

@info "Test: crit_t()"
@test crit_t(20, 0.05, twotailed = false) == 1.7247182429207868
@test crit_t(20, 0.05, twotailed = true) == 2.0859634472658644

@info "Test: z2p()"
@test z2p(1.0) == 0.15865525393145702

@info "Test: bootstrap_ci()"
x = rand(10, 100)
s1, s2, s3 = bootstrap_ci(x)
@test length(s1) == 10
@test length(s2) == 10
@test length(s3) == 10

@info "Test: bootstrap_stat()"
x = rand(10, 100)
s = bootstrap_stat(x; f = "abs(maximum(obj))")
@test length(s) == 3000
s = bootstrap_stat(x; n1 = 100, f = v -> abs(maximum(v)))
@test length(s) == 100 && all(s .> 0)

@info "Test: cip()"
@test cip(0.5, 10) == (0.19010248384771866, 0.8098975161522813)

@info "Test: cl2z()"
@test cl2z(0.95) == 1.9599639845400576

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
@test r1r2_zscore(r1 = 0.3, r2 = 0.6, n1 = 50, n2 = 50) == -1.8597036746544668

@info "Test: varp()"
@test varp(0.5, 10) == 0.025

@info "Test: varc()"
@test varc([0, 1, 2, 3], [2, 8, 27, 45]) == 0.589129780186691

@info "Test: stdp()"
@test stdp(0.5, 10) == 0.15811388300841897

@info "Test: stdc()"
@test stdc([0, 1, 2, 3], [2, 8, 27, 45]) == 0.7675479009069669

@info "Test: sep()"
@test sep(0.5, 10) == 0.15811388300841897

@info "Test: sep_diff()"
@test sep_diff(0.5, 0.6, 10, 15) == 0.2024845673131659

@info "Test: sen()"
@test sen(25) == 5

@info "Test: sen_diff()"
@test sen_diff(10, 15) == 5

@info "Test: rmna()"
@test NeuroAnalyzer.rmna([1, NaN, 2.0, missing]) == [1.0, 2.0]

@info "Test: cim()"
@test cim([1, 2, 3, 4], d = :t) == (0.4457397432394794, 4.554260256760521)
@test cim([1, 2, 3, 4], d = :z) == (1.2348486881183376, 3.765151311881662)
# one-sided: lower / upper confidence bound with the one-tailed critical value
e1 = crit_t(3, 0.05; twotailed = false) * sem([1, 2, 3, 4])
@test cim([1, 2, 3, 4]; twotailed = false) == (2.5 - e1, Inf)
@test cim([1, 2, 3, 4]; twotailed = false, tail = :upper) == (-Inf, 2.5 + e1)
@test_throws ArgumentError cim([1, 2, 3, 4]; twotailed = false, tail = :both)

@info "Test: t2p()"
@test t2p(2, df = 2, twotailed = true) == 0.18350341907227394

@info "Test: chi2p()"
@test chi2p(2, df = 2) == 0.36787944117144233

@info "Test: f2p()"
@test f2p(2, df1 = 4, df2 = 3) == 0.2978022709324748

@info "Test: crit_chi()"
@test crit_chi(10) == 3.94029913611906

@info "Test: cis()"
@test cis([1, 2, 3, 4]) == (0.7313348599303633, 4.813533834942627)
# Base.cis() for numbers is unchanged
@test cis(0.0) == 1.0 + 0.0im

@info "Test: civ()"
@test civ([1, 2, 3, 4]) == (0.5348506773493641, 23.17010798013747)
@test cis([1, 2, 3, 4]) == sqrt.(civ([1, 2, 3, 4]))

@info "Test: p2o()"
@test p2o(0.5) == 1.0

@info "Test: o2p()"
@test o2p(1.0) == 0.5

@info "Test: pcacomp()"
m = rand(4, 5)
p = pcacomp(m)
@test size(p.pc, 1) == 4 && 1 <= size(p.pc, 2) <= 4
df = DataFrame(m, :auto)
p = pcacomp(df, names(df))
@test size(p.pc, 1) == 4 && 1 <= size(p.pc, 2) <= 4

@info "Test: biplot()"
@test biplot(df, names(df)) isa GLMakie.Figure

@info "Test: screeplot()"
@test screeplot(df, names(df)) isa GLMakie.Figure

@info "Test: npca()"
@test npca(m, type = :var, value = 0.9) <= size(m, 2)
@test npca(m, type = :eig, value = 1) <= size(m, 2)

@info "Test: arf()"
df = DataFrame(:sex=>["F", "M", "F", "M", "F"], :group=>[1, 2, 1, 2, 2])
@test arf(df, :sex) == [3.0 2.0 5.0; 0.6 0.4 1.0; 60.0 40.0 100.0]
@test arf(df, :group) == [2.0 3.0 5.0; 0.4 0.6 1.0; 40.0 60.0 100.0]


##
## values: SciPy 1.17 references (scipy.stats) and formulas
##

@info "Test: cor_test() values"
x = [1.0, 2, 3, 4, 5, 6, 7, 8]
y = [2.0, 1, 4, 3, 7, 8, 6, 9]
# scipy.stats.pearsonr: r, p, 95% CI (Fisher z); t = r √(df / (1 - r²))
ct = cor_test(x, y)
@test ct.r ≈ 0.8964214570007951
@test ct.p ≈ 0.002566766096253243 rtol = 1e-6
@test collect(ct.rc) ≈ [0.520608645286509, 0.9812527738012913] rtol = 1e-6
@test ct.ts[1] ≈ 4.95433694306862
@test ct.df == 6
@test_throws ArgumentError cor_test(x, y[1:7])
@test_throws ArgumentError cor_test(x[1:3], y[1:3])

@info "Test: linreg() values"
# scipy.stats.linregress
lr = NeuroAnalyzer.linreg(x, y)
@test lr.c ≈ [0.17857142857142883, 1.0714285714285714]
@test lr.se ≈ [1.0920632036729039, 0.21626073957839226]
@test lr.R2 ≈ 0.8035714285714286
@test lr.R2adj ≈ 0.7708333333333334
@test lr.lf ≈ lr.c[1] .+ lr.c[2] .* x
@test_throws ArgumentError NeuroAnalyzer.linreg(x, y[1:7])

@info "Test: rfz() / r1r2_zscore() values"
@test rfz(0.5) ≈ atanh(0.5)
@test isfinite(rfz(1.0)) && rfz(1.0) > 10
@test r1r2_zscore(r1 = 0.5, r2 = 0.5, n1 = 20, n2 = 30) == 0.0
@test_throws ArgumentError rfz(1.5)
@test_throws ArgumentError r1r2_zscore(r1 = 0.3, r2 = 0.6, n1 = 3, n2 = 50)

@info "Test: ses() / sek() / semd() values"
# SE of skewness and kurtosis, n = 10: 0.687, 1.334 (standard tables)
@test ses(10) ≈ 0.6870429186215167
@test sek(10) ≈ 1.334248769989982
@test ses(collect(1:10)) == ses(10)
@test sek(collect(1:10)) == sek(10)
# SEK = 2 SES √((n² - 1) / ((n - 3)(n + 5)))
@test sek(25) ≈ 2 * ses(25) * sqrt((25^2 - 1) / (22 * 30))
@test_throws ArgumentError ses(2)
@test_throws ArgumentError sek(3)
@test_throws ArgumentError ses([1.0, 2.0])
@test NeuroAnalyzer.semd(x) ≈ 1.253 * std(x) / sqrt(8)
@test sep(0.2, 100) ≈ 0.04
@test sen_diff(9, 16) == 5.0
@test_throws ArgumentError sep(1.2, 10)

@info "Test: cmp_stat()"
d = collect(1.0:10.0)
@test cmp_stat(d, 7.5) == 0.3
@test cmp_stat(d, 7.5; type = :l) == 0.7
@test cmp_stat(d, 0.0) == 1.0
@test_throws ArgumentError cmp_stat(d, 1.0; type = :x)
@test_throws ArgumentError cmp_stat(Float64[], 1.0)

@info "Test: outlier_detect() / grubbs() values"
x = [1.0, 2, 3, 4, 5, 6, 30]
# G = (max - mean) / SD = 2.2355 > critical 2.0200 (two-sided, α = 0.05, n = 7)
@test grubbs(x)
@test grubbs(x; t = 1)
@test !grubbs(x; t = -1)
@test !grubbs([1.0, 2, 3, 4, 5, 6, 7])
@test_throws ArgumentError grubbs(x[1:6])
@test_throws ArgumentError grubbs(x; t = 2)
@test outlier_detect(x; method = :iqr) == [false, false, false, false, false, false, true]
@test outlier_detect(x; method = :g) == [false, false, false, false, false, false, true]
@test outlier_detect(vcat(-30.0, x[1:6]); method = :g) == [true, false, false, false, false, false, false]
@test_throws ArgumentError outlier_detect(x[1:5]; method = :g)
@test_throws ArgumentError outlier_detect(x; method = :x)
# indices of removed elements map back to the original vector
@test NeuroAnalyzer._original_index(2, [1], 5) == 3
@test NeuroAnalyzer._original_index(1, Int64[], 5) == 1

@info "Test: pcacomp() / npca() values"
m = [2.5 2.4 1.0; 0.5 0.7 2.1; 2.2 2.9 0.3; 1.9 2.2 1.7; 3.1 3.0 0.2; 2.3 2.7 1.1; 2.0 1.6 1.9; 1.0 1.1 2.6]
# eigenvalues of the correlation matrix (numpy.linalg.eigvalsh of the covariance of z-scored data)
p = pcacomp(m)
@test p.pc_model.prinvars ≈ [2.7704108324518177, 0.1721504588627501, 0.05743870868543334]
@test p.pcv ≈ [92.34702774839388, 5.7383486287583345, 1.914623622847777]
@test size(p.pc) == (8, 3)
# PC scores are uncorrelated, with variances = eigenvalues
@test abs(cor(p.pc[!, :PC1], p.pc[!, :PC2])) < 1e-10
@test var(p.pc[!, :PC1]) ≈ 2.7704108324518177
p = pcacomp(m; n = 1)
@test size(p.pc) == (8, 1)
@test_throws ArgumentError pcacomp(m; n = 4)
@test_throws ArgumentError pcacomp(m[1:1, :])
@test npca(m; type = :var, value = 0.9) == 1
@test npca(m; type = :var, value = 0.95) == 2
@test npca(m; type = :eig, value = 1) == 1
@test_throws ArgumentError npca(m; type = :var, value = 2)
@test_throws ArgumentError npca(m; type = :x, value = 1)
df = DataFrame(m, ["a", "b", "c"])
@test pcacomp(df, ["a", "b", "c"]).pcv ≈ pcacomp(m).pcv
@test_throws ArgumentError pcacomp(df, ["a", "x"])
mz = copy(m)
NeuroAnalyzer._zstd_columns!(mz)
@test vec(mean(mz; dims = 1)) ≈ zeros(3) atol = 1e-12
@test vec(std(mz; dims = 1)) ≈ ones(3)
mz = [1.0 2.0; 1.0 3.0]
NeuroAnalyzer._zstd_columns!(mz)
@test mz[:, 1] == [0.0, 0.0]
@test NeuroAnalyzer._fit_pca(m, 2, false) isa NeuroAnalyzer.MultivariateStats.PCA

@info "Test: bootstrap_ci() / bootstrap_stat() values"
# the bootstrap CI of the mean of a constant signal is the constant; noise: the CI contains the mean
s = repeat([1.0, 2.0, 3.0], 1, 50)
gm, ll, ul = bootstrap_ci(s; n1 = 100, n2 = 20)
@test gm ≈ [1.0, 2.0, 3.0] && ll ≈ gm && ul ≈ gm
s = randn(NeuroAnalyzer.Random.MersenneTwister(1), 5, 200) .+ 3
gm, ll, ul = bootstrap_ci(s; n1 = 500, n2 = 200)
@test all(ll .< gm .< ul)
@test all(abs.(gm .- 3) .< 0.3)
# width ≈ 2 × 1.96 × SD / √n2
@test all(abs.((ul .- ll) .- 2 * 1.96 / sqrt(200)) .< 0.1)
@test_throws ArgumentError bootstrap_ci(s; cl = 1.5)
@test_throws ArgumentError bootstrap_ci(s; n1 = 0)
@test all(bootstrap_stat(repeat([1.0, -4.0], 1, 10); n1 = 20, n2 = 5, f = v -> maximum(abs, v)) .== 4.0)
@test_throws ArgumentError bootstrap_stat(s; f = "undefined_fn(obj)")
bm = NeuroAnalyzer._bootstrap_mean(repeat([1.0, 2.0], 1, 4), 3)
@test bm == [1.0, 2.0]

@info "Test: permute() values"
# cyclic shifts: every row is a rotation of the input
s = collect(1.0:5.0)
pm = NeuroAnalyzer.permute(s, 20)
@test all(any(circshift(s, k) == pm[i, :] for k = 1:4) for i = 1:20)
pm = NeuroAnalyzer.permute([1.0 2 3; 4 5 6], 5)
@test all(sort(pm[i, 2, :]) == [4.0, 5, 6] for i = 1:5)
@test_throws ArgumentError NeuroAnalyzer.permute([1.0], 5)
@test_throws ArgumentError NeuroAnalyzer.permute(s, 0)

@info "Test: cimd() values"
# order statistics: ranks n/2 ∓ 1.96 √n / 2 (n = 100 → 40, 61)
x = collect(1.0:100.0)
@test NeuroAnalyzer._cimd_sorted(x, 0.95) == (40.0, 61.0)
@test NeuroAnalyzer._cimd_sorted(x, 0.99) == (37.0, 64.0)

@info "Test: summary_stats() helpers"
st = NeuroAnalyzer._summary_stats([1.0, 2.0, 3.0], 4)
@test st.ms == 1 && st.m == 2.0 && st.v == 1.0 && st.me == 2.0
@test isnan(NeuroAnalyzer._summary_stats([5.0], 1).v)
gd = [[1.0, 2.0, 3.0], [4.0, 6.0]]
df = NeuroAnalyzer._build_summary_df(idx -> (gd[idx], length(gd[idx])), ["a", "b"])
@test df[!, :group] == ["a", "b"]
@test df[!, :mean] == [2.0, 5.0]
@test df[!, :n] == [3, 2]
@test_throws ArgumentError NeuroAnalyzer._build_summary_df(idx -> (Float64[], 1), ["a"])
out = mktemp() do path, io
    redirect_stdout(io) do
        NeuroAnalyzer._print_summary_table(["a", "b"], df, 2)
    end
    flush(io)
    read(path, String)
end
@test occursin("median", out) && occursin("5.00", out)

@info "Test: make_table()"
out = mktemp() do path, io
    redirect_stdout(io) do
        make_table(header = ["Group" "A" "B"], data = Any["var1" 1 2.5; "var2" 3.25 4])
    end
    flush(io)
    read(path, String)
end
@test occursin("Group", out) && occursin("2.500", out) && occursin("3.250", out)
# integers are shown without decimals
@test occursin(" 1 ", out)
@test_throws ArgumentError make_table(header = ["A" "B"], data = Any[1 2 3])

@info "Test: misc values"
@test zscore(5, 3, 2) == 1.0
@test moe(400) == 0.05
@test rng([3.0, -1.0, 7.0]) == 8.0
@test mrng([3.0, -1.0, 7.0]) == 3.0
@test p2o(0.75) == 3.0
@test o2p(3.0) == 0.75
@test norminv(0.975) ≈ 1.959963984540054
@test z2p(1.959963984540054; twotailed = true) ≈ 0.05
@test p2z(z2p(1.3)) ≈ 1.3 || p2z(z2p(1.3); twotailed = false) ≈ 1.3
@test distance((1, 2, 3), (1, 2, 3)) == 0.0

true
