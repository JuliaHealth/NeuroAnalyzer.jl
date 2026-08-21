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
@test NeuroAnalyzer.zscore([1, 2, 3]) == [-1.0, 0.0, 1.0]
@test NeuroAnalyzer.zscore(1, 2.0, 3) == -0.3333333333333333

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
@test grubbs([1, 2, 3, 4, 5, 6, 100]) == true
@test outlier_detect(ones(10)) == zeros(10)

@info "Test: cor_test()"
_, _, _, _, df, _ = cor_test(ones(5), zeros(5))
@test df == 8

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

@info "Test: sem_diff()"
@test sem_diff(1:5, 2:6) == 1.0

@info "Test: linreg()"
_, _, c, _, _, _, _ = NeuroAnalyzer.linreg(ones(100), zeros(100))
@test c[1] == 0.0
@test isnan(c[2])

@info "Test: norminv()"
@test norminv(0.5) == 0.0

@info "Test: summary()"
@test NeuroAnalyzer.summary([1, 2.5, 3, NaN, 4, missing, 5.1]) == (
    n = 7,
    ms = 2,
    m = 3.12,
    v = 2.397,
    s = 1.5482,
    mn = 1.0,
    q1 = 2.5,
    me = 3.0,
    q3 = 4.0,
    mx = 5.1,
    mo = 1.0,
)
@test NeuroAnalyzer.summary(rand(10, 3), g = ["g1", "g2", "g3"], d = 2) isa DataFrame
@test NeuroAnalyzer.summary(rand(10), rand(11), rand(12), g = ["g1", "g2", "g3"], d = 2) isa
      DataFrame

@info "Test: cimd()"
@test cimd(collect(1:100)) == (41.0, 60.0)

@info "Test: cir()"
@test cir(r = 0.5, n = 50) == (0.2574878607682835, 0.6832563020988138)
@test cir(r = -0.5, n = 50) == (-0.6832563020988138, -0.2574878607682835)
@test cir([1, 2, 3, 4], [1, 2, 3.1, 4]) == (0.9660280689153002, 0.9999863933445201)

@info "Test: p2z()"
@test p2z(0.05, twotailed = false) == 1.6448536269514717
@test p2z(0.05, twotailed = true) == 1.9599639845400576

@info "Test: distance()"
@test NeuroAnalyzer.distance((0, 0), (1, 1)) == 1.4142135623730951

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

@info "Test: t2p()"
@test t2p(2, df = 2, twotailed = true) == 0.18350341907227394

@info "Test: chi2p()"
@test chi2p(2, df = 2) == 0.36787944117144233

@info "Test: f2p()"
@test f2p(2, df1 = 4, df2 = 3) == 0.2978022709324748

@info "Test: crit_chi()"
@test crit_chi(10) == 3.94029913611906

@info "Test: cis()"
@test NeuroAnalyzer.cis([1, 2, 3, 4]) == (0.7313348599303633, 4.813533834942627)

@info "Test: civ()"
@test NeuroAnalyzer.civ([1, 2, 3, 4]) == (0.5348506773493641, 23.17010798013747)
@test NeuroAnalyzer.cis([1, 2, 3, 4]) == sqrt.(NeuroAnalyzer.civ([1, 2, 3, 4]))

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

true
