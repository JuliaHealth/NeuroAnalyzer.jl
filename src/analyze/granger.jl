export granger

# granger(s; order, maxorder, conditional) works on a channels × samples × epochs array, and gc[i,j] means i → j.
# Model: it is a time-domain VAR with BIC order selection. Lags never cross epoch boundaries, and epochs are demeaned and pooled.
# Statistics: it reports ln(RSS ratio), an F-test p-value and Benjamini–Hochberg-adjusted p-values. conditional=true conditions on all other channels.

# Known methodological problems

# Convention: input array is channels × samples × epochs (assumed to match
# NeuroAnalyzer's signal layout - verify). Result[i, j] = influence of channel i ON channel j.

"""
    _gc_design(x, p)

Build the lagged VAR design for `x` (channels × samples × epochs). Lags never cross epoch boundaries. Column `(l-1)*ch + c` holds channel `c` at lag `l`.

# Returns

- `Y::Matrix{Float64}`: targets, (epochs·(n-p)) × ch
- `Z::Matrix{Float64}`: lagged regressors, (epochs·(n-p)) × (ch·p)
"""
function _gc_design(x::AbstractArray{<:Real, 3}, p::Int)
    ch, n, ep = size(x)
    n > p || throw(ArgumentError("model order ($p) must be smaller than epoch length ($n)"))
    m = n - p
    Y = Matrix{Float64}(undef, m * ep, ch)
    Z = Matrix{Float64}(undef, m * ep, ch * p)
    for e in 1:ep
        r = ((e - 1) * m + 1):(e * m)
        Y[r, :] = transpose(@view x[:, (p + 1):n, e])
        for l in 1:p, c in 1:ch
            Z[r, (l - 1) * ch + c] = @view x[c, (p - l + 1):(n - l), e]
        end
    end
    return Y, Z
end

# residual sum of squares of y regressed on Z (QR least squares)
_gc_rss(y::AbstractVector, Z::AbstractMatrix) = sum(abs2, y .- Z * (Z \ y))

"""
    _gc_order(x, maxorder)

Select common VAR order by minimizing BIC over 1:`maxorder`.
"""
function _gc_order(x::AbstractArray{<:Real, 3}, maxorder::Int)
    ch = size(x, 1)
    best_p, best_bic = 1, Inf
    for p in 1:maxorder
        Y, Z = _gc_design(x, p)
        N = size(Y, 1)
        R = Y .- Z * (Z \ Y)
        bic = logdet(Symmetric(R' * R / N)) + log(N) * ch^2 * p / N
        if bic < best_bic
            best_bic, best_p = bic, p
        end
    end
    return best_p
end

# Benjamini-Hochberg FDR over the off-diagonal entries of a p-value matrix
function _gc_fdr(pv::AbstractMatrix{<:Real})
    ch = size(pv, 1)
    idx = [CartesianIndex(i, j) for i in 1:ch for j in 1:ch if i != j]
    p = [pv[k] for k in idx]
    m = length(p)
    o = sortperm(p)
    adj = similar(p)
    running = 1.0
    for k in m:-1:1
        running = min(running, p[o[k]] * m / k)
        adj[o[k]] = running
    end
    out = ones(size(pv))
    for (k, c) in enumerate(idx)
        out[c] = adj[k]
    end
    return out
end

"""
    granger(s; order, maxorder, conditional)

Time-domain Granger causality between all channel pairs.

# Arguments

- `s::AbstractArray{<:Real, 3}`: signal, channels × samples × epochs
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC
- `maxorder::Int=20`: upper bound for BIC search
- `conditional::Bool=false`: if `true`, condition on all remaining channels (full multivariate GC; needs ch·order ≪ samples·epochs), otherwise pairwise

# Returns

Named tuple:
- `gc`: ch × ch matrix, `gc[i, j]` = ln(RSS_restricted / RSS_full), i → j
- `f`, `p`, `p_adj`: F statistic, p-value, BH-FDR adjusted p-value (same layout)
- `order`: model order used

# Notes

- Each epoch is demeaned per channel (no intercept in the model).
- Epochs are pooled as independent observations. Assumes stationarity and Gaussian residuals.
"""
function granger(s::AbstractArray{<:Real, 3}; order::Union{Int, Nothing} = nothing, maxorder::Int = 20, conditional::Bool = false)
    # validate
    ch = size(s, 1)
    ch >= 2 || throw(ArgumentError("at least 2 channels required"))

    x = Float64.(s)
    x .-= mean(x, dims = 2)

    p = isnothing(order) ? _gc_order(x, maxorder) : order
    p >= 1 || throw(ArgumentError("order must be ≥ 1"))

    Y, Z = _gc_design(x, p)
    N = size(Y, 1)
    chcols(c) = c:ch:(ch * p)

    gc = zeros(ch, ch)
    f = zeros(ch, ch)
    pv = ones(ch, ch)

    rss_full = conditional ? [_gc_rss(Y[:, j], Z) for j in 1:ch] : zeros(ch)
    Zpair = Matrix{Float64}(undef, N, 2p)

    for j in 1:ch, i in 1:ch
        i == j && continue
        y = @view Y[:, j]
        if conditional
            keep = setdiff(1:(ch * p), chcols(i))
            rss_r = _gc_rss(y, Z[:, keep])
            rss_u = rss_full[j]
            df2 = N - ch * p
        else
            Zpair[:, 1:p] = Z[:, chcols(j)]
            Zpair[:, (p + 1):(2p)] = Z[:, chcols(i)]
            rss_r = _gc_rss(y, @view Zpair[:, 1:p])
            rss_u = _gc_rss(y, Zpair)
            df2 = N - 2p
        end
        gc[i, j] = log(rss_r / rss_u)
        f[i, j] = ((rss_r - rss_u) / p) / (rss_u / df2)
        pv[i, j] = ccdf(FDist(p, df2), f[i, j])
    end

    return (gc = gc, f = f, p = pv, p_adj = _gc_fdr(pv), order = p)
end

"""
    _gc_var_fit(x, p)

Least-squares VAR(`p`) fit. `x` is channels × samples × epochs.

# Returns

- `A::Array{Float64, 3}`: ch × ch × p, `A[t, c, l]` = effect of channel `c` at lag `l` on target `t`
- `Σ::Matrix{Float64}`: ch × ch residual covariance
"""
function _gc_var_fit(x::AbstractArray{<:Real, 3}, p::Int)
    ch = size(x, 1)
    Y, Z = _gc_design(x, p)
    # (ch·p) × ch, B[(l-1)*ch + c, t]
    B = Z \ Y
    R = Y .- Z * B
    Σ = (R' * R) / size(Y, 1)
    A = Array{Float64, 3}(undef, ch, ch, p)
    for l in 1:p, c in 1:ch, t in 1:ch
        A[t, c, l] = B[(l - 1) * ch + c, t]
    end
    return A, Σ
end

# I added granger_spectral(...; conditional=true) to granger.jl; it uses Wilson spectral factorisation. It is still untested in Julia, because there is no Julia here; I checked it with a line-for-line NumPy port.

# How it works: one VAR is fitted to all channels, which gives the full spectrum S(ω). For each source channel, Wilson’s algorithm factorises the spectrum of all the other channels, which gives the reduced model. GC at each frequency then uses the Geweke 1984 / Ding et al. 2006 formula. That needs one factorisation per source channel, so ch in total.

# I corrected my own formula. My first version only removed the source’s contribution, and for z→y it integrated to 0.16 instead of 0.34. The Geweke/Ding normalisation makes the target’s innovation independent of the others first. That version integrates exactly to the time-domain conditional GC, so it’s the one in the file.

# Checks:

# Test    Result
# Factorising a known VAR spectrum    H and Σ recovered to about 1e-15, in 6–8 iterations
# Common driver: z→x at lag 1, z→y at lag 2   Spurious x→y drops from 0.247 (pairwise) to 0.0006 (conditional). z→x 0.59 and z→y 0.34 match high-order regression (0.59, 0.34)
# Spectrum integrated vs time-domain  Agrees to below 1e-15 on the 3-channel system, a random 5-channel VAR(3), and at 100× amplitude
# 2 channels  Conditional equals pairwise (difference 1.6e-15)

# Things to check:

# @warn on non-convergence: swap it for the package’s own warning helper if it has one.
# Grid size: the factorisation uses 2·(nfrq−1) points on the full circle. Narrow spectral peaks (high-Q alpha, line noise) need a larger nfrq, otherwise the autocovariance wraps around and biases the reduced model. There’s no automatic check for this.
# Conditioning set: conditioning on all other channels needs ch·order to be much smaller than the total number of samples. With dense EEG montages and volume conduction, the full VAR becomes ill-conditioned; use a channel subset or source space.
# No significance test: use permutation or surrogate testing.

"""
    granger_spectral(s; fs, order, maxorder, nfrq, conditional, maxiter, tol)

Spectral Granger causality (Geweke 1982, 1984; formulation of Ding, Chen & Bressler 2006).

# Arguments

- `s::AbstractArray{<:Real, 3}`: signal, channels × samples × epochs
- `fs::Int64`: sampling rate in Hz
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC on the full multichannel model
- `maxorder::Int=20`: upper bound for BIC search
- `nfrq::Int=256`: number of frequencies in 0:fs/2 (the factorisation uses 2·(nfrq - 1) points on the full circle)
- `conditional::Bool=false`
    - `conditional=false`: pairwise; a bivariate VAR is fitted for every channel pair, using a common model order.
    - `conditional=true`: conditional on all remaining channels. One VAR is fitted to all channels; for each source channel the reduced model (all channels except the source) is obtained by Wilson spectral factorisation of the corresponding sub-block of the full-model spectrum.
- `maxiter::Int=500`, `tol::Real=1e-10`: Wilson factorisation settings (`conditional=true` only)

# Returns

Named tuple:
- `gc::Array{Float64, 3}`: ch × ch × nfrq, `gc[i, j, k]` = i → j at `frq[k]`
- `frq::Vector{Float64}`: frequencies in Hz
- `order`: model order used

# Notes
- Averaging `gc[i, j, :]` over 0:fs/2 (trapezoid) reproduces the time-domain GC of the same model: pairwise ln(σ²_reduced / σ²_full) or, for `conditional=true`, ln(Σ_reduced[j, j] / Σ[j, j]).
- No significance test is provided (use permutation/surrogate testing).
- The reduced model is of infinite order; its accuracy depends on `nfrq` being large enough that the autocovariance has
decayed within 2·(nfrq - 1) lags (narrow spectral peaks need more points).
"""
function granger_spectral(s::AbstractArray{<:Real, 3}; fs::Int64, order::Union{Int, Nothing} = nothing, maxorder::Int = 20, nfrq::Int = 256, conditional::Bool = false, maxiter::Int = 500, tol::Real = 1e-10)
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1"))
    nfrq >= 2 || throw(ArgumentError("nfrq must be ≥ 2"))
    ch = size(s, 1)
    ch >= 2 || throw(ArgumentError("at least 2 channels required"))

    x = Float64.(s)
    x .-= mean(x, dims = 2)

    p = isnothing(order) ? _gc_order(x, maxorder) : order
    p >= 1 || throw(ArgumentError("order must be ≥ 1"))

    frq = collect(range(0, fs / 2, length = nfrq))
    # e^{-iωl}, ω = 2πf/fs
    E = [cis(-2π * frq[k] * l / fs) for l in 1:p, k in 1:nfrq]

    gc = zeros(ch, ch, nfrq)

    if conditional
        A, Σ = _gc_var_fit(x, p)
        N2 = 2 * (nfrq - 1)
        H = Array{ComplexF64, 3}(undef, ch, ch, N2)
        S = similar(H)
        for k in 1:N2
            Af = Matrix{ComplexF64}(I, ch, ch)
            for l in 1:p
                Af .-= @view(A[:, :, l]) .* cis(-2π * (k - 1) * l / N2)
            end
            H[:, :, k] = inv(Af)
            S[:, :, k] = H[:, :, k] * Σ * H[:, :, k]'
        end
        # source
        for i in 1:ch
            # reduced model: everything except the source
            r = setdiff(1:ch, i)
            Hr, Σr = _gc_wilson(S[r, r, :]; maxiter = maxiter, tol = tol)
            for k in 1:nfrq
                Hrinv = inv(Hr[:, :, k])
                Hk = H[r, :, k]
                # target
                for (t, j) in enumerate(r)
                    # full model normalised so that ε_j ⟂ other innovations:
                    # column j of H·P⁻¹ equals H·Σ[:, j] / Σ[j, j]
                    q = sum(Hrinv[t, :] .* (Hk * Σ[:, j]))
                    gc[i, j, k] = log(Σr[t, t] * Σ[j, j] / abs2(q))
                end
            end
        end
        return (gc = gc, frq = frq, order = p)
    end

    for a in 1:(ch - 1), b in (a + 1):ch
        A, Σ = _gc_var_fit(@view(x[[a, b], :, :]), p)
        # partial (conditional) variances of the innovations
        # var of ε_a given ε_b
        σa = Σ[1, 1] - Σ[1, 2]^2 / Σ[2, 2]
        # var of ε_b given ε_a
        σb = Σ[2, 2] - Σ[1, 2]^2 / Σ[1, 1]
        for k in 1:nfrq
            Af = Matrix{ComplexF64}(I, 2, 2)
            for l in 1:p
                Af .-= @view(A[:, :, l]) .* E[l, k]
            end
            H = inv(Af)
            S = H * Σ * H'
            # H[target, source]; a → b uses H[2,1], b → a uses H[1,2]
            gc[a, b, k] = log(real(S[2, 2]) / (real(S[2, 2]) - σa * abs2(H[2, 1])))
            gc[b, a, k] = log(real(S[1, 1]) / (real(S[1, 1]) - σb * abs2(H[1, 2])))
        end
    end

    return (gc = gc, frq = frq, order = p)
end

# Volume conduction and common reference: both create spurious GC in scalp EEG. Run it on source-space or Laplacian data, or interpret with care.
# Stationarity: nothing checks it. Short epochs and heavy filtering (high-order or notch filters) distort the VAR fit.
# Conditional mode: it needs ch·order ≪ samples·epochs. With many channels, downsample or restrict channels first.
# BIC order: it is a single common order, and the sample count differs slightly per order, so it is not strictly like-for-like.
# Residuals: the F-test assumes Gaussian residuals. A permutation test is the more defensible option for EEG.
function granger(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, kwargs...)
    # ch_idx = <repo helper that resolves `ch` to indices, honouring the global exclude_bads>
    s = obj.data[ch_idx, :, :]
    return granger(s; kwargs...)

end
