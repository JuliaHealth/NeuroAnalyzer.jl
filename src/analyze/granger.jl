export granger
export granger_spectral

"""
    _gc_design(x, p)

Build the lagged VAR design for `x` (channels × samples × epochs). Lags never cross epoch boundaries. Column `(l-1)*ch + c` holds channel `c` at lag `l`.

# Returns

- `Y::Matrix{Float64}`: targets, (epochs·(n-p)) × ch
- `Z::Matrix{Float64}`: lagged regressors, (epochs·(n-p)) × (ch·p)
"""
function _gc_design(x::AbstractArray{<:Real, 3}, p::Int)
    ch, n, ep = size(x)
    n > p ||
        throw(ArgumentError("Model order ($p) must be smaller than epoch length ($n)."))
    m = n - p
    Y = Matrix{Float64}(undef, m * ep, ch)
    Z = Matrix{Float64}(undef, m * ep, ch * p)
    for e = 1:ep
        r = ((e - 1) * m + 1):(e * m)
        Y[r, :] = transpose(@view x[:, (p + 1):n, e])
        for l = 1:p, c = 1:ch
            Z[r, (l - 1) * ch + c] = @view x[c, (p - l + 1):(n - l), e]
        end
    end
    return Y, Z
end

# residual sum of squares of y regressed on Z (QR least squares)
_gc_rss(y::AbstractVector, Z::AbstractMatrix) = sum(abs2, y .- Z * (Z \ y))

"""
    _gc_gram(x, P; <keyword arguments>)

Accumulate the Gram matrices of the lagged VAR design of order `P` (see `_gc_design()`) block by block, without building the whole design. Every model of order ≤ `P` and every channel subset is then fitted from sub-blocks of these matrices; all of them use the same observations (samples `P+1:n` of every epoch).

# Arguments

- `x::AbstractArray{<:Real, 3}`: signal array, shape (channels, samples, epochs)
- `P::Int`: maximum model order
- `block::Int=4096`: number of rows of the design built at a time

# Returns

Named tuple:

- `G::Matrix{Float64}`: Zᵀ Z, (ch·P) × (ch·P)
- `C::Matrix{Float64}`: Zᵀ Y, (ch·P) × ch
- `YY::Matrix{Float64}`: Yᵀ Y, ch × ch
- `N::Int64`: number of observations
- `ch::Int64`: number of channels
- `P::Int64`: maximum model order
"""
function _gc_gram(x::AbstractArray{<:Real, 3}, P::Int; block::Int = 4096)
    ch, n, ep = size(x)
    P >= 1 || throw(ArgumentError("Model order must be ≥ 1."))
    n > P ||
        throw(ArgumentError("Model order ($P) must be smaller than epoch length ($n)."))
    K = ch * P
    G = zeros(K, K)
    C = zeros(K, ch)
    YY = zeros(ch, ch)
    bl = min(block, n - P)
    Zb = Matrix{Float64}(undef, bl, K)
    Yb = Matrix{Float64}(undef, bl, ch)
    for e = 1:ep, t0 = (P + 1):bl:n
        m = min(bl, n - t0 + 1)
        @inbounds for c = 1:ch, r = 1:m
            Yb[r, c] = x[c, t0 + r - 1, e]
        end
        @inbounds for l = 1:P, c = 1:ch
            col = (l - 1) * ch + c
            for r = 1:m
                Zb[r, col] = x[c, t0 + r - 1 - l, e]
            end
        end
        Z = view(Zb, 1:m, :)
        Y = view(Yb, 1:m, :)
        mul!(G, Z', Z, 1.0, 1.0)
        mul!(C, Z', Y, 1.0, 1.0)
        mul!(YY, Y', Y, 1.0, 1.0)
    end
    return (G = (G + G') ./ 2, C = C, YY = (YY + YY') ./ 2, N = ep * (n - P), ch = ch, P = P)
end

# columns of the Gram design for channels `chs` at lags 1:p (lag-major, `chs` order kept)
_gc_cols(gr::NamedTuple, chs::AbstractVector{<:Integer}, p::Int) =
    [(l - 1) * gr.ch + c for l = 1:p for c in chs]

# residual sum of squares of channel j regressed on the design columns S
function _gc_rss_gram(gr::NamedTuple, j::Int, S::AbstractVector{<:Integer})
    c = gr.C[S, j]
    return gr.YY[j, j] - dot(c, Symmetric(gr.G[S, S]) \ c)
end

"""
    _gc_var_gram(gr, chs, p)

Least-squares VAR(`p`) fit of the channels `chs` from the Gram matrices returned by `_gc_gram()`.

# Returns

- `A::Array{Float64, 3}`: k × k × p, `A[t, c, l]` = effect of channel `chs[c]` at lag `l` on target `chs[t]`
- `Σ::Matrix{Float64}`: k × k residual covariance
"""
function _gc_var_gram(gr::NamedTuple, chs::AbstractVector{<:Integer}, p::Int)
    p <= gr.P || throw(ArgumentError("Model order ($p) must be ≤ $(gr.P)."))
    k = length(chs)
    S = _gc_cols(gr, chs, p)
    Cs = gr.C[S, chs]
    # (k·p) × k, B[(l-1)*k + c, t]
    B = Symmetric(gr.G[S, S]) \ Cs
    Σ = (gr.YY[chs, chs] .- Cs' * B) ./ gr.N
    Σ = (Σ + Σ') ./ 2
    A = Array{Float64, 3}(undef, k, k, p)
    for l = 1:p, c = 1:k, t = 1:k
        A[t, c, l] = B[(l - 1) * k + c, t]
    end
    return A, Σ
end

# BIC-optimal order (1:maxorder) of the VAR model of the channels `chs`
function _gc_bic(gr::NamedTuple, chs::AbstractVector{<:Integer}, maxorder::Int)
    k = length(chs)
    best_p, best_bic = 1, Inf
    for p = 1:min(maxorder, gr.P)
        _, Σ = _gc_var_gram(gr, chs, p)
        F = cholesky(Symmetric(Σ); check = false)
        issuccess(F) || continue
        bic = logdet(F) + log(gr.N) * k^2 * p / gr.N
        if bic < best_bic
            best_bic, best_p = bic, p
        end
    end
    return best_p
end

"""
    _gc_order(x, maxorder)

Select a common VAR order of all channels by minimizing BIC over 1:`maxorder` (all orders are fitted to the same observations).
"""
function _gc_order(x::AbstractArray{<:Real, 3}, maxorder::Int)
    maxorder >= 1 || throw(ArgumentError("maxorder must be ≥ 1."))
    gr = _gc_gram(x, min(maxorder, size(x, 2) - 1))
    return _gc_bic(gr, 1:size(x, 1), maxorder)
end

"""
    _gc_var_fit(x, p)

Least-squares VAR(`p`) fit. `x` is channels × samples × epochs.

# Returns

- `A::Array{Float64, 3}`: ch × ch × p, `A[t, c, l]` = effect of channel `c` at lag `l` on target `t`
- `Σ::Matrix{Float64}`: ch × ch residual covariance
"""
_gc_var_fit(x::AbstractArray{<:Real, 3}, p::Int) =
    _gc_var_gram(_gc_gram(x, p), 1:size(x, 1), p)

# Benjamini-Hochberg FDR over the off-diagonal entries of a p-value matrix
function _gc_fdr(pv::AbstractMatrix{<:Real})
    ch = size(pv, 1)
    idx = [CartesianIndex(i, j) for i = 1:ch for j = 1:ch if i != j]
    adj = p_adjust([pv[k] for k in idx]; method = :fdr)
    out = ones(size(pv))
    for (k, c) in enumerate(idx)
        out[c] = adj[k]
    end
    return out
end

# [g]₊ operator of the Wilson factorization: causal part (positive lags and half of the
# zero lag, upper triangle) of a matrix function on the full unit circle (n × n × N2)
function _gc_plus(g::AbstractArray{ComplexF64, 3})::Array{ComplexF64, 3}
    N2 = size(g, 3)
    gam = ifft(g, 3)
    gam[:, :, 1] = triu(0.5 .* gam[:, :, 1])
    # negative lags
    gam[:, :, (N2 ÷ 2 + 2):N2] .= 0
    return fft(gam, 3)
end

"""
    _gc_wilson(S; <keyword arguments>)

Wilson (1972) spectral matrix factorization `S(ω) = H(ω) Σ H(ω)ᴴ` of a cross-spectral matrix sampled on the full unit circle (n × n × N2, ω_k = 2π(k-1)/N2), as implemented in FieldTrip (`sfactorization_wilson`, Dhamala et al. 2008).

# Arguments

- `S::AbstractArray{<:Complex, 3}`: cross-spectral matrix
- `maxiter::Int=500`: maximum number of iterations
- `tol::Real=1e-10`: convergence tolerance (maximum change of the innovation covariance Σ between iterations, relative to its maximum absolute value)

# Returns

- `H::Array{ComplexF64, 3}`: minimum-phase transfer function (n × n × N2), `H(ω) → I` normalized
- `Σ::Matrix{Float64}`: innovation covariance (n × n)
"""
function _gc_wilson(
    S::AbstractArray{<:Complex, 3};
    maxiter::Int = 500,
    tol::Real = 1e-10,
)::Tuple{Array{ComplexF64, 3}, Matrix{Float64}}
    n, _, N2 = size(S)
    n == 1 && return _gc_cepstral(S)
    S = ComplexF64.(S)
    # initial factor: Cholesky factor of the zero-lag auto-covariance
    gam0 = real(ifft(S, 3)[:, :, 1])
    psi = repeat(Matrix{ComplexF64}(cholesky(Symmetric((gam0 + gam0') / 2)).U), 1, 1, N2)
    psi_new = similar(psi)
    g = similar(psi)
    pf = plan_fft!(g, 3)
    pi = plan_ifft!(g, 3)
    tmp = Matrix{ComplexF64}(undef, n, n)
    Σ_prev = fill(Inf, n, n)
    converged = false
    for _ = 1:maxiter
        for k = 1:N2
            pinv_k = inv(@view psi[:, :, k])
            gk = @view g[:, :, k]
            mul!(tmp, pinv_k, @view S[:, :, k])
            mul!(gk, tmp, pinv_k')
            for d = 1:n
                gk[d, d] += 1
            end
        end
        # [g]₊ in place (see _gc_plus())
        pi * g
        for b = 1:n, a = 1:n
            g[a, b, 1] = a <= b ? 0.5 * g[a, b, 1] : 0.0im
        end
        g[:, :, (N2 ÷ 2 + 2):N2] .= 0
        pf * g
        for k = 1:N2
            mul!(view(psi_new, :, :, k), view(psi, :, :, k), view(g, :, :, k))
        end
        psi, psi_new = psi_new, psi
        # convergence of the innovation covariance Σ = A0 A0ᵀ (A0: zero-lag coefficient);
        # the factor itself stalls at a round-off floor (~1e-7 relative) before Σ converges
        A0 = real(dropdims(sum(psi; dims = 3); dims = 3)) ./ N2
        Σ = A0 * A0'
        d = maximum(abs, Σ .- Σ_prev)
        Σ_prev = Σ
        # relative change: the data may be in any units (e.g. μV)
        if d < tol * maximum(abs, Σ)
            converged = true
            break
        end
    end
    converged ||
        _warn("Wilson spectral factorization did not converge in $maxiter iterations.")
    # zero-lag coefficient A0: Σ = A0 A0ᵀ, H = ψ A0⁻¹
    A0 = real(ifft(psi, 3)[:, :, 1])
    A0inv = inv(A0)
    H = similar(psi)
    for k = 1:N2
        H[:, :, k] = psi[:, :, k] * A0inv
    end
    return H, A0 * A0'
end

"""
    _gc_cepstral(S)

Exact minimum-phase factorization `S(ω) = H(ω) σ² H(ω)ᴴ` of a univariate spectrum (1 × 1 × N2, full unit circle) by the cepstral (Kolmogorov) method: `σ² = exp(mean(log S))`; same output as `_gc_wilson()`.
"""
function _gc_cepstral(S::AbstractArray{<:Complex, 3})::Tuple{Array{ComplexF64, 3}, Matrix{Float64}}
    N2 = size(S, 3)
    c = real(ifft(log.(max.(real(vec(S)), floatmin(Float64)))))
    # causal part of the cepstrum (log of the minimum-phase factor), zero lag excluded:
    # H is normalized to 1 at lag 0
    L = zeros(ComplexF64, N2)
    h = N2 ÷ 2
    L[2:h] = c[2:h]
    # lag N2/2 (even N2) is shared by the causal and anti-causal parts
    L[h + 1] = iseven(N2) ? c[h + 1] / 2 : c[h + 1]
    H = reshape(exp.(fft(L)), 1, 1, N2)
    return H, fill(exp(c[1]), 1, 1)
end

"""
    _gc_spectrum(A, Σ, N2)

Transfer function `H(ω) = (I - Σₗ Aₗ e^{-iωl})⁻¹` and spectral matrix `S(ω) = H Σ Hᴴ` of a VAR model at `N2` frequencies on the full unit circle (ω_k = 2π(k-1)/N2).
"""
function _gc_spectrum(A::AbstractArray{<:Real, 3}, Σ::AbstractMatrix{<:Real}, N2::Int)
    k, _, p = size(A)
    H = Array{ComplexF64, 3}(undef, k, k, N2)
    S = similar(H)
    Af = Matrix{ComplexF64}(undef, k, k)
    for f = 1:N2
        copyto!(Af, I(k))
        for l = 1:p
            Af .-= @view(A[:, :, l]) .* cis(-2π * (f - 1) * l / N2)
        end
        Hf = inv(Af)
        H[:, :, f] = Hf
        S[:, :, f] = Hf * Σ * Hf'
    end
    return H, S
end

"""
    _gc_model(A, Σ, N2; <keyword arguments>)

Granger causality between all channels of one VAR model, conditional on the remaining channels of the model. For every source channel `i` the reduced model (all channels except `i`) is obtained from the full model by Wilson factorization of the corresponding sub-block of its spectral matrix, i.e. without refitting (the reduced model of a VAR process is of infinite order; a refitted finite-order model is biased).

# Arguments

- `A::AbstractArray{<:Real, 3}`, `Σ::AbstractMatrix{<:Real}`: full VAR model (see `_gc_var_gram()`)
- `N2::Int`: number of frequencies on the full unit circle
- `spectral::Bool=true`: if `false`, only the time-domain values are computed
- `maxiter::Int=500`, `tol::Real=1e-10`: Wilson factorization settings

# Returns

- `t::Matrix{Float64}`: k × k, `t[i, j] = ln(Σʳ_jj / Σ_jj)`, i → j
- `s::Array{Float64, 3}`: k × k × (N2÷2 + 1), spectral GC at ω in [0, π] (empty if `spectral=false`); its frequency average equals `t`
"""
function _gc_model(
    A::AbstractArray{<:Real, 3},
    Σ::AbstractMatrix{<:Real},
    N2::Int;
    spectral::Bool = true,
    maxiter::Int = 500,
    tol::Real = 1e-10,
)
    k = size(Σ, 1)
    nf = N2 ÷ 2 + 1
    H, S = _gc_spectrum(A, Σ, N2)
    t = zeros(k, k)
    s = spectral ? zeros(k, k, nf) : zeros(0, 0, 0)
    # source channel
    for i = 1:k
        # reduced model: all channels except the source
        r = setdiff(1:k, i)
        Hr, Σr = _gc_wilson(S[r, r, :]; maxiter = maxiter, tol = tol)
        for (u, j) in enumerate(r)
            t[i, j] = log(Σr[u, u] / Σ[j, j])
        end
        spectral || continue
        for f = 1:nf
            Hrinv = inv(Hr[:, :, f])
            hs = H[r, :, f] * Σ
            # target channel
            for (u, j) in enumerate(r)
                # full model normalized so that ε_j ⟂ other innovations:
                # column j of H·P⁻¹ equals H·Σ[:, j] / Σ[j, j]
                q = sum(Hrinv[u, v] * hs[v, j] for v in eachindex(r))
                s[i, j, f] = log(Σr[u, u] * Σ[j, j] / abs2(q))
            end
        end
    end
    return t, s
end

# model orders: ch × ch matrix (diagonal 0); a common order (conditional) or one order per
# channel pair (pairwise); BIC if `order` is nothing
function _gc_orders(
    gr::NamedTuple,
    order::Union{Int, Nothing},
    maxorder::Int,
    conditional::Bool,
)::Matrix{Int64}
    ch = gr.ch
    o = zeros(Int64, ch, ch)
    if !isnothing(order) || conditional
        p = isnothing(order) ? _gc_bic(gr, 1:ch, maxorder) : order
        for i = 1:ch, j = 1:ch
            i != j && (o[i, j] = p)
        end
    else
        for a = 1:(ch - 1), b = (a + 1):ch
            o[a, b] = o[b, a] = _gc_bic(gr, [a, b], maxorder)
        end
    end
    return o
end

# Granger causality of all channel pairs: time domain (ch × ch) and, if `spectral`,
# frequency domain (ch × ch × (N2÷2 + 1))
function _gc_stats(
    gr::NamedTuple,
    orders::Matrix{Int64},
    conditional::Bool,
    N2::Int;
    spectral::Bool = false,
    maxiter::Int = 500,
    tol::Real = 1e-10,
)
    ch = gr.ch
    gt = zeros(ch, ch)
    gs = spectral ? zeros(ch, ch, N2 ÷ 2 + 1) : zeros(0, 0, 0)
    if conditional
        A, Σ = _gc_var_gram(gr, 1:ch, orders[1, 2])
        t, s = _gc_model(A, Σ, N2; spectral = spectral, maxiter = maxiter, tol = tol)
        gt .= t
        spectral && (gs .= s)
    else
        for a = 1:(ch - 1), b = (a + 1):ch
            A, Σ = _gc_var_gram(gr, [a, b], orders[a, b])
            t, s = _gc_model(A, Σ, N2; spectral = spectral, maxiter = maxiter, tol = tol)
            gt[a, b] = t[1, 2]
            gt[b, a] = t[2, 1]
            if spectral
                gs[a, b, :] = s[1, 2, :]
                gs[b, a, :] = s[2, 1, :]
            end
        end
    end
    return gt, gs
end

# F-test of the nested finite-order regressions (exact restricted model under H₀)
function _gc_ftest(gr::NamedTuple, orders::Matrix{Int64}, conditional::Bool)
    ch = gr.ch
    f = zeros(ch, ch)
    pv = ones(ch, ch)
    if conditional
        p = orders[1, 2]
        df2 = gr.N - ch * p
        Sall = _gc_cols(gr, 1:ch, p)
        rss_u = [_gc_rss_gram(gr, j, Sall) for j = 1:ch]
    end
    for j = 1:ch, i = 1:ch
        i == j && continue
        p = orders[i, j]
        if conditional
            rss_r = _gc_rss_gram(gr, j, _gc_cols(gr, setdiff(1:ch, i), p))
            ru = rss_u[j]
        else
            df2 = gr.N - 2p
            rss_r = _gc_rss_gram(gr, j, _gc_cols(gr, [j], p))
            ru = _gc_rss_gram(gr, j, _gc_cols(gr, [j, i], p))
        end
        f[i, j] = ((rss_r - ru) / p) / (ru / df2)
        pv[i, j] = ccdf(FDist(p, df2), max(f[i, j], 0.0))
    end
    return f, pv
end

"""
    _gc_perm(x, P, orders, conditional, N2, nperm, rng, obs_t, obs_s; <keyword arguments>)

Time-shift permutation test: every channel of every epoch is circularly shifted by an independent random number of samples (at least `max(P + 1, n ÷ 10)` from either end), which destroys the dependencies between channels and keeps the spectrum of each channel. The statistics are recomputed with the observed model orders.

# Returns

- `pt::Matrix{Float64}`: p-values of the time-domain GC
- `ps::Array{Float64, 3}`: pointwise p-values of the spectral GC (empty if `obs_s` is empty)
- `pm::Array{Float64, 3}`: p-values of the spectral GC corrected over frequencies (maximum statistic over frequencies of each channel pair)
"""
function _gc_perm(
    x::Array{Float64, 3},
    P::Int,
    orders::Matrix{Int64},
    conditional::Bool,
    N2::Int,
    nperm::Int,
    rng::Int,
    obs_t::Matrix{Float64},
    obs_s::Array{Float64, 3};
    maxiter::Int = 500,
    tol::Real = 1e-10,
)
    ch, n, ep = size(x)
    spectral = !isempty(obs_s)
    m = max(P + 1, n ÷ 10)
    n - m >= m || throw(
        ArgumentError(
            "Epochs ($n samples) are too short for the time-shift permutation test with model order $P.",
        ),
    )
    # all shifts are drawn first: the result does not depend on the number of threads
    sh = rand(Random.MersenneTwister(rng), m:(n - m), ch, ep, nperm)
    nt = min(Threads.nthreads(), nperm)
    chunks = collect(Iterators.partition(1:nperm, cld(nperm, nt)))
    # BLAS runs single-threaded inside the threaded loop (no oversubscription)
    nb = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    res = try
        tasks = map(chunks) do idx
            Threads.@spawn begin
                y = similar(x)
                ct = zeros(Int, ch, ch)
                cs = zeros(Int, size(obs_s))
                cm = zeros(Int, size(obs_s))
                for k in idx
                    @inbounds for e = 1:ep, c = 1:ch
                        s = sh[c, e, k]
                        for t = 1:n
                            y[c, t, e] = x[c, mod1(t - s, n), e]
                        end
                    end
                    gt, gs = _gc_stats(
                        _gc_gram(y, P),
                        orders,
                        conditional,
                        N2;
                        spectral = spectral,
                        maxiter = maxiter,
                        tol = tol,
                    )
                    ct .+= gt .>= obs_t
                    if spectral
                        cs .+= gs .>= obs_s
                        cm .+= maximum(gs; dims = 3) .>= obs_s
                    end
                end
                (ct, cs, cm)
            end
        end
        fetch.(tasks)
    finally
        BLAS.set_num_threads(nb)
    end
    pt = (1 .+ sum(r[1] for r in res)) ./ (nperm + 1)
    ps = (1 .+ sum(r[2] for r in res)) ./ (nperm + 1)
    pm = (1 .+ sum(r[3] for r in res)) ./ (nperm + 1)
    for i = 1:ch
        pt[i, i] = 1.0
        if spectral
            ps[i, i, :] .= 1.0
            pm[i, i, :] .= 1.0
        end
    end
    return pt, ps, pm
end

# shared setup of _granger() and _granger_spectral(): demeaned data, Gram matrices, orders
function _gc_setup(
    s::AbstractArray{<:Real, 3},
    order::Union{Int, Nothing},
    maxorder::Int,
    conditional::Bool,
    nfrq::Int,
    nperm::Int,
)
    ch = size(s, 1)
    ch >= 2 || throw(ArgumentError("At least 2 channels are required."))
    isnothing(order) || order >= 1 || throw(ArgumentError("order must be ≥ 1."))
    maxorder >= 1 || throw(ArgumentError("maxorder must be ≥ 1."))
    nfrq >= 2 || throw(ArgumentError("nfrq must be ≥ 2."))
    nperm >= 0 || throw(ArgumentError("nperm must be ≥ 0."))
    x = Float64.(s)
    # each epoch is demeaned per channel (no intercept in the model)
    x .-= mean(x; dims = 2)
    P = isnothing(order) ? min(maxorder, size(x, 2) - 1) : order
    gr = _gc_gram(x, P)
    orders = _gc_orders(gr, order, maxorder, conditional)
    # number of regressors of the largest (full) model
    k = (conditional ? ch : 2) * maximum(orders)
    gr.N > k || throw(
        ArgumentError(
            "Too few samples ($(gr.N) observations) for $k model parameters: reduce the model order or the number of channels.",
        ),
    )
    return x, gr, orders
end

"""
    _granger(s; <keyword arguments>)

Time-domain Granger causality between all channel pairs of a 3-D signal array.

# Arguments

- `s::AbstractArray{<:Real, 3}`: signal array, shape (channels, samples, epochs)
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC (for every channel pair if `conditional=false`, for all channels if `conditional=true`)
- `maxorder::Int=20`: upper bound for BIC search
- `conditional::Bool=false`: if `true`, condition on all remaining channels (full multivariate GC; needs channels·order ≪ samples·epochs), otherwise pairwise
- `nfrq::Int=256`: frequency resolution of the spectral factorization used to derive the reduced models (2·(nfrq - 1) points on the full circle)
- `nperm::Int=0`: number of time-shift permutations; 0: no permutation test
- `rng::Int=42`: random seed of the permutation test
- `maxiter::Int=500`, `tol::Real=1e-10`: Wilson factorization settings

# Returns

Named tuple:

- `gc::Matrix{Float64}`: channels × channels, `gc[i, j]` = ln(Σ_reduced / Σ_full) of the innovation variance of j, i → j
- `f::Matrix{Float64}`, `p::Matrix{Float64}`, `p_adj::Matrix{Float64}`: F statistic, p-value, Benjamini-Hochberg FDR-adjusted p-value (same layout)
- `p_perm::Matrix{Float64}`, `p_perm_adj::Matrix{Float64}`: permutation p-value and its FDR-adjusted value (`NaN` if `nperm=0`)
- `order::Matrix{Int64}`: model order used for each channel pair (diagonal: 0)
"""
function _granger(
    s::AbstractArray{<:Real, 3};
    order::Union{Int, Nothing} = nothing,
    maxorder::Int = 20,
    conditional::Bool = false,
    nfrq::Int = 256,
    nperm::Int = 0,
    rng::Int = 42,
    maxiter::Int = 500,
    tol::Real = 1e-10,
)
    x, gr, orders = _gc_setup(s, order, maxorder, conditional, nfrq, nperm)
    ch = gr.ch
    N2 = 2 * (nfrq - 1)
    gc, _ = _gc_stats(gr, orders, conditional, N2; maxiter = maxiter, tol = tol)
    f, pv = _gc_ftest(gr, orders, conditional)
    p_perm = fill(NaN, ch, ch)
    p_perm_adj = fill(NaN, ch, ch)
    if nperm > 0
        p_perm, _, _ = _gc_perm(
            x,
            gr.P,
            orders,
            conditional,
            N2,
            nperm,
            rng,
            gc,
            zeros(0, 0, 0);
            maxiter = maxiter,
            tol = tol,
        )
        p_perm_adj = _gc_fdr(p_perm)
    end
    return (
        gc = gc,
        f = f,
        p = pv,
        p_adj = _gc_fdr(pv),
        p_perm = p_perm,
        p_perm_adj = p_perm_adj,
        order = orders,
    )
end

"""
    _granger_spectral(s; <keyword arguments>)

Spectral Granger causality (Geweke 1982, 1984; formulation of Ding, Chen & Bressler 2006) for a 3-D signal array.

# Arguments

- `s::AbstractArray{<:Real, 3}`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC (for every channel pair if `conditional=false`, for all channels if `conditional=true`)
- `maxorder::Int=20`: upper bound for BIC search
- `nfrq::Int=256`: number of frequencies in 0:fs/2 (the factorization uses 2·(nfrq - 1) points on the full circle)
- `conditional::Bool=false`:
    - `false`: pairwise; a bivariate VAR is fitted for every channel pair
    - `true`: conditional on all remaining channels; one VAR is fitted to all channels and for each source channel the reduced model (all channels except the source) is obtained by Wilson spectral factorization of the corresponding sub-block of the full-model spectrum
- `nperm::Int=0`: number of time-shift permutations; 0: no permutation test
- `rng::Int=42`: random seed of the permutation test
- `maxiter::Int=500`, `tol::Real=1e-10`: Wilson factorization settings

# Returns

Named tuple:

- `gc::Array{Float64, 3}`: channels × channels × nfrq, `gc[i, j, k]` = i → j at `frq[k]`
- `frq::Vector{Float64}`: frequencies in Hz
- `p_perm::Array{Float64, 3}`: pointwise permutation p-value (same layout; `NaN` if `nperm=0`)
- `p_perm_max::Array{Float64, 3}`: permutation p-value corrected over frequencies (maximum statistic; `NaN` if `nperm=0`)
- `order::Matrix{Int64}`: model order used for each channel pair (diagonal: 0)
"""
function _granger_spectral(
    s::AbstractArray{<:Real, 3};
    fs::Int64,
    order::Union{Int, Nothing} = nothing,
    maxorder::Int = 20,
    nfrq::Int = 256,
    conditional::Bool = false,
    nperm::Int = 0,
    rng::Int = 42,
    maxiter::Int = 500,
    tol::Real = 1e-10,
)
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    x, gr, orders = _gc_setup(s, order, maxorder, conditional, nfrq, nperm)
    ch = gr.ch
    N2 = 2 * (nfrq - 1)
    frq = collect(range(0, fs / 2; length = nfrq))
    gt, gc = _gc_stats(
        gr,
        orders,
        conditional,
        N2;
        spectral = true,
        maxiter = maxiter,
        tol = tol,
    )
    p_perm = fill(NaN, ch, ch, nfrq)
    p_perm_max = fill(NaN, ch, ch, nfrq)
    if nperm > 0
        _, p_perm, p_perm_max = _gc_perm(
            x,
            gr.P,
            orders,
            conditional,
            N2,
            nperm,
            rng,
            gt,
            gc;
            maxiter = maxiter,
            tol = tol,
        )
    end
    return (gc = gc, frq = frq, p_perm = p_perm, p_perm_max = p_perm_max, order = orders)
end

"""
    granger(obj; <keyword arguments>)

Time-domain Granger causality between all pairs of selected channels of a NEURO object.

A vector autoregressive (VAR) model is fitted to every channel pair (`conditional=false`) or to all selected channels (`conditional=true`); lags never cross epoch boundaries, epochs are demeaned per channel and pooled as independent observations. The model order is chosen by BIC (for every channel pair or for all channels). The reduced model (without the source channel) is derived from the full model by spectral factorization, not refitted: a refitted finite-order reduced model is misspecified and biases GC (Barnett & Seth 2014). Therefore the time-domain GC equals the frequency average of the spectral GC of `granger_spectral()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), at least 2
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC
- `maxorder::Int=20`: upper bound for BIC search
- `conditional::Bool=false`: if `true`, condition on all remaining channels (full multivariate GC), otherwise pairwise
- `nfrq::Int=256`: frequency resolution of the spectral factorization (increase for narrow spectral peaks)
- `nperm::Int=0`: number of time-shift permutations; 0: no permutation test
- `rng::Int=42`: random seed of the permutation test
- `maxiter::Int=500`, `tol::Real=1e-10`: Wilson factorization settings

# Returns

Named tuple:

- `gc::Matrix{Float64}`: channels × channels, `gc[i, j]` = ln(Σ_reduced / Σ_full) of the innovation variance of j, i → j
- `f::Matrix{Float64}`, `p::Matrix{Float64}`, `p_adj::Matrix{Float64}`: F statistic, p-value, Benjamini-Hochberg FDR-adjusted p-value (same layout)
- `p_perm::Matrix{Float64}`, `p_perm_adj::Matrix{Float64}`: permutation p-value and its FDR-adjusted value (`NaN` if `nperm=0`)
- `order::Matrix{Int64}`: model order used for each channel pair (diagonal: 0)
- `labels::Vector{String}`: channel labels

# Notes

- The F-test compares nested finite-order regressions; it assumes stationarity and Gaussian residuals and with long recordings even negligible GC is significant.
- The permutation test circularly shifts every channel of every epoch by an independent random lag (at least 10% of the epoch length). This keeps the spectrum of each channel and destroys all dependencies between channels; the null hypothesis is "no coupling between any channels". Permutation p-values cannot be smaller than 1 / (nperm + 1).
- Volume conduction and a common reference create spurious GC in scalp EEG: use source-space or Laplacian (CSD) data, or interpret with care.
- Short epochs and heavy filtering (high-order or notch filters) distort the VAR fit.
- Conditional GC needs channels·order ≪ samples·epochs; with many channels downsample or select fewer channels.

# References

1. Barnett, L., & Seth, A. K. (2014). The MVGC multivariate Granger causality toolbox: a new approach to Granger-causal inference. Journal of Neuroscience Methods, 223, 50–68.
"""
function granger(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    order::Union{Int, Nothing} = nothing,
    maxorder::Int = 20,
    conditional::Bool = false,
    nfrq::Int = 256,
    nperm::Int = 0,
    rng::Int = 42,
    maxiter::Int = 500,
    tol::Real = 1e-10,
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)
    # resolve channel names to integer indices, optionally skipping bad channels
    ch = _get_ch(obj, ch, exclude_bad_channels)
    length(ch) >= 2 || throw(ArgumentError("At least 2 channels must be selected."))

    r = _granger(
        @view(obj.data[ch, :, :]);
        order = order,
        maxorder = maxorder,
        conditional = conditional,
        nfrq = nfrq,
        nperm = nperm,
        rng = rng,
        maxiter = maxiter,
        tol = tol,
    )
    return merge(r, (labels = labels(obj)[ch],))
end

"""
    granger_spectral(obj; <keyword arguments>)

Spectral Granger causality (Geweke 1982, 1984; formulation of Ding, Chen & Bressler 2006) between all pairs of selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), at least 2
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC (for every channel pair if `conditional=false`, for all channels if `conditional=true`)
- `maxorder::Int=20`: upper bound for BIC search
- `nfrq::Int=256`: number of frequencies in 0:fs/2
- `conditional::Bool=false`: `false`: pairwise (bivariate VAR for every channel pair); `true`: conditional on all remaining channels (Wilson spectral factorization of the reduced models)
- `nperm::Int=0`: number of time-shift permutations; 0: no permutation test
- `rng::Int=42`: random seed of the permutation test
- `maxiter::Int=500`, `tol::Real=1e-10`: Wilson factorization settings

# Returns

Named tuple:

- `gc::Array{Float64, 3}`: channels × channels × nfrq, `gc[i, j, k]` = i → j at `frq[k]`
- `frq::Vector{Float64}`: frequencies in Hz
- `p_perm::Array{Float64, 3}`: pointwise permutation p-value (same layout; `NaN` if `nperm=0`)
- `p_perm_max::Array{Float64, 3}`: permutation p-value corrected over frequencies: the observed value at each frequency is compared with the maximum over frequencies of every permutation (`NaN` if `nperm=0`)
- `order::Matrix{Int64}`: model order used for each channel pair (diagonal: 0)
- `labels::Vector{String}`: channel labels

# Notes

- Averaging `gc[i, j, :]` over 0:fs/2 (trapezoid rule) reproduces the time-domain GC of `granger()`.
- The reduced models are of infinite order; their accuracy depends on `nfrq` being large enough that the auto-covariance has decayed within 2·(nfrq - 1) lags (narrow spectral peaks, e.g. high-Q alpha or line noise, need more points).
- `p_perm_max` controls the family-wise error over frequencies for each channel pair (not over channel pairs); see also the notes of `granger()`.
"""
function granger_spectral(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    order::Union{Int, Nothing} = nothing,
    maxorder::Int = 20,
    nfrq::Int = 256,
    conditional::Bool = false,
    nperm::Int = 0,
    rng::Int = 42,
    maxiter::Int = 500,
    tol::Real = 1e-10,
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)
    # resolve channel names to integer indices, optionally skipping bad channels
    ch = _get_ch(obj, ch, exclude_bad_channels)
    length(ch) >= 2 || throw(ArgumentError("At least 2 channels must be selected."))

    r = _granger_spectral(
        @view(obj.data[ch, :, :]);
        fs = sr(obj),
        order = order,
        maxorder = maxorder,
        nfrq = nfrq,
        conditional = conditional,
        nperm = nperm,
        rng = rng,
        maxiter = maxiter,
        tol = tol,
    )
    return merge(r, (labels = labels(obj)[ch],))
end
