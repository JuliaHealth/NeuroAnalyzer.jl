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
    n > p || throw(ArgumentError("Model order ($p) must be smaller than epoch length ($n)."))
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

Select a common VAR order by minimizing BIC over 1:`maxorder`.
"""
function _gc_order(x::AbstractArray{<:Real, 3}, maxorder::Int)
    maxorder >= 1 || throw(ArgumentError("maxorder must be ≥ 1."))
    ch = size(x, 1)
    best_p, best_bic = 1, Inf
    for p in 1:min(maxorder, size(x, 2) - 1)
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
    _granger(s; <keyword arguments>)

Time-domain Granger causality between all channel pairs of a 3-D signal array.

# Arguments

- `s::AbstractArray{<:Real, 3}`: signal array, shape (channels, samples, epochs)
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC
- `maxorder::Int=20`: upper bound for BIC search
- `conditional::Bool=false`: if `true`, condition on all remaining channels (full multivariate GC; needs channels·order ≪ samples·epochs), otherwise pairwise

# Returns

Named tuple:

- `gc::Matrix{Float64}`: channels × channels, `gc[i, j]` = ln(RSS_restricted / RSS_full), i → j
- `f::Matrix{Float64}`, `p::Matrix{Float64}`, `p_adj::Matrix{Float64}`: F statistic, p-value, Benjamini-Hochberg FDR-adjusted p-value (same layout)
- `order::Int64`: model order used
"""
function _granger(
    s::AbstractArray{<:Real, 3};
    order::Union{Int, Nothing} = nothing,
    maxorder::Int = 20,
    conditional::Bool = false,
)
    # validate
    ch = size(s, 1)
    ch >= 2 || throw(ArgumentError("At least 2 channels are required."))
    x = Float64.(s)
    # each epoch is demeaned per channel (no intercept in the model)
    x .-= mean(x, dims = 2)
    p = isnothing(order) ? _gc_order(x, maxorder) : order
    p >= 1 || throw(ArgumentError("order must be ≥ 1."))
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
- `tol::Real=1e-10`: convergence tolerance (maximum absolute change of the factor)

# Returns

- `H::Array{ComplexF64, 3}`: minimum-phase transfer function (n × n × N2), `H(ω) → I` normalized
- `Σ::Matrix{Float64}`: innovation covariance (n × n)
"""
function _gc_wilson(
    S::AbstractArray{<:Complex, 3};
    maxiter::Int = 500,
    tol::Real = 1e-10,
)::Tuple{Array{ComplexF64, 3}, Matrix{Float64}}
    N2 = size(S, 3)
    S = ComplexF64.(S)
    # initial factor: Cholesky factor of the zero-lag auto-covariance
    gam0 = real(ifft(S, 3)[:, :, 1])
    psi = repeat(Matrix{ComplexF64}(cholesky(Symmetric((gam0 + gam0') / 2)).U), 1, 1, N2)
    g = similar(psi)
    converged = false
    for _ in 1:maxiter
        for k in 1:N2
            pinv_k = inv(psi[:, :, k])
            g[:, :, k] = pinv_k * S[:, :, k] * pinv_k' + I
        end
        gp = _gc_plus(g)
        psi_new = similar(psi)
        for k in 1:N2
            psi_new[:, :, k] = psi[:, :, k] * gp[:, :, k]
        end
        d = maximum(abs, psi_new .- psi)
        psi = psi_new
        if d < tol
            converged = true
            break
        end
    end
    converged || _warn("Wilson spectral factorization did not converge in $maxiter iterations.")
    # zero-lag coefficient A0: Σ = A0 A0ᵀ, H = ψ A0⁻¹
    A0 = real(ifft(psi, 3)[:, :, 1])
    A0inv = inv(A0)
    H = similar(psi)
    for k in 1:N2
        H[:, :, k] = psi[:, :, k] * A0inv
    end
    return H, A0 * A0'
end

"""
    _granger_spectral(s; <keyword arguments>)

Spectral Granger causality (Geweke 1982, 1984; formulation of Ding, Chen & Bressler 2006) for a 3-D signal array.

# Arguments

- `s::AbstractArray{<:Real, 3}`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC on the full multichannel model
- `maxorder::Int=20`: upper bound for BIC search
- `nfrq::Int=256`: number of frequencies in 0:fs/2 (the factorization uses 2·(nfrq - 1) points on the full circle)
- `conditional::Bool=false`:
    - `false`: pairwise; a bivariate VAR is fitted for every channel pair, using a common model order
    - `true`: conditional on all remaining channels; one VAR is fitted to all channels and for each source channel the reduced model (all channels except the source) is obtained by Wilson spectral factorization of the corresponding sub-block of the full-model spectrum
- `maxiter::Int=500`, `tol::Real=1e-10`: Wilson factorization settings (`conditional=true` only)

# Returns

Named tuple:

- `gc::Array{Float64, 3}`: channels × channels × nfrq, `gc[i, j, k]` = i → j at `frq[k]`
- `frq::Vector{Float64}`: frequencies in Hz
- `order::Int64`: model order used
"""
function _granger_spectral(
    s::AbstractArray{<:Real, 3};
    fs::Int64,
    order::Union{Int, Nothing} = nothing,
    maxorder::Int = 20,
    nfrq::Int = 256,
    conditional::Bool = false,
    maxiter::Int = 500,
    tol::Real = 1e-10,
)
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    nfrq >= 2 || throw(ArgumentError("nfrq must be ≥ 2."))
    ch = size(s, 1)
    ch >= 2 || throw(ArgumentError("At least 2 channels are required."))
    x = Float64.(s)
    x .-= mean(x, dims = 2)
    p = isnothing(order) ? _gc_order(x, maxorder) : order
    p >= 1 || throw(ArgumentError("order must be ≥ 1."))
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
        # source channel
        for i in 1:ch
            # reduced model: all channels except the source
            r = setdiff(1:ch, i)
            Hr, Σr = _gc_wilson(S[r, r, :]; maxiter = maxiter, tol = tol)
            for k in 1:nfrq
                Hrinv = inv(Hr[:, :, k])
                Hk = H[r, :, k]
                # target channel
                for (t, j) in enumerate(r)
                    # full model normalized so that ε_j ⟂ other innovations:
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

"""
    granger(obj; <keyword arguments>)

Time-domain Granger causality between all pairs of selected channels of a NEURO object.

A vector autoregressive (VAR) model with a common order (BIC) is fitted; lags never cross epoch boundaries, epochs are demeaned per channel and pooled as independent observations.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), at least 2
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC
- `maxorder::Int=20`: upper bound for BIC search
- `conditional::Bool=false`: if `true`, condition on all remaining channels (full multivariate GC), otherwise pairwise

# Returns

Named tuple:

- `gc::Matrix{Float64}`: channels × channels, `gc[i, j]` = ln(RSS_restricted / RSS_full), i → j
- `f::Matrix{Float64}`, `p::Matrix{Float64}`, `p_adj::Matrix{Float64}`: F statistic, p-value, Benjamini-Hochberg FDR-adjusted p-value (same layout)
- `order::Int64`: model order used

# Notes

- Stationarity and Gaussian residuals are assumed (the F-test relies on them); permutation or surrogate testing is more defensible for EEG.
- Volume conduction and a common reference create spurious GC in scalp EEG: use source-space or Laplacian (CSD) data, or interpret with care.
- Short epochs and heavy filtering (high-order or notch filters) distort the VAR fit.
- Conditional GC needs channels·order ≪ samples·epochs; with many channels downsample or select fewer channels.
"""
function granger(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    order::Union{Int, Nothing} = nothing,
    maxorder::Int = 20,
    conditional::Bool = false,
)
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    length(ch) >= 2 || throw(ArgumentError("At least 2 channels must be selected."))

    return _granger(
        @view(obj.data[ch, :, :]);
        order = order,
        maxorder = maxorder,
        conditional = conditional,
    )
end

"""
    granger_spectral(obj; <keyword arguments>)

Spectral Granger causality (Geweke 1982, 1984; formulation of Ding, Chen & Bressler 2006) between all pairs of selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), at least 2
- `order::Union{Int, Nothing}=nothing`: VAR order; if `nothing`, chosen by BIC on the full multichannel model
- `maxorder::Int=20`: upper bound for BIC search
- `nfrq::Int=256`: number of frequencies in 0:fs/2
- `conditional::Bool=false`: `false`: pairwise (bivariate VAR for every channel pair); `true`: conditional on all remaining channels (Wilson spectral factorization of the reduced models)
- `maxiter::Int=500`, `tol::Real=1e-10`: Wilson factorization settings (`conditional=true` only)

# Returns

Named tuple:

- `gc::Array{Float64, 3}`: channels × channels × nfrq, `gc[i, j, k]` = i → j at `frq[k]`
- `frq::Vector{Float64}`: frequencies in Hz
- `order::Int64`: model order used

# Notes

- Averaging `gc[i, j, :]` over 0:fs/2 (trapezoid rule) reproduces the time-domain GC of the same model.
- The reduced model of the conditional GC is of infinite order; its accuracy depends on `nfrq` being large enough that the auto-covariance has decayed within 2·(nfrq - 1) lags (narrow spectral peaks, e.g. high-Q alpha or line noise, need more points).
- No significance test is provided (use permutation/surrogate testing); see also the notes of `granger()`.
"""
function granger_spectral(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    order::Union{Int, Nothing} = nothing,
    maxorder::Int = 20,
    nfrq::Int = 256,
    conditional::Bool = false,
    maxiter::Int = 500,
    tol::Real = 1e-10,
)
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    length(ch) >= 2 || throw(ArgumentError("At least 2 channels must be selected."))

    return _granger_spectral(
        @view(obj.data[ch, :, :]);
        fs = sr(obj),
        order = order,
        maxorder = maxorder,
        nfrq = nfrq,
        conditional = conditional,
        maxiter = maxiter,
        tol = tol,
    )
end
