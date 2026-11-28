# Zapline-plus (method = :zapline of remove_powerline()): removal of power line noise by spatial
# filtering
#
# 1. de Cheveigné A. ZapLine: A simple and effective method to remove power line artifacts.
#    NeuroImage. 2020;207:116356.
# 2. Klug M, Kloosterman NA. Zapline-plus: A Zapline extension for automatic and adaptive removal
#    of frequency-specific noise artifacts in M/EEG. Hum Brain Mapp. 2022;43(9):2743-2758.

# Zapline-plus of channels `ch` (indices) of a continuous object, in-place; returns one row per
# chunk: `:start`, `:stop` (seconds), `:nremove`, `:reduction` (dB), `:prominence` (dB, after
# cleaning)
function _zapline!(
    obj::NeuroAnalyzer.NEURO;
    ch::Vector{Int64},
    fline::Real,
    nremove::Union{Int64, Symbol},
    chunk::Symbol,
    chunk_len::Real,
    min_chunk::Real,
    sigma::Real,
    max_prop::Real,
    nfft::Int64,
    pr::Real,
)::DataFrame
    # validate
    fs = sr(obj)
    _check_var(chunk, [:adaptive, :fixed, :none], "chunk")
    nremove isa Symbol && _check_var(nremove, [:auto], "nremove")
    nremove isa Int64 && nremove < 1 && throw(ArgumentError("nremove must be ≥ 1."))
    chunk_len > 0 || throw(ArgumentError("chunk_len must be > 0."))
    min_chunk > 0 || throw(ArgumentError("min_chunk must be > 0."))
    sigma > 0 || throw(ArgumentError("sigma must be > 0."))
    0 < max_prop <= 1 || throw(ArgumentError("max_prop must be in (0, 1]."))
    nfft >= 16 || throw(ArgumentError("nfft must be ≥ 16."))
    length(ch) >= 2 || throw(ArgumentError("method = :zapline requires at least 2 channels."))
    nremove isa Int64 && nremove >= length(ch) &&
        throw(ArgumentError("nremove must be < number of channels ($(length(ch)))."))

    # samples × channels
    X = Matrix{Float64}(permutedims(obj.data[ch, :, 1]))
    n = size(X, 1)
    period = fs / fline
    n > 2 * period || throw(ArgumentError("Signal is too short for pl_frq = $fline Hz."))

    # chunks (sample ranges)
    rngs = if chunk === :none
        [1:n]
    elseif chunk === :fixed
        _zl_fixed_chunks(n, round(Int64, chunk_len * fs))
    else
        _zl_adaptive_chunks(X, fs, fline, round(Int64, min_chunk * fs))
    end

    nmax = max(1, floor(Int64, max_prop * length(ch)))
    fl_n = fline / fs
    df = DataFrame(
        :start => Float64[], :stop => Float64[], :nremove => Int64[], :reduction => Float64[],
        :prominence => Float64[],
    )

    _info("Removing power line noise at $fline Hz and its harmonics (Zapline-plus, $(length(rngs)) chunk(s))")

    # components needed in heavily contaminated chunks (nremove = :auto)
    k_need = Tuple{Int64, Float64}[]

    for r in rngs
        Xc = X[r, :]
        XX = _zl_smooth(Xc, period)
        R = Xc .- XX
        nf = min(nfft, length(r))
        c1 = _zl_bias(R, fl_n, nf)
        todss, scores = _zl_dss(R' * R, c1)
        clean(k) = (Z = R * todss[:, 1:k]; R .- Z * (Z \ R))
        if nremove === :auto
            # outliers of the DSS scores; then (Zapline-plus) more components while the line
            # peak remains (prominence ≥ pr dB) and decreases, up to nmax (removing components
            # without line noise lowers the spectrum around the peak, so the prominence rises)
            k = _zl_nremove(scores, sigma, nmax)
            Rc = clean(k)
            prom = _zl_prominence(XX .+ Rc, fs, fline)
            while prom >= pr && k < min(nmax, length(scores) - 1)
                Rc_n = clean(k + 1)
                prom_n = _zl_prominence(XX .+ Rc_n, fs, fline)
                prom_n < prom || break
                k, Rc, prom = k + 1, Rc_n, prom_n
            end
            # heavily contaminated chunk (the peak remains, limited by nmax): components that
            # would remove (or reduce most) the peak, along the same path (not applied)
            if prom >= pr && k == nmax
                kn, pn = k, prom
                while pn >= pr && kn < length(scores) - 1
                    p_next = _zl_prominence(XX .+ clean(kn + 1), fs, fline)
                    p_next < pn || break
                    kn, pn = kn + 1, p_next
                end
                kn > k && push!(k_need, (kn, pn))
            end
        else
            k = min(nremove, length(scores))
            Rc = clean(k)
            prom = _zl_prominence(XX .+ Rc, fs, fline)
        end
        X[r, :] = XX .+ Rc
        red = 10 * log10(LinearAlgebra.tr(c1) / LinearAlgebra.tr(_zl_bias(Rc, fl_n, nf)))
        push!(df, ((first(r) - 1) / fs, last(r) / fs, k, red, prom))
    end

    n_left = count(>=(pr), df.prominence)
    if n_left > 0
        msg = "Line noise peak (≥ $pr dB) remains in $n_left of $(nrow(df)) chunk(s)"
        if nremove !== :auto
            _warn("$msg with nremove = $nremove; consider nremove = :auto or a larger nremove.")
        elseif isempty(k_need)
            _warn(
                "$msg: removing more components does not reduce it (line noise in many components or not stationary); consider method = :iir.",
            )
        else
            # smallest max_prop (rounded up to 0.05) for the components needed in all chunks
            kmax = maximum(first, k_need)
            mp = min(ceil(kmax / length(ch) * 20) / 20, 1.0)
            pmax = maximum(last, k_need)
            _warn(
                "$msg: heavily contaminated recording (up to $kmax of $(length(ch)) components needed, max_prop = $max_prop allows $nmax); consider max_prop = $mp" *
                (pmax >= pr ? " (a peak of up to $(round(pmax; digits = 1)) dB would remain)" : "") *
                ", or method = :iir.",
            )
        end
    end

    obj.data[ch, :, 1] = permutedims(X)

    return df
end

# prominence (dB) of the power line peak: mean power within fline ± 0.25 Hz over the median
# power 2–5 Hz away, in the channel-averaged power spectrum of `Y` (samples × channels)
# (the mean, not the maximum, of the peak bins: no bias from the variance of the spectrum)
function _zl_prominence(Y::AbstractMatrix, fs::Int64, fline::Real)::Float64
    wl = min(size(Y, 1), 4 * fs)
    ps = [_psd(@view(Y[:, j]); fs = fs, db = false, wlen = wl, woverlap = wl ÷ 2) for j in axes(Y, 2)]
    p = mean(q.p for q in ps)
    f = ps[1].f
    pk = findall(x -> abs(x - fline) <= 0.25, f)
    nb = findall(x -> 2 <= abs(x - fline) <= 5, f)
    (isempty(pk) || isempty(nb)) && return 0.0
    return 10 * log10(mean(p[pk]) / median(p[nb]))
end

# zero-phase moving average of `T` samples (T may be fractional) along the rows of `X`
# (samples × channels): the average of the piecewise-constant signal over [t - T / 2, t + T / 2]
# around the centre of each sample, from the linearly interpolated cumulative sum (exact zeros of
# the response at the multiples of 1 / T cycles/sample)
function _zl_smooth(X::AbstractMatrix, T::Real)::Matrix{Float64}
    n, c = size(X)
    C = vcat(zeros(1, c), cumsum(X; dims = 1)) # C[k + 1]: integral over [0, k]
    function integral!(out, p)
        k = min(floor(Int64, p), n - 1)
        fr = p - k
        @views out .= C[k + 1, :] .+ fr .* (C[k + 2, :] .- C[k + 1, :])
        return out
    end
    Y = zeros(n, c)
    ia = zeros(c)
    ib = zeros(c)
    h = T / 2
    # at the edges the window (one full period) is shifted inside the signal, so that the line
    # frequency is still cancelled
    Tc = min(Float64(T), Float64(n))
    for t in 1:n
        a = t - 0.5 - h
        a = clamp(a, 0.0, n - Tc)
        b = a + Tc
        integral!(ia, a)
        integral!(ib, b)
        @views Y[t, :] .= (ib .- ia) ./ (b - a)
    end
    return Y
end

# cross-spectral matrix of `R` (samples × channels) at the line frequency `fl` (cycles/sample)
# and its harmonics: Hann-windowed segments of `nfft` samples, 50% overlap
function _zl_bias(R::AbstractMatrix, fl::Real, nfft::Int64)::Matrix{Float64}
    n, c = size(R)
    nfft = min(nfft, n)
    w = DSP.hanning(nfft)
    bins = unique([round(Int64, k * fl * nfft) + 1 for k in 1:floor(Int64, 0.5 / fl) if k * fl < 0.5])
    bins = Base.filter(b -> 2 <= b <= nfft ÷ 2 + 1, bins)
    c1 = zeros(c, c)
    for s in 1:max(1, nfft ÷ 2):(n - nfft + 1)
        Zf = FFTW.fft(w .* @view(R[s:(s + nfft - 1), :]), 1)
        for b in bins
            z = Zf[b, :]
            c1 .+= real.(z * z')
        end
    end
    return c1
end

# DSS: components maximizing the ratio of the bias covariance `c1` to the covariance `c0`;
# returns the unmixing matrix (columns: components, by decreasing score) and the scores
function _zl_dss(c0::AbstractMatrix, c1::AbstractMatrix)::Tuple{Matrix{Float64}, Vector{Float64}}
    e0 = LinearAlgebra.eigen(LinearAlgebra.Symmetric(c0))
    # whitening (directions of negligible variance, e.g. after re-referencing, are dropped)
    keep = e0.values .> 1e-9 * maximum(e0.values)
    W = e0.vectors[:, keep] * LinearAlgebra.Diagonal(1 ./ sqrt.(e0.values[keep]))
    e1 = LinearAlgebra.eigen(LinearAlgebra.Symmetric(W' * c1 * W))
    idx = sortperm(e1.values; rev = true)
    return W * e1.vectors[:, idx], e1.values[idx]
end

# number of components to remove: iterative outlier detection of the (decreasing) DSS scores
# (mean + sigma × SD of the remaining scores), between 1 and nmax
function _zl_nremove(scores::AbstractVector, sigma::Real, nmax::Int64)::Int64
    s = collect(scores)
    k = 0
    while length(s) > 2
        thr = mean(s) + sigma * std(s)
        m = count(>(thr), s)
        m == 0 && break
        k += m
        s = s[(m + 1):end]
    end
    return clamp(k, 1, nmax)
end

# fixed chunks of `len` samples; a last chunk shorter than len / 2 is merged with the previous one
function _zl_fixed_chunks(n::Int64, len::Int64)::Vector{UnitRange{Int64}}
    len >= n && return [1:n]
    starts = collect(1:len:n)
    if length(starts) > 1 && n - starts[end] + 1 < len / 2
        pop!(starts)
    end
    return [starts[i]:(i < length(starts) ? starts[i + 1] - 1 : n) for i in eachindex(starts)]
end

# adaptive chunks (after Zapline-plus): boundaries at the largest changes (above the 95th
# percentile) of the covariance of the signal band-passed at fline ± 3 Hz between consecutive
# 1-second segments, at least `n_min` samples apart and from the edges
function _zl_adaptive_chunks(
    X::AbstractMatrix,
    fs::Int64,
    fline::Real,
    n_min::Int64,
)::Vector{UnitRange{Int64}}
    n, c = size(X)
    n < 2 * n_min && return [1:n]
    # messages of the internal filtering are not shown
    verbose_tmp = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false
    Xf = try
        flt = filter_create(;
            fprototype = :butterworth,
            ftype = :bp,
            cutoff = (max(fline - 3, 0.5), min(fline + 3, fs / 2 - 0.5)),
            fs = fs,
            order = 4,
        )
        reduce(hcat, [_filter_apply(@view(X[:, j]); flt = flt) for j in 1:c])
    finally
        NeuroAnalyzer.verbose = verbose_tmp
    end
    seg = fs
    nseg = n ÷ seg
    nseg < 3 && return [1:n]
    covs = [cov(@view(Xf[((k - 1) * seg + 1):(k * seg), :])) for k in 1:nseg]
    # boundary k: between the segments k and k + 1 (sample k × seg + 1)
    d = [sum(abs2, covs[k + 1] .- covs[k]) for k in 1:(nseg - 1)]
    thr = quantile(d, 0.95)
    bnd = Int64[]
    for k in sortperm(d; rev = true)
        d[k] > thr || break
        p = k * seg + 1
        (p - 1 >= n_min && n - p + 1 >= n_min) || continue
        all(b -> abs(b - p) >= n_min, bnd) || continue
        push!(bnd, p)
    end
    sort!(bnd)
    starts = vcat(1, bnd)
    return [starts[i]:(i < length(starts) ? starts[i + 1] - 1 : n) for i in eachindex(starts)]
end
