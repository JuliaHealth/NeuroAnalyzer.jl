# CleanLine (method = :cleanline of remove_powerline()): removal of power line noise by
# multi-taper estimation and subtraction of the line sinusoids in sliding windows, channel by
# channel
#
# 1. Mullen T. CleanLine EEGLAB plugin. San Diego, CA: Neuroimaging Informatics Tools and
#    Resources Clearinghouse (NITRC); 2012.
# 2. Mitra PP, Bokil H. Observed Brain Dynamics. Oxford University Press; 2008 (Chronux:
#    rmlinesmovingwinc).

# CleanLine of channels `ch` (indices) of a continuous object, in-place; returns one row per
# channel and line frequency: `:channel`, `:frequency` (Hz), `:windows` (proportion of windows
# with a significant sinusoid), `:reduction` (reduction of the power at the frequency, dB)
function _cleanline!(
    obj::NeuroAnalyzer.NEURO;
    ch::Vector{Int64},
    fline::Real,
    p::Real,
    wlen::Real,
    wstep::Real,
    bw::Real,
    maxiter::Int64,
)::DataFrame
    # validate
    fs = sr(obj)
    0 < p < 1 || throw(ArgumentError("p must be in (0, 1)."))
    wlen > 0 || throw(ArgumentError("wlen must be > 0."))
    0 < wstep <= wlen || throw(ArgumentError("wstep must be in (0, wlen]."))
    bw > 0 || throw(ArgumentError("bw must be > 0."))
    maxiter >= 1 || throw(ArgumentError("maxiter must be ≥ 1."))
    n_win = round(Int64, wlen * fs)
    n_step = max(1, round(Int64, wstep * fs))
    n = epoch_len(obj)
    n_win <= n || throw(ArgumentError("wlen must be ≤ the signal length ($(n / fs) s)."))
    # tapers: time-bandwidth product NW = wlen × bw / 2, K = 2NW - 1 tapers
    nw = wlen * bw / 2
    k = floor(Int64, 2 * nw - 1)
    k >= 2 || throw(ArgumentError("wlen × bw must be ≥ 3 (at least 2 tapers)."))

    # the fundamental and its harmonics below the Nyquist frequency (scan band within it)
    frqs = [h * fline for h in 1:floor(Int64, fs / 2 / fline) if h * fline + bw / 2 < fs / 2]
    isempty(frqs) && throw(ArgumentError("pl_frq + bw / 2 must be < $(fs / 2) Hz (Nyquist)."))

    clabels = labels(obj)
    df = DataFrame(
        :channel => String[], :frequency => Float64[], :windows => Float64[],
        :reduction => Float64[],
    )

    _info("Removing power line noise at $fline Hz and its harmonics (CleanLine, $(length(ch)) channel(s))")

    res = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, length(ch))
    progbar = Progress(length(ch); dt = 1, barlen = 20, color = :white, enabled = progress_bar)
    Threads.@threads :static for i in eachindex(ch)
        res[i] = _cleanline_channel(
            obj.data[ch[i], :, 1], fs, frqs, p, n_win, n_step, nw, k, maxiter,
        )
        progress_bar && next!(progbar)
    end
    progress_bar && finish!(progbar)

    for (i, ch_idx) in enumerate(ch)
        x0 = obj.data[ch_idx, :, 1]
        y, sig = res[i]
        p0, f = _psd(x0; fs = fs, db = true)
        p1, _ = _psd(y; fs = fs, db = true)
        for (j, f0) in enumerate(frqs)
            b = vsearch(f0, f)
            push!(df, (clabels[ch_idx], f0, sig[j], p0[b] - p1[b]))
        end
        obj.data[ch_idx, :, 1] = y
    end

    return df
end

# CleanLine of one signal: in windows of `n_win` samples (step `n_step`), the sinusoids at the
# line frequencies `frqs` (each scanned within ± bw / 2 for the largest Thomson F statistic)
# that are significant (p < `p`, corrected for the scanned frequencies) are estimated with `k` Slepian tapers (time-bandwidth `nw`) and
# subtracted, repeated up to `maxiter` times; the fitted sinusoids of overlapping windows are
# combined with smooth (Hann) weights and subtracted from the signal; returns the cleaned
# signal and, for each line frequency, the proportion of windows with a significant sinusoid
function _cleanline_channel(
    x::Vector{Float64},
    fs::Int64,
    frqs::Vector{<:Real},
    p::Real,
    n_win::Int64,
    n_step::Int64,
    nw::Real,
    k::Int64,
    maxiter::Int64,
)::Tuple{Vector{Float64}, Vector{Float64}}
    n = length(x)
    V = DSP.dpss(n_win, nw, k)
    V0 = vec(sum(V; dims = 1))
    sV0 = sum(abs2, V0)
    # zero-padded FFT: frequency resolution fs / nfft
    nfft = 2 * nextpow(2, n_win)
    f = (0:(nfft ÷ 2)) .* (fs / nfft)
    half_bw = nw / (n_win / fs) # bw / 2
    bins = [findall(x -> abs(x - f0) <= half_bw, f) for f0 in frqs]
    # the largest F of the scanned frequencies is tested: Bonferroni correction for the number
    # of independent frequencies in the scanned band (bw × wlen, resolution 1 / wlen)
    n_ind = max(1, round(Int64, 2 * half_bw * n_win / fs))
    f_crit = quantile(FDist(2, 2 * k - 2), 1 - p / n_ind)
    t = (0:(n_win - 1)) ./ fs
    # overlap-add weights (Hann without zeros)
    w = [0.5 - 0.5 * cos(2π * i / (n_win + 1)) for i in 1:n_win]

    starts = collect(1:n_step:(n - n_win + 1))
    last(starts) == n - n_win + 1 || push!(starts, n - n_win + 1)

    fit = zeros(n)
    wsum = zeros(n)
    nsig = zeros(Int64, length(frqs))
    buf = zeros(nfft, k)
    plan = FFTW.plan_rfft(buf, 1)
    for s in starts
        r = x[s:(s + n_win - 1)]
        seg_fit = zeros(n_win)
        for it in 1:maxiter
            buf .= 0
            buf[1:n_win, :] .= V .* r
            J = plan * buf
            found = false
            for (j, bj) in enumerate(bins)
                isempty(bj) && continue
                # Thomson F statistic of a sinusoid at each scanned frequency
                best_f, best_b, best_mu = -Inf, 0, 0.0im
                for b in bj
                    Jb = @view J[b, :]
                    mu = sum(Jb .* V0) / sV0
                    F = (k - 1) * abs2(mu) * sV0 / sum(abs2, Jb .- mu .* V0)
                    F > best_f && ((best_f, best_b, best_mu) = (F, b, mu))
                end
                best_f > f_crit || continue
                s_fit = 2 .* real.(best_mu .* cis.(2π * f[best_b] .* t))
                r .-= s_fit
                seg_fit .+= s_fit
                found = true
                it == 1 && (nsig[j] += 1)
            end
            found || break
        end
        fit[s:(s + n_win - 1)] .+= w .* seg_fit
        wsum[s:(s + n_win - 1)] .+= w
    end

    return x .- fit ./ wsum, nsig ./ length(starts)
end
