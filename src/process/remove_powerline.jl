export remove_powerline
export remove_powerline!
export detect_powerline
export detect_powerline!

# power line frequencies in use worldwide
const _LINE_FREQUENCIES = (50, 60)

# zero-phase IIR notch filter of a signal vector (as filter(obj; fprototype = :iirnotch), two-pass)
_notch(x::AbstractVector, fs::Int64, frq::Real, bw::Real)::Vector{Float64} =
    _filter_apply(
        x;
        flt = filter_create(; fprototype = :iirnotch, cutoff = frq, fs = fs, bw = bw),
    )

# for every notch bandwidth in `bws`: variance of the PSD (dB) within [f1, f2] of the notch-filtered
# signal (the candidates are independent and evaluated in parallel)
function _notch_band_var(
    x::AbstractVector,
    fs::Int64,
    frq::Real,
    bws::AbstractVector{<:Real},
    f1::Real,
    f2::Real,
)::Vector{Float64}
    v = zeros(length(bws))
    Threads.@threads :static for idx in eachindex(bws)
        p, f = _psd(_notch(x, fs, frq, bws[idx]); fs = fs, db = true)
        v[idx] = var(@view(p[vsearch(f1, f):vsearch(f2, f)]))
    end
    return v
end

# prominence (dB) of the PSD at f0 over the median PSD 2–5 Hz away from f0
function _line_prominence(p::AbstractVector, f::AbstractVector, f0::Real)::Float64
    nb = findall(x -> 2 <= abs(x - f0) <= 5, f)
    pk = findall(x -> abs(x - f0) <= 0.5, f)
    (isempty(nb) || isempty(pk)) && return -Inf
    return maximum(p[pk]) - median(p[nb])
end

"""
    remove_powerline(obj; <keyword arguments>)

Remove power line noise and its harmonics from a continuous NEURO object.

`method = :iir` (notch filters), for each selected channel:

1. the power spectrum (Welch, dB) is checked for a peak (prominence ≥ `pr` dB) at the power line frequency `pl_frq` (within `pl_frq ± d`); channels without a peak are not modified
2. zero-phase IIR notch filters of bandwidths `q`, `2q`, …, 10 Hz are applied to the signal and the bandwidth giving the flattest spectrum (lowest variance of the PSD within `pl_frq ± d`) is used
3. the procedure is repeated for each harmonic (`2 × pl_frq`, `3 × pl_frq`, … below the Nyquist frequency) that shows a peak in this channel (bandwidths up to 5 Hz, window `± d / 2`)

`method = :zapline` (Zapline-plus): spatial filtering that removes only the line noise components of the selected channels together (at least 2 channels, e.g. all EEG channels), so that the signal at the line frequency and the spectrum around it are preserved, unlike notch filtering. For each chunk of the signal:

1. the signal is smoothed with a moving average of one period of the line frequency (`sr / pl_frq` samples), which cancels the line frequency and all its harmonics
2. the residual (signal minus smoothed signal) contains the line noise and high-frequency activity; denoising source separation (DSS) finds its spatial components ordered by the ratio of their power at the line frequency and its harmonics to their total power
3. the first `nremove` components are regressed out of the residual, and the residual is added back to the smoothed signal

with the Zapline-plus additions:

- `nremove = :auto`: the number of components is detected per chunk by iterative outlier detection of the DSS scores (mean + `sigma` × SD); then components are added while the line peak remains in the cleaned chunk (prominence of the channel-averaged spectrum ≥ `pr` dB); at least 1 and at most `max_prop` × number of channels (if the peak remains, e.g. in heavily contaminated recordings, a warning suggests a `max_prop` value that would remove it, or `method = :iir`)
- `chunk = :adaptive`: the signal is split where the spatial pattern of the line noise changes (peaks of the distance between the covariance matrices of consecutive 1-second segments of the signal band-passed at `pl_frq` ± 3 Hz), chunks at least `min_chunk` seconds long; `:fixed`: chunks of `chunk_len` seconds; `:none`: the whole signal

`method = :cleanline` (CleanLine): estimation and subtraction of the line sinusoids, channel by channel (any number of channels), in sliding windows (`wlen` seconds, step `wstep`):

1. the window is multiplied by Slepian tapers (time-bandwidth `wlen × bw / 2`) and Fourier transformed
2. for the line frequency and each harmonic, the frequency within ± `bw / 2` with the largest Thomson F statistic (a sinusoid in the window) is found; if the sinusoid is significant (p < `p`, Bonferroni-corrected for the `bw × wlen` independent frequencies scanned), its amplitude and phase are estimated and it is subtracted from the window; repeated (up to `maxiter` times) until no significant sinusoid remains
3. the fitted sinusoids of the overlapping windows are combined with smooth (Hann) weights and subtracted from the signal

The rest of the spectrum is not changed; non-stationary noise (amplitude or frequency changes within a window) is removed incompletely.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be continuous (1 epoch)
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pl_frq::Real=obj.header.recording[:line_frequency]`: power line frequency in Hz; default is read from the object header
- `method::Symbol=:iir`: method: `:iir` (notch filters), `:zapline` (Zapline-plus) or `:cleanline` (CleanLine)
- `pr::Real=2.0`: minimum peak prominence in dB (`:iir`: peaks to filter; `:zapline` with `nremove = :auto`: components are added while the remaining line peak is ≥ `pr` dB)
- `d::Real=5.0`: half-width (Hz) of the frequency window around the power line frequency for the peak detection and the bandwidth optimization (`d / 2` for harmonics) (`:iir`)
- `q::Real=0.1`: bandwidth optimization step size in Hz; must be in `[0.01, 5)` (`:iir`)
- `nremove::Union{Int64, Symbol}=:auto`: number of components to remove in each chunk, or `:auto` (`:zapline`)
- `chunk::Symbol=:adaptive`: chunking: `:adaptive`, `:fixed` or `:none` (`:zapline`)
- `chunk_len::Real=150`: chunk length in seconds for `chunk = :fixed` (`:zapline`)
- `min_chunk::Real=30`: minimum chunk length in seconds for `chunk = :adaptive` (`:zapline`)
- `sigma::Real=3`: outlier threshold (number of SDs) for `nremove = :auto` (`:zapline`)
- `max_prop::Real=0.2`: maximum proportion of components removed for `nremove = :auto` (`:zapline`)
- `nfft::Int64=1024`: FFT length for the DSS bias (the power at the line frequency and its harmonics) (`:zapline`)
- `p::Real=0.01`: significance level of the Thomson F-test, corrected for the number of scanned frequencies (`:cleanline`)
- `wlen::Real=4`: window length in seconds (`:cleanline`)
- `wstep::Real=1`: window step in seconds (`:cleanline`)
- `bw::Real=2`: bandwidth (Hz) of the Slepian tapers; the sinusoids are searched within ± `bw / 2` of the line frequency and its harmonics (`:cleanline`)
- `maxiter::Int64=10`: maximum number of fits in each window (`:cleanline`)

# Returns

- `NeuroAnalyzer.NEURO`: new object with power line noise removed
- `DataFrame`:
    - `:iir`: filtered peaks: `:channel` (label), `:frequency` (Hz), `:bandwidth` (optimized notch bandwidth, Hz)
    - `:zapline`: one row per chunk: `:start`, `:stop` (seconds), `:nremove` (number of removed components), `:reduction` (reduction of the power at the line frequency and its harmonics, dB), `:prominence` (remaining line peak, dB)
    - `:cleanline`: one row per channel and line frequency: `:channel` (label), `:frequency` (Hz), `:windows` (proportion of windows with a significant sinusoid), `:reduction` (reduction of the power at the frequency, dB)

# References

1. de Cheveigné A. ZapLine: A simple and effective method to remove power line artifacts. NeuroImage. 2020;207:116356.
2. Klug M, Kloosterman NA. Zapline-plus: A Zapline extension for automatic and adaptive removal of frequency-specific noise artifacts in M/EEG. Hum Brain Mapp. 2022;43(9):2743-2758.
3. Mullen T. CleanLine EEGLAB plugin. San Diego, CA: Neuroimaging Informatics Tools and Resources Clearinghouse (NITRC); 2012.
"""
function remove_powerline(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pl_frq::Real = obj.header.recording[:line_frequency],
    method::Symbol = :iir,
    pr::Real = 2.0,
    d::Real = 5.0,
    q::Real = 0.1,
    nremove::Union{Int64, Symbol} = :auto,
    chunk::Symbol = :adaptive,
    chunk_len::Real = 150,
    min_chunk::Real = 30,
    sigma::Real = 3,
    max_prop::Real = 0.2,
    nfft::Int64 = 1024,
    p::Real = 0.01,
    wlen::Real = 4,
    wstep::Real = 1,
    bw::Real = 2,
    maxiter::Int64 = 10,
)::Tuple{NeuroAnalyzer.NEURO, DataFrame}
    obj_new = deepcopy(obj)
    df = remove_powerline!(
        obj_new;
        ch = ch,
        pl_frq = pl_frq,
        method = method,
        pr = pr,
        d = d,
        q = q,
        nremove = nremove,
        chunk = chunk,
        chunk_len = chunk_len,
        min_chunk = min_chunk,
        sigma = sigma,
        max_prop = max_prop,
        nfft = nfft,
        p = p,
        wlen = wlen,
        wstep = wstep,
        bw = bw,
        maxiter = maxiter,
    )

    return obj_new, df
end

"""
    remove_powerline!(obj; <keyword arguments>)

Remove power line noise and its harmonics in-place (see `remove_powerline()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pl_frq::Real=obj.header.recording[:line_frequency]`: power line frequency in Hz; default is read from the object header
- `method::Symbol=:iir`: method: `:iir` (notch filters), `:zapline` (Zapline-plus) or `:cleanline` (CleanLine)
- `pr::Real=2.0`: minimum peak prominence in dB (`:iir`: peaks to filter; `:zapline` with `nremove = :auto`: components are added while the remaining line peak is ≥ `pr` dB)
- `d::Real=5.0`: half-width (Hz) of the frequency window around the power line frequency for the peak detection and the bandwidth optimization (`d / 2` for harmonics) (`:iir`)
- `q::Real=0.1`: bandwidth optimization step size in Hz; must be in `[0.01, 5)` (`:iir`)
- `nremove::Union{Int64, Symbol}=:auto`: number of components to remove in each chunk, or `:auto` (`:zapline`)
- `chunk::Symbol=:adaptive`: chunking: `:adaptive`, `:fixed` or `:none` (`:zapline`)
- `chunk_len::Real=150`: chunk length in seconds for `chunk = :fixed` (`:zapline`)
- `min_chunk::Real=30`: minimum chunk length in seconds for `chunk = :adaptive` (`:zapline`)
- `sigma::Real=3`: outlier threshold (number of SDs) for `nremove = :auto` (`:zapline`)
- `max_prop::Real=0.2`: maximum proportion of components removed for `nremove = :auto` (`:zapline`)
- `nfft::Int64=1024`: FFT length for the DSS bias (`:zapline`)
- `p::Real=0.01`: significance level of the Thomson F-test, corrected for the number of scanned frequencies (`:cleanline`)
- `wlen::Real=4`: window length in seconds (`:cleanline`)
- `wstep::Real=1`: window step in seconds (`:cleanline`)
- `bw::Real=2`: bandwidth (Hz) of the Slepian tapers; the sinusoids are searched within ± `bw / 2` of the line frequency and its harmonics (`:cleanline`)
- `maxiter::Int64=10`: maximum number of fits in each window (`:cleanline`)

# Returns

- `DataFrame`:
    - `:iir`: filtered peaks: `:channel` (label), `:frequency` (Hz), `:bandwidth` (optimized notch bandwidth, Hz)
    - `:zapline`: one row per chunk: `:start`, `:stop` (seconds), `:nremove` (number of removed components), `:reduction` (dB), `:prominence` (dB)
    - `:cleanline`: one row per channel and line frequency: `:channel`, `:frequency` (Hz), `:windows` (proportion of windows with a significant sinusoid), `:reduction` (dB)
"""
function remove_powerline!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pl_frq::Real = obj.header.recording[:line_frequency],
    method::Symbol = :iir,
    pr::Real = 2.0,
    d::Real = 5.0,
    q::Real = 0.1,
    nremove::Union{Int64, Symbol} = :auto,
    chunk::Symbol = :adaptive,
    chunk_len::Real = 150,
    min_chunk::Real = 30,
    sigma::Real = 3,
    max_prop::Real = 0.2,
    nfft::Int64 = 1024,
    p::Real = 0.01,
    wlen::Real = 4,
    wstep::Real = 1,
    bw::Real = 2,
    maxiter::Int64 = 10,
)::DataFrame
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("remove_powerline() requires a continuous (1-epoch) object."))
    fs = sr(obj)
    pl_frq > 0 || throw(ArgumentError("pl_frq must be > 0."))
    pl_frq < fs / 2 || throw(ArgumentError("pl_frq must be < $(fs / 2) Hz (Nyquist)."))
    _check_var(method, [:iir, :zapline, :cleanline], "method")

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    if method === :zapline
        df = _zapline!(
            obj;
            ch = ch,
            fline = pl_frq,
            nremove = nremove,
            chunk = chunk,
            chunk_len = chunk_len,
            min_chunk = min_chunk,
            sigma = sigma,
            max_prop = max_prop,
            nfft = nfft,
            pr = pr,
        )
        _history!(
            obj, "remove_powerline",
            (; ch, pl_frq, method, pr, nremove, chunk, chunk_len, min_chunk, sigma, max_prop, nfft),
        )
        return df
    end

    if method === :cleanline
        df = _cleanline!(
            obj;
            ch = ch,
            fline = pl_frq,
            p = p,
            wlen = wlen,
            wstep = wstep,
            bw = bw,
            maxiter = maxiter,
        )
        _history!(obj, "remove_powerline", (; ch, pl_frq, method, p, wlen, wstep, bw, maxiter))
        return df
    end

    q >= 0.01 || throw(ArgumentError("q must be ≥ 0.01."))
    q < 5 || throw(ArgumentError("q must be < 5."))
    d > 0 || throw(ArgumentError("d must be > 0."))
    clabels = labels(obj)

    # the fundamental and its harmonics below the Nyquist frequency
    frqs = [k * pl_frq for k = 1:floor(Int64, fs / 2 / pl_frq) if k * pl_frq < fs / 2]
    bw_fund = collect(q:q:10.0)
    bw_harm = collect(q:q:5.0)

    df = DataFrame(:channel => String[], :frequency => Float64[], :bandwidth => Float64[])
    no_peak = String[]

    _info("Removing power line noise at $pl_frq Hz and its harmonics")

    # messages of the internal filtering are not shown; verbose is restored even if an error occurs
    verbose_tmp = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false
    try
        progbar = Progress(
            length(ch);
            dt = 1,
            barlen = 20,
            color = :white,
            enabled = progress_bar,
        )

        for ch_idx in ch
            x = obj.data[ch_idx, :, 1]
            for (h_idx, h) in enumerate(frqs)
                hw = h_idx == 1 ? d : d / 2
                f1, f2 = max(h - hw, 0), min(h + hw, fs / 2)

                # is there a peak at this frequency (within 1 Hz) in this channel?
                p, f = _psd(x; fs = fs, db = true)
                band = vsearch(f1, f):vsearch(f2, f)
                pks, _ = findpeaks1d(p[band]; prominence = pr)
                if !any(i -> abs(f[band][i] - h) <= 1, pks)
                    h_idx == 1 && push!(no_peak, clabels[ch_idx])
                    continue
                end

                # the notch band must lie within (0, Nyquist)
                bws = Base.filter(
                    b -> h - b / 2 > 0 && h + b / 2 < fs / 2,
                    h_idx == 1 ? bw_fund : bw_harm,
                )
                isempty(bws) && continue
                v = _notch_band_var(x, fs, h, bws, f1, f2)
                bw = bws[argmin(v)]
                x = _notch(x, fs, h, bw)
                push!(df, (clabels[ch_idx], h, bw))
                _notch_add!(obj, [ch_idx], h)
            end
            obj.data[ch_idx, :, 1] = x

            progress_bar && next!(progbar)
        end
        progress_bar && finish!(progbar)
    finally
        NeuroAnalyzer.verbose = verbose_tmp
    end

    isempty(no_peak) || _warn(
        "No power line peak at $pl_frq Hz in channel(s) $(join(no_peak, ", ")); power line frequency not filtered (check pl_frq or whether the signal is already filtered).",
    )
    _history!(obj, "remove_powerline", (; ch, pl_frq, method, pr, d, q))

    return df
end

"""
    _detect_powerline(s; <keyword arguments>)

Detect the power line noise frequency (50 or 60 Hz) in a signal.

The prominence of the power spectrum (Welch, dB) at 50 and 60 Hz over the median power 2–5 Hz away is calculated; the frequency with the larger prominence is returned if the prominence is ≥ `pr` dB.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 2
- `pr::Real=3`: minimum prominence in dB

# Returns

- `Float64`: power line frequency in Hz (50 or 60), 0 if not detected
"""
function _detect_powerline(s::AbstractVector; fs::Int64, pr::Real = 3)::Float64
    # validate
    fs >= 2 || throw(ArgumentError("fs must be ≥ 2."))

    p, f = _psd(s; fs = fs, db = true)
    best, best_pr = 0.0, -Inf
    for f0 in _LINE_FREQUENCIES
        f0 + 5 < fs / 2 || continue
        pr_f0 = _line_prominence(p, f, f0)
        pr_f0 > best_pr && ((best, best_pr) = (Float64(f0), pr_f0))
    end

    return best_pr >= pr ? best : 0.0
end

"""
    detect_powerline(obj; <keyword arguments>)

Detect the power line noise frequency (50 or 60 Hz).

For each channel and epoch, the prominence of the power spectrum at 50 and 60 Hz over the neighbouring frequencies (2–5 Hz away) is calculated; the frequency detected in most channels/epochs is returned.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `pr::Real=3`: minimum prominence in dB

# Returns

- `Float64`: power line frequency in Hz (50 or 60), 0 if not detected
"""
function detect_powerline(obj::NeuroAnalyzer.NEURO; pr::Real = 3)::Float64
    sr(obj) > 2 * (minimum(_LINE_FREQUENCIES) + 5) ||
        throw(
            ArgumentError(
                "Sampling rate must be > $(2 * (minimum(_LINE_FREQUENCIES) + 5)) Hz to detect power line noise.",
            ),
        )

    # number of channels
    ch_n = size(obj, 1)
    # number of epochs
    ep_n = size(obj, 3)

    # pre-allocate output
    noise_frq = zeros(ch_n, ep_n)

    # initialize progress bar
    progbar =
        Progress(ch_n * ep_n; dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        noise_frq[ch_idx, ep_idx] =
            _detect_powerline(@view(obj.data[ch_idx, :, ep_idx]); fs = sr(obj), pr = pr)

        # update progress bar
        progress_bar && next!(progbar)
    end
    progress_bar && finish!(progbar)

    # the frequency detected most often
    n = [count(==(f0), noise_frq) for f0 in _LINE_FREQUENCIES]
    maximum(n) == 0 && return 0.0

    return Float64(_LINE_FREQUENCIES[argmax(n)])
end

"""
    detect_powerline!(obj; <keyword arguments>)

Detect the power line noise frequency (see `detect_powerline()`) and store it in the object header (`obj.header.recording[:line_frequency]`); the header is not changed if no power line noise is detected.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `pr::Real=3`: minimum prominence in dB

# Returns

- `Nothing`
"""
function detect_powerline!(obj::NeuroAnalyzer.NEURO; pr::Real = 3)::Nothing
    noise_frq = detect_powerline(obj; pr = pr)
    if noise_frq == 0
        _warn(
            "No power line noise detected; line frequency not changed ($(obj.header.recording[:line_frequency]) Hz).",
        )
    else
        _info("Power line noise detected at $noise_frq Hz")
        obj.header.recording[:line_frequency] = noise_frq
    end
    _history!(obj, "detect_powerline!", (; pr))

    return nothing
end
