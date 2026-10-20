export remove_powerline
export remove_powerline!
export detect_powerline
export detect_powerline!

# power line frequencies in use worldwide
const _LINE_FREQUENCIES = (50, 60)

# zero-phase IIR notch filter of a signal vector (as filter(obj; fprototype = :iirnotch), two-pass)
_notch(x::AbstractVector, fs::Int64, frq::Real, bw::Real)::Vector{Float64} =
    _filter_apply(x; flt = filter_create(; fprototype = :iirnotch, cutoff = frq, fs = fs, bw = bw))

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

For each selected channel:

1. the power spectrum (Welch, dB) is checked for a peak (prominence ≥ `pr` dB) at the power line frequency `pl_frq` (within `pl_frq ± d`); channels without a peak are not modified
2. zero-phase IIR notch filters of bandwidths `q`, `2q`, …, 10 Hz are applied to the signal and the bandwidth giving the flattest spectrum (lowest variance of the PSD within `pl_frq ± d`) is used
3. the procedure is repeated for each harmonic (`2 × pl_frq`, `3 × pl_frq`, … below the Nyquist frequency) that shows a peak in this channel (bandwidths up to 5 Hz, window `± d / 2`)

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be continuous (1 epoch)
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pl_frq::Real=obj.header.recording[:line_frequency]`: power line frequency in Hz; default is read from the object header
- `method::Symbol=:iir`: filtering method (currently only `:iir` is supported)
- `pr::Real=2.0`: minimum peak prominence in dB
- `d::Real=5.0`: half-width (Hz) of the frequency window around the power line frequency for the peak detection and the bandwidth optimization (`d / 2` for harmonics)
- `q::Real=0.1`: bandwidth optimization step size in Hz; must be in `[0.01, 5)`

# Returns

- `NeuroAnalyzer.NEURO`: new object with power line noise removed
- `DataFrame`: filtered peaks: `:channel` (label), `:frequency` (Hz), `:bandwidth` (optimized notch bandwidth, Hz)
"""
function remove_powerline(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pl_frq::Real = obj.header.recording[:line_frequency],
    method::Symbol = :iir,
    pr::Real = 2.0,
    d::Real = 5.0,
    q::Real = 0.1,
)::Tuple{NeuroAnalyzer.NEURO, DataFrame}
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("remove_powerline() requires a continuous (1-epoch) object."))
    fs = sr(obj)
    pl_frq > 0 || throw(ArgumentError("pl_frq must be > 0."))
    pl_frq < fs / 2 || throw(ArgumentError("pl_frq must be < $(fs / 2) Hz (Nyquist)."))
    q >= 0.01 || throw(ArgumentError("q must be ≥ 0.01."))
    q < 5 || throw(ArgumentError("q must be < 5."))
    d > 0 || throw(ArgumentError("d must be > 0."))
    _check_var(method, [:iir], "method")

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    clabels = labels(obj)

    # create new dataset
    obj_tmp = deepcopy(obj)

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
        progbar = Progress(length(ch); dt = 1, barlen = 20, color = :white, enabled = progress_bar)

        for ch_idx in ch
            x = obj_tmp.data[ch_idx, :, 1]
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
                bws = Base.filter(b -> h - b / 2 > 0 && h + b / 2 < fs / 2, h_idx == 1 ? bw_fund : bw_harm)
                isempty(bws) && continue
                v = _notch_band_var(x, fs, h, bws, f1, f2)
                bw = bws[argmin(v)]
                x = _notch(x, fs, h, bw)
                push!(df, (clabels[ch_idx], h, bw))
            end
            obj_tmp.data[ch_idx, :, 1] = x

            progress_bar && next!(progbar)
        end
        progress_bar && finish!(progbar)
    finally
        NeuroAnalyzer.verbose = verbose_tmp
    end

    isempty(no_peak) || _warn(
        "No power line peak at $pl_frq Hz in channel(s) $(join(no_peak, ", ")); power line frequency not filtered (check pl_frq or whether the signal is already filtered).",
    )
    push!(obj_tmp.history, "remove_powerline(OBJ, ch=$ch, pl_frq=$pl_frq, method=:$method, pr=$pr, d=$d, q=$q)")

    return obj_tmp, df
end

"""
    remove_powerline!(obj; <keyword arguments>)

Remove power line noise and its harmonics in-place (see `remove_powerline()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pl_frq::Real=obj.header.recording[:line_frequency]`: power line frequency in Hz; default is read from the object header
- `method::Symbol=:iir`: filtering method (currently only `:iir` is supported)
- `pr::Real=2.0`: minimum peak prominence in dB
- `d::Real=5.0`: half-width (Hz) of the frequency window around the power line frequency for the peak detection and the bandwidth optimization (`d / 2` for harmonics)
- `q::Real=0.1`: bandwidth optimization step size in Hz; must be in `[0.01, 5)`

# Returns

- `DataFrame`: filtered peaks: `:channel` (label), `:frequency` (Hz), `:bandwidth` (optimized notch bandwidth, Hz)
"""
function remove_powerline!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pl_frq::Real = obj.header.recording[:line_frequency],
    method::Symbol = :iir,
    pr::Real = 2.0,
    d::Real = 5.0,
    q::Real = 0.1,
)::DataFrame
    obj_tmp, df = remove_powerline(obj; ch = ch, pl_frq = pl_frq, method = method, pr = pr, d = d, q = q)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history

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
        throw(ArgumentError("Sampling rate must be > $(2 * (minimum(_LINE_FREQUENCIES) + 5)) Hz to detect power line noise."))

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
        noise_frq[ch_idx, ep_idx] = _detect_powerline(@view(obj.data[ch_idx, :, ep_idx]); fs = sr(obj), pr = pr)

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
        _warn("No power line noise detected; line frequency not changed ($(obj.header.recording[:line_frequency]) Hz).")
    else
        _info("Power line noise detected at $noise_frq Hz")
        obj.header.recording[:line_frequency] = noise_frq
    end

    return nothing
end
