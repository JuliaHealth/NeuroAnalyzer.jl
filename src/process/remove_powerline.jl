export remove_powerline
export remove_powerline!
export detect_powerline
export detect_powerline!

"""
    remove_powerline(obj; <keyword arguments>)

Remove power line noise and its harmonics from a continuous NEURO object.

For each selected channel, an IIR notch filter is optimised by scanning a range of bandwidth values and choosing the one that minimises spectral variance around the noise frequency. The same process is repeated for any harmonic peaks detected above the fundamental.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be continuous (1 epoch)
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pl_frq::Real=obj.header.recording[:line_frequency]`: power line frequency in Hz; default is read from the object header
- `method::Symbol=:iir`: filtering method (currently only `:iir` is supported)
- `pr::Real=2.0`: minimum peak prominence in dB for harmonic detection
- `d::Real=5.0`: minimum distance between detected peaks in Hz
- `q::Real=0.1`: bandwidth optimization step size in Hz; must be in `[0.01, 5)`

# Returns

- `NeuroAnalyzer.NEURO`: new object with power line noise removed
- `DataFrame`: detected peaks with their optimized notch bandwidths

# Throws

- `ArgumentError`: if the object has more than one epoch, `pl_frq` is out of range, `q` is out of range, or no power line peak is found near `pl_frq`
"""
function remove_powerline(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pl_frq::Real = obj.header.recording[:line_frequency],
    method::Symbol = :iir,
    pr::Real = 2.0,
    d::Real = 5.0,
    q::Real = 0.1
)::Tuple{NeuroAnalyzer.NEURO, DataFrame}

    nepochs(obj) == 1 || throw(ArgumentError("remove_powerline() requires a continuous (1-epoch) object."))
    pl_frq >= 0 || throw(ArgumentError("pl_frq must be ≥ 0."))
    pl_frq <= sr(obj)/2 || throw(ArgumentError("pl_frq must be ≤ $(sr(obj)/2) Hz (Nyquist)."))
    q >= 0.01 || throw(ArgumentError("q must be ≥ 0.01."))
    q < 5 || throw(ArgumentError("q must be < 5."))

    _check_var(method, [:iir], "method")

    ch_idx_vec = get_channel(obj, ch=ch)
    clabels = labels(obj)

    obj_new = deepcopy(obj)
    pl_best_bw = Float64[]
    pks_frq = Float64[]
    pks_best_bw = Vector{Float64}[]

    verbose_tmp = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false

    _info("Removing power line noise at $pl_frq Hz and its harmonics")

    if method === :iir
        bw_values = collect(q:q:10.0)

        progbar   = Progress(length(ch_idx_vec), dt=1, barlen=20, color=:white, enabled=progress_bar)

        # use enumerate so we index pl_best_bw by loop position
        for (loop_idx, ch_idx) in enumerate(ch_idx_vec)
            ch_label = clabels[ch_idx]

            # --- detect and validate the fundamental power line peak ---
            psd_data = psd(obj_new, ch=ch_label, db=true)
            f = psd_data.f
            p = psd_data.p[:]
            f_band = vsearch(pl_frq - d, f):vsearch(pl_frq + d, f)
            p_band = p[f_band]
            f_band_vals = f[f_band]
            peaks, _ = findpeaks1d(p_band; prominence=pr)
            length(peaks) > 0 || throw(ArgumentError(
                "No power line peak found near $pl_frq Hz in channel $ch_label. Check pl_frq or whether the signal is already filtered."))
            pl_amp = vsearch(maximum(p_band[peaks]), p_band)
            pl_frq_detected = f_band_vals[pl_amp]
            (pl_frq_detected < pl_frq - d || pl_frq_detected > pl_frq + d) &&
                throw(ArgumentError("Channel $ch_label: power line peak at $pl_frq_detected Hz is outside the expected range."))

            # --- optimize notch bandwidth for the fundamental ---
            v = [begin
                    obj_tmp = NeuroAnalyzer.filter(obj_new, ch=ch_label, fprototype=:iirnotch, cutoff=pl_frq, bw=bw)
                    p2, f2 = psd(obj_tmp; ch=ch_label, db=true)
                    seg = p2[vsearch(pl_frq - d, f2):vsearch(pl_frq + d, f2)]
                    var(seg)
                 end for bw in bw_values]
            best_bw = bw_values[vsearch(minimum(v), v)]
            push!(pl_best_bw, best_bw)

            NeuroAnalyzer.filter!(obj_new, ch=ch_label, fprototype=:iirnotch, cutoff=pl_frq, bw=pl_best_bw[loop_idx])

            # --- detect harmonics above the fundamental ---
            p, f = psd(obj_new, ch=ch_label, db=true)
            p = p[:]
            f_above = vsearch(2 * pl_frq - 2 * d, f)
            p_above = p[f_above:end]
            f_above_vals = f[f_above:end]

            if isempty(pks_frq)
                pr_tmp, d_tmp = pr, d
                peaks, _ = findpeaks1d(p_above, prominence=pr_tmp, distance=vsearch(d_tmp, f))
                while isempty(peaks)
                    pr_tmp -= 0.1
                    d_tmp  -= 0.1
                    peaks, _ = findpeaks1d(p_above, prominence=pr_tmp, distance=vsearch(d_tmp, f))
                end
                for idx in eachindex(peaks)
                    push!(pks_frq, f_above_vals[peaks[idx]])
                end
            end

            n = length(pks_frq)
            if n > 0
                bw_harm  = collect(q:q:5.0)
                best_bw_harm = zeros(n)
                for peak_idx in eachindex(pks_frq)
                    vh = [begin
                             obj_tmp = NeuroAnalyzer.filter(
                                 obj_new, ch=ch_label, fprototype=:iirnotch,
                                 cutoff=pks_frq[peak_idx], bw=bw)
                             p2, f2 = psd(obj_tmp, ch=ch_label, db=true)
                             seg = p2[vsearch(pks_frq[peak_idx] - d/2, f2):vsearch(pks_frq[peak_idx] + d/2, f2)]
                             var(seg)
                          end for bw in bw_harm]
                    best_bw_harm[peak_idx] = bw_harm[vsearch(minimum(vh), vh)]
                end
                for peak_idx in eachindex(pks_frq)
                    NeuroAnalyzer.filter!(obj_new, ch=ch_label, fprototype=:iirnotch,
                                          cutoff=pks_frq[peak_idx], bw=best_bw_harm[peak_idx])
                end
                push!(pks_best_bw, best_bw_harm)
            end

            progress_bar && next!(progbar)
        end

        # --- build results DataFrame ---
        df = DataFrame("channel" => clabels[ch_idx_vec], "power line bandwidth" => pl_best_bw)
        if !isempty(pks_frq)
            n_ch = length(ch_idx_vec)
            n_pk = length(pks_frq)
            pks_frq_m    = zeros(n_ch, n_pk)
            pks_bw_m     = zeros(n_ch, n_pk)
            peak_frq_names = ["peak $i frequency"  for i in 1:n_pk]
            peak_bw_names  = ["peak $i bandwidth"  for i in 1:n_pk]
            for i in 1:n_ch
                pks_frq_m[i, :] = pks_frq
                pks_bw_m[i, :]  = pks_best_bw[i]
            end
            df = hcat(df, DataFrame(pks_frq_m, peak_frq_names))
            df = hcat(df, DataFrame(pks_bw_m,  peak_bw_names))
        end
    end

    NeuroAnalyzer.verbose = verbose_tmp
    push!(obj_new.history, "remove_powerline(OBJ, pl_frq=$pl_frq, method=:$method, pr=$pr, d=$d, q=$q)")

    return obj_new, df

end

"""
    remove_powerline!(obj; <keyword arguments>)

Remove power line noise in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pl_frq::Real=obj.header.recording[:line_frequency]`: power line frequency in Hz; default is read from the object header
- `method::Symbol=:iir`: filtering method (currently only `:iir` is supported)
- `pr::Real=2.0`: minimum peak prominence in dB for harmonic detection
- `d::Real=5.0`: minimum distance between detected peaks in Hz
- `q::Real=0.1`: bandwidth optimization step size in Hz; must be in `[0.01, 5)`

# Returns

- `DataFrame`: detected peaks with their optimized notch bandwidths

# Throws

- `ArgumentError`: if the object has more than one epoch, `pl_frq` is out of range, `q` is out of range, or no power line peak is found near `pl_frq`
"""
function remove_powerline!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pl_frq::Real = obj.header.recording[:line_frequency],
    method::Symbol = :iir,
    pr::Real = 2.0,
    d::Real = 5.0,
    q::Real = 0.1
)::DataFrame

    obj_new, df = remove_powerline(
        obj,
        ch = ch,
        pl_frq = pl_frq,
        method = method,
        pr = pr,
        d = d,
        q = q
    )
    obj.data = obj_new.data
    obj.history = obj_new.history

    return df

end

"""
    detect_powerline(s; <keyword arguments>)

Detect the dominant power line noise frequency in a signal.

Fits a sine + cosine model at each integer frequency from 1 to `fs÷2` Hz and returns the frequency with maximum fitted power.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 2

# Returns

- `Float64`: dominant noise frequency in Hz

# Throws

- `ArgumentError`: if `fs < 2`
"""
function detect_powerline(s::AbstractVector; fs::Int64)::Float64

    fs >= 2 || throw(ArgumentError("fs must be ≥ 2."))

    n = length(s)
    t = range(0, n / fs, length=n)

    noise_power = zeros(fs ÷ 2)
    for freq in 1:(fs ÷ 2)
        df = DataFrame(
            :signal => s,
            :b0 => ones(n),
            :b1 => sin.(2π * freq .* t),
            :b2 => cos.(2π * freq .* t),
        )
        lr = GLM.lm(@formula(signal ~ b0 + b1 + b2), df)
        b  = GLM.coef(lr)
        noise_power[freq] = b[3]^2 + b[4]^2
    end

    return Float64(vsearch(maximum(noise_power), noise_power))

end

"""
    detect_powerline(obj)

Detect the dominant power line noise frequency for each channel and epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Array{Float64, 2}`: noise frequency in Hz; shape `(channels,  epochs)`
"""
function detect_powerline(obj::NeuroAnalyzer.NEURO)::Array{Float64, 2}

    ch_n = size(obj, 1)
    ep_n = size(obj, 3)

    noise_frq = zeros(ch_n, ep_n)

    # initialize progress bar
    progbar = Progress(ch_n * ep_n, dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        noise_frq[ch_idx, ep_idx] = detect_powerline(@view(obj.data[ch_idx, :, ep_idx]); fs=sr(obj))

        # update progress bar
        progress_bar && next!(progbar)
    end

    return noise_frq

end

"""
    detect_powerline!(obj)

Detect power line noise and store the median detected frequency in the object header in-place.
# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; `header.recording[:line_frequency]` is updated in-place

# Returns

- `Nothing`
"""
function detect_powerline!(obj::NeuroAnalyzer.NEURO)::Nothing

    noise_frq = detect_powerline(obj)
    med_frq = median(noise_frq)
    _info("Power line noise detected at $med_frq Hz")
    obj.header.recording[:line_frequency] = med_frq

    return nothing

end
