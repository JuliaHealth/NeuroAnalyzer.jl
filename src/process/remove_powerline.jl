export remove_powerline
export remove_powerline!
export detect_powerline
export detect_powerline!

"""
    remove_powerline(obj; <keyword arguments>)

Remove power line noise and its peaks above power line frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `pl_frq::Real=obj.header.recording[:line_frequency]`: power line frequency, default is read from the OBJ header
- `method::Symbol=:iir`:
    - `:iir`: use IIR filter
- `pr::Real=2.0`: prominence of noise peaks in dB
- `d::Real=5.0`: minimum distance between peaks in Hz
- `q::Real=0.1`: optimization step size

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
- `df::DataFrame`: list of peaks detected
"""
function remove_powerline(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pl_frq::Real=obj.header.recording[:line_frequency], method::Symbol=:iir, pr::Real=2.0, d::Float64=5.0, q::Real=0.1)::Tuple{NeuroAnalyzer.NEURO, DataFrame}

    @assert nepochs(obj) == 1 "remove_powerline() must be applied to a continuous signal."
    @assert pl_frq >= 0 "pl_freq must be ≥ 0."
    @assert pl_frq <= sr(obj) / 2 "pl_freq must be ≤ $(sr(obj) / 2)."
    @assert q >= 0.01 "q must be ≥ 0.01."
    @assert q < 5 "q must be < 5."

    ch = get_channel(obj, ch=ch)
    clabels = labels(obj)

    _check_var(method, [:iir], "method")

    obj_new = deepcopy(obj)

    _info("Removing power line noise at $pl_frq Hz and its harmonics")

    verbose_tmp = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false

    # number of peaks
    pl_best_bw = Float64[]
    pks_frq = Float64[]
    pks_best_bw = Vector{Float64}[]

    if method === :iir
        bw_values = collect(0:q:20.0)[2:end]

        # initialize progress bar
        progress_bar && (progbar = Progress(length(ch), dt=1, barlen=20, color=:white))

        @inbounds for ch_idx in ch
            # detect power line peak
            p, f = psd(obj_new, ch=clabels[ch_idx], db=true)
            p = p[:]
            f_pl = vsearch(pl_frq - d, f):vsearch(pl_frq + d, f)
            p_tmp = p[f_pl]
            f_tmp = f[f_pl]
            peaks, _ = findpeaks1d(p_tmp, prominence=pr)
            @assert length(peaks) > 0 "No power line noise peak detected, check pl_frq value (perhaps the signal has already been filtered?)."
            pl_amp = vsearch(maximum(p_tmp[peaks]), p_tmp)
            pl_frq_detected = f_tmp[vsearch(pl_amp, p_tmp)]
            @assert !(pl_frq_detected < pl_frq - d || pl_frq_detected > pl_frq + d) "Channel $(clabels[ch_idx]): power line noise peak detected at $pl_frq_detected, check pl_frq value."
            # pl_wdth = peakwidths1d(p_tmp, peaks)[1][vsearch(pl_amp, peaks)]

            v = zeros(length(bw_values))
            Threads.@threads for bw_idx in eachindex(bw_values)
                obj_tmp = NeuroAnalyzer.filter(obj_new, ch=clabels[ch_idx], fprototype=:iirnotch, cutoff=pl_frq, bw=bw_values[bw_idx])
                p, f = psd(obj_tmp, ch=clabels[ch_idx], db=true)
                f1 = vsearch(pl_frq - d, f)
                f2 = vsearch(pl_frq + d, f)
                seg = p[f1:f2]
                v[bw_idx] = var(seg)
            end

            push!(pl_best_bw, bw_values[vsearch(minimum(v), v)])
            NeuroAnalyzer.filter!(obj_new, ch=clabels[ch_idx], fprototype=:iirnotch, cutoff=pl_frq, bw=pl_best_bw[ch_idx])

            # detect peaks
            p, f = psd(obj_new, ch=clabels[ch_idx], db=true)
            p = p[:]
            f_pl = vsearch(2 * pl_frq - 2 * d, f)
            p_tmp = p[f_pl:end]
            f_tmp = f[f_pl:end]

            # prominence of 2 dB, distance of 5 Hz
            if length(pks_frq) == 0
                peaks, _ = findpeaks1d(p_tmp, prominence=pr, distance=vsearch(d, f))
                if length(peaks) > 0
                    for idx in eachindex(peaks)
                        push!(pks_frq, f_tmp[peaks[idx]])
                    end
                    pks_wdth = peakwidths1d(p_tmp, peaks)[1]
                end
            end

            n = length(pks_frq)
            bw_values = collect(0:q:5.0)[2:end]
            if n > 0
                best_bw = zeros(n)
                @inbounds for peak_idx in eachindex(pks_frq)
                    v = zeros(length(bw_values))
                    Threads.@threads for bw_idx in eachindex(bw_values)
                        obj_tmp = NeuroAnalyzer.filter(obj_new, ch=clabels[ch_idx], fprototype=:iirnotch, cutoff=pks_frq[peak_idx], bw=bw_values[bw_idx])
                        p, f = psd(obj_tmp, ch=clabels[ch_idx], db=true)
                        f1 = vsearch(pks_frq[peak_idx] - d / 2, f)
                        f2 = vsearch(pks_frq[peak_idx] + d / 2, f)
                        seg = p[f1:f2]
                        seg_f = f[f1:f2]
                        v[bw_idx] = var(seg)
                    end
                    best_bw[peak_idx] = bw_values[vsearch(minimum(v), v)]
                end
                @inbounds for peak_idx in eachindex(pks_frq)
                    NeuroAnalyzer.filter!(obj_new, ch=clabels[ch_idx], fprototype=:iirnotch, cutoff=pks_frq[peak_idx], bw=best_bw[peak_idx])
                end
                push!(pks_best_bw, best_bw)
            end

            # update progress bar
            progress_bar && next!(progbar)
        end

        df = DataFrame("channel"=>clabels[ch], "power line bandwidth"=>pl_best_bw)
        if length(pks_frq) > 0
            pks_best_bw_m = zeros(length(ch), length(pks_frq))
            pks_frq_m = zeros(length(ch), length(pks_frq))
            peak_names = repeat(["peak "], length(pks_frq))
            peak_frqs = repeat(["peak "], length(pks_frq))
            for idx in eachindex(pks_frq)
                peak_names[idx] = peak_names[idx] * "$(idx) bandwidth"
                peak_frqs[idx] = peak_frqs[idx] * "$(idx) frequency"
            end
            for idx in eachindex(ch)
                pks_best_bw_m[idx, :] = pks_best_bw[idx]
                pks_frq_m[idx, :] = pks_frq
            end
            df1 = DataFrame(pks_frq_m, peak_frqs)
            df2 = DataFrame(pks_best_bw_m, peak_names)
            df = hcat(df, df1)
            df = hcat(df, df2)
        end

        NeuroAnalyzer.verbose = verbose_tmp

        reset_components!(obj_new)
        push!(obj_new.history, "remove_powerline(OBJ, pl_frq=$pl_frq, method=:$method, pr=$pr, d=$d, q=$q)")

        return obj_new, df

    end

end

"""
    remove_powerline!(obj; <keyword arguments>)

Remove power line noise and harmonics.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `pl_frq::Real=obj.header.recording[:line_frequency]`: power line frequency, default is read from the OBJ header
- `method::Symbol=:iir`: use IIR filter
- `pr::Real=2.0`: prominence of noise peaks in dB
- `d::Real=5.0`: minimum distance between peaks in Hz
- `q::Real=0.1`: optimization step size

# Returns

- `df::DataFrame`: list of peaks detected
"""
function remove_powerline!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pl_frq::Real=obj.header.recording[:line_frequency], method::Symbol=:iir, pr::Real=2.0, d::Float64=5.0, q::Real=0.1)::DataFrame

    obj_new, df = remove_powerline(obj, ch=ch, pl_frq=pl_frq, method=method, pr=pr, d=d, q=q)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return df

end

"""
    detect_powerline(s; <keyword arguments>)

Detect power line noise frequency.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate

# Returns

- `noise_frq::Float64`: peak noise frequency in Hz
"""
function detect_powerline(s::AbstractVector; fs::Int64)::Float64

    n = length(s)
    t = linspace(0, n * 1/fs, n)
    noise_power = zeros(fs ÷ 2)

    for noise_idx in 1:(fs ÷ 2)
        df = DataFrame(:signal=>s, :b0=>ones(n), :b1=>sin.(2 * pi * noise_idx .* t), :b2=>cos.(2 * pi * noise_idx .* t))
        lr = GLM.lm(@formula(signal ~ b0 + b1 + b2), df)
        b = StatsKit.coef(lr)
        noise_power[noise_idx] = b[3]^2 + b[4]^2
    end

    # maximum noise frequency in Hz
    noise_frq = vsearch(maximum(noise_power), noise_power)

    return noise_frq

end

"""
    detect_powerline(obj)

Detect power line noise frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `noise_frq::Array{Float64, 3}`: peak noise frequency in Hz
"""
function detect_powerline(obj::NeuroAnalyzer.NEURO)::Array{Float64, 3}

    ch_n = size(obj, 1)
    ep_n = size(obj, 3)

    noise_frq = zeros(ch_n, ep_n)

    # initialize progress bar
    progress_bar && (progbar = Progress(ch_n * ep_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        @Threads.threads for ch_idx in 1:ch_n
            noise_frq[ch_idx, ep_idx] = @views detect_powerline(obj.data[ch_n, :, ep_n], fs=sr(obj))
        end

        # update progress bar
        progress_bar && next!(progbar)
    end

    return noise_frq

end

"""
    detect_powerline!(obj)

Detect power line noise frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Nothing
"""
function detect_powerline!(obj::NeuroAnalyzer.NEURO)::Nothing

    ch_n = size(obj, 1)
    ep_n = size(obj, 3)

    noise_frq = zeros(ch_n, ep_n)

    # initialize progress bar
    progress_bar && (progbar = Progress(ch_n * ep_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        @Threads.threads for ch_idx in 1:ch_n
            noise_frq[ch_idx, ep_idx] = @views detect_powerline(obj.data[ch_n, :, ep_n], fs=sr(obj))
        end

        # update progress bar
        progress_bar && next!(progbar)
    end

    _info("Power line noise detected at $(median(noise_frq)) Hz")
    obj.header.recording[:line_frequency] = median(noise_frq)

    return nothing

end
