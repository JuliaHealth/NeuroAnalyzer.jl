export remove_powerline
export remove_powerline!

"""
    remove_powerline(obj; <keyword arguments>)

Remove power line noise and its peaks above power line frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pl_frq::Real=50`: power line frequency
- `method::Symbol=:iir`: use IIR filter
- `pr::Real=2.0`: prominence of noise peaks in dB
- `d::Real=5.0`: minimum distance between peaks in Hz
- `q::Real=0.1`: optimization step size

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function remove_powerline(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pl_frq::Real=50, method::Symbol=:iir, pr::Real=2.0, d::Float64=5.0, q::Real=0.1)

    @assert epoch_n(obj) == 1 "remove_powerline() should be applied to a continuous signal."

    _check_channels(obj, ch)    
    _check_var(method, [:iir], "method")

    obj_new = deepcopy(obj)

    @assert pl_frq >= 0 "pl_freq must be ≥ 0."
    @assert pl_frq <= sr(obj) / 2 "pl_freq must be ≤ $(sr(obj) / 2)."
    @assert q >= 0.01 "q must be ≥ 0.01."
    @assert q < 5 "q must be < 5."

    # number of peaks
    pl_best_bw = zeros(length(ch))
    pl_frq_detected = zeros(length(ch))
    pks_frq = Vector{Float64}()
    pks_best_bw = Vector{Float64}[]

    if method === :iir

        _info("Removing power line noise at $pl_frq Hz and its peaks")

        bw_values = collect(0:q:20.0)[2:end]

        # initialize progress bar
        progress_bar == true && (progbar = Progress(length(ch), dt=1, barlen=20, color=:white))

        # Threads.@threads for ch_idx in eachindex(ch)
        @inbounds @simd for ch_idx in eachindex(ch)
            # detect power line peak
            p, f = psd(obj_new, ch=ch[ch_idx], norm=true)
            p = p[:]
            p[end] = p[end - 1]
            f_pl = vsearch(pl_frq - 5 * d, f):vsearch(pl_frq + 5 * d, f)
            p_tmp = p[f_pl]
            f_tmp = f[f_pl]
            peaks, _ = findpeaks1d(p_tmp, prominence=pr)
            @assert length(peaks) > 0 "No power line peak detected, check pl_frq value (perhaps the signal has already been filtered?)."
            pl_amp = vsearch(maximum(p_tmp[peaks]), p_tmp)
            pl_frq_detected = f_tmp[vsearch(pl_amp, p_tmp)]
            @assert !(pl_frq_detected < pl_frq - d || pl_frq_detected > pl_frq + d) "Power line peak detected at $pl_frq_detected, check pl_frq value."
            # pl_wdth = peakwidths1d(p_tmp, peaks)[1][vsearch(pl_amp, peaks)]

            v = zeros(length(bw_values))
            Threads.@threads for bw_idx in eachindex(bw_values)
                obj_tmp = NeuroAnalyzer.filter(obj_new, ch=ch[ch_idx], fprototype=:iirnotch, cutoff=pl_frq, bw=bw_values[bw_idx])
                p, f = psd(obj_tmp, ch=ch[ch_idx], norm=true)
                f1 = vsearch(pl_frq - 5 * d, f)
                f2 = vsearch(pl_frq + 5 * d, f)
                seg = p[f1:f2]
                v[bw_idx] = var(seg)
            end

            pl_best_bw[ch_idx] = bw_values[vsearch(minimum(v), v)]
            NeuroAnalyzer.filter!(obj_new, ch=ch[ch_idx], fprototype=:iirnotch, cutoff=pl_frq, bw=pl_best_bw[ch_idx])

            # detect peaks
            p, f = psd(obj_new, ch=ch[ch_idx], norm=true)
            p = p[:]
            p[end] = p[end - 1]
            f_pl = vsearch(2 * pl_frq - 2 * d, f)
            p_tmp = p[f_pl:end]
            f_tmp = f[f_pl:end]

            # prominence of 2 dB, distance of 5 Hz
            if length(pks_frq) == 0
                peaks, properties = findpeaks1d(p_tmp, prominence=pr, distance=vsearch(d, f))
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
                @inbounds @simd for peak_idx in 1:n
                    v = zeros(length(bw_values))
                    Threads.@threads for bw_idx in eachindex(bw_values)
                        obj_tmp = NeuroAnalyzer.filter(obj_new, ch=ch[ch_idx], fprototype=:iirnotch, cutoff=pks_frq[peak_idx], bw=bw_values[bw_idx])
                        p, f = psd(obj_tmp, ch=ch[ch_idx], norm=true)
                        f1 = vsearch(pks_frq[peak_idx] - d / 2, f)
                        f2 = vsearch(pks_frq[peak_idx] + d / 2, f)
                        seg = p[f1:f2]
                        seg_f = f[f1:f2]
                        v[bw_idx] = var(seg)
                    end
                    best_bw[peak_idx] = bw_values[vsearch(minimum(v), v)]
                end
                @inbounds @simd for peak_idx in 1:n
                    NeuroAnalyzer.filter!(obj_new, ch=ch[ch_idx], fprototype=:iirnotch, cutoff=pks_frq[peak_idx], bw=best_bw[peak_idx])
                end
                push!(pks_best_bw, best_bw)
            end

            # update progress bar
            progress_bar == true && next!(progbar)
        end

        df = DataFrame("channel"=>ch, "power line bandwidth"=>pl_best_bw)
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

        verbose && println(df)

        reset_components!(obj_new)
        push!(obj_new.history, "remove_powerline(OBJ, pl_frq=$pl_frq, method=:$method, pr=$pr, d=$d, q=$q)")

        return obj_new

    end

end

"""
    remove_powerline!(obj; <keyword arguments>)

Remove power line noise and harmonics.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pl_frq::Real=50`: power line frequency
- `method::Symbol=:iir`: use IIR filter
- `pr::Real=2.0`: prominence of noise peaks in dB
- `d::Real=5.0`: minimum distance between peaks in Hz
- `q::Real=0.1`: optimization step size
"""
function remove_powerline!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pl_frq::Real=50, method::Symbol=:iir, pr::Real=2.0, d::Float64=5.0, q::Real=0.1)

    obj_new = remove_powerline(obj, ch=ch, pl_frq=pl_frq, method=method, pr=pr, d=d, q=q)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end