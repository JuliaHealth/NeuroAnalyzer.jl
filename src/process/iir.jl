# notch filters (method = :iir of remove_powerline()): zero-phase IIR notch filters at the power
# line frequency and its harmonics, bandwidth optimized for every channel and peak

# notch filtering of channels `ch` (indices) of a continuous object, in-place; returns one row per
# filtered peak: `:channel` (label), `:frequency` (Hz), `:bandwidth` (optimized notch bandwidth,
# Hz); the filtered frequencies are added to the `:notch` header field
function _iir!(
    obj::NeuroAnalyzer.NEURO;
    ch::Vector{Int64},
    fline::Real,
    pr::Real,
    d::Real,
    q::Real,
)::DataFrame
    # validate
    fs = sr(obj)
    q >= 0.01 || throw(ArgumentError("q must be ≥ 0.01."))
    q < 5 || throw(ArgumentError("q must be < 5."))
    d > 0 || throw(ArgumentError("d must be > 0."))
    clabels = labels(obj)

    # the fundamental and its harmonics below the Nyquist frequency
    frqs = [k * fline for k = 1:floor(Int64, fs / 2 / fline) if k * fline < fs / 2]
    bw_fund = collect(q:q:10.0)
    bw_harm = collect(q:q:5.0)

    df = DataFrame(:channel => String[], :frequency => Float64[], :bandwidth => Float64[])
    no_peak = String[]

    _info("Removing power line noise at $fline Hz and its harmonics")

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
        "No power line peak at $fline Hz in channel(s) $(join(no_peak, ", ")); power line frequency not filtered (check pl_frq or whether the signal is already filtered).",
    )

    return df
end
