export lat_idx

"""
    lat_idx(obj; <keyword arguments>)

Calculate Lateralization Index. LI is the log-ratio of right vs left hemisphere power at a given frequency (or frequency range):

LI = log( P_right / P_left )

LI > 0 → right-dominant
LI < 0 → left-dominant
LI = 0 → symmetric

Hemisphere channels are detected automatically via `channel_pick()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `frq::Union{Real, Tuple{<:Real, <:Real}}`: frequency at which the index is calculated; if range is provided, than averaged index across the range is calculated
- `method::Symbol=:welch`: PSD method:
- `:welch`: Welch's periodogram
- `:fft`: fast Fourier transform
- `:mt`: multi-tapered periodogram
- `:stft`: short-time Fourier transform
- `:mw`: Morlet wavelet convolution
- `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `lidx::Float64`: lateralization index (positive = right-dominant)
"""
function lat_idx(
    obj::NeuroAnalyzer.NEURO;
    frq::Union{Real, Tuple{<:Real, <:Real}},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::Float64

    _check_datatype(obj, ["meg", "eeg", "erp", "erf"])

    @assert length(channel_pick(obj, p = :l)) > 0 "Could not detect left hemisphere channels, check OBJ labels."
    @assert length(channel_pick(obj, p = :r)) > 0 "Could not detect right hemisphere channels, check OBJ labels."

    ch_l = get_channel(obj, ch = channel_pick(obj, p = :l))
    ch_r = get_channel(obj, ch = channel_pick(obj, p = :r))

    # shared PSD keyword arguments — avoids repeating them four times
    psd_kwargs = (
        fs = sr(obj),
        db = false,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )

    _log_off()

    # compute PSDs for left and right hemispheres
    # for ERP/ERF use epoch 1 (the trial-averaged waveform); keep it 3-D with
    # a singleton epoch dimension so _chk3d() inside psd() does not fail
    if datatype(obj) in ["erp", "erf"]
        psd_l = psd(@view(obj.data[ch_l, :, 1:1]); psd_kwargs...)
        psd_r = psd(@view(obj.data[ch_r, :, 1:1]); psd_kwargs...)
    else
        psd_l = psd(@view(obj.data[ch_l, :, :]); psd_kwargs...)
        psd_r = psd(@view(obj.data[ch_r, :, :]); psd_kwargs...)
    end

    _log_on()

    f = psd_l.f
    p_left = psd_l.p  # (ch_l, freq, ep_n)
    p_right = psd_r.p # (ch_r, freq, ep_n)

    # average across epochs (no-op for ERP/ERF which already has ep_n=1),
    # then average across channels to get a single (freq,) power vector
    if size(p_left, 3) > 1
        p_left  = dropdims(mean(p_left,  dims = 3), dims = 3)
        p_right = dropdims(mean(p_right, dims = 3), dims = 3)
    else
        p_left  = dropdims(p_left,  dims = 3)
        p_right = dropdims(p_right, dims = 3)
    end
    # now (ch, freq) — average over channels to get (freq,)
    p_left  = vec(mean(p_left,  dims = 1))
    p_right = vec(mean(p_right, dims = 1))

    # average across epochs
    size(p_left, 3) > 1 && (p_left = mean(p_left, dims = 3))
    size(p_right, 3) > 1 && (p_right = mean(p_right, dims = 3))
    # average across channels
    p_left = mean(p_left, dims = 1)
    p_right = mean(p_right, dims = 1)

    # extract power at the requested frequency or frequency range
    if frq isa Real
        frq_idx = vsearch(frq, f)
        p_left_val = p_left[frq_idx]
        p_right_val = p_right[frq_idx]
    else
        _check_tuple(frq, (0, sr(obj) / 2), "frq")
        frq_idx1 = vsearch(frq[1], f)
        frq_idx2 = vsearch(frq[2], f)
        p_left_val = mean(@view(p_left[frq_idx1:frq_idx2]))
        p_right_val = mean(@view(p_right[frq_idx1:frq_idx2]))
    end

    lidx = log(p_right / p_left)

    return lidx

end
