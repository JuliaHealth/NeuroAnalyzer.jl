# Example NeuroAnalyzer pipelines for typical EEG analyses.
#
# Usage (from the NeuroAnalyzer.jl folder):
#
#   using NeuroAnalyzer
#   include("pipelines/test_pipeline.jl")
#
#   eeg = import_recording("recording.edf")
#   eeg = preprocess_pipeline(eeg; locs_file = "standard-10-20.ced")
#   save_neuro(eeg; file_name = "recording_clean.hdf")
#
#   # resting state: band powers, relative powers, alpha peak, aperiodic exponent (per channel)
#   r = spectral_pipeline(eeg)
#
#   # task: ERP time-locked to marker "S  1" (200 ms before to 800 ms after the stimulus)
#   r = erp_pipeline(eeg; marker = "S  1", offset = 0.2, ep_len = 1.0)
#   plot_erp(r.erp; ch = "Cz")

using DataFrames

"""
    test_pipeline(eeg; <keyword arguments>)

Basic EEG processing pipeline: power line noise removal, band-pass filtering and common average reference.

# Arguments

- `eeg::NeuroAnalyzer.NEURO`: continuous EEG recording
- `pl_frq::Real=50`: power line frequency in Hz
- `hp::Real=0.5`: high-pass filter cutoff in Hz
- `lp::Real=45`: low-pass filter cutoff in Hz

# Returns

- `NeuroAnalyzer.NEURO`: processed copy of `eeg`
"""
function test_pipeline(
    eeg::NeuroAnalyzer.NEURO;
    pl_frq::Real = 50,
    hp::Real = 0.5,
    lp::Real = 45,
)::NeuroAnalyzer.NEURO
    obj = deepcopy(eeg)
    remove_powerline!(obj; ch = "eeg", pl_frq = pl_frq)
    _bandpass!(obj; hp = hp, lp = lp)
    reference_avg!(obj)
    return obj
end

# zero-phase band-pass: a 4th-order Butterworth high-pass (a FIR filter for low cutoffs would need
# thousands of taps) and a Hamming-window FIR low-pass with a 2 Hz transition band
function _bandpass!(obj::NeuroAnalyzer.NEURO; hp::Real, lp::Real, ch = "eeg")::Nothing
    NeuroAnalyzer.filter!(obj; ch = ch, fprototype = :butterworth, ftype = :hp, cutoff = hp, order = 4, report = false)
    NeuroAnalyzer.filter!(obj; ch = ch, fprototype = :fir, ftype = :lp, cutoff = lp, bw = 2, report = false)
    return nothing
end

"""
    preprocess_pipeline(eeg; <keyword arguments>)

Typical preprocessing of a continuous EEG recording:

1. channel locations (optional)
2. DC removal, power line noise removal (notch filters at the line frequency and its harmonics)
3. zero-phase band-pass filtering
4. bad channel detection (flat, extreme amplitude, outlying variance) and interpolation from the neighboring channels (requires locations)
5. common average reference
6. optional ICA-based artifact removal (e.g. eye blinks): components to remove must be chosen after inspection, e.g. with `plot_icatopo()`

# Arguments

- `eeg::NeuroAnalyzer.NEURO`: continuous EEG recording
- `locs_file::Union{Nothing, String}=nothing`: channel locations file (e.g. `.ced`, `.locs`, `.elc`)
- `pl_frq::Real=eeg.header.recording[:line_frequency]`: power line frequency in Hz
- `hp::Real=0.5`: high-pass filter cutoff in Hz
- `lp::Real=45`: low-pass filter cutoff in Hz
- `bad_channels::Bool=true`: detect and interpolate bad channels
- `amp_t::Real=500`: amplitude threshold (μV) for bad channels
- `ica_n::Union{Nothing, Int64}=nothing`: number of ICA components (`nothing`: no ICA)
- `ica_remove::Vector{Int64}=Int64[]`: ICA components to remove

# Returns

- `NeuroAnalyzer.NEURO`: preprocessed copy of `eeg`; detected bad channels are marked in `obj.header.recording[:bad_channel]` (and interpolated if locations are available)
"""
function preprocess_pipeline(
    eeg::NeuroAnalyzer.NEURO;
    locs_file::Union{Nothing, String} = nothing,
    pl_frq::Real = eeg.header.recording[:line_frequency],
    hp::Real = 0.5,
    lp::Real = 45,
    bad_channels::Bool = true,
    amp_t::Real = 500,
    ica_n::Union{Nothing, Int64} = nothing,
    ica_remove::Vector{Int64} = Int64[],
)::NeuroAnalyzer.NEURO
    nepochs(eeg) == 1 || throw(ArgumentError("preprocess_pipeline() expects a continuous (non-epoched) recording."))
    obj = deepcopy(eeg)

    # 1. channel locations
    !isnothing(locs_file) && load_locs!(obj; file_name = locs_file)

    # 2. DC offset and power line noise
    remove_dc!(obj; ch = "eeg")
    remove_powerline!(obj; ch = "eeg", pl_frq = pl_frq)

    # 3. band-pass filtering (before bad channel detection: slow drifts would dominate the statistics)
    _bandpass!(obj; hp = hp, lp = lp)

    # 4. bad channels: flat signal, extreme amplitude, outlying variance
    if bad_channels
        channel_reject!(obj; ch = "eeg", method = [:flat, :amp, :var], amp_t = amp_t)
        bad = labels(obj)[obj.header.recording[:bad_channel]]
        if !isempty(bad) && !isempty(obj.locs)
            for ch in bad
                # spherical-spline-like interpolation from the neighboring channels
                plinterpolate_channel!(obj; ch = ch, ep = 1)
            end
            obj.header.recording[:bad_channel] .= false
        end
    end

    # 5. common average reference (bad channels are excluded from the average if not interpolated)
    reference_avg!(obj)

    # 6. ICA
    if !isnothing(ica_n)
        ica = ica_decompose(obj; ch = "eeg", n = ica_n)
        isempty(ica_remove) ||
            ica_remove!(obj; ch = "eeg", ic = ica.ic, ic_mw = ica.ic_mw, ic_idx = ica_remove)
    end

    return obj
end

"""
    spectral_pipeline(eeg; <keyword arguments>)

Resting-state spectral analysis of a (preprocessed) continuous EEG recording. The signal is split into epochs, epochs with artifacts are rejected and, for every channel, the band powers (averaged over epochs), relative band powers, the individual alpha peak frequency and the aperiodic (1/f) exponent (FOOOF / specparam) are calculated.

# Arguments

- `eeg::NeuroAnalyzer.NEURO`: continuous EEG recording
- `ep_len::Real=4`: epoch length in seconds
- `amp_t::Real=150`: amplitude threshold (μV) for epoch rejection
- `bands::Vector{Symbol}=[:delta, :theta, :alpha, :beta, :gamma_1]`: frequency bands (see `band_frq()`)
- `flim::Tuple{Real, Real}=(1, 40)`: frequency range of the total power and of the aperiodic fit

# Returns

Named tuple:

- `df::DataFrame`: one row per channel: absolute (`<band>`) and relative (`<band>_rel`) band powers, alpha peak frequency (`apf`), aperiodic offset and exponent
- `rejected::Vector{Int64}`: rejected epochs
"""
function spectral_pipeline(
    eeg::NeuroAnalyzer.NEURO;
    ep_len::Real = 4,
    amp_t::Real = 150,
    bands::Vector{Symbol} = [:delta, :theta, :alpha, :beta, :gamma_1],
    flim::Tuple{Real, Real} = (1, 40),
)::@NamedTuple{df::DataFrame, rejected::Vector{Int64}}
    obj = epoch(eeg; ep_len = ep_len)

    # epochs with extreme amplitudes in any EEG channel
    rejected = epoch_reject(obj; ch = "eeg", method = :amp, amp_t = amp_t)
    length(rejected) < nepochs(obj) || throw(ArgumentError("All epochs were rejected; check amp_t."))
    isempty(rejected) || delete_epoch!(obj; ep = rejected)

    ch = get_channel(obj; type = "eeg")
    mean_ep(m) = vec(sum(m; dims = 2) ./ size(m, 2))
    df = DataFrame(:channel => ch)

    # absolute and relative band powers (Welch PSD, 2 s windows)
    wlen = min(2 * sr(obj), epoch_len(obj))
    total = mean_ep(band_power(obj; ch = ch, flim = flim, wlen = wlen))
    for b in bands
        bp = mean_ep(band_power(obj; ch = ch, flim = band_frq(obj; band = b), wlen = wlen))
        df[!, b] = bp
        df[!, Symbol(string(b) * "_rel")] = bp ./ total
    end

    # individual alpha peak frequency
    df[!, :apf] = mean_ep(peak_frq(obj; ch = ch, flim = (7, 14), wlen = wlen))

    # aperiodic component: offset and exponent of the 1/f fit (FOOOF, fixed mode)
    ff = fooof(obj; ch = ch, flim = flim, wlen = wlen)
    df[!, :ap_offset] = vec(sum(ff.ap[:, 1, :]; dims = 2) ./ size(ff.ap, 3))
    df[!, :ap_exponent] = vec(sum(ff.ap[:, 3, :]; dims = 2) ./ size(ff.ap, 3))

    return (; df, rejected)
end

"""
    erp_pipeline(eeg; <keyword arguments>)

Event-related potential (ERP) analysis of a (preprocessed) continuous EEG recording: epochs time-locked to a marker, artifact rejection, baseline correction and averaging, and the latencies and amplitudes of the largest positive and negative deflections.

# Arguments

- `eeg::NeuroAnalyzer.NEURO`: continuous EEG recording with markers
- `marker::String`: marker value (see `eeg.markers`)
- `offset::Real=0.2`: epoch start before the marker in seconds (also the baseline)
- `ep_len::Real=1.0`: epoch length in seconds
- `amp_t::Real=100`: amplitude threshold (μV) for epoch rejection
- `twin::Tuple{Real, Real}=(0, ep_len - offset)`: time window (s, from the marker) for the peaks

# Returns

Named tuple:

- `erp::NeuroAnalyzer.NEURO`: ERP object (epoch 1: the average, epochs 2…: the accepted epochs)
- `rejected::Vector{Int64}`: rejected epochs
- `peaks::DataFrame`: one row per channel: latency (s) and amplitude of the positive and the negative peak within `twin`
"""
function erp_pipeline(
    eeg::NeuroAnalyzer.NEURO;
    marker::String,
    offset::Real = 0.2,
    ep_len::Real = 1.0,
    amp_t::Real = 100,
    twin::Tuple{Real, Real} = (0, ep_len - offset),
)::@NamedTuple{erp::NeuroAnalyzer.NEURO, rejected::Vector{Int64}, peaks::DataFrame}
    obj = epoch(eeg; marker = marker, offset = offset, ep_len = ep_len)

    rejected = epoch_reject(obj; ch = "eeg", method = :amp, amp_t = amp_t)
    length(rejected) < nepochs(obj) || throw(ArgumentError("All epochs were rejected; check amp_t."))
    isempty(rejected) || delete_epoch!(obj; ep = rejected)

    # baseline: the pre-stimulus interval
    erp = average_epochs(obj; bl = (-offset, 0))

    # peaks of the average within the time window
    ch = get_channel(erp; ch = get_channel(erp; type = datatype(erp) == "erp" ? "eeg" : datatype(erp)))
    t = erp.epoch_time
    w = findall(x -> twin[1] <= x <= twin[2], t)
    peaks = DataFrame(:channel => labels(erp)[ch])
    s = erp.data[ch, w, 1]
    ipos = [argmax(s[i, :]) for i in axes(s, 1)]
    ineg = [argmin(s[i, :]) for i in axes(s, 1)]
    peaks[!, :pos_latency] = t[w][ipos]
    peaks[!, :pos_amplitude] = [s[i, ipos[i]] for i in axes(s, 1)]
    peaks[!, :neg_latency] = t[w][ineg]
    peaks[!, :neg_amplitude] = [s[i, ineg[i]] for i in axes(s, 1)]

    return (; erp, rejected, peaks)
end
