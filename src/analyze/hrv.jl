export hrv_detect
export hrv_analyze

"""
    hrv_detect(obj)

Detect heart rate variability (HRV) R-peaks from the ECG channel.

The ECG channel is located automatically from the `channel_type` field of `obj.header.recording`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

Named tuple:

- `nn_seg::Vector{Float64}`: NN (normal-to-normal) interval durations [ms]
- `r_idx::Vector{Float64}`: sample indices of detected R-peaks
"""
function hrv_detect(obj::NeuroAnalyzer.NEURO)::@NamedTuple{nn_seg::Vector{Float64}, r_idx::Vector{Float64}}

    @assert "ecg" in obj.header.recording[:channel_type] "OBJ does not contain ECG channel."

    # locate the ECG channel index
    ch = get_channel(obj, type = "ecg")
    _info("ECG channel found: $(ch[1])")
    ch = get_channel(obj, ch = ch)

    # flatten the ECG channel across all epochs into a single continuous vector
    # this treats a multi-epoch recording as one unbroken signal for peak detection
    ecg = vec(obj.data[ch, :, :])

    # detect R-peaks: threshold at mean + 2σ to select prominent peaks only
    r_idx, _ = findpeaks1d(ecg, height = mean(ecg) + 2 * std(ecg))

    # convert inter-peak sample distances to milliseconds
    nn_seg = diff(r_idx) ./ sr(obj) * 1000

    _info("Detected NN segments: $(length(nn_seg))")

    return (nn_seg = nn_seg, r_idx = r_idx)

end

"""
    hrv_analyze(nn_seg)

Analyze time-domain heart rate variability (HRV) statistics from NN intervals. NN interval = interval between successive R-peaks (equivalent to RR interval when ectopic beats have been removed, hence "Normal-to-Normal").

# Arguments

- `nn_seg::Vector{Float64}`: list of NN segments [ms]

# Returns

Named tuple:

- `menn::Float64`: mean NN interval [ms]
- `mdnn::Float64`: median NN interval [ms]
- `vnn::Float64`: variance of NN intervals [ms²]
- `sdnn::Float64`: standard deviation of NN intervals [ms]
- `rmssd::Float64`: root mean square of successive differences between adjacent NN intervals [ms] — √(mean(diff(NN)²))
- `sdsd::Float64`: standard deviation of successive differences between adjacent NN intervals [ms] — std(diff(NN))
- `nn50::Int64`: number of pairs of successive NNs differing by > 50 ms
- `pnn50::Float64`: proportion of NN50 to total number of NN intervals
- `nn20::Int64`: number of pairs of successive NNs differing by > 20 ms
- `pnn20::Float64`: proportion of NN20 to total number of NN intervals
"""
function hrv_analyze(
        nn_seg::Vector{Float64}
    )::@NamedTuple{
        menn::Float64,
        mdnn::Float64,
        vnn::Float64,
        sdnn::Float64,
        rmssd::Float64,
        sdsd::Float64,
        nn50::Int64,
        pnn50::Float64,
        nn20::Int64,
        pnn20::Float64,
    }

    nn_diff = diff(nn_seg)
    nn_total = length(nn_seg)

    menn = round(mean(nn_seg), digits = 3)
    mdnn = round(median(nn_seg), digits = 3)
    vnn = round(var(nn_seg), digits = 3)
    sdnn = round(std(nn_seg), digits = 3)

    # RMSSD: root mean square of successive differences.
    rmssd = round(sqrt(sum(abs2, nn_diff) / length(nn_diff)), digits = 3)
    sdsd  = round(std(nn_diff), digits = 3)

    # NN50 / NN20: count pairs of successive intervals differing by > threshold.
    nn50  = count(d -> abs(d) > 50, nn_diff)
    pnn50 = round(nn50 / nn_total, digits = 3)

    nn20  = count(d -> abs(d) > 20, nn_diff)
    pnn20 = round(nn20 / nn_total, digits = 3)

    return (
        menn = menn,
        mdnn = mdnn,
        vnn = vnn,
        sdnn = sdnn,
        rmssd = rmssd,
        sdsd = sdsd,
        nn50 = nn50,
        pnn50 = pnn50,
        nn20 = nn20,
        pnn20 = pnn20,
    )

end
