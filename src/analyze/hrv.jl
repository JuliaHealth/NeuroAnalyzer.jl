export hrv_detect
export hrv_analyze

"""
    hrv_detect(obj)

Detect heart rate variability (HRV). Requires ECG channel (which will be automatically detected based on `channel_type` field).

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple containing:
- `nn_seg::Vector{Float64}`: list of NN segments [msec]
- `r_idx::Vector{Float64}`: index of R peaks
"""
function hrv_detect(obj::NeuroAnalyzer.NEURO)::NamedTuple{(:nn_seg, :r_idx), Tuple{Vector{Float64}, Vector{Float64}}}

    @assert "ecg" in obj.header.recording[:channel_type] "OBJ does not contain ECG channel."
    ch = get_channel(obj, type="ecg")
    _info("ECG channel found: $(labels(obj)[ch])")
    ecg = eeg.data[ch, :, :][:]
    r_idx, _ = findpeaks1d(ecg, height=mean(ecg) + 2*std(ecg))

    # convert to ms
    nn_seg = diff(r_idx) ./ sr(eeg) * 1000

    _info("Detected NN segments: $(length(nn_seg))")

    return (nn_seg=nn_seg, r_idx=r_idx)

end

"""
    hrv_analyze(nn_seg)

Analyze heart rate variability (HRV).

# Arguments

- `nn_seg::Vector{Float64}`: list of NN segments [msec]

# Returns

Named tuple containing:
- `menn::Float64`: the mean of NN segments
- `mdnn::Float64`: the median of NN segments
- `vnn::Float64`: the variance of NN segments
- `sdnn::Float64`: the standard deviation of NN segments
- `rmssd::Float64`: ("root mean square of successive differences"), the square root of the mean of the squares of the successive differences between adjacent NNs
- `sdsd::Float64`: ("standard deviation of successive differences"), the standard deviation of the successive differences between adjacent NNs
- `nn50::Float64`: the number of pairs of successive NNs that differ by more than 50 ms
- `pnn50::Float64`, the proportion of NN50 divided by total number of NNs
- `nn20::Float64`, the number of pairs of successive NNs that differ by more than 20 ms
- `pnn20::Float64`, the proportion of NN20 divided by total number of NNs
"""
function hrv_analyze(nn_seg::Vector{Float64})::NamedTuple{(:menn, :mdnn, :vnn, :sdnn, :rmssd, :sdsd, :nn50, :pnn50, :nn20, :pnn20), Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}

    nn_diff = diff(nn_seg)

    menn = round(mean(nn_seg), digits=3)
    mdnn =  round(median(nn_seg), digits=3)
    vnn =  round(var(nn_seg), digits=3)
    sdnn =  round(std(nn_seg), digits=3)
    rmssd = round(sqrt(mean(nn_diff .^ 2)), digits=3)
    sdsd = round(std(nn_diff), digits=3)
    nn50 = length(findall(abs.(nn_diff) .> 50))
    pnn50 = round(nn50 / length(nn_seg), digits=3)
    nn20 = length(findall(abs.(nn_diff) .> 20))
    pnn20 = round(nn20 / length(nn_seg), digits=3)

    return(menn=menn, mdnn=mdnn, vnn=vnn, sdnn=sdnn, rmssd=rmssd, sdsd=sdsd, nn50=nn50, pnn50=pnn50, nn20=nn20, pnn20=pnn20)

end