export tpt_detect
export tpt_analyze

"""
    tpt_detect(obj)

Detect pinches in TPT recording.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `p_idx::Vector{Int64}`: index of pinches locations
"""
function tpt_detect(obj::NeuroAnalyzer.NEURO)::Vector{Int64}

    _check_datatype(obj, "tpt")

    pos_x = obj.data[1, :, :][:]
    pos_y = obj.data[2, :, :][:]
    pos_z = obj.data[3, :, :][:]
    acc_x = obj.data[4, :, :][:]
    acc_y = obj.data[5, :, :][:]
    acc_z = obj.data[6, :, :][:]

    p_idx_x, _ = findpeaks1d(pos_x, distance=(sr(obj) ÷ 2), height=mean(pos_x) + 2*std(pos_x))
    p_idx_y, _ = findpeaks1d(pos_y, distance=(sr(obj) ÷ 2), height=mean(pos_y) + 2*std(pos_y))
    p_idx_z, _ = findpeaks1d(pos_z, distance=(sr(obj) ÷ 2), height=mean(pos_z) + 2*std(pos_z))
    p_idx_acc_x, _ = findpeaks1d(acc_x, distance=(sr(obj) ÷ 2), height=mean(acc_x) + 2*std(acc_x))
    p_idx_acc_y, _ = findpeaks1d(acc_y, distance=(sr(obj) ÷ 2), height=mean(acc_y) + 2*std(acc_y))
    p_idx_acc_z, _ = findpeaks1d(acc_z, distance=(sr(obj) ÷ 2), height=mean(acc_z) + 2*std(acc_z))

    p_idx = round.(Int64, p_idx_x)

    _info("Detected pinches: $(length(p_idx))")

    return p_idx

end

"""
    tpt_analyze(nn_seg)

Analyze pinches in TPT recording.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple containing:
- `n::Int64`: number of pinches
- `t_mean::Float64`: mean interval between pinches [s]
- `t_median::Float64`: median interval between pinches [s]
- `t_rmssd::Float64`: ("root mean square of successive differences"), the square root of the mean of the squares of the successive differences between adjacent pinches
- `t_sdsd::Float64`: ("standard deviation of successive differences"), the standard deviation of the successive differences between adjacent pinches
"""
function tpt_analyze(obj::NeuroAnalyzer.NEURO)::@NamedTuple{n::Int64, t_mean::Float64, t_median::Float64, t_rmssd::Float64, t_sdsd::Float64}

    p_idx = tpt_detect(obj)
    t = obj.time_pts[p_idx]

    n = length(p_idx)
    t_diff = round.(diff(t), digits=2)
    t_mean = round(mean(t_diff), digits=2)
    t_median = round(median(t_diff), digits=2)
    t_rmssd = round(sqrt(mean(t_diff .^ 2)), digits=2)
    t_sdsd = round(std(t_diff), digits=2)

    return(n=n, t_mean=t_mean, t_median=t_median, t_rmssd=t_rmssd, t_sdsd=t_sdsd)

end