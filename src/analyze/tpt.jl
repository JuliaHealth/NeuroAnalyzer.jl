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

    p_idx1 = NeuroAnalyzer._tpt_peaks(obj.data[4, :, 1])
    p_idx2 = NeuroAnalyzer._tpt_peaks(obj.data[5, :, 1])
    p_idx3 = NeuroAnalyzer._tpt_peaks(.-obj.data[6, :, 1])
    p_idx = sort(union(p_idx1, p_idx2, p_idx3))
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
- `t_mean::Float64`: mean interval between pinches [ms]
- `t_median::Float64`: median interval between pinches [ms]
- `t_rmssd::Float64`: ("root mean square of successive differences"), the square root of the mean of the squares of the successive differences between adjacent pinches [ms]
- `t_sdsd::Float64`: ("standard deviation of successive differences"), the standard deviation of the successive differences between adjacent pinches [ms]

# Note

Return nothing if no pinches are detected.
"""
function tpt_analyze(obj::NeuroAnalyzer.NEURO)::Union{@NamedTuple{n::Int64, t_mean::Float64, t_median::Float64, t_rmssd::Float64, t_sdsd::Float64}, Nothing}

    p_idx = tpt_detect(obj)
    t = obj.time_pts[p_idx] .* 1000

    n = length(p_idx)
    if n == 1
        _warn("Only 1 pinch was detected, intervals cannot be calculated.")
        return nothing
    elseif n > 0
        t_diff = round.(diff(t), digits=1)
        t_mean = round(mean(t_diff), digits=1)
        t_median = round(median(t_diff), digits=1)
        t_rmssd = round(sqrt(mean(t_diff .^ 2)), digits=1)
        t_sdsd = round(std(t_diff), digits=1)
        return (n=n, t_mean=t_mean, t_median=t_median, t_rmssd=t_rmssd, t_sdsd=t_sdsd)
    else
        return nothing
    end


end
