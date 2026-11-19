export tpt_detect
export tpt_analyze

"""
    tpt_detect(obj)

Detect pinches in TPT recording: peaks of the movement envelope (3-D acceleration after removing single-sample glitches and the posture of the hand, smoothed so that the several acceleration peaks of one pinch are merged), at least 0.25 s apart. The detection threshold adapts to the noise and to the size of the pinches in the recording, so the acceleration units do not matter.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Vector{Int64}`: index of pinches locations
"""
function tpt_detect(obj::NeuroAnalyzer.NEURO)::Vector{Int64}
    _check_datatype(obj, "tpt")

    # acceleration channels: X, Y, Z (channels 4-6 of a TPT recording)
    p_idx = _tpt_peaks(obj.data[4:6, :, 1], sr(obj))
    _info("Detected pinches: $(length(p_idx))")

    return p_idx
end

"""
    tpt_analyze(obj)

Analyze pinches in TPT recording.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

Named tuple:

- `n::Int64`: number of pinches
- `t_mean::Float64`: mean interval between pinches [ms]
- `t_median::Float64`: median interval between pinches [ms]
- `t_rmssd::Float64`: ("root mean square of successive differences"), the square root of the mean of the squares of the successive differences between adjacent intervals [ms]; `NaN` for < 3 pinches
- `t_sdsd::Float64`: ("standard deviation of successive differences"), the standard deviation of the successive differences between adjacent intervals [ms]; `NaN` for < 4 pinches

# Note

Return nothing if no pinches are detected.
"""
function tpt_analyze(
    obj::NeuroAnalyzer.NEURO,
)::Union{
    @NamedTuple{
        n::Int64,
        t_mean::Float64,
        t_median::Float64,
        t_rmssd::Float64,
        t_sdsd::Float64,
    },
    Nothing,
}
    p_idx = tpt_detect(obj)
    t = obj.time_pts[p_idx] .* 1000

    n = length(p_idx)
    if n == 1
        _warn("Only 1 pinch was detected, intervals cannot be calculated.")
        return nothing
    elseif n > 0
        t_diff = diff(t)
        t_mean = mean(t_diff)
        t_median = median(t_diff)
        # successive differences of the intervals
        sd = diff(t_diff)
        t_rmssd = isempty(sd) ? NaN : sqrt(mean(abs2, sd))
        t_sdsd = length(sd) < 2 ? NaN : std(sd)
        return (; n, t_mean, t_median, t_rmssd, t_sdsd)
    else
        return nothing
    end
end
