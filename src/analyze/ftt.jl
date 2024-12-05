export ftt_analyze

"""
    ftt_analyze(t)

Analyze taps in TPT recording.

# Arguments

Named tuple containing:
- `taps::Vector{Int64}`: number of taps per trial
- `tap_t::Vector{Vector{Float64}}`: taps time point [ms]
- `tap_d::Vector{Vector{Float64}}`: taps duration [ms]
- `taps_int::Vector{Int64}`: number of taps per trial during intervals
- `tap_t_int::Vector{Vector{Float64}}`: taps time point [ms] during intervals
- `tap_d_int::Vector{Vector{Float64}}`: taps duration [ms] during intervals`

# Returns

Named tuple containing:
- `n::Int64`: number of taps
- `t_mean::Float64`: mean interval between taps [s]
- `t_median::Float64`: median interval between taps [s]
- `t_rmssd::Float64`: ("root mean square of successive differences"), the square root of the mean of the squares of the successive differences between adjacent taps
- `t_sdsd::Float64`: ("standard deviation of successive differences"), the standard deviation of the successive differences between adjacent taps
"""
function ftt_analyze(t::@NamedTuple{taps::Vector{Int64}, tap_t::Vector{Vector{Float64}}, tap_d::Vector{Vector{Float64}}, taps_int::Vector{Int64}, tap_t_int::Vector{Vector{Float64}}, tap_d_int::Vector{Vector{Float64}}})::@NamedTuple{n::Int64, t_mean::Float64, t_median::Float64, t_rmssd::Float64, t_sdsd::Float64}

    n = sum(t.taps)
    t_diff = Float64[]
    for idx in axes(t.tap_t, 1)
        if length(t.tap_t[idx]) == 1
            push!(t_diff, t.tap_t[idx][1])
        elseif length(t.tap_t[idx]) > 1
            append!(t_diff, round.(diff(t.tap_t[idx]), digits=2))
        end
    end
    t_mean = round(mean(t_diff), digits=2)
    t_median = round(median(t_diff), digits=2)
    t_rmssd = round(sqrt(mean(t_diff .^ 2)), digits=2)
    t_sdsd = round(std(t_diff), digits=2)

    return(n=n, t_mean=t_mean, t_median=t_median, t_rmssd=t_rmssd, t_sdsd=t_sdsd)

end