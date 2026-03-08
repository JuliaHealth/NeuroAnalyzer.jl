export ftt_analyze

"""
    ftt_analyze(t)

Analyze taps in a TPT (Tapping Performance Task) recording. Summarizes inter-tap interval statistics from a TPT (Tapping Performance Task) recording.

# Arguments

Named tuple containing:

- `taps::Vector{Int64}`: number of taps per trial
- `tap_t::Vector{Vector{Float64}}`: tap time points per trial [ms]
- `tap_d::Vector{Vector{Float64}}`: tap durations per trial [ms]
- `taps_int::Vector{Int64}`: number of taps per trial during intervals
- `tap_t_int::Vector{Vector{Float64}}`: tap time points during intervals [ms]
- `tap_d_int::Vector{Vector{Float64}}`: tap durations during intervals [ms]

# Returns

Named tuple containing:

- `n::Int64`: total number of taps across all trials
- `t_mean::Float64`: mean inter-tap interval [ms]
- `t_median::Float64`: median inter-tap interval [ms]
- `t_rmssd::Float64`: root mean square of successive differences between adjacent inter-tap intervals [ms]: `√(mean( diff(ITI)²))`
- `t_sdsd::Float64`: standard deviation of successive differences [ms]: `std(diff(ITI))`
"""
function ftt_analyze(
    t::@NamedTuple{
        taps::Vector{Int64},
        tap_t::Vector{Vector{Float64}},
        tap_d::Vector{Vector{Float64}},
        taps_int::Vector{Int64},
        tap_t_int::Vector{Vector{Float64}},
        tap_d_int::Vector{Vector{Float64}},
    }
)::@NamedTuple{
    n::Int64,
    t_mean::Float64,
    t_median::Float64,
    t_rmssd::Float64,
    t_sdsd::Float64
}

    # total tap count across all trials
    n = sum(t.taps)

    # accumulate inter-tap intervals (ITI) across all trials
    # a trial with 1 tap contributes that tap's timestamp (no ITI computable)
    # trials with ≥ 2 taps contribute their successive differences.
    t_iti = Float64[]
    for idx in axes(t.tap_t, 1)
        if length(t.tap_t[idx]) == 1
            push!(t_iti, t.tap_t[idx][1])
        elseif length(t.tap_t[idx]) > 1
            append!(t_iti, round.(diff(t.tap_t[idx]), digits = 1))
        end
    end

    t_mean = round(mean(t_iti), digits = 1)
    t_median = round(median(t_iti), digits = 1)

    # successive differences of the ITI vector (difference-of-differences)
    sd = diff(t_iti)

    # sum(abs2, x) = Σ xᵢ² — avoids allocating x .^ 2 array.
    t_rmssd = round(sqrt(sum(abs2, sd) / length(sd)), digits = 1)
    t_sdsd  = round(std(sd),                          digits = 1)

    return (n = n, t_mean = t_mean, t_median = t_median, t_rmssd = t_rmssd, t_sdsd = t_sdsd)

end
