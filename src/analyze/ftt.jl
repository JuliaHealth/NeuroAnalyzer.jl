export ftt_analyze

"""
    ftt_analyze(t)

Analyze taps in a TPT (Tapping Performance Task) recording. Summarizes inter-tap interval statistics from a TPT (Tapping Performance Task) recording.

# Arguments

Named tuple:

- `taps::Vector{Int64}`: number of taps per trial
- `tap_t::Vector{Vector{Float64}}`: tap time points per trial [ms]
- `tap_d::Vector{Vector{Float64}}`: tap durations per trial [ms]
- `taps_int::Vector{Int64}`: number of taps per trial during intervals
- `tap_t_int::Vector{Vector{Float64}}`: tap time points during intervals [ms]
- `tap_d_int::Vector{Vector{Float64}}`: tap durations during intervals [ms]

# Returns

Named tuple:

- `n::Int64`: total number of taps across all trials
- `t_mean::Float64`: mean inter-tap interval [ms]
- `t_median::Float64`: median inter-tap interval [ms]
- `t_rmssd::Float64`: root mean square of successive differences between adjacent inter-tap intervals [ms]: `√(mean( diff(ITI)²))`
- `t_sdsd::Float64`: standard deviation of successive differences [ms]: `std(diff(ITI))`

Intervals and their successive differences are computed within trials; values are `NaN` when there are too few taps.
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
    t_sdsd::Float64,
}
    # total tap count across all trials
    n = sum(t.taps)

    # inter-tap intervals (ITI) within each trial (a trial with < 2 taps has none) and
    # their successive differences within each trial (not across trial boundaries)
    t_iti = Float64[]
    sd = Float64[]
    for idx in axes(t.tap_t, 1)
        length(t.tap_t[idx]) > 1 || continue
        iti = diff(t.tap_t[idx])
        append!(t_iti, iti)
        append!(sd, diff(iti))
    end

    t_mean = isempty(t_iti) ? NaN : mean(t_iti)
    t_median = isempty(t_iti) ? NaN : median(t_iti)
    t_rmssd = isempty(sd) ? NaN : sqrt(mean(abs2, sd))
    t_sdsd = length(sd) < 2 ? NaN : std(sd)

    return (; n, t_mean, t_median, t_rmssd, t_sdsd)
end
