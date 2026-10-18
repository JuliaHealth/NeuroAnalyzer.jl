# zscore(x) for arrays is StatsBase.zscore() (`(x .- mean(x)) ./ std(x)`), re-exported
export zscore

"""
    zscore(x, m, sd)

Calculate the Z-score of a single value given a known mean and standard deviation.

Computed as `(x − m) / sd`. This is a method of `StatsBase.zscore()`.

# Arguments

- `x::Real`: observed value
- `m::Real`: population or sample mean
- `sd::Real`: population or sample standard deviation; must not be zero

# Returns

- `Float64`: Z-score
"""
function StatsBase.zscore(x::Real, m::Real, sd::Real)::Float64
    # validate
    sd != 0 || throw(ArgumentError("sd must not be zero."))

    return (x - m) / sd
end
