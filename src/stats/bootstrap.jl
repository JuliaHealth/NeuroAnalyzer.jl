export bootstrap_ci
export bootstrap_stat

"""
    bootstrap_ci(s; <keyword arguments>)

Calculate a bootstrap confidence interval for a signal matrix.

Algorithm:
1. Stage 1 (`n1` iterations): randomly draw `n2` epochs (with replacement), average them → produces `n1` bootstrap mean traces.
2. Stage 2: for each time point, sort the `n1` bootstrap means and take the `cl/2` and `1 − cl/2` quantiles as the lower and upper CI bounds.

# Arguments

- `s::AbstractMatrix`: signal matrix, shape `(time_points, epochs)`
- `n1::Int64=3000`: number of bootstrap resamples (outer loop); must be ≥ 1
- `n2::Int64=1000`: number of epochs drawn per resample (inner loop); must be ≥ 1
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

Named tuple:

- `s_avg::Vector{Float64}`: bootstrap grand mean (averaged across all `n1` resamples)
- `s_ci_l::Vector{Float64}`: lower CI bound at each time point
- `s_ci_h::Vector{Float64}`: upper CI bound at each time point

# Throws

- `ArgumentError`: if `cl ∉ (0, 1)`, `n1 < 1`, or `n2 < 1`

# See also

[`bootstrap_stat`](@ref)
"""
function bootstrap_ci(
    s::AbstractMatrix;
    n1::Int64 = 3000,
    n2::Int64 = 1000,
    cl::Float64 = 0.95
)::@NamedTuple{
    s_avg::Vector{Float64},
    s_ci_l::Vector{Float64},
    s_ci_h::Vector{Float64}
}

    _bin(cl, (0.0, 1.0), "cl")
    @assert n1 > 0 "n1 must be > 0."
    @assert n2 > 0 "n2 must be > 0."

    tp_n = size(s, 1) # number of time points
    ep_n = size(s, 2) # number of epochs

    # stage 1: build n1 bootstrap mean traces of shape (n1 × tp_n)
    progbar = Progress(n1, dt=1, barlen=20, color=:white, enabled=progress_bar)

    s_boot = zeros(n1, tp_n)
    @inbounds for idx1 in 1:n1

        s_tmp = zeros(tp_n, n2)

        for idx2 in 1:n2
            s_tmp[:, idx2] = @view s[:, rand(1:ep_n)]
        end

        s_boot[idx1, :] = vec(mean(s_tmp; dims=2))
        progress_bar && next!(progbar)

    end

    # stage 2: derive CI bounds from the empirical quantiles of the bootstrap distribution
    s_ci_l = zeros(tp_n)
    s_ci_h = zeros(tp_n)
    # lower tail probability
    ci_l = (1.0 - cl) / 2
    # upper tail probability
    ci_h = 1.0 - ci_l

    # guard: ensure indices are at least 1 and at most n1
    ci_l_idx = max(1, round(Int64, ci_l * n1))
    ci_h_idx = min(n1, round(Int64, ci_h * n1))

    @inbounds for tp in axes(s_boot, 2)
        tpt = sort(@view s_boot[:, tp])
        s_ci_l[tp] = tpt[ci_l_idx]
        s_ci_h[tp] = tpt[ci_h_idx]
    end

    s_avg = vec(mean(s_boot, dims=1))

    return (s_avg=s_avg, s_ci_l=s_ci_l, s_ci_h=s_ci_h)

end

"""
    bootstrap_stat(s; <keyword arguments>)

Calculate a scalar statistic of a signal matrix using bootstrapping.

At each of the `n1` iterations, `n2` epochs are drawn at random (with replacement), averaged into a single trace, and the expression `f` is evaluated on that trace. The result is a distribution of the statistic across bootstrap resamples.

The formula string `f` must reference the current signal trace using the placeholder variable `obj`. For example: `f = "abs(maximum(obj))"`.

# Arguments

- `s::AbstractMatrix`: signal matrix, shape `(time_points, epochs)`
- `n1::Int64=3000`: number of bootstrap resamples; must be ≥ 1
- `n2::Int64=1000`: number of epochs drawn per resample; must be ≥ 1
- `f::String`: Julia expression to evaluate on each bootstrap mean trace; use `obj` as the placeholder for the current trace vector

# Returns

- `Vector`: bootstrap distribution of the statistic; length `n1`

# Throws
- `ArgumentError`: if `n1 < 1`, `n2 < 1`, or the formula `f` fails the dry-run evaluation

# See also

[`bootstrap_ci`](@ref)
"""
function bootstrap_stat(
    s::AbstractMatrix;
    n1::Int64 = 3000,
    n2::Int64 = 1000,
    f::String
)::AbstractVector

    @assert n1 > 0 "n1 must be > 0."
    @assert n2 > 0 "n2 must be > 0."

    tp_n = size(s, 1)
    ep_n = size(s, 2)

    # dry run on the first epoch to infer the output element type and validate f
    f_dry = replace(f, "obj" => "$(s[:, 1])")
    local out_tmp
    try
        out_tmp = eval(Meta.parse(f_dry))
    catch err
        throw(ArgumentError("Formula dry-run failed. Check expression `f`. Error: $err"))
    end

    out     = zeros(typeof(out_tmp), n1)
    s_boot  = zeros(n1, tp_n)

    # initialize progress bar
    progbar = Progress(n1, dt=1, barlen=20, color=:white, enabled=progress_bar)

    @inbounds for idx1 in 1:n1
        # draw n2 epochs with replacement and average to one trace
        s_tmp = zeros(tp_n, n2)
        for idx2 in 1:n2
            s_tmp[:, idx2] = @view s[:, rand(1:ep_n)]
        end
        s_boot[idx1, :] = vec(mean(s_tmp; dims=2))

        # evaluate the user formula on this bootstrap mean trace
        f_tmp = replace(f, "obj" => "$(s_boot[idx1, :])")
        try
            out[idx1] = eval(Meta.parse(f_tmp))
        catch err
            error("Formula failed at resample $idx1. Check expression `f`. Error: $err")
        end

        progress_bar && next!(progbar)
    end

    return out

end
