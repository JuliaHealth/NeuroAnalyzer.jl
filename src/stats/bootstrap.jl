export bootstrap_ci
export bootstrap_stat

# mean trace of n epochs drawn at random (with replacement) from s (samples × epochs)
_bootstrap_mean(s::AbstractMatrix, n::Int64)::Vector{Float64} =
    vec(mean(@view(s[:, rand(axes(s, 2), n)]); dims = 2))

"""
    bootstrap_ci(s; <keyword arguments>)

Calculate a bootstrap confidence interval for a signal matrix.

Algorithm:
1. Stage 1 (`n1` iterations): randomly draw `n2` epochs (with replacement), average them → produces `n1` bootstrap mean traces.
2. Stage 2: for each time point, sort the `n1` bootstrap means and take the `cl/2` and `1 − cl/2` quantiles as the lower and upper CI bounds.

# Arguments

- `s::AbstractMatrix`: signal matrix, shape (samples, epochs)
- `n1::Int64=3000`: number of bootstrap resamples (outer loop); must be ≥ 1
- `n2::Int64=1000`: number of epochs drawn per resample (inner loop); must be ≥ 1
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

Named tuple:

- `gm::Vector{Float64}`: bootstrap grand mean (averaged across all `n1` resamples)
- `ll::Vector{Float64}`: lower CI bound at each time point
- `ul::Vector{Float64}`: upper CI bound at each time point
"""
function bootstrap_ci(
    s::AbstractMatrix;
    n1::Int64 = 3000,
    n2::Int64 = 1000,
    cl::Float64 = 0.95,
)::@NamedTuple{
    gm::Vector{Float64},
    ll::Vector{Float64},
    ul::Vector{Float64},
}
    # validate
    _bin(cl, (0.0, 1.0), "cl")
    n1 > 0 || throw(ArgumentError("n1 must be > 0."))
    n2 > 0 || throw(ArgumentError("n2 must be > 0."))

    # number of time points
    tp_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 2)

    # stage 1: build n1 bootstrap mean traces of shape (n1 × tp_n)
    progbar = Progress(n1; dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    s_boot = zeros(n1, tp_n)
    @inbounds for idx1 = 1:n1
        s_boot[idx1, :] = _bootstrap_mean(s, n2)
        progress_bar && next!(progbar)
    end

    # stage 2: derive CI bounds from the empirical quantiles of the bootstrap distribution
    ci_l = (1.0 - cl) / 2
    ll = [quantile(@view(s_boot[:, tp]), ci_l) for tp in axes(s_boot, 2)]
    ul = [quantile(@view(s_boot[:, tp]), 1.0 - ci_l) for tp in axes(s_boot, 2)]

    gm = vec(mean(s_boot; dims = 1))

    return (; gm, ll, ul)
end

"""
    bootstrap_stat(s; <keyword arguments>)

Calculate a scalar statistic of a signal matrix using bootstrapping.

At each of the `n1` iterations, `n2` epochs are drawn at random (with replacement), averaged into a single trace, and the expression `f` is evaluated on that trace. The result is a distribution of the statistic across bootstrap resamples.

`f` is a function of the trace (e.g. `f = x -> abs(maximum(x))`) or a string holding a Julia expression in which the trace is the variable `obj` (e.g. `f = "abs(maximum(obj))"`).

# Arguments

- `s::AbstractMatrix`: signal matrix, shape (samples, epochs)
- `n1::Int64=3000`: number of bootstrap resamples; must be ≥ 1
- `n2::Int64=1000`: number of epochs drawn per resample; must be ≥ 1
- `f::Union{Function, String}`: statistic evaluated on each bootstrap mean trace: a function, or a Julia expression using `obj` as the current trace vector

# Returns

- `AbstractVector`: bootstrap distribution of the statistic; length `n1`
"""
function bootstrap_stat(
    s::AbstractMatrix;
    n1::Int64 = 3000,
    n2::Int64 = 1000,
    f::Union{Function, String},
)::AbstractVector
    # validate
    n1 > 0 || throw(ArgumentError("n1 must be > 0."))
    n2 > 0 || throw(ArgumentError("n2 must be > 0."))

    # a string expression is compiled once into a function of `obj`
    fn = if f isa String
        ex = try
            Meta.parse(f)
        catch err
            throw(ArgumentError("Cannot parse expression `f`: $err"))
        end
        eval(:(obj -> $ex))
    else
        f
    end
    apply_f(x) = Base.invokelatest(fn, x)

    # dry run on the first epoch to infer the output type and validate f
    s_boot_tmp = try
        apply_f(s[:, 1])
    catch err
        throw(ArgumentError("Formula dry-run failed. Check expression `f`. Error: $err"))
    end
    result = Vector{typeof(s_boot_tmp)}(undef, n1)

    # initialize progress bar
    progbar = Progress(n1; dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    for idx1 = 1:n1
        # evaluate the statistic on the mean of n2 epochs drawn with replacement
        result[idx1] = apply_f(_bootstrap_mean(s, n2))
        progress_bar && next!(progbar)
    end

    return result
end
