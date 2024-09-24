export bootstrap_ci
export bootstrap_stat

"""
    bootstrap_ci(s; <keyword arguments>)

Calculate Confidence Interval using bootstrapping.

# Arguments

- `s::AbstractMatrix`: signal (time points × epochs)
- `n1::Int64=3000`: number of stage 1 resamplings -- number of samples in the resampled signal
- `n2::Int64=1000`: number of stage 2 resamplings -- number of samples sampled from the signal
- `ci::Float64=0.95`: confidence interval

# Returns

Named tuple containing:
- `s_avg::Vector{Float64}`: averaged signal
- `s_ci_l::Vector{Float64}`: lower bound of the confidence interval
- `s_ci_h::Vector{Float64}`: upper bound of the confidence interval
"""
function bootstrap_ci(s::AbstractMatrix; n1::Int64=3000, n2::Int64=1000, ci::Float64=0.95)

    @assert ci < 1.0 "ci must be < 1.0."

    # initialize progress bar
    progress_bar && (progbar = Progress(n1, dt=1, barlen=20, color=:white))

    s_avg = zeros(n1, size(s, 1))
    @inbounds for idx1 in 1:n1
        s_tmp = zeros(size(s, 1), n2)
        @inbounds for idx2 in 1:n2
            ep_idx = rand(1:size(s, 2))
            s_tmp[:, idx2] = @views s[:, ep_idx]
        end
        s_avg[idx1, :] = vec(mean(s_tmp, dims=2))

        # update progress bar
        progress_bar && next!(progbar)
    end

    s_ci_l = zeros(size(s_avg, 2))
    s_ci_h = zeros(size(s_avg, 2))

    ci_l = round((1.0 - 0.95) / 2, digits=3)
    ci_h = 1.0 - ci_l

    @inbounds for idx in axes(s_avg, 2)
        tpt = sort(@views s_avg[:, idx])
        ci_l_idx = round(Int64, ci_l * size(s_avg, 1))
        ci_h_idx = round(Int64, ci_h * size(s_avg, 1))
        s_ci_l[idx] = tpt[ci_l_idx]
        s_ci_h[idx] = tpt[ci_h_idx]
    end

    s_avg = vec(mean(s_avg, dims=1))

    return (s_avg=s_avg, s_ci_l=s_ci_l, s_ci_h=s_ci_h)

end

"""
    bootstrap_stat(s; <keyword arguments>)

Calculate signal statistic using bootstrapping.

# Arguments

- `s::AbstractMatrix`: signal (time points × epochs)
- `n1::Int64=3000`: number of stage 1 resamplings -- number of samples in the resampled signal
- `n2::Int64=1000`: number of stage 2 resamplings -- number of samples sampled from the signal
- `f::String`: statistic function to be applied, e.g. `f="abs(maximum(OBJ))"; signal is given using variable `OBJ` here.

# Returns

- `out::AbstractVector`
"""
function bootstrap_stat(s::AbstractMatrix; n1::Int64=3000, n2::Int64=1000, f::String)::AbstractVector

    f_tmp = replace(f, "OBJ" => "$(s[:, 1])")
    out_tmp = nothing
    try
        out_tmp = eval(Meta.parse(f_tmp))
    catch
        @error "Formula is incorrect."
    end
    out = zeros(typeof(out_tmp), n1)

    # initialize progress bar
    progress_bar && (progbar = Progress(n1, dt=1, barlen=20, color=:white))

    s_avg = zeros(n1, size(s, 1))
    @inbounds for idx1 in 1:n1
        s_tmp = zeros(size(s, 1), n2)
        @inbounds for idx2 in 1:n2
            ep_idx = rand(1:size(s, 2))
            s_tmp[:, idx2] = @views s[:, ep_idx]
        end
        s_avg[idx1, :] = mean(s_tmp, dims=2)
        f_tmp = replace(f, "OBJ" => "$(s_avg[idx1, :])")
        out[idx1] = eval(Meta.parse(f_tmp))

        # update progress bar
        progress_bar && next!(progbar)
    end

    return out

end
