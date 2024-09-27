export epoch_stats
export channel_stats

"""
    epoch_stats(obj)

Calculate epochs statistics.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple containing:
- `e_mean::Vector{Float64}`: mean
- `e_median::Vector{Float64}`: median
- `e_std::Vector{Float64}`: standard deviation
- `e_var::Vector{Float64}`: variance
- `e_kurt::Vector{Float64}`: kurtosis
- `e_skew::Vector{Float64}`: skewness
- `e_mean_diff::Vector{Float64}`: mean diff value
- `e_median_diff::Vector{Float64}`: median diff value
- `e_max_dif::Vector{Float64}`: max difference
- `e_dev_mean::Vector{Float64}`: deviation from channel mean
"""
function epoch_stats(obj::NeuroAnalyzer.NEURO)::NamedTuple{(:e_mean, :e_median, :e_std, :e_var, :e_kurt, :e_skew, :e_mean_diff, :e_median_diff, :e_max_dif, :e_dev_mean), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}}

    ep_n = nepochs(obj)

    e_mean = zeros(ep_n)
    e_median = zeros(ep_n)
    e_std = zeros(ep_n)
    e_var = zeros(ep_n)
    e_kurt = zeros(ep_n)
    e_skew = zeros(ep_n)
    e_mean_diff = zeros(ep_n)
    e_median_diff = zeros(ep_n)
    e_max_dif = zeros(ep_n)
    e_dev_mean = zeros(ep_n)

    @inbounds for ep_idx in 1:ep_n
        e_mean[ep_idx] = @views mean(obj.data[:, :, ep_idx])
        e_median[ep_idx] = @views median(obj.data[:, :, ep_idx])
        e_std[ep_idx] = @views std(obj.data[:, :, ep_idx])
        e_var[ep_idx] = @views var(obj.data[:, :, ep_idx])
        e_kurt[ep_idx] = @views kurtosis(obj.data[:, :, ep_idx])
        e_skew[ep_idx] = @views skewness(obj.data[:, :, ep_idx])
        e_mean_diff[ep_idx] = @views mean(diff(obj.data[:, :, ep_idx], dims=2))
        e_median_diff[ep_idx] = @views median(diff(obj.data[:, :, ep_idx], dims=2))
        e_max_dif[ep_idx] = @views maximum(obj.data[:, :, ep_idx]) - minimum(obj.data[:, :, ep_idx])
        e_dev_mean[ep_idx] = @views abs(mean(obj.data[:, :, ep_idx])) - mean(obj.data[:, :, ep_idx])
    end

    return (e_mean=e_mean, e_median=e_median, e_std=e_std, e_var=e_var, e_kurt=e_kurt, e_skew=e_skew, e_mean_diff=e_mean_diff, e_median_diff=e_median_diff, e_max_dif=e_max_dif, e_dev_mean=e_dev_mean)

end

"""
    channel_stats(obj)

Calculate channels statistics per epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple containing:
- `c_mean::Matrix{Float64}`: mean
- `c_median::Matrix{Float64}`: median
- `c_std::Matrix{Float64}`: standard deviation
- `c_var::Matrix{Float64}`: variance
- `c_kurt::Matrix{Float64}`: kurtosis
- `c_skew::Matrix{Float64}`: skewness
- `c_mean_diff::Matrix{Float64}`: mean diff value
- `c_median_diff::Matrix{Float64}`: median diff value
- `c_max_dif::Matrix{Float64}`: max difference
- `c_dev_mean::Matrix{Float64}`: deviation from channel mean
"""
function channel_stats(obj::NeuroAnalyzer.NEURO)::NamedTuple{(:c_mean, :c_median, :c_std, :c_var, :c_kurt, :c_skew, :c_mean_diff, :c_median_diff, :c_max_dif, :c_dev_mean), Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}

    ch_n = nchannels(obj)
    ep_n = nepochs(obj)

    c_mean = zeros(ch_n, ep_n)
    c_median = zeros(ch_n, ep_n)
    c_std = zeros(ch_n, ep_n)
    c_var = zeros(ch_n, ep_n)
    c_kurt = zeros(ch_n, ep_n)
    c_skew = zeros(ch_n, ep_n)
    c_mean_diff = zeros(ch_n, ep_n)
    c_median_diff = zeros(ch_n, ep_n)
    c_max_dif = zeros(ch_n, ep_n)
    c_dev_mean = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            c_mean[ch_idx, ep_idx] = @views mean(obj.data[ch_idx, :, ep_idx])
            c_median[ch_idx, ep_idx] = @views median(obj.data[ch_idx, :, ep_idx])
            c_std[ch_idx, ep_idx] = @views std(obj.data[ch_idx, :, ep_idx])
            c_var[ch_idx, ep_idx] = @views var(obj.data[ch_idx, :, ep_idx])
            c_kurt[ch_idx, ep_idx] = @views kurtosis(obj.data[ch_idx, :, ep_idx])
            c_skew[ch_idx, ep_idx] = @views skewness(obj.data[ch_idx, :, ep_idx])
            c_mean_diff[ch_idx, ep_idx] = @views mean(diff(obj.data[ch_idx, :, ep_idx]))
            c_median_diff[ch_idx, ep_idx] = @views median(diff(obj.data[ch_idx, :, ep_idx]))
            c_max_dif[ch_idx, ep_idx] = @views maximum(obj.data[ch_idx, :, ep_idx]) - minimum(obj.data[ch_idx, :, ep_idx])
            c_dev_mean[ch_idx, ep_idx] = @views abs(mean(obj.data[ch_idx, :, ep_idx])) - mean(obj.data[ch_idx, :, ep_idx])
        end
    end

    return (c_mean=c_mean, c_median=c_median, c_std=c_std, c_var=c_var, c_kurt=c_kurt, c_skew=c_skew, c_mean_diff=c_mean_diff, c_median_diff=c_median_diff, c_max_dif=c_max_dif, c_dev_mean=c_dev_mean)

end
