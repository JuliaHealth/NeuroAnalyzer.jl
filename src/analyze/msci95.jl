export msci95

"""
    msci95(s)

Calculate mean, std and 95% confidence interval.

# Arguments

- `s::Vector{Float64}`

# Returns

Named tuple containing:
- `s_m::Float64`: mean
- `s_s::Float64`: standard deviation
- `s_u::Float64`: upper 95% CI
- `s_l::Float64`: lower 95% CI
"""
function msci95(s::Vector{Float64})
    s_m = mean(s)
    s_s = std(s) / sqrt(length(s))
    s_u = s_m + 1.96 * s_s
    s_l = s_m - 1.96 * s_s
    return (s_m=s_m, s_s=s_s, s_u=s_u, s_l=s_l)
end

"""
    msci95(s; n, method)

Calculate mean, std and 95% confidence interval for each the channel.

# Arguments

- `s::AbstractArray`
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal method (`:normal`) or `n`-times boostrapping (`:boot`)

# Returns

Named tuple containing:
- `s_m::Vector{Float64}`: mean
- `s_s::Vector{Float64}`: standard deviation
- `s_u::Vector{Float64}`: upper 95% CI
- `s_l::Vector{Float64}`: lower 95% CI
"""
function msci95(s::AbstractArray; n::Int64=3, method::Symbol=:normal)

    _check_var(method, [:normal, :boot], "method")
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    if method === :normal
        s_m = mean(s, dims=1)'
        s_s = std(s, dims=1)' / sqrt(size(s, 1))
        s_u = s_m + 1.96 * s_s
        s_l = s_m - 1.96 * s_s
    else
        s_tmp1 = zeros(size(s, 1) * n, size(s, 2))
        Threads.@threads for idx1 in 1:size(s, 1) * n
            s_tmp2 = zeros(size(s))
            sample_idx = rand(1:size(s, 1), size(s, 1))
            @inbounds @simd for idx2 in 1:size(s, 1)
                s_tmp2[idx2, :] = s[sample_idx[idx2], :]'
            end
            s_tmp1[idx1, :] = mean(s_tmp2, dims=1)
        end

        s_m = mean(s_tmp1, dims=1)'
        s_s = std(s_tmp1, dims=1)' / sqrt(size(s_tmp1, 1))
        s_sorted = sort(s_tmp1, dims=1)
        s_l = s_sorted[round(Int, 0.025 * size(s_tmp1, 1)), :]
        s_u = s_sorted[round(Int, 0.975 * size(s_tmp1, 1)), :]
    end

    return (s_m=vec(s_m[:, 1]), s_s=vec(s_s[:, 1]), s_u=vec(s_u[:, 1]), s_l=vec(s_l[:, 1]))
end

"""
    msci95(s1, s2)

Calculate mean and 95% confidence interval for 2 signals.

# Arguments

- `s1::Vector{Float64}`
- `s2:Vector{Float64}`

# Returns

Named tuple containing:
- `s_m::Float64`: mean
- `s_s::Float64`: standard deviation
- `s_u::Float64`: upper 95% CI
- `s_l::Float64`: lower 95% CI
"""
function msci95(s1::Vector{Float64}, s2::Vector{Float64})

    length(s1) == length(s2) || throw(ArgumentError("s1 and s2 must have the same length."))

    s_m = zeros(length(s1))
    s_s = zeros(length(s1))
    s_u = zeros(length(s1))
    s_l = zeros(length(s1))

    s1_mean = mean(s1)
    s2_mean = mean(s2)
    s_m = s1_mean - s2_mean
    s1_sd = std(s1) / sqrt(length(s1))
    s2_sd = std(s2) / sqrt(length(s2))
    s_s = sqrt(s1_sd^2 + s2_sd^2)
    s_u = s_m + 1.96 * s_s
    s_l = s_m - 1.96 * s_s

    return (s_m=s_m, s_s=s_s, s_u=s_u, s_l=s_l)
end

"""
    msci95(obj; ch, n, method)

Calculate mean, standard deviation and 95% confidence interval for channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal (`:normal`) method or `n`-times bootstrapping (`:boot`)

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean
- `s_s::Matrix{Float64}`: standard deviation
- `s_u::Matrix{Float64}`: upper 95% CI
- `s_l::Matrix{Float64}`: lower 95% CI
"""
function msci95(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), n::Int64=3, method::Symbol=:normal)

    _check_var(method, [:normal, :boot], "method")

    _check_channels(obj, ch)
    ep_len = epoch_len(obj)
    ep_n = epoch_n(obj)

    s_m = zeros(ep_n, ep_len)
    s_s = zeros(ep_n, ep_len)
    s_u = zeros(ep_n, ep_len)
    s_l = zeros(ep_n, ep_len)

    Threads.@threads for ep_idx in 1:ep_n
        s_m[ep_idx, :], s_s[ep_idx, :], s_u[ep_idx, :], s_l[ep_idx, :] = @views msci95(obj.data[ch, :, ep_idx], n=n, method=method)
    end

    return (mean=s_m, sd=s_s, upper=s_u, lower=s_l)
end

"""
    msci95(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculates mean, standard deviation and 95% confidence interval for two channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2:NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean by epochs
- `s_s::Matrix{Float64}`: standard deviation by epochs
- `s_u::Matrix{Float64}`: upper 95% CI bound by epochs
- `s_l::Matrix{Float64}`: lower 95% CI bound by epochs
"""
function msci95(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    # check channels
    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    # check epochs
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    ep_n = length(epoch1)
    ep_len = epoch_len(obj1)

    s_m = zeros(ep_n, ep_len)
    s_s = zeros(ep_n, ep_len)
    s_u = zeros(ep_n, ep_len)
    s_l = zeros(ep_n, ep_len)

    Threads.@threads for ep_idx in 1:ep_n
        s1_mean = @views mean(obj1.data[channel1, :, epoch1[ep_idx]], dims=1)
        s2_mean = @views mean(obj2.data[channel2, :, epoch2[ep_idx]], dims=1)
        s_m[ep_idx, :] = s1_mean - s2_mean
        s1_sd = @views std(obj1.data[channel1, :, epoch1[ep_idx]], dims=1) / sqrt(size(obj1.data[channel1, :, epoch1[ep_idx]], 2))
        s2_sd = @views std(obj2.data[channel2, :, epoch2[ep_idx]], dims=1) / sqrt(size(obj2.data[channel2, :, epoch2[ep_idx]], 2))
        s_s[ep_idx, :] = sqrt.(s1_sd.^2 .+ s2_sd.^2)
        s_u[ep_idx, :] = @. s_m[ep_idx, :] + 1.96 * s_s[ep_idx, :]
        s_l[ep_idx, :] = @. s_m[ep_idx, :] - 1.96 * s_s[ep_idx, :]
    end

    return (s_m=s_m, s_s=s_s, s_u=s_u, s_l=s_l)
end
