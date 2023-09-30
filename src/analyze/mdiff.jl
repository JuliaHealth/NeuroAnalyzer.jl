export mdiff

"""
    mdiff(s1, s2; n, method)

Calculate mean difference and 95% confidence interval for 2 signals.

# Arguments

- `s1::AbstractMatrix`
- `s2::AbstractMatrix`
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:absdiff`:
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `st::Vector{Float64}`
- `sts::Float64`
- `p::Float64`
"""
function mdiff(s1::AbstractMatrix, s2::AbstractMatrix; n::Int64=3, method::Symbol=:absdiff)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _check_var(method, [:absdiff, :diff2int], "method")
    @assert n >= 1 "n must be ≥ 1."

    s1_mean = vec(mean(s1, dims=1))
    s2_mean = vec(mean(s2, dims=1))

    if method === :absdiff
        # statistic: maximum difference
        s_diff = s1_mean - s2_mean
        sts = maximum(abs.(s_diff))
    else
        # statistic: integrated area of the squared difference
        s_diff_squared = (s1_mean - s2_mean).^2
        sts = simpson(s_diff_squared)
    end

    ss = [s1; s2]
    st = zeros(size(s1, 1) * n)

    Threads.@threads for idx1 in 1:(size(s1, 1) * n)
        s_tmp1 = zeros(size(s1, 1), size(s1, 2))
        sample_idx = rand(1:size(ss, 1), size(ss, 1))
        @inbounds @simd for idx2 in 1:size(s1, 1)
            s_tmp1[idx2, :] = @views ss[sample_idx[idx2], :]'
        end
        s1_mean = vec(mean(s_tmp1, dims=1))
        s_tmp1 = zeros(size(s1, 1), size(s1, 2))
        sample_idx = rand(1:size(ss, 1), size(ss, 1))
        @inbounds @simd for idx2 in 1:size(s1, 1)
            s_tmp1[idx2, :] = @views ss[sample_idx[idx2], :]'
        end
        s2_mean = vec(mean(s_tmp1, dims=1))
        if method === :absdiff
            # statistic: maximum difference
            s_diff = s1_mean - s2_mean
            @inbounds st[idx1] = maximum(abs.(s_diff))
        else
            # statistic: integrated area of the squared difference
            s_diff_squared = (s1_mean - s2_mean).^2
            @inbounds st[idx1] = simpson(s_diff_squared)
        end
    end

    p = length(st[st .> sts]) / size(s1, 1) * n
    p > 1 && (p = 1.0)

    return (st=st, sts=sts, p=p)

end

"""
    mdiff(s1, s2; n, method)

Calculate mean difference and its 95% CI between channels.

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:absdiff`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `st::Matrix{Float64}`
- `sts::Vector{Float64}`
- `p::Vector{Float64}`
"""
function mdiff(s1::AbstractArray, s2::AbstractArray; n::Int64=3, method::Symbol=:absdiff)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    st = zeros(ep_n, ch_n * n)
    sts = zeros(ep_n)
    p = zeros(ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        st[ep_idx, :], sts[ep_idx], p[ep_idx] = mdiff(s1[:, :, ep_idx], s2[:, :, ep_idx], n=n, method=method)
    end

    return (st=st, sts=sts, p=p)
end

"""
    mdiff(obj; ch, n, method)

Calculate mean difference and its 95% CI between channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:absdiff`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `st::Matrix{Float64}`
- `sts::Vector{Float64}`
- `p::Vector{Float64}`
"""
function mdiff(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), n::Int64=3, method::Symbol=:absdiff)

    _check_channels(obj, ch)

    st, sts, p = mdiff(obj.data[ch, :, :], n=n, method=method)

    return (st=st, sts=sts, p=p)

end

"""
    mdiff(obj1, obj2; channel1, channel2, epoch1, epoch2, n, method)

Calculates mean difference and 95% confidence interval for two channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2:NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `channel2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `epoch1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `n::Int64`: number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `st::Matrix{Float64}`
- `sts::Vector{Float64}`
- `p::Vector{Float64}`
"""
function mdiff(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1), channel2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2), epoch1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1)), epoch2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2)), n::Int64=3, method::Symbol=:absdiff)

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    @assert length(channel1) == length(channel2) "ch1 and ch2 must have the same length."
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    @assert length(epoch1) == length(epoch2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    st, sts, p = @views mdiff(obj1.data[channel1, :, epoch1], obj2.data[channel2, :, epoch2], n=n, method=method)

    return (st=st, sts=sts, p=p)
end
