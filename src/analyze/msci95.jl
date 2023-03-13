export msci95

"""
    msci95(s)

Calculate mean, standard deviation and 95% confidence interval.

# Arguments

- `s::AbstractVector`

# Returns

Named tuple containing:
- `sm::Float64`: mean
- `ss::Float64`: standard deviation
- `su::Float64`: upper 95% CI
- `sl::Float64`: lower 95% CI
"""
function msci95(s::AbstractVector)

    sm = mean(s)
    ss = std(s) / sqrt(length(s))
    su = sm + 1.96 * ss
    sl = sm - 1.96 * ss

    return (sm=sm, ss=ss, su=su, sl=sl)

end

"""
    msci95(s; n, method)

Calculate mean, standard deviation and 95% confidence interval.

# Arguments

- `s::AbstractMatrix`
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal method (`:normal`) or `n`-times boostrapping (`:boot`)

# Returns

Named tuple containing:
- `sm::Vector{Float64}`: mean
- `ss::Vector{Float64}`: standard deviation
- `su::Vector{Float64}`: upper 95% CI
- `sl::Vector{Float64}`: lower 95% CI
"""
function msci95(s::AbstractMatrix; n::Int64=3, method::Symbol=:normal)

    _check_var(method, [:normal, :boot], "method")
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    if method === :normal
        sm = mean(s, dims=1)'
        ss = std(s, dims=1)' / sqrt(size(s, 1))
        su = sm + 1.96 * ss
        sl = sm - 1.96 * ss
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

        sm = mean(s_tmp1, dims=1)'
        ss = std(s_tmp1, dims=1)' / sqrt(size(s_tmp1, 1))
        ssorted = sort(s_tmp1, dims=1)
        sl = ssorted[round(Int, 0.025 * size(s_tmp1, 1)), :]
        su = ssorted[round(Int, 0.975 * size(s_tmp1, 1)), :]
    end

    return (sm=vec(sm[:, 1]), ss=vec(ss[:, 1]), su=vec(su[:, 1]), sl=vec(sl[:, 1]))

end

"""
    msci95(s; n, method)

Calculate mean, standard deviation and 95% confidence interval.

# Arguments

- `s::AbstractArray`
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal method (`:normal`) or `n`-times boostrapping (`:boot`)

# Returns

Named tuple containing:
- `sm::Array{Float64}`: mean
- `ss::Array{Float64}`: standard deviation
- `su::Array{Float64}`: upper 95% CI
- `sl::Array{Float64}`: lower 95% CI
"""
function msci95(s::AbstractArray; n::Int64=3, method::Symbol=:normal)

    _check_var(method, [:normal, :boot], "method")

    ep_len = size(s, 2)
    ep_n = size(s, 3)
    
    sm = zeros(ep_n, ep_len)
    ss = zeros(ep_n, ep_len)
    su = zeros(ep_n, ep_len)
    sl = zeros(ep_n, ep_len)

    Threads.@threads for ep_idx in 1:ep_n
        @inbounds sm[ep_idx, :], ss[ep_idx, :], su[ep_idx, :], sl[ep_idx, :] = @views msci95(s[:, :, ep_idx], n=n, method=method)
    end

    return (sm=sm, ss=ss, su=su, sl=sl)

end

"""
    msci95(obj; ch, n, method)

Calculate mean, standard deviation and 95% confidence interval.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal (`:normal`) method or `n`-times bootstrapping (`:boot`)

# Returns

Named tuple containing:
- `sm::Matrix{Float64}`: mean
- `ss::Matrix{Float64}`: standard deviation
- `su::Matrix{Float64}`: upper 95% CI
- `sl::Matrix{Float64}`: lower 95% CI
"""
function msci95(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), n::Int64=3, method::Symbol=:normal)

    _check_channels(obj, ch)

    sm, ss, su, sl = @views msci95(obj.data, n=n, method=method)

    return (sm=sm, ss=ss, su=su, sl=sl)

end

"""
    msci95(s1, s2)

Calculate mean difference, standard deviation and 95% confidence interval.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

Named tuple containing:
- `sm::Float64`: mean
- `ss::Float64`: standard deviation
- `su::Float64`: upper 95% CI
- `sl::Float64`: lower 95% CI
"""
function msci95(s1::AbstractVector, s2::AbstractVector)

    length(s1) == length(s2) || throw(ArgumentError("s1 and s2 must have the same length."))

    sm = zeros(length(s1))
    ss = zeros(length(s1))
    su = zeros(length(s1))
    sl = zeros(length(s1))

    s1_mean = mean(s1)
    s2_mean = mean(s2)
    sm = s1_mean - s2_mean
    s1_sd = std(s1) / sqrt(length(s1))
    s2_sd = std(s2) / sqrt(length(s2))
    ss = sqrt(s1_sd^2 + s2_sd^2)
    su = sm + 1.96 * ss
    sl = sm - 1.96 * ss

    return (sm=sm, ss=ss, su=su, sl=sl)

end

"""
    msci95(s1, s2)

Calculate mean difference, standard deviation and 95% confidence interval.

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`

# Returns

Named tuple containing:
- `sm::Array{Float64}`: mean
- `ss::Array{Float64}`: standard deviation
- `su::Array{Float64}`: upper 95% CI
- `sl::Array{Float64}`: lower 95% CI
"""
function msci95(s1::AbstractArray, s2::AbstractArray; n::Int64=3, method::Symbol=:normal)

    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))
    _check_var(method, [:normal, :boot], "method")

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    sm = zeros(ch_n, ep_n)
    ss = zeros(ch_n, ep_n)
    su = zeros(ch_n, ep_n)
    sl = zeros(ch_n, ep_n)

    Threads.@threads for ep_idx in 1:ep_n
        @inbounds @simd for ch_idx in 1:ch_n
            sm[ch_idx, ep_idx], ss[ch_idx, ep_idx], su[ch_idx, ep_idx], sl[ch_idx, ep_idx] = @views msci95(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx])
        end
    end

    return (sm=sm, ss=ss, su=su, sl=sl)

end

"""
    msci95(obj1, obj2; ch1, ch2, ep1, ep2)

Calculate mean difference, standard deviation and 95% confidence interval.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2:NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `sm::Matrix{Float64}`: mean
- `ss::Matrix{Float64}`: standard deviation
- `su::Matrix{Float64}`: upper 95% CI bound
- `sl::Matrix{Float64}`: lower 95% CI bound
"""
function msci95(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj2.header.recording[:data_type])), ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj2)))

    # check channels
    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    length(ch1) == length(ch2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    sm, ss, su, sl = @views msci95(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return (sm=sm, ss=ss, su=su, sl=sl)

end
