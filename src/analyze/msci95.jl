export msci95

"""
    msci95(s; <keyword arguments>)

Calculate mean, standard deviation and 95% confidence interval.

# Arguments

- `s::AbstractVector`
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal method (`:normal`) or `n`-times boostrapping (`:boot`)

# Returns

Named tuple containing:
- `sm::Float64`: mean
- `ss::Float64`: standard deviation
- `su::Float64`: upper 95% CI
- `sl::Float64`: lower 95% CI
"""
function msci95(s::AbstractVector; n::Int64=3, method::Symbol=:normal)::@NamedTuple{sm::Float64, ss::Float64, su::Float64, sl::Float64}

    _check_var(method, [:normal, :boot], "method")
    @assert n >= 1 "n must be ≥ 1."

    if method === :normal
        sm = mean(s)
        ss = std(s) / sqrt(length(s))
        su = sm + 1.96 * ss
        sl = sm - 1.96 * ss
    else
        s_tmp1 = zeros(length(s) * n)
        Threads.@threads for idx1 in eachindex(s) * n
            s_tmp2 = zeros(length(s))
            sample_idx = rand(1:length(s), length(s))
            @inbounds for idx2 in eachindex(s)
                s_tmp2[idx2] = s[sample_idx[idx2]]
            end
            s_tmp1[idx1] = mean(s_tmp2)
        end

        sm = mean(s_tmp1)'
        ss = std(s_tmp1)' / sqrt(length(s_tmp1))
        ssorted = sort(s_tmp1)
        sl = ssorted[round(Int, 0.025 * length(s_tmp1)), :]
        su = ssorted[round(Int, 0.975 * length(s_tmp1)), :]
    end

    return (sm=sm, ss=ss, su=su, sl=sl)

end

"""
    msci95(s; <keyword arguments>)

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
function msci95(s::AbstractMatrix; n::Int64=3, method::Symbol=:normal)::@NamedTuple{sm::Vector{Float64}, ss::Vector{Float64}, su::Vector{Float64}, sl::Vector{Float64}}

    _check_var(method, [:normal, :boot], "method")
    @assert n >= 1 "n must be ≥ 1."

    if method === :normal
        sm = mean(s, dims=1)'
        ss = std(s, dims=1)' / sqrt(size(s, 1))
        su = sm + 1.96 * ss
        sl = sm - 1.96 * ss
    else
        s_tmp1 = zeros(size(s, 1) * n, size(s, 2))
        Threads.@threads for idx1 in axes(s, 1) * n
            s_tmp2 = zeros(size(s))
            sample_idx = rand(axes(s, 1), size(s, 1))
            @inbounds for idx2 in axes(s, 1)
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
    msci95(s; <keyword arguments>)

Calculate mean, standard deviation and 95% confidence interval.

# Arguments

- `s::AbstractArray`
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal method (`:normal`) or `n`-times boostrapping (`:boot`)

# Returns

Named tuple containing:
- `sm::Matrix{Float64}`: mean
- `ss::Matrix{Float64}`: standard deviation
- `su::Matrix{Float64}`: upper 95% CI
- `sl::Matrix{Float64}`: lower 95% CI
"""
function msci95(s::AbstractArray; n::Int64=3, method::Symbol=:normal)::@NamedTuple{sm::Matrix{Float64}, ss::Matrix{Float64}, su::Matrix{Float64}, sl::Matrix{Float64}}

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
function msci95(s1::AbstractVector, s2::AbstractVector)::@NamedTuple{sm::Float64, ss::Float64, su::Float64, sl::Float64}

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

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
- `sm::Matrix{Float64}`: mean
- `ss::Matrix{Float64}`: standard deviation
- `su::Matrix{Float64}`: upper 95% CI
- `sl::Matrix{Float64}`: lower 95% CI
"""
function msci95(s1::AbstractArray, s2::AbstractArray)::@NamedTuple{sm::Matrix{Float64}, ss::Matrix{Float64}, su::Matrix{Float64}, sl::Matrix{Float64}}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    sm = zeros(ch_n, ep_n)
    ss = zeros(ch_n, ep_n)
    su = zeros(ch_n, ep_n)
    sl = zeros(ch_n, ep_n)

    Threads.@threads for ep_idx in 1:ep_n
        @inbounds for ch_idx in 1:ch_n
            sm[ch_idx, ep_idx], ss[ch_idx, ep_idx], su[ch_idx, ep_idx], sl[ch_idx, ep_idx] = @views msci95(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx])
        end
    end

    return (sm=sm, ss=ss, su=su, sl=sl)

end

"""
    msci95(obj; <keyword arguments>)

Calculate mean, standard deviation and 95% confidence interval.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal (`:normal`) method or `n`-times bootstrapping (`:boot`)

# Returns

Named tuple containing:
- `sm::Matrix{Float64}`: mean
- `ss::Matrix{Float64}`: standard deviation
- `su::Matrix{Float64}`: upper 95% CI
- `sl::Matrix{Float64}`: lower 95% CI
"""
function msci95(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, n::Int64=3, method::Symbol=:normal)::@NamedTuple{sm::Matrix{Float64}, ss::Matrix{Float64}, su::Matrix{Float64}, sl::Matrix{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    sm, ss, su, sl = @views NeuroAnalyzer.msci95(obj.data[ch, :, :], n=n, method=method)

    return (sm=sm, ss=ss, su=su, sl=sl)

end

"""
    msci95(obj1, obj2; <keyword arguments>)

Calculate mean difference, standard deviation and 95% confidence interval.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2:NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}`: list of channels
- `ch2::Union{String, Vector{String}}`: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `sm::Matrix{Float64}`: mean
- `ss::Matrix{Float64}`: standard deviation
- `su::Matrix{Float64}`: upper 95% CI bound
- `sl::Matrix{Float64}`: lower 95% CI bound
"""
function msci95(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))::@NamedTuple{sm::Matrix{Float64}, ss::Matrix{Float64}, su::Matrix{Float64}, sl::Matrix{Float64}}

    # check channels
    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal."

    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1)) and ep2 ($(length(ep2)) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    sm, ss, su, sl = @views NeuroAnalyzer.msci95(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return (sm=sm, ss=ss, su=su, sl=sl)

end
