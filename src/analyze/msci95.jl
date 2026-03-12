export msci95

"""
    msci95(s; <keyword arguments>)

Calculate mean, standard error and 95% CI.

Two methods:

`:normal` – analytical CI using the normal approximation: CI = mean ± 1.96 * (std / √N)
`:boot` – empirical CI from n*length(s) bootstrap resamples: CI = [2.5th percentile, 97.5th percentile] of bootstrap means

# Arguments

- `s::AbstractVector`: signal vector
- `n::Int64=3`: bootstrap resampling factor (total resamples = `n * length(s)`)
- `method::Symbol=:normal`: `:normal` (analytical) or `:boot` (bootstrap)

# Returns

Named tuple:

- `sm::Float64`: mean
- `se::Float64`: standard error
- `su::Float64`: upper 95% CI
- `sl::Float64`: lower 95% CI
"""
function msci95(
    s::AbstractVector;
    n::Int64 = 3,
    method::Symbol = :normal
)::@NamedTuple{
    sm::Float64,
    se::Float64,
    su::Float64,
    sl::Float64
}

    _check_var(method, [:normal, :boot], "method")
    @assert n >= 1 "n must be ≥ 1."

    if method === :normal

        sm = mean(s)
        se = std(s) / sqrt(length(s))
        su = sm + 1.96 * se
        sl = sm - 1.96 * se

    else

        n_boot  = length(s) * n
        s_tmp1  = zeros(n_boot)

        @inbounds Threads.@threads :dynamic for idx1 in 1:n_boot
            s_tmp2 = zeros(length(s))
            sample_idx = rand(1:length(s), length(s))
            @inbounds for idx2 in eachindex(s)
                s_tmp2[idx2] = s[sample_idx[idx2]]
            end
            s_tmp1[idx1] = mean(s_tmp2)
        end

        sm = mean(s_tmp1)
        se = std(s_tmp1) / sqrt(n_boot)

        ssorted = sort(s_tmp1)
        sl = ssorted[round(Int, 0.025 * n_boot)]
        su = ssorted[round(Int, 0.975 * n_boot)]

    end

    return (
        sm = sm,
        se = se,
        su = su,
        sl = sl
    )

end

"""
    msci95(s; <keyword arguments>)

Calculate mean, standard error and 95% CI across rows.

Two methods:

`:normal` – analytical CI using the normal approximation: CI = mean ± 1.96 * (std / √N)
`:boot` – empirical CI from n*length(s) bootstrap resamples: CI = [2.5th percentile, 97.5th percentile] of bootstrap means

# Arguments

- `s::AbstractMatrix`: signal matrix (channels/trials × samples)
- `n::Int64=3`: bootstrap resampling factor
- `method::Symbol=:normal`: `:normal` (analytical) or `:boot` (bootstrap)

# Returns

Named tuple:

- `sm::Vector{Float64}`: column-wise mean
- `se::Vector{Float64}`: column-wise standard error
- `su::Vector{Float64}`: upper 95% CI per sample
- `sl::Vector{Float64}`: lower 95% CI per sample
"""
function msci95(
    s::AbstractMatrix;
    n::Int64 = 3,
    method::Symbol = :normal
)::@NamedTuple{
    sm::Vector{Float64},
    se::Vector{Float64},
    su::Vector{Float64},
    sl::Vector{Float64}
}

    _check_var(method, [:normal, :boot], "method")
    @assert n >= 1 "n must be ≥ 1."

    if method === :normal

        sm = dropdims(mean(s, dims = 1), dims = 1)
        se = dropdims(std(s,  dims = 1), dims = 1) ./ sqrt(size(s, 1))
        su = sm .+ 1.96 .* se
        sl = sm .- 1.96 .* se

    else

        n_boot = size(s, 1) * n
        s_tmp1 = zeros(n_boot, size(s, 2))

        @inbounds Threads.@threads :dynamic for idx1 in 1:n_boot
            s_tmp2 = zeros(size(s))
            sample_idx = rand(axes(s, 1), size(s, 1))
            @inbounds for idx2 in axes(s, 1)
                s_tmp2[idx2, :] = s[sample_idx[idx2], :]
            end
            s_tmp1[idx1, :] = mean(s_tmp2, dims = 1)
        end

        sm = dropdims(mean(s_tmp1, dims = 1), dims = 1)
        se = dropdims(std(s_tmp1,  dims = 1), dims = 1) ./ sqrt(n_boot)
        ssorted = sort(s_tmp1, dims = 1)
        sl = ssorted[round(Int, 0.025 * n_boot), :]
        su = ssorted[round(Int, 0.975 * n_boot), :]

    end

    return (
        sm = vec(sm),
        se = vec(se),
        su = vec(su),
        sl = vec(sl)
    )

end

"""
    msci95(s; <keyword arguments>)

Calculate mean, standard error and 95% CI.

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)
- `n::Int64=3`: bootstrap resampling factor
- `method::Symbol=:normal`: `:normal` (analytical) or `:boot` (bootstrap)

# Returns

Named tuple:

- `sm::Matrix{Float64}`: mean, shape `(epochs, samples)`
- `se::Matrix{Float64}`: standard error, shape `(epochs, samples)`
- `su::Matrix{Float64}`: upper 95% CI, shape `(epochs, samples)`
- `sl::Matrix{Float64}`: lower 95% CI, shape `(epochs, samples)`
"""
function msci95(
    s::AbstractArray;
    n::Int64 = 3,
    method::Symbol = :normal
)::@NamedTuple{
    sm::Matrix{Float64},
    se::Matrix{Float64},
    su::Matrix{Float64},
    sl::Matrix{Float64}
}

    _check_var(method, [:normal, :boot], "method")

    # epoch lengths
    ep_len = size(s, 2)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate outputs
    sm = zeros(ep_n, ep_len)
    se = zeros(ep_n, ep_len)
    su = zeros(ep_n, ep_len)
    sl = zeros(ep_n, ep_len)

    # calculate over epochs
    @inbounds Threads.@threads :dynamic for ep_idx in 1:ep_n
            msci_data = msci95(@view(s[:, :, ep_idx]), n = n, method = method)
            sm[ep_idx, :] = msci_data.sm
            se[ep_idx, :] = msci_data.se
            su[ep_idx, :] = msci_data.su
            sl[ep_idx, :] = msci_data.sl
    end

    return (sm = sm, se = se, su = su, sl = sl)

end

"""
    msci95(s1, s2)

Calculate mean difference, standard error and 95% CI between two vectors.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

Named tuple:

- `sm::Float64`: mean
- `se::Float64`: pooled standard error
- `su::Float64`: upper 95% CI
- `sl::Float64`: lower 95% CI
"""
function msci95(
    s1::AbstractVector,
    s2::AbstractVector
)::@NamedTuple{
    sm::Float64,
    se::Float64,
    su::Float64,
    sl::Float64
}

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    sm = mean(s1) - mean(s2)
    s1_se = std(s1) / sqrt(length(s1))
    s2_se = std(s2) / sqrt(length(s2))
    ss = sqrt(s1_se^2 + s2_se^2)
    su = sm + 1.96 * ss
    sl = sm - 1.96 * ss

    return (
        sm = sm,
        se = se,
        su = su,
        sl = sl
    )

end

"""
    msci95(s1, s2)

Calculate mean difference, standard error and 95% CI per channel and epoch.

# Arguments

- `s1::AbstractArray`: signal array (channels, samples, epochs)
- `s2::AbstractArray`: signal array (channels, samples, epochs)

# Returns

Named tuple:

- `sm::Matrix{Float64}`: mean difference, shape `(channels, epochs)`
- `se::Matrix{Float64}`: pooled SE, shape `(channels, epochs)`
- `su::Matrix{Float64}`: upper 95% CI, shape `(channels, epochs)`
- `sl::Matrix{Float64}`: lower 95% CI, shape `(channels, epochs)`
"""
function msci95(
        s1::AbstractArray, s2::AbstractArray
    )::@NamedTuple{sm::Matrix{Float64}, se::Matrix{Float64}, su::Matrix{Float64}, sl::Matrix{Float64}}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    # number of channels
    ch_n = size(s1, 1)
    # number of epochs
    ep_n = size(s1, 3)

    # pre-allocate outputs
    sm = zeros(ch_n, ep_n)
    se = zeros(ch_n, ep_n)
    su = zeros(ch_n, ep_n)
    sl = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        result = msci95(
            @view(s1[ch_idx, :, ep_idx]),
            @view(s2[ch_idx, :, ep_idx]),
        )
        sm[ch_idx, ep_idx] = result.sm
        ss[ch_idx, ep_idx] = result.ss
        su[ch_idx, ep_idx] = result.su
        sl[ch_idx, ep_idx] = result.sl
    end

    return (
        sm = sm,
        se = se,
        su = su,
        sl = sl
    )

end

"""
    msci95(obj; <keyword arguments>)

Calculate mean, standard error and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `n::Int64=3`: bootstrap resampling factor
- `method::Symbol=:normal`: `:normal` (analytical) or `:boot` (bootstrap)

# Returns

Named tuple:

- `sm::Matrix{Float64}`: mean
- `se::Matrix{Float64}`: standard error
- `su::Matrix{Float64}`: upper 95% CI
- `sl::Matrix{Float64}`: lower 95% CI
"""
function msci95(
        obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, n::Int64 = 3, method::Symbol = :normal
    )::@NamedTuple{sm::Matrix{Float64}, se::Matrix{Float64}, su::Matrix{Float64}, sl::Matrix{Float64}}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    msci_data = NeuroAnalyzer.msci95(@view(obj.data[ch, :, :]), n = n, method = method)

    return msci_data

end

"""
    msci95(obj1, obj2; <keyword arguments>)

Calculate mean difference, standard ERROR and 95% CI between two objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2:NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

Named tuple:

- `sm::Matrix{Float64}`: mean difference, shape `(channels, epochs)`
- `se::Matrix{Float64}`: pooled SE, shape `(channels, epochs)`
- `su::Matrix{Float64}`: upper 95% CI bound, shape `(channels, epochs)`
- `sl::Matrix{Float64}`: lower 95% CI bound, shape `(channels, epochs)`
"""
function msci95(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
)::@NamedTuple{
    sm::Matrix{Float64},
    se::Matrix{Float64},
    su::Matrix{Float64},
    sl::Matrix{Float64}
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    msci_data = NeuroAnalyzer.msci95(
        @view(obj1.data[ch1, :, ep1]),
        @view(obj2.data[ch2, :, ep2]),
    )

    return msci_data

end
