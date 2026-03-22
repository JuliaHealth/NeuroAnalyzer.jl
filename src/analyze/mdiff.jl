export mdiff

"""
    mdiff(s1, s2; <keyword arguments>)

Calculate the mean difference and its bootstrap p-value for two signal matrices. Estimate whether the mean difference between two signal matrices is statistically significant using a permutation/bootstrap approach:

1. Compute the observed statistic (max absolute difference or integrated squared difference) between the mean waveforms of s1 and s2.
2. Pool s1 and s2 rows; repeatedly resample with replacement and compute the same statistic on random pairs of sub-groups.
3. p = proportion of bootstrap statistics exceeding the observed statistic.

# Arguments

- `s1::AbstractMatrix`: signal matrix (channels, samples)
- `s2::AbstractMatrix`: signal matrix (channels, samples)
- `n::Int64=3`: number of bootstrap iterations per channel
- `method::Symbol=:absdiff`: test statistic:
    - `:absdiff`: maximum absolute difference between mean waveforms
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple:

- `st::Vector{Float64}`: bootstrap distribution of the test statistic
- `sts::Float64`: observed test statistic
- `p::Float64`: proportion of bootstrap statistics exceeding `sts`
"""
function mdiff(
    s1::AbstractMatrix,
    s2::AbstractMatrix;
    n::Int64 = 3,
    method::Symbol = :absdiff
)::@NamedTuple{
    st::Vector{Float64},
    sts::Float64,
    p::Float64
}

    _check_var(method, [:absdiff, :diff2int], "method")
    !(size(s1) == size(s2)) && throw(ArgumentError("s1 and s2 must have the same size."))
    !(n >= 1) && throw(ArgumentError("n must be ≥ 1."))

    # calculate means over channels
    s1_mean = vec(mean(s1, dims = 1))
    s2_mean = vec(mean(s2, dims = 1))

    # observed test statistic.
    if method === :absdiff
        # statistic: maximum absolute difference between mean waveforms
        sts = maximum(abs, s1_mean .- s2_mean)
    else
        # statistic: integrated area of the squared difference
        sts = simpson((s1_mean .- s2_mean) .^ 2)
    end

    # pooled sample for bootstrapping
    ss = [s1; s2]
    n_boot = size(s1, 1) * n
    st = zeros(n_boot)

    # pre-allocate bootstrap buffers outside the thread loop 
    # each thread needs its own buffers to avoid races
    # allocate inside the loop per iteration since Threads.@threads
    # does not provide thread-local storage here
    Threads.@threads :dynamic for idx in 1:n_boot

        # sample two independent bootstrap groups from the pooled data
        s_tmp1 = zeros(size(s1, 1), size(s1, 2))
        idx_a  = rand(axes(ss, 1), size(s1, 1))
        @inbounds for idx2 in axes(s1, 1)
            s_tmp1[idx2, :] = ss[idx_a[idx2], :]
        end
        bs1_mean = vec(mean(s_tmp1, dims = 1))

        s_tmp2 = zeros(size(s1, 1), size(s1, 2))
        idx_b  = rand(axes(ss, 1), size(s1, 1))
        @inbounds for idx2 in axes(s1, 1)
            s_tmp2[idx2, :] = ss[idx_b[idx2], :]
        end
        bs2_mean = vec(mean(s_tmp2, dims = 1))

        if method === :absdiff
            # statistic: maximum absolute difference between mean waveforms
            @inbounds st[idx] = maximum(abs, bs1_mean .- bs2_mean)
        else
            # statistic: integrated area of the squared difference
            @inbounds st[idx] = simpson((bs1_mean .- bs2_mean) .^ 2)
        end

    end

    p = count(x -> x > sts, st) / n_boot
    p > 1.0 && (p = 1.0)

    return (st = st, sts = sts, p = p)

end

"""
    mdiff(s1, s2; <keyword arguments>)

Calculate the mean difference and its bootstrap p-value for each epoch.

# Arguments

- `s1::AbstractArray`: signal array (channels, samples, epochs)
- `s2::AbstractArray`: signal array (channels, samples, epochs)
- `n::Int64=3`: number of bootstrap iterations per channel
- `method::Symbol=:absdiff`: test statistic:
    - `:absdiff`: maximum absolute difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple:

- `st::Matrix{Float64}`: bootstrap distributions, shape `(epochs, channels*n)`
- `sts::Vector{Float64}`: observed statistics per epoch
- `p::Vector{Float64}`: p-values per epoch
"""
function mdiff(
    s1::AbstractArray,
    s2::AbstractArray;
    n::Int64 = 3,
    method::Symbol = :absdiff
)::@NamedTuple{
    st::Matrix{Float64},
    sts::Vector{Float64},
    p::Vector{Float64}
}

    !(size(s1) == size(s2)) && throw(ArgumentError("s1 and s2 must have the same size."))
    !(n >= 1) && throw(ArgumentError("n must be ≥ 1."))

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s1)
    _chk3d(s2)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate outputs
    st = zeros(ep_n, ch_n * n)
    sts = zeros(ep_n)
    p = zeros(ep_n)

    # calculate over epochs
    @inbounds @Threads.threads :dynamic for ep_idx in 1:ep_n
        mdriff_data = mdiff(
            @view(s1[:, :, ep_idx]),
            @view(s2[:, :, ep_idx]),
            n = n,
            method = method,
        )
        st[ep_idx, :] = mdriff_data.st
        sts[ep_idx] = mdriff_data.sts
        p[ep_idx] = mdriff_data.p
    end

    return (st = st, sts = sts, p = p)
end

"""
    mdiff(obj1, obj2; <keyword arguments>)

Calculate the mean difference and its bootstrap p-value for two objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2:NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)
- `n::Int64`: number of bootstraps
- `method::Symbol=:absdiff`: test statistic:
    - `:absdiff`: maximum absolute difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple:

- `st::Matrix{Float64}`: bootstrap distributions, shape `(epochs, channels*n)`
- `sts::Vector{Float64}`: observed statistics per epoch
- `p::Vector{Float64}`: p-values per epoch
"""
function mdiff(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
    n::Int64 = 3,
    method::Symbol = :absdiff
)::@NamedTuple{
    st::Matrix{Float64},
    sts::Vector{Float64},
    p::Vector{Float64}
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    (length(ch1) == length(ch2)) || throw(ArgumentError("Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."))

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    (length(ep1) == length(ep2)) || throw(ArgumentError("Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."))
    (epoch_len(obj1) == epoch_len(obj2)) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    mdiff_data = NeuroAnalyzer.mdiff(
        @view(obj1.data[ch1, :, ep1]),
        @view(obj2.data[ch2, :, ep2]),
        n = n,
        method = method,
    )

    return mdiff_data

end
