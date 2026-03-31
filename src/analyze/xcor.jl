export xcor

"""
    xcor(s1, s2; <keyword arguments>)

Calculate cross-correlation between two 1-D signal vectors.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector (must be the same length as `s1`)
- `l::Int64=round(Int64, min(length(s1) - 1, 10 * log10(length(s1))))`: maximum lag in samples; lags range is `−l : l`
- `demean::Bool=true`: subtract the mean before computing cross-correlation
- `biased::Bool=true`: use biased (÷ n) or unbiased (÷ n−lag) estimator
- `method::Symbol=:sum`: computation method:
    - `:sum`: manual lag-shifted dot product
    - `:cov`: use Julia's `cor()` (`biased` is ignored)
    - `:stat`: use `StatsBase.crosscor()` (`biased` is ignored)

# Returns

- `Array{Float64, 3}`: cross-correlation at lags `−l:l`
"""
function xcor(
    s1::AbstractVector,
    s2::AbstractVector;
    l::Int64 = round(Int64, min(length(s1) - 1, 10 * log10(length(s1)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::Array{Float64, 3}

    # validate
    _check_var(method, [:sum, :cov, :stat], "method")
    length(s1) == length(s2) ||
        throw(ArgumentError("Both signals must have the same length."))

    # optionally remove the DC component (mean) from both signals before
    # computing cross-correlation to eliminate offset bias
    if demean
        s1_tmp = remove_dc(s1)
        s2_tmp = remove_dc(s2)
    else
        s1_tmp = s1
        s2_tmp = s2
    end

    # pre-allocate outputs
    xc = zeros(l + 1)
    xc_neg = zeros(l + 1)

    if method === :sum
        # ---- positive lags: s1 leads s2 by idx samples ---------------
        for idx in 0:l
            # shift s1 forward by idx: align s1[1+idx:end] with s2[1:end-idx]
            xc[idx + 1] = @views sum(s1_tmp[(1 + idx):end] .* s2_tmp[1:(end - idx)])
            # normalize: biased divides by n; unbiased by (n − lag) to
            # correct for the reduced number of overlapping samples
            xc[idx + 1] /= biased ? length(s1) : (length(s1) - idx)
        end
        # ---- negative lags: s2 leads s1 by idx samples ---------------
        for idx in 0:l
            xc_neg[idx + 1] = @views sum(s1_tmp[1:(end - idx)] .* s2_tmp[(1 + idx):end])
            xc_neg[idx + 1] /= biased ? length(s1) : (length(s1) - idx)
        end
    elseif method === :cor
        for idx in 0:l
            xc[idx + 1] = @views cor(s1_tmp[(1 + idx):end], s2_tmp[1:(end - idx)])
        end
        for idx in 0:l
            xc_neg[idx + 1] = @views cor(s1_tmp[1:(end - idx)], s2_tmp[(1 + idx):end])
        end
    elseif method === :stat
        # StatsBase crosscor handles demeaning internally; `biased` is ignored.
        xc = crosscor(s1, s2, 0:l; demean = demean)
        xc_neg = crosscor(s2, s1, 0:l; demean = demean)
    end

    # concatenate negative lags (reversed) with positive lags to produce a
    # symmetric lag vector from −l to +l. xc_neg[1] is lag 0 (same as xc[1])
    # so drop the duplicate when concatenating
    xc = vcat(reverse(xc_neg), xc[2:end])
    if method === :sum
        xc = xc ./ (std(s1) * std(s2))
    end
    xc = round.(xc; digits = 3)

    return reshape(xc, 1, :, 1)
end

"""
    xcor(s1, s2; <keyword arguments>)

Calculate cross-correlation for a pair of 2-D arrays.

# Arguments

- `s1::AbstractMatrix`: signal matrix (channels, epochs)
- `s2::AbstractMatrix`: signal matrix, same size as `s1`
- `l::Int64=round(Int64, min(size(s1[1, :, 1], 1) - 1, 10 * log10(size(s1[1, :, 1], 1))))`: maximum lag in samples; lags range is `−l : l`
- `demean::Bool=true`: subtract the mean before computing cross-correlation
- `biased::Bool=true`: use biased (÷ n) or unbiased (÷ n−lag) estimator
- `method::Symbol=:sum`: computation method:
    - `:sum`: manual lag-shifted dot product
    - `:cov`: use Julia's `cor()` (`biased` is ignored)
    - `:stat`: use `StatsBase.crosscor()` (`biased` is ignored)

# Returns

- `Array{Float64, 3}`
"""
function xcor(
    s1::AbstractMatrix,
    s2::AbstractMatrix;
    l::Int64 = round(Int64, min(size(s1, 1), 10 * log10(size(s1, 1)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::Array{Float64, 3}

    # validate
    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))

    # number of epochs
    ep_n = size(s1, 2)

    # pre-allocate output
    xc = zeros(1, length((-l):l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        xc[1, :, ep_idx] = @views xcor(
            s1[:, ep_idx],
            s2[:, ep_idx];
            l = l,
            demean = demean,
            biased = biased,
            method = method,
        )
    end

    return xc
end

"""
    xcor(s1, s2; <keyword arguments>)

Calculate cross-correlation for two 3-D signal arrays.

# Arguments

- `s1::AbstractArray`: signal array, shape (channels, samples, epochs)
- `s2::AbstractArray`: signal array, shape (channels, samples, epochs)
- `l::Int64=round(Int64, min(size(s1[1, :, 1], 1) - 1, 10 * log10(size(s1[1, :, 1], 1))))`: maximum lag in samples; lags range is `−l : l`
- `demean::Bool=true`: subtract the mean before computing cross-correlation
- `biased::Bool=true`: use biased (÷ n) or unbiased (÷ n−lag) estimator
- `method::Symbol=:sum`: computation method:
    - `:sum`: manual lag-shifted dot product
    - `:cov`: use Julia's `cor()` (`biased` is ignored)
    - `:stat`: use `StatsBase.crosscor()` (`biased` is ignored)

# Returns

- `Array{Float64, 3}`
"""
function xcor(
    s1::AbstractArray,
    s2::AbstractArray;
    l::Int64 = round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s1)
    _chk3d(s2)
    # validate
    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))

    # number of channels
    ch_n = size(s1, 1)
    # number of epochs
    ep_n = size(s1, 3)

    # pre-allocate output
    xc = zeros(ch_n, length((-l):l), ep_n)

    # calculate over channels and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        xc[ch_idx, :, ep_idx] = @views xcor(
            s1[ch_idx, :, ep_idx],
            s2[ch_idx, :, ep_idx],
            l = l,
            demean = demean,
            biased = biased,
            method = method,
        )
    end

    return xc
end

"""
    xcor(obj1, obj2; <keyword arguments>)

Calculate cross-correlation between selected channels of two NEURO objects.

For ERP/ERF objects the trial-averaged cross-correlation is prepended as epoch 1.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s) in `obj1`
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s) in `obj2`
- `ep1::Union{Int64, Vector{Int64}, UnitRange{Int64}}=_c(nepochs(obj1))`: epoch number(s) in `obj1`
- `ep2::Union{Int64, Vector{Int64}, UnitRange{Int64}}=_c(nepochs(obj2))`: epoch number(s) in `obj2`
- `l::Real=1`: maximum lag in samples; lags range is `−l : l`
- `demean::Bool=true`: subtract the mean before computing cross-correlation
- `biased::Bool=true`: use biased (÷ n) or unbiased (÷ n−lag) estimator
- `method::Symbol=:sum`: computation method:
    - `:sum`: manual lag-shifted dot product
    - `:cov`: use Julia's `cor()` (`biased` is ignored)
    - `:stat`: use `StatsBase.crosscor()` (`biased` is ignored)

# Returns

Named tuple:

- `xc::Array{Float64, 3}`: cross-correlation
- `lags::Vector{Float64}`: lag values in seconds
"""
function xcor(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, UnitRange{Int64}} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, UnitRange{Int64}} = _c(nepochs(obj2)),
    l::Real = 1,
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::@NamedTuple{
    xc::Array{Float64, 3},
    lags::Vector{Float64},
}

    # validate
    sr(obj1) == sr(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same sampling rate."))
    length(ch1) == length(ch2) ||
        throw(
            ArgumentError(
                "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal.",
            ),
        )
    length(ep1) == length(ep2) ||
        throw(
            ArgumentError(
                "Lengths of ep1 ($(length(ep1)) and ep2 ($(length(ep2)) must be equal.",
            ),
        )
    epoch_len(obj1) == epoch_len(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 =
        exclude_bads ? get_channel(obj1; ch = ch1, exclude = "bad") :
        get_channel(obj1; ch = ch1, exclude = "")
    ch2 =
        exclude_bads ? get_channel(obj2; ch = ch2, exclude = "bad") :
        get_channel(obj2; ch = ch2, exclude = "")
    isempty(ch1) && throw(ArgumentError("No channels selected."))
    isempty(ch2) && throw(ArgumentError("No channels selected."))
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    # validate lag bound against the epoch length (both in seconds)
    max_l = epoch_len(obj1) / sr(obj1)
    (0 <= l <= max_l) ||
        throw(ArgumentError("l must be in [0, $max_l] seconds."))

    l_samples = round(Int64, l * sr(obj1))

    if datatype(obj1) == "erp" && datatype(obj2) == "erp"
        xc = @views xcor(
            reshape(obj1.data[ch1, :, 2:end], length(ch1), :, (nepochs(obj1) - 1)),
            reshape(obj2.data[ch2, :, 2:end], length(ch2), :, (nepochs(obj2) - 1)),
            l = l,
            demean = demean,
            biased = biased,
            method = method,
        )
        xc = cat(mean(xc; dims = 3), xc; dims = 3)
    else
        xc = @views xcor(
            obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2], l = l, demean = demean,
            biased = biased, method = method,
        )
    end

    # convert lag indices back to seconds for the returned lag axis
    lags = collect((-l_samples):l_samples) ./ sr(obj1)

    return (; xc, lags)
end
