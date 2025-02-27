export xcor

"""
    xcor(s1, s2; <keyword arguments>)

Calculate cross-correlation.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `l::Int64=round(Int64, min(length(s1) - 1, 10 * log10(length(s1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing cross-correlation
- `biased::Bool=true`: calculate biased or unbiased cross-correlation
- `method::Symbol=:sum`: method of calculating cross-correlation:
    - `:sum`: `acf = Σ(s1[1:end - l] .* s1[1+l:end]) ./ (std(s1) × std(s2))`
    - `:cor`: `acf = cor(s1[1:end - l], s2[1+l:end])`, `biased` value is ignored
    - `:stat`: use StatsBase `crosscor()`, `biased` value is ignored

# Returns

- `xc::Array{Float64, 3}`
"""
function xcor(s1::AbstractVector, s2::AbstractVector; l::Int64=round(Int64, min(length(s1) - 1, 10 * log10(length(s1)))), demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::Array{Float64, 3}

    _check_var(method, [:sum, :cor, :stat], "method")

    @assert length(s1) == length(s2) "Both signals must have the same length."

    xc = zeros(l + 1)
    xc_neg = zeros(l + 1)

    if demean
        s1_tmp = delmean(s1)
        s2_tmp = delmean(s2)
    else
        s1_tmp = s1
        s2_tmp = s2
    end

    if method === :sum
        for idx in 0:l
            xc[idx + 1] = @views sum(s1_tmp[(1 + idx):end] .* s2_tmp[1:(end - idx)])
            if biased
                xc[idx + 1] /= length(s1)
            else
                xc[idx + 1] /= (length(s1) - idx)
            end
        end
        for idx in 0:l
            xc_neg[idx + 1] = @views sum(s1_tmp[1:(end - idx)] .* s2_tmp[(1 + idx):end])
            if biased
                xc_neg[idx + 1] /= length(s1)
            else
                xc_neg[idx + 1] /= (length(s1) - idx)
            end
        end
    elseif method === :cor
        for idx in 0:l
            xc[idx + 1] = @views cor(s1_tmp[(1 + idx):end], s2_tmp[1:(end - idx)])
        end
        for idx in 0:l
            xc_neg[idx + 1] = @views cor(s1_tmp[1:(end - idx)], s2_tmp[(1 + idx):end])
        end
    elseif method === :stat
        xc = crosscor(s1, s2, 0:l, demean=demean)
        xc_neg = crosscor(s2, s1, 0:l, demean=demean)
    end

    xc = vcat(reverse(xc_neg), xc[2:end])
    if method === :sum
        xc = xc ./ (std(s1) * std(s2))
    end
    xc = round.(xc, digits=3)

    return reshape(xc, 1, :, 1)

end

"""
    xcor(s1, s2; <keyword arguments>)

Calculate cross-correlation.

# Arguments

- `s1::AbstractMatrix`
- `s2::AbstractMatrix`
- `l::Int64=round(Int64, min(size(s1[1, :, 1], 1) - 1, 10 * log10(size(s1[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing cross-correlation
- `biased::Bool=true`: calculate biased or unbiased cross-correlation
- `method::Symbol=:sum`: method of calculating cross-correlation:
    - `:sum`: `acf = Σ(s1[1:end - l] .* s1[1+l:end]) ./ var(s)`
    - `:cor`: `acf = cor(s1[1:end - l], s2[1+l:end])`
    - `:stat`: use StatsBase `crosscor()`, `biased` value is ignored

# Returns

- `xc::Array{Float64, 3}`
"""
function xcor(s1::AbstractMatrix, s2::AbstractMatrix; l::Int64=round(Int64, min(size(s1, 1), 10 * log10(size(s1, 1)))), demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::Array{Float64, 3}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ep_n = size(s1, 2)
    xc = zeros(1, length(-l:l), ep_n)
    @inbounds for ep_idx in 1:ep_n
        xc[1, :, ep_idx] = @views xcor(s1[1, ep_idx], s2[1, ep_idx], l=l, demean=demean, biased=biased, method=method)
    end

    return xc

end

"""
    xcor(s1, s2; <keyword arguments>)

Calculate cross-correlation.

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `l::Int64=round(Int64, min(size(s1[1, :, 1], 1) - 1, 10 * log10(size(s1[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing cross-correlation
- `biased::Bool=true`: calculate biased or unbiased cross-correlation
- `method::Symbol=:sum`: method of calculating cross-correlation:
    - `:sum`: `acf = Σ(s1[1:end - l] .* s1[1+l:end]) ./ var(s)`
    - `:cor`: `acf = cor(s1[1:end - l], s2[1+l:end])`
    - `:stat`: use StatsBase `crosscor()`, `biased` value is ignored

# Returns

- `xc::Array{Float64, 3}`
"""
function xcor(s1::AbstractArray, s2::AbstractArray; l::Int64=round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2)))), demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::Array{Float64, 3}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _chk3d(s1)
    _chk3d(s2)

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    xc = zeros(ch_n, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :static for ch_idx in 1:ch_n
            xc[ch_idx, :, ep_idx] = @views xcor(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], l=l, demean=demean, biased=biased, method=method)
        end
    end

    return xc

end

"""
    xcor(obj1, obj2; <keyword arguments>)

Calculate cross-correlation. For ERP return trial-averaged cross-correlation.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}`: list of channels
- `ch2::Union{String, Vector{String}}`: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `l::Real=1`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing cross-correlation
- `biased::Bool=true`: calculate biased or unbiased cross-correlation
- `method::Symbol=:sum`: method of calculating cross-correlation:
    - `:sum`: `acf = Σ(s1[1:end - l] .* s1[1+l:end]) ./ var(s)`
    - `:cor`: `acf = cor(s1[1:end - l], s2[1+l:end])`
    - `:stat`: use StatsBase `crosscor()`, `biased` value is ignored

# Returns

Named tuple containing:
- `xc::Array{Float64, 3}`: cross-correlation
- `l::Vector{Float64}`: lags [s]
"""
function xcor(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)), l::Real=1, demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::@NamedTuple{xc::Array{Float64, 3}, l::Vector{Float64}}

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    ch1 = get_channel(obj1, ch=ch1)
    ch2 = get_channel(obj2, ch=ch2)
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    @assert l <= size(obj1, 2) "l must be ≤ $(size(obj1, 2))."
    @assert l >= 0 "l must be ≥ 0."

    if datatype(obj1) == "erp" && datatype(obj2) == "erp"
        xc = @views xcor(reshape(obj1.data[ch1, :, 2:end], length(ch1), :, (nepochs(obj1) - 1)), reshape(obj2.data[ch2, :, 2:end], length(ch2), :, (nepochs(obj2) - 1)), l=l, demean=demean, biased=biased, method=method)
        xc = cat(mean(xc, dims=3), xc, dims=3)
    else
        xc = @views xcor(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2], l=l, demean=demean, biased=biased, method=method)
    end

    return (xc=xc, l=collect(-l:l) .* 1/sr(obj1))

end
