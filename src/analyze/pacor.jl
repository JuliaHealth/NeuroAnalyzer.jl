export pacor

"""
   pacor(s; l, demean, method)

Calculate partial auto-correlation.

# Arguments

- `s::AbstractVector`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `method::Symbol=:yw`: method of calculating auto-correlation:
    - `:yw`: computes the partial autocorrelations using the Yule-Walker equations
    - `:reg`: computes the partial autocorrelations via successive regression models

# Returns

- `pac::Matrix{Float64}`

# Notes

If you get `ERROR: PosDefException: matrix is not positive definite; Cholesky factorization failed.`, try lowering `l` value or change method to `:yw`.
"""
function pacor(s::AbstractVector; l::Int64=round(Int64, min(length(s) - 1, 10 * log10(length(s)))), demean::Bool=true, method::Symbol=:yw)

    _check_var(method, [:reg, :yw], "method")

    method === :reg && (method = :regression)
    method === :yw && (method = :yulewalker)

    pac = zeros(l + 1)
    pac_neg = zeros(l + 1)

    if demean
        s_tmp = delmean(s)
    else
        s_tmp = s
    end

    pac = pacf(s_tmp, collect(0:l), method=method)
    pac_neg = pacf(reverse(s_tmp), collect(0:l), method=method)

    pac = vcat(reverse(pac_neg), pac[2:end])
    pac = round.(pac, digits=3)

    return reshape(pac, 1, :, 1)

end

"""
   pacor(s; l, demean, method)

Calculate partial auto-correlation.

# Arguments

- `s::AbstractMatrix`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `method::Symbol=:yw`: method of calculating auto-correlation:
    - `:yw`: computes the partial autocorrelations using the Yule-Walker equations
    - `:reg`: computes the partial autocorrelations via successive regression models

# Returns

- `pac::Matrix{Float64}`
"""
function pacor(s::AbstractMatrix; l::Int64=round(Int64, min(size(s[:, 1], 1) - 1, 10 * log10(size(s[:, 1], 1)))), demean::Bool=true, method::Symbol=:yw)

    ep_n = size(s, 2)

    pac = zeros(1, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        pac[1, :, ep_idx] = @views reshape(pacor(s[:, ep_idx], l=l, demean=demean, method=method), 1, :, ep_n)
    end

    return pac

end

"""
   pacor(s; l, demean, method)

Calculate partial auto-correlation.

# Arguments

- `s::AbstractArray`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `method::Symbol=:yw`: method of calculating auto-correlation:
    - `:yw`: computes the partial autocorrelations using the Yule-Walker equations
    - `:reg`: computes the partial autocorrelations via successive regression models

# Returns

- `pac::Matrix{Float64}`
"""
function pacor(s::AbstractArray; l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1)))), demean::Bool=true, method::Symbol=:yw)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    pac = zeros(ch_n, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pac[ch_idx, :, ep_idx] = @views pacor(s[ch_idx, :, ep_idx], l=l, demean=demean, method=method)
        end
    end

    return pac

end

"""
   pacor(obj; ch, lag, demean, method)

Calculate partial auto-correlation. For ERP return trial-averaged auto-correlation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `l::Real=1`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `method::Symbol=:yw`: method of calculating auto-correlation:
    - `:yw`: computes the partial autocorrelations using the Yule-Walker equations
    - `:reg`: computes the partial autocorrelations via successive regression models

# Returns

Named tuple containing:
- `pac::Array{Float64, 3}`
- `l::Vector{Float64}`: lags [s]
"""
function pacor(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), l::Real=1, demean::Bool=true, method::Symbol=:yw)

    @assert (l > 1 && method === :yw) "For :yw method, l must be > 1."

    _check_channels(obj, ch)
    @assert l <= size(obj, 2) "l must be ≤ $(size(obj, 2))."
    @assert l >= 0 "l must be ≥ 0."

    if datatype(obj) == "erp"
        pac = @views pacor(obj.data[ch, :, 2:end], l=l, demean=demean, method=method)
        pac = cat(mean(pac, dims=3), pac, dims=3)
    else
        pac = @views pacor(obj.data[ch, :, :], l=l, demean=demean, method=method)
    end

    return (pac=pac, l=collect(-l:l) .* 1/sr(obj))

end
