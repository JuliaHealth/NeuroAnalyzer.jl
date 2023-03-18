export acov

"""
   acov(s; lag, norm)

Calculate auto-covariance.

# Arguments

- `s::AbstractVector`
- `lag::Int64=1`: lags range is `-lag:lag`
- `norm::Bool=false`: normalize auto-covariance

# Returns

Named tuple containing:
- `ac::Vector{Float64}`: auto-covariance
- `l::Vector{Int64}`: lags
"""
function acov(s::AbstractVector; lag::Int64=1, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))
    l = collect(-lag:lag)

    ss2 = sum(s.^2)

    ac = zeros(length(l))

    @inbounds @simd for idx in eachindex(l)
        # no lag
        l[idx] == 0 && (ac[idx] = ss2)
        # positive lag
        l[idx] > 0 && (ac[idx] = @views sum(s[(1 + l[idx]):end] .* s[1:(end - l[idx])]))
        # negative lag
        l[idx] < 0 && (ac[idx] = @views sum(s[1:(end - abs(l[idx]))] .* s[(1 + abs(l[idx])):end]))
    end
    norm == true && (ac ./ length(s))

    return (ac=ac, l=l)

end

"""
   acov(s; lag, norm)

Calculate auto-covariance.

# Arguments

- `s::AbstractArray`
- `lag::Int64=1`: lags range is `-lag:lag`
- `norm::Bool=false`: normalize auto-covariance

# Returns

Named tuple containing:
- `ac::Matrix{Float64}`
- `l::Vector{Float64}`
"""
function acov(s::AbstractArray; lag::Int64=1, norm::Bool=false)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    l = collect(-lag:lag)

    ac = zeros(ch_n, length(-lag:lag), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ac[ch_idx, :, ep_idx], _ = @views acov(s[ch_idx, :, ep_idx], lag=lag, norm=norm)
        end
    end

    return (ac=ac, l=l)

end

"""
   acov(obj; ch, lag, norm)

Calculate auto-covariance.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `lag::Int64=1`: lags range is `-lag:lag`
- `norm::Bool=false`: normalize auto-covariance

# Returns

Named tuple containing:
- `ac::Matrix{Float64}`
- `l::Vector{Float64}`: lags in ms
"""
function acov(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), lag::Int64=1, norm::Bool=false)

    _check_channels(obj, ch)

    ac, _ = @views acov(obj.data[ch, :, :], lag=lag, norm=norm)
    l = 1/sr(obj) .* collect(-lag:lag) .* 1000

    return (ac=ac, l=l)

end
