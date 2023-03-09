export detrend
export detrend!

"""
    detrend(signal; type, offset, order, span, fs)

Perform piecewise detrending.

# Arguments

- `signal::AbstractVector`
- `type::Symbol=:ls`:
    - `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:linear`: linear trend is subtracted from `signal`
    - `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` is subtracted
    - `:loess`: fit and subtract loess approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`
- `fs::Int64=0`: sampling frequency

# Returns
- `s_det::Vector{Float64}`
"""
function detrend(signal::AbstractVector; type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0, fs::Int64=0)

    _check_var(type, [:ls, :linear, :constant, :poly, :loess, :hp], "type")
    f <= 0 && throw(ArgumentError("f must be > 0."))
    order < 1 && throw(ArgumentError("order must be ≥ 1."))

    if type === :loess
        t = collect(1.0:1:length(signal))
        model = loess(t, signal, span=f)
        trend = Loess.predict(model, t)
        return signal .- trend
    elseif type === :poly
        t = collect(1:1:length(signal))        
        p = Polynomials.fit(t, signal, order)
        trend = zeros(length(signal))
        for idx in eachindex(signal)
            trend[idx] = p(t[idx])
        end
        return signal .- trend
    elseif type === :constant
        offset == 0 && (offset = mean(signal))
        return signal .- mean(signal)
    elseif type === :ls
        T = eltype(signal)
        N = size(signal, 1)
        # create linear trend matrix
        A = similar(signal, T, N, 2)
        A[:,2] .= T(1)
        A[:,1] .= range(T(0),T(1),length=N)
        # create linear trend matrix
        R = transpose(A) * A
        # do the matrix inverse for 2×2 matrix
        Rinv = inv(Array(R)) |> typeof(R)
        factor = Rinv * transpose(A)
        return signal .- A * (factor * signal)
    elseif type === :linear
        trend = linspace(signal[1], signal[end], length(signal))
        return signal .- trend
    elseif type === :hp
        fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
        flt = filter_create(fprototype=:fir, ftype=:hp, cutoff=f, fs=fs, n=length(signal))
        return filter_apply(signal, flt=flt)
    end
end

"""
    detrend(obj; channel, type, offset, order, f)

Perform piecewise detrending.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `type::Symbol=:linear`: detrending method
    - `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:linear`: linear trend is subtracted from `signal`
    - `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` is subtracted
    - `:loess`: fit and subtract LOESS approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for `:constant` detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function detrend(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)

    _check_var(type, [:ls, :linear, :constant, :poly, :loess, :hp], "type")

    ep_n = epoch_n(obj)
    fs = sr(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
            @views obj_new.data[channel[ch_idx], :, ep_idx] = detrend(obj_new.data[channel[ch_idx], :, ep_idx], type=type, offset=offset, order=order, f=f, fs=fs)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "detrend(OBJ, channel=$channel, type=$type, offset=$offset, order=$order, f=$f)")

    return obj_new
end

"""
    detrend!(obj; channel, type, offset, order, span)

Perform piecewise detrending.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `type::Symbol=:linear`: detrending method
    - `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:linear`: linear trend is subtracted from `signal`
    - `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` order is subtracted
    - `:loess`: fit and subtract LOESS approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`
"""
function detrend!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)

    obj_tmp = detrend(obj, channel=channel, type=type, offset=offset, order=order, f=f)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing
end

