export detrend
export detrend!

"""
    detrend(s; type, offset, order, span, fs)

Perform piecewise detrending.

# Arguments

- `s::AbstractVector`
- `type::Symbol=:ls`:
    - `:ls`: the result of a linear least-squares fit to `s` is subtracted from `s`
    - `:linear`: linear trend is subtracted from `s`
    - `:constant`: `offset` or the mean of `s` (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` is subtracted
    - `:loess`: fit and subtract loess approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`
- `fs::Int64=0`: sampling frequency

# Returns
- `s_new::Vector{Float64}`
"""
function detrend(s::AbstractVector; type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0, fs::Int64=0)

    _check_var(type, [:ls, :linear, :constant, :poly, :loess, :hp], "type")
    f <= 0 && throw(ArgumentError("f must be > 0."))
    order < 1 && throw(ArgumentError("order must be ≥ 1."))

    if type === :loess
        t = collect(1.0:1:length(s))
        model = loess(t, s, span=f)
        trend = Loess.predict(model, t)
        s_new = s .- trend
    elseif type === :poly
        t = collect(1:1:length(s))        
        p = Polynomials.fit(t, s, order)
        trend = zeros(length(s))
        for idx in eachindex(s)
            trend[idx] = p(t[idx])
        end
        s_new = s .- trend
    elseif type === :constant
        offset == 0 && (offset = mean(s))
        s_new = s .- mean(s)
    elseif type === :ls
        T = eltype(s)
        N = size(s, 1)
        # create linear trend matrix
        A = similar(s, T, N, 2)
        A[:,2] .= T(1)
        A[:,1] .= range(T(0),T(1),length=N)
        # create linear trend matrix
        R = transpose(A) * A
        # do the matrix inverse for 2×2 matrix
        Rinv = inv(Array(R)) |> typeof(R)
        factor = Rinv * transpose(A)
        s_new = s .- A * (factor * s)
    elseif type === :linear
        trend = linspace(s[1], s[end], length(s))
        s_new = s .- trend
    elseif type === :hp
        fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
        flt = filter_create(fprototype=:fir, ftype=:hp, cutoff=f, fs=fs, n=length(s))
        s_new = filter_apply(s, flt=flt)
    end

    return s_new

end

"""
    detrend(s; type, offset, order, f)

Perform piecewise detrending.

# Arguments

- `s::AbstractArray`
- `type::Symbol=:linear`: detrending method
    - `:ls`: the result of a linear least-squares fit to signal is subtracted from signal
    - `:linear`: linear trend is subtracted from signal
    - `:constant`: `offset` or the mean of signal (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` is subtracted
    - `:loess`: fit and subtract LOESS approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for `:constant` detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`
- `fs::Int64=0`: sampling frequency

# Returns

- `s_new::Array`{Float64, 3}
"""
function detrend(s::AbstractArray; type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0, fs::Int64=0)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views detrend(s_new[ch_idx, :, ep_idx], type=type, offset=offset, order=order, f=f, fs=fs)
        end
    end

    return s_new
end

"""
    detrend(obj; ch, type, offset, order, f)

Perform piecewise detrending.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `type::Symbol=:linear`: detrending method
    - `:ls`: the result of a linear least-squares fit to signal is subtracted from signal
    - `:linear`: linear trend is subtracted from signal
    - `:constant`: `offset` or the mean of signal (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` is subtracted
    - `:loess`: fit and subtract LOESS approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for `:constant` detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function detrend(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = detrend(obj.data[ch, :, :], type=type, offset=offset, order=order, f=f, fs=sr(obj))
    reset_components!(obj_new)
    push!(obj_new.header.history, "detrend(OBJ, ch=$ch, type=$type, offset=$offset, order=$order, f=$f)")

    return obj_new

end

"""
    detrend!(obj; ch, type, offset, order, span)

Perform piecewise detrending.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `type::Symbol=:linear`: detrending method
    - `:ls`: the result of a linear least-squares fit to signal is subtracted from signal
    - `:linear`: linear trend is subtracted from signal
    - `:constant`: `offset` or the mean of signal (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` order is subtracted
    - `:loess`: fit and subtract LOESS approximation
    - `:hp`: use HP filter
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`
"""
function detrend!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)

    obj_tmp = detrend(obj, ch=ch, type=type, offset=offset, order=order, f=f)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end

