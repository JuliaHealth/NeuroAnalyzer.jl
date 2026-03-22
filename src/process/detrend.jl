export detrend
export detrend!

"""
    detrend(s; <keyword arguments>)

Remove a trend.

# Arguments

- `s::AbstractVector`: signal vector; must contain at least 2 elements
- `type::Symbol=:linear`: detrending method:
    - `:linear`: subtract a 1st-order polynomial (linear trend); equivalent to `:poly` with `order=1`
    - `:ls`: subtract a linear least-squares fit using an efficient closed-form solution
    - `:poly`: subtract a polynomial of degree `order` fitted with `Polynomials.fit`
    - `:loess`: subtract a Loess smooth fitted with smoothing factor `f`
    - `:mean`: subtract the signal mean
    - `:constant`: subtract `offset`
- `offset::Real=0`: constant subtracted when `type = :constant`
- `order::Int64=1`: polynomial degree for `type = :poly`; must be ≥ 1 and < `length(s)`
- `f::Float64=1.0`: loess smoothing span ∈ `(0, 1]` for `type = :loess`

# Returns

- `Vector{Float64}`: detrended signal of the same length as `s`

# Throws

- `ArgumentError`: if `type` is invalid, `f ≤ 0`, or `order < 1`

# See also

[`detrend(::AbstractArray)`](@ref), [`detrend(::NeuroAnalyzer.NEURO)`](@ref)
"""
function detrend(
    s::AbstractVector;
    type::Symbol = :linear,
    offset::Real = 0,
    order::Int64 = 1,
    f::Float64 = 1.0
)::Vector{Float64}

    _check_var(type, [:ls, :linear, :mean, :constant, :poly, :loess], "type")
    f > 0 || throw(ArgumentError("f must be > 0."))
    order >= 1 || throw(ArgumentError("order must be ≥ 1."))

    if type === :loess

        t = collect(1.0:length(s))
        model = Loess.loess(t, Vector{Float64}(s); span=f)
        return s .- Loess.predict(model, t)

    elseif type === :poly

        t = collect(1:length(s))
        p = Polynomials.fit(t, s, order)
        trend = [p(ti) for ti in t]
        return s .- trend

    elseif type === :mean

        return s .- mean(s)

    elseif type === :constant

        return s .- offset

    elseif type === :ls

        T = eltype(s)
        N = length(s)
        # build design matrix A = [t  1] with t ∈ [0, 1]
        A      = Matrix{T}(undef, N, 2)
        A[:, 1] = range(T(0), T(1); length=N)
        A[:, 2] .= one(T)
        # closed-form OLS: coefficients = (AᵀA)⁻¹ Aᵀ s
        R = A' * A
        Rinv = inv(R)
        return s .- A * (Rinv * (A' * s))

    elseif type === :linear

        # Least-squares linear fit using the backslash operator
        t = collect(1.0:length(s))
        A = hcat(t, ones(length(s)))
        coef = A \ s
        trend = A * coef
        return s .- trend

    end

end

"""
    detrend(s; <keyword arguments>)

Remove a trend.

# Arguments

- `s::AbstractArray`: signal array, shape `(channels, samples, epochs)`
- `type::Symbol=:linear`: detrending method:
    - `:linear`: subtract a 1st-order polynomial (linear trend); equivalent to `:poly` with `order=1`
    - `:ls`: subtract a linear least-squares fit using an efficient closed-form solution
    - `:poly`: subtract a polynomial of degree `order` fitted with `Polynomials.fit`
    - `:loess`: subtract a Loess smooth fitted with smoothing factor `f`
    - `:mean`: subtract the signal mean
    - `:constant`: subtract `offset`
- `offset::Real=0`: constant subtracted when `type = :constant`
- `order::Int64=1`: polynomial degree for `type = :poly`; must be ≥ 1 and < `length(s)`
- `f::Float64=1.0`: loess smoothing span ∈ `(0, 1]` for `type = :loess`

# Returns

- `Array{Float64, 3}`: detrended array of the same shape as `s`

# Throws

- `ArgumentError`: if `s` is not 3-dimensional

# See also

[`detrend(::AbstractVector)`](@ref), [`detrend(::NeuroAnalyzer.NEURO)`](@ref)
"""
function detrend(
        s::AbstractArray; type::Symbol = :linear, offset::Real = 0, order::Int64 = 1, f::Float64 = 1.0
    )::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = detrend(
            @view(s[ch_idx, :, ep_idx]),
            type=type,
            offset=offset,
            order=order,
            f=f
        )
    end

    return s_new

end

"""
    detrend(obj; <keyword arguments>)

Remove a trend from selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `type::Symbol=:linear`: detrending method:
    - `:linear`: subtract a 1st-order polynomial (linear trend); equivalent to `:poly` with `order=1`
    - `:ls`: subtract a linear least-squares fit using an efficient closed-form solution
    - `:poly`: subtract a polynomial of degree `order` fitted with `Polynomials.fit`
    - `:loess`: subtract a Loess smooth fitted with smoothing factor `f`
    - `:mean`: subtract the signal mean
    - `:constant`: subtract `offset`
- `offset::Real=0`: constant subtracted when `type = :constant`
- `order::Int64=1`: polynomial degree for `type = :poly`; must be ≥ 1 and < `length(s)`
- `f::Float64=1.0`: loess smoothing span ∈ `(0, 1]` for `type = :loess`

# Returns

- `NeuroAnalyzer.NEURO`: new object with detrended channels

# See also

[`detrend!`](@ref), [`detrend(::AbstractArray)`](@ref)
"""
function detrend(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    type::Symbol = :linear,
    offset::Real = 0,
    order::Int64 = 1,
    f::Float64 = 1.0
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = detrend(
        @view(obj.data[ch, :, :]),
        type=type,
        offset=offset,
        order=order,
        f=f
    )
    push!(obj_new.history, "detrend(OBJ, ch=$ch, type=$type, offset=$offset, order=$order, f=$f)")

    return obj_new

end

"""
    detrend!(obj; <keyword arguments>)

Remove a trend in-place from selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `type::Symbol=:linear`: detrending method:
    - `:linear`: subtract a 1st-order polynomial (linear trend); equivalent to `:poly` with `order=1`
    - `:ls`: subtract a linear least-squares fit using an efficient closed-form solution
    - `:poly`: subtract a polynomial of degree `order` fitted with `Polynomials.fit`
    - `:loess`: subtract a Loess smooth fitted with smoothing factor `f`
    - `:mean`: subtract the signal mean
    - `:constant`: subtract `offset`
- `offset::Real=0`: constant subtracted when `type = :constant`
- `order::Int64=1`: polynomial degree for `type = :poly`; must be ≥ 1 and < `length(s)`
- `f::Float64=1.0`: loess smoothing span ∈ `(0, 1]` for `type = :loess`

# Returns

- `Nothing`
"""
function detrend!(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        type::Symbol = :linear,
        offset::Real = 0,
        order::Int64 = 1,
        f::Float64 = 1.0
    )::Nothing

    obj_new = detrend(obj, ch = ch, type = type, offset = offset, order = order, f = f)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
