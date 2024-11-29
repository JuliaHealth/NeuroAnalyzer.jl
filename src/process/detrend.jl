export detrend
export detrend!

"""
    detrend(s; <keyword arguments>)

Perform piecewise detrending.

# Arguments

- `s::AbstractVector`
- `type::Symbol=:linear`:
    - `:loess`: fit loess approximation and subtract it from `s`
    - `:poly`: polynomial of `order` is subtracted from `s`
    - `:mean`: the mean of `s` is subtracted from `s`
    - `:constant`: `offset` is subtracted from `s`
    - `:ls`: the result of a linear least-squares fit to `s` is subtracted from `s`
    - `:linear`: linear trend (1st order polynomial) is subtracted from `s`
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

# Returns
- `s_new::Vector{Float64}`
"""
function detrend(s::AbstractVector; type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)::Vector{Float64}

    _check_var(type, [:ls, :linear, :mean, :constant, :poly, :loess, :hp], "type")
    @assert f > 0 "f must be > 0."
    @assert order >= 1 "order must be ≥ 1."

    if type === :loess
        t = collect(1.0:1:length(s))
        model = Loess.loess(t, Vector(s), span=f)
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
    elseif type === :mean
        s_new = s .- mean(s)
    elseif type === :constant
        s_new = s .- offset
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
        t = collect(1:1:length(s))
        p = Polynomials.fit(t, s, 1)
        trend = zeros(length(s))
        for idx in eachindex(s)
            trend[idx] = p(t[idx])
        end
        s_new = s .- trend
    end

    return s_new

end

"""
    detrend(s; <keyword arguments>)

Perform piecewise detrending.

# Arguments

- `s::AbstractArray`
- `type::Symbol=:linear`: detrending method
    - `:loess`: fit loess approximation and subtract it from `s`
    - `:poly`: polynomial of `order` is subtracted from `s`
    - `:mean`: the mean of `s` is subtracted from `s`
    - `:constant`: `offset` is subtracted from `s`
    - `:ls`: the result of a linear least-squares fit to `s` is subtracted from `s`
    - `:linear`: linear trend (1st order polynomial) is subtracted from `s`
- `offset::Real=0`: constant for `:constant` detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

# Returns

- `s_new::Array{Float64, 3}`
"""
function detrend(s::AbstractArray; type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views detrend(s[ch_idx, :, ep_idx], type=type, offset=offset, order=order, f=f)
        end
    end

    return s_new

end

"""
    detrend(obj; <keyword arguments>)

Perform piecewise detrending.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `type::Symbol=:linear`: detrending method
    - `:loess`: fit loess approximation and subtract it from `s`
    - `:poly`: polynomial of `order` is subtracted from `s`
    - `:mean`: the mean of `s` is subtracted from `s`
    - `:constant`: `offset` is subtracted from `s`
    - `:ls`: the result of a linear least-squares fit to `s` is subtracted from `s`
    - `:linear`: linear trend (1st order polynomial) is subtracted from `s`
- `offset::Real=0`: constant for `:constant` detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function detrend(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = detrend(obj.data[ch, :, :], type=type, offset=offset, order=order, f=f)
    reset_components!(obj_new)
    push!(obj_new.history, "detrend(OBJ, ch=$ch, type=$type, offset=$offset, order=$order, f=$f)")

    return obj_new

end

"""
    detrend!(obj; <keyword arguments>)

Perform piecewise detrending.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `type::Symbol=:linear`: detrending method
    - `:loess`: fit loess approximation and subtract it from `s`
    - `:poly`: polynomial of `order` is subtracted from `s`
    - `:mean`: the mean of `s` is subtracted from `s`
    - `:constant`: `offset` is subtracted from `s`
    - `:ls`: the result of a linear least-squares fit to `s` is subtracted from `s`
    - `:linear`: linear trend (1st order polynomial) is subtracted from `s`
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

# Returns

Nothing
"""
function detrend!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, type::Symbol=:linear, offset::Real=0, order::Int64=1, f::Float64=1.0)::Nothing

    obj_new = detrend(obj, ch=ch, type=type, offset=offset, order=order, f=f)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
