export linspace
export logspace
export cmax
export cmin
export tuple_order
export cums
export f_nearest

"""
    linspace(start, stop, length)

Generates a sequence of evenly spaced numbers between `start` and `stop`.

# Arguments

- `start::Number`
- `stop::Number`
- `n::Int64`: sequence length

# Returns

- `range::Vector`
"""
function linspace(start::Number, stop::Number, n::Int64)

    n < 2 && throw(ArgumentError("n must be ≥ 2."))

    return collect(range(start, stop, n))

end

"""
    logspace(start, stop, n)

Generates a sequence of log10-spaced numbers between `start` and `stop`.

# Arguments

- `start::Number`
- `stop::Number`
- `n::Int64`: sequence length

# Returns

- `range::Vector{<:Number}`
"""
function logspace(start::Number, stop::Number, n::Int64)

    n < 2 && throw(ArgumentError("n must be ≥ 2."))

    return collect(exp10.(range(start, stop, n)))

end

"""
    cmax(x)

Return maximum value of the complex vector.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmax::ComplexF64`
"""
function cmax(x::Vector{<:Complex})

    return argmax(abs, x)

end

"""
    cmin(x)

Return minimum value of the complex vector.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmin::ComplexF64`
"""
function cmin(x::Vector{<:Complex})

    return argmin(abs, x)

end

"""
    tuple_order(t, rev)

Order tuple elements in ascending or descending (`rev=true`) order.

# Arguments

- `t::Tuple{Real, Real}`
- `rev::Bool=false`

# Returns

- `t::Tuple{Real, Real}`
"""
function tuple_order(t::Tuple{Real, Real}, rev::Bool=false)

    (rev == false && t[1] > t[2]) && (t = (t[2], t[1]))
    (rev == true && t[1] < t[2]) && (t = (t[2], t[1]))

    return t

end

"""
    cums(signal)

Calculate cumulative sum.

# Arguments

- `signal::Array{<:Real, 3}`

# Returns

- `signal_cs::Array{Float64, 3}`
"""
function cums(signal::Array{<:Real, 3})
    
    ch_n, _, ep_n = size(signal)
    signal_cs = similar(signal)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            signal_cs[ch_idx, :, ep_idx] = @views cumsum(signal[ch_idx, :, ep_idx])
        end
    end

    return signal_cs

end

"""
    f_nearest(m, pos)

Find nearest position tuple `pos` in vector of positions `m`.

# Arguments

- `m::Matrix{Tuple{Float64, Float64}}`
- `p::Tuple{Float64, Float64}`

# Returns

- `pos::Tuple{Int64, Int64}`: row and column in m
"""
function f_nearest(m::Matrix{Tuple{Float64, Float64}}, p::Tuple{Float64, Float64})

    d = zeros(size(m))

    @inbounds @simd for idx1 in 1:size(m, 1)
        for idx2 in 1:size(m, 2)
            d[idx1, idx2] = euclidean(m[idx1, idx2], p)
        end
    end

    return (findmin(d)[2][1], findmin(d)[2][2])

end
