export linspace
export logspace
export cmax
export cmin
export tuple_order

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

Return maximum value of the complex vector`x`.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmax::ComplexF64`
"""
function cmax(x::Vector{ComplexF64})
    return argmax(abs, x)
end

"""
    cmin(x)

Return minimum value of the complex vector`x`.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmin::ComplexF64`
"""
function cmin(x::Vector{ComplexF64})
    return argmin(abs, x)
end

"""
    tuple_order(t, rev)

Order tuple elements in ascending or descending (rev=true) order.

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

