export vsearch
export vsplit
export minat
export maxat
export vreduce

"""
    vsearch(y, x)

Return the positions of the string in the vector of strings.

# Arguments

- `y::String`: value of interest
- `x::Vector{String}`: vector to search within

# Returns

- `idx::Union{Int64, Nothing}`
"""
function vsearch(y::String, x::Vector{String})::Union{Int64, Nothing}

    idx = nothing
    if y in x
        idx = findfirst(isequal(y), x)
    end

    return idx

end

"""
    vsearch(y, x; <keyword arguments>)

Return the positions of the value in the vector.

# Arguments

- `y::Real`: value of interest
- `x::AbstractVector`: vector to search within
- `acc::Bool=false`: if true, return the difference between `y` and `x[idx]`

# Returns

- `idx::Int64`
- `d::Real`: the difference between `y` and `x[idx]`
"""
function vsearch(y::Real, x::AbstractVector; acc::Bool=false)::Union{Int64, Tuple{Int64, Real}}

    @assert length(x) > 0 "Length of x must be > 0."
    d, idx = findmin(abs.(x .- y))

    return acc ? (idx, d) : idx

end

"""
    vsearch(y, x; <keyword arguments>)

Return the positions of the value in the vector.

# Arguments

- `y::AbstractVector`: vector of interest
- `x::AbstractVector`: vector to search within
- `acc::Bool=false`: if true, return the difference between `y` and `x[idx:idx + length(y)]`

# Returns

- `idx::Int64`
- `d::Real`: the difference between `y` and `x[idx:idx + length(y)]`
"""
function vsearch(y::AbstractVector, x::AbstractVector; acc::Bool=false)::Union{AbstractVector, Tuple{AbstractVector, AbstractVector}}

    @assert length(x) > 0 "Length of x must be > 0."
    @assert length(y) <= length(x) "Length of y must be ≤ length 'x'"

    idx = zeros(length(y))
    d = zeros(length(y))

    @inbounds for y_idx in eachindex(y)
        d[y_idx], idx[y_idx] = findmin(abs.(x .- y[y_idx]))
    end

    return acc ? (convert.(Int64, idx), d) : idx

end

"""
    vsplit(x, n)

Splits vector into pieces.

# Argument

- `x::AbstractVector`
- `n::Int64`: length of one piece

# Returns

- `x::Vector{AbstractVector}`
"""
function vsplit(x::AbstractVector, n::Int64=1)::Vector{AbstractVector}

    @assert length(x) > 0 "Length of x must be > 0."
    @assert n >= 1 "n must be ≥ 1."
    @assert length(x) % n == 0 "Length of x must be a multiple of n."

    x_m = reshape(x, length(x) ÷ n, n)
    result = [x_m[1, :]]
    @inbounds for idx in 2:size(x_m, 1)
        result = vcat(result, [x_m[idx, :]])
    end

    return result

end

"""
    minat(x, y)

Find minimum value of one vector and return value at its index from another vector.

# Argument

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `value::Real`
- `idx::Int64`
"""
function minat(x::AbstractVector, y::AbstractVector)::Tuple{Real, Int64}

    @assert length(x) > 0 "Length of x must be > 0."
    @assert length(y) > 0 "Length of y must be > 0."
    @assert length(x) == length(x) "x and y length must be equal."

    idx = vsearch(minimum(x), x)
    value = y[idx]

    return value, idx

end

"""
    maxat(x, y)

Find maximum value of one vector and return value at its index from another vector.

# Argument

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `value::Real`
- `idx::Int64`
"""
function maxat(x::AbstractVector, y::AbstractVector)::Tuple{Real, Int64}

    @assert length(x) > 0 "Length of x must be > 0."
    @assert length(y) > 0 "Length of y must be > 0."
    @assert length(x) == length(x) "Lengths of x and y must be equal."

    idx = vsearch(maximum(x), x)
    value = y[idx]

    return value, idx

end

"""
    vreduce(x, f; <keyword arguments>)

Reduce two vectors at indices of the second vector being multiplications of a constant. Useful e.g. for simplifying values across frequencies, when the number of frequencies (and thus values) is high.

# Arguments

- `x::AbstractVector`: e.g. signal data
- `f::AbstractVector`: e.g. frequencies
- `n::Float64=0.5`: reduce at multiplications of this value

# Returns

- `x_new::AbstractVector`
- `f_new::AbstractVector`
"""
function vreduce(x::AbstractVector, f::AbstractVector; n::Float64=0.5)::Tuple{AbstractVector, AbstractVector}

    @assert length(x) > 0 "Length of x must be > 0."
    @assert length(f) > 0 "Length of f must be > 0."
    @assert length(x) == length(f) "Lengths of x and f must be equal."

    f1_idx = vsearch(round(f[1]), f)
    f2_idx = vsearch(round(f[end]), f)
    f1 = round(f[f1_idx])
    f2 = round(f[f2_idx])

    f_new = collect(f1:n:f2)
    x_new = zeros(length(f_new))
    @inbounds for idx in eachindex(f_new)
        f_idx = vsearch(f_new[idx], f)
        x_new[idx] = x[f_idx]
    end

    return x_new, f_new

end
