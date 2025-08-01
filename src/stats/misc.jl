export k_categories
export slope
export distance
export count_thresh
export cmp_stat
export permute
export logit
export ss
export na
export df
export center

"""
    k_categories(n)

Calculate number of categories for a given sample size `n`.

# Arguments

- `n::Int64`: sample size

# Returns

Named tuple containing:
- `k1::Float64`: sqrt(n)
- `k2::Float64`: 1 + 3.222 * log10(n)
"""
function k_categories(n::Int64)::@NamedTuple{k1::Float64, k2::Float64}

    k1 = sqrt(n)
    k2 = 1 + 3.222 * log10(n)

    return (k1=k1, k2=k2)

end

"""
    slope(p1, p2)

Calculate slope of the line crossing two points.

# Arguments

- `p1::Tuple{Real, Real}`
- `p2::Tuple{Real, Real}`

# Returns

- `s::Float64`: slope
"""
function slope(p1::Tuple{Real, Real}, p2::Tuple{Real, Real})::Float64

    @assert p2[1] - p1[1] != 0 "p2[1] - p1[1] must not be 0."

    s = (p2[2] - p1[2]) / (p2[1] - p1[1])

    return s

end

"""
    distance(p1, p2)

Calculate distance between two points.

# Arguments

- `p1::Tuple{Real, Real}`
- `p2::Tuple{Real, Real}`

# Returns

- `d::Float64`: distance
"""
function distance(p1::Tuple{Real, Real}, p2::Tuple{Real, Real})::Float64

    d = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)

    return d

end

"""
    count_thresh(x; <keyword arguments>)

Collect thresholded elements, e.g. in a topographical map.

# Arguments

- `x::AbstractMatrix`
- `t::Real`: threshold value
- `t_type::Symbol=:g`: rule for thresholding:
    - `:eq`: =
    - `:geq`: ≥
    - `:leq`: ≤
    - `:g`: >
    - `:l`: <

# Returns

Named tuple containing:
- `x_t::Matrix{Bool}`: thresholded matrix
- `n::Int64`: number of elements
"""
function count_thresh(x::AbstractMatrix; t::Real, t_type::Symbol=:g)::@NamedTuple{x_t::Matrix{Bool}, n::Int64}

    _check_var(t_type, [:eq, :geq, :leq, :g, :l], "t_type")

    x_t = zeros(Bool, size(x))

    if t_type === :eq
        x_t[x .== t] .= true
    elseif t_type === :g
        x_t[x .> t] .= true
    elseif t_type === :geq
        x_t[x .>= t] .= true
    elseif t_type === :l
        x_t[x .< t] .= true
    elseif t_type === :leq
        x_t[x .<= t] .= true
    end

    n = count(==(true), x_t)

    return (x_t=x_t, n=n)

end

"""
    cmp_stat(stat_dist, v)

Calculate proportion of elements below or above a given statistic value.

# Arguments

- `stat_dist::AbstractVector`: statistic values distribution
- `v::Real`: statistic value
- `type::Symbol=:g`: calculation proportion of elements greater (`:g`) or lesser (`:l`) than `v`

# Returns

- `p::Float64`
"""
function cmp_stat(stat_dist::AbstractVector, v::Real; type::Symbol=:g)::Float64

    _check_var(type, [:g, :l], "type")
    @assert length(stat_dist) > 0 "Length of stat_dist is 0, cannot compute."

    type === :g && return count(stat_dist .> v) / length(stat_dist)
    type === :l && return count(stat_dist .< v) / length(stat_dist)

end

"""
    permute(s, n)

Permute signal data.

# Arguments

- `s::AbstractVector`
- `n::Int64`: number of permutations

# Returns

- `s_new::Matrix{Float64}`
"""
function permute(s::AbstractVector, n::Int64)::Matrix{Float64}

    @assert n > 0 "n must be > 0."

    s_new = zeros(n, length(s))
    for idx in 1:n
        x = rand(2:length(s))
        s1 = s[x:end]
        s2 = s[1:(x - 1)]
        s_new[idx, :] = @views vcat(s1, s2)
    end

    return s_new

end

"""
    permute(s, n)

Permute signal data.

# Arguments

- `s::AbstractArray`
- `n::Int64`: number of permutations

# Returns

- `s_new::Union{Array{Float64, 3}, Array{Float64, 4}}`
"""
function permute(s::AbstractArray, n::Int64)::Union{Array{Float64, 3}, Array{Float64, 4}}

    @assert n > 0 "n must be > 0."
    @assert ndims(s) <= 3 "permute() only works for arrays of ≤ 3 dimensions."

    if ndims(s) == 2
        s_new = zeros(n, size(s,1 ), size(s,2 ))
        @inbounds for idx1 in 1:n
            Threads.@threads :greedy for idx2 in axes(s, 1)
                x = rand(2:size(s, 2))
                s1 = s[idx2, x:end]
                s2 = s[idx2, 1:(x - 1)]
                s_new[idx1, idx2, :] = @views vcat(s1, s2)
            end
        end
    else
        s_new = zeros(n, size(s, 1), size(s, 2), size(s, 3))
        for idx1 in 1:n
            @inbounds for idx2 in axes(s, 1)
                Threads.@threads :greedy for idx3 in axes(s, 3)
                    x = rand(2:size(s, 2))
                    s1 = s[idx2, x:end, idx3]
                    s2 = s[idx2, 1:(x - 1), idx3]
                    s_new[idx1, idx2, :, idx3] = @views vcat(s1, s2)
                end
            end
        end
    end

    return s_new

end

"""
    logit(p)

Convert proportion to logit.

# Arguments

- `p::Float64`: proportion

# Returns

- `l::Float64`
"""
function logit(p::Float64)::Float64

    _in(p, (0.0, 1.0), "p")

    l = log(p / (1 - p))

    return l

end

"""
    ss(x)

Calculate sum of squares.

# Arguments

- `x::AbstractVector`

# Returns

- `s::Float64`
"""
function ss(x::AbstractVector)::Float64

    m = mean(x)
    s = sum((x .- m).^2)

    return s

end

"""
    na(x)

Return values of x ignoring NaNs and Missing values.

# Arguments

- `x::AbstractVector`

# Returns

- `x::Vector{Float64}`
"""
function na(x::AbstractVector)::Vector{Float64}

    x = x[.!ismissing.(x)]
    x = x[.!isnan.(x)]

    return Float64.(x)

end


"""
    df(x)

Calculate degrees of freedom.

# Arguments

- `x::AbstractVector`

# Returns

- `df::Int64`
"""
function df(x::AbstractVector)::Int64

    return length(x) - 1

end

"""
    center(x)

Center values by subtracting mean.

# Arguments

- `x::AbstractVector`

# Returns

- `x::Vector{Float64}`
"""
function center(x::AbstractVector)::Vector{Float64}

    m = mean(x)

    return x .- m

end
