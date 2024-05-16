export outlier_detect
export grubbs

"""
    outlier_detect(x; method)

Detect outliers.

# Arguments

- `x::AbstractVector`
- `method::Symbol=iqr`: detecting methods:
    - `:iqr`: interquartile range
    - `:z`: z-score
    - `:g`: Grubbs test

# Returns

- `o::Vector{Bool}`: index of outliers
"""
function outlier_detect(x::AbstractVector; method::Symbol=:iqr)

    _check_var(method, [:iqr, :z, :g], "method")
    o = zeros(Bool, length(x))

    if method === :iqr
        m1 = quantile(x, 0.25) - 1.5 * iqr(x)
        m2 = quantile(x, 0.75) + 1.5 * iqr(x)
        o[x .< m1] .= true
        o[x .> m2] .= true
    elseif method === :z
        z = z_score(x)
        o[z .< -3] .= true
        o[z .> 3] .= true
    else
        @assert length(x) > 6 "For :g method length(x) must be > 6."
        x_tmp = deepcopy(x)
        for _ in length(x_tmp):-1:6
            _, m_idx = findmax(x_tmp)
            if grubbs(x_tmp, t=1)
                o[m_idx] = true
                deleteat!(x_tmp, m_idx)
            end
        end
        x_tmp = deepcopy(x)
        for _ in length(x_tmp):-1:6
            _, m_idx = findmin(x_tmp)
            if grubbs(x_tmp, t=-1)
                o[m_idx] = true
                deleteat!(x_tmp, m_idx)
            end
        end
    end

    return o

end

"""
    grubbs(x; alpha, t)

Perform Grubbs test for outlier.

# Arguments

- `x::AbstractVector`
- `alpha::Float64=0.95`
- `t::Int64=0`: test type:
    - `-1`: test whether the minimum value is an outlier
    - `0`: two-sided test
    - `1`: test whether the maximum value is an outlier

# Returns

- `g::Bool`: true: outlier exists, false: there is no outlier
"""
function grubbs(x::AbstractVector; alpha::Float64=0.95, t::Int64=0)

    n = length(x)
    df = n - 2

    @assert t in [-1, 0, 1] "t must be -1, 0 or 1."

    if t == 0
        two_sided = true
        g = maximum(abs.(x .- mean(x))) / std(x)
    elseif t == -1
        two_sided = false
        g = (mean(x) - minimum(x)) / std(x)
    elseif t == 1
        two_sided = false
        g = (maximum(x) - mean(x)) / std(x)
    end

    p = two_sided ? (1 - alpha) / (2 * n) : (1 - alpha) / n
    t_critical = quantile(TDist(df), 1 - p)
    h = (n - 1) * t_critical / sqrt(n * (df + t_critical^2))

    return g < h ? false : true

end

