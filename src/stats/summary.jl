export summary

"""
    summary(x; g)

Return summary statistics.

# Arguments

- `x::AbstractVector`
- `g::String=""`: group name

# Returns

Named tuple containing:
- `n::Int64`: number of observations
- `ms::Int64`: missings
- `mm::Float64`: mean
- `v::Float64`: variance
- `s::Float64`: standard deviation
- `min::Float64`: minimum
- `q1::Float64`: 1st quartile
- `md::Float64`: median
- `q3::Float64`: 3rd quartile
- `max::Float64`: maximum
- `mo::Float64`: mode
"""
function summary(x::AbstractVector; g::String="")::@NamedTuple{n::Int64, ms::Int64, m::Float64, v::Float64, s::Float64, min::Float64, q1::Float64, me::Float64, q3::Float64, max::Float64, mo::Float64}

    x_tmp = na(x)
    n = length(x)
    ms = n - length(x_tmp)
    m = mean(x_tmp)
    v = var(x_tmp)
    s = std(x_tmp)
    min, max = extrema(x_tmp)
    me = median(x_tmp)
    q1 = quantile(x_tmp, 0.25)
    q3 = quantile!(x_tmp, 0.75)
    mo = mode(x_tmp)

    make_table(header=["" g];
               data=["n" n;
                     "missing" ms;
                     "mean" m;
                     "var" v;
                     "std" s;
                     "min" min;
                     "Q1" q1;
                     "median" me;
                     "Q3" q3;
                     "max" max;
                     "mode" mo])

    return (n=n, ms=ms,  m=m, v=v, s=s, min=min, q1=q1, me=me, q3=q3, max=max, mo=mo)

end

"""
    summary(x; g, d)

Return summary statistics.

# Arguments

- `x::AbstractMatrix`
- `g::Vector{String}`: group names
- `d::Int64`: round to `d` digits

# Returns

`- df::DataFrame` containing:
    - `:group`: group name
    - `:n`: number of observations
    - `:missing`: missings
    - `:m`: mean
    - `:v`: variance
    - `:s`: standard deviation
    - `:min`: minimum
    - `:Q1`: 1st quartile
    - `:median`: median
    - `:Q3`: 3rd quartile
    - `:max`: maximum
    - `:mode`: mode
"""
function summary(x::AbstractMatrix; g::Vector{String}, d::Int64=3)::DataFrame

    @assert length(g) == size(x, 2) "Number of group names ($length(g)) must be equal to the number of groups $(size(x, 2))."

    cnames = [:group, :n, :missing, :mean, :var, :std, :min, :Q1, :median, :Q3, :max, :mode]
    ctypes = [String, Int64, Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64]
    cols = (; zip(cnames, type[] for type in ctypes )...)
    df = DataFrame(cols)

    for idx in axes(x, 2)
        x_tmp = na(x[:, idx])
        n = length(x[:, idx])
        ms = n - length(x_tmp)
        m = round(mean(x_tmp), digits=d)
        v = round(var(x_tmp), digits=d)
        s = round(std(x_tmp), digits=d)
        min, max = round.(extrema(x_tmp), digits=d)
        me = round(median(x_tmp), digits=d)
        q1 = round(quantile(x_tmp, 0.25), digits=d)
        q3 = round(quantile!(x_tmp, 0.75), digits=d)
        mo = round(mode(x_tmp), digits=d)
        push!(df, [g[idx], n, ms, m, v, s, min, q1, me, q3, max, mo])
    end

    make_table(header=["group" reshape(g[:, :], 1, :)];
               data=["n" reshape(df[:, :n], 1, :);
                     "missing" reshape(df[:, :missing], 1, :);
                     "mean" reshape(df[:, :mean], 1, :);
                     "var" reshape(df[:, :var], 1, :);
                     "std" reshape(df[:, :std], 1, :);
                     "min" reshape(df[:, :min], 1, :);
                     "Q1" reshape(df[:, :Q1], 1, :);
                     "median" reshape(df[:, :median], 1, :);
                     "Q3" reshape(df[:, :Q3], 1, :);
                     "max" reshape(df[:, :max], 1, :);
                     "mode" reshape(df[:, :mode], 1, :)])

    return df

end

"""
    summary(x; g, d)

Return summary statistics.

# Arguments

- `x::AbstractMatrix`
- `g::Vector{String}`: group names
- `d::Int64`: round to `d` digits

# Returns

`- df::DataFrame` containing:
    - `:group`: group name
    - `:n`: number of observations
    - `:missing`: missings
    - `:m`: mean
    - `:v`: variance
    - `:s`: standard deviation
    - `:min`: minimum
    - `:Q1`: 1st quartile
    - `:median`: median
    - `:Q3`: 3rd quartile
    - `:max`: maximum
    - `:mode`: mode
"""
function summary(x::AbstractArray...; g::Vector{String}, d::Int64=3)::DataFrame

    @assert length(g) == length(x) "Number of group names ($length(g)) must be equal to the number of groups $(size(x, 2))."

    cnames = [:group, :n, :missing, :mean, :var, :std, :min, :Q1, :median, :Q3, :max, :mode]
    ctypes = [String, Int64, Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64]
    cols = (; zip(cnames, type[] for type in ctypes )...)
    df = DataFrame(cols)

    for idx in eachindex(x)
        x_tmp = na(x[idx])
        n = length(x[idx])
        ms = n - length(x_tmp)
        m = round(mean(x_tmp), digits=d)
        v = round(var(x_tmp), digits=d)
        s = round(std(x_tmp), digits=d)
        min, max = round.(extrema(x_tmp), digits=d)
        me = round(median(x_tmp), digits=d)
        q1 = round(quantile(x_tmp, 0.25), digits=d)
        q3 = round(quantile!(x_tmp, 0.75), digits=d)
        mo = round(mode(x_tmp), digits=d)
        push!(df, [g[idx], n, ms, m, v, s, min, q1, me, q3, max, mo])
    end

    make_table(header=["group" reshape(g[:, :], 1, :)];
               data=["n" reshape(df[:, :n], 1, :);
                     "missing" reshape(df[:, :missing], 1, :);
                     "mean" reshape(df[:, :mean], 1, :);
                     "var" reshape(df[:, :var], 1, :);
                     "std" reshape(df[:, :std], 1, :);
                     "min" reshape(df[:, :min], 1, :);
                     "Q1" reshape(df[:, :Q1], 1, :);
                     "median" reshape(df[:, :median], 1, :);
                     "Q3" reshape(df[:, :Q3], 1, :);
                     "max" reshape(df[:, :max], 1, :);
                     "mode" reshape(df[:, :mode], 1, :)])

    return df

end
