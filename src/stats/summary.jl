export summary

# ---------------------------------------------------------------------------
# internal helper: compute all summary statistics for one clean vector
# returns a NamedTuple matching the shared column layout used by all methods
# ---------------------------------------------------------------------------
function _summary_stats(x_clean::Vector{Float64}, n_total::Int64, d::Int64)
    ms = n_total - length(x_clean)
    m = round(mean(x_clean), digits=d)
    v = round(var(x_clean), digits=d)
    s = round(std(x_clean), digits=d)
    mn, mx= round.(extrema(x_clean), digits=d)
    q1 = round(quantile(x_clean, 0.25), digits=d)
    me = round(median(x_clean), digits=d)
    q3 = round(quantile(x_clean, 0.75), digits=d)
    mo = round(Float64(mode(x_clean)), digits=d)
    return (; ms, m, v, s, mn, q1, me, q3, mx, mo)
end

# ---------------------------------------------------------------------------
# internal helpers shared by the matrix and varargs methods
# ---------------------------------------------------------------------------

"""Build the shared summary DataFrame, calling `get_data(idx)` per group."""
function _build_summary_df(get_data::Function, g::Vector{String}, d::Int64)::DataFrame
    cnames = [:group, :n, :missing, :mean, :var, :std, :min,
              :Q1, :median, :Q3, :max, :mode]
    ctypes = [String, Int64, Int64, Float64, Float64, Float64,
              Float64, Float64, Float64, Float64, Float64, Float64]
    df = DataFrame([name => type[] for (name, type) in zip(cnames, ctypes)])

    for idx in eachindex(g)
        x_clean, n_total = get_data(idx)
        st = _summary_stats(x_clean, n_total, d)
        push!(df, [g[idx], n_total, st.ms, st.m, st.v, st.s,
                   st.mn, st.q1, st.me, st.q3, st.mx, st.mo])
    end
    return df
end

"""Print the multi-group summary table."""
function _print_summary_table(g::Vector{String}, df::DataFrame)
    make_table(
        header = Matrix{String}(reshape(["group"; g], 1, :)),
        data = Any[
            "n" reshape(df[!, :n], 1, :);
            "missing" reshape(df[!, :missing], 1, :);
            "mean" reshape(df[!, :mean], 1, :);
            "var" reshape(df[!, :var], 1, :);
            "std" reshape(df[!, :std], 1, :);
            "min" reshape(df[!, :min], 1, :);
            "Q1" reshape(df[!, :Q1], 1, :);
            "median" reshape(df[!, :median],  1, :);
            "Q3" reshape(df[!, :Q3], 1, :);
            "max" reshape(df[!, :max], 1, :);
            "mode" reshape(df[!, :mode], 1, :);
        ],
    )
end

"""
    summary(x; g)

Return summary statistics for a single vector, printing a formatted table.

# Arguments
- `x::AbstractVector`: input data; may contain `Missing` or `NaN` values (removed before computation)
- `g::String=""`: optional group label shown in the table header

# Returns

Named tuple:

- `n::Int64`: total number of observations (including missing)
- `ms::Int64`: number of missing / NaN observations
- `m::Float64`: mean
- `v::Float64`: variance
- `s::Float64`: standard deviation
- `min::Float64`: minimum
- `q1::Float64`: first quartile (25th percentile)
- `me::Float64`: median
- `q3::Float64`: third quartile (75th percentile)
- `max::Float64`: maximum
- `mo::Float64`: mode

# Throws

- `ArgumentError`: if all values are missing/NaN (no data to summarise)

# See also

[`summary(::AbstractMatrix)`](@ref)
"""
function summary(
    x::AbstractVector;
    g::String = ""
)::@NamedTuple{
    n::Int64,
    ms::Int64,
    m::Float64,
    v::Float64,
    s::Float64,
    mn::Float64,
    q1::Float64,
    me::Float64,
    q3::Float64,
    mx::Float64,
    mo::Float64,
}

    x_clean = rmna(x)
    !(length(x_clean) > 0) && throw(ArgumentError("No non-missing observations remain in x."))

    n  = length(x)
    st = _summary_stats(x_clean, n, 4)   # use 4 d.p. for the scalar method

    make_table(
        header = Matrix{String}(["" g]),   # ensure correct type for make_table
        data = Any[
            "n" n;
            "missing" st.ms;
            "mean" st.m;
            "var" st.v;
            "std" st.s;
            "min" st.mn;
            "Q1" st.q1;
            "median" st.me;
            "Q3" st.q3;
            "max" st.mx;
            "mode" st.mo;
        ],
    )

    return (; ms, m, v, s, mn, q1, me, q3, mx, mo)

end

"""
    summary(x; g, d)

Return summary statistics for each column of a matrix, printing a formatted table.

# Arguments

- `x::AbstractMatrix`: data matrix; each column is treated as one group. may contain `Missing` or `NaN` values
- `g::Vector{String}`: group labels; length must equal `size(x, 2)`
- `d::Int64=3`: number of decimal places for rounding

# Returns

- `DataFrame`: one row per group with columns `:group`, `:n`, `:missing`, `:mean`, `:var`, `:std`, `:min`, `:Q1`,
  `:median`, `:Q3`, `:max`, `:mode`

# Throws

- `ArgumentError`: if `length(g) ≠ size(x, 2)`

# See also

[`summary(::AbstractVector)`](@ref), [`summary(::AbstractArray...)`](@ref)
"""
function summary(
    x::AbstractMatrix;
    g::Vector{String},
    d::Int64 = 3
)::DataFrame

    !(length(g) == size(x, 2)) && throw(ArgumentError("length(g) ($(length(g))) must equal size(x, 2) ($(size(x, 2)))."))

    df = _build_summary_df(g, d) do idx
        rmna(x[:, idx]), length(x[:, idx])
    end

    _print_summary_table(g, df)

    return df

end

"""
    summary(x; g, d)

Return summary statistics for each vector in a varargs list, printing a formatted table.

# Arguments

- `x::AbstractArray...`: one or more data vectors (passed as positional arguments)
- `g::Vector{String}`: group labels; length must equal the number of vectors
- `d::Int64=3`: number of decimal places for rounding

# Returns

- `DataFrame`: one row per group with columns `:group`, `:n`, `:missing`, `:mean`, `:var`, `:std`, `:min`, `:Q1`,
  `:median`, `:Q3`, `:max`, `:mode`

# Throws
- `ArgumentError`: if `length(g) ≠ length(x)`

# See also

[`summary(::AbstractMatrix)`](@ref)
"""
function summary(
    x::AbstractArray...;
    g::Vector{String},
    d::Int64 = 3
)::DataFrame

    !(length(g) == length(x)) && throw(ArgumentError("Number of group names ($length(g)) must be equal to the number of groups $(size(x, 2))."))

    !(length(g) == length(x)) && throw(ArgumentError("length(g) ($(length(g))) must equal the number of arrays ($(length(x)))."))

    df = _build_summary_df(g, d) do idx
        rmna(x[idx]), length(x[idx])
    end

    _print_summary_table(g, df)

    return df

end
