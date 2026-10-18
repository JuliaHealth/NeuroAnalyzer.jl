
# ---------------------------------------------------------------------------
# internal helper: compute all summary statistics for one clean vector
# returns a NamedTuple matching the shared column layout used by all methods
# ---------------------------------------------------------------------------
# (exact values; rounding is done only when the table is printed)
function _summary_stats(x_clean::Vector{Float64}, n_total::Int64)
    ms = n_total - length(x_clean)
    m = mean(x_clean)
    # variance / SD are undefined for a single observation
    v = length(x_clean) > 1 ? var(x_clean) : NaN
    s = sqrt(v)
    mn, mx = extrema(x_clean)
    q1 = quantile(x_clean, 0.25)
    me = median(x_clean)
    q3 = quantile(x_clean, 0.75)
    mo = Float64(StatsKit.mode(x_clean))
    return (; ms, m, v, s, mn, q1, me, q3, mx, mo)
end

# ---------------------------------------------------------------------------
# internal helpers shared by the matrix and varargs methods
# ---------------------------------------------------------------------------

"""Build the shared summary DataFrame, calling `get_data(idx)` per group."""
function _build_summary_df(get_data::Function, g::Vector{String})::DataFrame
    cnames = [
        :group, :n, :missing, :mean, :var, :std, :min,
        :Q1, :median, :Q3, :max, :mode,
    ]
    ctypes = [
        String, Int64, Int64, Float64, Float64, Float64,
        Float64, Float64, Float64, Float64, Float64, Float64,
    ]
    df = DataFrame([name => type[] for (name, type) in zip(cnames, ctypes)])

    for idx in eachindex(g)
        x_clean, n_total = get_data(idx)
        isempty(x_clean) &&
            throw(ArgumentError("Group $(g[idx]) has no non-missing observations."))
        st = _summary_stats(x_clean, n_total)
        push!(
            df,
            [
                g[idx], n_total, st.ms, st.m, st.v, st.s,
                st.mn, st.q1, st.me, st.q3, st.mx, st.mo,
            ],
        )
    end
    return df
end

"""Print the multi-group summary table."""
function _print_summary_table(g::Vector{String}, df::DataFrame, d::Int64)
    return make_table(;
        digits = d,
        header = Matrix{String}(reshape(["group"; g], 1, :)),
        data = Any[
            "n" reshape(df[!, :n], 1, :);
            "missing" reshape(df[!, :missing], 1, :);
            "mean" reshape(df[!, :mean], 1, :);
            "var" reshape(df[!, :var], 1, :);
            "std" reshape(df[!, :std], 1, :);
            "min" reshape(df[!, :min], 1, :);
            "Q1" reshape(df[!, :Q1], 1, :);
            "median" reshape(df[!, :median], 1, :);
            "Q3" reshape(df[!, :Q3], 1, :);
            "max" reshape(df[!, :max], 1, :);
            "mode" reshape(df[!, :mode], 1, :);
        ],
    )
end

export summary_stats

"""
    summary_stats(x; g, d)

Return summary statistics for a single vector, printing a formatted table.

# Arguments

- `x::AbstractVector`: input data; may contain `Missing` or `NaN` values (removed before computation)
- `g::String=""`: optional group label shown in the table header
- `d::Int64=4`: number of decimal places shown in the table (returned values are not rounded)

# Returns

Named tuple:

- `n::Int64`: total number of observations (including missing)
- `ms::Int64`: number of missing / NaN observations
- `m::Float64`: mean
- `v::Float64`: variance
- `s::Float64`: standard deviation
- `mn::Float64`: minimum
- `q1::Float64`: first quartile (25th percentile)
- `me::Float64`: median
- `q3::Float64`: third quartile (75th percentile)
- `mx::Float64`: maximum
- `mo::Float64`: mode
"""
function summary_stats(
    x::AbstractVector;
    g::String = "",
    d::Int64 = 4,
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
    # remove NaNs
    x_clean = rmna(x)
    # validate
    length(x_clean) > 0 || throw(ArgumentError("No non-missing observations remain in x."))

    n = length(x)
    st = _summary_stats(x_clean, n)

    make_table(;
        digits = d,
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

    return (
        n = n,
        ms = st.ms,
        m = st.m,
        v = st.v,
        s = st.s,
        mn = st.mn,
        q1 = st.q1,
        me = st.me,
        q3 = st.q3,
        mx = st.mx,
        mo = st.mo,
    )
end

"""
    summary_stats(x; g, d)

Return summary statistics for each column of a matrix, printing a formatted table.

# Arguments

- `x::AbstractMatrix`: data matrix; each column is treated as one group. may contain `Missing` or `NaN` values
- `g::Vector{String}`: group labels; length must equal `size(x, 2)`
- `d::Int64=3`: number of decimal places shown in the table (returned values are not rounded)

# Returns

- `DataFrame`: one row per group with columns `:group`, `:n`, `:missing`, `:mean`, `:var`, `:std`, `:min`, `:Q1`,
  `:median`, `:Q3`, `:max`, `:mode`
"""
function summary_stats(
    x::AbstractMatrix;
    g::Vector{String},
    d::Int64 = 3,
)::DataFrame
    # validate
    length(g) == size(x, 2) || throw(
        ArgumentError("length(g) ($(length(g))) must equal size(x, 2) ($(size(x, 2)))."),
    )

    df = _build_summary_df(g) do idx
        return rmna(x[:, idx]), length(x[:, idx])
    end

    _print_summary_table(g, df, d)

    return df
end

"""
    summary(x1, x2, xs...; g, d)

Return summary statistics for two or more data vectors (groups), printing a formatted table.

This is a method of `Base.summary()` (it takes at least two arrays, so `summary(x)` and `summary(io, x)` are unchanged). For a single vector or the columns of a matrix use [`summary_stats`](@ref).

# Arguments

- `x1, x2, xs...::AbstractArray`: two or more data vectors (passed as positional arguments)
- `g::Vector{String}`: group labels; length must equal the number of vectors
- `d::Int64=3`: number of decimal places shown in the table (returned values are not rounded)

# Returns

- `DataFrame`: one row per group with columns `:group`, `:n`, `:missing`, `:mean`, `:var`, `:std`, `:min`, `:Q1`,
  `:median`, `:Q3`, `:max`, `:mode`
"""
function Base.summary(
    x1::AbstractArray,
    x2::AbstractArray,
    xs::AbstractArray...;
    g::Vector{String},
    d::Int64 = 3,
)::DataFrame
    x = (x1, x2, xs...)
    # validate
    length(g) == length(x) || throw(
        ArgumentError(
            "length(g) ($(length(g))) must equal the number of arrays ($(length(x))).",
        ),
    )

    df = _build_summary_df(g) do idx
        return rmna(x[idx]), length(x[idx])
    end

    _print_summary_table(g, df, d)

    return df
end
