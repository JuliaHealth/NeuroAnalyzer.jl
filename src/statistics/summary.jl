export summary

"""
    summary(x)

Return summary statistics.

# Arguments

- `x::AbstractVector`

# Returns

Named tuple containing:
- `mm::Float64`: mean
- `s::Float64`: standard deviation
- `v::Float64`: variance
- `me::Float64`: median
- `mo::Float64`: mode
"""
function summary(x::AbstractVector; g::String="")::@NamedTuple{mm::Float64, s::Float64, v::Float64, me::Float64, mo::Float64}

    mm = mean(x)
    v = var(x)
    s = std(x)
    me = median(x)
    mo = mode(x)

    if verbose
        make_table(header=["" g];
                   data=["N" length(x);
                         "mean" mm;
                         "std dev" s;
                         "variance" v;
                         "median" me;
                         "mode" mo])
    end

    return (mm=mm, s=s, v=v, me=me, mo=mo)

end

"""
    summary(x, y)

Return summary statistics.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`
- `g1::String="1"`: group 1 name
- `g2::String="2"`: group 2 name

# Returns

Named tuple containing:
- `mm1::Float64`: mean
- `mm2::Float64`: mean
- `s1::Float64`: standard deviation
- `s2::Float64`: standard deviation
- `v1::Float64`: variance
- `v2::Float64`: variance
- `me1::Float64`: median
- `me2::Float64`: median
- `mo1::Float64`: mode
- `mo2::Float64`: mode
"""
function summary(x::AbstractVector, y::AbstractVector; g1::String="1", g2::String="2")::@NamedTuple{mm1::Float64, mm2::Float64, s1::Float64, s2::Float64, v1::Float64, v2::Float64, me1::Float64, me2::Float64, mo1::Float64, mo2::Float64}

    mm1 = mean(x)
    mm2 = mean(y)
    v1 = var(x)
    v2 = var(y)
    s1 = std(x)
    s2 = std(y)
    me1 = median(x)
    me2 = median(y)
    mo1 = mode(x)
    mo2 = mode(y)

    if verbose
        make_table(header=["Group" g1 g2];
                   data=["N" length(x) length(y);
                         "mean" mm1 mm2;
                         "std dev" s1 s2;
                         "variance" v1 v2;
                         "median" me1 me2;
                         "mode" mo1 mo2])
    end

    return (mm1=mm1, mm2=mm2, s1=s1, s2=s2, v1=v1, v2=v2, me1=me1, me2=me2, mo1=mo1, mo2=mo2)

end
