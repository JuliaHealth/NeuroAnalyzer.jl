export spec_seg

"""
    spec_seg(sp, sf, st; t, f)

Return spectrogram segment.

# Arguments

- `sp::Matrix{Float64}`: spectrogram powers
- `st::Vector{Float64}`: spectrogram time
- `sf::Vector{Float64}`: spectrogram frequencies
- `t::Tuple{Real, Real}`: time bounds
- `f::Tuple{Real, Real}`: frequency bounds

# Returns

Named tuple containing:
- `segp::Matrix{Float64}`: powers
- `segs::Shape{Real, Int64}`: shape for plotting
- `tidx::Tuple{Real, Real}`: time indices
- `fidx::Tuple{Real, Real}`: frequency indices
"""
function spec_seg(sp::Matrix{Float64}, st::Vector{Float64}, sf::Vector{Float64}; t::Tuple{Real, Real}, f::Tuple{Real, Real})

    t = tuple_order(t)
    f = tuple_order(f)

    @assert t[1] >= st[1] "t[1] must be ≥ $(st[1])."
    @assert t[2] <= st[end] "t[2] must be ≤ $(st[end])."
    @assert f[1] >= sf[1] "f[1] must be ≥ $(sf[1])."
    @assert f[2] <= sf[end] "f[2] must be ≤ $(sf[end])."

    fidx1 = vsearch(f[1], sf)
    fidx2 = vsearch(f[2], sf)
    tidx1 = vsearch(t[1], st)
    tidx2 = vsearch(t[2], st)

    segp = sp[fidx1:fidx2, tidx1:tidx2]
    segs = Shape([(st[tidx1], sf[fidx1]), (st[tidx2], sf[fidx1]), (st[tidx2], sf[fidx2]), (st[tidx1], sf[fidx2])])

    return (segp=segp, segs=segs, tidx=(tidx1, tidx2), fidx=(fidx1, fidx2))

end

"""
    spec_seg(sp, sf, st; ch, t, f)

Return spectrogram segment.

# Arguments

- `sp::AbstractArray`: spectrogram powers
- `st::AbstractVector`: spectrogram time
- `sf::AbstractVector`: spectrogram frequencies
- `ch::Int64`: channel
- `t::Tuple{Real, Real}`: time bounds
- `f::Tuple{Real, Real}`: frequency bounds

# Returns

Named tuple containing:
- `segp::Array{Float64, 3}`: segment of powers
- `segs::Shape{Real, Int64}`: segment coordinates (shape for plotting)
- `tidx::Tuple{Real, Real}`: time indices
- `fidx::Tuple{Real, Real}`: frequency indices
"""
function spec_seg(sp::AbstractArray, st::AbstractVector, sf::AbstractVector; ch::Int64, t::Tuple{Real, Real}, f::Tuple{Real, Real})

    t = tuple_order(t)
    f = tuple_order(f)

    @assert ch >= 1 "ch must be ≥ 1."
    @assert ch <= size(sp, 3) "ch must be ≤ $(size(sp, 3))."

    @assert t[1] >= st[1] "t[1] must be ≥ $(st[1])."
    @assert t[2] <= st[end] "t[2] must be ≤ $(st[end])."
    @assert f[1] >= sf[1] "f[1] must be ≥ $(sf[1])."
    @assert f[2] <= sf[end] "f[2] must be ≤ $(sf[end])."

    fidx1 = vsearch(f[1], sf)
    fidx2 = vsearch(f[2], sf)
    tidx1 = vsearch(t[1], st)
    tidx2 = vsearch(t[2], st)
    segp = sp[fidx1:fidx2, tidx1:tidx2, ch, :]
    segs = Shape([(st[tidx1], sf[fidx1]), (st[tidx2], sf[fidx1]), (st[tidx2], sf[fidx2]), (st[tidx1], sf[fidx2])])

    return (segp=segp, segs=segs, tidx=(tidx1, tidx2), fidx=(fidx1, fidx2))

end
