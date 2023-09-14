export spec_seg

"""
    spec_seg(sp, st, sf; t, f)

Return spectrogram segment.

# Arguments

- `sp::Matrix{Float64}`: spectrogram powers
- `sf::Vector{Float64}`: spectrogram frequencies
- `st::Vector{Float64}`: spectrogram time
- `t::Tuple{Real, Real}`: time bounds
- `f::Tuple{Real, Real}`: frequency bounds

# Returns

Named tuple containing:
- `segp::Matrix{Float64}`: powers
- `segs::Vector{Tuple{Float64, Float64}}`: segment coordinates, for plotting should be converted by `Plots.Shape(segs)`
- `tidx::Tuple{Real, Real}`: time indices
- `fidx::Tuple{Real, Real}`: frequency indices
"""
function spec_seg(sp::Matrix{Float64}, sf::Vector{Float64}, st::Vector{Float64}; t::Tuple{Real, Real}, f::Tuple{Real, Real})

    _check_tuple(t, "t", (st[1], st[end]))
    _check_tuple(f, "f", (sf[1], sf[end]))

    fidx1 = vsearch(f[1], sf)
    fidx2 = vsearch(f[2], sf)
    tidx1 = vsearch(t[1], st)
    tidx2 = vsearch(t[2], st)

    segp = sp[fidx1:fidx2, tidx1:tidx2]
    segs = ([(st[tidx1], sf[fidx1]), (st[tidx2], sf[fidx1]), (st[tidx2], sf[fidx2]), (st[tidx1], sf[fidx2])])

    return (segp=segp, segs=segs, tidx=(tidx1, tidx2), fidx=(fidx1, fidx2))

end

"""
    spec_seg(sp, sf, st; ch, t, f)

Return spectrogram segment.

# Arguments

- `sp::AbstractArray`: spectrogram powers
- `sf::AbstractVector`: spectrogram frequencies
- `st::AbstractVector`: spectrogram time
- `ch::Int64`: channel
- `t::Tuple{Real, Real}`: time bounds
- `f::Tuple{Real, Real}`: frequency bounds

# Returns

Named tuple containing:
- `segp::Array{Float64, 3}`: segment of powers
- `segs::Vector{Tuple{Float64, Float64}}`: segment coordinates, for plotting should be converted by `Plots.Shape(segs)`
- `tidx::Tuple{Real, Real}`: time indices
- `fidx::Tuple{Real, Real}`: frequency indices
"""
function spec_seg(sp::AbstractArray, sf::AbstractVector, st::AbstractVector; ch::Int64, t::Tuple{Real, Real}, f::Tuple{Real, Real})

    _check_tuple(t, "t", (st[1], st[end]))
    _check_tuple(f, "f", (sf[1], sf[end]))
    @assert ch in 1:size(sp, 3) "ch must be in [1, $(size(sp, 3))]."

    fidx1 = vsearch(f[1], sf)
    fidx2 = vsearch(f[2], sf)
    tidx1 = vsearch(t[1], st)
    tidx2 = vsearch(t[2], st)
    segp = sp[fidx1:fidx2, tidx1:tidx2, ch, :]
    segs = ([(st[tidx1], sf[fidx1]), (st[tidx2], sf[fidx1]), (st[tidx2], sf[fidx2]), (st[tidx1], sf[fidx2])])

    return (segp=segp, segs=segs, tidx=(tidx1, tidx2), fidx=(fidx1, fidx2))

end
