export spec_seg
export flim
export tlim

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

"""
    flim(p, f; frq_lim)

Trim power spectrum or spectrogram array to a range of frequencies.

# Arguments

- `p::AbstractArray`: powers
- `f::AbstractVector`: frequencies
- `frq_lim::Tuple{Real, Real}`: frequency bounds

# Returns

Named tuple containing:
- `p::AbstractArray`: powers
- `f::AbstractVector`: frequencies
"""
function flim(p::AbstractArray, f::AbstractVector; frq_lim::Tuple{Real, Real})

    @assert ndims(p) in [3, 4] "Input array must have 3 (power spectrum) or 4 (spectrogram) dimensions."

    _check_tuple(frq_lim, "frq_lim", (f[1], f[end]))

    f1_idx = vsearch(frq_lim[1], f)
    f2_idx = vsearch(frq_lim[2], f)
    f_new = @views f[f1_idx:f2_idx]

    if ndims(p) == 3
        # power spectrum
        p_new = @views p[:, f1_idx:f2_idx, :]
    else
        # spectrogram
        p_new = @views p[f1_idx:f2_idx, :, :, :]
    end

    return (p=p_new, f=f_new)

end

"""
    tlim(p, f; frq_lim)

Trim spectrogram array to a range of time points.

# Arguments

- `p::AbstractArray`: powers
- `t::AbstractVector`: time points
- `seg::Tuple{Real, Real}`: time segment

# Returns

Named tuple containing:
- `p::AbstractArray`: powers
- `t::AbstractVector`: time points
"""
function tlim(p::AbstractArray, t::AbstractVector; seg::Tuple{Real, Real})

    @assert ndims(p) == 4 "Input array must have 4 dimensions."

    _check_tuple(seg, "seg", (t[1], t[end]))

    t1_idx = vsearch(seg[1], t)
    t2_idx = vsearch(seg[2], t)
    t_new = @views t[t1_idx:t2_idx]
    p_new = @views p[:, t1_idx:t2_idx, :, :]

    return (p=p_new, t=t_new)

end
