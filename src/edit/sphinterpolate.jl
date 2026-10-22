export sphinterpolate_channel
export sphinterpolate_channel!

"""
    _sph_g(a, b; <keyword arguments>)

Spherical spline matrix G between two sets of electrode positions on the unit sphere (Perrin et al. 1989): `G[i, j] = 1/4π Σₖ (2k + 1) Pₖ(cos θᵢⱼ) / (k(k + 1))ᵐ`, `k = 1:n`, where `θᵢⱼ` is the angle between `a[i, :]` and `b[j, :]` and `Pₖ` is the Legendre polynomial of order `k`.

# Arguments

- `a::AbstractMatrix`: positions (rows), unit vectors, shape (points, 3)
- `b::AbstractMatrix`: positions (rows), unit vectors, shape (points, 3)
- `m::Int64`: spline flexibility constant
- `n::Int64`: number of Legendre terms

# Returns

- `Matrix{Float64}`: G matrix, shape (rows of `a`, rows of `b`)
"""
function _sph_g(a::AbstractMatrix, b::AbstractMatrix; m::Int64, n::Int64)::Matrix{Float64}
    cosang = clamp.(a * b', -1.0, 1.0)
    g = zeros(size(cosang))
    @inbounds for idx in eachindex(cosang)
        x = cosang[idx]
        # Legendre polynomials by the three-term recurrence, starting with P₁
        p0, p1 = 1.0, x
        s = 3 * p1 / 2.0^m
        for k = 2:n
            p0, p1 = p1, ((2k - 1) * x * p1 - (k - 1) * p0) / k
            s += (2k + 1) * p1 / (Float64(k) * (k + 1))^m
        end
        g[idx] = s / (4π)
    end
    return g
end

"""
    sphinterpolate_channel(obj; <keyword arguments>)

Interpolate EEG channel(s) using spherical splines (Perrin et al. 1989).

The signal on the scalp is modelled as a spherical spline fitted to the good EEG channels (bad channels and the interpolated channels are not used) and evaluated at the positions of the interpolated channels. The electrode positions are projected onto the unit sphere. The interpolation is linear: each interpolated channel is a fixed weighted sum of the good channels, computed once and applied to all samples.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"` with channel locations (Cartesian coordinates)
- `ch::Union{String, Vector{String}, Regex}`: channel(s) to interpolate
- `ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj))`: epoch number(s) within to interpolate; default: all epochs
- `seg::Union{Nothing, Tuple{Real, Real}}=nothing`: time segment (from, to) in seconds of the epoch time to interpolate, e.g. a short artifact in a continuous recording; default: the whole epoch(s)
- `m::Int64=4`: spline flexibility constant; must be in `[2, 10]`; higher values give a smoother (stiffer) spline
- `n::Int64=50`: number of Legendre terms; must be ≥ 1
- `lambda::Real=1e-5`: regularization (smoothing) of the spline fit; must be ≥ 0
- `reset_bads::Bool=false`: if `true`, unmark the interpolated channel(s) as bad (`obj.header.recording[:bad_channel]`)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object

# References

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184–187.
"""
function sphinterpolate_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj)),
    seg::Union{Nothing, Tuple{Real, Real}} = nothing,
    m::Int64 = 4,
    n::Int64 = 50,
    lambda::Real = 1e-5,
    reset_bads::Bool = false,
)::NeuroAnalyzer.NEURO
    # validate
    _check_datatype(obj, "eeg")
    _has_locs(obj)
    (2 <= m <= 10) || throw(ArgumentError("m must be in [2, 10]."))
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    lambda >= 0 || throw(ArgumentError("lambda must be ≥ 0."))
    _check_epochs(obj, ep)
    ep = _n2v(ep)
    # samples to interpolate within each epoch
    t_rng = isnothing(seg) ? (1:epoch_len(obj)) : _epoch_seg(obj, seg)

    # channels to interpolate and the good EEG channels used for the spline
    ch_int = get_channel(obj; ch = ch)
    isempty(ch_int) && throw(ArgumentError("No channels selected."))
    signal = get_channel(obj; ch = get_channel(obj; type = "eeg"))
    all(in(signal), ch_int) ||
        throw(ArgumentError("ch must be EEG channels; cannot interpolate other channels."))
    ch_good = setdiff(signal, ch_int, findall(obj.header.recording[:bad_channel]))
    length(ch_good) >= 3 ||
        throw(ArgumentError("At least 3 good EEG channels are needed for interpolation."))

    # positions (matched by label), projected onto the unit sphere
    loc_labels = lowercase.(String.(obj.locs[!, :label]))
    function _xyz(chs::Vector{Int64})::Matrix{Float64}
        idx = [findfirst(==(lowercase(l)), loc_labels) for l in labels(obj)[chs]]
        any(isnothing, idx) &&
            throw(ArgumentError("Locations of all interpolated and good channels are required."))
        p = Matrix{Float64}(obj.locs[Vector{Int64}(idx), [:loc_x, :loc_y, :loc_z]])
        r = sqrt.(sum(abs2, p; dims = 2))
        any(iszero, r) && throw(
            ArgumentError("Some channels have no Cartesian coordinates, use locs_sph2cart() first."),
        )
        return p ./ r
    end
    xyz_good = _xyz(ch_good)
    xyz_int = _xyz(ch_int)

    # spline with a constant term: [G + λI 1; 1ᵀ 0] [c; c0] = [V; 0]; the values at the
    # interpolated positions are [G_int 1] [c; c0], i.e. a fixed linear map of V
    ng = length(ch_good)
    A = [(_sph_g(xyz_good, xyz_good; m = m, n = n) + lambda * I(ng)) ones(ng);
         ones(1, ng) 0.0]
    B = [_sph_g(xyz_int, xyz_good; m = m, n = n) ones(length(ch_int))]
    W = (B / A)[:, 1:ng]

    # create new dataset
    obj_tmp = deepcopy(obj)
    @inbounds for ep_idx in ep
        obj_tmp.data[ch_int, t_rng, ep_idx] = W * @view(obj.data[ch_good, t_rng, ep_idx])
    end

    reset_bads && (obj_tmp.header.recording[:bad_channel][ch_int] .= false)

    push!(
        obj_tmp.history,
        "sphinterpolate_channel(obj; ch=$ch, ep=$ep, seg=$seg, m=$m, n=$n, lambda=$lambda, reset_bads=$reset_bads)",
    )

    return obj_tmp
end

"""
    sphinterpolate_channel!(obj; <keyword arguments>)

Interpolate EEG channel(s) using spherical splines in-place (see `sphinterpolate_channel()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"` with channel locations
- `ch::Union{String, Vector{String}, Regex}`: channel(s) to interpolate
- `ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj))`: epoch number(s) within to interpolate; default: all epochs
- `seg::Union{Nothing, Tuple{Real, Real}}=nothing`: time segment (from, to) in seconds of the epoch time to interpolate; default: the whole epoch(s)
- `m::Int64=4`: spline flexibility constant; must be in `[2, 10]`
- `n::Int64=50`: number of Legendre terms; must be ≥ 1
- `lambda::Real=1e-5`: regularization of the spline fit; must be ≥ 0
- `reset_bads::Bool=false`: if `true`, unmark the interpolated channel(s) as bad (`obj.header.recording[:bad_channel]`)

# Returns

- `Nothing`
"""
function sphinterpolate_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj)),
    seg::Union{Nothing, Tuple{Real, Real}} = nothing,
    m::Int64 = 4,
    n::Int64 = 50,
    lambda::Real = 1e-5,
    reset_bads::Bool = false,
)::Nothing
    obj_tmp = sphinterpolate_channel(
        obj; ch = ch, ep = ep, seg = seg, m = m, n = n, lambda = lambda, reset_bads = reset_bads,
    )
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
