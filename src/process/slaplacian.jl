export slaplacian

"""
    slaplacian(obj; m, n, s)

Transform signal channels using surface Laplacian.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `m::Int64=8`: constant positive integer for smoothness
- `n::Int64=8`: Legendre polynomial order
- `s::Float64=10^-5`: smoothing factor

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
- `G::Matrix{Float64}`: transformation matrix
- `H::Matrix{Float64}`: transformation matrix

# Source

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-7
"""
function slaplacian(obj::NeuroAnalyzer.NEURO; m::Int64=4, n::Int64=8, s::Float64=10^-5)

    obj.header.has_locs == false && throw(ArgumentError("Channel locations not available, use load_locs() or add_locs() first."))

    m < 1 && throw(ArgumentError("m must be ≥ 1."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))
    s <= 0 && throw(ArgumentError("s must be > 0."))

    channels = signal_channels(obj)
    locs = obj.locs
    ch_n = nrow(locs)
    ep_n = epoch_n(obj)

    length(channels) > nrow(locs) && throw(ArgumentError("Some channels do not have locations."))

    G = zeros(ch_n, ch_n)
    H = zeros(ch_n, ch_n)
    cosdist = zeros(ch_n, ch_n)

    # r = locs[!, :loc_radius_sph]
    # x = locs[!, :loc_x] ./ maximum(r)
    # y = locs[!, :loc_y] ./ maximum(r)
    # z = locs[!, :loc_z] ./ maximum(r)

    x = locs[!, :loc_x]
    y = locs[!, :loc_y]
    z = locs[!, :loc_z]
    x, y, z = _locnorm(x, y, z)

    # cosine distance
    Threads.@threads for i = 1:ch_n
        @inbounds @simd  for j = 1:ch_n
            cosdist[i, j]  =  1 - (( (x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2 ) / 2 )
        end
    end

    # compute Legendre polynomial
    legpoly = zeros(n, ch_n, ch_n)
    Threads.@threads for idx1 in 1:n
        @inbounds @simd for idx2 in 1:ch_n
            legpoly[idx1, idx2, :] = @views legendre.(cosdist[idx2, :], idx1)
        end
    end

    # compute G and H
    Threads.@threads for i = 1:ch_n
        @inbounds @simd for j = 1:ch_n
            g = 0
            h = 0
            @inbounds @simd for idx_n = 1:n
                g = g + ((2 * idx_n + 1) * legpoly[idx_n, i, j] / ((idx_n * (idx_n + 1))^m))
                h = h - ((2 * (idx_n + 1) * legpoly[idx_n, i, j]) / ((idx_n * (idx_n + 1))^(m - 1)))
            end
            G[i, j] =  g / (4 * pi)
            H[i, j] = -h / (4 * pi)
        end
    end

    # add smoothing factor to the diagonal
    Gs = G + I(ch_n) * s
    GsinvS = sum(inv(Gs))

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        data = @views obj.data[channels, :, ep_idx]
        # dataGs = data[channels, :]' / Gs
        dataGs = Gs / data'
        # C = dataGs .- (sum(dataGs,dims=2)/sum(GsinvS))*GsinvS
        C = data .- (sum(dataGs, dims=2) / sum(GsinvS)) * GsinvS
        # compute surface Laplacian
        obj_new.data[channels, :, ep_idx] = (C'*H)'
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "slaplacian(OBJ, m=m, n=n, s=s)")

    return obj_new, G, H
end

"""
    slaplacian!(obj; m, n, s)

Transform signal channels using surface Laplacian.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `m::Int64=8`: constant positive integer for smoothness
- `n::Int64=8`: Legendre polynomial order
- `s::Float64=10^-5`: smoothing factor

# Returns

- `G::Matrix{Float64}`: transformation matrix
- `H::Matrix{Float64}`: transformation matrix

# Source

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-7
"""
function slaplacian!(obj::NeuroAnalyzer.NEURO; m::Int64=4, n::Int64=8, s::Float64=10^-5)

    obj_tmp, G, H = slaplacian(obj, m=m, n=n, s=s)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return G, H
end