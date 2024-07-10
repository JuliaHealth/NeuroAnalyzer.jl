export csd
export gh

"""
    csd(obj; <keyword arguments>)

Transform data using Current Source Density (CSD) transformation based on spherical spline surface Laplacian.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `m::Int64=4`: positive integer constant that affects spherical spline flexibility, higher `m` values mean increased rigidity
- `n::Int64=8`: Legendre polynomial order
- `lambda::Float64=10^-5`: smoothing factor

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: with `csd` channel types and `µV/m²` units

# Source

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-187
Kayser J, Tenke CE. Principal components analysis of Laplacian waveforms as a generic method for identifying ERP generator patterns: I. Evaluation with auditory oddball tasks. Clin Neurophysiol 2006;117(2):348-368
"""
function csd(obj::NeuroAnalyzer.NEURO; m::Int64=4, n::Int64=8, lambda::Float64=10^-5)

    _check_datatype(obj, "eeg")
    @assert _has_locs(obj) "Channel locations not available, use load_locs() or add_locs() first."
    @assert !(m < 2 || m > 10) "m must be in [2, 10]."
    @assert n >= 1 "n must be ≥ 1."
    @assert lambda > 0 "lambda must be > 0."

    chs = signal_channels(obj)
    locs = obj.locs
    ch_n = nrow(locs)
    ep_n = nepochs(obj)

    @assert length(chs) <= nrow(locs) "Some channels do not have locations."

    G, H = gh(locs, m=m, n=n)

    # add smoothing factor to the diagonal
    Gs = G + I(ch_n) * lambda

    # calculate inverse
    Gs_inv = inv(Gs)

    # sum per row
    Gs_rs = sum(Gs_inv, dims=2)
    Gs_inv_sum = sum(Gs_rs)

    obj_new = deepcopy(obj)
    @inbounds for ep_idx in 1:ep_n
        data = @views obj.data[chs, :, ep_idx]
        # dataGs = data[chs, :]' / Gs
        dataGs = Gs / data'
        # C = dataGs .- (sum(dataGs,dims=2)/sum(Gs_inv_sum))*Gs_inv_sum
        C = data .- (sum(dataGs, dims=2) / sum(Gs_inv_sum)) * Gs_inv_sum
        # compute surface Laplacian
        obj_new.data[chs, :, ep_idx] = (C'*H)'
    end

    obj_new.header.recording[:data_type] = "csd"
    obj_new.header.recording[:channel_type][chs] .= "csd"
    obj_new.header.recording[:units][chs] .= "µV/m²"

    reset_components!(obj_new)
    push!(obj_new.history, "csd(OBJ, m=$m, n=$n, lambda=$lambda)")

    return obj_new

end

"""
    csd!(obj; <keyword arguments>)

Transform data using Current Source Density (CSD) transformation based on spherical spline surface Laplacian.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `m::Int64=4`: positive integer constant that affects spherical spline flexibility, higher `m` values mean increased rigidity
- `n::Int64=8`: Legendre polynomial order
- `lambda::Float64=10^-5`: smoothing factor

# Returns

- `G::Matrix{Float64}`: transformation matrix (SP spline)
- `H::Matrix{Float64}`: transformation matrix (CSD spline)

# Source

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-7
"""
function csd!(obj::NeuroAnalyzer.NEURO; m::Int64=4, n::Int64=8, lambda::Float64=10^-5)

    obj_new = csd(obj, m=m, n=n, lambda=lambda)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    gh(locs; <keyword arguments>)

Generate G and H matrices.

# Arguments

- `locs::DataFrame`
- `m::Int64=4`: positive integer constant that affects spherical spline flexibility, higher `m` values mean increased rigidity
- `n::Int64=8`: Legendre polynomial order

# Returns

Named tuple containing:
- `G::Matrix{Float64}`: transformation matrix (SP spline)
- `H::Matrix{Float64}`: transformation matrix (CSD spline)

# Source

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-7
"""
function gh(locs::DataFrame; m::Int64=4, n::Int64=8)

    @assert !(m < 2 || m > 10) "m must be in [2, 10]."
    @assert n >= 1 "n must be ≥ 1."

    ch_n = nrow(locs)

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
    x, y, z = _locs_norm(x, y, z)

    # compute all cosine distances
    Threads.@threads for i = 1:ch_n
        @inbounds  for j = 1:ch_n
            cosdist[i, j]  =  1 - (( (x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2 ) / 2 )
        end
    end

    # compute Legendre polynomial
    legpoly = zeros(n, ch_n, ch_n)
    Threads.@threads for idx1 in 1:n
        @inbounds for idx2 in 1:ch_n
            legpoly[idx1, idx2, :] = @views legendre.(cosdist[idx2, :], idx1)
        end
    end

    # compute G and H
    Threads.@threads for i = 1:ch_n
        @inbounds for j = 1:ch_n
            g = 0
            h = 0
            @inbounds for idx_n = 1:n
                g = g + ((2.0 * idx_n + 1.0) * legpoly[idx_n, i, j] / ((idx_n * (idx_n + idx_n))^m))
                h = h - ((2.0 * (idx_n + 1.0) * legpoly[idx_n, i, j]) / ((idx_n * (idx_n + idx_n))^(m - 1)))
            end
            G[i, j] =  g / 4 / pi
            H[i, j] = -h / 4 / pi
        end
    end

    return (G=G, H=H)

end
