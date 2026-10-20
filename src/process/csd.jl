export csd
export csd!
export gh

"""
    csd(obj; <keyword arguments>)

Transform EEG data using the Current Source Density (CSD) transformation based on a spherical spline surface Laplacian.

The algorithm (Perrin et al. 1989):

1. Compute the G (potential) and H (Laplacian) spline matrices from electrode locations.
2. Add a Tikhonov regularisation term `λI` to G.
3. Solve for the constrained spline coefficients C.
4. Compute the surface Laplacian as `H × C`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"` with channel locations defined
- `m::Int64=4`: spline flexibility constant; must be in `[2, 10]`; higher values increase rigidity
- `n::Int64=8`: Legendre polynomial order; must be ≥ 1
- `lambda::Float64=10^-5`: Tikhonov smoothing factor; must be > 0

# Returns

- `NeuroAnalyzer.NEURO`: new object with CSD-transformed data, channel types set to `"csd"`, and units set to `"μV/m²"`

# References

1. Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184–187.
2. Kayser J, Tenke CE. Principal components analysis of Laplacian waveforms as a generic method for identifying ERP generator patterns: I. Evaluation with auditory oddball tasks. Clinical Neurophysiology. 2006;117(2):348–368.
"""
function csd(
    obj::NeuroAnalyzer.NEURO;
    m::Int64 = 4,
    n::Int64 = 8,
    lambda::Float64 = 10^-5,
)::NeuroAnalyzer.NEURO
    # validate
    _check_datatype(obj, "eeg")
    _has_locs(obj)
    (m >= 2 && m <= 10) || throw(ArgumentError("m must be in [2, 10]."))
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    lambda > 0 || throw(ArgumentError("lambda must be > 0."))

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = get_channel(obj; type = datatype(obj)))
    isempty(ch) && throw(ArgumentError("No channels selected."))
    _check_ch_locs(ch, labels(obj), obj.locs.label)
    # locations in the order of the channels
    locs = obj.locs[[findfirst(==(l), obj.locs.label) for l in labels(obj)[ch]], :]

    # number of channels
    ch_n = DataFrames.nrow(locs)
    # number of epochs
    ep_n = nepochs(obj)

    G, H = gh(locs; m = m, n = n)

    # regularized G matrix and its inverse
    Gs = G + I(ch_n) * lambda
    Gs_inv = inv(Gs)

    # row sums of the inverse and their total - used for the zero-sum constraint
    Gs_rs = vec(sum(Gs_inv; dims = 2))
    Gs_inv_sum = sum(Gs_rs)

    # create new dataset
    obj_tmp = deepcopy(obj)

    @inbounds for ep_idx = 1:ep_n
        # data: (ch_n × samples)
        data = @view obj.data[ch, :, ep_idx]

        # solve Gs * C_unconstrained = data for each time point
        # (samples × ch_n): spline coefficients (unconstrained)
        dataGs = (Gs \ data)'
        # enforce the zero-sum constraint on the coefficients (Gs C + c0 = V, Σ C = 0):
        # c0 = 1ᵀGs⁻¹V / 1ᵀGs⁻¹1, C = Gs⁻¹V - c0 Gs⁻¹1
        # C[t, i] = dataGs[t, i] - (Σ_j dataGs[t, j]) / Gs_inv_sum * Gs_rs[i]
        # (samples,)
        correction = vec(sum(dataGs; dims = 2)) / Gs_inv_sum
        # (samples × ch_n)
        C = dataGs .- correction .* Gs_rs'

        # surface Laplacian via H matrix
        # back to (ch_n × samples)
        obj_tmp.data[ch, :, ep_idx] = (C * H')'
    end

    obj_tmp.header.recording[:data_type] = "csd"
    obj_tmp.header.recording[:channel_type][ch] .= "csd"
    obj_tmp.header.recording[:unit][ch] .= "μV/m²"

    push!(obj_tmp.history, "csd(obj, m=$m, n=$n, lambda=$lambda)")

    return obj_tmp
end

"""
    csd!(obj; <keyword arguments>)

Transform EEG data using the CSD transformation in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"` with channel locations defined; modified in-place
- `m::Int64=4`: spline flexibility constant; must be in `[2, 10]`; higher values increase rigidity
- `n::Int64=8`: Legendre polynomial order; must be ≥ 1
- `lambda::Float64=10^-5`: Tikhonov smoothing factor; must be > 0

# Returns

- `Nothing`

# References

1. Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184–187.
2. Kayser J, Tenke CE. Principal components analysis of Laplacian waveforms as a generic method for identifying ERP generator patterns: I. Evaluation with auditory oddball tasks. Clinical Neurophysiology. 2006;117(2):348–368.
"""
function csd!(
    obj::NeuroAnalyzer.NEURO;
    m::Int64 = 4,
    n::Int64 = 8,
    lambda::Float64 = 10^-5,
)::Nothing
    obj_tmp = csd(obj; m = m, n = n, lambda = lambda)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end

"""
    gh(locs; <keyword arguments>)

Compute the G (potential spline) and H (CSD spline) matrices for a set of electrode locations on the unit sphere.

# Arguments

- `locs::DataFrame`: electrode locations; must contain columns `:loc_x`, `:loc_y`, `:loc_z` (Cartesian coordinates, normalised to the unit sphere)
- `m::Int64=4`: spline flexibility constant; must be in `[2, 10]`; higher values increase rigidity
- `n::Int64=8`: Legendre polynomial order; must be ≥ 1

# Returns

Named tuple:

- `G::Matrix{Float64}`: potential spline matrix `(ch_n, ch_n)`
- `H::Matrix{Float64}`: CSD spline matrix `(ch_n, ch_n)`

# References

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184–187.
"""
function gh(
    locs::DataFrame;
    m::Int64 = 4,
    n::Int64 = 8,
)::@NamedTuple{
    G::Matrix{Float64},
    H::Matrix{Float64},
}
    (m >= 2 && m <= 10) || throw(ArgumentError("m must be in [2, 10]."))
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))

    ch_n = DataFrames.nrow(locs)

    G = zeros(ch_n, ch_n)
    H = zeros(ch_n, ch_n)

    # project the electrodes onto the unit sphere (cos(angle) = 1 - |ri - rj|² / 2 holds for unit vectors only)
    r = sqrt.(locs.loc_x .^ 2 .+ locs.loc_y .^ 2 .+ locs.loc_z .^ 2)
    any(iszero, r) && throw(ArgumentError("Electrode location at the origin; CSD requires 3-D locations."))
    x, y, z = locs.loc_x ./ r, locs.loc_y ./ r, locs.loc_z ./ r

    # --- cosine distances between all electrode pairs ---
    cosdist = zeros(ch_n, ch_n)
    # thread over rows; each thread writes to a unique row - no contention
    Threads.@threads :static for i = 1:ch_n
        @inbounds for j = 1:ch_n
            cosdist[i, j] = 1 - ((x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2) / 2
        end
    end

    # --- Legendre polynomials for each order up to n ---
    legpoly = zeros(n, ch_n, ch_n)
    Threads.@threads :static for idx1 = 1:n
        @inbounds for idx2 = 1:ch_n
            legpoly[idx1, idx2, :] = legendre.(cosdist[idx2, :], idx1)
        end
    end

    # --- G and H matrices from the spline formula (Perrin et al. 1989, eq. 3–4) ---
    G = zeros(ch_n, ch_n)
    H = zeros(ch_n, ch_n)
    Threads.@threads :static for i = 1:ch_n
        @inbounds for j = 1:ch_n
            g = 0.0
            h = 0.0
            for k = 1:n
                denom_g = (k * (k + 1))^m
                denom_h = (k * (k + 1))^(m - 1)
                lp = legpoly[k, i, j]
                g += (2k + 1) * lp / denom_g
                h -= (2k + 1) * lp / denom_h
            end
            G[i, j] = g / (4π)
            H[i, j] = -h / (4π)
        end
    end

    return (; G, H)
end
