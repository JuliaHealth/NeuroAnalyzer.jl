export generate_ssp_projectors
export apply_ssp_projectors
export apply_ssp_projectors!

# scaling factor from SI units (T, T/m, V) to the units stored in NEURO data
function _ssp_unit_scale(u::String)::Float64
    u in ("T", "T/m", "V") && return 1.0
    u == "fT" && return 1.0e15
    # T/m → fT/cm: 1e15 / 100
    u == "fT/cm" && return 1.0e13
    u == "μV" && return 1.0e6
    u == "nV" && return 1.0e9
    u == "mV" && return 1.0e3
    u == "pT" && return 1.0e12
    throw(ArgumentError("Unsupported unit for SSP channel: $u."))
end

"""
    generate_ssp_projectors(obj; <keyword arguments>)

Generate SSP (Signal-Space Projection) projectors from the SSP data embedded in a MEG recording.

Following MNE-Python (`_make_projector`), the selected projection vectors are restricted to good channels (bad channel entries are set to zero), normalized to unit length, re-orthogonalised via SVD, and linearly dependent vectors are discarded (relative singular-value threshold 0.01). The projector `P = I - U Uᵀ` is built in SI units (T, T/m, V), in which the vectors are defined, and is returned converted to the units of the data: `D P D⁻¹`, where `D` holds the per-channel unit scaling (e.g. 10¹⁵ for fT, 10¹³ for fT/cm). This keeps the projection correct for projectors that span sensors of different types.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `pidx::Union{Int64, Vector{Int64}}=0`: projection index/indices to use; `0` (default) selects all available projections

# Returns

Named tuple:

- `ssp_projectors::Matrix{Float64}`: projection operator to be applied to the data of SSP channels (in data units)
- `U::Matrix{Float64}`: orthonormal basis of the noise subspace (in SI units)
"""
function generate_ssp_projectors(
    obj::NeuroAnalyzer.NEURO;
    pidx::Union{Int64, Vector{Int64}} = 0,
)::@NamedTuple{
    ssp_projectors::Matrix{Float64},
    U::Matrix{Float64},
}
    # validate
    _check_datatype(obj, "meg")
    rec = obj.header.recording
    (:ssp_data in keys(rec) && :ssp_channels in keys(rec)) ||
        throw(ArgumentError("OBJ does not contain SSP projections."))
    ssp_data = rec[:ssp_data]
    ssp_mask = rec[:ssp_channels]
    n_proj = size(ssp_data, 1)
    n_proj > 0 || throw(ArgumentError("OBJ does not contain SSP projections."))
 
    # consistency between SSP data and the current channel list
    length(ssp_mask) == size(obj.data, 1) || throw(
        ArgumentError(
            "SSP channel mask ($(length(ssp_mask))) does not match the number of channels ($(size(obj.data, 1))); channels were probably added or removed after import.",
        ),
    )
    n_ssp_ch = count(ssp_mask)
    size(ssp_data, 2) == n_ssp_ch || throw(
        ArgumentError(
            "SSP vectors have $(size(ssp_data, 2)) elements but $n_ssp_ch SSP channels are marked.",
        ),
    )

    # resolve the projection selection
    if pidx isa Int64
        if pidx == 0
            pidx = collect(1:n_proj)
        else
            (1 <= pidx <= n_proj) || throw(ArgumentError("pidx must be in [1, $n_proj]."))
            pidx = [pidx]
        end
    else
        isempty(pidx) && throw(ArgumentError("pidx must not be empty."))
        pidx = sort(unique(pidx))
        (pidx[1] >= 1 && pidx[end] <= n_proj) ||
            throw(ArgumentError("pidx must be in [1, $n_proj]."))
    end

    # n_ch × n_sel: each column is one projection vector (SI units)
    ssp_vecs = Matrix{Float64}(permutedims(ssp_data[pidx, :]))

    # exclude bad channels: their entries are set to zero (as in MNE-Python)
    bads = rec[:bad_channel][ssp_mask]
    norm_orig = [norm(ssp_vecs[:, k]) for k in axes(ssp_vecs, 2)]
    ssp_vecs[bads, :] .= 0

    # normalize each vector to unit length
    for k in axes(ssp_vecs, 2)
        n = norm(ssp_vecs[:, k])
        n > 0 || continue
        # warn when excluding bad channels removed a sizeable part of the vector
        norm_orig[k] > 0 && n / norm_orig[k] < 0.9 && _warn(
            "Projection vector $(pidx[k]) lost $(round(100 * (1 - n / norm_orig[k]); digits = 1))% of its norm after excluding bad channels.",
        )
        ssp_vecs[:, k] ./= n
    end

    # re-orthogonalise via SVD and discard linearly dependent vectors
    U, S, _ = svd(ssp_vecs)
    (isempty(S) || S[1] == 0) &&
        throw(ArgumentError("All selected projection vectors are zero."))
    n_keep = count(S ./ S[1] .> 0.01)
    n_keep < n_ssp_ch || throw(
        ArgumentError(
            "Number of independent projection vectors ($n_keep) equals the number of channels; the projector would remove all signal.",
        ),
    )
    U = U[:, 1:n_keep]

    # projector in SI units, then converted to data units: D P D⁻¹
    P = Matrix{Float64}(I(n_ssp_ch)) .- (U * U')
    d = [_ssp_unit_scale(u) for u in rec[:unit][ssp_mask]]
    ssp_projectors = d .* P ./ d'

    return (; ssp_projectors, U)
end

"""
    apply_ssp_projectors(obj; <keyword arguments>)

Apply SSP projectors generated from embedded projections to a MEG object. All epochs are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `pidx::Union{Int64, Vector{Int64}}=0`: projection index/indices to use; `0` (default) selects all available projections

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object with SSP projections applied
"""
function apply_ssp_projectors(
    obj::NeuroAnalyzer.NEURO;
    pidx::Union{Int64, Vector{Int64}} = 0,
)::NeuroAnalyzer.NEURO
    # validate
    _check_datatype(obj, "meg")

    # generate the projector matrix and the noise-subspace basis U
    ssp_projectors, U = generate_ssp_projectors(obj; pidx = pidx)
    _info("Applying $(size(U, 2)) SSP projection$(_pl(size(U, 2)))")

    # create new dataset
    obj_tmp = deepcopy(obj)

    # project all samples of all epochs at once: (n_ssp_ch × n_samples·n_epochs)
    ssp_idx = findall(obj.header.recording[:ssp_channels])
    n_s = size(obj.data, 2)
    n_e = size(obj.data, 3)
    x = reshape(obj.data[ssp_idx, :, :], length(ssp_idx), n_s * n_e)
    obj_tmp.data[ssp_idx, :, :] = reshape(ssp_projectors * x, length(ssp_idx), n_s, n_e)

    push!(obj_tmp.history, "apply_ssp_projectors(obj, pidx=$pidx)")

    return obj_tmp
end

"""
    apply_ssp_projectors!(obj; <keyword arguments>)

Apply SSP projectors from embedded projections.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `pidx::Union{Int64, Vector{Int64}}=0`: projection index/indices to use; `0` (default) selects all available projections

# Returns

- `Nothing`
"""
function apply_ssp_projectors!(
    obj::NeuroAnalyzer.NEURO;
    pidx::Union{Int64, Vector{Int64}} = 0,
)::Nothing
    obj_tmp = apply_ssp_projectors(obj; pidx = pidx)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
