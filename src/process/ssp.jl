export generate_ssp_projectors
export apply_ssp_projectors
export apply_ssp_projectors!

"""
    generate_ssp_projectors(obj; <keyword arguments>)

Generate SSP (Signal-Space Projection) projectors from the SSP data embedded in a MEG recording.

The projectors are constructed by extracting the selected projection vectors, re-orthogonalising them via SVD, discarding linearly dependent vectors (relative singular-value threshold 0.01, following MNE-Python), and forming `I - U Uᵀ` to project the data onto the space orthogonal to the noise subspace.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `pidx::Union{Int64, Vector{Int64}}=0`: projection index/indices to use; `0` (default) selects all available projections

# Returns

Named tuple:

- `ssp_projectors::Matrix{Float64}`: projection operator `I − U Uᵀ`
- `U::Matrix{Float64}}`: SVD U matrix (orthonormal basis of the noise subspace)
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
    :ssp_data in keys(obj.header.recording) ||
        throw(ArgumentError("OBJ does not contain SSP projections."))
    n_proj = size(obj.header.recording[:ssp_data], 1)
    n_proj > 0 ||
        throw(ArgumentError("OBJ does not contain SSP projections."))

    # resolve the projection selection
    if pidx isa Int64 && pidx == 0
        # default: use all available projections
        pidx = collect(1:n_proj)
    elseif pidx isa Int64
        # single projection index — validate range then wrap in a vector
        # so the rest of the function works uniformly on Vector{Int64}
        (1 <= pidx <= n_proj) ||
            throw(ArgumentError("pidx must be in [1, $n_proj]."))
        pidx = [pidx]
    else
        # multiple projection indices — sort ascending, then validate range
        pidx = sort(pidx)
        (pidx[1] >= 1 && pidx[end] <= n_proj) ||
            throw(ArgumentError("pidx must be in [1, $n_proj]."))
    end

    # Extract the selected projection vectors.
    # ssp_data is stored as (n_projections × n_channels)
    # transpose to (n_channels × n_selected) so each column is one projection vector
    ssp_vecs = obj.header.recording[:ssp_data][pidx, :]' # n_ch × n_sel

    # re-orthogonalise the projection vectors via SVD
    # U contains orthonormal columns spanning the same subspace as ssp_vecs
    # S contains the singular values ordered largest to smallest
    U, S, _ = svd(ssp_vecs)

    # discard linearly dependent vectors using a relative singular-value
    # threshold of 0.01 (1 % of the largest singular value)
    # this matches the implementation in proj.py of the MNE-Python project
    n_keep = sum(S ./ S[1] .> 0.01)
    U = U[:, 1:n_keep]

    # build the orthogonal projector: I − UUᵀ
    # multiplying a signal by this matrix removes its component in the
    # noise subspace spanned by U (i.e. the SSP projections)
    n_ssp_ch = size(U, 1)
    ssp_projectors = Matrix{Float64}(I(n_ssp_ch)) .- (U * U')

    return (; ssp_projectors, U)
end

"""
    apply_ssp_projectors(obj; <keyword arguments>)

Apply SSP projectors generated from embedded projections to a MEG object.

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

    # create new dataset
    obj_new = deepcopy(obj)

    # generate the projector matrix and the noise-subspace basis U
    ssp_projectors, U = generate_ssp_projectors(obj; pidx = pidx)
    _info("Applying $(size(U, 2)) SSP projection$(_pl(size(U, 2)))")

    ssp_mask = obj.header.recording[:ssp_channels]
    obj_new.data[ssp_mask, :, 1] = ssp_projectors * obj.data[ssp_mask, :, 1]
    push!(obj_new.history, "apply_ssp_projectors(obj, pidx=$pidx)")

    return obj_new
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
    obj_new = apply_ssp_projectors(obj; pidx = pidx)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
