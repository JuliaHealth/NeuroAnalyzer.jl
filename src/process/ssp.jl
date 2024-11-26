export generate_ssp_projectors
export apply_ssp_projectors
export apply_ssp_projectors!

"""
    generate_ssp_projectors(obj; <keyword arguments>)

Generate SSP projectors from embedded projections.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `projectors::Union{Int64, Vector{Int64}}=0`: list of projections used for generating projectors

# Returns

Named tuple containing:
- `ssp_projectors::Matrix{Float64}` : projectors
- `U::Matrix{Float64}}`: SVD U orthogonal matrix
"""
function generate_ssp_projectors(obj::NeuroAnalyzer.NEURO; projectors::Union{Int64, Vector{Int64}}=0)::@NamedTuple{ssp_projectors::Matrix{Float64}, U::Matrix{Float64}}

    @assert "ssp_data" in keys(obj.header.recording) "OBJ does not contain SSP projectors."
    @assert size(obj.header.recording[:ssp_data], 1) > 0 "OBJ does not contain SSP projectors."

    if projectors == 0
        projectors = 1:size(obj.header.recording[:ssp_data], 1)
    end

    if isa(projectors, Int64)
        @assert projectors >= 1 && projectors <= size(obj.header.recording[:ssp_data], 1) "projectors must be in [1, $(size(obj.header.recording[:ssp_data], 1))]."
    else
        projectors = sort(projectors)
        @assert projectors[1] >= 1 && projectors[end] <= size(obj.header.recording[:ssp_data], 1) "projectors must be in [1, $(size(obj.header.recording[:ssp_data], 1))]."
    end

    # extract projectors
    ssp_projectors = obj.header.recording[:ssp_data][projectors, :]'

    # reorthogonalize the vectors
    U, S, _ = svd(ssp_projectors)

    # remove linearly dependent vectors -- this code comes from proj.py of the mne-python project 
    nproj = sum(S ./ S[1] .> 0.01)
    U = U[:, 1:nproj]

    # create projectors
    ssp_projectors = I(count(obj.header.recording[:ssp_channels])) .- (U * U')

    return (ssp_projectors=ssp_projectors, U=U)

end


"""
    apply_ssp_projectors(obj; <keyword arguments>)

Apply SSP projectors from embedded projections.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `projectors::Union{Int64, Vector{Int64}}=0`: list of projections used for generating projectors

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function apply_ssp_projectors(obj::NeuroAnalyzer.NEURO; projectors::Union{Int64, Vector{Int64}}=0)::NeuroAnalyzer.NEURO

    obj_new = deepcopy(obj)

    # generate
    ssp_projectors, U = generate_ssp_projectors(obj, projectors=projectors)

    # apply
    _info("Applying $(size(U, 2)) SSP projectors")
    obj_new.data[obj.header.recording[:ssp_channels], :, 1] = ssp_projectors * obj.data[obj.header.recording[:ssp_channels], :, 1]

    reset_components!(obj_new)
    push!(obj_new.history, "apply_ssp_projectors(OBJ, projectors=$projectors)")

    return obj_new

end

"""
    apply_ssp_projectors!(obj; <keyword arguments>)

Apply SSP projectors from embedded projections.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `projectors::Union{Int64, Vector{Int64}}=0`: list of projections used for generating projectors

# Returns

Nothing
"""
function apply_ssp_projectors!(obj::NeuroAnalyzer.NEURO; projectors::Union{Int64, Vector{Int64}}=0)::Nothing

    obj_new = apply_ssp_projectors(obj, projectors=projectors)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return obj_new

end
