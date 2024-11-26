export generate_ssp_projectors
export apply_ssp_projectors
export apply_ssp_projectors!

"""
    generate_ssp_projectors(obj; <keyword arguments>)

Generate SSP projectors from embedded projections.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `proj::Union{Int64, Vector{Int64}}=0`: list of projections used for generating projectors, by default use all available projections

# Returns

Named tuple containing:
- `ssp_projectors::Matrix{Float64}` : projectors
- `U::Matrix{Float64}}`: SVD U orthogonal matrix
"""
function generate_ssp_projectors(obj::NeuroAnalyzer.NEURO; proj::Union{Int64, Vector{Int64}}=0)::@NamedTuple{ssp_projectors::Matrix{Float64}, U::Matrix{Float64}}

    @assert "ssp_data" in keys(obj.header.recording) "OBJ does not contain SSP projections."
    @assert size(obj.header.recording[:ssp_data], 1) > 0 "OBJ does not contain SSP projections."

    # by default use all available projections
    if proj == 0
        proj = 1:size(obj.header.recording[:ssp_data], 1)
    end

    if isa(proj, Int64)
        @assert proj >= 1 && proj <= size(obj.header.recording[:ssp_data], 1) "proj must be in [1, $(size(obj.header.recording[:ssp_data], 1))]."
    else
        proj = sort(proj)
        @assert proj[1] >= 1 && proj[end] <= size(obj.header.recording[:ssp_data], 1) "proj must be in [1, $(size(obj.header.recording[:ssp_data], 1))]."
    end

    # extract projections
    ssp_projectors = obj.header.recording[:ssp_data][proj, :]'

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
- `proj::Union{Int64, Vector{Int64}}=0`: list of projections used for generating projectors, by default use all available projections

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function apply_ssp_projectors(obj::NeuroAnalyzer.NEURO; proj::Union{Int64, Vector{Int64}}=0)::NeuroAnalyzer.NEURO

    obj_new = deepcopy(obj)

    # generate
    ssp_projectors, U = generate_ssp_projectors(obj, proj=proj)

    # apply
    _info("Applying $(size(U, 2)) SSP projections")
    obj_new.data[obj.header.recording[:ssp_channels], :, 1] = ssp_projectors * obj.data[obj.header.recording[:ssp_channels], :, 1]

    reset_components!(obj_new)
    push!(obj_new.history, "apply_ssp_projectors(OBJ, proj=$proj)")

    return obj_new

end

"""
    apply_ssp_projectors!(obj; <keyword arguments>)

Apply SSP projectors from embedded projections.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `proj::Union{Int64, Vector{Int64}}=0`: list of projections used for generating projectors, by default use all available projections

# Returns

Nothing
"""
function apply_ssp_projectors!(obj::NeuroAnalyzer.NEURO; proj::Union{Int64, Vector{Int64}}=0)::Nothing

    obj_new = apply_ssp_projectors(obj, proj=proj)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return obj_new

end
