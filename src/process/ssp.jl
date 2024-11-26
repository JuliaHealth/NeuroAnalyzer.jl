export generate_ssp_projectors
export apply_ssp_projectors

function generate_ssp_projectors(obj::NeuroAnalyzer.NEURO; projectors::Union{Int64, Vector{Int64}}=0)::@NamedTuple{ssp_projectors::Array{Float64, 2}, U::Array{Float64, 2}}

    @assert size(obj.header.recording[:ssp_data], 1) > 0 "OBJ does not have SSP projectors."

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

function apply_ssp_projectors(obj::NeuroAnalyzer.NEURO; projectors::Union{Int64, Vector{Int64}}=0)::NeuroAnalyzer.NEURO

    obj_new = deepcopy(obj)

    # generate
    ssp_projectors, U = generate_ssp_projectors(obj, projectors=projectors)

    # apply
    _info("Applying $(size(U, 2)) SSP projectors")
    obj_new.data[obj.header.recording[:ssp_channels], :, 1] = ssp_projectors * obj.data[obj.header.recording[:ssp_channels], :, 1]

    return obj_new

end