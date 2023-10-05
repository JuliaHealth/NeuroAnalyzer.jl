export intensity2od
export intensity2od!

"""
    intensity2od(s)

Convert NIRS intensity (RAW data) to optical density (OD).

# Arguments

- `s::AbstractArray`

# Returns

- `od::AbstractArray`
"""
function intensity2od(s::AbstractArray)

    sm = mean(abs.(s), dims=2)
    od = -log.(abs.(s) ./ (ones(size(s)) .* sm))

    return od

end

"""
    intensity2od(obj; ch)

Convert NIRS intensity (RAW data) to optical density (OD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_int"))`: index of channels, default is NIRS intensity channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function intensity2od(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_int"))

    _check_channels(obj, ch)
    _check_datatype(obj, "nirs")
    @assert length(get_channel_bytype(obj, type="nirs_int")) > 0 "OBJ does not contain NIRS intensity channels."
    _check_channels(get_channel_bytype(obj, type="nirs_int"), ch)

    obj_new = deepcopy(obj)
    
    # add channels
    obj_new.data = vcat(obj.data[ch, :, :], reshape(intensity2od(obj.data[ch, :, :]), length(ch), epoch_len(obj), nepochs(obj)), obj.data[setdiff(collect(1:size(obj.data, 1)), ch), :, :])
    
    # update header
    obj_new.header.recording[:wavelength_index] = vcat(obj.header.recording[:wavelength_index][ch], obj.header.recording[:wavelength_index][ch])
    obj_new.header.recording[:channel_pairs] = vcat(obj.header.recording[:channel_pairs][ch, :], obj.header.recording[:channel_pairs][ch, :])
    obj_new.header.recording[:channel_type] = vcat(obj.header.recording[:channel_type][ch], repeat(["nirs_od"], length(ch)), obj.header.recording[:channel_type][setdiff(collect(1:size(obj.data, 1)), ch)])
    obj_new.header.recording[:labels] = vcat(obj.header.recording[:labels][ch], obj.header.recording[:labels][ch], obj.header.recording[:labels][setdiff(collect(1:size(obj.data, 1)), ch)])
    obj_new.header.recording[:units] = vcat(obj.header.recording[:units][ch], obj.header.recording[:units][ch], obj.header.recording[:units][setdiff(collect(1:size(obj.data, 1)), ch)])

    reset_components!(obj_new)
    push!(obj_new.history, "intensity2od(OBJ, ch=$ch)")

    return obj_new

end

"""
    intensity2od!(obj; ch)

Convert NIRS intensity (RAW data) to optical density (OD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_int"))`: index of channels, default is NIRS intensity channels
"""
function intensity2od!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_int"))

    obj_new = intensity2od(obj, ch=ch)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
