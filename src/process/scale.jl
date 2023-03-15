export scale
export scale!

"""
    scale(obj; ch, factor)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `factor::Real`: signal is multiplied by `factor`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function scale(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), factor::Real)

    _check_channels(obj, ch)
    
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views obj_new.data[ch, :, :] .* factor
    reset_components!(obj_new)
    push!(obj_new.header.history, "scale(OBJ, ch=$ch, factor=$factor)")

    return obj_new

end

"""
    scale!(obj; ch, factor)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `factor::Real`: signal is multiplied by `factor`
"""
function scale!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), factor::Real)

    obj_tmp = scale(obj, ch=ch, factor=factor)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end
