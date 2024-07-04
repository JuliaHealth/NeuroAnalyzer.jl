export scale
export scale!

"""
    scale(obj; ch, factor)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `factor::Real`: signal is multiplied by `factor`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function scale(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), factor::Real)

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views obj_new.data[ch, :, :] .* factor
    reset_components!(obj_new)
    push!(obj_new.history, "scale(OBJ, ch=$ch, factor=$factor)")

    return obj_new

end

"""
    scale!(obj; ch, factor)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `factor::Real`: signal is multiplied by `factor`
"""
function scale!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), factor::Real)

    obj_new = scale(obj, ch=ch, factor=factor)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
