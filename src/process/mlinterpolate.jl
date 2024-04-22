export mlinterpolate_channel
export mlinterpolate_channel!

"""
    mlinterpolate_channel(obj; ch, ep, ep_ref, model)

Interpolate channel using a machine-learning model.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel number to interpolate
- `ep::Int64`: epoch number within to interpolate
- `ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep)`: reference epoch(s), default is all epochs except the interpolated one
- `model<:MLJ.Model`: MLJ regressor model

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function mlinterpolate_channel(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Int64, ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep), model::T) where {T <: MLJ.Model}

    channels = get_channel_bytype(obj, type=obj.header.recording[:data_type])
    @assert length(channels) > 1 "signal must contain > 1 signal channel."
    @assert ch in channels "ch must be a signal channel; cannot interpolate non-signal channels."
    @assert nepochs(obj) > 1 "Training the model requires the signal to have > 1 epoch."

    _check_epochs(obj, ep_ref)
    @assert !(ep in ep_ref) "ep must not be in ep_rep."

    ch_ref = setdiff(channels, ch)
    signal_ref = _make_epochs(obj.data[channels, :, ep_ref], ep_n=1)[:, :]

    # train
    y = signal_ref[ch, :, 1]
    x = table(signal_ref[ch_ref, :, 1]')
    mach = MLJ.machine(model, x, y)
    MLJ.fit!(mach)
    # reconstruct
    yhat = MLJ.predict(mach, x)

    _info("Accuracy report:")
    m = MLJ.RSquared()
    _info(" R²: $(round(m(yhat, y), digits=4))")
    m = MLJ.RootMeanSquaredError()
    _info(" RMSE: $(round(m(yhat, y), digits=4))")

    # predict
    obj_new = deepcopy(obj)
    x = table(obj.data[ch_ref, :, ep]')
    obj_new.data[ch, :, ep] = MLJ.predict(mach, x)

    reset_components!(obj_new)
    push!(obj_new.history, "mlinterpolate_channel(OBJ, ch=$ch, ep=$ep, ep_ref=$ep_ref, model)")

    return obj_new
    
end

"""
    mlinterpolate_channel!(obj; ch, ep, ep_ref, model)

Interpolate channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel number to interpolate
- `ep::Int64`: epoch number within to interpolate
- `ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep)`: reference epoch(s), default is all epochs except the interpolated one
- `model::T where T <: DataType`: MLJ regressor model

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function mlinterpolate_channel!(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Int64, ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep), model::T) where {T <: MLJ.Model}

    obj_new = mlinterpolate_channel(obj, ch=ch, ep=ep, ep_ref=ep_ref, model=model)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
