export mlinterpolate_channel
export mlinterpolate_channel!

"""
    mlinterpolate_channel(obj; ch, ep, model)

Interpolate channel using a machine-learning model.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel number to interpolate
- `ep::Int64`: epoch number within to interpolate
- `model<:MLJ.Model`: MLJ regressor model

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function mlinterpolate_channel(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Int64, model::T) where {T <: MLJ.Model}

    channels = get_channel_bytype(obj, type=obj.header.recording[:data_type])
    @assert length(channels) > 1 "signal must contain > 1 signal channel."
    @assert ch in channels "ch must be a signal channel; cannot interpolate non-signal channels."
    @assert nepochs(obj) > 1 "Training the model requires the signal to have > 1 epoch."

    good_chs = setdiff(channels, ch)
    good_eps = setdiff(1:nepochs(obj), ep)
    good_signal = _make_epochs(obj.data[:, :, good_eps], ep_n=1)[:, :]

    # train
    y = good_signal[ch, :, 1]
    x = table(good_signal[good_chs, :, 1]')
    mach = MLJ.machine(model, x, y)
    MLJ.fit!(mach)
    yhat = MLJ.predict(mach, x)

    _info("Accuracy report:")
    m = MLJ.RSquared()
    _info(" R²: $(round(m(yhat, y), digits=4))")
    m = MLJ.RootMeanSquaredError()
    _info(" RMSE: $(round(m(yhat, y), digits=4))")

    # predict
    obj_new = deepcopy(obj)
    x = table(obj.data[good_chs, :, ep]')
    obj_new.data[ch, :, ep] = MLJ.predict(mach, x)

    reset_components!(obj_new)
    push!(obj_new.history, "mlinterpolate_channel(OBJ, ch=$ch, ep=$ep, model)")

    return obj_new
    
end

"""
    mlinterpolate_channel!(obj; ch, ep, model)

Interpolate channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel number to interpolate
- `ep::Int64`: epoch number within to interpolate
- `model::T where T <: DataType`: MLJ regressor model

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function mlinterpolate_channel!(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Int64, model::T) where {T <: MLJ.Model}

    obj_new = mlinterpolate_channel(obj, ch=ch, ep=ep, model=model)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
