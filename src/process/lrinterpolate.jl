export lrinterpolate_channel
export lrinterpolate_channel!

"""
    lrinterpolate_channel(obj; <keyword arguments>)

Interpolate channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel to interpolate
- `ep::Int64`: epoch number(s) within to interpolate
- `ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep)`: reference epoch(s), default is all epochs except the interpolated one

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function lrinterpolate_channel(obj::NeuroAnalyzer.NEURO; ch::String, ep::Int64, ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep))::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)[1]
    channels = get_channel(obj, ch=get_channel(obj, type=datatype(obj)))
    @assert length(channels) > 1 "signal must contain > 1 signal channel."
    @assert ch in channels "ch must be a signal channel; cannot interpolate non-signal channels."
    @assert nepochs(obj) > 1 "Training the model requires the signal to have > 1 epoch."

    _check_epochs(obj, ep_ref)
    @assert !(ep in ep_ref) "ep must not be in ep_rep."

    signal_src = @views obj.data[:, :, ep]
    ch_ref = setdiff(channels, ch)
    signal_ref = @views _make_epochs(obj.data[:, :, ep_ref], ep_n=1)

    # train
    df = @views DataFrame(hcat(signal_ref[ch, :, 1], signal_ref[ch_ref, :, 1]'), :auto)
    train, test = _split(df, 0.80)
    fm = Term(:x1) ~ sum(Term.(Symbol.(names(df[!, Not(:x1)]))))
    linear_regressor = GLM.lm(fm, train)
    prediction = GLM.predict(linear_regressor, test)
    accuracy_testdf = DataFrame(signal_actual = test[!, :x1], signal_predicted = prediction)
    accuracy_testdf.error = accuracy_testdf[!, :signal_actual]
    acc_rmse = sqrt(sum((accuracy_testdf.error).^2)) / length(accuracy_testdf.error)
    acc_mae = mean(abs.(accuracy_testdf.error))
    R2, R2adj, aic, bic = infcrit(linear_regressor)

    _info("Accuracy report:")
    _info(" R²: $(round(R2, digits=3))")
    _info(" R² adj: $(round(R2adj, digits=3))")
    _info(" AIC: $(round(aic, digits=3))")
    _info(" BIC: $(round(bic, digits=3))")
    _info(" RMSE: $(round(acc_rmse, digits=3))")
    _info(" MAE: $(round(acc_mae, digits=3))")

    # predict
    obj_new = deepcopy(obj)
    df = @views DataFrame(hcat(signal_src[ch, :], signal_src[ch_ref, :]'), :auto)
    obj_new.data[ch, :, ep] = GLM.predict(linear_regressor, df)

    reset_components!(obj_new)
    push!(obj_new.history, "lrinterpolate_channel(OBJ, ch=$ch, ep=$ep, ep_ref=$ep_ref)")

    return obj_new

end

"""
    lrinterpolate_channel!(obj; <keyword arguments>)

Interpolate channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel to interpolate
- `ep::Int64`: epoch number(s) within to interpolate
- `ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep)`: reference epoch(s), default is all epochs except the interpolated one

# Returns

Nothing
"""
function lrinterpolate_channel!(obj::NeuroAnalyzer.NEURO; ch::String, ep::Int64, ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep))::Nothing

    obj_new = lrinterpolate_channel(obj, ch=ch, ep=ep, ep_ref=ep_ref)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
