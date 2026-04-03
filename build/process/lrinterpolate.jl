export lrinterpolate_channel
export lrinterpolate_channel!

"""
    lrinterpolate_channel(obj; <keyword arguments>)

Interpolate channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel to interpolate
- `ep::Int64`: epoch index to interpolate
- `ep_ref::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=setdiff(_c(nepochs(obj)), ep)`: reference epochs used for training; default is all epochs except `ep`

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function lrinterpolate_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Int64,
    ep_ref::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = setdiff(_c(nepochs(obj)), ep),
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)[1]
    isempty(ch) && throw(ArgumentError("No channels selected."))
    channels = get_channel(obj; ch = get_channel(obj; type = datatype(obj)))
    length(channels) > 1 ||
        throw(ArgumentError("signal must contain > 1 signal channel."))
    ch in channels ||
        throw(
            ArgumentError(
                "ch must be a signal channel; cannot interpolate non-signal channels.",
            ),
        )
    nepochs(obj) > 1 ||
        throw(ArgumentError("Training the model requires the signal to have > 1 epoch."))
    _check_epochs(obj, ep_ref)
    ep in ep_ref && throw(ArgumentError("ep must not be in ep_rep."))

    # source data
    signal_src = @view(obj.data[:, :, ep])

    # reference channels and data
    ch_ref = setdiff(channels, ch)
    signal_ref =
        reshape(obj.data, size(obj.data, 1), (size(obj.data, 2) * size(obj.data, 3)), 1)

    # train
    df = DataFrame(
        hcat(@view(signal_ref[ch, :, 1]), @view(signal_ref[ch_ref, :, 1])'),
        :auto,
    )
    train, test = _split(df, 0.8)
    fm = Term(:x1) ~ sum(Term.(Symbol.(names(df[!, Not(:x1)]))))
    linear_regressor = GLM.lm(fm, train)
    prediction = GLM.predict(linear_regressor, test)
    accuracy_testdf =
        DataFrame(; signal_actual = test[!, :x1], signal_predicted = prediction)
    accuracy_testdf.error = accuracy_testdf[!, :signal_actual]
    acc_rmse = sqrt(sum((accuracy_testdf.error) .^ 2)) / length(accuracy_testdf.error)
    acc_mae = mean(abs.(accuracy_testdf.error))
    R2, R2adj, aic, bic = infcrit(linear_regressor)

    _info("Accuracy report:")
    _info(" R²: $(round(R2, digits = 3))")
    _info(" R² adj: $(round(R2adj, digits = 3))")
    _info(" AIC: $(round(aic, digits = 3))")
    _info(" BIC: $(round(bic, digits = 3))")
    _info(" RMSE: $(round(acc_rmse, digits = 3))")
    _info(" MAE: $(round(acc_mae, digits = 3))")

    # predict

    # create new dataset
    obj_new = deepcopy(obj)

    df = DataFrame(
        hcat(
            @view(signal_src[ch, :]),
            @view(signal_src[ch_ref, :])',
        ),
        :auto,
    )
    obj_new.data[ch, :, ep] = GLM.predict(linear_regressor, df)

    push!(obj_new.history, "lrinterpolate_channel(obj; ch=$ch, ep=$ep, ep_ref=$ep_ref)")

    return obj_new
end

"""
    lrinterpolate_channel!(obj; <keyword arguments>)

Interpolate channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel to interpolate
- `ep::Int64`: epoch index(s) within to interpolate
- `ep_ref::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=setdiff(_c(nepochs(obj)), ep)`: reference epochs used for training; default is all epochs except `ep`

# Returns

- `Nothing`
"""
function lrinterpolate_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Int64,
    ep_ref::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = setdiff(_c(nepochs(obj)), ep),
)::Nothing
    obj_new = lrinterpolate_channel(obj; ch = ch, ep = ep, ep_ref = ep_ref)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
