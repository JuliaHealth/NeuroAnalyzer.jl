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
    ep_ref::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = setdiff(
        _c(nepochs(obj)),
        ep,
    ),
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

    # training data: the reference epochs only (all their samples), flattened
    ep_ref = _n2v(ep_ref)
    ch_ref = setdiff(channels, ch)
    y = vec(obj.data[ch, :, ep_ref])
    x = reshape(obj.data[ch_ref, :, ep_ref], length(ch_ref), :)'
    df = DataFrame(hcat(y, x), :auto)

    # linear regression of ch on the other signal channels
    fm = Term(:x1) ~ sum(Term.(Symbol.(names(df)[2:end])))
    linear_regressor = GLM.lm(fm, df)

    # accuracy on the reference epochs
    err = y .- GLM.predict(linear_regressor)
    R2, R2adj, aic, bic = infcrit(linear_regressor)
    _info("Accuracy report (reference epochs):")
    _info(" R²: $(round(R2, digits = 3))")
    _info(" R² adj: $(round(R2adj, digits = 3))")
    _info(" AIC: $(round(aic, digits = 3))")
    _info(" BIC: $(round(bic, digits = 3))")
    _info(" RMSE: $(round(sqrt(mean(abs2, err)), digits = 3))")
    _info(" MAE: $(round(mean(abs, err), digits = 3))")

    # predict ch in epoch ep from the other channels of that epoch
    obj_tmp = deepcopy(obj)
    df_ep = DataFrame(hcat(obj.data[ch, :, ep], obj.data[ch_ref, :, ep]'), :auto)
    obj_tmp.data[ch, :, ep] = GLM.predict(linear_regressor, df_ep)

    push!(obj_tmp.history, "lrinterpolate_channel(obj; ch=$ch, ep=$ep, ep_ref=$ep_ref)")

    return obj_tmp
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
    ep_ref::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = setdiff(
        _c(nepochs(obj)),
        ep,
    ),
)::Nothing
    obj_tmp = lrinterpolate_channel(obj; ch = ch, ep = ep, ep_ref = ep_ref)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
