export lrinterpolate_channel
export lrinterpolate_channel!

"""
    lrinterpolate_channel(obj; ch, ep)

Interpolate channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number to interpolate
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function lrinterpolate_channel(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Union{Int64, Vector{Int64}, <:AbstractRange})

    channels = get_channel_bytype(obj, type=obj.header.recording[:data_type])
    @assert ch in channels "channel must be signal channel; cannot interpolate non-signal channels."

    bad_signal = obj.data[:, :, ep]
    good_eps = setdiff(1:nepochs(obj), ep)
    good_chs = setdiff(channels, ch)
    good_signal = _make_epochs(obj.data[:, :, good_eps], ep_n=1)

    # train
    df = @views DataFrame(hcat(good_signal[ch, :, 1], good_signal[good_chs, :, 1]'), :auto)
    train, test = _split(df, 0.80)
    fm = Term(:x1) ~ sum(Term.(Symbol.(names(df[!, Not(:x1)]))))
    linear_regressor = lm(fm, train)
    acc_r2 = r2(linear_regressor)
    prediction = GLM.predict(linear_regressor, test)
    accuracy_testdf = DataFrame(signal_actual = test[!, :x1], signal_predicted = prediction)
    accuracy_testdf.error = accuracy_testdf[!, :signal_actual] 
    accuracy_testdf[!, :signal_predicted]
    acc_mae = mean(abs.(accuracy_testdf.error))
    aic, bic = infcrit(linear_regressor)

    _info("R² for the linear regressor: $(round(acc_r2, digits=3))")
    _info("MAE (test dataset): $(round(acc_mae, digits=3))")
    _info("AIC: $(round(aic, digits=3))")
    _info("BIC: $(round(bic, digits=3))")

    # predict
    obj_new = deepcopy(obj) 
    @inbounds for ep_idx in ep
        df = @views DataFrame(hcat(bad_signal[ch, :, ep_idx], bad_signal[good_chs, :, ep_idx]'), :auto)
        obj_new.data[ch, :, ep_idx] = GLM.predict(linear_regressor, df)
    end

    reset_components!(obj_new)
    push!(obj_new.history, "lrinterpolate_channel(OBJ, ch=$ch, ep=$ep)")

    return obj_new
    
end

"""
    lrinterpolate_channel!(obj; ch, ep)

Interpolate channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number to interpolate
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function lrinterpolate_channel!(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Union{Int64, Vector{Int64}, <:AbstractRange})

    obj_new = lrinterpolate_channel(obj, ch=ch, ep=ep)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
