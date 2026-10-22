export lrinterpolate_channel
export lrinterpolate_channel!

"""
    lrinterpolate_channel(obj; <keyword arguments>)

Interpolate channel using linear regression: the channel is modelled as a linear combination of the other good signal channels (bad channels are not used), the model is fitted to reference data and used to predict the channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel to interpolate
- `ep::Int64`: epoch index to interpolate (`1` for a continuous signal)
- `ep_ref::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=setdiff(_c(nepochs(obj)), ep)`: reference epochs used for training; default is all epochs except `ep`
- `seg::Union{Nothing, Tuple{Real, Real}}=nothing`: time segment (from, to) in seconds of the epoch time to interpolate; the samples of epoch `ep` outside the segment are added to the training data, so that a segment of a continuous signal can be repaired; default: the whole epoch (requires reference epochs)

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
    seg::Union{Nothing, Tuple{Real, Real}} = nothing,
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch_idx = get_channel(obj; ch = ch)
    isempty(ch_idx) && throw(ArgumentError("No channels selected."))
    ch_idx = ch_idx[1]
    channels = get_channel(obj; ch = get_channel(obj; type = datatype(obj)))
    ch_idx in channels ||
        throw(
            ArgumentError(
                "ch must be a signal channel; cannot interpolate non-signal channels.",
            ),
        )
    _check_epochs(obj, ep)
    ep_ref = _n2v(ep_ref)
    _check_epochs(obj, ep_ref)
    ep in ep_ref && throw(ArgumentError("ep must not be in ep_ref."))
    # predictors: the other good signal channels
    ch_ref = setdiff(channels, ch_idx, findall(obj.header.recording[:bad_channel]))
    isempty(ch_ref) && throw(ArgumentError("No good signal channels to predict ch from."))

    # samples to predict; training data: the reference epochs (all samples) and the
    # samples of epoch ep outside the segment
    if isnothing(seg)
        isempty(ep_ref) && throw(
            ArgumentError("Training the model requires reference epochs (or a segment, seg)."),
        )
        t_rng = 1:epoch_len(obj)
        t_train = Int64[]
    else
        t_rng = _epoch_seg(obj, seg)
        t_train = setdiff(1:epoch_len(obj), t_rng)
    end
    y = vcat(vec(obj.data[ch_idx, :, ep_ref]), obj.data[ch_idx, t_train, ep])
    x = vcat(
        reshape(obj.data[ch_ref, :, ep_ref], length(ch_ref), :)',
        obj.data[ch_ref, t_train, ep]',
    )
    length(y) > length(ch_ref) + 1 ||
        throw(ArgumentError("Not enough training samples for $(length(ch_ref)) predictors."))
    df = DataFrame(hcat(y, x), :auto)

    # linear regression of ch on the other signal channels
    fm = Term(:x1) ~ sum(Term.(Symbol.(names(df)[2:end])))
    linear_regressor = GLM.lm(fm, df)

    # accuracy on the training data
    err = y .- GLM.predict(linear_regressor)
    R2, R2adj, aic, bic = infcrit(linear_regressor)
    _info("Accuracy report (training data):")
    _info(" R²: $(round(R2, digits = 3))")
    _info(" R² adj: $(round(R2adj, digits = 3))")
    _info(" AIC: $(round(aic, digits = 3))")
    _info(" BIC: $(round(bic, digits = 3))")
    _info(" RMSE: $(round(sqrt(mean(abs2, err)), digits = 3))")
    _info(" MAE: $(round(mean(abs, err), digits = 3))")

    # predict ch in epoch ep (within the segment) from the other channels
    obj_tmp = deepcopy(obj)
    df_ep = DataFrame(hcat(obj.data[ch_idx, t_rng, ep], obj.data[ch_ref, t_rng, ep]'), :auto)
    obj_tmp.data[ch_idx, t_rng, ep] = GLM.predict(linear_regressor, df_ep)

    push!(
        obj_tmp.history,
        "lrinterpolate_channel(obj; ch=$ch, ep=$ep, ep_ref=$ep_ref, seg=$seg)",
    )

    return obj_tmp
end

"""
    lrinterpolate_channel!(obj; <keyword arguments>)

Interpolate channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel to interpolate
- `ep::Int64`: epoch index to interpolate (`1` for a continuous signal)
- `ep_ref::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=setdiff(_c(nepochs(obj)), ep)`: reference epochs used for training; default is all epochs except `ep`
- `seg::Union{Nothing, Tuple{Real, Real}}=nothing`: time segment (from, to) in seconds of the epoch time to interpolate; default: the whole epoch

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
    seg::Union{Nothing, Tuple{Real, Real}} = nothing,
)::Nothing
    obj_tmp = lrinterpolate_channel(obj; ch = ch, ep = ep, ep_ref = ep_ref, seg = seg)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
