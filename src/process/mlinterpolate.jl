export mlinterpolate_channel
export mlinterpolate_channel!

"""
    mlinterpolate_channel(obj; <keyword arguments>)

Interpolate a single channel in a specified epoch using a machine-learning regression model trained on the remaining epochs.

The model is trained with all other signal channels as features and the target channel as the response. After training, the fitted model predicts the target channel's values for the epoch to interpolate.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: name of the channel to interpolate
- `ep::Int64`: index of the epoch to interpolate
- `ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep)`: reference epochs used
  for training; default is all epochs except `ep`
- `model::T where T <: MLJ.Model`: any MLJ regressor (e.g. `RandomForestRegressor`)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object with the specified channel/epoch replaced by the model's prediction
"""
function mlinterpolate_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Int64,
    ep_ref::Union{Int64, Vector{Int64}, AbstractRange} = setdiff(_c(nepochs(obj)), ep),
    model::T,
)::NeuroAnalyzer.NEURO where {T <: MLJ.Model}

    # resolve channel names to integer indices
    channels = get_channel(obj; type = datatype(obj))
    channel_labels = labels(obj)[channels]

    # validate
    length(channels) > 1 ||
        throw(ArgumentError("Signal must contain > 1 signal channel."))
    ch in channel_labels ||
        throw(
            ArgumentError(
                "\"$ch\" is not a signal channel; cannot interpolate non-signal channels.",
            ),
        )
    nepochs(obj) > 1 ||
        throw(ArgumentError("Training the model requires > 1 epoch."))
    _check_epochs(obj, ep_ref)
    ep in ep_ref && throw(ArgumentError("ep ($ep) must not be included in ep_ref."))

    ch_idx = get_channel(obj; ch = ch)[1]     # String → Int64
    ch_ref = setdiff(channel_indices, ch_idx) # all signal channels except ch

    # ------------------------------------------------------------------ #
    # build training data from reference epochs                          #
    # flatten the selected epochs into a single time axis:               #
    #   (ch × samples × n_ref_epochs) → (ch × samples*n_ref_epochs)      #
    # ------------------------------------------------------------------ #
    # extract only the signal channels from the reference epochs
    ref_data = obj.data[channel_indices, :, ep_ref] # signal ch × samples × ep_ref
    n_ch_sig = length(channel_indices)
    n_samples = epoch_len(obj)
    n_ref = length(ep_ref)

    # reshape to (signal_channels × total_samples)
    ref_flat = reshape(ref_data, n_ch_sig, n_samples * n_ref)

    # find the position of ch_idx within channel_indices (1-based into ref_flat)
    ch_pos = findfirst(==(ch_idx), channel_indices)
    ch_ref_pos = setdiff(1:n_ch_sig, ch_pos)

    # target vector: the channel to interpolate, flattened over reference epochs
    y = ref_flat[ch_pos, :]

    # feature matrix: all other signal channels, transposed so rows = samples
    x = MLJBase.table(ref_flat[ch_ref_pos, :]')

    # ------------------------------------------------------------------ #
    # Train the model                                                    #
    # ------------------------------------------------------------------ #
    mach = MLJ.machine(model, x, y)
    MLJ.fit!(mach)

    # report training accuracy metrics for diagnostics
    yhat = MLJ.predict(mach, x)
    _info("Training accuracy on reference epochs:")
    _info("  R²:   $(round(MLJ.RSquared()(yhat, y); digits = 4))")
    _info("  RMSE: $(round(MLJ.RootMeanSquaredError()(yhat, y); digits = 4))")

    # ------------------------------------------------------------------ #
    # Predict the interpolated epoch                                     #
    # ------------------------------------------------------------------ #
    obj_new = deepcopy(obj)

    # feature matrix for the epoch to reconstruct: same other-channel layout
    x_pred = MLJBase.table(obj.data[ch_ref, :, ep]')
    obj_new.data[ch_idx, :, ep] = MLJ.predict(mach, x_pred)

    push!(
        obj_new.history,
        "mlinterpolate_channel(obj; ch=$ch, ep=$ep, ep_ref=$ep_ref, model=$(typeof(model)))",
    )

    return obj_new
end

"""
    mlinterpolate_channel!(obj; <keyword arguments>)

Interpolate a channel using an MLJ regression model, modifying `obj` in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: name of the channel to interpolate
- `ep::Int64`: index of the epoch to interpolate
- `ep_ref::Union{Int64, Vector{Int64}, AbstractRange}=setdiff(_c(nepochs(obj)), ep)`: reference epochs for training; default is all epochs except `ep`
- `model::T where T <: MLJ.Model`: any MLJ regressor (e.g. `RandomForestRegressor`)

# Returns

- `Nothing`
"""
function mlinterpolate_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Int64,
    ep_ref::Union{Int64, Vector{Int64}, AbstractRange} = setdiff(_c(nepochs(obj)), ep),
    model::T,
)::Nothing where {T <: MLJ.Model}
    obj_new = mlinterpolate_channel(obj; ch = ch, ep = ep, ep_ref = ep_ref, model = model)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
