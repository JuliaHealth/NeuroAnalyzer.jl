export intensity2od
export intensity2od!

"""
    intensity2od(s)

Convert NIRS raw intensity signal to optical density (OD).

The optical density is defined as:

    OD = −log(I / Ī)

where `Ī` is the mean absolute intensity across the sample (time) dimension. The mean is computed per channel and per epoch so that the reference level adapts to each recording segment independently.

# Arguments

- `s::AbstractArray`: intensity array with axes, shape (channels, samples, epochs)

# Returns

- `od::AbstractArray`: optical density array, same shape as `s`
"""
function intensity2od(s::AbstractArray)::AbstractArray

    # compute the reference level: mean absolute intensity over the sample dimension (dim 2)
    # result shape is (ch × 1 × epochs)
    sm = mean(abs.(s); dims = 2)

    od = -log.(abs.(s) ./ sm)

    return od
end

"""
    intensity2od(obj; <keyword arguments>)

Convert NIRS intensity (RAW) channels to optical density (OD) and append the OD channels to the object. Both the original intensity channels and the new OD channels are retained.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj, type="nirs_int"))`: list of channels, default
  is all NIRS intensity channels

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object with OD channels appended after the selected intensity channels
"""
function intensity2od(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex} = get_channel(obj, type = "nirs_int"),
    )::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)

    # validate
    length(get_channel(obj; type = "nirs_int")) > 0 ||
        throw(ArgumentError("OBJ does not contain NIRS intensity channels."))
    _check_datatype(obj, "nirs")
    _check_channels(get_channel(obj; type = "nirs_int"), ch)

    # create new dataset
    obj_new = deepcopy(obj)

    # index of non-selected channels (kept as-is in the output)
    other_ch = setdiff(axes(obj.data, 1), ch)

    # ------------------------------------------------------------------ #
    # Signal data: [original intensity | new OD | remaining channels]    #
    # ------------------------------------------------------------------ #
    obj_new.data = vcat(
        obj.data[ch, :, :],
        reshape(
            intensity2od(@view(obj.data[ch, :, :])),
            length(ch), epoch_len(obj), nepochs(obj)
        ),
        obj.data[other_ch, :, :],
    )

    # ------------------------------------------------------------------ #
    # header updates: all three sections must cover the same channel     #
    # order as the data array: [ch | OD(ch) | other_ch].                 #
    # ------------------------------------------------------------------ #

    # wavelength index: replicate for OD channels; pass through for others
    obj_new.header.recording[:wavelength_index] = vcat(
        obj.header.recording[:wavelength_index][ch],        # intensity
        obj.header.recording[:wavelength_index][ch],        # OD (same wavelength)
        obj.header.recording[:wavelength_index][other_ch],
    )

    # optode pairs: same fix as wavelength_index
    obj_new.header.recording[:optode_pairs] = vcat(
        obj.header.recording[:optode_pairs][ch, :],   # intensity
        obj.header.recording[:optode_pairs][ch, :],   # OD (same pairs)
        obj.header.recording[:optode_pairs][other_ch, :],
    )

    # channel type: OD channels get the "nirs_od" type string
    obj_new.header.recording[:channel_type] = vcat(
        obj.header.recording[:channel_type][ch],
        repeat(["nirs_od"], length(ch)),
        obj.header.recording[:channel_type][other_ch],
    )

    # channel labels: append " OD" suffix to distinguish from the raw channels
    obj_new.header.recording[:label] = vcat(
        obj.header.recording[:label][ch],              # intensity
        obj.header.recording[:label][ch] .* " OD",     # OD (same pairs)
        obj.header.recording[:label][other_ch],
    )

    # units: OD channels inherit the same unit string as their source channels
    obj_new.header.recording[:unit] = vcat(
        obj.header.recording[:unit][ch],              # intensity
        obj.header.recording[:unit][ch],              # OD (same pairs)
        obj.header.recording[:unit][other_ch],
    )

    # reset bad-channel flags to match the new channel count
    obj_new.header.recording[:bad_channel] = zeros(Bool, size(obj_new.data, 1))

    push!(obj_new.history, "intensity2od(obj; ch=$ch)")

    return obj_new
end

"""
    intensity2od!(obj; <keyword arguments>)

Convert NIRS intensity channels to optical density in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj, type="nirs_int"))`: list of channels, default is NIRS intensity channels

# Returns

- `Nothing`
"""
function intensity2od!(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex} = get_channel(obj, type = "nirs_int"),
    )::Nothing
    obj_new = intensity2od(obj; ch = ch)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing
end
