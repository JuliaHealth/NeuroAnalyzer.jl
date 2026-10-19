export od2conc
export od2conc!

"""
    od2conc(obj; <keyword arguments>)

Convert NIRS optical density (OD) to concentration (HbO, HbR, HbT).

`ppf` parameter is the partial path length factors for each wavelength. This is a vector of factors per wavelength. Typical value is ~6 for each wavelength if the absorption change is uniform over the volume of tissue measured. To approximate the partial volume effect of a small localized absorption change within an adult human head, this value could be as small as 0.1. Convention is becoming to set `ppf=1` and to not divide by the source-detector separation such that the resultant "concentration" is in units of Molar mm (or Molar cm if those are the spatial units). This is becoming wide spread in the literature but there is no fixed citation. Use a value of 1 to choose this option.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj, type="nirs_od"))`: list of channels, default is NIRS intensity channels
- `ppf::Vector{Real}=ones(length(obj.header.recording[:wavelengths]))`: partial path length factors for each wavelength

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function od2conc(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = get_channel(obj; type = "nirs_od"),
    ppf::Vector{<:Real} = ones(length(obj.header.recording[:wavelengths])),
)::NeuroAnalyzer.NEURO
    # validate
    length(get_channel(obj; type = "nirs_od")) > 0 || throw(
        ArgumentError("OBJ does not contain NIRS OD channels, use intensity2od() first."),
    )
    length(ppf) == length(obj.header.recording[:wavelengths]) ||
        throw(ArgumentError("ppf length does not correspond to the number of wavelengths."))
    _check_datatype(obj, "nirs")
    _check_channels(get_channel(obj; type = "nirs_od"), ch)

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_tmp = deepcopy(obj)

    ep_len = epoch_len(obj)
    ep_n = nepochs(obj)

    wl = obj.header.recording[:wavelengths]
    wl_idx = obj.header.recording[:wavelength_index][ch]
    chp = obj.header.recording[:optode_pairs][ch, :]

    e = zeros(length(wl), 2)
    for w_idx in eachindex(wl)
        e[w_idx, :] = _wl2ext(wl[w_idx])
    end
    e = e ./ 10 # convert from /cm to /mm
    einv = inv(e' * e) * e'

    # one concentration triplet per optode pair, from its first-wavelength channel
    lst = findall(wl_idx .== 1)
    dc = zeros(3, ep_len, length(lst), ep_n)

    # source-detector distances (only used for ppf != 1): optode positions from obj.locs
    use_rho = !all(≈(1), ppf)
    rho = ones(length(lst))
    if use_rho
        rec = obj.header.recording
        for (idx, idx1) in enumerate(lst)
            sl = rec[:src_labels][chp[idx1, 1]]
            dl = rec[:det_labels][chp[idx1, 2]]
            si = findfirst(==(sl), obj.locs[!, :label])
            di = findfirst(==(dl), obj.locs[!, :label])
            (isnothing(si) || isnothing(di)) &&
                throw(ArgumentError("Optode locations of $sl / $dl not found in OBJ locs."))
            rho[idx] = sqrt(
                (obj.locs[si, :loc_x] - obj.locs[di, :loc_x])^2 +
                (obj.locs[si, :loc_y] - obj.locs[di, :loc_y])^2,
            )
        end
        _warn(
            "Source-detector distances are computed from OBJ locs coordinates; concentrations are in μM only if these are in mm.",
        )
    end

    @inbounds for ep_idx = 1:ep_n
        dod = @view(obj_tmp.data[ch, :, ep_idx])
        for (idx, idx1) in enumerate(lst)
            # channels of the same optode pair at the other wavelengths
            idx2 = findall(
                wl_idx .> 1 .&& chp[:, 1] .== chp[idx1, 1] .&& chp[:, 2] .== chp[idx1, 2],
            )
            # (Homer3 hmrR_OD2Conc) ppf = 1: Molar·mm, no division by the distance
            dc[1:2, :, idx, ep_idx] = einv * (dod[vcat(idx1, idx2), :] ./ (rho[idx] .* ppf))
        end
        dc[3, :, :, ep_idx] = dc[1, :, :, ep_idx] + dc[2, :, :, ep_idx]
    end
    # convert to μM
    dc .*= 10^6

    # add channels: HbO, HbR, HbT of each optode pair
    for dc_idx in axes(dc, 3)
        obj_tmp.data = vcat(obj_tmp.data, dc[:, :, dc_idx, :])
    end

    # update header
    r = obj_tmp.header.recording
    r[:channel_type] = vcat(r[:channel_type], repeat(["nirs_hbo", "nirs_hbr", "nirs_hbt"], size(dc, 3)))
    r[:unit] = vcat(r[:unit], repeat([use_rho ? "μM" : "μM·mm"], 3 * size(dc, 3)))
    for idx1 in lst
        l = split(obj.header.recording[:label][ch[idx1]], ' ')[1]
        r[:label] = vcat(r[:label], ["$l HbO", "$l HbR", "$l HbT"])
        # a concentration channel is bad if any of its OD channels is bad
        idx2 = findall(wl_idx .> 1 .&& chp[:, 1] .== chp[idx1, 1] .&& chp[:, 2] .== chp[idx1, 2])
        r[:bad_channel] = vcat(r[:bad_channel], fill(any(obj.header.recording[:bad_channel][ch[vcat(idx1, idx2)]]), 3))
    end
    r[:channel_order] = _sort_channels(r[:channel_type])

    push!(obj_tmp.history, "od2conc(obj; ch=$ch)")

    return obj_tmp
end

"""
    od2conc!(obj; <keyword arguments>)

Convert NIRS optical density (OD) to concentration (HbO, HbR, HbT).

`ppf` parameter is the partial path length factors for each wavelength. This is a vector of factors per wavelength. Typical value is ~6 for each wavelength if the absorption change is uniform over the volume of tissue measured. To approximate the partial volume effect of a small localized absorption change within an adult human head, this value could be as small as 0.1. Convention is becoming to set `ppf=1` and to not divide by the source-detector separation such that the resultant "concentration" is in units of Molar mm (or Molar cm if those are the spatial units). This is becoming wide spread in the literature but there is no fixed citation. Use a value of 1 to choose this option.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj, type="nirs_od"))`: list of channels, default is NIRS intensity channels
- `ppf::Vector{Real}=ones(length(obj.header.recording[:wavelengths]))`: partial path length factors for each wavelength

# Returns

- `Nothing`
"""
function od2conc!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = get_channel(obj; type = "nirs_od"),
    ppf::Vector{<:Real} = ones(length(obj.header.recording[:wavelengths])),
)::Nothing
    obj_tmp = od2conc(obj; ch = ch, ppf = ppf)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
