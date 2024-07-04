export od2conc
export od2conc!

"""
    od2conc(obj; ch, ppf)

Convert NIRS optical density (OD) to concentration (HbO, HbR, HbT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_od"))`: index of channels, default is NIRS intensity channels
- `ppf::Vector{Real}=ones(length(obj.header.recording[:wavelengths]))`: Partial path length factors for each wavelength. This is a vector of factors per wavelength. Typical value is ~6 for each wavelength if the absorption change is uniform over the volume of tissue measured. To approximate the partial volume effect of a small localized absorption change within an adult human head, this value could be as small as 0.1. Convention is becoming to set `ppf=1` and to not divide by the source-detector separation such that the resultant "concentration" is in units of Molar mm (or Molar cm if those are the spatial units). This is becoming wide spread in the literature but there is no fixed citation. Use a value of 1 to choose this option.

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function od2conc(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_od"), ppf::Vector{<:Real}=ones(length(obj.header.recording[:wavelengths])))

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])
    _check_datatype(obj, "nirs")
    @assert length(get_channel_bytype(obj, type="nirs_od")) > 0 "OBJ does not contain NIRS OD channels, use intensity2od() first."
    _check_channels(get_channel_bytype(obj, type="nirs_od"), ch)
    @assert length(ppf) == length(obj.header.recording[:wavelengths]) "ppf length does not correspond to the number of wavelengths."

    obj_new = deepcopy(obj)

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

    lst = findall(wl_idx .== 1)
    dc = zeros(3, ep_len, length(lst), ep_n)

    @inbounds for ep_idx in 1:ep_n
        dod = @views obj_new.data[ch, :, ep_idx]

        for idx=eachindex(lst)

            idx1 = lst[idx]
            idx2 = findall(wl_idx .> 1 .&& chp[:, 1] .== chp[idx1, 1] .&& chp[:, 2] .== chp[idx1, 2])

            x = obj.locs[:, :loc_x]
            y = obj.locs[:, :loc_y]

            src_pos = hcat(x, y)'[:, eachindex(unique(chp[:, 1]))]
            det_pos = hcat(x, y)'[:, 1 + (end - length(unique(chp[:, 2]))):end]
            rho = norm(src_pos[:, chp[idx1, 1]] - det_pos[:, chp[idx1, 2]])

            if ppf[1] ≈ 1
                dc[1:2, :, idx, ep_idx] = einv * (dod[vcat(idx1, idx2), :] ./ (ones(ep_len, length(ppf)) .* rho .* ppf')')
            else
                dc[1:2, :, idx, ep_idx] = einv * (dod[vcat(idx1, idx2), :] ./ (ones(ep_len, length(ppf)))')
            end

        end

        dc[3, :, :, ep_idx] = dc[1, :, :, ep_idx] + dc[2, :, :, ep_idx]
        # convert to μM
        dc .*= 10^6
    end

    # add channels
    for dc_idx in 1:size(dc, 3)
        obj_new.data = vcat(obj_new.data, dc[:, :, dc_idx, :])
    end

    # update header
    obj_new.header.recording[:channel_type] = vcat(obj.header.recording[:channel_type], repeat(["nirs_hbo", "nirs_hbr", "nirs_hbt"], size(dc, 3)))
    obj_new.header.recording[:units] = vcat(obj.header.recording[:units], repeat(["μM/mm"], 3 * size(dc, 3)))
    for idx in 1:size(dc, 3)
        obj_new.header.recording[:labels] = vcat(obj_new.header.recording[:labels], ["$(split((obj.header.recording[:labels][idx]), ' ')[1]) HbO", "$(split((obj.header.recording[:labels][idx]), ' ')[1]) HbR", "$(split((obj.header.recording[:labels][idx]), ' ')[1]) HbT"])
    end
    obj_new.header.recording[:labels] = replace.(obj_new.header.recording[:labels], ".0"=>"")

    #=
    for idx in 1:size(dc, 3)
        obj_new.header.recording[:optode_pairs] = vcat(obj_new.header.recording[:optode_pairs], repeat(chp[unique(chp), :][idx, :]', 3))
    end
    obj_new.header.recording[:wavelength_index] = vcat(obj_new.header.recording[:wavelength_index], repeat([-1], 3 * size(dc, 3)))
    =#

    reset_components!(obj_new)
    push!(obj_new.history, "od2conc(OBJ, ch=$ch)")

    return obj_new

end

"""
    od2conc(obj; ch, ppf)

Convert NIRS optical density (OD) to concentration (HbO, HbR, HbT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_od"))`: index of channels, default is NIRS intensity channels
- `ppf::Vector{Real}=ones(length(obj.header.recording[:wavelengths]))`: Partial path length factors for each wavelength. This is a vector of factors per wavelength. Typical value is ~6 for each wavelength if the absorption change is uniform over the volume of tissue measured. To approximate the partial volume effect of a small localized absorption change within an adult human head, this value could be as small as 0.1. Convention is becoming to set `ppf=1` and to not divide by the source-detector separation such that the resultant "concentration" is in units of Molar mm (or Molar cm if those are the spatial units). This is becoming wide spread in the literature but there is no fixed citation. Use a value of 1 to choose this option.
"""
function od2conc!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_od"), ppf::Vector{<:Real}=ones(length(obj.header.recording[:wavelengths])))

    obj_new = od2conc(obj, ch=ch, ppf=ppf)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
