export od2conc
export od2conc!
export nirs_dpf

"""
    nirs_dpf(; <keyword arguments>)

Calculate the differential pathlength factor (DPF) of the adult and child head for a given age and wavelength, using the general equation of Scholkmann & Wolf (2013):

`DPF(λ, A) = 223.3 + 0.05624 · A^0.8493 − 5.723 × 10⁻⁷ · λ³ + 0.001245 · λ² − 0.9025 · λ`

The equation was derived for ages 0–72 years and wavelengths 690–832 nm; values outside these ranges are extrapolated (with a warning).

# Arguments

- `age::Real`: age in years; must be ≥ 0
- `wl::Union{Real, AbstractVector{<:Real}}`: wavelength(s) in nm

# Returns

- `Union{Float64, Vector{Float64}}`: DPF value(s)

# References

1. Scholkmann F, Wolf M. General equation for the differential pathlength factor of the frontal human head depending on wavelength and age. Journal of Biomedical Optics. 2013;18(10):105004.
"""
function nirs_dpf(;
    age::Real,
    wl::Union{Real, AbstractVector{<:Real}},
)::Union{Float64, Vector{Float64}}
    # validate
    age >= 0 || throw(ArgumentError("age must be ≥ 0."))
    all(>(0), wl) || throw(ArgumentError("wl must be > 0."))
    age > 72 && _warn("DPF equation was derived for ages 0-72 years; age = $age is extrapolated.")
    any(w -> w < 690 || w > 832, wl) &&
        _warn("DPF equation was derived for wavelengths 690-832 nm; other wavelengths are extrapolated.")

    f(λ) = 223.3 + 0.05624 * age^0.8493 - 5.723e-7 * λ^3 + 0.001245 * λ^2 - 0.9025 * λ

    return wl isa Real ? f(wl) : f.(wl)
end

# scale factor converting a length in `unit` to mm; `nothing` if the unit is unknown
function _length_to_mm(unit::String)::Union{Nothing, Float64}
    u = lowercase(strip(unit))
    u == "mm" && return 1.0
    u == "cm" && return 10.0
    u == "m" && return 1000.0
    return nothing
end

"""
    od2conc(obj; <keyword arguments>)

Convert NIRS optical density (OD) to concentration (HbO, HbR, HbT) using the modified Beer-Lambert law.

For each optode pair, the OD changes at all wavelengths are converted to the concentration changes of oxygenated (HbO) and deoxygenated (HbR) hemoglobin (least-squares solution using the extinction coefficients of Wray et al. 1988); HbT = HbO + HbR.

The path length is set using one of:

- default (`ppf` and `age` not given): partial path length factor 1 and no division by the source-detector distance, so the concentrations are in μM·mm (the convention of Homer3)
- `age`: the differential pathlength factor (DPF) of each wavelength is calculated for the given age (see `nirs_dpf()`), and the OD is divided by DPF × source-detector distance: concentrations in μM
- `ppf`: partial path length factors given for each wavelength (e.g. ~6 for a uniform absorption change, as small as 0.1 for a small localized change); with `ppf ≠ 1` the OD is divided by `ppf` × source-detector distance: concentrations in μM

The source-detector distances are calculated from the 3-D optode positions in `obj.locs` and converted to mm using the length unit of the recording (`obj.header.recording[:length_unit]`, read from SNIRF and Homer `.nirs` files); if the unit is unknown, mm are assumed (with a warning).

The new channels (types `nirs_hbo`, `nirs_hbr`, `nirs_hbt`) are inserted after the optical channels; they keep the optode pair of their OD channels (wavelength index 0).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj, type="nirs_od"))`: list of channels, default is all NIRS OD channels
- `ppf::Union{Nothing, Vector{<:Real}}=nothing`: partial path length factors for each wavelength
- `age::Union{Nothing, Real}=nothing`: age of the subject in years, used to calculate the DPF; cannot be used together with `ppf`

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function od2conc(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = get_channel(obj; type = "nirs_od"),
    ppf::Union{Nothing, Vector{<:Real}} = nothing,
    age::Union{Nothing, Real} = nothing,
)::NeuroAnalyzer.NEURO
    # validate
    _check_datatype(obj, "nirs")
    length(get_channel(obj; type = "nirs_od")) > 0 || throw(
        ArgumentError("OBJ does not contain NIRS OD channels, use intensity2od() first."),
    )
    rec = obj.header.recording
    wl = rec[:wavelengths]
    !isnothing(ppf) && !isnothing(age) &&
        throw(ArgumentError("ppf and age cannot be specified together."))
    if !isnothing(ppf)
        length(ppf) == length(wl) || throw(
            ArgumentError("ppf length does not correspond to the number of wavelengths."),
        )
        all(>(0), ppf) || throw(ArgumentError("ppf values must be > 0."))
    end
    _check_channels(get_channel(obj; type = "nirs_od"), ch)

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # path length factors
    pf = if !isnothing(age)
        nirs_dpf(; age = age, wl = wl)
    elseif !isnothing(ppf)
        Float64.(ppf)
    else
        ones(length(wl))
    end
    # divide by the source-detector distance unless all factors are 1 (μM·mm convention)
    use_rho = !all(≈(1), pf)

    ep_len = epoch_len(obj)
    ep_n = nepochs(obj)
    wl_idx = rec[:wavelength_index][ch]
    chp = rec[:optode_pairs][ch, :]

    # extinction coefficients (1/(M·cm) → 1/(M·mm)) and their least-squares inverse
    e = zeros(length(wl), 2)
    for w_idx in eachindex(wl)
        e[w_idx, :] = _wl2ext(wl[w_idx])
    end
    e = e ./ 10
    einv = inv(e' * e) * e'

    # one concentration triplet per optode pair, from its first-wavelength channel;
    # OD channels of the same pair at the other wavelengths
    lst = findall(wl_idx .== 1)
    isempty(lst) && throw(ArgumentError("No OD channels of the first wavelength selected."))
    other_wl(idx1) =
        findall(wl_idx .> 1 .&& chp[:, 1] .== chp[idx1, 1] .&& chp[:, 2] .== chp[idx1, 2])

    # source-detector distances in mm (3-D optode positions)
    rho = ones(length(lst))
    if use_rho
        scale = _length_to_mm(string(get(rec, :length_unit, "")))
        if isnothing(scale)
            _warn(
                "Length unit of the optode positions is unknown; positions are assumed to be in mm.",
            )
            scale = 1.0
        end
        for (idx, idx1) in enumerate(lst)
            sl = rec[:src_labels][chp[idx1, 1]]
            dl = rec[:det_labels][chp[idx1, 2]]
            si = findfirst(==(sl), obj.locs[!, :label])
            di = findfirst(==(dl), obj.locs[!, :label])
            (isnothing(si) || isnothing(di)) &&
                throw(ArgumentError("Optode locations of $sl / $dl not found in OBJ locs."))
            d = sqrt(
                (obj.locs[si, :loc_x] - obj.locs[di, :loc_x])^2 +
                (obj.locs[si, :loc_y] - obj.locs[di, :loc_y])^2 +
                (obj.locs[si, :loc_z] - obj.locs[di, :loc_z])^2,
            )
            d > 0 || throw(ArgumentError("Source-detector distance of $sl-$dl is 0."))
            rho[idx] = d * scale
        end
    end

    # concentration changes: (HbO, HbR, HbT) × samples × pairs × epochs
    dc = zeros(3, ep_len, length(lst), ep_n)
    @inbounds for ep_idx = 1:ep_n
        dod = @view(obj.data[ch, :, ep_idx])
        for (idx, idx1) in enumerate(lst)
            idx_all = vcat(idx1, other_wl(idx1))
            # path length factor of each OD channel (by its wavelength)
            l = rho[idx] .* pf[wl_idx[idx_all]]
            dc[1:2, :, idx, ep_idx] = einv * (dod[idx_all, :] ./ l)
        end
        dc[3, :, :, ep_idx] = dc[1, :, :, ep_idx] + dc[2, :, :, ep_idx]
    end
    # M → μM
    dc .*= 10^6

    # new channels: HbO, HbR, HbT of each optode pair
    hb_data = reduce(vcat, [dc[:, :, idx, :] for idx in axes(dc, 3)])
    n_hb = size(hb_data, 1)
    hb_type = repeat(["nirs_hbo", "nirs_hbr", "nirs_hbt"], length(lst))
    hb_unit = fill(use_rho ? "μM" : "μM·mm", n_hb)
    hb_label = String[]
    hb_bad = Bool[]
    hb_pairs = zeros(Int64, n_hb, 2)
    for (idx, idx1) in enumerate(lst)
        l = split(rec[:label][ch[idx1]], ' ')[1]
        append!(hb_label, ["$l HbO", "$l HbR", "$l HbT"])
        # a concentration channel is bad if any of its OD channels is bad
        append!(hb_bad, fill(any(rec[:bad_channel][ch[vcat(idx1, other_wl(idx1))]]), 3))
        hb_pairs[(3 * idx - 2):(3 * idx), :] .= chp[idx1, :]'
    end

    # insert the new channels after the optical channels (described by :wavelength_index
    # and :optode_pairs), before the remaining (e.g. auxiliary) channels
    n_opt = length(rec[:wavelength_index])
    opt = 1:n_opt
    rest = (n_opt + 1):size(obj.data, 1)
    obj_tmp = deepcopy(obj)
    obj_tmp.data = vcat(obj.data[opt, :, :], hb_data, obj.data[rest, :, :])
    r = obj_tmp.header.recording
    r[:channel_type] = vcat(rec[:channel_type][opt], hb_type, rec[:channel_type][rest])
    r[:unit] = vcat(rec[:unit][opt], hb_unit, rec[:unit][rest])
    r[:label] = vcat(rec[:label][opt], hb_label, rec[:label][rest])
    r[:bad_channel] = vcat(rec[:bad_channel][opt], hb_bad, rec[:bad_channel][rest])
    # concentration channels belong to their optode pair, but not to a wavelength (index 0)
    r[:wavelength_index] = vcat(rec[:wavelength_index], zeros(Int64, n_hb))
    r[:optode_pairs] = vcat(rec[:optode_pairs], hb_pairs)
    r[:channel_order] = _sort_channels(r[:channel_type])

    push!(
        obj_tmp.history,
        "od2conc(obj; ch=$ch, ppf=$(isnothing(ppf) ? "nothing" : ppf), age=$(isnothing(age) ? "nothing" : age))",
    )

    return obj_tmp
end

"""
    od2conc!(obj; <keyword arguments>)

Convert NIRS optical density (OD) to concentration (HbO, HbR, HbT) using the modified Beer-Lambert law. See `od2conc()` for details.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj, type="nirs_od"))`: list of channels, default is all NIRS OD channels
- `ppf::Union{Nothing, Vector{<:Real}}=nothing`: partial path length factors for each wavelength
- `age::Union{Nothing, Real}=nothing`: age of the subject in years, used to calculate the DPF; cannot be used together with `ppf`

# Returns

- `Nothing`
"""
function od2conc!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = get_channel(obj; type = "nirs_od"),
    ppf::Union{Nothing, Vector{<:Real}} = nothing,
    age::Union{Nothing, Real} = nothing,
)::Nothing
    obj_tmp = od2conc(obj; ch = ch, ppf = ppf, age = age)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
