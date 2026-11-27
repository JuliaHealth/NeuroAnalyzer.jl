export nirs_distance
export nirs_sci
export nirs_cv
export nirs_quality

# NIRS: optical channels (described by :wavelength_index / :optode_pairs) come first; for
# each source-detector pair: its label ("S1_D1"), its row in :optode_pairs and the channels
# of the pair of a given type
function _nirs_pairs(obj::NeuroAnalyzer.NEURO)::Vector{NamedTuple}
    rec = obj.header.recording
    chp = rec[:optode_pairs]
    n_opt = size(chp, 1)
    out = NamedTuple[]
    seen = Tuple{Int64, Int64}[]
    for j = 1:n_opt
        p = (chp[j, 1], chp[j, 2])
        (p[1] == 0 || p[2] == 0 || p in seen) && continue
        push!(seen, p)
        ch = [k for k = 1:n_opt if (chp[k, 1], chp[k, 2]) == p]
        label = rec[:src_labels][p[1]] * "_" * rec[:det_labels][p[2]]
        push!(out, (; label, src = p[1], det = p[2], ch))
    end

    return out
end

# channels of `pair` of type `type`, ordered by wavelength index
function _nirs_pair_ch(obj::NeuroAnalyzer.NEURO, pair::NamedTuple, type::String)::Vector{Int64}
    rec = obj.header.recording
    ch = [k for k in pair.ch if rec[:channel_type][k] == type]

    return ch[sortperm(rec[:wavelength_index][ch])]
end

"""
    nirs_distance(obj)

Calculate the source-detector distance of each NIRS optical channel (the channels described by the optode pairs, i.e. intensity, OD and concentration channels).

The distances are calculated from the 3-D optode positions in `obj.locs` and converted to mm using the unit of the locations (`locs_coord()`, or guessed from the head size if it is not known); if it cannot be determined, mm are assumed (with a warning).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

Named tuple:

- `ch::Vector{String}`: channel labels
- `d::Vector{Float64}`: source-detector distances in mm
"""
function nirs_distance(obj::NeuroAnalyzer.NEURO)::@NamedTuple{
    ch::Vector{String},
    d::Vector{Float64},
}
    # validate
    _check_chtype(obj, "nirs")
    rec = obj.header.recording
    chp = rec[:optode_pairs]
    n_opt = size(chp, 1)
    n_opt > 0 || throw(ArgumentError("OBJ does not contain NIRS optical channels."))

    lc = locs_coord(obj)
    # (a guessed unit is not used for distances)
    scale = lc.guessed ? nothing : _length_to_mm(lc.unit)
    if isnothing(scale)
        _warn("Length unit of the optode positions is unknown; positions are assumed to be in mm.")
        scale = 1.0
    end

    d = zeros(n_opt)
    for j = 1:n_opt
        (chp[j, 1] == 0 || chp[j, 2] == 0) && (d[j] = NaN; continue)
        sl = rec[:src_labels][chp[j, 1]]
        dl = rec[:det_labels][chp[j, 2]]
        si = findfirst(==(sl), obj.locs[!, :label])
        di = findfirst(==(dl), obj.locs[!, :label])
        (isnothing(si) || isnothing(di)) &&
            throw(ArgumentError("Optode locations of $sl / $dl not found in OBJ locs."))
        d[j] = scale * sqrt(
            (obj.locs[si, :loc_x] - obj.locs[di, :loc_x])^2 +
            (obj.locs[si, :loc_y] - obj.locs[di, :loc_y])^2 +
            (obj.locs[si, :loc_z] - obj.locs[di, :loc_z])^2,
        )
    end

    return (; ch = rec[:label][1:n_opt], d)
end

"""
    nirs_sci(obj; <keyword arguments>)

Calculate the scalp coupling index (SCI) of each NIRS source-detector pair.

The SCI is the correlation between the signals of the two wavelengths of a pair, band-pass filtered in the range of the cardiac pulsation (`flim`). With good optode-scalp coupling, the heartbeat is visible in both wavelengths and the SCI is close to 1; values below about 0.75 indicate poor coupling. The intensity channels are used (or the OD channels, if there are no intensity channels). For pairs with more than two wavelengths, the first two are used.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `flim::Tuple{Real, Real}=(0.5, 2.5)`: frequency range of the cardiac pulsation in Hz

# Returns

Named tuple:

- `pair::Vector{String}`: source-detector pairs (e.g. `"S1_D1"`)
- `sci::Matrix{Float64}`: SCI, shape (pairs, epochs)

# References

1. Pollonini L, Olds C, Abaya H, Bortfeld H, Beauchamp MS, Oghalai JS. Auditory cortex activation to natural speech and simulated cochlear implant speech measured with functional near-infrared spectroscopy. Hearing Research. 2014;309:84–93.
"""
function nirs_sci(
    obj::NeuroAnalyzer.NEURO;
    flim::Tuple{Real, Real} = (0.5, 2.5),
)::@NamedTuple{
    pair::Vector{String},
    sci::Matrix{Float64},
}
    # validate
    _check_chtype(obj, "nirs")
    fs = sr(obj)
    _check_tuple(flim, (0, fs / 2), "flim")
    type = !isempty(get_channel(obj; type = "nirs_int")) ? "nirs_int" : "nirs_od"
    isempty(get_channel(obj; type = type)) &&
        throw(ArgumentError("OBJ does not contain NIRS intensity or OD channels."))

    flt = digitalfilter(Bandpass(flim[1], flim[2]), Butterworth(4); fs = fs)
    pairs = [p for p in _nirs_pairs(obj) if length(_nirs_pair_ch(obj, p, type)) >= 2]
    isempty(pairs) && throw(ArgumentError("No source-detector pairs with two wavelengths."))
    ep_n = nepochs(obj)
    sci = zeros(length(pairs), ep_n)
    for (idx, p) in enumerate(pairs), ep_idx = 1:ep_n
        c = _nirs_pair_ch(obj, p, type)
        s1 = filtfilt(flt, obj.data[c[1], :, ep_idx] .- mean(obj.data[c[1], :, ep_idx]))
        s2 = filtfilt(flt, obj.data[c[2], :, ep_idx] .- mean(obj.data[c[2], :, ep_idx]))
        sci[idx, ep_idx] = (std(s1) > 0 && std(s2) > 0) ? cor(s1, s2) : 0.0
    end

    return (; pair = [p.label for p in pairs], sci)
end

"""
    nirs_cv(obj)

Calculate the coefficient of variation (CV) of the NIRS intensity channels: `100 × SD / mean` (in %). High CV (e.g. > 15%) indicates a noisy channel (poor coupling, motion, low light level).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

Named tuple:

- `ch::Vector{String}`: channel labels
- `cv::Matrix{Float64}`: CV in %, shape (channels, epochs)
"""
function nirs_cv(obj::NeuroAnalyzer.NEURO)::@NamedTuple{
    ch::Vector{String},
    cv::Matrix{Float64},
}
    # validate
    _check_chtype(obj, "nirs")
    ch = get_channel(obj; ch = get_channel(obj; type = "nirs_int"))
    isempty(ch) && throw(ArgumentError("OBJ does not contain NIRS intensity channels."))

    cv = zeros(length(ch), nepochs(obj))
    for (idx, c) in enumerate(ch), ep_idx = 1:nepochs(obj)
        s = @view(obj.data[c, :, ep_idx])
        m = mean(abs, s)
        cv[idx, ep_idx] = m > 0 ? 100 * std(s) / m : Inf
    end

    return (; ch = labels(obj)[ch], cv)
end

"""
    nirs_quality(obj; <keyword arguments>)

Assess the signal quality of each NIRS source-detector pair: scalp coupling index (SCI, see `nirs_sci()`), coefficient of variation of the intensity channels (CV, see `nirs_cv()`; the largest CV of the pair's wavelengths) and the source-detector distance (see `nirs_distance()`). A pair is marked as bad if its SCI is below `sci_t` or its CV is above `cv_t`. For epoched objects, the SCI and CV are averaged over epochs.

Bad pairs can be removed using `delete_pair()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `sci_t::Real=0.75`: SCI threshold
- `cv_t::Real=15`: CV threshold in % (ignored if there are no intensity channels)
- `flim::Tuple{Real, Real}=(0.5, 2.5)`: frequency range of the cardiac pulsation in Hz (for the SCI)

# Returns

- `DataFrame`: columns `pair`, `distance` (mm), `sci`, `cv` (`NaN` without intensity channels), `bad`
"""
function nirs_quality(
    obj::NeuroAnalyzer.NEURO;
    sci_t::Real = 0.75,
    cv_t::Real = 15,
    flim::Tuple{Real, Real} = (0.5, 2.5),
)::DataFrame
    # validate
    _check_chtype(obj, "nirs")
    -1 <= sci_t <= 1 || throw(ArgumentError("sci_t must be in [-1, 1]."))
    cv_t > 0 || throw(ArgumentError("cv_t must be > 0."))

    s = nirs_sci(obj; flim = flim)
    sci = vec(mean(s.sci; dims = 2))
    pairs = Dict(p.label => p for p in _nirs_pairs(obj))
    dist = nirs_distance(obj).d

    cv = fill(NaN, length(s.pair))
    has_int = !isempty(get_channel(obj; type = "nirs_int"))
    if has_int
        c = nirs_cv(obj)
        cv_ch = Dict(zip(c.ch, vec(mean(c.cv; dims = 2))))
        for (idx, p) in enumerate(s.pair)
            ch = _nirs_pair_ch(obj, pairs[p], "nirs_int")
            cv[idx] = maximum(cv_ch[labels(obj)[k]] for k in ch)
        end
    end
    d = [dist[first(pairs[p].ch)] for p in s.pair]
    bad = (sci .< sci_t) .|| (has_int .&& (cv .> cv_t))

    return DataFrame(; pair = s.pair, distance = d, sci = sci, cv = cv, bad = bad)
end

# NIRS channel locations: midpoints of the source and detector of each optical channel
# (labels = channel labels), e.g. for topographical maps of channel values
function _nirs_chlocs(obj::NeuroAnalyzer.NEURO)::DataFrame
    rec = obj.header.recording
    chp = rec[:optode_pairs]
    lbl = String[]
    xyz = Vector{Float64}[]
    for j in axes(chp, 1)
        (chp[j, 1] == 0 || chp[j, 2] == 0) && continue
        si = findfirst(==(rec[:src_labels][chp[j, 1]]), obj.locs[!, :label])
        di = findfirst(==(rec[:det_labels][chp[j, 2]]), obj.locs[!, :label])
        (isnothing(si) || isnothing(di)) && continue
        push!(lbl, rec[:label][j])
        push!(
            xyz,
            [(obj.locs[si, c] + obj.locs[di, c]) / 2 for c in (:loc_x, :loc_y, :loc_z)],
        )
    end
    n = length(lbl)
    locs = DataFrame(
        :label => lbl,
        :loc_radius => zeros(n),
        :loc_theta => zeros(n),
        :loc_x => [p[1] for p in xyz],
        :loc_y => [p[2] for p in xyz],
        :loc_z => [p[3] for p in xyz],
        :loc_radius_sph => zeros(n),
        :loc_theta_sph => zeros(n),
        :loc_phi_sph => zeros(n),
    )
    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    return locs
end
