export nirs_ssr
export nirs_ssr!

"""
    nirs_ssr(obj; <keyword arguments>)

Short-separation regression of NIRS signals.

Short-separation channels (source-detector distance below `short_t`, typically ~8 mm) measure mostly the extracerebral (scalp) hemodynamics: systemic physiology such as the heartbeat, respiration and Mayer waves, and blood flow changes in the skin. For each long channel, the nearest short channel (distance between the midpoints of the two pairs) of the same type (and, for OD channels, the same wavelength) is regressed out, separately in each epoch:

`y_corrected = y − β · (x − mean(x))`, where `β = cov(x, y) / var(x)`

The short channels are not changed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj; type = ["nirs_hbo", "nirs_hbr"])`: long channels to correct, default is all HbO and HbR channels
- `short_t::Real=15`: channels with source-detector distance < `short_t` mm are short-separation channels

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object

# References

1. Saager RB, Berger AJ. Direct characterization and removal of interfering absorption trends in two-layer turbid media. Journal of the Optical Society of America A. 2005;22(9):1874–1882.
2. Scholkmann F, Kleiser S, Metz AJ, et al. A review on continuous wave functional near-infrared spectroscopy and imaging instrumentation and methodology. NeuroImage. 2014;85:6–27.
"""
function nirs_ssr(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = get_channel(obj; type = ["nirs_hbo", "nirs_hbr"]),
    short_t::Real = 15,
)::NeuroAnalyzer.NEURO
    # validate
    _check_datatype(obj, "nirs")
    short_t > 0 || throw(ArgumentError("short_t must be > 0."))
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    rec = obj.header.recording
    n_opt = size(rec[:optode_pairs], 1)
    all(<=(n_opt), ch) ||
        throw(ArgumentError("Only NIRS optical channels (with optode pairs) can be corrected."))

    d = nirs_distance(obj).d
    ct = rec[:channel_type]
    wl = rec[:wavelength_index]
    short = [j for j = 1:n_opt if d[j] < short_t]
    isempty(short) && throw(ArgumentError("No short-separation channels (distance < $short_t mm)."))

    # channel positions: midpoints of the source and detector (in the units of OBJ locs)
    function _mid(j)
        sl = rec[:src_labels][rec[:optode_pairs][j, 1]]
        dl = rec[:det_labels][rec[:optode_pairs][j, 2]]
        si = findfirst(==(sl), obj.locs[!, :label])
        di = findfirst(==(dl), obj.locs[!, :label])
        return [
            (obj.locs[si, c] + obj.locs[di, c]) / 2 for c in (:loc_x, :loc_y, :loc_z)
        ]
    end

    obj_tmp = deepcopy(obj)
    pairs_used = String[]
    for c in ch
        d[c] < short_t && continue
        # short channels of the same type (and wavelength for intensity / OD channels)
        cand = [
            s for s in short if ct[s] == ct[c] && (wl[c] == 0 || wl[s] == wl[c])
        ]
        if isempty(cand)
            _warn("No short-separation channel of type $(ct[c]) for $(labels(obj)[c]).")
            continue
        end
        mc = _mid(c)
        s = cand[argmin([sum(abs2, _mid(k) .- mc) for k in cand])]
        push!(pairs_used, "$(labels(obj)[c]) ← $(labels(obj)[s])")
        for ep = 1:nepochs(obj)
            y = @view(obj.data[c, :, ep])
            x = obj.data[s, :, ep] .- mean(obj.data[s, :, ep])
            vx = sum(abs2, x)
            vx > 0 || continue
            b = sum(x .* (y .- mean(y))) / vx
            obj_tmp.data[c, :, ep] = y .- b .* x
        end
    end
    isempty(pairs_used) || _info("Short-separation regression: $(join(pairs_used, ", "))")
    push!(obj_tmp.history, "nirs_ssr(obj; ch=$(labels(obj)[ch]), short_t=$short_t)")

    return obj_tmp
end

"""
    nirs_ssr!(obj; <keyword arguments>)

Short-separation regression of NIRS signals. See `nirs_ssr()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj; type = ["nirs_hbo", "nirs_hbr"])`: long channels to correct, default is all HbO and HbR channels
- `short_t::Real=15`: channels with source-detector distance < `short_t` mm are short-separation channels

# Returns

- `Nothing`
"""
function nirs_ssr!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = get_channel(obj; type = ["nirs_hbo", "nirs_hbr"]),
    short_t::Real = 15,
)::Nothing
    obj_tmp = nirs_ssr(obj; ch = ch, short_t = short_t)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history

    return nothing
end
