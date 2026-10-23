export reference_ce
export reference_ce!
export reference_avg
export reference_avg!
export reference_a
export reference_a!
export reference_m
export reference_m!
export reference_plap
export reference_plap!
export reference_slap
export reference_slap!
export reference_custom
export reference_custom!

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

# indices (into obj.data) of EEG channels, in data order
function _eeg_ch_idx(obj::NeuroAnalyzer.NEURO)::Vector{Int64}
    ch = get_channel(obj; ch = get_channel(obj; type = "eeg"))
    isempty(ch) && throw(ArgumentError("OBJ does not contain EEG channels."))
    return ch
end

# channel locations for channels `ch`, row i corresponding to ch[i]
function _locs_ch_order(obj::NeuroAnalyzer.NEURO, ch::Vector{Int64})::DataFrame
    _has_locs(obj)
    _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
    rows = _find_bylabel(obj.locs, labels(obj)[ch])
    length(rows) == length(ch) ||
        throw(ArgumentError("Could not match all channels to channel locations."))
    return obj.locs[rows, :]
end

# 2-D (topographic) electrode positions from polar coordinates, as used for plotting
function _locs_polar_xy(locs::DataFrame)::Tuple{Vector{Float64}, Vector{Float64}}
    r = Float64.(locs[!, :loc_radius])
    θ = Float64.(locs[!, :loc_theta])
    all(iszero, r) &&
        throw(ArgumentError("Polar coordinates are not available, use locs_cart2pol() first."))
    return r .* cosd.(θ), r .* sind.(θ)
end

# great-circle distances (radians) between electrode positions projected onto the unit sphere
function _locs_gc_dist(locs::DataFrame)::Matrix{Float64}
    xyz = hcat(Float64.(locs[!, :loc_x]), Float64.(locs[!, :loc_y]), Float64.(locs[!, :loc_z]))
    r = sqrt.(sum(abs2, xyz; dims = 2))
    any(iszero, r) && throw(
        ArgumentError("Some channels have no Cartesian coordinates, use locs_sph2cart() first."),
    )
    xyz ./= r
    n = size(xyz, 1)
    # atan(|a×b|, a⋅b) is exact for identical and nearby positions (acos(a⋅b) is not)
    return [atan(norm(cross(xyz[i, :], xyz[j, :])), dot(xyz[i, :], xyz[j, :])) for i = 1:n, j = 1:n]
end

# inverse-distance weights; zero/non-finite distances get weight 0 (excluded)
function _inv_dist_weights(d::AbstractMatrix{Float64})::Matrix{Float64}
    return [(isfinite(x) && x > 0) ? 1 / x : 0.0 for x in d]
end

# append suffix (one String for all, or one per channel) to labels in header and locs
function _append_label_suffix!(
    obj_tmp::NeuroAnalyzer.NEURO,
    obj::NeuroAnalyzer.NEURO,
    ch::Vector{Int64},
    suffix::Union{String, Vector{String}},
)::Nothing
    lbl = labels(obj)[ch]
    sfx = suffix isa String ? fill(suffix, length(ch)) : suffix
    length(sfx) == length(ch) || throw(ArgumentError("Suffix length mismatch."))
    obj_tmp.header.recording[:label][ch] .*= sfx
    for (i, l) in enumerate(lbl)
        r = _find_bylabel(obj_tmp.locs, l)
        r isa Int64 || continue  # channel has no location
        obj_tmp.locs[r, :label] = obj_tmp.locs[r, :label] * sfx[i]
    end
    return nothing
end

# column-wise weighted median of X (rows = channels, columns = samples), weights w > 0;
# when the cumulative weight hits exactly half, the two middle values are averaged
# (so equal weights reproduce the ordinary median)
function _wmedian_cols(
    X::AbstractMatrix{Float64},
    w::AbstractVector{Float64},
)::Vector{Float64}
    k, n = size(X)
    length(w) == k || throw(ArgumentError("Length of weights must equal number of rows."))
    out = Vector{Float64}(undef, n)
    v = Vector{Float64}(undef, k)
    p = Vector{Int64}(undef, k)
    half = sum(w) / 2
    tol = k * eps(half)
    @inbounds for j = 1:n
        for i = 1:k
            v[i] = X[i, j]
        end
        sortperm!(p, v)
        c = 0.0
        for m = 1:k
            c += w[p[m]]
            if c >= half - tol
                out[j] =
                    (abs(c - half) <= tol && m < k) ? (v[p[m]] + v[p[m + 1]]) / 2 :
                    v[p[m]]
                break
            end
        end
    end
    return out
end

function _copy_into!(obj::NeuroAnalyzer.NEURO, obj_tmp::NeuroAnalyzer.NEURO)::Nothing
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj.locs = obj_tmp.locs
    return nothing
end

"""
    reference_ce(obj; <keyword arguments>)

Re-reference EEG channels to a common electrode or the average of multiple electrodes. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`
- `ch::Union{String, Vector{String}, Regex}`: reference channel name(s); if multiple channels are given their mean (or median) is used as the reference
- `med::Bool=false`: if `true`, use median instead of mean for multi-channel reference

# Returns

- `NeuroAnalyzer.NEURO`: new object with re-referenced EEG channels

# Notes

- Reference channels that are EEG channels are re-referenced too; with a single reference electrode that channel becomes flat (zero).
"""
function reference_ce(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    med::Bool = false,
)::NeuroAnalyzer.NEURO
    # validate
    _check_datatype(obj, "eeg")

    # reference channel indices
    ch == "" && throw(ArgumentError("No reference channels selected."))
    ref_ch_idx = get_channel(obj; ch = ch)
    isempty(ref_ch_idx) && throw(ArgumentError("No reference channels selected."))
    single = length(ref_ch_idx) == 1

    # channels to re-reference
    sig_ch_idx = _eeg_ch_idx(obj)
    ep_n = nepochs(obj)

    # create new dataset
    obj_tmp = deepcopy(obj)

    # the reference is the same for all channels within an epoch
    Threads.@threads :static for ep_idx = 1:ep_n
        ref = if single
            @view obj.data[ref_ch_idx[1], :, ep_idx]
        else
            src = @view obj.data[ref_ch_idx, :, ep_idx]
            med ? vec(median(src; dims = 1)) : vec(mean(src; dims = 1))
        end
        for s_ch in sig_ch_idx
            @views obj_tmp.data[s_ch, :, ep_idx] .= obj.data[s_ch, :, ep_idx] .- ref
        end
    end

    suffix = single ? "-$(labels(obj)[ref_ch_idx[1]])" : (med ? "-cmed" : "-cavg")
    _append_label_suffix!(obj_tmp, obj, sig_ch_idx, suffix)

    obj_tmp.header.recording[:reference] = if single
        "common ($(labels(obj)[ref_ch_idx[1]]))"
    else
        "common ($(Base.join(labels(obj)[ref_ch_idx], ", ")) $(med ? "median" : "averaged"))"
    end
    push!(obj_tmp.history, "reference_ce(obj; ch=$ch, med=$med)")

    return obj_tmp
end

"""
    reference_ce!(obj; <keyword arguments>)

Re-reference EEG channels to a common electrode in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: reference channel name(s); if multiple channels are given their mean (or median) is used as the reference
- `med::Bool=false`: if `true`, use median instead of mean for multi-channel reference

# Returns

- `Nothing`

# Notes

- Reference channels that are EEG channels are re-referenced too; with a single reference electrode that channel becomes flat (zero).
"""
function reference_ce!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    med::Bool = false,
)::Nothing
    _copy_into!(obj, reference_ce(obj; ch = ch, med = med))
    return nothing
end

"""
    reference_avg(obj; <keyword arguments>)

Re-reference EEG channels to the common average reference (CAR). Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`
- `exclude_fpo::Bool=false`: exclude Fp1, Fp2 (due to eye blinks), O1, O2 (due to head movements) from the CAR calculation
- `exclude_current::Bool=false`: exclude the current channel from its own CAR calculation
- `average::Bool=true`: subtract the mean (or median) reference; if `false`, subtract the sum
- `med::Bool=false`: use median instead of mean (weighted median when `weighted=true`); ignored when `average=false`
- `weighted::Bool=false`: weight reference channels by inverse planar distance in the 2-D topographic (polar) projection of electrode positions (requires polar channel locations)

# Returns

- `NeuroAnalyzer.NEURO`: new object with re-referenced EEG channels

# Notes

- Weighted reference: `Σ wᵢxᵢ / Σ wᵢ` with `wᵢ = 1/dᵢ`; with `average=false` this is multiplied by the number of reference channels (weighted analogue of the sum).
- With `weighted=true` the current channel (distance 0) and co-located channels always get weight 0.
"""
function reference_avg(
    obj::NeuroAnalyzer.NEURO;
    exclude_fpo::Bool = false,
    exclude_current::Bool = false,
    average::Bool = true,
    med::Bool = false,
    weighted::Bool = false,
)::NeuroAnalyzer.NEURO
    # validate
    _check_datatype(obj, "eeg")

    # channels that will be referenced
    sig_ch_idx = _eeg_ch_idx(obj)
    sig_lbl = labels(obj)[sig_ch_idx]
    ch_n = length(sig_ch_idx)
    ep_n = nepochs(obj)

    # rows (positions within sig_ch_idx) of Fp/O channels to optionally exclude
    fpo_rows = Int64[]
    if exclude_fpo
        l = lowercase.(sig_lbl)
        for lbl in ("fp1", "fp2", "o1", "o2")
            i = findfirst(isequal(lbl), l)
            !isnothing(i) && push!(fpo_rows, i)
        end
    end

    # inverse-distance weight matrix (row i: weights for channel i)
    W = Matrix{Float64}(undef, 0, 0)
    if weighted
        x, y = _locs_polar_xy(_locs_ch_order(obj, sig_ch_idx))
        d = [hypot(x[i] - x[j], y[i] - y[j]) for i = 1:ch_n, j = 1:ch_n]
        W = _inv_dist_weights(d)
    end

    # reference rows per channel (do not depend on epoch)
    ref_rows = Vector{Vector{Int64}}(undef, ch_n)
    for i = 1:ch_n
        rr = setdiff(1:ch_n, exclude_current ? vcat(fpo_rows, i) : fpo_rows)
        weighted && (rr = Base.filter(j -> W[i, j] > 0, rr))
        isempty(rr) && throw(
            ArgumentError("No reference channels left for $(sig_lbl[i]) after exclusions."),
        )
        ref_rows[i] = rr
    end
    # same reference for all channels → compute once per epoch
    shared = !weighted && !exclude_current

    # read-only source / output buffer
    src = obj.data[sig_ch_idx, :, :]
    dst = similar(src)

    _ref(X) =
        average ? (med ? vec(median(X; dims = 1)) : vec(mean(X; dims = 1))) :
        vec(sum(X; dims = 1))

    Threads.@threads :static for ep_idx = 1:ep_n
        ref_shared = shared ? _ref(@view src[ref_rows[1], :, ep_idx]) : nothing
        for i = 1:ch_n
            ref = if shared
                ref_shared
            elseif weighted
                rr = ref_rows[i]
                w = W[i, rr]
                X = @view src[rr, :, ep_idx]
                if average && med
                    _wmedian_cols(X, w)
                else
                    r = vec(w' * X) ./ sum(w)
                    average ? r : r .* length(rr)
                end
            else
                _ref(@view src[ref_rows[i], :, ep_idx])
            end
            @views dst[i, :, ep_idx] .= src[i, :, ep_idx] .- ref
        end
    end

    # create new dataset
    obj_tmp = deepcopy(obj)
    obj_tmp.data[sig_ch_idx, :, :] = dst

    # med is ignored when average=false, so it only affects the averaged suffixes
    suffix = (weighted ? "-w" : "-") * (average ? (med ? "med" : "avg") : "sum")
    _append_label_suffix!(obj_tmp, obj, sig_ch_idx, suffix)

    obj_tmp.header.recording[:reference] = if average
        string(med ? "median" : "average", weighted ? " (weighted)" : "")
    else
        weighted ? "sum (weighted)" : "sum"
    end
    push!(
        obj_tmp.history,
        "reference_avg(obj, exclude_fpo=$exclude_fpo, exclude_current=$exclude_current, average=$average, med=$med, weighted=$weighted)",
    )

    return obj_tmp
end

"""
    reference_avg!(obj; <keyword arguments>)

Re-reference EEG channels to the common average reference in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`; modified in-place
- `exclude_fpo::Bool=false`: exclude Fp1, Fp2 (due to eye blinks), O1, O2 (due to head movements) from the CAR calculation
- `exclude_current::Bool=false`: exclude the current channel from its own CAR calculation
- `average::Bool=true`: subtract the mean (or median) reference; if `false`, subtract the sum
- `med::Bool=false`: use median instead of mean (weighted median when `weighted=true`); ignored when `average=false`
- `weighted::Bool=false`: weight reference channels by inverse planar distance in the 2-D topographic (polar) projection of electrode positions (requires polar channel locations)

# Returns

- `Nothing`

# Notes

- Weighted reference: `Σ wᵢxᵢ / Σ wᵢ` with `wᵢ = 1/dᵢ`; with `average=false` this is multiplied by the number of reference channels (weighted analogue of the sum).
- With `weighted=true` the current channel (distance 0) and co-located channels always get weight 0.
"""
function reference_avg!(
    obj::NeuroAnalyzer.NEURO;
    exclude_fpo::Bool = false,
    exclude_current::Bool = false,
    average::Bool = true,
    med::Bool = false,
    weighted::Bool = false,
)::Nothing
    _copy_into!(
        obj,
        reference_avg(
            obj;
            exclude_fpo = exclude_fpo,
            exclude_current = exclude_current,
            average = average,
            med = med,
            weighted = weighted,
        ),
    )
    return nothing
end

# ---------------------------------------------------------------------------
# internal helper: auricular/mastoid reference, shared by reference_a and reference_m
# ---------------------------------------------------------------------------
function _reference_paired(
    obj::NeuroAnalyzer.NEURO,
    r1::String,
    r2::String,
    type::Symbol,
    name::String,
)::NeuroAnalyzer.NEURO
    # validate
    _check_datatype(obj, "eeg")
    _check_var(type, [:l, :i, :c], "type")
    r1 in labels(obj) || throw(ArgumentError("OBJ does not contain $r1 channel."))
    r2 in labels(obj) || throw(ArgumentError("OBJ does not contain $r2 channel."))

    ch = _eeg_ch_idx(obj)
    lbl = labels(obj)[ch]
    ch_n = length(ch)
    ep_n = nepochs(obj)
    i1 = get_channel(obj; ch = r1)[1]
    i2 = get_channel(obj; ch = r2)[1]

    # reference per row of `ch`: 0 = linked, 1 = r1, 2 = r2
    assign = zeros(Int64, ch_n)
    if type !== :l
        rows_of = p -> findall(in(Set(channel_pick(obj; pick = p))), lbl)
        assign .= -1
        assign[rows_of(:central)] .= 0
        assign[rows_of(:left)] .= type === :i ? 1 : 2
        assign[rows_of(:right)] .= type === :i ? 2 : 1
        unc = findall(==(-1), assign)
        if !isempty(unc)
            _warn(
                "Channels not classified as left/central/right, using linked reference: $(Base.join(lbl[unc], ", ")).",
            )
            assign[unc] .= 0
        end
    end

    s_ref = similar(obj.data, ch_n, size(obj.data, 2), ep_n)
    Threads.@threads :static for ep_idx = 1:ep_n
        a = @view obj.data[i1, :, ep_idx]
        b = @view obj.data[i2, :, ep_idx]
        linked = (a .+ b) ./ 2
        for i = 1:ch_n
            ref = assign[i] == 0 ? linked : (assign[i] == 1 ? a : b)
            @views s_ref[i, :, ep_idx] .= obj.data[ch[i], :, ep_idx] .- ref
        end
    end

    # create new dataset
    obj_tmp = deepcopy(obj)
    obj_tmp.data[ch, :, :] = s_ref
    sfx = ["-$r1$r2", "-$r1", "-$r2"]
    _append_label_suffix!(obj_tmp, obj, ch, [sfx[k + 1] for k in assign])
    obj_tmp.header.recording[:reference] = Dict(
        :l => "$name (linked)",
        :i => "$name (ipsilateral)",
        :c => "$name (contralateral)",
    )[type]

    return obj_tmp
end

"""
    reference_a(obj; <keyword arguments>)

Re-reference EEG channels to auricular electrodes (A1, A2). Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must contain A1 and A2 channels; must be of type `"eeg"`
- `type::Symbol=:l`: reference type:
    - `:l`: linked - average of A1 and A2
    - `:i`: ipsilateral - A1 for left channels, A2 for right channels
    - `:c`: contralateral - A2 for left channels, A1 for right channels
- `med::Bool=false`: kept for API compatibility; has no effect (the median of two electrodes equals their mean)

# Returns

- `NeuroAnalyzer.NEURO`: new re-referenced object

# Notes

- For `:i`/`:c`, midline channels and channels not classified as left/right use the linked reference (a warning lists the unclassified ones).
"""
function reference_a(
    obj::NeuroAnalyzer.NEURO;
    type::Symbol = :l,
    med::Bool = false,
)::NeuroAnalyzer.NEURO
    obj_tmp = _reference_paired(obj, "A1", "A2", type, "auricular")
    push!(obj_tmp.history, "reference_a(obj, type=$type, med=$med)")
    return obj_tmp
end

"""
    reference_a!(obj; <keyword arguments>)

Re-reference EEG channels to auricular electrodes in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must contain A1 and A2 channels; must be of type `"eeg"`; modified in-place
- `type::Symbol=:l`: reference type:
    - `:l`: linked - average of A1 and A2
    - `:i`: ipsilateral - A1 for left channels, A2 for right channels
    - `:c`: contralateral - A2 for left channels, A1 for right channels
- `med::Bool=false`: kept for API compatibility; has no effect (the median of two electrodes equals their mean)

# Returns

- `Nothing`

# Notes

- For `:i`/`:c`, midline channels and channels not classified as left/right use the linked reference (a warning lists the unclassified ones).
"""
function reference_a!(
    obj::NeuroAnalyzer.NEURO;
    type::Symbol = :l,
    med::Bool = false,
)::Nothing
    _copy_into!(obj, reference_a(obj; type = type, med = med))
    return nothing
end

"""
    reference_m(obj; <keyword arguments>)

Re-reference EEG channels to mastoid electrodes (M1, M2). Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must contain M1 and M2 channels; must be of type `"eeg"`
- `type::Symbol=:l`: reference type:
    - `:l`: linked - average of M1 and M2
    - `:i`: ipsilateral - M1 for left channels, M2 for right channels
    - `:c`: contralateral - M2 for left channels, M1 for right channels
- `med::Bool=false`: kept for API compatibility; has no effect (the median of two electrodes equals their mean)

# Returns

- `NeuroAnalyzer.NEURO`: new re-referenced object

# Notes

- For `:i`/`:c`, midline channels and channels not classified as left/right use the linked reference (a warning lists the unclassified ones).
"""
function reference_m(
    obj::NeuroAnalyzer.NEURO;
    type::Symbol = :l,
    med::Bool = false,
)::NeuroAnalyzer.NEURO
    obj_tmp = _reference_paired(obj, "M1", "M2", type, "mastoid")
    push!(obj_tmp.history, "reference_m(obj, type=$type, med=$med)")
    return obj_tmp
end

"""
    reference_m!(obj; <keyword arguments>)

Re-reference EEG channels to mastoid electrodes in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must contain M1 and M2 channels; must be of type `"eeg"`; modified in-place
- `type::Symbol=:l`: reference type:
    - `:l`: linked - average of M1 and M2
    - `:i`: ipsilateral - M1 for left channels, M2 for right channels
    - `:c`: contralateral - M2 for left channels, M1 for right channels
- `med::Bool=false`: kept for API compatibility; has no effect (the median of two electrodes equals their mean)

# Returns

- `Nothing`

# Notes

- For `:i`/`:c`, midline channels and channels not classified as left/right use the linked reference (a warning lists the unclassified ones).
"""
function reference_m!(
    obj::NeuroAnalyzer.NEURO;
    type::Symbol = :l,
    med::Bool = false,
)::Nothing
    _copy_into!(obj, reference_m(obj; type = type, med = med))
    return nothing
end

# ---------------------------------------------------------------------------
# internal helper: shared nearest-neighbor Laplacian core
# used by both reference_plap and reference_slap
# ---------------------------------------------------------------------------
function _laplacian_reference(
    obj::NeuroAnalyzer.NEURO;
    nn::Int64,
    weighted::Bool,
    med::Bool,
    spherical::Bool,
)::Tuple{Vector{Int64}, Array{Float64, 3}}
    # validate
    _check_datatype(obj, "eeg")

    ch = _eeg_ch_idx(obj)
    ch_n = length(ch)
    ep_n = nepochs(obj)
    nn >= 1 || throw(ArgumentError("nn must be ≥ 1."))
    nn <= ch_n - 1 || throw(ArgumentError("nn must be ≤ $(ch_n - 1)."))

    # locations in channel order
    locs = _locs_ch_order(obj, ch)
    d = if spherical
        _locs_gc_dist(locs)
    else
        x, y = _locs_polar_xy(locs)
        [hypot(x[i] - x[j], y[i] - y[j]) for i = 1:ch_n, j = 1:ch_n]
    end
    # eliminate auto-referencing (and co-located electrodes)
    for i = 1:ch_n
        d[i, i] = Inf
    end
    d[d .== 0] .= Inf

    lbl = labels(obj)[ch]
    nn_idx = Vector{Vector{Int64}}(undef, ch_n)
    for i = 1:ch_n
        nb = sortperm(d[i, :])[1:nn]
        any(isinf, d[i, nb]) && throw(
            ArgumentError("Channel $(lbl[i]) has fewer than $nn neighbours at distinct locations."),
        )
        nn_idx[i] = nb
    end
    W = weighted ? [1 ./ d[i, nn_idx[i]] for i = 1:ch_n] : Vector{Vector{Float64}}()

    s_ref = similar(obj.data, ch_n, size(obj.data, 2), ep_n)
    Threads.@threads :static for ep_idx = 1:ep_n
        for i = 1:ch_n
            X = @view obj.data[ch[nn_idx[i]], :, ep_idx]
            ref = if weighted
                med ? _wmedian_cols(X, W[i]) : vec(W[i]' * X) ./ sum(W[i])
            else
                med ? vec(median(X; dims = 1)) : vec(mean(X; dims = 1))
            end
            @views s_ref[i, :, ep_idx] .= obj.data[ch[i], :, ep_idx] .- ref
        end
    end
    return ch, s_ref
end

# label suffix / header name for Laplacian references, e.g. "-wplapmed",
# "weighted planar Laplacian (median, 4)"
_lap_suffix(kind::String, weighted::Bool, med::Bool)::String =
    "-" * (weighted ? "w" : "") * kind * (med ? "med" : "")
_lap_name(kind::String, nn::Int64, weighted::Bool, med::Bool)::String =
    (weighted ? "weighted " : "") * "$kind Laplacian (" * (med ? "median, " : "") * "$nn)"

"""
    reference_plap(obj; <keyword arguments>)

Re-reference EEG channels using the planar Laplacian. Neighbours and distances are calculated in the 2-D topographic (polar) projection of electrode positions, as used for plotting. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must have channel locations (polar coordinates); must be of type `"eeg"`
- `nn::Int64=4`: number of nearest-neighbour electrodes; must satisfy `1 ≤ nn ≤ nch-1`
- `weighted::Bool=false`: weight neighbours by inverse distance (`Σ wᵢxᵢ / Σ wᵢ`, `wᵢ = 1/dᵢ`)
- `med::Bool=false`: use median instead of mean (weighted median when `weighted=true`)

# Returns

- `NeuroAnalyzer.NEURO`: new re-referenced object
"""
function reference_plap(
    obj::NeuroAnalyzer.NEURO;
    nn::Int64 = 4,
    weighted::Bool = false,
    med::Bool = false,
)::NeuroAnalyzer.NEURO
    ch, s_ref =
        _laplacian_reference(obj; nn = nn, weighted = weighted, med = med, spherical = false)

    # create new dataset
    obj_tmp = deepcopy(obj)
    obj_tmp.data[ch, :, :] = s_ref
    _append_label_suffix!(obj_tmp, obj, ch, _lap_suffix("plap", weighted, med))
    obj_tmp.header.recording[:reference] = _lap_name("planar", nn, weighted, med)
    push!(obj_tmp.history, "reference_plap(obj, nn=$nn, weighted=$weighted, med=$med)")

    return obj_tmp
end

"""
    reference_plap!(obj; <keyword arguments>)

Re-reference EEG channels using the planar Laplacian in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must have channel locations; must be of type `"eeg"`; modified in-place
- `nn::Int64=4`: number of nearest-neighbour electrodes; must satisfy `1 ≤ nn ≤ nch-1`
- `weighted::Bool=false`: weight neighbours by inverse distance
- `med::Bool=false`: use median instead of mean (weighted median when `weighted=true`)

# Returns

- `Nothing`
"""
function reference_plap!(
    obj::NeuroAnalyzer.NEURO;
    nn::Int64 = 4,
    weighted::Bool = false,
    med::Bool = false,
)::Nothing
    _copy_into!(obj, reference_plap(obj; nn = nn, weighted = weighted, med = med))
    return nothing
end

"""
    reference_slap(obj; <keyword arguments>)

Re-reference EEG channels using the spherical Laplacian. Neighbours and distances are calculated as great-circle distances between electrode positions projected onto the unit sphere. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must have channel locations (Cartesian coordinates); must be of type `"eeg"`
- `nn::Int64=4`: number of nearest-neighbour electrodes; must satisfy `1 ≤ nn ≤ nch-1`
- `weighted::Bool=false`: weight neighbours by inverse distance (`Σ wᵢxᵢ / Σ wᵢ`, `wᵢ = 1/dᵢ`)
- `med::Bool=false`: use median instead of mean (weighted median when `weighted=true`)

# Returns

- `NeuroAnalyzer.NEURO`: new re-referenced object
"""
function reference_slap(
    obj::NeuroAnalyzer.NEURO;
    nn::Int64 = 4,
    weighted::Bool = false,
    med::Bool = false,
)::NeuroAnalyzer.NEURO
    ch, s_ref =
        _laplacian_reference(obj; nn = nn, weighted = weighted, med = med, spherical = true)

    # create new dataset
    obj_tmp = deepcopy(obj)
    obj_tmp.data[ch, :, :] = s_ref
    _append_label_suffix!(obj_tmp, obj, ch, _lap_suffix("slap", weighted, med))
    obj_tmp.header.recording[:reference] = _lap_name("spherical", nn, weighted, med)
    push!(obj_tmp.history, "reference_slap(obj, nn=$nn, weighted=$weighted, med=$med)")

    return obj_tmp
end

"""
    reference_slap!(obj; <keyword arguments>)

Re-reference EEG channels using the spherical Laplacian in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must have channel locations; must be of type `"eeg"`; modified in-place
- `nn::Int64=4`: number of nearest-neighbour electrodes; must satisfy `1 ≤ nn ≤ nch-1`
- `weighted::Bool=false`: weight neighbours by inverse distance
- `med::Bool=false`: use median instead of mean (weighted median when `weighted=true`)

# Returns

- `Nothing`
"""
function reference_slap!(
    obj::NeuroAnalyzer.NEURO;
    nn::Int64 = 4,
    weighted::Bool = false,
    med::Bool = false,
)::Nothing
    _copy_into!(obj, reference_slap(obj; nn = nn, weighted = weighted, med = med))
    return nothing
end

# default montage: standard longitudinal bipolar ("double banana")
const _BIP_LONGITUDINAL = [
    "Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1",
    "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2",
    "F4-C4", "C4-P4", "P4-O2",
]

# resolve one montage entry into (channel, nothing) or (ch1, ch2);
# labels containing '-' are handled by testing every split point
function _parse_montage_entry(
    ref::String,
    lbls::Vector{String},
)::Tuple{String, Union{Nothing, String}}
    ref in lbls && return (ref, nothing)
    hits = Tuple{String, String}[]
    for i in findall(==('-'), ref)
        a = ref[1:prevind(ref, i)]
        b = ref[nextind(ref, i):end]
        a in lbls && b in lbls && push!(hits, (a, b))
    end
    isempty(hits) &&
        throw(ArgumentError("Montage entry \"$ref\": label(s) not found among EEG channels."))
    length(hits) > 1 && throw(ArgumentError("Montage entry \"$ref\" is ambiguous."))
    return hits[1]
end

# location of one montage channel: `:first` copies the location of the first electrode,
# `:mid` uses the midpoint between the two electrodes; the Cartesian midpoint is projected
# back onto the head surface (mean radius of the pair), the polar midpoint is taken between
# the 2-D topographic positions (so it follows whatever projection the polar coordinates use)
function _montage_loc(
    locs::DataFrame,
    lbl::String,
    a::String,
    b::Union{Nothing, String},
    loc::Symbol,
)::Union{Nothing, DataFrame}
    ra = _find_bylabel(locs, a)
    ra isa Int64 || return nothing
    row = locs[[ra], :]
    row[1, :label] = lbl
    (isnothing(b) || loc === :first) && return row
    rb = _find_bylabel(locs, b)
    rb isa Int64 || return nothing

    # Cartesian and spherical
    p1 = [Float64(locs[ra, c]) for c in (:loc_x, :loc_y, :loc_z)]
    p2 = [Float64(locs[rb, c]) for c in (:loc_x, :loc_y, :loc_z)]
    p = (p1 .+ p2) ./ 2
    np = norm(p)
    np > 0 && (p .*= ((norm(p1) + norm(p2)) / 2) / np)
    row[1, :loc_x], row[1, :loc_y], row[1, :loc_z] = p
    row[1, :loc_radius_sph], row[1, :loc_theta_sph], row[1, :loc_phi_sph] = cart2sph(p...)

    # polar
    x1, y1 = pol2cart(locs[ra, :loc_radius], locs[ra, :loc_theta])
    x2, y2 = pol2cart(locs[rb, :loc_radius], locs[rb, :loc_theta])
    row[1, :loc_radius], row[1, :loc_theta] = cart2pol((x1 + x2) / 2, (y1 + y2) / 2)

    return row
end

"""
    reference_custom(obj; <keyword arguments>)

Re-reference EEG channels using a custom bipolar or referential montage.

Each entry in `ref_list` is either a single channel name (copied as-is) or a `"Ch1-Ch2"` pair (output = Ch1 − Ch2). All EEG channels are replaced by the new montage channels; non-EEG channels are preserved.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`
- `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: montage definition; default is a standard longitudinal bipolar montage
- `ref_name::String="longitudinal-BIP"`: montage name stored in the header
- `loc::Symbol=:mid`: location of bipolar montage channels:
    - `:mid`: midpoint between the two electrodes
    - `:first`: location of the first electrode of the pair

# Returns

- `NeuroAnalyzer.NEURO`: new object with montage-referenced channels

# Notes

- `"Ch1-Ch2"` → output = amplitude(Ch1) − amplitude(Ch2).
- A montage channel is marked bad if any of its source channels was marked bad.
- Single-channel entries keep the location of their channel; for `"Ch1-Ch2"` entries the location is set according to `loc`. Locations of the replaced EEG channels are removed.
- With `loc=:mid` the Cartesian midpoint is projected back onto the head surface and the spherical coordinates are recalculated from it; the polar midpoint is taken between the 2-D topographic positions.
- A montage channel whose source channel(s) have no location gets no location.

Examples of montages:

- bipolar transverse: ["Fp2-Fp1", "F8-Fp2", "F8-F4", "F4-Fz", "Fz-F3", "F3-F7", "Fp1-F7", "T4-C4", "C4-Cz", "Cz-C3", "C3-T3", "T6-P4", "P4-Pz", "Pz-P3", "P3-T5", "O2-O1"], "transverse-BIP"
- bipolar longitudinal: ["Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["Fz", "Cz", "Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["FPz-Fz", "Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
"""
function reference_custom(
    obj::NeuroAnalyzer.NEURO;
    ref_list::Vector{String} = _BIP_LONGITUDINAL,
    ref_name::String = "longitudinal-BIP",
    loc::Symbol = :mid,
)::NeuroAnalyzer.NEURO
    _check_datatype(obj, "eeg")
    isempty(ref_list) && throw(ArgumentError("ref_list must not be empty."))
    _check_var(loc, [:mid, :first], "loc")

    ch = _eeg_ch_idx(obj)
    ch_lbls = labels(obj)[ch]
    lbl2idx = Dict(l => ch[i] for (i, l) in enumerate(ch_lbls))
    entries = [_parse_montage_entry(ref, ch_lbls) for ref in ref_list]

    n_new = length(ref_list)
    ep_n = nepochs(obj)
    rec0 = obj.header.recording
    s = zeros(n_new, epoch_len(obj), ep_n)
    bad_new = zeros(Bool, n_new)
    for (k, (a, b)) in enumerate(entries)
        ia = lbl2idx[a]
        if isnothing(b)
            @views s[k, :, :] .= obj.data[ia, :, :]
            bad_new[k] = rec0[:bad_channel][ia]
        else
            ib = lbl2idx[b]
            @views s[k, :, :] .= obj.data[ia, :, :] .- obj.data[ib, :, :]
            bad_new[k] = rec0[:bad_channel][ia] || rec0[:bad_channel][ib]
        end
    end

    # rebuild the object: montage channels first, then the non-EEG channels
    # (done directly instead of delete_channel(), which refuses to delete all channels)
    keep = setdiff(1:nchannels(obj), ch)
    obj_tmp = deepcopy(obj)
    rec = obj_tmp.header.recording
    obj_tmp.data = vcat(s, obj.data[keep, :, :])
    rec[:label] = vcat(ref_list, rec0[:label][keep])
    rec[:channel_type] = vcat(fill("eeg", n_new), rec0[:channel_type][keep])
    for f in (:unit, :prefiltering, :transducers, :gain)
        rec[f] = vcat(fill(rec0[f][ch[1]], n_new), rec0[f][keep])
    end
    rec[:bad_channel] = vcat(bad_new, collect(Bool, rec0[:bad_channel][keep]))
    rec[:channel_order] = _sort_channels(rec[:channel_type])
    rec[:reference] = ref_name

    # montage channel locations replace the EEG channel locations
    eeg_loc_rows = _find_bylabel(obj.locs, ch_lbls)
    new_locs = filter(
        !isnothing,
        [_montage_loc(obj.locs, ref_list[k], a, b, loc) for (k, (a, b)) in enumerate(entries)],
    )
    obj_tmp.locs = vcat(
        new_locs...,
        obj.locs[setdiff(1:size(obj.locs, 1), eeg_loc_rows), :];
        cols = :setequal,
    )

    push!(
        obj_tmp.history,
        "reference_custom(obj, ref_list=$ref_list, ref_name=$ref_name, loc=:$loc)",
    )

    return obj_tmp
end

"""
    reference_custom!(obj; <keyword arguments>)

Re-reference EEG channels using a custom montage in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`; modified in-place
- `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: montage definition; default is a standard longitudinal bipolar montage
- `ref_name::String="longitudinal-BIP"`: montage name stored in the header
- `loc::Symbol=:mid`: location of bipolar montage channels:
    - `:mid`: midpoint between the two electrodes
    - `:first`: location of the first electrode of the pair

# Returns

- `Nothing`

# Notes

- `"Ch1-Ch2"` → output = amplitude(Ch1) − amplitude(Ch2).
- A montage channel is marked bad if any of its source channels was marked bad.
- Single-channel entries keep the location of their channel; for `"Ch1-Ch2"` entries the location is set according to `loc`. Locations of the replaced EEG channels are removed.
- With `loc=:mid` the Cartesian midpoint is projected back onto the head surface and the spherical coordinates are recalculated from it; the polar midpoint is taken between the 2-D topographic positions.
- A montage channel whose source channel(s) have no location gets no location.

Examples of montages:

- bipolar transverse: ["Fp2-Fp1", "F8-Fp2", "F8-F4", "F4-Fz", "Fz-F3", "F3-F7", "Fp1-F7", "T4-C4", "C4-Cz", "Cz-C3", "C3-T3", "T6-P4", "P4-Pz", "Pz-P3", "P3-T5", "O2-O1"], "transverse-BIP"
- bipolar longitudinal: ["Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["Fz", "Cz", "Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["FPz-Fz", "Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
"""
function reference_custom!(
    obj::NeuroAnalyzer.NEURO;
    ref_list::Vector{String} = _BIP_LONGITUDINAL,
    ref_name::String = "longitudinal-BIP",
    loc::Symbol = :mid,
)::Nothing
    _copy_into!(obj, reference_custom(obj; ref_list = ref_list, ref_name = ref_name, loc = loc))
    return nothing
end
