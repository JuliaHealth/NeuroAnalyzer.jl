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

"""
    reference_ce(obj; <keyword arguments>)

Re-reference EEG channels to a common electrode or the average of multiple electrodes. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`
- `ch::Union{String, Vector{String}, Regex}`: reference channel name(s); if multiple channels are given their mean (or median) is used as the reference
- `med::Bool=false`: if `true`, use median instead of mean for multi-channel reference

# Returns

- `NeuroAnalyzer.NEURO`: new object with re-referenced EEG channels

# See also

[`reference_ce!`](@ref), [`reference_avg`](@ref)
"""
function reference_ce(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    med::Bool = false
)::NeuroAnalyzer.NEURO

    _check_datatype(obj, "eeg")

    # reference channel indices
    ref_ch_idx  = get_channel(obj, ch=ch)
    # channels to re-reference
    sig_ch_idx  = get_channel(obj, ch=get_channel(obj; type="eeg"))

    # number of channels
    ch_n = length(sig_ch_idx)
    # number of epochs
    ep_n = nepochs(obj)

    obj_new = deepcopy(obj)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]

        s_ch  = sig_ch_idx[ch_idx]

        # skip when the current channel IS the reference (single-channel case)
        length(ref_ch_idx) == 1 && s_ch == ref_ch_idx[1] && continue

        ref = if length(ref_ch_idx) == 1
            @view obj.data[ref_ch_idx[1], :, ep_idx]
        else
            src = @view obj.data[ref_ch_idx, :, ep_idx]
            med ? vec(median(src; dims=1)) : vec(mean(src; dims=1))
        end
        obj_new.data[s_ch, :, ep_idx] .-= ref
    end

    suffix = length(ref_ch_idx) == 1 ? "-$(labels(obj)[ref_ch_idx[1]])" : "-cavg"
    obj_new.header.recording[:label][sig_ch_idx] .*= suffix
    obj_new.locs[:, :label][sig_ch_idx] .*= suffix

    obj_new.header.recording[:reference] = if length(ref_ch_idx) == 1
        "common ($(labels(obj)[ref_ch_idx[1]]))"
    else
        "common ($(join(labels(obj)[ref_ch_idx], ", ")) averaged)"
    end
    push!(obj_new.history, "reference_ce(OBJ, ch=$ch, med=$med)")

    return obj_new

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

# See also

[`reference_ce`](@ref), [`reference_avg`](@ref)
"""
function reference_ce!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    med::Bool = false
)::Nothing

    obj_new = reference_ce(obj, ch = ch, med = med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

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
- `med::Bool=false`: use median instead of mean
- `weighted::Bool=false`: weight reference channels by inverse distance (requires channel locations)

# Returns

- `NeuroAnalyzer.NEURO`: new object with re-referenced EEG channels

# See also

[`reference_avg!`](@ref), [`reference_ce`](@ref)
"""
function reference_avg(
    obj::NeuroAnalyzer.NEURO;
    exclude_fpo::Bool = false,
    exclude_current::Bool = false,
    average::Bool = true,
    med::Bool = false,
    weighted::Bool = false
)::NeuroAnalyzer.NEURO

    _check_datatype(obj, "eeg")

    # channels that will be referenced
    sig_ch_idx = get_channel(obj, ch=get_channel(obj; type="eeg"))
    ch_n       = length(sig_ch_idx)
    ep_n       = nepochs(obj)

    # pre-compute channel locations and distance matrix when weighted
    if weighted
        _has_locs(obj)
        chs = intersect(obj.locs[!, :label], labels(obj)[sig_ch_idx])
        locs = Base.filter(:label => in(chs), obj.locs)
        _check_ch_locs(sig_ch_idx, labels(obj), obj.locs[!, :label])
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
    end

    # pre-compute indices of Fp/O channels to optionally exclude
    fpo_idx = Int64[]
    if exclude_fpo
        l = lowercase.(labels(obj))
        for lbl in ("fp1", "fp2", "o1", "o2")
            i = findfirst(isequal(lbl), l)
            !isnothing(i) && push!(fpo_idx, i)
        end
    end

    # read-only source
    src = obj.data[sig_ch_idx, :, :]
    # output buffer
    dst = copy(src)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]

        chs2exclude = copy(fpo_idx)
        exclude_current && push!(chs2exclude, ch_idx)
        ref_rows = setdiff(1:ch_n, unique(chs2exclude))

        if weighted
            w = zeros(ch_n)
            for w_idx in 1:ch_n
                w[w_idx] = euclidean(
                    [loc_x[ch_idx], loc_y[ch_idx]],
                    [loc_x[w_idx],  loc_y[w_idx]],
                )
            end
            w = 1 .- normalize_n(w)
            ref_chs = @view src[ref_rows, :, ep_idx]
            w_sub   = w[ref_rows]
            ref_ch  = med ? vec(median(w_sub .* ref_chs; dims=1)) :
                             vec(mean(w_sub  .* ref_chs; dims=1))
        else
            ref_chs = @view src[ref_rows, :, ep_idx]
            ref_ch  = if average
                med ? vec(median(ref_chs; dims=1)) : vec(mean(ref_chs; dims=1))
            else
                vec(sum(ref_chs; dims=1))
            end
        end

        @inbounds dst[ch_idx, :, ep_idx] = src[ch_idx, :, ep_idx] .- ref_ch
    end

    obj_new = deepcopy(obj)
    obj_new.data[sig_ch_idx, :, :] = dst

    suffix = average ? (weighted ? "-wavg" : "-avg") : (weighted ? "-wsum" : "-sum")
    obj_new.header.recording[:label][sig_ch_idx] .*= suffix
    ch_locs = _find_bylabel(obj.locs, labels(obj)[sig_ch_idx])
    obj_new.locs[ch_locs, :label] .*= suffix

    obj_new.header.recording[:reference] = if average
        weighted ? "average (weighted)" : "average"
    else
        weighted ? "sum (weighted)" : "sum"
    end
    push!(
        obj_new.history,
        "reference_avg(OBJ, exclude_fpo=$exclude_fpo, exclude_current=$exclude_current, average=$average, med=$med, weighted=$weighted)",
    )

    return obj_new

end

"""
    reference_avg!(obj; <keyword arguments>)

Re-reference EEG channels to the common average reference in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`; modified in-place
- `exclude_fpo::Bool=false`: exclude Fp1, Fp2 (due to eye blinks), O1, O2 (due to head movements) from the CAR calculation
- `exclude_current::Bool=false`: exclude the current channel from its own CAR calculation
- `average::Bool=true`: subtract the mean (or median) reference; if `false`, subtract the sum
- `med::Bool=false`: use median instead of mean
- `weighted::Bool=false`: weight reference channels by inverse distance (requires channel locations)
# Returns

- `Nothing`

# See also

[`reference_avg`](@ref), [`reference_ce`](@ref)
"""
function reference_avg!(
    obj::NeuroAnalyzer.NEURO;
    exclude_fpo::Bool = false,
    exclude_current::Bool = false,
    average::Bool = true,
    med::Bool = false,
    weighted::Bool = false
)::Nothing

    obj_new = reference_avg(
        obj,
        exclude_fpo = exclude_fpo,
        exclude_current = exclude_current,
        average = average,
        med = med,
        weighted = weighted,
    )
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end

# ---------------------------------------------------------------------------
# internal helper: apply auricular/mastoid reference for one electrode pair.
# used by reference_a and reference_m to avoid code duplication.
# ---------------------------------------------------------------------------
function _apply_paired_reference!(
    s_ref::AbstractArray,
    s::AbstractArray,
    ref_data::AbstractArray,
    picks::Vector{Int64},
    ref_label::Vector{String},
    suffix::String,
    ep_n::Int64,
    med::Bool,
)
    for ep_idx in 1:ep_n
        ref_ch = med ? vec(median(ref_data[:, :, ep_idx]; dims=1)) :
                       vec(mean(ref_data[:, :, ep_idx];   dims=1))
        # thread over channels for this epoch
        Threads.@threads :dynamic for ch_idx in picks
            @inbounds s_ref[ch_idx, :, ep_idx] = s[ch_idx, :, ep_idx] .- ref_ch
            ref_label[ch_idx] = suffix
        end
    end
end

function _apply_single_reference!(
    s_ref::AbstractArray,
    s::AbstractArray,
    ref_data::AbstractArray,
    picks::Vector{Int64},
    ref_label::Vector{String},
    suffix::String,
    ep_n::Int64,
)
    for ep_idx in 1:ep_n
        ref_ch = vec(ref_data[:, :, ep_idx])
        Threads.@threads :dynamic for ch_idx in picks
            @inbounds s_ref[ch_idx, :, ep_idx] = s[ch_idx, :, ep_idx] .- ref_ch
            ref_label[ch_idx] = suffix
        end
    end
end

"""
    reference_a(obj; <keyword arguments>)

Re-reference EEG channels to auricular electrodes (A1, A2). Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must contain A1 and A2 channels; must be of type `"eeg"`
- `type::Symbol=:l`: reference type:
    - `:l`: linked — average of A1 and A2
    - `:i`: ipsilateral — A1 for left channels, A2 for right channels
    - `:c`: contralateral — A2 for left channels, A1 for right channels
- `med::Bool=false`: use median instead of mean

# Returns

- `NeuroAnalyzer.NEURO`: new re-referenced object

# See also

[`reference_a!`](@ref), [`reference_m`](@ref)
"""
function reference_a(
    obj::NeuroAnalyzer.NEURO;
    type::Symbol = :l,
    med::Bool = false
)::NeuroAnalyzer.NEURO

    _check_datatype(obj, "eeg")
    _check_var(type, [:l, :i, :c], "type")
    "A1" in labels(obj) || throw(ArgumentError("OBJ does not contain A1 channel."))
    "A2" in labels(obj) || throw(ArgumentError("OBJ does not contain A2 channel."))

    ch = get_channel(obj, ch=get_channel(obj; type="eeg"))
    obj_new = deepcopy(obj)
    s = obj_new.data[ch, :, :]
    a1 = extract_channel(obj, ch="A1")
    a2 = extract_channel(obj, ch="A2")
    ch_n = size(s, 1)
    ep_n = size(s, 3)
    s_ref = copy(s)
    ref_label = fill("", length(ch))

    linked = vcat(a1, a2)

    if type === :l
        _apply_paired_reference!(s_ref, s, linked, collect(1:ch_n), ref_label, "-A1A2", ep_n, med)
    elseif type === :i
        _apply_paired_reference!(s_ref, s, linked,
            get_channel(obj, ch=channel_pick(obj, pick=:central)), ref_label, "-A1A2", ep_n, med)
        _apply_single_reference!(s_ref, s, a1,
            get_channel(obj, ch=channel_pick(obj, pick=:left)),  ref_label, "-A1", ep_n)
        _apply_single_reference!(s_ref, s, a2,
            get_channel(obj, ch=channel_pick(obj, pick=:right)), ref_label, "-A2", ep_n)
    elseif type === :c
        _apply_paired_reference!(s_ref, s, linked,
            get_channel(obj, ch=channel_pick(obj, pick=:central)), ref_label, "-A1A2", ep_n, med)
        _apply_single_reference!(s_ref, s, a2,
            get_channel(obj, ch=channel_pick(obj, pick=:left)),  ref_label, "-A2", ep_n)
        _apply_single_reference!(s_ref, s, a1,
            get_channel(obj, ch=channel_pick(obj, pick=:right)), ref_label, "-A1", ep_n)
    end

    obj_new.data[ch, :, :] = s_ref
    obj_new.header.recording[:label][ch] .*= ref_label
    ch_locs = _find_bylabel(obj.locs, labels(obj)[ch])
    obj_new.locs[ch_locs, :label] .*= ref_label
    obj_new.header.recording[:reference] = Dict(
        :l => "auricular (linked)",
        :i => "auricular (ipsilateral)",
        :c => "auricular (contralateral)",
    )[type]
    push!(obj_new.history, "reference_a(OBJ, type=$type, med=$med)")

    return obj_new

end

"""
    reference_a!(obj; <keyword arguments>)

Re-reference EEG channels to auricular electrodes in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must contain A1 and A2 channels; must be of type `"eeg"`; modified in-place
- `type::Symbol=:l`: reference type:
    - `:l`: linked — average of A1 and A2
    - `:i`: ipsilateral — A1 for left channels, A2 for right channels
    - `:c`: contralateral — A2 for left channels, A1 for right channels
- `med::Bool=false`: use median instead of mean

# Returns

- `Nothing`

# See also

[`reference_a`](@ref), [`reference_m`](@ref)
"""
function reference_a!(
    obj::NeuroAnalyzer.NEURO;
    type::Symbol = :l,
    med::Bool = false
)::Nothing

    obj_new = reference_a(obj, type = type, med = med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end

"""
    reference_m(obj; <keyword arguments>)

Re-reference EEG channels to mastoid electrodes (M1, M2). Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must contain M1 and M2 channels; must be of type `"eeg"`; modified in-place
- `type::Symbol=:l`: reference type:
    - `:l`: linked — average of M1 and M2
    - `:i`: ipsilateral — M1 for left channels, M2 for right channels
    - `:c`: contralateral — M2 for left channels, M1 for right channels
- `med::Bool=false`: use median instead of mean

# Returns

- `NeuroAnalyzer.NEURO`: new re-referenced object

# See also

[`reference_m!`](@ref), [`reference_a`](@ref)
"""
function reference_m(obj::NeuroAnalyzer.NEURO; type::Symbol = :l, med::Bool = false)::NeuroAnalyzer.NEURO

    _check_datatype(obj, "eeg")
    _check_var(type, [:l, :i, :c], "type")
    "M1" in labels(obj) || throw(ArgumentError("OBJ does not contain M1 channel."))
    "M2" in labels(obj) || throw(ArgumentError("OBJ does not contain M2 channel."))

    ch = get_channel(obj; ch=get_channel(obj; type="eeg"))
    obj_new = deepcopy(obj)
    s = obj_new.data[ch, :, :]
    m1 = extract_channel(obj; ch="M1")
    m2 = extract_channel(obj; ch="M2")
    ch_n = size(s, 1)
    ep_n = size(s, 3)
    s_ref = copy(s)
    ref_label = fill("", length(ch))

    linked = vcat(m1, m2)

    if type === :l
        _apply_paired_reference!(s_ref, s, linked, collect(1:ch_n), ref_label, "-M1M2", ep_n, med)
    elseif type === :i
        _apply_paired_reference!(s_ref, s, linked,
            get_channel(obj; ch=channel_pick(obj, pick=:central)), ref_label, "-M1M2", ep_n, med)
        _apply_single_reference!(s_ref, s, m1,
            get_channel(obj; ch=channel_pick(obj, pick=:left)),  ref_label, "-M1", ep_n)
        _apply_single_reference!(s_ref, s, m2,
            get_channel(obj; ch=channel_pick(obj, pick=:right)), ref_label, "-M2", ep_n)
    elseif type === :c
        _apply_paired_reference!(s_ref, s, linked,
            get_channel(obj; ch=channel_pick(obj, pick=:central)), ref_label, "-M1M2", ep_n, med)
        _apply_single_reference!(s_ref, s, m2,
            get_channel(obj; ch=channel_pick(obj, pick=:left)),  ref_label, "-M2", ep_n)
        _apply_single_reference!(s_ref, s, m1,
            get_channel(obj; ch=channel_pick(obj, pick=:right)), ref_label, "-M1", ep_n)
    end

    obj_new.data[ch, :, :] = s_ref
    obj_new.header.recording[:label][ch] .*= ref_label
    ch_locs = _find_bylabel(obj.locs, labels(obj)[ch])
    obj_new.locs[ch_locs, :label] .*= ref_label
    obj_new.header.recording[:reference] = Dict(
        :l => "mastoid (linked)",
        :i => "mastoid (ipsilateral)",
        :c => "mastoid (contralateral)",
    )[type]
    push!(obj_new.history, "reference_m(OBJ, type=$type, med=$med)")

    return obj_new

end

"""
    reference_m!(obj; <keyword arguments>)

Re-reference EEG channels to mastoid electrodes in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must contain M1 and M2 channels; must be of type `"eeg"`; modified in-place
- `type::Symbol=:l`: reference type:
    - `:l`: linked — average of M1 and M2
    - `:i`: ipsilateral — M1 for left channels, M2 for right channels
    - `:c`: contralateral — M2 for left channels, M1 for right channels
- `med::Bool=false`: use median instead of mean

# Returns

- `Nothing`
"""
function reference_m!(
    obj::NeuroAnalyzer.NEURO;
    type::Symbol = :l,
    med::Bool = false
)::Nothing

    obj_new = reference_m(obj, type = type, med = med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end

# ---------------------------------------------------------------------------
# internal helper: shared nearest-neighbour Laplacian core
# used by both reference_plap and reference_slap
# ---------------------------------------------------------------------------
function _laplacian_reference(
    obj::NeuroAnalyzer.NEURO,
    d::Matrix{Float64},
    nn::Int64,
    weighted::Bool,
    med::Bool,
    loc_x::Vector{Float64},
    loc_y::Vector{Float64},
    loc_z::Union{Nothing, Vector{Float64}}=nothing
)::Matrix{Float64}

    ch = get_channel(obj; ch=get_channel(obj; type="eeg"))
    s = @view obj.data[ch, :, :]
    ch_n, ep_n = size(s, 1), size(s, 3)

    nn_idx = zeros(Int64, ch_n, nn)
    for idx in 1:ch_n
        nn_idx[idx, :] = sortperm(d[idx, :])[1:nn]
    end

    s_ref = zeros(ch_n, size(s, 2), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ref_chs = @view s[nn_idx[ch_idx, :], :, ep_idx]

        if weighted
            w = zeros(nn)
            for w_idx in 1:nn
                ni = nn_idx[ch_idx, w_idx]
                w[w_idx] = isnothing(loc_z) ?
                    euclidean([loc_x[ch_idx], loc_y[ch_idx]], [loc_x[ni], loc_y[ni]]) :
                    _sph_distance_cart(loc_x[ch_idx], loc_y[ch_idx], loc_z[ch_idx],
                                       loc_x[ni],     loc_y[ni],     loc_z[ni])
            end
            w       = 1 .- normalize_n(w)
            ref_ch  = med ? vec(median(w .* ref_chs; dims=1)) :
                             vec(mean(w  .* ref_chs; dims=1))
        else
            ref_ch = med ? vec(median(ref_chs; dims=1)) : vec(mean(ref_chs; dims=1))
        end
        @inbounds s_ref[ch_idx, :, ep_idx] = s[ch_idx, :, ep_idx] .- ref_ch
    end
    return s_ref
end

"""
    reference_plap(obj; <keyword arguments>)

Re-reference EEG channels using the planar (2-D Euclidean) Laplacian. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must have channel locations; must be of type `"eeg"`
- `nn::Int64=4`: number of nearest-neighbour electrodes; must satisfy `1 ≤ nn < nch-1`
- `weighted::Bool=false`: weight neighbours by inverse distance
- `med::Bool=false`: use median instead of mean.

# Returns
- `NeuroAnalyzer.NEURO`: new re-referenced object

# See also

[`reference_plap!`](@ref), [`reference_slap`](@ref)
"""
function reference_plap(
    obj::NeuroAnalyzer.NEURO;
    nn::Int64 = 4,
    weighted::Bool = false,
    med::Bool = false
)::NeuroAnalyzer.NEURO

    _check_datatype(obj, "eeg")
    _has_locs(obj)

    ch = get_channel(obj; ch=get_channel(obj; type="eeg"))
    chs = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
    ch_n = length(ch)

    nn >= 1 || throw(ArgumentError("nn must be ≥ 1."))
    nn < ch_n - 1 || throw(ArgumentError("nn must be < $(ch_n - 1)."))

    loc_x = locs[!, :loc_x]
    loc_y = locs[!, :loc_y]

    # Euclidean distances
    d = [euclidean([loc_x[i], loc_y[i]], [loc_x[j], loc_y[j]]) for i in 1:ch_n, j in 1:ch_n]
    # eliminate auto-referencing
    d[d .== 0] .= Inf

    s_ref = _laplacian_reference(obj, d, nn, weighted, med, loc_x, loc_y)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = s_ref
    suffix = weighted ? "-wplap" : "-plap"
    obj_new.header.recording[:label][ch] .*= suffix
    ch_locs = _find_bylabel(obj.locs, labels(obj)[ch])
    obj_new.locs[ch_locs, :label] .*= suffix
    obj_new.header.recording[:reference] = weighted ? "weighted planar Laplacian ($nn)" : "planar Laplacian ($nn)"
    push!(obj_new.history, "reference_plap(OBJ, nn=$nn, weighted=$weighted, med=$med)")

    return obj_new

end

"""
    reference_plap!(obj; <keyword arguments>)

Re-reference EEG channels using the planar Laplacian in-place. Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must have channel locations; must be of type `"eeg"`; modified in-place
- `nn::Int64=4`: number of nearest-neighbour electrodes; must satisfy `1 ≤ nn < nch-1`
- `weighted::Bool=false`: weight neighbours by inverse distance
- `med::Bool=false`: use median instead of mean.

# Returns

- `Nothing`

# See also

[`reference_plap`](@ref), [`reference_slap`](@ref)
"""
function reference_plap!(
    obj::NeuroAnalyzer.NEURO;
    nn::Int64 = 4,
    weighted::Bool = false,
    med::Bool = false
)::Nothing

    obj_new = reference_plap(obj, nn = nn, weighted = weighted, med = med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end


"""
    reference_slap(obj; <keyword arguments>)

Re-reference EEG channels using the spherical Laplacian (great-circle distance). Only EEG-type channels are modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must have channel locations; must be of type `"eeg"`
- `nn::Int64=4`: number of nearest-neighbour electrodes; must satisfy `1 ≤ nn < nch-1`
- `weighted::Bool=false`: weight neighbours by inverse distance
- `med::Bool=false`: use median instead of mean.

# Returns

- `NeuroAnalyzer.NEURO`: new re-referenced object

# See also

[`reference_slap!`](@ref), [`reference_plap`](@ref)
"""
function reference_slap(
    obj::NeuroAnalyzer.NEURO;
    nn::Int64 = 4,
    weighted::Bool = false,
    med::Bool = false
)::NeuroAnalyzer.NEURO

    _check_datatype(obj, "eeg")
    _has_locs(obj)

    ch = get_channel(obj; ch=get_channel(obj; type="eeg"))
    chs = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
    ch_n = length(ch)

    nn >= 1 || throw(ArgumentError("nn must be ≥ 1."))
    nn < ch_n - 1 || throw(ArgumentError("nn must be < $(ch_n - 1)."))

    loc_x = locs[!, :loc_x]
    loc_y = locs[!, :loc_y]
    loc_z = locs[!, :loc_z]

    # Euclidean distance
    d = [_sph_distance_cart(loc_x[i], loc_y[i], loc_z[i], loc_x[j], loc_y[j], loc_z[j])
         for i in 1:ch_n, j in 1:ch_n]
    # eliminate auto-referencing
    d[d .== 0] .= Inf

    s_ref   = _laplacian_reference(obj, d, nn, weighted, med, loc_x, loc_y, loc_z)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = s_ref
    suffix  = weighted ? "-wslap" : "-slap"
    obj_new.header.recording[:label][ch] .*= suffix
    ch_locs = _find_bylabel(obj.locs, labels(obj)[ch])
    obj_new.locs[ch_locs, :label] .*= suffix
    obj_new.header.recording[:reference] = weighted ? "weighted spherical Laplacian ($nn)" : "spherical Laplacian ($nn)"
    push!(obj_new.history, "reference_slap(OBJ, nn=$nn, weighted=$weighted, med=$med)")

    return obj_new

end

"""
    reference_slap!(obj; <keyword arguments>)

Re-reference EEG channels using the spherical Laplacian in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must have channel locations; must be of type `"eeg"`; modified in-place
- `nn::Int64=4`: number of nearest-neighbour electrodes; must satisfy `1 ≤ nn < nch-1`
- `weighted::Bool=false`: weight neighbours by inverse distance
- `med::Bool=false`: use median instead of mean.

# Returns

- `Nothing`

# See also

[`reference_slap`](@ref), [`reference_plap`](@ref)
"""
function reference_slap!(
    obj::NeuroAnalyzer.NEURO;
    nn::Int64 = 4,
    weighted::Bool = false,
    med::Bool = false
)::Nothing

    obj_new = reference_slap(obj, nn = nn, weighted = weighted, med = med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end

"""
    reference_custom(obj; <keyword arguments>)

Re-reference EEG channels using a custom bipolar or referential montage.

Each entry in `ref_list` is either a single channel name (copied as-is) or a `"Ch1-Ch2"` pair (output = Ch2 − Ch1). All referenced EEG channels are replaced by the new montage channels; non-EEG channels are preserved.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: montage definition; default is a standard longitudinal bipolar montage
- `ref_name::String="longitudinal-BIP"`: montage name stored in the header

# Returns

- `NeuroAnalyzer.NEURO`: new object with montage-referenced channels

# Notes

- `"Ch1-Ch2"` → output = amplitude(Ch2) − amplitude(Ch1).
- The channel location table is not updated (marked as TODO).

Examples of montages:

- bipolar transverse: ["Fp2-Fp1", "F8-Fp2", "F8-F4", "F4-Fz", "Fz-F3", "F3-F7", "Fp1-F7", "T4-C4", "C4-Cz", "Cz-C3", "C3-T3", "T6-P4", "P4-Pz", "Pz-P3", "P3-T5", "O2-O1"], "transverse-BIP"
- bipolar longitudinal: ["Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["Fz", "Cz", "Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["FPz-Fz", "Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
"""
function reference_custom(
    obj::NeuroAnalyzer.NEURO;
    ref_list::Vector{String} = [
        "Fz-Cz","Cz-Pz","Fp1-F7","Fp1-F3","F7-T3","T3-T5","T5-O1",
        "F3-C3","C3-P3","P3-O1","Fp2-F8","F8-T4","T4-T6","T6-O2",
        "Fp2-F4","F4-C4","C4-P4","P4-O2"
    ],
    ref_name::String = "longitudinal-BIP"
)::NeuroAnalyzer.NEURO

    _check_datatype(obj, "eeg")

    ch = get_channel(obj; ch=get_channel(obj; type="eeg"))
    ch_lbls = labels(obj)[ch]

    for ref in ref_list
        if '-' in ref
            m = match(r"(.+)-(.+)", ref)
            m[1] in ch_lbls || throw(ArgumentError("Label $(m[1]) not found in OBJ."))
            m[2] in ch_lbls || throw(ArgumentError("Label $(m[2]) not found in OBJ."))
        else
            ref in ch_lbls || throw(ArgumentError("Label $ref not found in OBJ."))
        end
    end

    ep_n = nepochs(obj)
    s = zeros(length(ref_list), epoch_len(obj), ep_n)

    @inbounds for ep_idx in 1:ep_n
        for (ref_idx, ref) in enumerate(ref_list)
            if '-' in ref
                m = match(r"(.+)-(.+)", ref)
                idx1 = get_channel(obj; ch=string(m[1]))
                idx2 = get_channel(obj; ch=string(m[2]))
                s[ref_idx, :, ep_idx] = obj.data[idx1, :, ep_idx] .- obj.data[idx2, :, ep_idx]
            else
                s[ref_idx, :, ep_idx] = obj.data[get_channel(obj; ch=ref), :, ep_idx]
            end
        end
    end

    obj_new = delete_channel(obj; ch=get_channel(obj; type="eeg"))
    obj_new.data = vcat(s, obj_new.data)
    rec = obj_new.header.recording
    rec[:label] = vcat(ref_list, labels(obj_new))
    rec[:reference] = ref_name
    rec[:channel_type] = vcat(fill("eeg", length(ref_list)), rec[:channel_type])
    rec[:unit] = vcat(fill("μV",  length(ref_list)), rec[:unit])
    rec[:prefiltering] = vcat(fill(obj.header.recording[:prefiltering][1], length(ref_list)), rec[:prefiltering])
    rec[:transducers]  = vcat(fill(obj.header.recording[:transducers][1],  length(ref_list)), rec[:transducers])
    rec[:gain] = vcat(fill(obj.header.recording[:gain][1],         length(ref_list)), rec[:gain])
    _info("Bad channels matrix will be reset.")
    rec[:bad_channel]  = falses(size(obj_new.data, 1))   # was: zeros — Bool is more appropriate

    # TODO: update obj_new.locs for the new montage channels

    push!(obj_new.history, "reference_custom(OBJ, ref_list=$ref_list, ref_name=$ref_name)")

    return obj_new

end

"""
    reference_custom!(obj; <keyword arguments>)

Re-reference EEG channels using a custom montage in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be of type `"eeg"`; modified in-place
- `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: montage definition; default is a standard longitudinal bipolar montage
- `ref_name::String="longitudinal-BIP"`: montage name stored in the header

# Returns

- `Nothing`

# Notes

- `"Ch1-Ch2"` → output = amplitude(Ch2) − amplitude(Ch1).
- The channel location table is not updated (marked as TODO).

Examples of montages:

- bipolar transverse: ["Fp2-Fp1", "F8-Fp2", "F8-F4", "F4-Fz", "Fz-F3", "F3-F7", "Fp1-F7", "T4-C4", "C4-Cz", "Cz-C3", "C3-T3", "T6-P4", "P4-Pz", "Pz-P3", "P3-T5", "O2-O1"], "transverse-BIP"
- bipolar longitudinal: ["Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["Fz", "Cz", "Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["FPz-Fz", "Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
"""
function reference_custom!(
    obj::NeuroAnalyzer.NEURO;
    ref_list::Vector{String}=[
        "Fz-Cz","Cz-Pz","Fp1-F7","Fp1-F3","F7-T3","T3-T5","T5-O1",
        "F3-C3","C3-P3","P3-O1","Fp2-F8","Fp2-F4","F8-T4","T4-T6","T6-O2",
        "F4-C4","C4-P4","P4-O2",
    ],
    ref_name::String="longitudinal-BIP"
)::Nothing

    obj_new = reference_custom(obj, ref_list = ref_list, ref_name = ref_name)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end
