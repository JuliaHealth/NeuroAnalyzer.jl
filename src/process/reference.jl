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
export reference_custom
export reference_custom!

"""
    reference_ce(obj; <keyword arguments>)

Reference to common electrode(s). Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels used as reference; if multiple channels are specified, their average is used as the reference
- `med::Bool=false`: use median instead of mean

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function reference_ce(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, med::Bool=false)

    _check_datatype(obj, "eeg")

    # keep signal channels
    ch = get_channel(obj, ch=ch)
    chs = [get_channel(obj, ch=idx)[1] for idx in get_channel(obj, type="eeg")]

    obj_new = deepcopy(obj)

    s = @view obj_new.data[chs, :, :]

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if ch isa Int64
                ref_ch = @views vec(s[ch, :, ep_idx])
                if ch_idx != ch
                    @views s[ch_idx, :, ep_idx] .-= ref_ch
                end
            else
                if !med
                    ref_ch = @views vec(mean(s[ch, :, ep_idx], dims=1))
                else
                    ref_ch = @views vec(median(s[ch, :, ep_idx], dims=1))
                end
                @views s[ch_idx, :, ep_idx] .-= ref_ch
            end
        end
    end

    obj_new.header.recording[:label][chs] .*= ch isa Int64 ? "-$(labels(obj)[ch])" : "-cavg"
    obj_new.data[chs, :, :] = s
    if ch isa Int64
        obj_new.header.recording[:reference] = "common ($(labels(obj)[ch]))"
    else
        obj_new.header.recording[:reference] = "common ($(labels(obj)[ch]) averaged)"
    end
    reset_components!(obj_new)
    push!(obj_new.history, "reference_ce(OBJ, ch=$ch, med=$med)")

    return obj_new

end

"""
    reference_ce!(obj; <keyword arguments>)

Reference to common electrode(s). Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels used as reference; if multiple channels are specified, their average is used as the reference
- `med::Bool=false`: use median instead of mean
"""
function reference_ce!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, med::Bool=false)

    obj_new = reference_ce(obj, ch=ch, med=med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    reference_avg(obj; <keyword arguments>)

Reference to averaged reference. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `exclude_fpo::Bool=false`: exclude Fp1, Fp2 (due to eye blinks), O1, O2 (due to head movements) from CAR calculation
- `exclude_current::Bool=false`: exclude current channel from CAR calculation
- `average::Bool=true`: average reference channels prior to subtracting, otherwise add all reference channels
- `med::Bool=false`: use median instead of mean
- `weighted::Bool=false`: use weighted reference channels (weights depend on the distance from the current electrode)

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function reference_avg(obj::NeuroAnalyzer.NEURO; exclude_fpo::Bool=false, exclude_current::Bool=false, average::Bool=true, med::Bool=false, weighted::Bool=false)

    _check_datatype(obj, "eeg")

    # keep signal channels
    chs = [get_channel(obj, ch=idx)[1] for idx in get_channel(obj, type="eeg")]

    # source signals
    obj_new = deepcopy(obj)
    src = @view deepcopy(obj_new).data[chs, :, :]
    ch_n = size(src, 1)
    ep_n = size(src, 3)
    # destination signals
    dst = @view deepcopy(obj_new).data[chs, :, :]

    if weighted
        @assert length(chs) == nrow(obj.locs) "Some channels do not have locations."
        _has_locs(obj)

        locs_idx = _find_bylabel(obj.locs, labels(obj)[chs])
        loc_x = obj.locs[locs_idx, :loc_x]
        loc_y = obj.locs[locs_idx, :loc_y]

        # Euclidean distance matrix
        d = zeros(ch_n, ch_n)
        for idx1 in 1:ch_n, idx2 in 1:ch_n
                d[idx1, idx2] = euclidean([loc_x[idx1], loc_y[idx1]], [loc_x[idx2], loc_y[idx2]])
        end
        # set weights not to reference to itself
        d[d .== 0] .= Inf
        w = sum(1 ./ d)
    end

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n

            if weighted
                src = @view deepcopy(obj).data[chs, :, :]
                w = zeros(ch_n)
                # calculate vector of weights - distances between the current electrode and each of the reference electrodes
                for w_idx in 1:ch_n
                    w[w_idx] = euclidean([loc_x[ch_idx], loc_y[ch_idx]], [loc_x[w_idx], loc_y[w_idx]])
                end
                # invert distances, so that closer electrodes have higher weights
                w = 1 .- normalize_n(w)
                # apply weights
                src .*= w
            end

            chs2exclude = Vector{Int64}()
            if exclude_fpo
                l = lowercase.(labels(obj))
                "fp1" in l && push!(chs2exclude, findfirst(isequal("fp1"), l))
                "fp2" in l && push!(chs2exclude, findfirst(isequal("fp2"), l))
                "o1" in l && push!(chs2exclude, findfirst(isequal("o1"), l))
                "o2" in l && push!(chs2exclude, findfirst(isequal("o2"), l))
            end
            exclude_current && push!(chs2exclude, ch_idx)
            ref_chs = @view src[setdiff(1:ch_n, unique(chs2exclude)), :, ep_idx]

            if average
                if !med
                    ref_ch = vec(mean(ref_chs, dims=1))
                else
                    ref_ch = vec(median(ref_chs, dims=1))
                end
            else
                ref_ch = vec(sum(ref_chs, dims=1))
            end
            dst[ch_idx, :, ep_idx] .-= ref_ch
        end
    end

    if average
        obj_new.header.recording[:label][chs] .*= weighted ? "-wavg" : "-avg"
    else
        obj_new.header.recording[:label][chs] .*= weighted ? "-wsum" : "-sum"
    end

    obj_new.data[chs, :, :] = dst
    if average
        obj_new.header.recording[:reference] = weighted ? "average (weighted)" : "average"
    else
        obj_new.header.recording[:reference] = weighted ? "sum (weighted)" : "sum"
    end
    reset_components!(obj_new)
    push!(obj_new.history, "reference_avg(OBJ, exclude_fpo=$exclude_fpo, exclude_current=$exclude_current, average=$average, med=$med, weighted=$weighted)")

    return obj_new

end

"""
    reference_avg!(obj; <keyword arguments>)

Reference to averaged reference. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `exclude_fpo::Bool=false`: exclude Fp1, Fp2 (due to eye blinks), O1, O2 (due to head movements) from CAR calculation
- `exclude_current::Bool=false`: exclude current channel from CAR mean calculation
- `average::Bool=true`: average reference channels prior to subtracting, otherwise add all reference channels
- `med::Bool=false`: use median instead of mean
- `weighted::Bool=false`: use weighted reference channels (weights depend on the distance from the current electrode)
"""
function reference_avg!(obj::NeuroAnalyzer.NEURO; exclude_fpo::Bool=false, exclude_current::Bool=false, average::Bool=true, med::Bool=false, weighted::Bool=false)

    obj_new = reference_avg(obj, exclude_fpo=exclude_fpo, exclude_current=exclude_current, average=average, med=med, weighted=weighted)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    reference_a(obj; <keyword arguments>)

Reference to auricular (A1, A2) channels. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:l`:
    - `:l`: linked - average of A1 and A2
    - `:i`: ipsilateral - A1 for left channels, A2 for right channels
    - `:c`: contraletral - A1 for right channels, A2 for left channels
- `med::Bool=false`: use median instead of mean

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function reference_a(obj::NeuroAnalyzer.NEURO; type::Symbol=:l, med::Bool=false)

    _check_datatype(obj, "eeg")
    _check_var(type, [:l, :i, :c], "type")
    @assert "a1" in lowercase.(labels(obj)) "OBJ does not contain A1 channel."
    @assert "a2" in lowercase.(labels(obj)) "OBJ does not contain A2 channel."

    # keep signal channels
    chs = [get_channel(obj, ch=idx)[1] for idx in get_channel(obj, type="eeg")]
    obj_new = deepcopy(obj)
    s = @view obj_new.data[chs, :, :]

    a1_idx = labels(obj)[findfirst(isequal("a1"), lowercase.(labels(obj)))]
    a2_idx = labels(obj)[findfirst(isequal("a2"), lowercase.(labels(obj)))]
    a1 = extract_channel(obj, ch=a1_idx)
    a2 = extract_channel(obj, ch=a2_idx)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    s_ref = similar(s)

    ref_label = repeat([""], length(chs))

    if type === :l
        @inbounds for ep_idx in 1:ep_n
            if !med
                ref_ch = @views vec(mean(vcat(a1[:, :, ep_idx], a2[:, :, ep_idx]), dims=1))
            else
                ref_ch = @views vec(median(vcat(a1[:, :, ep_idx], a2[:, :, ep_idx]), dims=1))
            end
            Threads.@threads for ch_idx in 1:ch_n
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-A1A2"
            end
        end
    elseif type === :i
        c_picks = channel_pick(obj, p=:central)
        @inbounds for ep_idx in 1:ep_n
            if !med
                ref_ch = @views vec(mean(vcat(a1[:, :, ep_idx], a2[:, :, ep_idx]), dims=1))
            else
                ref_ch = @views vec(median(vcat(a1[:, :, ep_idx], a2[:, :, ep_idx]), dims=1))
            end
            Threads.@threads for ch_idx in c_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-A1A2"
            end
        end
        l_picks = channel_pick(obj, p=:left)
        @inbounds for ep_idx in 1:ep_n
            ref_ch = @views vec(a1[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-A1"
            end
        end
        r_picks = channel_pick(obj, p=:right)
        @inbounds for ep_idx in 1:ep_n
            ref_ch = @views vec(a2[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-A2"
            end
        end
    elseif type === :c
        c_picks = channel_pick(obj, p=:central)
        @inbounds for ep_idx in 1:ep_n
            if !med
                ref_ch = @views vec(mean(vcat(a1[:, :, ep_idx], a2[:, :, ep_idx]), dims=1))
            else
                ref_ch = @views vec(median(vcat(a1[:, :, ep_idx], a2[:, :, ep_idx]), dims=1))
            end
            Threads.@threads for ch_idx in c_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-A1A2"
            end
        end
        l_picks = channel_pick(obj, p=:left)
        @inbounds for ep_idx in 1:ep_n
            ref_ch = @views vec(a2[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-A2"
            end
        end
        r_picks = channel_pick(obj, p=:right)
        @inbounds for ep_idx in 1:ep_n
            ref_ch = @views vec(a1[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-A1"
            end
        end
    end

    obj_new.data[chs, :, :] = s_ref
    obj_new.header.recording[:label][chs] .*= ref_label
    if type === :l
        obj_new.header.recording[:reference] = "auricular (linked)"
    elseif type === :i
        obj_new.header.recording[:reference] = "auricular (ipsilateral)"
    elseif type === :c
        obj_new.header.recording[:reference] = "auricular (contralateral)"
    end
    reset_components!(obj_new)
    push!(obj_new.history, "reference_a(OBJ, type=$type, med=$med)")

    return obj_new

end

"""
    reference_a!(obj; <keyword arguments>)

Reference to auricular (A1, A2) channels. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:l`:
    - `:l`: linked - average of A1 and A2
    - `:i`: ipsilateral - A1 for left channels, A2 for right channels
    - `:c`: contraletral - A1 for right channels, A2 for left channels
- `med::Bool=false`: use median instead of mean
"""
function reference_a!(obj::NeuroAnalyzer.NEURO; type::Symbol=:l, med::Bool=false)

    obj_new = reference_a(obj, type=type, med=med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    reference_m(obj; <keyword arguments>)

Reference to mastoid (M1, M2) channels. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:l`:
    - `:l`: linked - average of M1 and M2
    - `:i`: ipsilateral - M1 for left channels, M2 for right channels
    - `:c`: contraletral - M1 for right channels, M2 for left channels
- `med::Bool=false`: use median instead of mean

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function reference_m(obj::NeuroAnalyzer.NEURO; type::Symbol=:l, med::Bool=false)

    _check_datatype(obj, "eeg")
    _check_var(type, [:l, :i, :c], "type")
    @assert "m1" in lowercase.(labels(obj)) "OBJ does not contain M1 channel."
    @assert "m2" in lowercase.(labels(obj)) "OBJ does not contain M2 channel."

    # keep signal channels
    chs = [get_channel(obj, ch=idx)[1] for idx in get_channel(obj, type="eeg")]
    obj_new = deepcopy(obj)
    s = @view obj_new.data[chs, :, :]

    m1_idx = labels(obj)[findfirst(isequal("m1"), lowercase.(labels(obj)))]
    m2_idx = labels(obj)[findfirst(isequal("m2"), lowercase.(labels(obj)))]
    m1 = extract_channel(obj, ch=m1_idx)
    m2 = extract_channel(obj, ch=m2_idx)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    s_ref = similar(s)

    ref_label = repeat([""], length(chs))

    if type === :l
        @inbounds for ep_idx in 1:ep_n
            if !med
                ref_ch = @views vec(mean(vcat(m1[:, :, ep_idx], m2[:, :, ep_idx]), dims=1))
            else
                ref_ch = @views vec(median(vcat(m1[:, :, ep_idx], m2[:, :, ep_idx]), dims=1))
            end
            Threads.@threads for ch_idx in 1:ch_n
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-M1M2"
            end
        end
    elseif type === :i
        c_picks = channel_pick(obj, p=:central)
        @inbounds for ep_idx in 1:ep_n
            if !med
                ref_ch = @views vec(mean(vcat(m1[:, :, ep_idx], m2[:, :, ep_idx]), dims=1))
            else
                ref_ch = @views vec(median(vcat(m1[:, :, ep_idx], m2[:, :, ep_idx]), dims=1))
            end
            Threads.@threads for ch_idx in c_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-M1M2"
            end
        end
        l_picks = channel_pick(obj, p=:left)
        @inbounds for ep_idx in 1:ep_n
            ref_ch = @views vec(m1[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-M1"
            end
        end
        r_picks = channel_pick(obj, p=:right)
        @inbounds for ep_idx in 1:ep_n
            ref_ch = @views vec(m2[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-M2"
            end
        end
    elseif type === :c
        c_picks = channel_pick(obj, p=:central)
        @inbounds for ep_idx in 1:ep_n
            if !med
                ref_ch = @views vec(mean(vcat(m1[:, :, ep_idx], m2[:, :, ep_idx]), dims=1))
            else
                ref_ch = @views vec(median(vcat(m1[:, :, ep_idx], m2[:, :, ep_idx]), dims=1))
            end
            Threads.@threads for ch_idx in c_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-M1M2"
            end
        end
        l_picks = channel_pick(obj, p=:left)
        @inbounds for ep_idx in 1:ep_n
            ref_ch = @views vec(m2[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-M2"
            end
        end
        r_picks = channel_pick(obj, p=:right)
        @inbounds for ep_idx in 1:ep_n
            ref_ch = @views vec(m1[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
                ref_label[ch_idx] = "-M1"
            end
        end
    end

    obj_new.data[chs, :, :] = s_ref
    obj_new.header.recording[:label][chs] .*= ref_label
    if type === :l
        obj_new.header.recording[:reference] = "mastoid (linked)"
    elseif type === :i
        obj_new.header.recording[:reference] = "mastoid (ipsilateral)"
    elseif type === :c
        obj_new.header.recording[:reference] = "mastoid (contralateral)"
    end
    reset_components!(obj_new)
    push!(obj_new.history, "reference_m(OBJ, type=$type, med=$med)")

    return obj_new

end

"""
    reference_m!(obj; <keyword arguments>)

Reference to mastoid (M1, M2) channels. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:l`:
    - `:l`: linked - average of M1 and M2
    - `:i`: ipsilateral - M1 for left channels, M2 for right channels
    - `:c`: contraletral - M1 for right channels, M2 for left channels
- `med::Bool=false`: use median instead of mean
"""
function reference_m!(obj::NeuroAnalyzer.NEURO; type::Symbol=:l, med::Bool=false)

    obj_new = reference_m(obj, type=type, med=med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    reference_plap(obj; <keyword arguments>)

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `nn::Int64=4`: use `nn` adjacent electrodes
- `weighted::Bool=false`: use mean of `nn` nearest channels if false; if true, mean of `nn` nearest channels is weighted by distance to the referenced channel
- `med::Bool=false`: use median instead of mean

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function reference_plap(obj::NeuroAnalyzer.NEURO; nn::Int64=4, weighted::Bool=false, med::Bool=false)

    _check_datatype(obj, "eeg")
    _has_locs(obj)

    # keep signal channels
    chs = get_channel(obj, ch=get_channel(obj, type="eeg"))
    s = obj.data[chs, :, :]

    @assert length(chs) == nrow(obj.locs) "Some channels do not have locations."

    ch_n = size(s, 1)
    @assert nn >= 1 "nn must be ≥ 1"
    @assert nn < ch_n - 1 "nn must be < $(ch_n - 1)"
    ep_n = size(s, 3)

    loc_x = obj.locs[1:ch_n, :loc_x]
    loc_y = obj.locs[1:ch_n, :loc_y]

    # Euclidean distance matrix
    d = zeros(ch_n, ch_n)
    for idx1 in 1:ch_n
        for idx2 in 1:ch_n
            d[idx1, idx2] = euclidean([loc_x[idx1], loc_y[idx1]], [loc_x[idx2], loc_y[idx2]])
        end
    end
    # set weights not to reference to itself
    d[d .== 0] .= Inf

    # nn nearest neighbors index matrix
    nn_idx = zeros(Int64, ch_n, nn)
    for idx1 in 1:ch_n
        # nn_idx[idx1, :] = sortperm(d[idx1, :])[2:(nn + 1)] # 1st neighbor is the electrode itself
        nn_idx[idx1, :] = sortperm(d[idx1, :])[1:nn]
    end

    s_ref = zeros(size(s))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ref_chs = @view s[nn_idx[ch_idx, :], :, ep_idx]
            if !weighted
                if !med
                    ref_ch = vec(mean(ref_chs, dims=1))
                else
                    ref_ch = vec(median(ref_chs, dims=1))
                end
            else
                w = zeros(nn)
                for w_idx in 1:nn
                    w[w_idx] = euclidean([loc_x[ch_idx], loc_y[ch_idx]], [loc_x[nn_idx[ch_idx, w_idx]], loc_y[nn_idx[ch_idx, w_idx]]])
                end
                w = 1 .- w
                if !med
                    ref_ch = vec(mean(w .* ref_chs, dims=1))
                else
                    ref_ch = vec(median(w .* ref_chs, dims=1))
                end
            end
            s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
        end
    end

    obj_new = deepcopy(obj)
    obj_new.header.recording[:label][chs] .*= weighted ? "-wplap" : "-plap"
    obj_new.data[chs, :, :] = s_ref
    obj_new.header.recording[:reference] = weighted ? "weighted Laplacian ($nn)" : "Laplacian ($nn)"
    reset_components!(obj_new)
    push!(obj_new.history, "reference_plap(OBJ, nn=$nn, weighted=$weighted, med=$med)")

    return obj_new

end

"""
    reference_plap!(obj; <keyword arguments>)

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `nn::Int64=4`: use `nn` adjacent electrodes
- `weighted::Bool=false`: use distance weights; use mean of nearest channels if false
- `med::Bool=false`: use median instead of mean
"""
function reference_plap!(obj::NeuroAnalyzer.NEURO; nn::Int64=4, weighted::Bool=false, med::Bool=false)

    obj_new = reference_plap(obj, nn=nn, weighted=weighted, med=med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    reference_custom(obj; <keyword arguments>)

Reference using custom montage. Only signal channels are processed. Custom montage may be imported using `import_montage()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: list of channel pairs
- `ref_name::String="longitudinal-BIP"`: name of the montage

# Returns

- `obj_new::NeuroAnalyzer.NEURO`

# Notes

If the reference contains a single channel (e.g. "Fz"), than the channel is copied to the referenced data.
For each reference pair (e.g. "Fz-Cz"), the referenced channel is equal to the amplitude of channel 2 ("Cz") - amplitude of channel 1 ("Fz").

Examples of montages:
- bipolar transverse: ["Fp2-Fp1", "F8-Fp2", "F8-F4", "F4-Fz", "Fz-F3", "F3-F7", "Fp1-F7", "T4-C4", "C4-Cz", "Cz-C3", "C3-T3", "T6-P4", "P4-Pz", "Pz-P3", "P3-T5", "O2-O1"], "transverse-BIP"
- bipolar longitudinal: ["Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["Fz", "Cz", "Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
- bipolar longitudinal: ["FPz-Fz", "Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
"""
function reference_custom(obj::NeuroAnalyzer.NEURO; ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T4", "T4-T6", "T6-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], ref_name::String="longitudinal-BIP")

    _check_datatype(obj, "eeg")

    chs = get_channel(obj, ch=get_channel(obj, type="eeg"))

    for ref_idx in eachindex(ref_list)
        if '-' in ref_list[ref_idx]
            m = match(r"(.*)-(.+)", ref_list[ref_idx])
            @assert m[1] in labels(obj)[chs] "Label $(m[1]) does not match OBJ labels."
            @assert m[2] in labels(obj)[chs] "Label $(m[2]) does not match OBJ labels."
        else
            @assert ref_list[ref_idx] in labels(obj)[chs] "Label $(ref_list[ref_idx]) does not match OBJ labels."
        end
    end

    ep_n = nepochs(obj)
    s = zeros(length(ref_list), epoch_len(obj), ep_n)

    @inbounds for ep_idx in 1:ep_n
        for ref_idx in eachindex(ref_list)
            if '-' in ref_list[ref_idx]
                m = match(r"(.+)-(.+)", ref_list[ref_idx])
                ref1 = get_channel(obj, ch=string(m[1]))
                ref2 = get_channel(obj, ch=string(m[2]))
                s[ref_idx, :, ep_idx] = @views obj.data[ref1, :, ep_idx] - obj.data[ref2, :, ep_idx]
            else
                s[ref_idx, :, ep_idx] = @views obj.data[get_channel(obj, ch=ref_list[ref_idx]), :, ep_idx]
            end
        end
    end

    obj_new = delete_channel(obj, ch=get_channel(obj, type="eeg"))
    obj_new.data = vcat(s, obj_new.data)
    obj_new.header.recording[:label] = vcat(ref_list, labels(obj_new))
    obj_new.header.recording[:reference] = ref_name
    obj_new.header.recording[:channel_type] = vcat(repeat(["eeg"], length(ref_list)), obj_new.header.recording[:channel_type])
    obj_new.header.recording[:unit] = vcat(repeat(["μV"], length(ref_list)), obj_new.header.recording[:unit])
    obj_new.header.recording[:prefiltering] = vcat(repeat([obj.header.recording[:prefiltering][1]], length(ref_list)), obj_new.header.recording[:prefiltering])
    obj_new.header.recording[:transducers] = vcat(repeat([obj.header.recording[:transducers][1]], length(ref_list)), obj_new.header.recording[:transducers])
    obj_new.header.recording[:gain] = vcat(repeat([obj.header.recording[:gain][1]], length(ref_list)), obj_new.header.recording[:gain])

    reset_components!(obj_new)
    push!(obj_new.history, "reference_custom(OBJ, ref_list=$ref_list, ref_name=$ref_name)")

    return obj_new

end

"""
    reference_custom!(obj; <keyword arguments>)

Reference using custom montage. Only signal channels are processed. Custom montage may be imported using `import_montage()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: list of channel pairs
- `ref_name::String="BIP ||"`: name of the montage

# Notes

If the reference contains a single channel (e.g. "Fz"), than the channel is copied to the referenced data.
For each reference pair (e.g. "Fz-Cz"), the referenced channel is equal to the amplitude of channel 2 ("Cz") - amplitude of channel 1 ("Fz").

Examples of montages:
- bipolar transverse: ["Fp2-Fp1", "F8-Fp2", "F8-F4", "F4-Fz", "Fz-F3", "F3-F7", "Fp1-F7", "T4-C4", "C4-Cz", "Cz-C3", "C3-T3", "T6-P4", "P4-Pz", "Pz-P3", "P3-T5", "O2-O1"], "BIP ="
- bipolar longitudinal: ["Fz", "Cz", "Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], "BIP ||"
- bipolar longitudinal: ["Fp-Fz", "Fz-Cz", "Cz-Pz", "Pz-O", "Fp1-F7", "Fp1-F3", "F7-T7", "T7-P7", "P7-O1", "F3-C3", "C3-P3", "P3-O1", "Fp1-F7", "Fp2-F4", "F8-T8", "T8-P8", "P8-O2", "F4-C4", "C4-P4", "P4-O2"], "BIP ||"
"""
function reference_custom!(obj::NeuroAnalyzer.NEURO; ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], ref_name::String="BIP ||")

    obj_new = reference_custom(obj, ref_list=ref_list, ref_name=ref_name)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end