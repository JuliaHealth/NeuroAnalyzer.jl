export reference_ch
export reference_ch!
export reference_car
export reference_car!
export reference_a
export reference_a!
export reference_m
export reference_m!
export reference_plap
export reference_plap!

"""
    reference_ch(obj; ch, med)

Reference to selected channel(s). Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
- `med::Bool=false`: use median instead of mean

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function reference_ch(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}, med::Bool=false)

    # keep signal channels
    _check_channels(obj, ch)
    chs = signal_channels(obj)

    s = @view obj_new.data[chs, :, :]

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if length(ch) == 1
                ref_ch = @views vec(s[ch, :, ep_idx])
                if ch_idx != ch
                    @views s[ch_idx, :, ep_idx] .-= ref_ch
                end
            else
                if med == false
                    ref_ch = @views vec(mean(s[ch, :, ep_idx], dims=1))
                else
                    ref_ch = @views vec(median(s[ch, :, ep_idx], dims=1))
                end
                @views s[ch_idx, :, ep_idx] .-= ref_ch
            end
        end
    end

    obj_new = deepcopy(obj)
    obj_new.data[chs, :, :] = s
    obj_new.header.recording[:reference] = "channel: $ch"
    reset_components!(obj_new)
    push!(obj_new.header.history, "reference_ch(OBJ, ch=$ch, med=$med")

    return obj_new

end

"""
    reference_ch!(obj; ch, med)

Reference to selected channel(s). Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
- `med::Bool=false`: use median instead of mean
"""
function reference_ch!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}, med::Bool=false)

    obj_tmp = reference_ch(obj, ch=ch, med=med)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end

"""
    reference_car(obj; exclude_fpo, exclude_current, med)

Reference to common average reference. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `exclude_fpo::Bool=false`: exclude Fp1, Fp2, O1, O2 from CAR calculation
- `exclude_current::Bool=true`: exclude current channel from CAR calculation
- `med::Bool=false`: use median instead of mean

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function reference_car(obj::NeuroAnalyzer.NEURO; exclude_fpo::Bool=false, exclude_current::Bool=true, med::Bool=false)

    # keep signal channels
    chs = signal_channels(obj)
    s = @view obj_new.data[chs, :, :]

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            chs2exclude = Vector{Int64}()
            if exclude_fpo == true
                l = lowercase.(labels(obj_new))
                "fp1" in l && push!(chs2exclude, findfirst(isequal("fp1"), l))
                "fp2" in l && push!(chs2exclude, findfirst(isequal("fp2"), l))
                "o1" in l && push!(chs2exclude, findfirst(isequal("o1"), l))
                "o2" in l && push!(chs2exclude, findfirst(isequal("o2"), l))
            end
            exclude_current == true && push!(chs2exclude, ch_idx)
            ref_chs = @view s[setdiff(1:ch_n, unique(chs2exclude)), :, ep_idx]
            if med == false
                ref_ch = vec(mean(ref_chs, dims=1))
            else
                ref_ch = vec(median(ref_chs, dims=1))
            end
            @views s[ch_idx, :, ep_idx] .-= ref_ch
        end
    end

    obj_new = deepcopy(obj)
    obj_new.data[chs, :, :] = s
    obj_new.header.recording[:reference] = "CAR"
    reset_components!(obj_new)
    push!(obj_new.header.history, "reference_car(OBJ, exclude_fpo=$exclude_fpo, exclude_current=$exclude_current, med=$med))")

    return obj_new

end

"""
    reference_car!(obj; exclude_fpo, exclude_current, med)

Reference to common average reference. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `exclude_fpo::Bool=false`: exclude Fp1, Fp2, O1, O2 from CAR mean calculation
- `exclude_current::Bool=true`: exclude current channel from CAR mean calculation
- `med::Bool=false`: use median instead of mean
"""
function reference_car!(obj::NeuroAnalyzer.NEURO; exclude_fpo::Bool=false, exclude_current::Bool=true, med::Bool=false)

    obj_tmp = reference_car(obj, exclude_fpo=exclude_fpo, exclude_current=exclude_current, med=med)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end

"""
    reference_a(obj; type, med)

Reference to auricular (A1, A2) channels. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:l`:
    - `:l`: linked - average of A1 and A2
    - `:i`: ipsilateral - A1 for left channels, A2 for right channels
    - `:c`: contraletral - A1 for right channels, A2 for left channels
- `med::Bool=false`: use median instead of mean

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function reference_a(obj::NeuroAnalyzer.NEURO; type::Symbol=:l, med::Bool=false)

    _check_var(type, [:l, :i, :c], "type")
    all(iszero, occursin.("a1", lowercase.(obj.header.recording[:labels]))) == false || throw(ArgumentError("OBJ does not contain A1 channel."))
    all(iszero, occursin.("a2", lowercase.(obj.header.recording[:labels]))) == false || throw(ArgumentError("OBJ does not contain A2 channel."))

    # keep signal channels
    chs = signal_channels(obj)
    s = @view obj.data[chs, :, :]

    a1_idx = findfirst(isequal("A1"), obj.header.recording[:labels])
    a2_idx = findfirst(isequal("A2"), obj.header.recording[:labels])
    a1 = extract_channel(obj, ch=a1_idx)
    a2 = extract_channel(obj, ch=a2_idx)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    s_ref = similar(s)

    if type === :l
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(mean([a1[:, :, ep_idx], a2[:, :, ep_idx]]))
            Threads.@threads for ch_idx in 1:ch_n
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
    elseif type === :i
        c_picks = pick(obj, p=:central)
        @inbounds @simd for ep_idx in 1:ep_n
            if med == false
                ref_ch = @views vec(mean([a1[:, :, ep_idx], a2[:, :, ep_idx]]))
            else
                ref_ch = @views vec(median([a1[:, :, ep_idx], a2[:, :, ep_idx]]))
            end
            Threads.@threads for ch_idx in c_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        l_picks = pick(obj, p=:left)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(a1[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        r_picks = pick(obj, p=:right)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(a2[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
    elseif type === :c
        c_picks = pick(obj, p=:central)
        @inbounds @simd for ep_idx in 1:ep_n
            if med == false
                ref_ch = @views vec(mean([a1[:, :, ep_idx], a2[:, :, ep_idx]]))
            else
                ref_ch = @views vec(median([a1[:, :, ep_idx], a2[:, :, ep_idx]]))
            end
            Threads.@threads for ch_idx in c_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        l_picks = pick(obj, p=:left)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(a2[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        r_picks = pick(obj, p=:right)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(a1[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
    end

    obj_new = deepcopy(obj)
    obj_new.data[chs, :, :] = s_ref
    obj_new.header.recording[:reference] = "A ($type)"
    reset_components!(obj_new)
    push!(obj_new.header.history, "reference_a(OBJ, type=$type, med=$med)")

    return obj_new

end

"""
    reference_a!(obj; type, med)

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

    obj_tmp = reference_a(obj, type=type, med=med)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end

"""
    reference_m(obj; type, med)

Reference to mastoid (M1, M2) channels. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:l`:
    - `:l`: linked - average of M1 and M2
    - `:i`: ipsilateral - M1 for left channels, M2 for right channels
    - `:c`: contraletral - M1 for right channels, M2 for left channels
- `med::Bool=false`: use median instead of mean

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function reference_m(obj::NeuroAnalyzer.NEURO; type::Symbol=:l, med::Bool=false)

    _check_var(type, [:l, :i, :c], "type")
    all(iszero, occursin.("m1", lowercase.(obj.header.recording[:labels]))) == false || throw(ArgumentError("OBJ does not contain M1 channel."))
    all(iszero, occursin.("m2", lowercase.(obj.header.recording[:labels]))) == false || throw(ArgumentError("OBJ does not contain M2 channel."))

    # keep signal channels
    channels = signal_channels(obj)
    signal = @view obj.data[channels, :, :]

    m1_idx = findfirst(isequal("M1"), obj.header.recording[:labels])
    m2_idx = findfirst(isequal("M2"), obj.header.recording[:labels])
    m1 = extract_channel(obj, ch=m1_idx)
    m2 = extract_channel(obj, ch=m2_idx)

    ch_n = size(signal, 1)
    ep_n = size(signal, 3)
    s_ref = similar(signal)

    if type === :l
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(mean([m1[:, :, ep_idx], m2[:, :, ep_idx]]))
            Threads.@threads for ch_idx in 1:ch_n
                s_ref[ch_idx, :, ep_idx] = @views signal[ch_idx, :, ep_idx] .- ref_ch
            end
        end
    elseif type === :i
        c_picks = pick(obj, p=:central)
        @inbounds @simd for ep_idx in 1:ep_n
            if med == false
                ref_ch = @views vec(mean([m1[:, :, ep_idx], m2[:, :, ep_idx]]))
            else
                ref_ch = @views vec(median([m1[:, :, ep_idx], m2[:, :, ep_idx]]))
            end
            Threads.@threads for ch_idx in c_picks
                s_ref[ch_idx, :, ep_idx] = @views signal[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        l_picks = pick(obj, p=:left)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(m1[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views signal[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        r_picks = pick(obj, p=:right)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(m2[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views signal[ch_idx, :, ep_idx] .- ref_ch
            end
        end
    elseif type === :c
        c_picks = pick(obj, p=:central)
        @inbounds @simd for ep_idx in 1:ep_n
            if med == false
                ref_ch = @views vec(mean([m1[:, :, ep_idx], m2[:, :, ep_idx]]))
            else
                ref_ch = @views vec(median([m1[:, :, ep_idx], m2[:, :, ep_idx]]))
            end
            Threads.@threads for ch_idx in c_picks
                s_ref[ch_idx, :, ep_idx] = @views signal[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        l_picks = pick(obj, p=:left)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(m2[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views signal[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        r_picks = pick(obj, p=:right)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(m1[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views signal[ch_idx, :, ep_idx] .- ref_ch
            end
        end
    end

    obj_new = deepcopy(obj)
    obj_new.data[channels, :, :] = s_ref
    obj_new.header.recording[:reference] = "M ($type)"
    reset_components!(obj_new)
    push!(obj_new.header.history, "reference_m(OBJ, type=$type, med=$med)")

    return obj_new
end

"""
    reference_m!(obj; type, med)

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

    obj_tmp = reference_m(obj, type=type, med=med)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end

"""
    reference_plap(obj; nn, weights)

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `nn::Int64=4`: use `nn` adjacent electrodes
- `weights::Bool=false`: use mean of `nn` nearest channels if false; if true, mean of `nn` nearest channels is weighted by distance to the referenced channel
- `med::Bool=false`: use median instead of mean

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function reference_plap(obj::NeuroAnalyzer.NEURO; nn::Int64=4, weights::Bool=false, med::Bool=false)

    obj.header.has_locs == false && throw(ArgumentError("Electrode locations not available, use load_locs() or add_locs() first."))

    # keep signal channels
    obj_tmp = deepcopy(obj)
    channels = signal_channels(obj)
    signal = obj_tmp.data[channels, :, :]

    length(channels) > nrow(obj.locs) && throw(ArgumentError("Some channels do not have locations."))

    ch_n = size(signal, 1)
    nn < 1 && throw(ArgumentError("nn must be ≥ 1"))
    nn > ch_n - 1 && throw(ArgumentError("nn must be < $(ch_n - 1)"))
    ep_n = size(signal, 3)
    
    loc_x = zeros(ch_n)
    loc_y = zeros(ch_n)
    for idx in 1:ch_n
        loc_x[idx] = obj.locs[idx, :loc_x]
        loc_y[idx] = obj.locs[idx, :loc_y]
    end

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

    s_ref = zeros(size(signal))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            reference_channels = @view signal[nn_idx[ch_idx, :], :, ep_idx]
            if weights == false
                if med == false
                    reference_channel = vec(mean(reference_channels, dims=1))
                else
                    reference_channel = vec(median(reference_channels, dims=1))
                end
            else
                g = Vector{Float64}()
                for idx1 in 1:nn
                    push!(g, 1 / d[ch_idx, nn_idx[ch_idx, idx1]] / sum(1 / d[ch_idx, nn_idx[ch_idx, :]]))
                end
                if med == false
                    reference_channel = vec(mean(g .* reference_channels, dims=1))
                else
                    reference_channel = vec(median(g .* reference_channels, dims=1))
                end
            end
            s_ref[ch_idx, :, ep_idx] = @views signal[ch_idx, :, ep_idx] .- reference_channel
        end
    end

    obj_tmp.data[channels, :, :] = s_ref
    obj_tmp.header.recording[:reference] = "PLAP ($nn)"
    reset_components!(obj_tmp)
    push!(obj_tmp.header.history, "reference_plap(OBJ, nn=$nn, med=$med))")

    return obj_tmp
end

"""
    reference_plap!(obj; nn, weights)

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `nn::Int64=4`: use `nn` adjacent electrodes
- `weights::Bool=false`: use distance weights; use mean of nearest channels if false
- `med::Bool=false`: use median instead of mean
"""
function reference_plap!(obj::NeuroAnalyzer.NEURO; nn::Int64=4, weights::Bool=false, med::Bool=false)

    obj_tmp = reference_plap(obj, nn=nn, weights=weights, med=med)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end
