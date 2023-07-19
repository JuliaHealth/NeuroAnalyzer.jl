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
export reference_custom
export reference_custom!

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

    _check_datatype(obj, :eeg)

    # keep signal channels
    _check_channels(obj, ch)
    chs = signal_channels(obj)

    s = @view obj.data[chs, :, :]

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
    push!(obj_new.history, "reference_ch(OBJ, ch=$ch, med=$med")

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

    obj_new = reference_ch(obj, ch=ch, med=med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

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

    _check_datatype(obj, :eeg)

    # keep signal channels
    chs = signal_channels(obj)
    s = @view obj.data[chs, :, :]

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
    push!(obj_new.history, "reference_car(OBJ, exclude_fpo=$exclude_fpo, exclude_current=$exclude_current, med=$med))")

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

    obj_new = reference_car(obj, exclude_fpo=exclude_fpo, exclude_current=exclude_current, med=med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

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

    _check_datatype(obj, :eeg)
    _check_var(type, [:l, :i, :c], "type")
    all(iszero, occursin.("a1", lowercase.(labels(obj)))) == false || throw(ArgumentError("OBJ does not contain A1 channel."))
    all(iszero, occursin.("a2", lowercase.(labels(obj)))) == false || throw(ArgumentError("OBJ does not contain A2 channel."))

    # keep signal channels
    chs = signal_channels(obj)
    s = @view obj.data[chs, :, :]

    a1_idx = findfirst(isequal("a1"), lowercase.(labels(obj)))
    a2_idx = findfirst(isequal("a2"), lowercase.(labels(obj)))
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
    push!(obj_new.history, "reference_a(OBJ, type=$type, med=$med)")

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

    obj_new = reference_a(obj, type=type, med=med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

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

    _check_datatype(obj, :eeg)
    _check_var(type, [:l, :i, :c], "type")
    all(iszero, occursin.("m1", lowercase.(labels(obj)))) == false || throw(ArgumentError("OBJ does not contain M1 channel."))
    all(iszero, occursin.("m2", lowercase.(labels(obj)))) == false || throw(ArgumentError("OBJ does not contain M2 channel."))

    # keep signal channels
    chs = signal_channels(obj)
    s = @view obj.data[chs, :, :]

    m1_idx = findfirst(isequal("m1"), lowercase.(labels(obj)))
    m2_idx = findfirst(isequal("m2"), lowercase.(labels(obj)))
    m1 = extract_channel(obj, ch=m1_idx)
    m2 = extract_channel(obj, ch=m2_idx)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    s_ref = similar(s)

    if type === :l
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(mean([m1[:, :, ep_idx], m2[:, :, ep_idx]]))
            Threads.@threads for ch_idx in 1:ch_n
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
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
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        l_picks = pick(obj, p=:left)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(m1[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        r_picks = pick(obj, p=:right)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(m2[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
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
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        l_picks = pick(obj, p=:left)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(m2[:, :, ep_idx])
            Threads.@threads for ch_idx in l_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
        r_picks = pick(obj, p=:right)
        @inbounds @simd for ep_idx in 1:ep_n
            ref_ch = @views vec(m1[:, :, ep_idx])
            Threads.@threads for ch_idx in r_picks
                s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
            end
        end
    end

    obj_new = deepcopy(obj)
    obj_new.data[chs, :, :] = s_ref
    obj_new.header.recording[:reference] = "M ($type)"
    reset_components!(obj_new)
    push!(obj_new.history, "reference_m(OBJ, type=$type, med=$med)")

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

    obj_new = reference_m(obj, type=type, med=med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

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

    _check_datatype(obj, :eeg)
    _has_locs(obj) == false && throw(ArgumentError("Electrode locations not available, use load_locs() or add_locs() first."))

    # keep signal channels
    chs = signal_channels(obj)
    s = obj.data[chs, :, :]

    length(chs) > nrow(obj.locs) && throw(ArgumentError("Some channels do not have locations."))

    ch_n = size(s, 1)
    nn < 1 && throw(ArgumentError("nn must be ≥ 1"))
    nn > ch_n - 1 && throw(ArgumentError("nn must be < $(ch_n - 1)"))
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

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ref_chs = @view s[nn_idx[ch_idx, :], :, ep_idx]
            if weights == false
                if med == false
                    ref_ch = vec(mean(ref_chs, dims=1))
                else
                    ref_ch = vec(median(ref_chs, dims=1))
                end
            else
                g = Vector{Float64}()
                for idx1 in 1:nn
                    push!(g, 1 / d[ch_idx, nn_idx[ch_idx, idx1]] / sum(1 / d[ch_idx, nn_idx[ch_idx, :]]))
                end
                if med == false
                    ref_ch = vec(mean(g .* ref_chs, dims=1))
                else
                    ref_ch = vec(median(g .* ref_chs, dims=1))
                end
            end
            s_ref[ch_idx, :, ep_idx] = @views s[ch_idx, :, ep_idx] .- ref_ch
        end
    end

    obj_new = deepcopy(obj)
    obj_new.data[chs, :, :] = s_ref
    obj_new.header.recording[:reference] = "PLAP ($nn)"
    reset_components!(obj_new)
    push!(obj_new.history, "reference_plap(OBJ, nn=$nn, med=$med))")

    return obj_new

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

    obj_new = reference_plap(obj, nn=nn, weights=weights, med=med)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    reference_custom(obj; ref_list, ref_name)

Reference using custom montage. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: list of channel pairs
- `ref_name::String="BIP ||"`: name of the montage

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

If the reference contains a single channel (e.g. "Fz"), than the channel is copied to the referenced data.
For each reference pair (e.g. "Fz-Cz"), the referenced channel is equal to the amplitude of channel 1 ("Fz") - amplitude of channel 2 ("Cz").

Examples of montages:
- bipolar transverse: ["Fp2-Fp1", "F8-Fp2", "F8-F4", "F4-Fz", "Fz-F3", "F3-F7", "Fp1-F7", "T4-C4", "C4-Cz", "Cz-C3", "C3-T3", "T6-P4", "P4-Pz", "Pz-P3", "P3-T5", "O2-O1"], "BIP ="
- bipolar longitudinal: ["Fz", "Cz", "Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], "BIP ||"
- bipolar longitudinal: ["Fp-Fz", "Fz-Cz", "Cz-Pz", "Pz-O", "Fp1-F7", "Fp1-F3", "F7-T7", "T7-P7", "P7-O1", "F3-C3", "C3-P3", "P3-O1", "Fp1-F7", "Fp2-F4", "F8-T8", "T8-P8", "P8-O2", "F4-C4", "C4-P4", "P4-O2"], "BIP ||"
"""
function reference_custom(obj::NeuroAnalyzer.NEURO; ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], ref_name::String="BIP ||")

    _check_datatype(obj, :eeg)

    chs = signal_channels(obj)

    for ref_idx in 1:length(ref_list)
        if '-' in ref_list[ref_idx]
            m = match(r"(.+)-(.+)", ref_list[ref_idx])
            m[1] in labels(obj)[chs] || throw(ArgumentError("Label $(m[1]) does not match OBJ labels."))
            m[2] in labels(obj)[chs] || throw(ArgumentError("Label $(m[2]) does not match OBJ labels."))
        else
            ref_list[ref_idx] in labels(obj)[chs] || throw(ArgumentError("Label $(ref_list[ref_idx]) does not match OBJ labels."))
        end
    end

    ep_n = epoch_n(obj)
    s = zeros(length(ref_list), epoch_len(obj), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        for ref_idx in 1:length(ref_list)
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

    obj_new = delete_channel(obj, ch=signal_channels(obj))
    obj_new.data = vcat(s, obj_new.data)
    obj_new.header.recording[:labels] = vcat(ref_list, labels(obj_new))
    obj_new.header.recording[:reference] = ref_name
    obj_new.header.recording[:channel_type] = vcat(repeat(["eeg"], length(ref_list)), obj_new.header.recording[:channel_type])
    obj_new.header.recording[:units] = vcat(repeat(["μV"], length(ref_list)), obj_new.header.recording[:units])
    obj_new.header.recording[:prefiltering] = vcat(repeat([obj.header.recording[:prefiltering][1]], length(ref_list)), obj_new.header.recording[:prefiltering])
    obj_new.header.recording[:transducers] = vcat(repeat([obj.header.recording[:transducers][1]], length(ref_list)), obj_new.header.recording[:transducers])
    obj_new.header.recording[:gain] = vcat(repeat([obj.header.recording[:gain][1]], length(ref_list)), obj_new.header.recording[:gain])

    reset_components!(obj_new)
    push!(obj_new.history, "reference_custom(OBJ, ref_list=$ref_list, ref_name=$ref_name)")

    return obj_new

end

"""
    reference_custom!(obj; ref_list, ref_name)

Reference using custom montage. Only signal channels are processed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: list of channel pairs
- `ref_name::String="BIP ||"`: name of the montage

# Notes

If the reference contains a single channel (e.g. "Fz"), than the channel is copied to the referenced data.
For each reference pair (e.g. "Fz-Cz"), the referenced channel is equal to the amplitude of channel 1 ("Fz") - amplitude of channel 2 ("Cz").

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