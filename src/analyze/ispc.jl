export ispc

"""
    ispc(s1, s2)

Calculate ISPC (Inter-Site-Phase Clustering) between `s1` and `s2`.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

Named tuple containing:
- `ispc_val::Float64`: ISPC value
- `ispc_ang::Float64`: ISPC angle
- `s_diff::Vector{Float64}`: signal difference (s2 - s1)
- `ph_diff::Vector{Float64}`: phase difference (s2 - s1)
- `s1_phase::Vector{Float64}`: signal 1 phase
- `s2_phase::Vector{Float64}`: signal 2 phase
"""
function ispc(s1::AbstractVector, s2::AbstractVector)::NamedTuple{ispc_val::Float64, ispc_ang::Float64, s_diff::Vector{Float64}, ph_diff::Vector{Float64}, s1_phase::Vector{Float64}, s2_phase::Vector{Float64}}

    @assert length(s1) == length(s2) "Both signals must have the same length."

    _, _, _, s1_phase = hspectrum(s1)
    _, _, _, s2_phase = hspectrum(s2)

    s_diff = s2 - s1
    ph_diff = s2_phase - s1_phase

    ispc_val = abs(mean(exp.(1im .* ph_diff)))
    ispc_ang = angle(mean(exp.(1im .* ph_diff)))

    return (ispc_val=ispc_val, ispc_ang=ispc_ang, s_diff=s_diff, ph_diff=ph_diff, s1_phase=s1_phase, s2_phase=s2_phase)

end

"""
    ispc(obj; <keyword arguments>)

Calculate ISPCs (Inter-Site-Phase Clustering).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

Named tuple containing:
- `ispc_val::Array{Float64, 3}`: ISPC value matrices over epochs
- `ispc_ang::Array{Float64, 3}`: ISPC angle matrices over epochs
"""
function ispc(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})::NamedTuple{ispc_val::Array{Float64, 3}, ispc_ang::Array{Float64, 3}}

    ch = get_channel(obj, ch=ch)
    ch_n = length(ch)
    ep_n = nepochs(obj)

    ispc_val = zeros(ch_n, ch_n, ep_n)
    ispc_ang = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                ispc_val[ch_idx1, ch_idx2, ep_idx], ispc_ang[ch_idx1, ch_idx2, ep_idx], _, _, _, _ = @views ispc(obj.data[ch[ch_idx1], :, ep_idx], obj.data[ch[ch_idx2], :, ep_idx])
            end
        end
    end

    # copy lower triangle to upper triangle
    ispc_val = _copy_lt2ut(ispc_val)
    ispc_ang = _copy_lt2ut(ispc_ang)

    return (ispc_val=ispc_val, ispc_ang=ispc_ang)

end

"""
    ispc(obj1, obj2; <keyword arguments>)

Calculate ISPC (Inter-Site-Phase Clustering).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `ispc_val::Matrix{Float64}`: ISPC value
- `ispc_ang::Matrix{Float64}`: ISPC angle
- `s_diff::Array{Float64, 3}`: signal difference (s2 - s1)
- `ph_diff::Array{Float64, 3}`: phase difference (s2 - s1)
- `s1_phase::Array{Float64, 3}`: signal 1 phase
- `s2_phase::Array{Float64, 3}`: signal 2 phase
"""
function ispc(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2)))::NamedTuple{ispc_val::Matrix{Float64}, ispc_ang::Matrix{Float64}, s_diff::Array{Float64, 3}, ph_diff::Array{Float64, 3}, s1_phase::Array{Float64, 3}, s2_phase::Array{Float64, 3}}

    ch1 = get_channel(obj1, ch=ch1)
    ch2 = get_channel(obj2, ch=ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    ep_n = length(ep1)
    ch_n = length(ch1)

    ispc_val = zeros(ch_n, ep_n)
    ispc_ang = zeros(ch_n, ep_n)
    s_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    ph_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    s1_phase = zeros(ch_n, epoch_len(obj1), ep_n)
    s2_phase = zeros(ch_n, epoch_len(obj1), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ispc_val[ch_idx, ep_idx], ispc_ang[ch_idx, ep_idx], s_diff[ch_idx, :, ep_idx], ph_diff[ch_idx, :, ep_idx], s1_phase[ch_idx, :, ep_idx], s2_phase[ch_idx, :, ep_idx] = @views ispc(obj1.data[ch1[ch_idx], :, ep1[ep_idx]], obj2.data[ch2[ch_idx], :, ep2[ep_idx]])
        end
    end

    return (ispc_val=ispc_val, ispc_ang=ispc_ang, s_diff=s_diff, ph_diff=ph_diff, s1_phase=s1_phase, s2_phase=s2_phase)

end
