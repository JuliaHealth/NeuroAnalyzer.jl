export dpli

"""
    dpli(s1, s2)

Calculate Directed Phase Lag Index (dPLI).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

Named tuple containing:
- `pv::Float64`: dPLI value
- `sd::Vector{Float64}`: signal difference (s2 - s1)
- `phd::Vector{Float64}`: phase difference (s2 - s1)
- `s1ph::Vector{Float64}`: signal 1 phase
- `s2ph::Vector{Float64}`: signal 2 phase

# Source

1. Stam, C. J., & van Straaten, E. C. W. (2012). Go with the flow: Use of a directed phase lag index (dPLI) to characterize patterns of phase relations in a large-scale model of brain dynamics. NeuroImage, 62(3), 1415–1428.
"""
function dpli(s1::AbstractVector, s2::AbstractVector)::@NamedTuple{pv::Float64, sd::Vector{Float64}, phd::Vector{Float64}, s1ph::Vector{Float64}, s2ph::Vector{Float64}}

    @assert length(s1) == length(s2) "Both signals must have the same length."

    # get instatenous phases
    _, _, _, s1ph = htransform(s1)
    _, _, _, s2ph = htransform(s2)

    # signal difference
    sd = s1 - s2

    # phase differences
    phd = s1ph - s2ph

    rel_phase = mod.(phd, 2.0 * pi)

    r1 = intersect(findall(x -> 0 .<= x, rel_phase), findall(x -> x .< pi, rel_phase))
    r2 = intersect(findall(x -> x .>= -pi, rel_phase), findall(x -> x .< 0, rel_phase))
    r3 = findall(x -> x .== 0, rel_phase)

    pv = (length(r1) + length(r2) - 0.5 * length(r3)) / length(s1)

    return (pv=pv, sd=sd, phd=phd, s1ph=s1ph, s2ph=s2ph)

end

"""
    dpli(obj1, obj2; <keyword arguments>)

Calculate Directed Phase Lag Index (dPLI).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `pv::Matrix{Float64}`: dPLI value
- `sd::Array{Float64, 3}`: signal difference (s2 - s1)
- `phd::Array{Float64, 3}`: phase difference (s2 - s1)
- `s1ph::Array{Float64, 3}`: signal 1 phase
- `s2ph::Array{Float64, 3}`: signal 2 phase
"""
function dpli(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))::@NamedTuple{pv::Matrix{Float64}, sd::Array{Float64, 3}, phd::Array{Float64, 3}, s1ph::Array{Float64, 3}, s2ph::Array{Float64, 3}}

    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1)) and ep2 ($(length(ep2)) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    ep_n = length(ep1)
    ch_n = length(ch1)

    pv = zeros(ch_n, ep_n)
    sd = zeros(ch_n, epoch_len(obj1), ep_n)
    phd = zeros(ch_n, epoch_len(obj1), ep_n)
    s1ph = zeros(ch_n, epoch_len(obj1), ep_n)
    s2ph = zeros(ch_n, epoch_len(obj1), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pv[ch_idx, ep_idx], sd[ch_idx, :, ep_idx], phd[ch_idx, :, ep_idx], s1ph[ch_idx, :, ep_idx], s2ph[ch_idx, :, ep_idx] = @views dpli(obj1.data[ch1[ch_idx], :, ep1[ep_idx]], obj2.data[ch2[ch_idx], :, ep2[ep_idx]])
        end
    end

    return (pv=pv, sd=sd, phd=phd, s1ph=s1ph, s2ph=s2ph)

end

"""
    dpli(obj; <keyword arguments>)

Calculate Directed Phase Lag Index (dPLI).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `pv::Array{Float64, 3}`: dPLI value
"""
function dpli(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    isa(ch, Int64) && (ch = [ch])

    pv = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                pv[ch_idx1, ch_idx2, ep_idx], _, _, _, _ = @views dpli(obj.data[ch[ch_idx1], :, ep_idx], obj.data[ch[ch_idx2], :, ep_idx])
            end
        end
    end

    pv = _copy_lt2ut(pv)

    return pv

end
