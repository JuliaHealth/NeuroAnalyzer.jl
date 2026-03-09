export dpli

"""
    dpli(s1, s2)

Calculate Directed Phase Lag Index (dPLI). dPLI quantifies the fraction of time that s1 leads s2 in phase: `dPLI = ( {t : 0 < phd[t] < π} + 0.5 · {t : phd[t] = 0} ) / N`

where phd = s1_phase − s2_phase ∈ (−π, π].

- dPLI > 0.5 → s1 consistently leads s2
- dPLI < 0.5 → s2 consistently leads s1
- dPLI ≈ 0.5 → no consistent phase relationship

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

Named tuple containing:

- `pv::Float64`: dPLI value ∈ [0, 1]
- `sd::Vector{Float64}`: signal difference (s2 - s1)
- `phd::Vector{Float64}`: instantaneous phase difference (s1_phase - s2_phase) ∈ (−π, π]
- `s1ph::Vector{Float64}`: instantaneous phase of s1
- `s2ph::Vector{Float64}`: instantaneous phase of s2

# Reference

 1. Stam, C. J., & van Straaten, E. C. W. (2012). Go with the flow: Use of a directed phase lag index (dPLI) to characterize patterns of phase relations in a large-scale model of brain dynamics. NeuroImage, 62(3), 1415–1428.
"""
function dpli(
    s1::AbstractVector,
    s2::AbstractVector
)::@NamedTuple{
    pv::Float64,
    sd::Vector{Float64},
    phd::Vector{Float64},
    s1ph::Vector{Float64},
    s2ph::Vector{Float64}
}

    @assert length(s1) == length(s2) "Both signals must have the same length."

    # extract instantaneous phase via Hilbert transform
    h1 = htransform(s1)
    h2 = htransform(s2)
    s1ph = h1.ph
    s2ph = h2.ph

    # signal difference and instantaneous phase difference.
    sd = s1 - s2

    # instantaneous phase difference ∈ (−π, π] by convention of htransform output
    phd = s1ph - s2ph

    # dPLI calculation
    # r1: 0 < phd < π → s1 leads s2
    # r2: −π < phd < 0 → s2 leads s1
    # r3: phd == 0 → instantaneously synchronous (no direction)
    r1 = count(x ->  0 <  x <  π, phd)
    r2 = count(x -> -π <  x <  0, phd)   # kept for completeness; cancels in formula
    r3 = count(x ->  x == 0,      phd)

    # dPLI: fraction of time s1 leads s2, with ties counting as 0.5
    # r2 cancels out of the numerator — only r1 and the tie-break matter
    pv = (r1 + 0.5 * r3) / length(s1)

    return (pv = pv, sd = sd, phd = phd, s1ph = s1ph, s2ph = s2ph)

end

"""
    dpli(obj1, obj2; <keyword arguments>)

Calculate Directed Phase Lag Index (dPLI). dPLI quantifies the fraction of time that s1 leads s2 in phase: `dPLI = ( {t : 0 < phd[t] < π} + 0.5 · {t : phd[t] = 0} ) / N`

where phd = s1_phase − s2_phase ∈ (−π, π].

- dPLI > 0.5 → s1 consistently leads s2
- dPLI < 0.5 → s2 consistently leads s1
- dPLI ≈ 0.5 → no consistent phase relationship

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

Named tuple containing:

- `pv::Matrix{Float64}`: dPLI values of shape `(channels, epochs)`
- `sd::Array{Float64, 3}`: signal difference of shape `(channels, samples, epochs)`
- `phd::Array{Float64, 3}`: phase difference of shape `(channels, samples, epochs)`
- `s1ph::Array{Float64, 3}`: signal 1 instantaneous phase of shape `(channels, samples, epochs)`
- `s2ph::Array{Float64, 3}`: signal 2 instantaneous phase of shape `(channels, samples, epochs)`
"""
function dpli(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
)::@NamedTuple{
    pv::Matrix{Float64},
    sd::Array{Float64, 3},
    phd::Array{Float64, 3},
    s1ph::Array{Float64, 3},
    s2ph::Array{Float64, 3},
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    # number of channels
    ch_n = length(ch1)
    # number of epochs
    ep_n = length(ep1)
    # epoch length
    ep_len  = epoch_len(obj1)

    # pre-allocate outputs
    pv   = zeros(ch_n, ep_n)
    sd   = zeros(ch_n, ep_len, ep_n)
    phd  = zeros(ch_n, ep_len, ep_n)
    s1ph = zeros(ch_n, ep_len, ep_n)
    s2ph = zeros(ch_n, ep_len, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        dpli_data = dpli(
            @view(obj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(obj2.data[ch2[ch_idx], :, ep2[ep_idx]]),
        )
        pv[ch_idx, ep_idx] = dpli_data.pv
        sd[ch_idx, :, ep_idx] = dpli_data.sd
        phd[ch_idx, :, ep_idx] = dpli_data.phd
        s1ph[ch_idx, :, ep_idx] = dpli_data.s1ph
        s2ph[ch_idx, :, ep_idx] = dpli_data.s2ph
    end

    return (pv = pv, sd = sd, phd = phd, s1ph = s1ph, s2ph = s2ph)

end

"""
    dpli(obj; <keyword arguments>)

Calculate Directed Phase Lag Index (dPLI). dPLI quantifies the fraction of time that s1 leads s2 in phase: `dPLI = ( {t : 0 < phd[t] < π} + 0.5 · {t : phd[t] = 0} ) / N`

where phd = s1_phase − s2_phase ∈ (−π, π].

- dPLI > 0.5 → s1 consistently leads s2
- dPLI < 0.5 → s2 consistently leads s1
- dPLI ≈ 0.5 → no consistent phase relationship

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `pv::Array{Float64, 3}`: dPLI values of shape `(channels, channels, epochs)`
"""
function dpli(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = length(ep)

    # pre-allocate outpu
    pv = zeros(ch_n, ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 in 1:ch_idx1
            pv[ch_idx1, ch_idx2, ep_idx] = dpli(
                @view(obj.data[ch[ch_idx1], :, ep_idx]),
                @view(obj.data[ch[ch_idx2], :, ep_idx]),
            ).pv
        end
    end

    # mirror lower triangle to upper triangle
    pv = _copy_lt2ut(pv)

    return pv

end
