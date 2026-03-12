export ispc

"""
    ispc(s1, s2)

Calculate ISPC (Inter-Site-Phase Clustering) between two signals. ISPC measures the consistency of the phase difference between two signals across trials (or, here, across time within a single trial):

ISPC value = |mean( exp(i·Δφ) )| (0 = no clustering, 1 = perfect lock)

ISPC angle = angle( mean( exp(i·Δφ) ) ) (preferred phase difference)

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

Named tuple:

- `ispcv::Float64`: ISPC value (phase-locking magnitude)
- `ispca::Float64`: ISPC angle (preferred phase difference)
- `sd::Vector{Float64}`: signal difference (s1 - s2)
- `phd::Vector{Float64}`: phase difference (s1 - s2)
- `s1ph::Vector{Float64}`: instantaneous phase of s1
- `s2ph::Vector{Float64}`: instantaneous phase of s2
"""
function ispc(
    s1::AbstractVector,
    s2::AbstractVector
)::@NamedTuple{
    ispcv::Float64,
    ispca::Float64,
    sd::Vector{Float64},
    phd::Vector{Float64},
    s1ph::Vector{Float64},
    s2ph::Vector{Float64}
}

    @assert length(s1) == length(s2) "Both signals must have the same length."

    h1 = htransform(s1)
    h2 = htransform(s2)
    s1ph = h1.ph
    s2ph = h2.ph

    # signal and phase differences
    sd  = s1 .- s2
    phd = s1ph .- s2ph

    # compute the mean complex phase difference once and reuse it for both the magnitude (ispcv) and angle (ispca)
    mcphd = mean(Base.cis.(phd)) # cis.(phd) = exp.(im .* phd)
    ispcv = abs(mcphd)
    ispca = DSP.angle(mcphd)

    return (
        ispcv = ispcv,
        ispca = ispca,
        sd = sd,
        phd = phd,
        s1ph = s1ph,
        s2ph = s2ph,
    )

end

"""
    ispc(obj; <keyword arguments>)

Calculate ISPC (Inter-Site Phase Clustering) for all channel pairs. ISPC measures the consistency of the phase difference between two signals across trials (or, here, across time within a single trial):

ISPC value = |mean( exp(i·Δφ) )| (0 = no clustering, 1 = perfect lock)

ISPC angle = angle( mean( exp(i·Δφ) ) ) (preferred phase difference)

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

Named tuple:

- `ispcv::Array{Float64, 3}`: ISPC value matrices, shape `(channels, channels, epochs)`
- `ispca::Array{Float64, 3}`: ISPC angle matrices, shape `(channels, channels, epochs)`
"""
function ispc(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex}
)::@NamedTuple{
    ispcv::Array{Float64, 3},
    ispca::Array{Float64, 3}
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    # pre-allocate symmetric output: channels × channels × epochs
    # only the lower triangle is computed
    ispcv = zeros(ch_n, ch_n, ep_n)
    ispca = zeros(ch_n, ch_n, ep_n)

    @inbounds Threads.@threads :dynamic for ep_idx in 1:ep_n
        for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1 - 1
                ispc_data = ispc(
                    @view(obj.data[ch[ch_idx1], :, ep_idx]),
                    @view(obj.data[ch[ch_idx2], :, ep_idx])
                )
                ispcv[ch_idx1, ch_idx2, ep_idx] = ispc_data.ispcv
                ispca[ch_idx1, ch_idx2, ep_idx] = ispc_data.ispca
            end
        end
    end

    # copy lower triangle to upper triangle
    ispcv = _copy_lt2ut(ispcv)
    ispca = _copy_lt2ut(ispca)

    return (ispcv = ispcv, ispca = ispca)

end

"""
    ispc(obj1, obj2; <keyword arguments>)

Calculate ISPC (Inter-Site Phase Clustering) between channel-matched pairs.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

Named tuple:

- `ispcv::Matrix{Float64}`: ISPC value, shape `(channels, epochs)`
- `ispca::Matrix{Float64}`: ISPC angle, shape `(channels, epochs)`
- `sd::Array{Float64, 3}`: signal difference (s1 - s2), shape `(channels, samples, epochs)`
- `phd::Array{Float64, 3}`: phase difference (s1 - s2), shape `(channels, samples, epochs)`
- `s1ph::Array{Float64, 3}`: signal 1 phases, shape `(channels, samples, epochs)`
- `s2ph::Array{Float64, 3}`: signal 2 phases, shape `(channels, samples, epochs)`
"""
function ispc(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2))
)::@NamedTuple{
    ispcv::Matrix{Float64},
    ispca::Matrix{Float64},
    sd::Array{Float64, 3},
    phd::Array{Float64, 3},
    s1ph::Array{Float64, 3},
    s2ph::Array{Float64, 3}
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

    # pre-allocate output
    ispcv = zeros(ch_n, ep_n)
    ispca = zeros(ch_n, ep_n)
    sd = zeros(ch_n, epoch_len(obj1), ep_n)
    phd = zeros(ch_n, epoch_len(obj1), ep_n)
    s1ph = zeros(ch_n, epoch_len(obj1), ep_n)
    s2ph = zeros(ch_n, epoch_len(obj1), ep_n)

    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ispc_data = ispc(
            @view(obj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(obj2.data[ch2[ch_idx], :, ep2[ep_idx]])
        )
        ispcv[ch_idx, ep_idx] = ispc_data.ispcv
        ispca[ch_idx, ep_idx] = ispc_data.ispca
        sd[ch_idx, :, ep_idx] = ispc_data.sd
        phd[ch_idx, :, ep_idx] = ispc_data.phd
        s1ph[ch_idx, :, ep_idx] = ispc_data.s1ph
        s2ph[ch_idx, :, ep_idx] = ispc_data.s2ph
    end

    return (
        ispcv = ispcv,
        ispca = ispca,
        sd = sd,
        phd = phd,
        s1ph = s1ph,
        s2ph = s2ph
    )

end
