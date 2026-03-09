export iplv

"""
    iplv(s1, s2)

Calculate Imaginary Phase Locking Value (IPLV) for a single channel pair. The IPLV is the absolute imaginary part of the mean complex phase difference:

IPLV = |Im( mean( exp(i·Δφ) ) )|

Unlike the standard PLV, the imaginary component is insensitive to spurious zero-lag synchrony introduced by volume conduction, making it more specific to true neural coupling.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

Named tuple containing:

- `ipl::Float64`: PLV value
- `sd::Vector{Float64}`: signal difference (s1 - s2)
- `phd::Vector{Float64}`: phase difference (s1 - s2)
- `s1ph::Vector{Float64}`: instantaneous phase of s1
- `s2ph::Vector{Float64}`: instantaneous phase of s2

# Reference

Aydore S, Pantazis D, Leahy RM. A note on the phase locking value and its properties. NeuroImage. 2013 July;74:231–44.
"""
function iplv(
    s1::AbstractVector, s2::AbstractVector
)::@NamedTuple{
    ipl::Float64,
    sd::Vector{Float64},
    phd::Vector{Float64},
    s1ph::Vector{Float64},
    s2ph::Vector{Float64}
}

    @assert length(s1) == length(s2) "Both signals must have the same length."

    # instantaneous phases via Hilbert transform
    h1 = htransform(s1)
    h2 = htransform(s2)
    s1ph = h1.ph
    s2ph = h2.ph

    # signal and phase differences
    sd  = s1 .- s2
    phd = s1ph .- s2ph

    # IPLV
    ipl = abs(imag(mean(Base.cis.(phd)))) # cis.(phd) = exp.(im .* phd)

    return (ipl = ipl, sd = sd, phd = phd, s1ph = s1ph, s2ph = s2ph)

end

"""
    iplv(obj1, obj2; <keyword arguments>)

Calculate Imaginary Phase Locking Value (IPLV) between channel-matched pairs. The IPLV is the absolute imaginary part of the mean complex phase difference:

IPLV = |Im( mean( exp(i·Δφ) ) )|

Unlike the standard PLV, the imaginary component is insensitive to spurious zero-lag synchrony introduced by volume conduction, making it more specific to true neural coupling.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

Named tuple containing:

- `ipl::Matrix{Float64}`: IPLV value, shape `(channels, epochs)`
- `sd::Array{Float64, 3}`: signal difference (s1 - s2), shape `(channels, samples, epochs)`
- `phd::Array{Float64, 3}`: phase difference (s1 - s2), shape `(channels, samples, epochs)`
- `s1ph::Array{Float64, 3}`: signal 1 phases, shape `(channels, samples, epochs)`
- `s2ph::Array{Float64, 3}`: signal 2 phases, shape `(channels, samples, epochs)`

# Reference

Aydore S, Pantazis D, Leahy RM. A note on the phase locking value and its properties. NeuroImage. 2013 July;74:231–44.
"""
function iplv(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2))
)::@NamedTuple{
    ipl::Matrix{Float64},
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
    ipl = zeros(ch_n, ep_n)
    sd = zeros(ch_n, epoch_len(obj1), ep_n)
    phd = zeros(ch_n, epoch_len(obj1), ep_n)
    s1ph = zeros(ch_n, epoch_len(obj1), ep_n)
    s2ph = zeros(ch_n, epoch_len(obj1), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        iplv_data = iplv(
            @view(obj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(obj2.data[ch2[ch_idx], :, ep2[ep_idx]]),
        )
        ipl[ch_idx, ep_idx] = iplv_data.ipl
        sd[ch_idx, :, ep_idx] = iplv_data.sd
        phd[ch_idx, :, ep_idx]  = iplv_data.phd
        s1ph[ch_idx, :, ep_idx] = iplv_data.s1ph
        s2ph[ch_idx, :, ep_idx] = iplv_data.s2ph
    end

    return (ipl = ipl, sd = sd, phd = phd, s1ph = s1ph, s2ph = s2ph)

end

"""
    iplv(obj; <keyword arguments>)

Calculate Imaginary Phase Locking Value (IPLV) for all channel pairs. The IPLV is the absolute imaginary part of the mean complex phase difference:

IPLV = |Im( mean( exp(i·Δφ) ) )|

Unlike the standard PLV, the imaginary component is insensitive to spurious zero-lag synchrony introduced by volume conduction, making it more specific to true neural coupling.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `ipl::Array{Float64, 3}`: IPLV matrix, shape `(channels, channels, epochs)`; symmetric with zeros on the diagonal

# Reference

Aydore S, Pantazis D, Leahy RM. A note on the phase locking value and its properties. NeuroImage. 2013 July;74:231–44.
"""
function iplv(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    # pre-allocate symmetric output: channels × channels × epochs
    # only the lower triangle is computed
    ipl = zeros(ch_n, ch_n, ep_n)

    # compute lower triangle (ch_idx2 < ch_idx1); diagonal is zero by definition
    # the outer loop is parallelized over epochs
    # the inner two loops over channel pairs are not nested @threads (no nesting issue here)
    @inbounds Threads.@threads for ep_idx in 1:ep_n
        for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1 - 1
                ipl[ch_idx1, ch_idx2, ep_idx] = iplv(
                    @view(obj.data[ch[ch_idx1], :, ep_idx]),
                    @view(obj.data[ch[ch_idx2], :, ep_idx]),
                ).ipl
            end
        end
    end

    # mirror lower triangle to upper triangle
    ipl = _copy_lt2ut(ipl)

    return ipl

end
