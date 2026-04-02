export pli

"""
    pli(s1, s2)

Calculate Phase Locking Index (PLI) for two 1-D signal vectors.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

Named tuple:

- `pv::Float64`: PLI value
- `sd::Vector{Float64}`: signal difference (s1 - s2)
- `phd::Vector{Float64}`: phase difference (s1 - s2)
- `s1ph::Vector{Float64}`: signal 1 phase
- `s2ph::Vector{Float64}`: signal 2 phase

# References

 1. Stam CJ, Nolte G, Daffertshofer A. Phase lag index: assessment of functional connectivity from multi channel EEG and MEG with diminished bias from common sources. Hum Brain Mapp. 2007 Nov;28(11):1178-93.
 2. Aydore S, Pantazis D, Leahy RM. A note on the phase locking value and its properties. NeuroImage. 2013 July;74:231–44.
"""
function pli(
    s1::AbstractVector,
    s2::AbstractVector,
)::@NamedTuple{
    pv::Float64,
    sd::Vector{Float64},
    phd::Vector{Float64},
    s1ph::Vector{Float64},
    s2ph::Vector{Float64},
}

    # validate
    length(s1) == length(s2) ||
        throw(ArgumentError("Both signals must have the same length."))

    # get instatenous phases
    ht_data1 = htransform(s1)
    ht_data2 = htransform(s2)
    s1ph = ht_data1.ph
    s2ph = ht_data2.ph

    # signal difference
    sd = s1 - s2

    # phase differences
    phd = s1ph - s2ph

    # PLI
    pv = abs(mean(sign.(phd)))

    return (; pv, sd, phd, s1ph, s2ph)
end

"""
    pli(obj1, obj2; <keyword arguments>)

Calculate Phase Locking Index (PLI) for two NEURO objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s) in `obj1`
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s) in `obj2`
- `ep1::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj1))`: epoch number(s) in `obj1`
- `ep2::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj2))`: epoch number(s) in `obj2`

# Returns

Named tuple:

- `pv::Matrix{Float64}`: PLI value
- `sd::Array{Float64, 3}`: signal difference (s1 - s2)
- `phd::Array{Float64, 3}`: phase difference (s1 - s2)
- `s1ph::Array{Float64, 3}`: signal 1 phase
- `s2ph::Array{Float64, 3}`: signal 2 phase
"""
function pli(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj2)),
)::@NamedTuple{
    pv::Matrix{Float64},
    sd::Array{Float64, 3},
    phd::Array{Float64, 3},
    s1ph::Array{Float64, 3},
    s2ph::Array{Float64, 3},
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 =
        exclude_bads ? get_channel(obj1; ch = ch1, exclude = "bad") :
        get_channel(obj1; ch = ch1, exclude = "")
    ch2 =
        exclude_bads ? get_channel(obj2; ch = ch2, exclude = "bad") :
        get_channel(obj2; ch = ch2, exclude = "")
    isempty(ch1) && throw(ArgumentError("No channels selected."))
    isempty(ch2) && throw(ArgumentError("No channels selected."))
    length(ch1) == length(ch2) ||
        throw(
            ArgumentError(
                "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal.",
            ),
        )

    # validate
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) ||
        throw(
            ArgumentError(
                "Lengths of ep1 ($(length(ep1)) and ep2 ($(length(ep2)) must be equal.",
            ),
        )
    epoch_len(obj1) == epoch_len(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))
    ep1 = _n2v(ep1)
    ep2 = _n2v(ep2)

    # number of channels
    ch_n = length(ch1)
    # number of epochs
    ep_n = length(ep1)

    # pre-allocate output
    pv = zeros(ch_n, ep_n)
    sd = zeros(ch_n, epoch_len(obj1), ep_n)
    phd = zeros(ch_n, epoch_len(obj1), ep_n)
    s1ph = zeros(ch_n, epoch_len(obj1), ep_n)
    s2ph = zeros(ch_n, epoch_len(obj1), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        pli_data = pli(
            @view(viewobj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(viewobj2.data[ch2[ch_idx], :, ep2[ep_idx]])
        )
        pv[ch_idx, ep_idx] = pli_data.pv
        sd[ch_idx, :, ep_idx] = pli_data.sd
        phd[ch_idx, :, ep_idx] = pli_data.phd
        s1ph[ch_idx, :, ep_idx] = pli_data.s1ph
        s2ph[ch_idx, :, ep_idx] = pli_data.s2ph
    end

    return (; pv, sd, phd, s1ph, s2ph)
end

"""
    pli(obj; <keyword arguments>)

Calculate Phase Locking Index (PLI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Array{Float64, 3}`: PLI value
"""
function pli(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    # pre-allocate output
    pv = zeros(ch_n, ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 in 1:ch_idx1
            pv[ch_idx1, ch_idx2, ep_idx] = pli(
                @view(obj.data[ch[ch_idx1], :, ep_idx]),
                @view(obj.data[ch[ch_idx2], :, ep_idx])
            ).pv
        end
    end

    # mirror the lower triangle to the upper triangle to produce the full symmetric matrix
    return _copy_lt2ut(pv)
end
