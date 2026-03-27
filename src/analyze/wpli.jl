export wpli

"""
    wpli(s1, s2; <keyword arguments>)

Calculate weighted PLI (Phase Locking Index) for two 1-D signal vectors.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `debiased::Bool=false`: if `true`, calculate debiased wPLI

# Returns

Named tuple:

- `pv::Float64`: wPLI value
- `sd::Vector{Float64}`: signal difference (s1 - s2)
- `phd::Vector{Float64}`: phase difference (s1 - s2)
- `s1ph::Vector{Float64}`: signal 1 phase
- `s2ph::Vector{Float64}`: signal 2 phase
"""
function wpli(
    s1::AbstractVector,
    s2::AbstractVector;
    debiased::Bool = false,
)::@NamedTuple{
    pv::Float64,
    sd::Vector{Float64},
    phd::Vector{Float64},
    s1ph::Vector{Float64},
    s2ph::Vector{Float64},
}
    length(s1) == length(s2) ||
        throw(ArgumentError("Both signals must have the same length."))

    # CPSD
    n = length(s1)
    ss1 = fft(detrend(s1; type = :mean) .* DSP.hanning(n)) / n
    ss2 = fft(detrend(s2; type = :mean) .* DSP.hanning(n)) / n
    pxy = conj.(ss1) .* ss2
    im_pxy = imag.(pxy)

    pli_data = pli(s1, s2)
    sd = pli_data.sd
    phd = pli_data.phd
    s1ph = pli_data.s1ph
    s2ph = pli_data.s2ph

    # wPLI
    num = sum(abs.(im_pxy) .* sign.(im_pxy))
    denom = sum(abs.(im_pxy))
    sum_sq = sum(im_pxy .^ 2)

    if debiased
        pv = (num^2 - sum_sq) / (denom^2 - sum_sq)
    else
        pv = num / denom
    end

    return (; pv, sd, phd, s1ph, s2ph)
end

"""
    wpli(obj1, obj2; <keyword arguments>)

Calculate weighted PLI (Phase Locking Index) for two NEURO objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s) in `obj1`
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s) in `obj2`
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s) in `obj1`
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s) in `obj2`
- `debiased::Bool=false`: if `true`, calculate debiased wPLI

# Returns

Named tuple:

- `pv::Matrix{Float64}`: PLI value
- `sd::Array{Float64, 3}`: signal difference (s1 - s2)
- `phd::Array{Float64, 3}`: phase difference (s1 - s2)
- `s1ph::Array{Float64, 3}`: signal 1 phase
- `s2ph::Array{Float64, 3}`: signal 2 phase
"""
function wpli(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
    debiased::Bool = false,
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
    length(ch1) == length(ch2) ||
        throw(
            ArgumentError(
                "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal.",
            ),
        )

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
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    ch_n = length(ch1)
    ep_n = length(ep1)

    pv = zeros(ch_n, ep_n)
    sd = zeros(ch_n, epoch_len(obj1), ep_n)
    phd = zeros(ch_n, epoch_len(obj1), ep_n)
    s1ph = zeros(ch_n, epoch_len(obj1), ep_n)
    s2ph = zeros(ch_n, epoch_len(obj1), ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        wpli_data = wpli(
            @view(obj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(obj2.data[ch2[ch_idx], :, ep2[ep_idx]]),
            debiased = debiased,
        )
        pv[ch_idx, ep_idx] = wpli_data.pv
        sd[ch_idx, :, ep_idx] = wpli_data.sd
        phd[ch_idx, :, ep_idx] = wpli_data.phd
        s1ph[ch_idx, :, ep_idx] = wpli_data.s1ph
        s2ph[ch_idx, :, ep_idx] = wpli_data.s2ph
    end

    return (; pv, sd, phd, s1ph, s2ph)
end

"""
    wpli(obj; <keyword arguments>)

Calculate weighted PLI (Phase Locking Index) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `debiased::Bool=false`: if `true`, calculate debiased wPLI

# Returns

- `Array{Float64, 3}`: wPLI value
"""
function wpli(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    debiased::Bool = false,
)::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isa(ch, Int64) && (ch = [ch])

    ch_n = length(ch)
    ep_n = nepochs(obj)

    pv = zeros(ch_n, ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 in 1:ch_idx1
            pv[ch_idx1, ch_idx2, ep_idx] = wpli(
                @view(obj.data[ch[ch_idx1], :, ep_idx]),
                @view(obj.data[ch[ch_idx2], :, ep_idx]),
                debiased = debiased,
            ).pv
        end
    end

    # mirror the lower triangle to the upper triangle to produce the full symmetric matrix
    return _copy_lt2ut(pv)
end
