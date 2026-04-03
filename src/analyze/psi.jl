export psi

"""
    psi(s1, s2; <keyword arguments>)

Calculate Phase Slope Index (PSI) for two 1-D signal vectors.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `flim::Tuple{Real, Real}=(1, fs / 2 - 1))`: frequency limits

# Returns

- `Tuple{Float64, Float64}`: PSI value (signal1 -> signal2, signal2 -> signal1)

# References

 1. Nolte, G., Ziehe, A., Nikulin, V. V., Schlögl, A., Krämer, N., Brismar, T., & Müller, K.-R. (2008). Robustly Estimating the Flow Direction of Information in Complex Physical Systems. Physical Review Letters. 2008; 100(23).
"""
function psi(
    s1::AbstractVector,
    s2::AbstractVector;
    fs::Int64,
    flim::Tuple{Real, Real} = (1, fs / 2 - 1),
)::Tuple{Float64, Float64}
    # validate
    length(s1) == length(s2) ||
        throw(ArgumentError("Both signals must have the same length."))
    _check_tuple(flim, (1, fs / 2 - 1), "flim")

    if flim[1] != round(Int64, flim[1])
        _warn("Lower frequency bound rounded to: $(round(Int64, flim[1])) Hz")
        flim = (round(Int64, flim[1]), flim[end])
    end
    if flim[end] != round(Int64, flim[end])
        _warn("Upper frequency bound rounded to: $(round(Int64, flim[end])) Hz")
        flim = (flim[1], round(Int64, flim[end]))
    end

    fs > length(s1) && (fs = length(s1))
    seglen = fs             # segment length (1 second) or signal length
    nboot = 256             # number of bootstrap iterations
    method = "boostrap"     # standard deviation estimation method
    detrend = true          # performs a 0th-order detrend across raw segments

    pv, _ = PhaseSlopeIndex.data2psi(
        [s1 s2], seglen; nboot = nboot, method = method, detrend = detrend,
        freqlist = Int64(flim[1]):1:Int64(flim[end]),
    )
    pv = (pv[1, 2], pv[2, 1])

    return pv
end

"""
    psi(obj1, obj2; <keyword arguments>)

Calculate Phase Slope Index (PSI) for two NEURO objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s) in `obj1`
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s) in `obj2`
- `ep1::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj1))`: epoch number(s) in `obj1`
- `ep2::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj2))`: epoch number(s) in `obj2`
- `flim::Tuple{Real, Real}=(1, sr(obj1) / 2 - 1))`: frequency limits

# Returns

- `Matrix{Float64}`: PSI value
"""
function psi(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj2)),
    flim::Tuple{Real, Real} = (1, sr(obj1) / 2 - 1),
)::Matrix{Tuple{Float64, Float64}}
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

    ch_n = length(ch1)
    ep_n = length(ep1)

    # pre-allocate output
    pv = Matrix{Tuple{Float64, Float64}}(undef, ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        pv[ch_idx, ep_idx] = psi(
            @view(obj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(obj2.data[ch2[ch_idx], :, ep2[ep_idx]]),
            fs = sr(obj1),
            flim = flim,
        )
    end

    return pv
end

"""
    psi(obj; <keyword arguments>)

Calculate Phase Slope Index (PSI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flim::Tuple{Real, Real}=(1, sr(obj) / 2 - 1))`: frequency limits

# Returns

- `Array{Tuple{Float64, Float64}, 3}`: PSI value
"""
function psi(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    flim::Tuple{Real, Real} = (1, sr(obj) / 2 - 1),
)::Array{Tuple{Float64, Float64}, 3}
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
    pv = Array{Tuple{Float64, Float64}}(undef, ch_n, ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 in 1:ch_idx1
            pv[ch_idx1, ch_idx2, ep_idx] = psi(
                @view(obj.data[ch[ch_idx1], :, ep_idx]),
                @view(obj.data[ch[ch_idx2], :, ep_idx]),
                fs = sr(obj), flim = flim,
            )
        end
    end

    # mirror the lower triangle to the upper triangle to produce the full symmetric matrix
    return _copy_lt2ut(pv)
end
