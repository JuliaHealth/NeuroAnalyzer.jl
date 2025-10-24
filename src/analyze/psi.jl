export psi

"""
    psi(s1, s2; <keyword arguments>)

Calculate Phase Slope Index (PSI).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(1, fs / 2 - 1))`: frequency bounds

# Returns

- `pv::Tuple{Float64, Float64}`: PSI value (signal1 - signal2, signal2- signal1)

# Source

1. Nolte, G., Ziehe, A., Nikulin, V. V., Schlögl, A., Krämer, N., Brismar, T., & Müller, K.-R. (2008). Robustly Estimating the Flow Direction of Information in Complex Physical Systems. Physical Review Letters. 2008; 100(23).
"""
function psi(s1::AbstractVector, s2::AbstractVector; fs::Int64, frq_lim::Tuple{Real, Real}=(1, fs / 2 - 1))::Tuple{Float64, Float64}

    @assert length(s1) == length(s2) "Both signals must have the same length."
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))
    @assert frq_lim[1] >= 1.0 "Lower frequency bound must be ≥ 1."
    @assert frq_lim[2] < fs/2 "Upper frequency bound must be ≤ $(round(Int64, fs / 2 - 1))."

    if frq_lim[1] != round(Int64, frq_lim[1])
        _warn("Lower frequency bound rounded to: $(round(Int64, frq_lim[1])) Hz")
        frq_lim = (round(Int64, frq_lim[1]), frq_lim[end])
    end
    if frq_lim[end] != round(Int64, frq_lim[end])
        _warn("Upper frequency bound rounded to: $(round(Int64, frq_lim[end])) Hz")
        frq_lim = (frq_lim[1], round(Int64, frq_lim[end]))
    end

    fs > length(s1) && (fs = length(s1))
    seglen = fs             # segment length (1 second) or signal length
    nboot = 256             # number of bootstrap iterations
    method = "boostrap"     # standard deviation estimation method
    detrend = true          # performs a 0th-order detrend across raw segments

    pv, _ = data2psi([s1 s2], seglen; nboot=nboot, method=method, detrend=detrend, freqlist=Int64(frq_lim[1]):1:Int64(frq_lim[end]))
    pv = (pv[1, 2], pv[2, 1])

    return pv

end

"""
    psi(obj1, obj2; <keyword arguments>)

Calculate Phase Slope Index (PSI).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `frq_lim::Tuple{Real, Real}=(1, sr(obj1) / 2 - 1))`: frequency bounds

# Returns

- `pv::Matrix{Float64}`: PSI value
"""
function psi(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)), frq_lim::Tuple{Real, Real}=(1, sr(obj1) / 2 - 1))::Matrix{Tuple{Float64, Float64}}

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

    pv = Matrix{Tuple{Float64, Float64}}(undef, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pv[ch_idx, ep_idx] = @views psi(obj1.data[ch1[ch_idx], :, ep1[ep_idx]], obj2.data[ch2[ch_idx], :, ep2[ep_idx]], fs=sr(obj1), frq_lim=frq_lim)
        end
    end

    return pv

end

"""
    psi(obj; <keyword arguments>)

Calculate Phase Slope Index (PSI).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `frq_lim::Tuple{Real, Real}=(1, sr(obj) / 2 - 1))`: frequency bounds

# Returns

- `pv::Array{Tuple{Float64, Float64}, 3}`: PSI value
"""
function psi(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, frq_lim::Tuple{Real, Real}=(1, sr(obj) / 2 - 1))::Array{Tuple{Float64, Float64}, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    isa(ch, Int64) && (ch = [ch])

    pv = Array{Tuple{Float64, Float64}}(undef, ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                pv[ch_idx1, ch_idx2, ep_idx] = @views psi(obj.data[ch[ch_idx1], :, ep_idx], obj.data[ch[ch_idx2], :, ep_idx], fs=sr(obj), frq_lim=frq_lim)
            end
        end
    end

    pv = _copy_lt2ut(pv)

    return pv

end
