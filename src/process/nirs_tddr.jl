export tddr
export tddr!

# TDDR of a single signal (Fishburn et al. 2019)
function _tddr(s::AbstractVector, fs::Int64)::Vector{Float64}
    x = Float64.(s)
    x_mean = mean(x)
    x = x .- x_mean

    # low-pass (< 0.5 Hz) component is corrected; the high-frequency component is kept
    fc = 0.5
    x_low = fc < fs / 2 ?
        filtfilt(digitalfilter(Lowpass(fc), Butterworth(3); fs = fs), x) : copy(x)
    x_high = x .- x_low

    # robust (Tukey biweight) estimate of the mean of the temporal derivative: derivatives
    # far from the bulk of the distribution (motion artifacts) get small weights
    d = diff(x_low)
    w = ones(length(d))
    tune = 4.685
    tol = sqrt(eps(Float64))
    mu = Inf
    for _ = 1:50
        mu0 = mu
        mu = sum(w .* d) / sum(w)
        dev = abs.(d .- mu)
        sigma = 1.4826 * median(dev)
        sigma > 0 || break
        r = dev ./ (sigma * tune)
        w = ((1 .- r .^ 2) .* (r .< 1)) .^ 2
        abs(mu - mu0) < tol * max(abs(mu), abs(mu0)) && break
    end

    # weighted derivative, integrated back to the signal
    x_corr = cumsum(vcat(0.0, w .* (d .- mu)))
    x_corr .-= mean(x_corr)

    return x_corr .+ x_high .+ x_mean
end

"""
    tddr(obj; <keyword arguments>)

Correct motion artifacts in NIRS signals using temporal derivative distribution repair (TDDR).

TDDR assumes that the temporal derivative of the signal (below 0.5 Hz) is normally distributed, while motion artifacts produce outlying derivatives. The derivatives are weighted with a robust (Tukey biweight) estimator, so that outlying derivatives (spikes and baseline shifts) are down-weighted, and the signal is reconstructed from the weighted derivatives; the high-frequency component (> 0.5 Hz) is kept unchanged. TDDR is parameter-free and is applied to OD or concentration channels (not to the raw intensity), separately in each epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj; type = "nirs_od")`: channel name(s), default is all OD channels

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object

# References

1. Fishburn FA, Ludlum RS, Vaidya CJ, Medvedev AV. Temporal Derivative Distribution Repair (TDDR): A motion correction method for fNIRS. NeuroImage. 2019;184:171–179.
"""
function tddr(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = get_channel(obj; type = "nirs_od"),
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    tddr!(obj_new; ch = ch)

    return obj_new
end

"""
    tddr!(obj; <keyword arguments>)

Correct motion artifacts in NIRS signals using temporal derivative distribution repair (TDDR). See `tddr()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=get_channel(obj; type = "nirs_od")`: channel name(s), default is all OD channels

# Returns

- `Nothing`
"""
function tddr!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = get_channel(obj; type = "nirs_od"),
)::Nothing
    # validate
    _check_datatype(obj, "nirs")
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    any(==("nirs_int"), obj.header.recording[:channel_type][ch]) &&
        _warn("TDDR should be applied to OD or concentration channels, not to the raw intensity.")

    fs = sr(obj)
    Threads.@threads :static for idx in CartesianIndices((length(ch), nepochs(obj)))
        c, ep = ch[idx[1]], idx[2]
        obj.data[c, :, ep] = _tddr(@view(obj.data[c, :, ep]), fs)
    end
    push!(obj.history, "tddr(obj; ch=$(labels(obj)[ch]))")

    return nothing
end
