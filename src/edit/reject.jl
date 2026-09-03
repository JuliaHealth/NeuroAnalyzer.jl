export channel_reject
export channel_reject!
export epoch_reject

# -------------------------------------------------------------------------------------------------
# Internal detection helpers
# Each function takes a (ch_n × samples) matrix and returns a Bool vector where true = bad channel.
# -------------------------------------------------------------------------------------------------

"""
Return true for each channel whose RMSE vs the cross-channel median lies
outside the 95 % confidence interval of the per-channel RMSE distribution.
"""
function detect_rmse(s::AbstractMatrix)::Vector{Bool}
    ch_n    = size(s, 1)
    ch_m    = vec(median(s; dims = 1))   # cross-channel median reference signal
    rmse_ch = [rmse(@view(s[i, :]), ch_m) for i = 1:ch_n]

    # compute CI once; was called twice per iteration (wasted allocation)
    ci = HypothesisTests.confint(OneSampleTTest(rmse_ch))
    return [rmse_ch[i] < ci[1] || rmse_ch[i] > ci[2] for i = 1:ch_n]
end

"""
Return true for each channel whose RMSD vs the cross-channel median lies
outside the 95 % CI of the per-channel RMSD distribution.
"""
function detect_rmsd(s::AbstractMatrix)::Vector{Bool}
    ch_n    = size(s, 1)
    ch_m    = vec(median(s; dims = 1))
    rmsd_ch = [Distances.rmsd(@view(s[i, :]), ch_m) for i = 1:ch_n]
    ci      = HypothesisTests.confint(OneSampleTTest(rmsd_ch))
    return [rmsd_ch[i] < ci[1] || rmsd_ch[i] > ci[2] for i = 1:ch_n]
end

"""
Return true for each channel whose Euclidean distance to the cross-channel
median lies outside the 95 % CI of the per-channel distance distribution.
"""
function detect_euclid(s::AbstractMatrix)::Vector{Bool}
    ch_n  = size(s, 1)
    ch_m  = vec(median(s; dims = 1))
    ed_ch = [Distances.euclidean(@view(s[i, :]), ch_m) for i = 1:ch_n]
    ci    = HypothesisTests.confint(OneSampleTTest(ed_ch))
    return [ed_ch[i] < ci[1] || ed_ch[i] > ci[2] for i = 1:ch_n]
end

"""
Return true for each channel whose peak amplitude exceeds ±`amp_t`.
"""
function detect_amp(s::AbstractMatrix; amp_t::Real = 400.0)::Vector{Bool}
    return [maximum(s[i, :]) > amp_t || minimum(s[i, :]) < -amp_t for i = 1:size(s, 1)]
end

"""
Return true for each channel with excessive peak-to-peak variation.

Sliding-window means are computed, their differences are rounded to the
nearest 100, and a channel is flagged if any rounded difference falls outside
`[mean − p_z×std, mean + p_z×std]` where `p_z = Φ⁻¹(p)`.
"""
function detect_p2p(s::AbstractMatrix; w::Int64 = 10, p::Float64 = 0.95)::Vector{Bool}
    w < size(s, 2) || throw(ArgumentError("w must be < $(size(s, 2))."))

    ch_n    = size(s, 1)
    bad_chs = zeros(Bool, ch_n)
    p_z     = quantile(Distributions.Normal(), p)

    for ch_idx = 1:ch_n
        v = @view s[ch_idx, :]
        # sliding-window means
        sm = [mean(@view v[idx:(idx + w)]) for idx = 1:w:(length(v) - w)]
        p2p = round.(diff(sm); digits = -2)

        s_m = mean(v)
        s_s = std(v)
        bad_chs[ch_idx] = any(p2p .> s_m + p_z * s_s) || any(p2p .< s_m - p_z * s_s)
    end

    return bad_chs
end

"""
Return true for each channel with at least two 10-sample windows where the
absolute z-scored signal exceeds the absolute z-scored TKEO by more than the
z-score corresponding to `p`.
"""
function detect_tkeo(
    s::AbstractMatrix,
    t::AbstractVector;
    tkeo_method::Symbol = :pow,
    p::Float64 = 0.95,
)::Vector{Bool}
    ch_n    = size(s, 1)
    bad_chs = zeros(Bool, ch_n)
    thresh  = cl2z(p)

    for ch_idx = 1:ch_n
        stkeo    = tkeo(@view(s[ch_idx, :]), t; method = tkeo_method)
        z_signal = vec(NeuroAnalyzer.zscore(s[ch_idx, :]))
        z_tkeo   = vec(NeuroAnalyzer.zscore(stkeo))

        # scan in 10-sample windows; count windows where signal z-score
        # substantially exceeds TKEO z-score (indicative of an artefact)
        w = max(1, length(stkeo) ÷ 10)
        bad_windows = 0
        for idx = 1:w:length(stkeo)
            i2 = min(idx + w - 1, length(stkeo))
            count(abs.(z_signal[idx:i2]) .- abs.(z_tkeo[idx:i2]) .> thresh) > 1 &&
                (bad_windows += 1)
        end

        bad_chs[ch_idx] = bad_windows > 1
    end

    return bad_chs
end

"""
Return true for each channel whose correlation with its nearest spatial
neighbour (outliers removed via RANSAC) falls below `ransac_r` in more than
`ransac_tr` proportion of sliding windows.
"""
function detect_ransac(
    s::AbstractMatrix;
    loc_x::Vector{Float64},
    loc_y::Vector{Float64},
    w::Int64 = 10,
    ransac_r::Float64 = 0.8,
    ransac_tr::Float64 = 0.4,
    ransac_t::Float64 = 100.0,
)::Vector{Bool}   # was: missing return type
    ch_n    = size(s, 1)
    bad_chs = zeros(Bool, ch_n)

    # bild pairwise Euclidean distance matrix between electrode positions
    d = zeros(ch_n, ch_n)
    @inbounds for i = 1:ch_n, j = 1:ch_n
        d[i, j] = euclidean([loc_x[i], loc_y[i]], [loc_x[j], loc_y[j]])
    end
    # prevent a channel from being its own nearest neighbour
    d[d .== 0] .= Inf

    @inbounds for ch_idx = 1:ch_n
        _, nearest_idx = findmin(d[ch_idx, :])
        y = @view s[ch_idx, :]
        x = @view s[nearest_idx, :]

        # remove RANSAC outliers before computing windowed correlation
        df  = DataFrame(:y => remove_dc(y), :x => remove_dc(x))
        reg = createRegressionSetting(@formula(y ~ x), df)
        out = ransac(reg; t = ransac_t, k = 128)["outliers"]
        idx = setdiff(1:length(x), out)
        xc  = x[idx]
        yc  = y[idx]

        # windowed Pearson correlation
        c = [cor(@view(xc[i:(i + w)]), @view(yc[i:(i + w)])) for i = 1:w:(length(xc) - w)]
        bad_chs[ch_idx] = sum(c .< ransac_r) / length(c) > ransac_tr
    end

    return bad_chs
end

"""
    channel_reject(obj; <keyword arguments>)

Detect bad channels using one or more artifact-detection methods.

A channel is marked bad if any of the selected methods flags it in at least one epoch. The returned mask has one entry per *selected* channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]`: one or more detection methods (see below); default: all available methods.
- `w::Int64=sr(obj)`: sliding-window width in samples (default = 1 s)
- `flat_tol::Float64=0.1`: flatness tolerance; a segment is flat if `|diff(window_mean)| < flat_tol`.
- `flat_fr::Float64=0.3`: maximum acceptable fraction of flat windows before marking a channel as bad
- `p::Float64=0.99`: probability threshold for `:p2p`, `:tkeo`, and `:z` methods; interpreted as a percentile (converted internally to a z-score)
- `tc::Float64=0.2`: reserved threshold parameter (not currently used by all methods)
- `tkeo_method::Symbol=:pow`: TKEO variant; see `tkeo()` for options
- `z::Real=3`: z-score threshold for `:kurt` and `:z` methods; must be > 0
- `ransac_r::Float64=0.8`: minimum acceptable Pearson r between a channel and its nearest neighbor (`:ransac` method)
- `ransac_tr::Float64=0.4`: maximum acceptable fraction of windows below `ransac_r` before marking a channel bad
- `ransac_t::Float64=100.0`: RANSAC inlier distance threshold (in signal units)
- `amp_t::Real=400.0`: amplitude rejection threshold; channels with `max > +amp_t` or `min < −amp_t` are flagged

# Detection methods

- `:flat`: channels with a large proportion of flat windows
- `:rmse`: channels whose RMSE vs the median reference is outside the 95 % CI
- `:rmsd`: same as `:rmse` using RMSD
- `:euclid`: same using Euclidean distance
- `:var`: channels with IQR-outlier variance
- `:p2p`: channels with excessive peak-to-peak variation
- `:tkeo`: channels where z-scored TKEO diverges from z-scored signal
- `:kurt`: channels with z-scored kurtosis exceeding `z`
- `:z`: channels with a large proportion of samples above the z-score threshold
- `:ransac`: channels poorly correlated with their nearest spatial neighbor
- `:amp`: channels exceeding ±`amp_t`

# Returns

- `Vector{Bool}`: bad-channel mask (length = number of selected channels); `true` = bad
"""
function channel_reject(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Union{Symbol, Vector{Symbol}} = [
        :flat,
        :rmse,
        :rmsd,
        :euclid,
        :var,
        :p2p,
        :tkeo,
        :kurt,
        :z,
        :ransac,
        :amp,
    ],
    w::Int64 = sr(obj),
    flat_tol::Float64 = 0.1,
    flat_fr::Float64 = 0.3,
    p::Float64 = 0.99,
    tc::Float64 = 0.2,
    tkeo_method::Symbol = :pow,
    z::Real = 3,
    ransac_r::Float64 = 0.8,
    ransac_tr::Float64 = 0.4,
    ransac_t::Float64 = 100.0,
    amp_t::Real = 400.0,
)::Vector{Bool}
    # validate
    _in(p, (0, 1), "p")
    _in(tc, (0, 1), "tc")
    _in(ransac_r, (0, 1), "ransac_r")
    _in(ransac_tr, (0, 1), "ransac_tr")

    # validate
    typeof(method) != Vector{Symbol} && (method = [method])
    for m in method
        _check_var(
            m,
            [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp],
            "method",
        )
    end

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    ch_list = labels(obj)[ch]

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # epoch length
    ep_len = epoch_len(obj)

    # methods that need more than one channel
    for m in (:rmse, :rmsd, :euclid)
        m in method && ch_n < 2 && throw(ArgumentError(":$m requires > 1 channel."))
    end

    # pre-allocate output: bad-channel mask over all channels; we will index by ch at the end
    bc = zeros(Bool, nchannels(obj))

    # -----------------------------------------------------------------------
    if :flat in method
        w < ep_len || throw(ArgumentError("w must be < $ep_len."))
        _info("Using :flat method")

        bad_chs = zeros(Bool, ch_n, ep_n)   # was: 1-D, then indexed as 2-D
        n_samples = size(obj.data, 2)

        Threads.@threads :static for linear_idx = 1:(ch_n * ep_n)
            ci = (linear_idx - 1) % ch_n + 1
            ei = (linear_idx - 1) ÷ ch_n + 1
            sm = [
                mean(@view obj.data[ch[ci], idx_w:(idx_w + w), ei])
                for idx_w = 1:w:(n_samples - w)
            ]
            r = count(abs.(diff(sm)) .< flat_tol) / length(sm)
            @inbounds bad_chs[ci, ei] = r > flat_fr
        end

        # a channel is bad if it was flat in any epoch
        bc[ch] = bc[ch] .|| vec(any(bad_chs; dims = 2))
    end

    if :rmse in method
        _info("Using :rmse method")
        @inbounds for ep_idx = 1:ep_n
            bad_chs = detect_rmse(@view(obj.data[ch, :, ep_idx]))
            bc[ch] = bc[ch] .|| bad_chs
        end
    end

    # -----------------------------------------------------------------------
    if :rmse in method
        _info("Using :rmse method")
        for ep_idx = 1:ep_n
            bad_chs = detect_rmse(@view obj.data[ch, :, ep_idx])
            bc[ch]  = bc[ch] .|| bad_chs
        end
    end

    # -----------------------------------------------------------------------
    if :rmsd in method
        _info("Using :rmsd method")
        for ep_idx = 1:ep_n
            bad_chs = detect_rmsd(@view obj.data[ch, :, ep_idx])
            bc[ch]  = bc[ch] .|| bad_chs
        end
    end

    # -----------------------------------------------------------------------
    if :euclid in method
        _info("Using :euclid method")
        for ep_idx = 1:ep_n
            bad_chs = detect_euclid(@view obj.data[ch, :, ep_idx])
            bc[ch]  = bc[ch] .|| bad_chs
        end
    end

    # -----------------------------------------------------------------------
    if :var in method
        _info("Using :var method")
        s_v = var(obj.data[ch, :, :]; dims = 2) # (ch_n × 1 × ep_n)
        # IQR outlier detection over all (channel, epoch) variance values
        o = reshape(outlier_detect(vec(s_v); method = :iqr), ch_n, ep_n)
        # a channel is bad if it was an outlier in any epoch
        bc[ch] = bc[ch] .|| vec(any(o; dims = 2))
    end

    # -----------------------------------------------------------------------
    if :p2p in method
        _info("Using :p2p method")
        for ep_idx = 1:ep_n
            bad_chs = detect_p2p(@view(obj.data[ch, :, ep_idx]); w = w, p = p)
            bc[ch]  = bc[ch] .|| bad_chs
        end
    end

    # -----------------------------------------------------------------------
    if :tkeo in method
        _info("Using :tkeo method")
        for ep_idx = 1:ep_n
            bad_chs = detect_tkeo(
                obj.data[ch, :, ep_idx], obj.time_pts; tkeo_method = tkeo_method, p = p,
            )
            bc[ch] = bc[ch] .|| bad_chs
        end
    end

    # -----------------------------------------------------------------------
    if :kurt in method
        _info("Using :kurt method")
        z > 0 || throw(ArgumentError("z must be > 0."))

        # Kurtosis per (channel, epoch), then z-scored globally
        k = zeros(ch_n, ep_n)
        Threads.@threads :static for linear_idx = 1:(ch_n * ep_n)
            ci = (linear_idx - 1) % ch_n + 1
            ei = (linear_idx - 1) ÷ ch_n + 1
            @inbounds k[ci, ei] = kurtosis(@view obj.data[ch[ci], :, ei])
        end
        k       = normalize_zscore(k)
        bad_idx = abs.(k) .> z
        bc[ch]  = bc[ch] .|| vec(any(bad_idx; dims = 2))
    end

    # -----------------------------------------------------------------------
    if :z in method
        _info("Using :z method")
        z > 0 || throw(ArgumentError("z must be > 0."))

        # pass 1 - global z-score threshold
        s_z     = normalize_zscore(obj.data[ch, :, :]; bych = false)
        above_z = abs.(s_z) .> z
        k       = zeros(ch_n, ep_n)
        Threads.@threads :static for linear_idx = 1:(ch_n * ep_n)
            ci = (linear_idx - 1) % ch_n + 1
            ei = (linear_idx - 1) ÷ ch_n + 1
            @inbounds k[ci, ei] = count(above_z[ci, :, ei]) / size(above_z, 2)
        end
        bc[ch] = bc[ch] .|| vec(any(k .> p; dims = 2))

        # pass 2 - per-channel z-score threshold (z + 1 to be stricter)
        above_z2 = abs.(normalize_zscore(obj.data[ch, :, :]; bych = false)) .> (z + 1)
        k2       = zeros(ch_n, ep_n)
        Threads.@threads :static for linear_idx = 1:(ch_n * ep_n)
            ci = (linear_idx - 1) % ch_n + 1
            ei = (linear_idx - 1) ÷ ch_n + 1
            @inbounds k2[ci, ei] = count(above_z2[ci, :, ei]) / size(above_z2, 2)
        end
        bc[ch] = bc[ch] .|| vec(any(k2 .> p; dims = 2))
    end

    # -----------------------------------------------------------------------
    if :ransac in method
        _info("Using :ransac method")
        _check_datatype(obj, ["eeg", "seeg", "ecog", "meg"])

        chs = get_channel(obj; type = ["eeg", "seeg", "ecog", "meg", "mag", "grad"])
        length(setdiff(ch_list, chs)) == 0 ||
            throw(ArgumentError("ch must contain only signal channels."))

        chs  = intersect(obj.locs[!, :label], ch_list)
        locs = Base.filter(:label => in(chs), obj.locs)

        for ep_idx = 1:ep_n
            bad_chs = detect_ransac(
                obj.data[ch, :, ep_idx];
                loc_x = locs[!, :loc_x], loc_y = locs[!, :loc_y],
                w = w, ransac_t = ransac_t, ransac_r = ransac_r, ransac_tr = ransac_tr,
            )
            bc[ch] = bc[ch] .|| bad_chs
        end
    end

    # -----------------------------------------------------------------------
    if :amp in method
        _info("Using :amp method")
        for ep_idx = 1:ep_n
            bad_chs = detect_amp(@view(obj.data[ch, :, ep_idx]); amp_t = amp_t)
            bc[ch]  = bc[ch] .|| bad_chs
        end
    end

    return bc[ch]
end

"""
    channel_reject!(obj; <keyword arguments>)

Detect bad channels and update the `:bad_channel` field in the object header in-place.

A channel is marked bad if any of the selected methods flags it in at least one epoch. The returned mask has one entry per *selected* channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]`: one or more detection methods (see below); default: all available methods.
- `w::Int64=sr(obj)`: sliding-window width in samples (default = 1 s)
- `flat_tol::Float64=0.1`: flatness tolerance; a segment is flat if `|diff(window_mean)| < flat_tol`.
- `flat_fr::Float64=0.3`: maximum acceptable fraction of flat windows before marking a channel as bad
- `p::Float64=0.99`: probability threshold for `:p2p`, `:tkeo`, and `:z` methods; interpreted as a percentile (converted internally to a z-score)
- `tc::Float64=0.2`: reserved threshold parameter (not currently used by all methods)
- `tkeo_method::Symbol=:pow`: TKEO variant; see `tkeo()` for options
- `z::Real=3`: z-score threshold for `:kurt` and `:z` methods; must be > 0
- `ransac_r::Float64=0.8`: minimum acceptable Pearson r between a channel and its nearest neighbor (`:ransac` method)
- `ransac_tr::Float64=0.4`: maximum acceptable fraction of windows below `ransac_r` before marking a channel bad
- `ransac_t::Float64=100.0`: RANSAC inlier distance threshold (in signal units)
- `amp_t::Real=400.0`: amplitude rejection threshold; channels with `max > +amp_t` or `min < −amp_t` are flagged

# Detection methods

- `:flat`: channels with a large proportion of flat windows
- `:rmse`: channels whose RMSE vs the median reference is outside the 95 % CI
- `:rmsd`: same as `:rmse` using RMSD
- `:euclid`: same using Euclidean distance
- `:var`: channels with IQR-outlier variance
- `:p2p`: channels with excessive peak-to-peak variation
- `:tkeo`: channels where z-scored TKEO diverges from z-scored signal
- `:kurt`: channels with z-scored kurtosis exceeding `z`
- `:z`: channels with a large proportion of samples above the z-score threshold
- `:ransac`: channels poorly correlated with their nearest spatial neighbor
- `:amp`: channels exceeding ±`amp_t`

# Returns

- `Nothing`
"""
function channel_reject!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Union{Symbol, Vector{Symbol}} = [
        :flat,
        :rmse,
        :rmsd,
        :euclid,
        :var,
        :p2p,
        :tkeo,
        :kurt,
        :z,
        :ransac,
        :amp,
    ],
    w::Int64 = sr(obj),
    flat_tol::Float64 = 0.1,
    flat_fr::Float64 = 0.3,
    p::Float64 = 0.99,
    tc::Float64 = 0.2,
    tkeo_method::Symbol = :pow,
    z::Real = 3,
    ransac_r::Float64 = 0.8,
    ransac_tr::Float64 = 0.4,
    ransac_t::Float64 = 100.0,
    amp_t::Real = 400.0,
)::Nothing
    bc = channel_reject(
        obj;
        ch = ch,
        method = method,
        w = w,
        flat_tol = flat_tol,
        flat_fr = flat_fr,
        p = p,
        tc = tc,
        tkeo_method = tkeo_method,
        z = z,
        ransac_r = ransac_r,
        ransac_tr = ransac_tr,
        ransac_t = ransac_t,
        amp_t = amp_t,
    )
    obj.header.recording[:bad_channel][get_channel(obj; ch = ch)] = bc

    return nothing
end

"""
    epoch_reject(obj; <keyword arguments>)

Detect bad epochs using one or more artifact-detection methods.

An epoch is marked bad if it contains at least `nbad` bad channels according to any of the selected methods.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]`: one or more detection methods (see below); default: all available methods.
- `w::Int64=sr(obj)`: sliding-window width in samples (default = 1 s)
- `flat_tol::Float64=0.1`: flatness tolerance; a segment is flat if `|diff(window_mean)| < flat_tol`.
- `flat_fr::Float64=0.3`: maximum acceptable fraction of flat windows before marking a channel as bad
- `p::Float64=0.99`: probability threshold for `:p2p`, `:tkeo`, and `:z` methods; interpreted as a percentile (converted internally to a z-score)
- `tc::Float64=0.2`: reserved threshold parameter (not currently used by all methods)
- `tkeo_method::Symbol=:pow`: TKEO variant; see `tkeo()` for options
- `z::Real=3`: z-score threshold for `:kurt` and `:z` methods; must be > 0
- `ransac_r::Float64=0.8`: minimum acceptable Pearson r between a channel and its nearest neighbor (`:ransac` method)
- `ransac_tr::Float64=0.4`: maximum acceptable fraction of windows below `ransac_r` before marking a channel bad
- `ransac_t::Float64=100.0`: RANSAC inlier distance threshold (in signal units)
- `amp_t::Real=400.0`: amplitude rejection threshold; channels with `max > +amp_t` or `min < −amp_t` are flagged
- `nbad::Int64=1`: minimum number of bad channels to declare an epoch bad; must be in `[1, nchannels(obj)]`

# Detection methods

- `:flat`: channels with a large proportion of flat windows
- `:rmse`: channels whose RMSE vs the median reference is outside the 95 % CI
- `:rmsd`: same as `:rmse` using RMSD
- `:euclid`: same using Euclidean distance
- `:var`: channels with IQR-outlier variance
- `:p2p`: channels with excessive peak-to-peak variation
- `:tkeo`: channels where z-scored TKEO diverges from z-scored signal
- `:kurt`: channels with z-scored kurtosis exceeding `z`
- `:z`: channels with a large proportion of samples above the z-score threshold
- `:ransac`: channels poorly correlated with their nearest spatial neighbor
- `:amp`: channels exceeding ±`amp_t`

# Returns

- `Vector{Int64}`: sorted, unique indices of bad epochs
"""
function epoch_reject(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Union{Symbol, Vector{Symbol}} = [
        :flat,
        :rmse,
        :rmsd,
        :euclid,
        :var,
        :p2p,
        :tkeo,
        :kurt,
        :z,
        :ransac,
        :amp,
    ],
    w::Int64 = sr(obj),
    flat_tol::Float64 = 0.1,
    flat_fr::Float64 = 0.3,
    p::Float64 = 0.99,
    tc::Float64 = 0.2,
    tkeo_method::Symbol = :pow,
    z::Real = 3,
    ransac_r::Float64 = 0.8,
    ransac_tr::Float64 = 0.4,
    ransac_t::Float64 = 100.0,
    amp_t::Real = 400.0,
    nbad::Int64 = 1,
)::Vector{Int64}
    # validate
    _in(p, (0, 1), "p")
    _in(tc, (0, 1), "tc")
    _in(ransac_r, (0, 1), "ransac_r")
    _in(ransac_tr, (0, 1), "ransac_tr")
    nbad >= 1 || throw(ArgumentError("nbad must be ≥ 1."))
    nbad <= size(obj, 1) || throw(ArgumentError("nbad must be ≤ $(size(obj, 1))."))

    typeof(method) != Vector{Symbol} && (method = [method])
    for m in method
        _check_var(
            m,
            [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp],
            "method",
        )
    end

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    ch_list = labels(obj)[ch]

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # epoch length
    ep_len = epoch_len(obj)

    # validate
    for m in (:rmse, :rmsd, :euclid)
        m in method && ch_n < 2 && throw(ArgumentError(":$m requires > 1 channel."))
    end

    # pre-allocate outputs: bad-channel mask and bad-epoch accumulator
    bc = zeros(Bool, nchannels(obj))
    be = Int64[]

    # -----------------------------------------------------------------------
    if :flat in method
        w < ep_len || throw(ArgumentError("w must be < $ep_len."))
        _info("Using :flat method")

        bad_chs   = zeros(Bool, ch_n, ep_n)
        n_samples = size(obj.data, 2)

        Threads.@threads :static for linear_idx = 1:(ch_n * ep_n)
            ci = (linear_idx - 1) % ch_n + 1
            ei = (linear_idx - 1) ÷ ch_n + 1
            sm = [
                mean(@view obj.data[ch[ci], idx_w:(idx_w + w), ei])
                for idx_w = 1:w:(n_samples - w)
            ]
            r = count(abs.(diff(sm)) .< flat_tol) / length(sm)
            @inbounds bad_chs[ci, ei] = r > flat_fr
        end

        bc[ch] = bc[ch] .|| vec(any(bad_chs; dims = 2))
        append!(be, findall(vec(sum(bad_chs; dims = 1)) .>= nbad))
    end

    # -----------------------------------------------------------------------
    if :rmse in method
        _info("Using :rmse method")
        for ep_idx = 1:ep_n
            bad_chs = detect_rmse(@view obj.data[ch, :, ep_idx])
            bc[ch]  = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end
    end

    # -----------------------------------------------------------------------
    if :rmsd in method
        _info("Using :rmsd method")
        for ep_idx = 1:ep_n
            bad_chs = detect_rmsd(@view obj.data[ch, :, ep_idx])
            bc[ch]  = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end
    end

    # -----------------------------------------------------------------------
    if :euclid in method
        _info("Using :euclid method")
        # each column of bad_mat is written by exactly one thread - no race condition
        bad_mat = zeros(Bool, ch_n, ep_n)
        Threads.@threads :static for ep_idx = 1:ep_n
            @inbounds bad_mat[:, ep_idx] = detect_euclid(@view obj.data[ch, :, ep_idx])
        end
        bc[ch] = bc[ch] .|| vec(any(bad_mat; dims = 2))
        append!(be, findall(vec(sum(bad_mat; dims = 1)) .>= nbad))
    end

    # -----------------------------------------------------------------------
    if :var in method
        _info("Using :var method")
        s_v = var(obj.data[ch, :, :]; dims = 2)
        o = reshape(outlier_detect(vec(s_v); method = :iqr), ch_n, ep_n)
        bc[ch] = bc[ch] .|| vec(any(o; dims = 2))
        append!(be, findall(vec(sum(o; dims = 1)) .>= nbad))
    end

    # -----------------------------------------------------------------------
    if :p2p in method
        _info("Using :p2p method")
        # sequential: bc[ch] and be are shared - cannot parallelise safely here
        for ep_idx = 1:ep_n
            bad_chs = detect_p2p(@view(obj.data[ch, :, ep_idx]); w = w, p = p)
            bc[ch]  = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end
    end

    # -----------------------------------------------------------------------
    if :tkeo in method
        _info("Using :tkeo method")
        for ep_idx = 1:ep_n
            bad_chs = detect_tkeo(
                obj.data[ch, :, ep_idx], obj.time_pts; tkeo_method = tkeo_method, p = p,
            )
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end
    end

    # -----------------------------------------------------------------------
    if :kurt in method
        _info("Using :kurt method")
        z > 0 || throw(ArgumentError("z must be > 0."))

        k = zeros(ch_n, ep_n)
        Threads.@threads :static for linear_idx = 1:(ch_n * ep_n)
            ci = (linear_idx - 1) % ch_n + 1
            ei = (linear_idx - 1) ÷ ch_n + 1
            @inbounds k[ci, ei] = kurtosis(@view obj.data[ch[ci], :, ei])
        end
        k       = normalize_zscore(k)
        bad_idx = abs.(k) .> z
        bc[ch]  = bc[ch] .|| vec(any(bad_idx; dims = 2))
        append!(be, findall(vec(sum(bad_idx; dims = 1)) .>= nbad))
    end

    # -----------------------------------------------------------------------
    if :z in method
        _info("Using :z method")
        z > 0 || throw(ArgumentError("z must be > 0."))

        # pass 1 - global z-score threshold
        s_z     = normalize_zscore(obj.data[ch, :, :]; bych = false)
        above_z = abs.(s_z) .> z
        k       = zeros(ch_n, ep_n)
        Threads.@threads :static for linear_idx = 1:(ch_n * ep_n)
            ci = (linear_idx - 1) % ch_n + 1
            ei = (linear_idx - 1) ÷ ch_n + 1
            @inbounds k[ci, ei] = count(above_z[ci, :, ei]) / size(above_z, 2)
        end
        bad_idx = k .> p
        bc[ch]  = bc[ch] .|| vec(any(bad_idx; dims = 2))
        append!(be, findall(vec(sum(bad_idx; dims = 1)) .>= nbad))

        # pass 2 - per-channel threshold (z + 1)
        above_z2 = abs.(normalize_zscore(obj.data[ch, :, :]; bych = false)) .> (z + 1)
        k2       = zeros(ch_n, ep_n)
        Threads.@threads :static for linear_idx = 1:(ch_n * ep_n)
            ci = (linear_idx - 1) % ch_n + 1
            ei = (linear_idx - 1) ÷ ch_n + 1
            @inbounds k2[ci, ei] = count(above_z2[ci, :, ei]) / size(above_z2, 2)
        end
        bad_idx2 = k2 .> p
        bc[ch]   = bc[ch] .|| vec(any(bad_idx2; dims = 2))
        append!(be, findall(vec(sum(bad_idx2; dims = 1)) .>= nbad))
    end

    # -----------------------------------------------------------------------
    if :ransac in method
        _info("Using :ransac method")
        _check_datatype(obj, ["eeg", "seeg", "ecog", "meg"])

        chs = get_channel(obj; type = ["eeg", "seeg", "ecog", "meg", "mag", "grad"])
        length(setdiff(ch_list, chs)) == 0 ||
            throw(ArgumentError("ch must contain only signal channels."))

        chs  = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)

        for ep_idx = 1:ep_n
            bad_chs = detect_ransac(
                obj.data[ch, :, ep_idx];
                loc_x = locs[!, :loc_x], loc_y = locs[!, :loc_y],
                w = w, ransac_t = ransac_t, ransac_r = ransac_r, ransac_tr = ransac_tr,
            )
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end
    end

    # -----------------------------------------------------------------------
    if :amp in method
        _info("Using :amp method")
        for ep_idx = 1:ep_n
            bad_chs = detect_amp(@view(obj.data[ch, :, ep_idx]); amp_t = amp_t)
            bc[ch]  = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end
    end

    return sort(unique(be))
end
