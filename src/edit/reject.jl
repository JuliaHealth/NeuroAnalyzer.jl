export channel_reject
export channel_reject!
export epoch_reject

function detect_flat(
    s::AbstractMatrix;
    w::Int64 = 10,
    flat_tol::Float64 = 0.1,
    flat_fr::Float64 = 0.3,
)::Vector{Bool}

    @assert w < size(s, 2) "w must be < $(size(s, 2))."

    ch_n = size(s, 1)
    bad_chs = zeros(Bool, ch_n)

    for ch_idx in 1:ch_n
        v = @views s[ch_idx, :]
        sm = Vector{Float64}()
        for idx in 1:w:(length(v) - w)
            @views push!(sm, mean(v[idx:(idx + w)]))
        end
        r = count(abs.(diff(sm)) .< flat_tol) / length(sm)
        bad_chs[ch_idx] = r > flat_fr && bad_chs[ch_idx]
    end

    return bad_chs

end

function detect_rmse(
    s::AbstractMatrix
)::Vector{Bool}

    ch_n = size(s, 1)
    bad_chs = zeros(Bool, ch_n)

    ch_m = @views vec(median(s, dims = 1))
    rmse_ch = zeros(ch_n)
    for ch_idx in 1:ch_n
        rmse_ch[ch_idx] = @views rmse(s[ch_idx, :], ch_m)
    end
    for ch_idx in 1:ch_n
        bad_chs[ch_idx] = rmse_ch[ch_idx] < HypothesisTests.confint(OneSampleTTest(rmse_ch))[1] ||
            rmse_ch[ch_idx] > HypothesisTests.confint(OneSampleTTest(rmse_ch))[2]
    end

    return bad_chs

end

function detect_rmsd(
    s::AbstractMatrix
)::Vector{Bool}

    ch_n = size(s, 1)
    bad_chs = zeros(Bool, ch_n)

    ch_m = @views vec(median(s, dims = 1))
    rmsd_ch = zeros(ch_n)
    for ch_idx in 1:ch_n
        rmsd_ch[ch_idx] = @views Distances.rmsd(s[ch_idx, :], ch_m)
    end
    for ch_idx in 1:ch_n
        bad_chs[ch_idx] = rmsd_ch[ch_idx] < HypothesisTests.confint(OneSampleTTest(rmsd_ch))[1] ||
            rmsd_ch[ch_idx] > HypothesisTests.confint(OneSampleTTest(rmsd_ch))[2]
    end

    return bad_chs

end

function detect_euclid(
    s::AbstractMatrix
)::Vector{Bool}

    ch_n = size(s, 1)
    bad_chs = zeros(Bool, ch_n)

    ch_m = @views vec(median(s, dims = 1))
    ed_ch = zeros(ch_n)
    for ch_idx in 1:ch_n
        ed_ch[ch_idx] = @views Distances.euclidean(s[ch_idx, :], ch_m)
    end
    for ch_idx in 1:ch_n
        bad_chs[ch_idx] = ed_ch[ch_idx] < HypothesisTests.confint(OneSampleTTest(ed_ch))[1] ||
            ed_ch[ch_idx] > HypothesisTests.confint(OneSampleTTest(ed_ch))[2]
    end

    return bad_chs

end

function detect_amp(
    s::AbstractMatrix;
    amp_t::Real = 400.0,
)::Vector{Bool}

    ch_n = size(s, 1)
    bad_chs = zeros(Bool, ch_n)

    for ch_idx in 1:ch_n
        bad_chs[ch_idx] = maximum(s[ch_idx, :]) > amp_t || minimum(s[ch_idx, :]) < -amp_t
    end

    return bad_chs

end

function detect_p2p(
    s::AbstractMatrix;
    w::Int64 = 10,
    p::Float64=0.95,
)::Vector{Bool}

    @assert w < size(s, 2) "w must be < $(size(s, 2))."

    ch_n = size(s, 1)
    bad_chs = zeros(Bool, ch_n)

    for ch_idx in 1:ch_n
        v = @views s[ch_idx, :]
        sm = Vector{Float64}()
        for idx in 1:w:(length(v) - w)
            @views push!(sm, mean(v[idx:(idx + w)]))
        end
        p2p = @views round.(diff(sm), digits = -2)
        s_m = @views mean(s[ch_idx, :])
        s_s = @views std(s[ch_idx, :])
        s_u = s_m + quantile.(Distributions.Normal(), p) * s_s
        s_l = s_m - quantile.(Distributions.Normal(), p) * s_s
        p2p_p = zeros(Bool, length(p2p))
        p2p_p[p2p .> s_u] .= true
        p2p_p[p2p .< s_l] .= true
        bad_chs[ch_idx] = sum(p2p_p) > 0
    end

    return bad_chs

end

function detect_tkeo(
    s::AbstractMatrix,
    t::AbstractVector;
    tkeo_method::Symbol=:pow,
    p::Float64=0.95,
)::Vector{Bool}

    ch_n = size(s, 1)
    bad_chs = zeros(Bool, ch_n)

    for ch_idx in 1:ch_n
        stkeo = @views tkeo(s[ch_idx, :], t, method = tkeo_method)
        z_signal = vec(NeuroAnalyzer.zscore(s[ch_idx, :]))
        z_tkeo = vec(NeuroAnalyzer.zscore(stkeo))
        # scan in 10-sample windows
        w = length(stkeo) ÷ 10
        bad_windows = 0
        for idx in 1:w:length(stkeo)
            count(abs.(z_signal[idx:(idx + w - 1)]) - abs.(z_tkeo[idx:(idx + w - 1)]) .> cl2z(p)) > 1 &&
                (bad_windows += 1)
        end
        # mark channel as bad if there is at least one bad window
        bad_chs[ch_idx] = bad_windows > 1
    end

    return bad_chs

end

"""
    channel_reject(obj; <keyword arguments>)

Detect bad channels.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]`: detection method:
      + `:flat`: flat channel(s)
      + `:rmse`: RMSE vs average channel outside of 95% CI
      + `:rmsd`: RMSD
      + `:euclid`: Euclidean distance
      + `:var`: mean signal variance outside of 95% CI and variance inter-quartile outliers
      + `:p2p`: mark bad channels based on peak-to-peak amplitude; good for detecting transient artifacts
      + `:tkeo`: mark bad channels based on z-score TKEO value outside of 95% CI
      + `:kurt`: mark bad channels based on z-scores of kurtosis values
      + `:z`: mark bad channels based on their z-score of amplitude
      + `:ransac`: calculate each channel correlation to its nearest neighbor with outliers removed using random sample consensus (RANSAC) in successive 1-s segments; channel signals exhibiting low correlation to signals in neighboring scalp channels in individual windows (here, r < `ransac_r` at more than `ransac_tr` of the data points) are marked as bad
      + `:amp`: mark bad channels based on their amplitude
  - `w::Int64=sr(obj)`: window width in samples (signal is averaged within `w`-width window), default is 1 second
  - `flat_tol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
  - `flat_fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
  - `p::Float64=0.95`: probability threshold (0.0 to 1.0) for marking a channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326); also threshold for `:z` method: percentage of channel length per epoch for marking a channel as bad
  - `tkeo_method::Symbol=:pow`: method of calculating TKEO, see `tkeo()` for details
  - `z::Real=3`: threshold number of z-scores
  - `ransac_r::Float64=0.8`: threshold (0.0 to 1.0) correlation between channels
  - `ransac_tr::Float64=0.4`: threshold (0.0 to 1.0) ratio of uncorrelated channels
  - `ransac_t::Float64=100.0`: threshold distance of a sample point to the regression hyperplane to determine if it fits the model well
  - `amp_t::Real=400.0`: two-sided rejection amplitude threshold (maximum > +amp_t or minimum < -amp_t)

# Returns

  - `bc::Vector{Bool}`: vector of bad channels
"""
function channel_reject(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Union{Symbol, Vector{Symbol}} = [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp],
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

    @assert !(p < 0 || p > 1) "p must in [0.0, 1.0]."
    @assert !(tc < 0 || tc > 1) "tc must in [0.0, 1.0]."
    @assert !(ransac_r < 0 || ransac_r > 1) "ransac_r must in [0.0, 1.0]."
    @assert !(ransac_tr < 0 || ransac_tr > 1) "ransac_tr must in [0.0, 1.0]."

    typeof(method) != Vector{Symbol} && (method = [method])
    for idx in method
        _check_var(idx, [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp], "method")
    end

    ch = get_channel(obj; ch = ch)
    ch_list = labels(obj)[ch]
    ch_n = length(ch)
    ep_n = nepochs(obj)

    :rmse in method && @assert ch_n > 1 ":rmse method requires > 1 channel."
    :rmsd in method && @assert ch_n > 1 ":rmsd method requires > 1 channel."
    :euclid in method && @assert ch_n > 1 ":euclid method requires > 1 channel."

    bc = zeros(Bool, ch_n)

    if :flat in method

        _info("Using :flat method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_flat(obj.data[ch, :, ep_idx], w=w, flat_tol=flat_tol, flat_fr=flat_fr)
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :rmse in method

        _info("Using :rmse method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_rmse(obj.data[ch, :, ep_idx])
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :rmsd in method

        _info("Using :rmsd method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_rmsd(obj.data[ch, :, ep_idx])
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :euclid in method

        _info("Using :euclid method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_euclid(obj.data[ch, :, ep_idx])
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :var in method

        _info("Using :var method")
        s_v = @views var(obj.data[ch, :, :], dims = 2)
        # mean variance
        s_mv = @views vec(mean(s_v, dims = 3))
        # variance outliers
        o = reshape(outlier_detect(vec(s_v); method = :iqr), ch_n, ep_n)

        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            ch_v = @views vec(var(obj.data[ch, :, ep_idx], dims = 2))
            s_mv = vcat(s_mv, ch_v)
            for ch_idx in 1:ch_n
                #if ch_v[ch_idx] > HypothesisTests.confint(OneSampleTTest(s_mv))[2] || o[ch_idx, ep_idx]
                if o[ch_idx, ep_idx]
                    bad_chs[ch_idx] = true
                end
            end
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :p2p in method

        _info("Using :p2p method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_p2p(obj.data[ch, :, ep_idx], w=w, p=p)
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :tkeo in method

        _info("Using :tkeo method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_tkeo(obj.data[ch, :, ep_idx], obj.time_pts, tkeo_method = tkeo_method, p = p)
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :kurt in method

        _info("Using :kurt method")
        @assert z > 0 "z must be > 0."
        k = zeros(ch_n, ep_n)
        @inbounds for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                k[ch_idx, ep_idx] = @views kurtosis(obj.data[ch[ch_idx], :, ep_idx])
            end
        end
        k = normalize_zscore(k)
        bad_idx = abs.(k) .> z
        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            for ch_idx in 1:ch_n
                if bad_idx[ch_idx, ep_idx]
                    bad_chs[ch_idx] = true
                end
            end
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :z in method

        _info("Using :z method")
        @assert z > 0 "z must be > 0."

        # by global-channel threshold
        k = zeros(ch_n, ep_n)
        s = @views normalize_zscore(obj.data[ch, :, :], bych = false)
        s = abs.(s) .> z
        @inbounds for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                k[ch_idx, ep_idx] = @views count(s[ch_idx, :, ep_idx]) / length(s[ch_idx, :, ep_idx])
            end
        end
        bad_idx = k .> p
        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            for ch_idx in 1:ch_n
                if bad_idx[ch_idx, ep_idx]
                    bad_chs[ch_idx] = true
                end
            end
            bc[ch] = bc[ch] .|| bad_chs
        end

        # by individual-channel threshold
        k = zeros(ch_n, ep_n)
        s = @views normalize_zscore(obj.data[ch, :, :], bych = false)
        s = abs.(s) .> (z + 1)
        @inbounds for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                k[ch_idx, ep_idx] = @views count(s[ch_idx, :, ep_idx]) / length(s[ch_idx, :, ep_idx])
            end
        end
        bad_idx = k .> p
        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            for ch_idx in 1:ch_n
                if bad_idx[ch_idx, ep_idx]
                    bad_chs[ch_idx] = true
                end
            end
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :ransac in method

        _info("Using :ransac method")

        _check_datatype(obj, ["eeg", "seeg", "ecog", "meg"])
        chs = get_channel(obj; type = ["eeg", "seeg", "ecog", "meg", "mag", "grad"])
        @assert length(setdiff(ch_list, chs)) == 0 "ch must contain only signal channels."
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)

        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]

        # Euclidean distance matrix
        d = zeros(ch_n, ch_n)
        @inbounds for idx1 in 1:ch_n
            for idx2 in 1:ch_n
                d[idx1, idx2] = euclidean([loc_x[idx1], loc_y[idx1]], [loc_x[idx2], loc_y[idx2]])
            end
        end
        # set weights not to reference to itself
        d[d .== 0] .= Inf

        s = @views obj.data[ch, :, :]

        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            for ch_idx in 1:ch_n
                _, nearest_idx = findmin(d[ch_idx, :])
                y = @views s[ch_idx, :, ep_idx]
                x = @views s[nearest_idx, :, ep_idx]
                df = DataFrame(:y=>remove_dc(y), :x=>remove_dc(x))
                reg = createRegressionSetting(@formula(y ~ x), df)
                o = ransac(reg, t = ransac_t, k = 128)["outliers"]
                idx = setdiff(1:length(x), o)
                x = x[idx]
                y = y[idx]
                c = Float64[]
                for w_idx in 1:w:(length(x) - w)
                    xx = @views x[w_idx:(w_idx + w)]
                    yy = @views y[w_idx:(w_idx + w)]
                    push!(c, cor(xx, yy))
                end
                bad_chs[ch_idx] = sum(c .< ransac_r) / length(c) > ransac_tr
            end
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    if :amp in method

        _info("Using :amp method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_amp(obj.data[ch, :, ep_idx], amp_t=amp_t)
            bc[ch] = bc[ch] .|| bad_chs
        end

    end

    return bc

end

"""
    channel_reject!(obj; <keyword arguments>)

Detect bad channels and update the `:bad_channel` field in the OBJ header.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]`: detection method:
      + `:flat`: flat channel(s)
      + `:rmse`: RMSE vs average channel outside of 95% CI
      + `:rmsd`: RMSD
      + `:euclid`: Euclidean distance
      + `:var`: mean signal variance outside of 95% CI and variance inter-quartile outliers
      + `:p2p`: mark bad channels based on peak-to-peak amplitude; good for detecting transient artifacts
      + `:tkeo`: mark bad channels based on z-score TKEO value outside of 95% CI
      + `:kurt`: mark bad channels based on z-scores of kurtosis values
      + `:z`: mark bad channels based on their z-score of amplitude
      + `:ransac`: calculate each channel correlation to its nearest neighbor with outliers removed using random sample consensus (RANSAC) in successive 1-s segments; channel signals exhibiting low correlation to signals in neighboring scalp channels in individual windows (here, r < `ransac_r` at more than `ransac_tr` of the data points) are marked as bad
      + `:amp`: mark bad channels based on their amplitude
  - `w::Int64=10`: window width in samples (signal is averaged within `w`-width window), default is 1 second
  - `flat_tol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
  - `flat_fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
  - `p::Float64=0.95`: probability threshold (0.0 to 1.0) for marking a channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326); also threshold for `:z` method: percentage of channel length per epoch for marking a channel as bad
  - `tkeo_method::Symbol=:pow`: method of calculating TKEO, see `tkeo()` for details
  - `z::Real=3`: threshold number of z-scores
  - `ransac_r::Float64=0.8`: threshold (0.0 to 1.0) correlation between channels
  - `ransac_tr::Float64=0.4`: threshold (0.0 to 1.0) ratio of uncorrelated channels
  - `ransac_t::Float64=100.0`: threshold distance of a sample point to the regression hyperplane to determine if it fits the model well
  - `amp_t::Real=400.0`: two-sided rejection amplitude threshold (maximum > +amp_t or minimum < -amp_t)

# Returns

  - `Nothing`
"""
function channel_reject!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Union{Symbol, Vector{Symbol}} = [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp],
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
    obj.header.recording[:bad_channel] = bc

    return nothing

end

"""
    epoch_reject(obj; <keyword arguments>)

Detect bad epochs.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]`: detection method:
      + `:flat`: flat channel(s)
      + `:rmse`: RMSE vs average channel outside of 95% CI
      + `:rmsd`: RMSD
      + `:euclid`: Euclidean distance
      + `:var`: mean signal variance outside of 95% CI and variance inter-quartile outliers
      + `:p2p`: mark bad epoch based on peak-to-peak amplitude; good for detecting transient artifacts
      + `:tkeo`: mark bad epoch based on z-score TKEO value outside of 95% CI
      + `:kurt`: mark bad epoch based on z-scores of kurtosis values
      + `:z`: mark bad epoch based on their z-score of amplitude
      + `:ransac`: calculate each channel correlation to its nearest neighbor with outliers removed using random sample consensus (RANSAC) in successive 1-s segments; channel signals exhibiting low correlation to signals in neighboring scalp epoch in individual windows (here, r < `ransac_r` at more than `ransac_tr` of the data points) are marked as bad
      + `:amp`: mark bad epoch based on their amplitude
  - `w::Int64=10`: window width in samples (signal is averaged within `w`-width window), default is 1 second
  - `flat_tol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
  - `flat_fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
  - `p::Float64=0.95`: probability threshold (0.0 to 1.0) for marking a channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326); also threshold for `:z` method: percentage of channel length per epoch for marking a channel as bad
  - `tkeo_method::Symbol=:pow`: method of calculating TKEO, see `tkeo()` for details
  - `z::Real=3`: threshold number of z-scores
  - `ransac_r::Float64=0.8`: threshold (0.0 to 1.0) correlation between channels
  - `ransac_tr::Float64=0.4`: threshold (0.0 to 1.0) ratio of uncorrelated channels
  - `ransac_t::Float64=100.0`: threshold distance of a sample point to the regression hyperplane to determine if it fits the model well
  - `amp_t::Real=400.0`: two-sided rejection amplitude threshold (maximum > +amp_t or minimum < -amp_t)
  - `nbad::Int64=1`: minimum number of bad channels to mark the epoch as bad

# Returns

  - `be::Vector{Int64}`: bad epochs numbers
"""
function epoch_reject(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Union{Symbol, Vector{Symbol}} = [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp],
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
    nbad::Int64=1,
)::Vector{Int64}

    @assert !(p < 0 || p > 1) "p must in [0.0, 1.0]."
    @assert !(tc < 0 || tc > 1) "tc must in [0.0, 1.0]."
    @assert !(ransac_r < 0 || ransac_r > 1) "ransac_r must in [0.0, 1.0]."
    @assert !(ransac_tr < 0 || ransac_tr > 1) "ransac_tr must in [0.0, 1.0]."
    @assert nbad >= 1 "nbad must be ≥ 1."
    @assert nbad <= size(obj, 1) "nbad must be ≤ $(size(obj, 1))."
    typeof(method) != Vector{Symbol} && (method = [method])
    for idx in method
        _check_var(idx, [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp], "method")
    end

    ch = get_channel(obj; ch = ch)
    ch_list = labels(obj)[ch]
    ch_n = length(ch)
    ep_n = nepochs(obj)

    :rmse in method && @assert ch_n > 1 ":rmse method requires > 1 channel."
    :rmsd in method && @assert ch_n > 1 ":rmsd method requires > 1 channel."
    :euclid in method && @assert ch_n > 1 ":euclid method requires > 1 channel."

    bc = zeros(Bool, ch_n)
    be = Int64[]

    if :flat in method

        _info("Using :flat method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_flat(obj.data[ch, :, ep_idx], w=w, flat_tol=flat_tol, flat_fr=flat_fr)
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :rmse in method

        _info("Using :rmse method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_rmse(obj.data[ch, :, ep_idx])
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :rmsd in method

        _info("Using :rmsd method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_rmsd(obj.data[ch, :, ep_idx])
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :euclid in method

        _info("Using :euclid method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_euclid(obj.data[ch, :, ep_idx])
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :var in method

        _info("Using :var method")
        s_v = @views var(obj.data[ch, :, :], dims = 2)
        # mean variance
        s_mv = @views vec(mean(s_v, dims = 3))
        # variance outliers
        o = reshape(outlier_detect(vec(s_v); method = :iqr), ch_n, ep_n)

        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            ch_v = @views vec(var(obj.data[ch, :, ep_idx], dims = 2))
            s_mv = vcat(s_mv, ch_v)
            for ch_idx in 1:ch_n
                #if ch_v[ch_idx] > HypothesisTests.confint(OneSampleTTest(s_mv))[2] || o[ch_idx, ep_idx]
                if o[ch_idx, ep_idx]
                    bad_chs[ch_idx] = true
                end
            end
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :p2p in method

        _info("Using :p2p method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_p2p(obj.data[ch, :, ep_idx], w=w, p=p)
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :tkeo in method

        _info("Using :tkeo method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_tkeo(obj.data[ch, :, ep_idx], obj.time_pts, tkeo_method = tkeo_method, p = p)
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :kurt in method

        _info("Using :kurt method")
        @assert z > 0 "z must be > 0."
        k = zeros(ch_n, ep_n)
        @inbounds for ep_idx in 1:ep_n
            for ch_idx in 1:ch_n
                k[ch_idx, ep_idx] = @views kurtosis(obj.data[ch[ch_idx], :, ep_idx])
            end
        end
        k = normalize_zscore(k)
        bad_idx = abs.(k) .> z
        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            for ch_idx in 1:ch_n
                if bad_idx[ch_idx, ep_idx]
                    bad_chs[ch_idx] = true
                end
            end
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :z in method

        _info("Using :z method")
        @assert z > 0 "z must be > 0."

        # by global-channel threshold
        k = zeros(ch_n, ep_n)
        s = @views normalize_zscore(obj.data[ch, :, :], bych = false)
        s = abs.(s) .> z
        @inbounds for ep_idx in 1:ep_n
            for ch_idx in 1:ch_n
                k[ch_idx, ep_idx] = @views count(s[ch_idx, :, ep_idx]) / length(s[ch_idx, :, ep_idx])
            end
        end
        bad_idx = k .> p
        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            for ch_idx in 1:ch_n
                if bad_idx[ch_idx, ep_idx]
                    bad_chs[ch_idx] = true
                end
            end
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

        # by individual-channel threshold
        k = zeros(ch_n, ep_n)
        s = @views normalize_zscore(obj.data[ch, :, :], bych = false)
        s = abs.(s) .> (z + 1)
        Threads.@threads for ep_idx in 1:ep_n
            @inbounds for ch_idx in 1:ch_n
                k[ch_idx, ep_idx] = @views count(s[ch_idx, :, ep_idx]) / length(s[ch_idx, :, ep_idx])
            end
        end
        bad_idx = k .> p
        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            for ch_idx in 1:ch_n
                if bad_idx[ch_idx, ep_idx]
                    bad_chs[ch_idx] = true
                end
            end
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :ransac in method

        _info("Using :ransac method")

        _check_datatype(obj, ["eeg", "seeg", "ecog", "meg"])
        chs = get_channel(obj; type = ["eeg", "seeg", "ecog", "meg", "mag", "grad"])
        @assert length(setdiff(ch_list, chs)) == 0 "ch must contain only signal channels."
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)

        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]

        # Euclidean distance matrix
        d = zeros(ch_n, ch_n)
        for idx1 in 1:ch_n
            for idx2 in 1:ch_n
                d[idx1, idx2] = euclidean([loc_x[idx1], loc_y[idx1]], [loc_x[idx2], loc_y[idx2]])
            end
        end
        # set weights not to reference to itself
        d[d .== 0] .= Inf

        s = @views obj.data[ch, :, :]

        @inbounds for ep_idx in 1:ep_n
            bad_chs = zeros(Bool, ch_n)
            for ch_idx in 1:ch_n
                _, nearest_idx = findmin(d[ch_idx, :])
                y = @views s[ch_idx, :, ep_idx]
                x = @views s[nearest_idx, :, ep_idx]
                df = DataFrame(:y=>remove_dc(y), :x=>remove_dc(x))
                reg = createRegressionSetting(@formula(y ~ x), df)
                o = ransac(reg, t = ransac_t, k = 128)["outliers"]
                idx = setdiff(1:length(x), o)
                x = x[idx]
                y = y[idx]
                c = Float64[]
                for w_idx in 1:w:(length(x) - w)
                    xx = @views x[w_idx:(w_idx + w)]
                    yy = @views y[w_idx:(w_idx + w)]
                    push!(c, cor(xx, yy))
                end
                if sum(c .< ransac_r) / length(c) > ransac_tr
                    bad_chs[ch_idx] = true
                end
            end
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    if :amp in method

        _info("Using :amp method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs = @views detect_amp(obj.data[ch, :, ep_idx], amp_t=amp_t)
            bc[ch] = bc[ch] .|| bad_chs
            count(bad_chs) >= nbad && push!(be, ep_idx)
        end

    end

    return sort(unique(be))

end

