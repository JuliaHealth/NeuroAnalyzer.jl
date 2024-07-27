export detect_bad
export detect_bad!

"""
    detect_bad(obj; <keyword arguments>)

Detect bad channels and epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]`: detection method:
    - `:flat`: flat channel(s)
    - `:rmse`: RMSE vs average channel outside of 95% CI
    - `:rmsd`: RMSD
    - `:euclid`: Euclidean distance
    - `:var`: mean signal variance outside of 95% CI and variance inter-quartile outliers
    - `:p2p`: mark bad channels based on peak-to-peak amplitude; good for detecting transient artifacts
    - `:tkeo`: mark bad channels based on z-score TKEO value outside of 95% CI
    - `:kurt`: mark bad channels based on z-scores of kurtosis values
    - `:z`: mark bad channels based on their z-score of amplitude
    - `:ransac`: calculate each channel correlation to its nearest neighbor with outliers removed using random sample consensus (RANSAC) in successive 1-s segments; channel signals exhibiting low correlation to signals in neighboring scalp channels in individual windows (here, r < `ransac_r` at more than `ransac_tr` of the data points) are marked as bad
    - `:amp`: mark bad channels based on their amplitude
- `w::Int64=10`: window width in samples (signal is averaged within `w`-width window)
- `flat_tol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
- `flat_fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
- `p::Float64=0.99`: probability threshold (0.0 to 1.0) for marking a channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326); also threshold for `:z` method: percentage of channel length per epoch for marking a channel as bad
- `tc::Float64=0.3`: threshold (0.0 to 1.0) of bad channels ratio to mark the epoch as bad
- `tkeo_method::Symbol=:pow`: method of calculating TKEO, see `tkeo()` for details
- `z::Real=3`: threshold number of z-scores 
- `ransac_r::Float64=0.8`: threshold (0.0 to 1.0) correlation between channels
- `ransac_tr::Float64=0.4`: threshold (0.0 to 1.0) ratio of uncorrelated channels
- `ransac_t::Float64=100.0`: threshold distance of a sample point to the regression hyperplane to determine if it fits the model well
- `amp_t::Float64=400.0`: two-sided rejection amplitude threshold (± 400 μV)

# Returns

Named tuple containing:
- `bm::Matrix{Bool}`: matrix of bad channels × epochs
- `be::Vector{Int64}`: list of bad epochs
"""
function detect_bad(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp], w::Int64=10, flat_tol::Float64=0.1, flat_fr::Float64=0.3, p::Float64=0.99, tc::Float64=0.2, tkeo_method::Symbol=:pow, z::Real=3, ransac_r::Float64=0.8, ransac_tr::Float64=0.4, ransac_t::Float64=100.0, amp_t::Float64=400.0)

    @assert !(p < 0 || p > 1) "p must in [0.0, 1.0]."
    @assert !(tc < 0 || tc > 1) "tc must in [0.0, 1.0]."
    @assert !(ransac_r < 0 || ransac_r > 1) "ransac_r must in [0.0, 1.0]."
    @assert !(ransac_tr < 0 || ransac_tr > 1) "ransac_tr must in [0.0, 1.0]."

    typeof(method) != Vector{Symbol} && (method = [method])
    for idx in method
        _check_var(idx, [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp], "method")
    end

    ch_list = isa(ch, String) ? [ch] : ch
    ch = get_channel(obj, ch=ch)
    ch_n = length(ch)
    ep_n = nepochs(obj)

    :rmse in method && @assert ch_n > 1 ":rmse method requires > 1 channel."
    :rmsd in method && @assert ch_n > 1 ":rmsd method requires > 1 channel."
    :euclid in method && @assert ch_n > 1 ":euclid method requires > 1 channel."

    bm = zeros(Bool, ch_n, ep_n)
    be = Vector{Int64}()

    if :flat in method
        _info("Using :flat method")
        @assert w < size(obj.data, 2) "w must be < $(size(obj.data, 2))."
        @inbounds for ep_idx in 1:ep_n
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                s = @views obj.data[ch[ch_idx], :, ep_idx]
                sm = Vector{Float64}()
                for idx in 1:w:(length(s) - w)
                    @views push!(sm, mean(s[idx:(idx + w)]))
                end
                r = count(abs.(diff(sm)) .< flat_tol) / length(sm)
                if r > flat_fr
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end
    end

    if :rmse in method
        _info("Using :rmse method")
        @inbounds for ep_idx in 1:ep_n
            ch_m = @views vec(median(obj.data[ch, :, ep_idx], dims=1))
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            rmse_ch = zeros(ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                rmse_ch[ch_idx] = @views rmse(obj.data[ch[ch_idx], :, ep_idx], ch_m)
            end
            Threads.@threads for ch_idx in 1:ch_n
                if rmse_ch[ch_idx] < HypothesisTests.confint(OneSampleTTest(rmse_ch))[1] || rmse_ch[ch_idx] > HypothesisTests.confint(OneSampleTTest(rmse_ch))[2]
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end
    end

    if :rmsd in method
        _info("Using :rmsd method")
        @inbounds for ep_idx in 1:ep_n
            ch_m = @views vec(median(obj.data[ch, :, ep_idx], dims=1))
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)

            rmsd_ch = zeros(ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                rmsd_ch[ch_idx] = @views Distances.rmsd(obj.data[ch[ch_idx], :, ep_idx], ch_m)
            end
            Threads.@threads for ch_idx in 1:ch_n
                if rmsd_ch[ch_idx] < HypothesisTests.confint(OneSampleTTest(rmsd_ch))[1] || rmsd_ch[ch_idx] > HypothesisTests.confint(OneSampleTTest(rmsd_ch))[2]
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end
    end

    if :euclid in method
        _info("Using :euclid method")
        @inbounds for ep_idx in 1:ep_n
            ch_m = @views vec(median(obj.data[ch, :, ep_idx], dims=1))
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)

            ed_ch = zeros(ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                ed_ch[ch_idx] = @views Distances.euclidean(obj.data[ch[ch_idx], :, ep_idx], ch_m)
            end
            Threads.@threads for ch_idx in 1:ch_n
                if ed_ch[ch_idx] < HypothesisTests.confint(OneSampleTTest(ed_ch))[1] || ed_ch[ch_idx] > HypothesisTests.confint(OneSampleTTest(ed_ch))[2]
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end
    end

    if :var in method
        _info("Using :var method")
        s_v = @views var(obj.data[ch, :, :], dims=2)
        # mean variance
        s_mv = @views vec(mean(s_v, dims=3))
        # variance outliers
        o = reshape(outlier_detect(vec(s_v), method=:iqr), ch_n, ep_n)

        @inbounds for ep_idx in 1:ep_n
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            ch_v = @views vec(var(obj.data[ch, :, ep_idx], dims=2))
            s_mv = vcat(s_mv, ch_v)
            Threads.@threads for ch_idx in 1:ch_n
                #if ch_v[ch_idx] > HypothesisTests.confint(OneSampleTTest(s_mv))[2] || o[ch_idx, ep_idx]
                if o[ch_idx, ep_idx]
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end
    end

    if :p2p in method
        _info("Using :p2p method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                s = Vector{Float64}()
                for length_idx in 1:w:(length(obj.data[ch[ch_idx], :, ep_idx]) - w)
                    @views push!(s, median(obj.data[ch[ch_idx], length_idx:(length_idx + w), ep_idx]))
                end
                p2p = @views round.(diff(s), digits=-2)
                s_m = @views mean(obj.data[ch[ch_idx], :, ep_idx])
                s_s = @views std(obj.data[ch[ch_idx], :, ep_idx])
                s_u = s_m + quantile.(Distributions.Normal(), p) * s_s
                s_l = s_m - quantile.(Distributions.Normal(), p) * s_s
                p2p_p = zeros(Bool, length(p2p))
                p2p_p[p2p .> s_u] .= true
                p2p_p[p2p .< s_l] .= true
                if sum(p2p_p) > 0
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            for ch_idx in 1:ch_n
                bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true)
            end
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end
    end

    if :tkeo in method
        _info("Using :tkeo method")
        @inbounds for ep_idx in 1:ep_n
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                s = @views tkeo(obj.data[ch[ch_idx], :, ep_idx], obj.epoch_time, method=tkeo_method)
                z_signal = vec(zscore(obj.data[ch[ch_idx], :, ep_idx]))
                z_tkeo = vec(zscore(s))
                # scan in 10-sample windows
                w = length(s) ÷ 10
                bad_windows = 0
                for idx in 1:w:length(s)
                    count(abs.(z_signal[idx:(idx + w - 1)]) - abs.(z_tkeo[idx:(idx + w - 1)]) .> ci2z(p)) > 1 && (bad_windows += 1)
                end
                # mark channel as bad if there is at least one bad window per epoch
                if bad_windows > 1
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
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
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                if bad_idx[ch_idx, ep_idx]
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end
    end

    if :z in method
        _info("Using :z method")
        @assert z > 0 "z must be > 0."

        # by global-channel threshold
        k = zeros(ch_n, ep_n)
        s = @views normalize_zscore(obj.data[ch, :, :], bych=false)
        s = abs.(s) .> z
        @inbounds for ep_idx in 1:ep_n
            for ch_idx in 1:ch_n
                k[ch_idx, ep_idx] = @views count(s[ch_idx, :, ep_idx]) / length(s[ch_idx, :, ep_idx])
            end
        end
        bad_idx = k .> p
        @inbounds for ep_idx in 1:ep_n
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                if bad_idx[ch_idx, ep_idx]
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end

        # by individual-channel threshold
        k = zeros(ch_n, ep_n)
        s = @views normalize_zscore(obj.data[ch, :, :], bych=false)
        s = abs.(s) .> (z + 1)
        @inbounds for ep_idx in 1:ep_n
            for ch_idx in 1:ch_n
                k[ch_idx, ep_idx] = @views count(s[ch_idx, :, ep_idx]) / length(s[ch_idx, :, ep_idx])
            end
        end
        bad_idx = k .> p
        @inbounds for ep_idx in 1:ep_n
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                if bad_idx[ch_idx, ep_idx]
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end

    end

    if :ransac in method
        _info("Using :ransac method")

        _check_datatype(obj, ["eeg", "seeg", "ecog", "meg"])
        chs = get_channel(obj, type=["eeg", "seeg", "ecog", "meg", "mag", "grad"])
        @assert length(setdiff(ch_list, chs)) == 0 "ch must contain only signal channels."
        locs = _ch_locs(obj, ch_list)
        chs = intersect(ch_list, obj.locs[!, :label])
        locs = Base.filter(:label => in(chs), obj.locs)
        @assert length(ch_list) == nrow(locs) "Some channels do not have locations."
        ch_n = length(ch)
        ep_n = nepochs(obj)
        
        # scan in 1-second windows
        w = sr(obj)

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
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                _, nearest_idx = findmin(d[ch_idx, :])
                y = @views s[ch_idx, :, ep_idx]
                x = @views s[nearest_idx, :, ep_idx]
                df = DataFrame(:y=>remove_dc(y), :x=>remove_dc(x))
                reg = createRegressionSetting(@formula(y ~ x), df)
                o = ransac(reg, t=ransac_t, k=128)["outliers"]
                idx = setdiff(1:length(x), o)
                x = x[idx]
                y = y[idx]
                c = Float64[]
                for w_idx in 1:w:(length(x) - w)
                    xx = @views x[w_idx:w_idx + w]
                    yy = @views y[w_idx:w_idx + w]
                    push!(c, cor(xx, yy))
                end
                if sum(c .< ransac_r) / length(c) > ransac_tr
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end
    end

    if :amp in method
        _info("Using :amp method")

        ch_n = length(ch)
        ep_n = nepochs(obj)
        s = @views obj.data[ch, :, :]

        @inbounds for ep_idx in 1:ep_n
            bad_chs_score = 0
            bad_chs = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                a = (amp(s[ch_idx, :, ep_idx])).p
                if a > amp_t
                    bad_chs_score += 1
                    bad_chs[ch_idx] = true
                end
            end
            [bad_chs[ch_idx] && (bm[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_chs_score / ch_n) > tc && push!(be, ep_idx)
        end
    end

    return (bm=bm, be=sort(unique(be)))

end

"""
    detect_bad!(obj; <keyword arguments>)

Detect bad channels and epochs and update the `:bad_channel` field in the OBJ header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]`: detection method:
    - `:flat`: flat channel(s)
    - `:rmse`: RMSE vs average channel outside of 95% CI
    - `:rmsd`: RMSD
    - `:euclid`: Euclidean distance
    - `:var`: mean signal variance outside of 95% CI and variance inter-quartile outliers
    - `:p2p`: mark bad channels based on peak-to-peak amplitude; good for detecting transient artifacts
    - `:tkeo`: mark bad channels based on z-score TKEO value outside of 95% CI
    - `:kurt`: mark bad channels based on z-scores of kurtosis values
    - `:z`: mark bad channels based on their z-score of amplitude
    - `:ransac`: calculate each channel correlation to its nearest neighbor with outliers removed using random sample consensus (RANSAC) in successive 1-s segments; channel signals exhibiting low correlation to signals in neighboring scalp channels in individual windows (here, r < `ransac_r` at more than `ransac_tr` of the data points) are marked as bad
    - `:amp`: mark bad channels based on their amplitude
- `w::Int64=10`: window width in samples (signal is averaged within `w`-width window)
- `flat_tol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
- `flat_fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
- `p::Float64=0.99`: probability threshold (0.0 to 1.0) for marking a channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326); also threshold for `:z` method: percentage of channel length per epoch for marking a channel as bad
- `tc::Float64=0.3`: threshold (0.0 to 1.0) of bad channels ratio to mark the epoch as bad
- `tkeo_method::Symbol=:pow`: method of calculating TKEO, see `tkeo()` for details
- `z::Real=3`: threshold number of z-scores 
- `ransac_r::Float64=0.8`: threshold (0.0 to 1.0) correlation between channels
- `ransac_tr::Float64=0.4`: threshold (0.0 to 1.0) ratio of uncorrelated channels
- `ransac_t::Float64=100.0`: threshold distance of a sample point to the regression hyperplane to determine if it fits the model well
- `amp_t::Float64=400.0`: two-sided rejection amplitude threshold (± 400 μV)
"""
function detect_bad!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp], w::Int64=10, flat_tol::Float64=0.1, flat_fr::Float64=0.3, p::Float64=0.99, tc::Float64=0.2, tkeo_method::Symbol=:pow, z::Real=3, ransac_r::Float64=0.8, ransac_tr::Float64=0.4, ransac_t::Float64=100.0, amp_t::Float64=400.0)

    bm, _ = detect_bad(obj, ch=ch, method=method, w=w, flat_tol=flat_tol, flat_fr=flat_fr, p=p, tc=tc, tkeo_method=tkeo_method, z=z, ransac_r=ransac_r, ransac_tr=ransac_tr, ransac_t=ransac_t, amp_t=amp_t)
    obj.header.recording[:bad_channel] = bm

    return nothing

end
