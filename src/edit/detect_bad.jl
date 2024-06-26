export detect_bad

"""
    detect_bad(obj; method, ch_t)

Detect bad channels and epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z]`: detection method:
    - `:flat`: flat channel(s)
    - `:rmse`: RMSE vs average channel outside of 95% CI
    - `:rmsd`: RMSD
    - `:euclid`: Euclidean distance
    - `:var`: mean signal variance outside of 95% CI and variance inter-quartile outliers
    - `:p2p`: peak-to-peak amplitude; good for detecting transient artifacts
    - `:tkeo`: z-score TKEO value outside of 95% CI
    - `:kurt`: z-scores of kurtosis values
    - `:z`: z-score of amplitude
- `w::Int64=10`: window width in samples (signal is averaged within `w`-width window)
- `ftol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
- `fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
- `p::Float64=0.99`: probability threshold (0.0 to 1.0) for marking a channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326); also threshold for `:z` method: percentage of channel length per epoch for marking a channel as bad
- `tc::Float64=0.3`: threshold (0.0 to 1.0) of bad channels ratio to mark the epoch as bad
- `tkeo_method::Symbol=:pow`: method of calculating TKEO, see `tkeo()` for details
- `z::Real=3`: threshold number of z-scores 

# Returns

Named tuple containing:
- `bm::Matrix{Bool}`: matrix of bad channels × epochs
- `be::Vector{Int64}`: list of bad epochs
"""
function detect_bad(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z], w::Int64=10, ftol::Float64=0.1, fr::Float64=0.3, p::Float64=0.99, tc::Float64=0.2, tkeo_method::Symbol=:pow, z::Real=3)

    typeof(method) != Vector{Symbol} && (method = [method])
    for idx in method
        _check_var(idx, [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z], "method")
    end

    @assert !(p < 0 || p > 1) "p must in [0.0, 1.0]."
    @assert !(tc < 0 || tc > 1) "tc must in [0.0, 1.0]."

    _check_channels(obj, ch)
    ch_n = length(ch)
    ep_n = nepochs(obj)

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
                r = count(abs.(diff(sm)) .< ftol) / length(sm)
                if r > fr
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
    end

    return (bm=bm, be=sort(unique(be)))

end
