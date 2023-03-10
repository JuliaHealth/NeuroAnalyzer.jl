export detect_bad
export detect_channel_flat

"""
    detect_bad(obj; method, ch_t)

Detect bad OBJ channels and epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p, :var]`: detection method:
    - `:flat`: flat channel(s)
    - `:p2p`: peak-to-peak amplitude; good for detecting transient artifacts
    - `:var`: mean signal variance outside of 95%CI and variance inter-quartile outliers
    - `:rmse`: RMSE vs average channel outside of 95%CI
    - `:rmsd`: RMSD
    - `:euclid`: Euclidean distance
- `w::Int64=10`: window width in samples (signal is averaged within `w`-width window)
- `ftol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
- `fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
- `p::Float64=0.95`: probability threshold (0.0 to 1.0) for marking channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326) 
- `tc::Float64=0.3`: threshold (0.0 to 1.0) of bad channels ratio to mark the epoch as bad

# Returns

Named tuple containing:
- `bad_m::Matrix{Bool}`: matrix of bad channels × epochs
- `bad_epochs::Vector{Int64}`: list of bad epochs
"""
function detect_bad(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj, type=Symbol(obj.header.recording[:data_type])), method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p], w::Int64=10, ftol::Float64=0.1, fr::Float64=0.3, p::Float64=0.95, tc::Float64=0.2)

    for idx in method
        _check_var(idx, [:flat, :rmse, :rmsd, :euclid, :var, :p2p], "method")
    end

    p < 0 || p > 1 && throw(ArgumentError("p must be ≥ 0.0 and ≤ 1.0"))
    tc < 0 || tc > 1 && throw(ArgumentError("t must be ≥ 0.0 and ≤ 1.0"))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)    

    signal = obj.data

    bad_m = zeros(Bool, ch_n, ep_n)
    bad_epochs = Vector{Int64}()

    if :flat in method
        @inbounds @simd for ep_idx in 1:ep_n
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                r = @views detect_channel_flat(signal[channel[ch_idx], :, ep_idx], w=w, tol=ftol)
                if r > fr
                    bad_channels_score += 1
                    bad_channels[ch_idx] = true
                end
            end
            [bad_channels[ch_idx] == true && (bad_m[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, ep_idx)
        end
    end
    
    if :rmse in method
        @inbounds @simd for ep_idx in 1:ep_n
            ch_m = @views vec(median(signal[channel, :, ep_idx], dims=1))
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)

            rmse_ch = zeros(ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                rmse_ch[ch_idx] = @views rmse(signal[channel[ch_idx], :, ep_idx], ch_m)
            end
            Threads.@threads for ch_idx in 1:ch_n
                if rmse_ch[ch_idx] < HypothesisTests.confint(OneSampleTTest(rmse_ch))[1] || rmse_ch[ch_idx] > HypothesisTests.confint(OneSampleTTest(rmse_ch))[2]
                    bad_channels_score += 1
                    bad_channels[ch_idx] = true
                end
            end
            [bad_channels[ch_idx] == true && (bad_m[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, ep_idx)
        end
    end

    if :rmsd in method
        @inbounds @simd for ep_idx in 1:ep_n
            ch_m = @views vec(median(signal[channel, :, ep_idx], dims=1))
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)

            rmsd_ch = zeros(ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                rmsd_ch[ch_idx] = @views Distances.rmsd(signal[channel[ch_idx], :, ep_idx], ch_m)
            end
            Threads.@threads for ch_idx in 1:ch_n
                if rmsd_ch[ch_idx] < HypothesisTests.confint(OneSampleTTest(rmsd_ch))[1] || rmsd_ch[ch_idx] > HypothesisTests.confint(OneSampleTTest(rmsd_ch))[2]
                    bad_channels_score += 1
                    bad_channels[ch_idx] = true
                end
            end
            [bad_channels[ch_idx] == true && (bad_m[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, ep_idx)
        end
    end

    if :euclid in method
        @inbounds @simd for ep_idx in 1:ep_n
            ch_m = @views vec(median(signal[channel, :, ep_idx], dims=1))
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)

            ed_ch = zeros(ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                ed_ch[ch_idx] = @views Distances.euclidean(signal[channel[ch_idx], :, ep_idx], ch_m)
            end
            Threads.@threads for ch_idx in 1:ch_n
                if ed_ch[ch_idx] < HypothesisTests.confint(OneSampleTTest(ed_ch))[1] || ed_ch[ch_idx] > HypothesisTests.confint(OneSampleTTest(ed_ch))[2]
                    bad_channels_score += 1
                    bad_channels[ch_idx] = true
                end
            end
            [bad_channels[ch_idx] == true && (bad_m[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, ep_idx)
        end
    end

    if :var in method
        s_v = @views var(signal[channel, :, :], dims=2)
        # mean variance
        s_mv = @views vec(mean(s_v, dims=3))
        # variance outliers
        o = reshape(outlier_detect(vec(s_v), method=:iqr), ch_n, ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)
            ch_v = @views vec(var(signal[channel, :, ep_idx], dims=2))
            s_mv = vcat(s_mv, ch_v)
            Threads.@threads for ch_idx in 1:ch_n
                #if ch_v[ch_idx] > HypothesisTests.confint(OneSampleTTest(s_mv))[2] || o[ch_idx, ep_idx]
                if o[ch_idx, ep_idx]
                    bad_channels_score += 1
                    bad_channels[ch_idx] = true
                end
            end
            [bad_channels[ch_idx] == true && (bad_m[ch_idx, ep_idx] = true) for ch_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, ep_idx)
        end
    end

    if :p2p in method
        @inbounds @simd for ep_idx in 1:ep_n
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)
            Threads.@threads for ch_idx in 1:ch_n
                s = Vector{Float64}()
                for length_idx in 1:w:(length(signal[channel[ch_idx], :, ep_idx]) - w)
                    @views push!(s, median(signal[channel[ch_idx], length_idx:(length_idx + w), ep_idx]))
                end
                p2p = @views round.(diff(s), digits=-2)
                s_m = @views mean(signal[channel[ch_idx], :, ep_idx])
                s_s = @views std(signal[channel[ch_idx], :, ep_idx])
                s_u = s_m + quantile.(Distributions.Normal(), p) * s_s
                s_l = s_m - quantile.(Distributions.Normal(), p) * s_s
                p2p_p = zeros(Bool, length(p2p))
                p2p_p[p2p .> s_u] .= true
                p2p_p[p2p .< s_l] .= true
                if sum(p2p_p) > 0
                    bad_channels_score += 1
                    bad_channels[ch_idx] = true
                end
            end
            for ch_idx in 1:ch_n
                bad_channels[ch_idx] == true && (bad_m[ch_idx, ep_idx] = true)
            end 
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, ep_idx)
        end
    end

    return (bad_m=bad_m, bad_epochs=sort(unique(bad_epochs)))
end

"""
    detect_channel_flat(signal; w, tol)

Detect flat channel.

# Arguments

- `signal::AbstractVector`
- `w::Int64=4`: window width in samples (signal is averaged within `w`-width window)
- `tol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance

# Returns

Named tuple containing:
- `r::Float64`: flat to non-flat segments ratio
"""
function detect_channel_flat(signal::AbstractVector; w::Int64=4, tol::Float64=0.1)
    w < length(signal) || throw(ArgumentError("w must be < $(length(signal))"))
    sm = Vector{Float64}()
    for idx in 1:w:(length(signal) - w)
        @views push!(sm, mean(signal[idx:(idx + w)]))
    end
    return count(abs.(diff(sm)) .< tol) / length(sm)
end

