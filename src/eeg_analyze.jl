
"""
    mi(obj; channel)

Calculate mutual information between EEG channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

- `mi::Array{Float64, 3}`
"""
function mi(obj::NeuroAnalyzer.NEURO; channel::Union{Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    mi = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        
        # create half of the matrix
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                mi[ch_idx1, ch_idx2, ep_idx] = @views s2_mi(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx])
            end
        end

        # copy to the other half
        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                mi[ch_idx1, ch_idx2, ep_idx] = @views mi[ch_idx2, ch_idx1, ep_idx]
            end
        end
    end

    return mi
end

"""
    mi(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculate mutual information between two EEG channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

- `mi::Array{Float64, 3}`
"""
function mi(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    # check channels
    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # check epochs
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ch_n = length(channel1)
    ep_n = length(epoch1)

    mi = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_n
                mi[ch_idx1, ch_idx2, ep_idx] = @views s2_mi(obj1.signals[channel1[ch_idx1], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx2], :, epoch2[ep_idx]])
            end
        end
    end

    return mi
end

"""
    entropy(obj; channel)

Calculate entropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

Named tuple containing:
- `ent::Array{Float64, 2}`
- `s_ent::Array{Float64, 2}`: Shanon entropy
- `le_ent::Array{Float64, 2}`: log energy entropy
"""
function entropy(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    ent = zeros(ch_n, ep_n)
    sent = zeros(ch_n, ep_n)
    leent = zeros(ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ent[ch_idx, ep_idx], sent[ch_idx, ep_idx], leent[ch_idx, ep_idx] = @views s_entropy(obj.data[channel[ch_idx], :, ep_idx])
        end
    end

    return (ent=ent, s_ent=sent, le_ent=leent)
end

"""
    negentropy(obj; channel)

Calculate negentropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

- `ne::Matrix{Float64}`
"""
function negentropy(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    ne = zeros(ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ne[ch_idx, ep_idx] = @views s_negentropy(obj.data[channel[ch_idx], :, ep_idx])
        end
    end

    return ne
end

"""
    tcoherence(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculate coherence (mean over time) and MSC (magnitude-squared coherence).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `ic::Array{Float64, 3}`: imaginary part of coherence
"""
function tcoherence(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)), pad::Int64=0)

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))

    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    c = zeros(length(channel1), epoch_len(obj1), length(epoch1))
    msc = zeros(length(channel1), epoch_len(obj1), length(epoch1))
    ic = zeros(length(channel1), epoch_len(obj1), length(epoch1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            c[ch_idx, :, ep_idx], msc[ch_idx, :, ep_idx], ic[ch_idx, :, ep_idx] = @views s2_tcoherence(obj1.signals[channel1[ch_idx], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx], :, epoch2[ep_idx]], pad=pad)
        end
    end

    return (c=c, msc=msc, ic=ic)
end

"""
    freqs(obj)

Return vector of frequencies and Nyquist frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple containing:
- `hz::Vector{Float64}`
- `nyquist::Float64`
"""
function freqs(obj::NeuroAnalyzer.NEURO)
    hz, nyq = s_freqs(obj.data[1, :, 1], sr(obj))
    return (hz=hz, nyquist=nyq)
end

"""
    difference(obj; channel, n, method)

Calculate mean difference and its 95% CI between EEG channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:absdiff`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `signals_statistic::Matrix{Float64}`
- `signals_statistic_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function difference(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), n::Int64=3, method::Symbol=:absdiff)

    ep_n = epoch_n(obj)
    _check_channels(obj, channel)

    s_stat = zeros(ep_n, length(channel) * n)
    s_stat_single = zeros(ep_n)
    p = zeros(ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        s_stat[ep_idx, :], s_stat_single[ep_idx], p[ep_idx] = s2_difference(obj.data[channel, :, ep_idx], obj.data[channel, :, ep_idx], n=n, method=method)
    end

    return (s_stat=s_stat, s_stat_single=s_stat_single, p=p)
end

"""
    channel_pick(obj; pick)

Return set of channel indices corresponding with `pick` of electrodes

# Arguments

- `pick::Vector{Symbol}`: pick of electrodes; picks may be combined, e.g. `[:left, :frontal]`
    - `:list`
    - `:central` (or `:c`)
    - `:left` (or `:l`)
    - `:right` (or `:r`)
    - `:frontal` (or `:f`)
    - `:temporal` (or `:t`)
    - `:parietal` (or `:p`)
    - `:occipital` (or `:o`)

# Returns

- `channels::Vector{Int64}`: channel numbers
"""
function channel_pick(obj::NeuroAnalyzer.NEURO; pick::Union{Symbol, Vector{Symbol}})

    length(labels(obj)) == 0 && throw(ArgumentError("EEG does not contain channel labels."))

    if typeof(pick) == Vector{Symbol}
        for idx in pick
            _check_var(idx, [:list, :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "pick")
        end

        # convert picks to channel labels
        c = Vector{Char}()
        for idx in pick
            (idx === :central || idx === :c) && push!(c, 'z')
            (idx === :frontal || idx === :f) && push!(c, 'F')
            (idx === :temporal || idx === :t) && push!(c, 'T')
            (idx === :parietal || idx === :p) && push!(c, 'P')
            (idx === :occipital || idx === :o) && push!(c, 'O')
        end
        
        # check which channels are in the picks list
        labels = labels(obj)
        channels = Vector{Int64}()
        for idx1 in eachindex(labels)
            for idx2 in eachindex(c)
                in(c[idx2], labels[idx1]) && push!(channels, idx1)
            end
        end

        # check for both :l and :r
        for idx1 in eachindex(pick)
            if (pick[idx1] === :left || pick[idx1] === :l)
                for idx2 in eachindex(pick)
                    if (pick[idx2] === :right || pick[idx2] === :r)
                        return channels
                    end
                end
            end
            if (pick[idx1] === :right || pick[idx1] === :r)
                for idx2 in eachindex(pick)
                    if (pick[idx2] === :left || pick[idx2] === :l)
                        return channels
                    end
                end
            end
        end

        labels = labels(obj)
        labels = labels[channels]
        pat = nothing
        for idx in pick
            # for :right remove lefts
            (idx === :right || idx === :r) && (pat = r"[z13579]$")
            # for :left remove rights
            (idx === :left || idx === :l) && (pat = r"[z02468]$")
        end
        if typeof(pat) == Regex
            for idx in length(labels):-1:1
                typeof(match(pat, labels[idx])) == RegexMatch && deleteat!(channels, idx)
            end
        end

        return channels
    else
        _check_var(pick, [:central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "pick")

        c = Vector{Char}()
        (pick === :central || pick === :c) && (c = ['z'])
        (pick === :left || pick === :l) && (c = ['1', '3', '5', '7', '9'])
        (pick === :right || pick === :r) && (c = ['2', '4', '6', '8'])
        (pick === :frontal || pick === :f) && (c = ['F'])
        (pick === :temporal || pick === :t) && (c = ['T'])
        (pick === :parietal || pick === :p) && (c = ['P'])
        (pick === :occipital || pick === :o) && (c = ['O'])

        labels = labels(obj)
        channels = Vector{Int64}()
        for idx1 in eachindex(c)
            for idx2 in eachindex(labels)
                in(c[idx1], labels[idx2]) && push!(channels, idx2)
            end
        end

        return channels
    end
end

"""
    epoch_stats(obj)

Calculate epochs statistics.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple containing:
- `e_mean::Vector(Float64)`: mean
- `e_median::Vector(Float64)`: median
- `e_std::Vector(Float64)`: standard deviation
- `e_var::Vector(Float64)`: variance
- `e_kurt::Vector(Float64)`: kurtosis
- `e_skew::Vector(Float64)`: skewness
- `e_mean_diff::Vector(Float64)`: mean diff value
- `e_median_diff::Vector(Float64)`: median diff value
- `e_max_dif::Vector(Float64)`: max difference
- `e_dev_mean::Vector(Float64)`: deviation from channel mean
"""
function epoch_stats(obj::NeuroAnalyzer.NEURO)

    ep_n = epoch_n(obj)

    e_mean = zeros(ep_n)
    e_median = zeros(ep_n)
    e_std = zeros(ep_n)
    e_var = zeros(ep_n)
    e_kurt = zeros(ep_n)
    e_skew = zeros(ep_n)
    e_mean_diff = zeros(ep_n)
    e_median_diff = zeros(ep_n)
    e_max_dif = zeros(ep_n)
    e_dev_mean = zeros(ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        e_mean[ep_idx] = @views mean(obj.data[:, :, ep_idx])
        e_median[ep_idx] = @views median(obj.data[:, :, ep_idx])
        e_std[ep_idx] = @views std(obj.data[:, :, ep_idx])
        e_var[ep_idx] = @views var(obj.data[:, :, ep_idx])
        e_kurt[ep_idx] = @views kurtosis(obj.data[:, :, ep_idx])
        e_skew[ep_idx] = @views skewness(obj.data[:, :, ep_idx])
        e_mean_diff = @views mean(diff(obj.data[:, :, ep_idx], dims=2))
        e_median_diff = @views median(diff(obj.data[:, :, ep_idx], dims=2))
        e_max_dif = @views maximum(obj.data[:, :, ep_idx]) - minimum(obj.data[:, :, ep_idx])
        e_dev_mean = @views abs(mean(obj.data[:, :, ep_idx])) - mean(obj.data[:, :, ep_idx])
    end

    return (e_mean=e_mean, e_median=e_median, e_std=e_std, e_var=e_var, e_kurt=e_kurt, e_skew=e_skew, e_mean_diff=e_mean_diff, e_median_diff=e_median_diff, e_max_dif=e_max_dif, e_dev_mean=e_dev_mean)
end

"""
    spectrogram(obj; channel, norm, mt, st, demean)

Calculate spectrogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to add
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limits
- `frq_n::Int64=0`: number of frequencies, default is length(frq_lim[1]:frq_lim[2])
- `norm::Bool=true`: normalize powers to dB
- `demean::Bool=true`: demean signal prior to analysis
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin
- `wt<:CWT=wavelet(Morlet(π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `s_pow::Array{Float64, 3}`
- `s_frq::Vector{Float64}`
- `s_t::Vector{Float64}`
"""
function spectrogram(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), pad::Int64=0, frq_lim::Tuple{Real, Real}=(0, 0), frq_n::Int64=0, method::Symbol=:standard, norm::Bool=true, demean::Bool=true, frq::Symbol=:log, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=6, wt::T=wavelet(Morlet(π), β=2)) where {T <: CWT}

    _check_var(method, [:standard, :stft, :mt, :mw, :gh, :cwt], "method")
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    # get frequency range
    fs = sr(obj)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_n == 0 && (frq_n = length(frq_lim[1]:frq_lim[2]))

    if method === :standard
        p_tmp, s_frq, _ = @views s_spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, mt=false, st=false, demean=demean)
    elseif method === :mt
        p_tmp, s_frq, _ = @views s_spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, mt=true, st=false, demean=demean)
    elseif method === :stft
        p_tmp, s_frq, _ = @views s_spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, mt=false, st=true, demean=demean)
    elseif method === :mw
    _, p_tmp, _, s_frq = @views s_wspectrogram(obj.data[1, :, 1], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc, demean=demean)
    elseif method === :gh
        p_tmp, s_frq = @views s_ghspectrogram(obj.data[1, :, 1], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, demean=demean, gw=gw)
    elseif method === :cwt
        p_tmp, s_frq = @views s_cwtspectrogram(obj.data[1, :, 1], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm, demean=demean)
    end

    s_t = linspace(0, (epoch_len(obj) / fs), size(p_tmp, 2))
    s_pow = zeros(size(p_tmp, 1), size(p_tmp, 2), ch_n, ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if method === :standard
                s_pow[:, :, ch_idx, ep_idx], _, _ = @views s_spectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=false, st=false, demean=demean)
            elseif method === :mt
                s_pow[:, :, ch_idx, ep_idx], _, _ = @views s_spectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=true, st=false, demean=demean)
            elseif method === :stft
                s_pow[:, :, ch_idx, ep_idx], _, _ = @views s_spectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=false, st=true, demean=demean)
            elseif method === :mw
                _, s_pow[:, :, ch_idx, ep_idx], _, _ = @views s_wspectrogram(obj.data[channel[ch_idx], :, ep_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc, demean=demean)
            elseif method === :gh
                s_pow[:, :, ch_idx, ep_idx], _, _ = @views s_ghspectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, demean=demean, gw=gw)
            elseif method === :cwt
                s_pow[:, :, ch_idx, ep_idx], _ = @views s_cwtspectrogram(obj.data[channel[ch_idx], :, ep_idx], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm, demean=demean)
            end

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    s_frq = round.(s_frq, digits=2)
    s_t = round.(s_t, digits=2)
    s_t .+= obj.epoch_time[1]

    return (s_pow=s_pow, s_frq=s_frq, s_t=s_t)
end

"""
    spectrum(obj; channel, pad, h)

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `h::Bool=false`: use Hilbert transform for calculations instead of FFT
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
- `amp::Array{Float64, 3}`: amplitudes
- `pow::Array{Float64, 3}`: powers
- `pha::Array{Float64, 3}: phase angles
"""
function spectrum(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), pad::Int64=0, h::Bool=false, norm::Bool=false)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fft_size = epoch_len(obj) + pad

    s_c = zeros(ComplexF64, ch_n, fft_size, ep_n)
    s_pha = zeros(ch_n, fft_size, ep_n)
    if h == true
        s_amp = zeros(ch_n, fft_size, ep_n)
        s_pow = zeros(ch_n, fft_size, ep_n)
    else
        s_amp = zeros(ch_n, fft_size ÷ 2, ep_n)
        s_pow = zeros(ch_n, fft_size ÷ 2, ep_n)
    end        

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if h == true
                s_c[ch_idx, :, ep_idx], s_amp[ch_idx, :, ep_idx], s_pow[ch_idx, :, ep_idx], s_pha[ch_idx, :, ep_idx] = @views s_hspectrum(obj.data[channel[ch_idx], :, ep_idx], pad=pad, norm=norm)
            else
                s_c[ch_idx, :, ep_idx], s_amp[ch_idx, :, ep_idx], s_pow[ch_idx, :, ep_idx], s_pha[ch_idx, :, ep_idx] = @views s_spectrum(obj.data[channel[ch_idx], :, ep_idx], pad=pad, norm=norm)
            end
        end
    end

    return (c=s_c, amp=s_amp, pow=s_pow, pha=s_pha)
end

"""
    s2t(obj; t)

Convert time in samples to seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::Int64`: time in samples

# Returns

- `t_s::Float64`: time in seconds
"""
function s2t(obj::NeuroAnalyzer.NEURO; t::Int64)
    return round(t / sr(obj), digits=2)
end

"""
    t2s(obj; t)

Convert time in seconds to samples.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::Real`: time in seconds

# Returns

- `t_s::Int64`: time in samples
"""
function t2s(obj::NeuroAnalyzer.NEURO; t::Real)
    return floor(Int64, t * sr(obj)) + 1
end

"""
    channel_stats(obj)

Calculate channels statistics per epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple containing:
- `c_mean::Matrix(Float64)`: mean
- `c_median::Matrix(Float64)`: median
- `c_std::Matrix(Float64)`: standard deviation
- `c_var::Matrix(Float64)`: variance
- `c_kurt::Matrix(Float64)`: kurtosis
- `c_skew::Matrix(Float64)`: skewness
- `c_mean_diff::Matrix(Float64)`: mean diff value
- `c_median_diff::Matrix(Float64)`: median diff value
- `c_max_dif::Matrix(Float64)`: max difference
- `c_dev_mean::Matrix(Float64)`: deviation from channel mean
"""
function channel_stats(obj::NeuroAnalyzer.NEURO)

    ch_n = channel_n(obj)
    ep_n = epoch_n(obj)
    
    c_mean = zeros(ch_n, ep_n)
    c_median = zeros(ch_n, ep_n)
    c_std = zeros(ch_n, ep_n)
    c_var = zeros(ch_n, ep_n)
    c_kurt = zeros(ch_n, ep_n)
    c_skew = zeros(ch_n, ep_n)
    c_mean_diff = zeros(ch_n, ep_n)
    c_median_diff = zeros(ch_n, ep_n)
    c_max_dif = zeros(ch_n, ep_n)
    c_dev_mean = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            c_mean[ch_idx, ep_idx] = @views mean(obj.data[ch_idx, :, ep_idx])
            c_median[ch_idx, ep_idx] = @views median(obj.data[ch_idx, :, ep_idx])
            c_std[ch_idx, ep_idx] = @views std(obj.data[ch_idx, :, ep_idx])
            c_var[ch_idx, ep_idx] = @views var(obj.data[ch_idx, :, ep_idx])
            c_kurt[ch_idx, ep_idx] = @views kurtosis(obj.data[ch_idx, :, ep_idx])
            c_skew[ch_idx, ep_idx] = @views skewness(obj.data[ch_idx, :, ep_idx])
            c_mean_diff[ch_idx, ep_idx] = @views mean(diff(obj.data[ch_idx, :, ep_idx]))
            c_median_diff[ch_idx, ep_idx] = @views median(diff(obj.data[ch_idx, :, ep_idx]))
            c_max_dif[ch_idx, ep_idx] = @views maximum(obj.data[ch_idx, :, ep_idx]) - minimum(obj.data[ch_idx, :, ep_idx])
            c_dev_mean[ch_idx, ep_idx] = @views abs(mean(obj.data[ch_idx, :, ep_idx])) - mean(obj.data[ch_idx, :, ep_idx])
        end
    end

    return (c_mean=c_mean, c_median=c_median, c_std=c_std, c_var=c_var, c_kurt=c_kurt, c_skew=c_skew, c_mean_diff=c_mean_diff, c_median_diff=c_median_diff, c_max_dif=c_max_dif, c_dev_mean=c_dev_mean)
end

"""
    snr(obj; channel, type)

Calculate SNR.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple containing:
- `snr::Matrix(Float64)`: SNR for each channel over frequencies 1:Nyquist
- `hz::Vector(Float64)`: frequencies
"""
function snr(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), type::Symbol=:rms)

    _check_var(type, [:mean, :rms], "type")
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    ep_n == 1 && throw(ArgumentError("EEG must contain ≥ 2 epochs."))

    hz, _ = s_freqs(obj.epoch_time)
    amp = zeros(ch_n, length(hz), ep_n)
    snr = zeros(ch_n, length(hz))

    # create spectrum for each channel
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            _, amp[ch_idx, :, ep_idx], _, _ = @views s_spectrum(obj.data[channel[ch_idx], :, ep_idx])
        end
    end

    # calculate SNR for each channel spectrum
    @inbounds @simd for hz_idx in 1:length(hz)
        Threads.@threads for ch_idx in 1:ch_n
            if type === :mean
                snr[ch_idx, hz_idx] = @views s_snr(amp[ch_idx, hz_idx, :])
            else
                snr[ch_idx, hz_idx] = @views s_snr2(amp[ch_idx, hz_idx, :])
            end
        end
    end

    return (snr=snr, hz=hz)
end

"""
    standardize(obj; channel)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

- `new::NeuroAnalyzer.NEURO`: standardized EEG
- `scaler::Matrix{Float64}`: standardizing matrix
"""
function standardize(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))
    
    _check_channels(obj, channel)
    ep_n = epoch_n(obj)

    scaler = Vector{Any}()

    new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        @views push!(scaler, StatsBase.fit(ZScoreTransform, obj.data[channel, :, ep_idx], dims=2)) 
        @views new.signals[channel,:, ep_idx] = StatsBase.transform(scaler[ep_idx], obj.data[channel, :, ep_idx])
    end

    reset_components!(new)
    push!(new.header[:history], "standardize(EEG)")

    return new, scaler
end

"""
    standardize!(obj; channel)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

- `scaler::Matrix{Float64}`: standardizing matrix
"""
function standardize!(obj::NeuroAnalyzer.NEURO)

    tmp, scaler = standardize(obj, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return scaler
end

"""
    fconv(obj; channel, kernel, norm)

Perform convolution in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel for convolution
- `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

- `s_convoluted::Array{Float64, 3}`: convoluted signal
"""
function fconv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), kernel::Union{Vector{<:Real}, Vector{ComplexF64}}, norm::Bool=true, pad::Int64=0)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    s_convoluted = zeros(ch_n, epoch_len(obj), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_convoluted[ch_idx, :, ep_idx] = @views s_fconv(obj.data[channel[ch_idx], :, ep_idx], kernel=kernel, norm=norm, pad=pad)
            
            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    return s_convoluted
end

"""
    tconv(obj; channel, kernel)

Perform convolution in the time domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

# Returns

- `s_convoluted::Array{Float64, 3}`: convoluted signal
"""
function tconv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), kernel::Union{Vector{<:Real}, Vector{ComplexF64}})

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    s_convoluted = zeros(ch_n, epoch_len(obj), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_convoluted[ch_idx, :, ep_idx] = s_tconv(obj.data[channel[ch_idx], :, ep_idx], kernel=kernel)
        end
    end

    return s_convoluted
end

"""
    dft(obj; channel)

Return FFT and DFT sample frequencies for a DFT.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

Named tuple containing:
- `sfft::Array{ComplexF64, 3}`: FFT
- `sf::Vector{Float64}`: sample frequencies
"""
function dft(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), pad::Int64=0)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    sfft = zeros(ComplexF64, ch_n, epoch_len(obj), ep_n)
    sf = nothing

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            sfft[ch_idx, :, ep_idx], sf = @views s_dft(obj.data[channel[ch_idx], :, ep_idx], fs=fs, pad=pad)
        end
    end

    return (sfft=sfft, sf=sf)
end

"""
    msci95(obj; channel, n, method)

Calculate mean, standard deviation and 95% confidence interval for EEG channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal (`:normal`) method or `n`-times bootstrapping (`:boot`)

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean
- `s_s::Matrix{Float64}`: standard deviation
- `s_u::Matrix{Float64}`: upper 95% CI
- `s_l::Matrix{Float64}`: lower 95% CI
"""
function msci95(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), n::Int64=3, method::Symbol=:normal)

    _check_var(method, [:normal, :boot], "method")

    _check_channels(obj, channel)
    epoch_len = epoch_len(obj)
    ep_n = epoch_n(obj)

    s_m = zeros(ep_n, epoch_len)
    s_s = zeros(ep_n, epoch_len)
    s_u = zeros(ep_n, epoch_len)
    s_l = zeros(ep_n, epoch_len)

    Threads.@threads for ep_idx in 1:ep_n
        s_m[ep_idx, :], s_s[ep_idx, :], s_u[ep_idx, :], s_l[ep_idx, :] = @views s_msci95(obj.data[channel, :, ep_idx], n=n, method=method)
    end

    return (mean=s_m, sd=s_s, upper=s_u, lower=s_l)
end

"""
    mean(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculates mean, standard deviation and 95% confidence interval for two EEG channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2:NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean by epochs
- `s_s::Matrix{Float64}`: standard deviation by epochs
- `s_u::Matrix{Float64}`: upper 95% CI bound by epochs
- `s_l::Matrix{Float64}`: lower 95% CI bound by epochs
"""
function mean(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    # check channels
    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # check epochs
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    epoch_len = epoch_len(obj1)

    s_m = zeros(ep_n, epoch_len)
    s_s = zeros(ep_n, epoch_len)
    s_u = zeros(ep_n, epoch_len)
    s_l = zeros(ep_n, epoch_len)

    Threads.@threads for ep_idx in 1:ep_n
        s1_mean = @views mean(obj1.signals[channel1, :, epoch1[ep_idx]], dims=1)
        s2_mean = @views mean(obj2.signals[channel2, :, epoch2[ep_idx]], dims=1)
        s_m[ep_idx, :] = s1_mean - s2_mean
        s1_sd = @views std(obj1.signals[channel1, :, epoch1[ep_idx]], dims=1) / sqrt(size(obj1.signals[channel1, :, epoch1[ep_idx]], 2))
        s2_sd = @views std(obj2.signals[channel2, :, epoch2[ep_idx]], dims=1) / sqrt(size(obj2.signals[channel2, :, epoch2[ep_idx]], 2))
        s_s[ep_idx, :] = sqrt.(s1_sd.^2 .+ s2_sd.^2)
        s_u[ep_idx, :] = @. s_m[ep_idx, :] + 1.96 * s_s[ep_idx, :]
        s_l[ep_idx, :] = @. s_m[ep_idx, :] - 1.96 * s_s[ep_idx, :]
    end

    return (s_m=s_m, s_s=s_s, s_u=s_u, s_l=s_l)
end

"""
    difference(obj1, obj2; channel1, channel2, epoch1, epoch2, n, method)

Calculates mean difference and 95% confidence interval for two EEG channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2:NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `n::Int64`: number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `s_stat::Matrix{Float64}`
- `s_stat_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function difference(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)), n::Int64=3, method::Symbol=:absdiff)

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)

    s_stat = zeros(ep_n, length(channel1) * n)
    s_stat_single = zeros(ep_n)
    p = zeros(ep_n)

    Threads.@threads for ep_idx in 1:ep_n
        s_stat[ep_idx, :], s_stat_single[ep_idx], p[ep_idx] = @views s2_difference(obj1.signals[channel1, :, epoch1[ep_idx]], obj2.signals[channel2, :, epoch2[ep_idx]], n=n, method=method)
    end

    return (s_stat=s_stat, statsitic_single=s_stat_single, p=p)
end

"""
   acov(obj; channel, lag=1, demean=false, norm=false)

Calculate autocovariance.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean obj prior to analysis
- `norm::Bool`: normalize autocovariance

# Returns

Named tuple containing:
- `acov::Matrix{Float64}`
- `lags::Vector{Float64}`: lags in ms
"""
function acov(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    acov = zeros(ch_n, length(-lag:lag), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            acov[ch_idx, :, ep_idx], _ = @views s_acov(obj.data[channel[ch_idx], :, ep_idx], lag=lag, demean=demean, norm=norm)
        end
    end

    lags = 1/sr(obj) .* collect(-lag:lag) .* 1000

    return (acov=acov, acov_lags=lags)
end

"""
    tenv(obj; channel, d)

Calculate temporal envelope (amplitude).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env::Array{Float64, 3}`: temporal envelope
- `s_t::Vector{Float64}`: signal time
"""
function tenv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), d::Int64=32)
    
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    t_env = zeros(ch_n, epoch_len(obj), ep_n)
    s_t = obj.epoch_time

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s = @view obj.data[channel[ch_idx], :, ep_idx]
            # find peaks
            p_idx = s_findpeaks(s, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(s))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_t[p_idx], s[p_idx])
                try
                    t_env[ch_idx, :, ep_idx] = model(s_t)
                catch
                    @error "CubicSpline error, using Loess."
                    model = loess(s_t[p_idx], s[p_idx], span=0.5)
                    t_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
                end
            else
                _info("Less than 5 peaks detected, using Loess.")
                model = loess(s_t[p_idx], s[p_idx], span=0.5)
                t_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
            end
            t_env[ch_idx, 1, ep_idx] = t_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (t_env=t_env, s_t=s_t)
end

"""
    tenv_mean(obj; channel, dims, d)

Calculate temporal envelope (amplitude): mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function tenv_mean(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = tenv(obj, channel=channel, d=d)
    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # mean over channels

        t_env_m = zeros(length(s_t), ep_n)
        t_env_u = zeros(length(s_t), ep_n)
        t_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            t_env_m[:, ep_idx] = mean(s_a[:, :, ep_idx], dims=1)
            s = std(t_env_m[:, ep_idx]) / sqrt(length(t_env_m[:, ep_idx]))
            t_env_u[:, ep_idx] = @. t_env_m[:, ep_idx] + 1.96 * s
            t_env_l[:, ep_idx] = @. t_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        t_env_m = zeros(length(s_t), ch_n)
        t_env_u = zeros(length(s_t), ch_n)
        t_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            t_env_m[:, ch_idx] = mean(s_a[ch_idx, :, :], dims=2)
            s = std(t_env_m[:, ch_idx]) / sqrt(length(t_env_m[:, ch_idx]))
            t_env_u[:, ch_idx] = @. t_env_m[:, ch_idx] + 1.96 * s
            t_env_l[:, ch_idx] = @. t_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        t_env_m, t_env_u, t_env_l, _ = tenv_mean(obj, dims=1, d=d)

        t_env_m = mean(t_env_m, dims=2)
        t_env_u = mean(t_env_u, dims=2)
        t_env_l = mean(t_env_l, dims=2)

        t_env_m = reshape(t_env_m, size(t_env_m, 1))
        t_env_u = reshape(t_env_u, size(t_env_u, 1))
        t_env_l = reshape(t_env_l, size(t_env_l, 1))
    end

    return (t_env_m=t_env_m, t_env_u=t_env_u, t_env_l=t_env_l, s_t=s_t)
end

"""
    tenv_median(obj; channel, dims, d)

Calculate temporal envelope (amplitude): median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function tenv_median(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = tenv(obj, channel=channel, d=d)
    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # median over channels

        t_env_m = zeros(length(s_t), ep_n)
        t_env_u = zeros(length(s_t), ep_n)
        t_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            t_env_m[:, ep_idx] = median(s_a[:, :, ep_idx], dims=1)
            t_idx = s_findpeaks(t_env_m[:, ep_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(t_env_m[:, ep_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], t_env_m[t_idx])
                try
                    t_env_m[:, ep_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(t_env_m[:, ep_idx]) / sqrt(length(t_env_m[:, ep_idx]))
            t_env_u[:, ep_idx] = @. t_env_m[:, ep_idx] + 1.96 * s
            t_env_l[:, ep_idx] = @. t_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        t_env_m = zeros(length(s_t), ch_n)
        t_env_u = zeros(length(s_t), ch_n)
        t_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            t_env_m[:, idx] = median(s_a[ch_idx, :, :], dims=2)
            t_idx = s_findpeaks(t_env_m[:, ch_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(t_env_m[:, ch_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], t_env_m[t_idx])
                try
                    t_env_m[:, ch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(t_env_m[:, ch_idx]) / sqrt(length(t_env_m[:, ch_idx]))
            t_env_u[:, ch_idx] = @. t_env_m[:, ch_idx] + 1.96 * s
            t_env_l[:, ch_idx] = @. t_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

        t_env_m, t_env_u, t_env_l, _ = tenv_median(obj, dims=1, d=d)

        t_env_m = median(t_env_m, dims=2)
        t_env_u = median(t_env_u, dims=2)
        t_env_l = median(t_env_l, dims=2)

        t_env_m = reshape(t_env_m, size(t_env_m, 1))
        t_env_u = reshape(t_env_u, size(t_env_u, 1))
        t_env_l = reshape(t_env_l, size(t_env_l, 1))
    end

    return (t_env_m=t_env_m, t_env_u=t_env_u, t_env_l=t_env_l, s_t=s_t)
end

"""
    penv(obj; channel, d)

Calculate power (in dB) envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `p_env::Array{Float64, 3}`: power spectrum envelope
- `p_env_frq::Vector{Float64}`: frequencies for each envelope
"""
function penv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), d::Int64=8, mt::Bool=false, nt::Int64=8)
    
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)

    psd_tmp, frq = s_psd(obj.data[1, :, 1], fs=fs, mt=mt, nt=nt)
    p_env = zeros(ch_n, length(psd_tmp), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            psd_pow, _ = s_psd(obj.data[channel[ch_idx], :, ep_idx], fs=fs, mt=mt, norm=true, nt=nt)
            # find peaks
            p_idx = s_findpeaks(psd_pow, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(psd_pow))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(frq[p_idx], psd_pow[p_idx])
                try
                    p_env[ch_idx, :, ep_idx] = model(frq)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            else
                p_env[ch_idx, :, ep_idx] = psd_pow
            end
            p_env[ch_idx, 1, ep_idx] = p_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (p_env=p_env, p_env_frq=frq)
end

"""
    penv_mean(obj; channel, dims, d)

Calculate power (in dB) envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: mean
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_mean(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=8, mt::Bool=false)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_f = psd(obj, channel=channel, norm=true, mt=mt)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # mean over channels

        p_env_m = zeros(length(s_f), ep_n)
        p_env_u = zeros(length(s_f), ep_n)
        p_env_l = zeros(length(s_f), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            p_env_m[:, ep_idx] = mean(s_p[:, :, ep_idx], dims=1)
            # find peaks
            p_idx = s_findpeaks(p_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ep_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ep_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(p_env_m[:, ep_idx]) / sqrt(length(p_env_m[:, ep_idx]))
            p_env_u[:, ep_idx] = @. p_env_m[:, ep_idx] + 1.96 * s
            p_env_l[:, ep_idx] = @. p_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        p_env_m = zeros(length(s_f), ch_n)
        p_env_u = zeros(length(s_f), ch_n)
        p_env_l = zeros(length(s_f), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            p_env_m[:, ch_idx] = mean(s_p[ch_idx, :, :], dims=2)
            # find peaks
            p_idx = s_findpeaks(p_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ch_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(p_env_m[:, ch_idx]) / sqrt(length(p_env_m[:, ch_idx]))
            p_env_u[:, ch_idx] = @. p_env_m[:, ch_idx] + 1.96 * s
            p_env_l[:, ch_idx] = @. p_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        p_env_m, p_env_u, p_env_l, _ = penv_mean(obj, dims=1, d=d)
        p_env_m = mean(p_env_m, dims=2)
        p_env_u = mean(p_env_u, dims=2)
        p_env_l = mean(p_env_l, dims=2)
        p_env_m = reshape(p_env_m, size(p_env_m, 1))
        p_env_u = reshape(p_env_u, size(p_env_u, 1))
        p_env_l = reshape(p_env_l, size(p_env_l, 1))
    end
    
    return (p_env_m=p_env_m, p_env_u=p_env_u, p_env_l=p_env_l, p_env_frq=s_f)
end

"""
    penv_median(obj; channel, dims, d)

Calculate power (in dB) envelope: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: median
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_median(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=8, mt::Bool=false)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_f = psd(obj, channel=channel, norm=true, mt=mt)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # median over channels

        p_env_m = zeros(length(s_f), ep_n)
        p_env_u = zeros(length(s_f), ep_n)
        p_env_l = zeros(length(s_f), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            p_env_m[:, ep_idx] = median(s_p[:, :, ep_idx], dims=1)
            # find peaks
            p_idx = s_findpeaks(p_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ep_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ep_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(p_env_m[:, ep_idx]) / sqrt(length(p_env_m[:, ep_idx]))
            p_env_u[:, ep_idx] = @. p_env_m[:, ep_idx] + 1.96 * s
            p_env_l[:, ep_idx] = @. p_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        p_env_m = zeros(length(s_f), ch_n)
        p_env_u = zeros(length(s_f), ch_n)
        p_env_l = zeros(length(s_f), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            p_env_m[:, ch_idx] = median(s_p[ch_idx, :, :], dims=2)
            # find peaks
            p_idx = s_findpeaks(p_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ch_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(p_env_m[:, ch_idx]) / sqrt(length(p_env_m[:, ch_idx]))
            p_env_u[:, ch_idx] = @. p_env_m[:, ch_idx] + 1.96 * s
            p_env_l[:, ch_idx] = @. p_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # median over channels and epochs
        
        p_env_m, p_env_u, p_env_l, _ = penv_median(obj, dims=1, d=d)
        p_env_m = median(p_env_m, dims=2)
        p_env_u = median(p_env_u, dims=2)
        p_env_l = median(p_env_l, dims=2)
        p_env_m = reshape(p_env_m, size(p_env_m, 1))
        p_env_u = reshape(p_env_u, size(p_env_u, 1))
        p_env_l = reshape(p_env_l, size(p_env_l, 1))
    end
    
    return (p_env_m=p_env_m, p_env_u=p_env_u, p_env_l=p_env_l, p_env_frq=s_f)
end

"""
    senv(obj; channel, d, mt, t)

Calculate spectral envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env::Array{Float64, 3}`: spectral envelope
- `s_env_t::Vector{Float64}`: spectrogram time
"""
function senv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)
    
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    s_tmp = @view obj.data[1, :, 1]
    interval = fs
    overlap = round(Int64, fs * 0.75)
    # for short signals always use multi-taper
    length(s_tmp) < 4 * fs && (mt = true)
    if mt == true
        spec_tmp = mt_spectrogram(s_tmp, fs=fs)
    else
        spec_tmp = spectrogram(s_tmp, interval, overlap, nfft=length(s_tmp), fs=fs, window=hanning)
    end
    sp_t = collect(spec_tmp.time)
    sp_t .+= obj.epoch_time[1]

    s_env = zeros(ch_n, length(sp_t), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            # prepare spectrogram
            if mt == true
                spec = @views mt_spectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs)
            else
                spec = @views spectrogram(obj.data[channel[ch_idx], :, ep_idx], interval, overlap, nfft=length(s_tmp), fs=fs, window=hanning)
            end
            s_frq = Vector(spec.freq)
            s_p = pow2db.(spec.power)

            # maximize all powers above threshold (t)
            if t !== nothing
                s_p[s_p .> t] .= 0
                reverse!(s_p)
                reverse!(s_frq)
            end
            
            f_idx = zeros(length(spec.time))
            m = maximum(s_p, dims=1)
            for idx2 in eachindex(m)
                f_idx[idx2] = s_frq[vsearch(m[idx2], s_p[:, idx2])]
            end
            # find peaks
            p_idx = s_findpeaks(f_idx, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(spec.time))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(sp_t[p_idx], f_idx[p_idx])
                try
                    s_env[ch_idx, :, ep_idx] = model(sp_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            else
                s_env[ch_idx, :, ep_idx] = f_idx
            end
            s_env[ch_idx, 1, ep_idx] = s_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (s_env=s_env, senv_t=sp_t)
end

"""
    senv_mean(obj; channel, dims, d, mt, t)

Calculate spectral envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env_m::Array{Float64, 3}`: spectral envelope: mean
- `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
- `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
- `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)
"""
function senv_mean(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)

    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_t = senv(obj, channel=channel, d=d, mt=mt, t=t)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # mean over channels

        s_env_m = zeros(length(s_t), ep_n)
        s_env_u = zeros(length(s_t), ep_n)
        s_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            s_env_m[:, ep_idx] = mean(s_p[:, :, ep_idx], dims=1)
            # find peaks
            s_idx = s_findpeaks(s_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # interpolate peaks using cubic spline or loess
            push!(s_idx, length(s_env_m[:, ep_idx]))
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ep_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(s_env_m[:, ep_idx]) / sqrt(length(s_env_m[:, ep_idx]))
            s_env_u[:, ep_idx] = @. s_env_m[:, ep_idx] + 1.96 * s
            s_env_l[:, ep_idx] = @. s_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        s_env_m = zeros(length(s_t), ch_n)
        s_env_u = zeros(length(s_t), ch_n)
        s_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            s_env_m[:, ch_idx] = mean(s_p[ch_idx, :, :], dims=2)
            # find peaks
            s_idx = s_findpeaks(s_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(s_env_m[:, ch_idx]) / sqrt(length(s_env_m[:, ch_idx]))
            s_env_u[:, ch_idx] = @. s_env_m[:, ch_idx] + 1.96 * s
            s_env_l[:, ch_idx] = @. s_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        s_env_m, s_env_u, s_env_l, _ = senv_mean(obj, dims=1, d=d, mt=mt)
        s_env_m = mean(s_env_m, dims=2)
        s_env_u = mean(s_env_u, dims=2)
        s_env_l = mean(s_env_l, dims=2)
        s_env_m = reshape(s_env_m, size(s_env_m, 1))
        s_env_u = reshape(s_env_u, size(s_env_u, 1))
        s_env_l = reshape(s_env_l, size(s_env_l, 1))
    end
    
    return (s_env_m=s_env_m, s_env_u=s_env_u, s_env_l=s_env_l, s_env_t=s_t)
end

"""
    senv_median(obj; channel, dims, d, mt)

Calculate spectral envelope: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env_m::Array{Float64, 3}`: spectral envelope: median
- `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
- `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
- `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)
"""
function senv_median(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_t = senv(obj, channel=channel, d=d, mt=mt, t=t)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # median over channels

        s_env_m = zeros(length(s_t), ep_n)
        s_env_u = zeros(length(s_t), ep_n)
        s_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            s_env_m[:, ep_idx] = median(s_p[:, :, ep_idx], dims=1)
            # find peaks
            s_idx = s_findpeaks(s_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, ep_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ep_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(s_env_m[:, ep_idx]) / sqrt(length(s_env_m[:, ep_idx]))
            s_env_u[:, ep_idx] = @. s_env_m[:, ep_idx] + 1.96 * s
            s_env_l[:, ep_idx] = @. s_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        s_env_m = zeros(length(s_t), ch_n)
        s_env_u = zeros(length(s_t), ch_n)
        s_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            s_env_m[:, ch_idx] = median(s_p[ch_idx, :, :], dims=2)
            # find peaks
            s_idx = s_findpeaks(s_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(s_env_m[:, ch_idx]) / sqrt(length(s_env_m[:, ch_idx]))
            s_env_u[:, ch_idx] = @. s_env_m[:, ch_idx] + 1.96 * s
            s_env_l[:, ch_idx] = @. s_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

        s_env_m, s_env_u, s_env_l, _ = senv_median(obj, dims=1, d=d, mt=mt)
        s_env_m = median(s_env_m, dims=2)
        s_env_u = median(s_env_u, dims=2)
        s_env_l = median(s_env_l, dims=2)
        s_env_m = reshape(s_env_m, size(s_env_m, 1))
        s_env_u = reshape(s_env_u, size(s_env_u, 1))
        s_env_l = reshape(s_env_l, size(s_env_l, 1))
    end
    
    return (s_env_m=s_env_m, s_env_u=s_env_u, s_env_l=s_env_l, s_env_t=s_t)
end

"""
    ispc(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculate ISPC (Inter-Site-Phase Clustering).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `ispc::Array{Float64, 2}`: ISPC value
- `ispc_angle::Array{Float64, 2}`: ISPC angle
- `signal_diff::Array{Float64, 3}`: signal difference (signal2 - signal1)
- `phase_diff::Array{Float64, 3}`: phase difference (signal2 - signal1)
- `s1_phase::Array{Float64, 3}`: signal 1 phase
- `s2_phase::Array{Float64, 3}`: signal 2 phase
"""
function ispc(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    ispc = zeros(ch_n, ep_n)
    ispc_angle = zeros(ch_n, ep_n)
    signal_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    phase_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    s1_phase = zeros(ch_n, epoch_len(obj1), ep_n)
    s2_phase = zeros(ch_n, epoch_len(obj1), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ispc[ch_idx, ep_idx], ispc_angle[ch_idx, ep_idx], signal_diff[ch_idx, :, ep_idx], phase_diff[ch_idx, :, ep_idx], s1_phase[ch_idx, :, ep_idx], s2_phase[ch_idx, :, ep_idx] = @views s2_ispc(obj1.signals[channel1[ch_idx], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx], :, epoch2[ep_idx]])
        end
    end

    return (ispc=ispc, ispc_angle=ispc_angle, signal_diff=signal_diff, phase_diff=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    itpc(obj; channel, t, w)

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs/trials.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc::Vector{Float64}`: ITPC or wITPC value
- `itpcz::Vector{Float64}`: Rayleigh's ITPC Z value
- `itpc_angle::Vector{Float64}`: ITPC angle
- `itpc_phases::Array{Float64, 2}`: phase difference (channel2 - channel1)
"""
function itpc(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), t::Int64, w::Union{Vector{<:Real}, Nothing}=nothing)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)
    t < 1 && throw(ArgumentError("t must be ≥ 1."))
    t > epoch_len(obj) && throw(ArgumentError("t must be ≤ $(epoch_len(obj))."))
    ep_n < 2 && throw(ArgumentError("EEG must contain ≥ 2 epochs."))

    itpc = zeros(ch_n)
    itpcz = zeros(ch_n)
    itpc_angle = zeros(ch_n)
    itpc_phases = zeros(ch_n, ep_n)

    Threads.@threads for ch_idx in 1:ch_n
        @inbounds itpc[ch_idx], itpcz[ch_idx], itpc_angle[ch_idx], itpc_phases[ch_idx, :] = @views s_itpc(reshape(obj.data[channel[ch_idx], :, :], 1, :, ep_n), t=t, w=w)
    end
    return (itpc=itpc, itpcz=itpcz, itpc_angle=itpc_angle, itpc_phases=itpc_phases)
end

"""
    pli(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculate PLI (Phase Lag Index) between `obj1` and `obj2`.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `pli::Array{Float64, 2}`: PLI value
- `signal_diff::Array{Float64, 3}`: signal difference (signal2 - signal1)
- `phase_diff::Array{Float64, 3}`: phase difference (signal2 - signal1)
- `s1_phase::Array{Float64, 3}`: signal 1 phase
- `s2_phase::Array{Float64, 3}`: signal 2 phase
"""
function pli(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    pli = zeros(ch_n, ep_n)
    signal_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    phase_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    s1_phase = zeros(ch_n, epoch_len(obj1), ep_n)
    s2_phase = zeros(ch_n, epoch_len(obj1), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pli[ch_idx, ep_idx], signal_diff[ch_idx, :, ep_idx], phase_diff[ch_idx, :, ep_idx], s1_phase[ch_idx, :, ep_idx], s2_phase[ch_idx, :, ep_idx] = @views s2_pli(obj1.signals[channel1[ch_idx], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx], :, epoch2[ep_idx]])
        end
    end

    return (pli=pli, signal_diff=signal_diff, phase_dif=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    pli(obj; channel)

Calculate PLIs (Phase Lag Index).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

- `pli_m::Array{Float64, 3}`: PLI value matrices over epochs
"""
function pli(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    pli_m = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                pli_m[ch_idx1, ch_idx2, ep_idx], _, _, _, _ = @views s2_pli(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx])
            end
        end
        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                pli_m[ch_idx1, ch_idx2, ep_idx] = @views pli_m[ch_idx2, ch_idx1, ep_idx]
            end
        end
    end

    return pli_m
end

"""
    ispc(obj; channel)

Calculate ISPCs (Inter-Site-Phase Clustering).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

- `ispc_m::Array{Float64, 3}`: ISPC value matrices over epochs
"""
function ispc(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    ispc_m = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                ispc_m[ch_idx1, ch_idx2, ep_idx], _, _, _, _, _ = @views s2_ispc(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx])
            end
        end

        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                ispc_m[ch_idx1, ch_idx2, ep_idx] = @views ispc_m[ch_idx2, ch_idx1, ep_idx]
            end
        end
    end

    return ispc_m
end

"""
    ec(obj1, obj2; type, channel1, channel2, epoch1, epoch2)

Calculate envelope correlation.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `type::Symbol=:amp`: envelope type:
    - `:amp`: amplitude
    - `:pow`: power
    - `:spec`: spectrogram
    - `:hamp`: Hilbert spectrum amplitude
- `channel1::Int64`
- `channel2::Int64`
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `ec::Vector{Float64}`: power correlation value
- `ec_p::Vector{Float64}`: power correlation p-value
"""
function ec(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; type::Symbol=:amp, channel1::Int64, channel2::Int64, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_var(type, [:amp, :pow, :spec, :hamp], "type")

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))

    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    
    ec_r = zeros(ep_n)
    ec_p = zeros(ep_n)

    # calculate envelopes
    if type === :amp
        s1, _ = tenv(obj1, channel=channel1)
        s2, _ = tenv(obj2, channel=channel2)
    elseif type === :pow
        s1, _ = penv(obj1, channel=channel1)
        s2, _ = penv(obj2, channel=channel2)
    elseif type === :spec
        s1, _ = senv(obj1, channel=channel1)
        s2, _ = senv(obj2, channel=channel2)
    elseif type === :hamp
        s1, _ = henv(obj1, channel=channel1)
        s2, _ = henv(obj2, channel=channel2)
    end
    s1 = s1[:, :, epoch1]
    s2 = s2[:, :, epoch2]
    
    # compare envelopes per epochs
    Threads.@threads for ep_idx in 1:ep_n
        ec = CorrelationTest(vec(s1[:, :, ep_idx]), vec(s2[:, :, ep_idx]))
        @inbounds ec_r[ep_idx] = ec.r
        @inbounds ec_p[ep_idx] = pvalue(ec)
    end

    return (ec=ec_r, ec_p=ec_p)
end

"""
    ged(obj1, obj2; channel1, channel2, epoch1, epoch2)

Perform generalized eigendecomposition.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: signal data to be analyzed
- `obj2::NeuroAnalyzer.NEURO`: original signal data
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

- `sged::Array{Float64, 3}`
- `ress::Matrix{Float64}`
- `resnormalized::Matrix{Float64}`
"""
function ged(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    sged = zeros(ch_n, epoch_len(obj1), ep_n)
    ress = zeros(ch_n, ep_n)
    resnormalized = zeros(ch_n, ep_n)

    Threads.@threads for ep_idx in 1:ep_n
        sged[:, :, ep_idx], ress[:, ep_idx], resnormalized[:, ep_idx] = @views s2_ged(obj1.signals[channel1, :, epoch1[ep_idx]], obj2.signals[channel2, :, epoch2[ep_idx]])
    end

    return (sged=sged, ress=ress, resnormalized=resnormalized)
end

"""
    frqinst(obj; channel)

Calculate instantaneous frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

- `frqinst::Array{Float64, 3}`
"""
function frqinst(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    frqinst = zeros(ch_n, epoch_len(obj), ep_n)
    fs = sr(obj)

    _info("frqinst() uses Hilbert transform, the signal should be narrowband for best results.")

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            frqinst[ch_idx, :, ep_idx] = @views s_frqinst(obj.data[channel[ch_idx], :, ep_idx], fs=fs)
        end
        # update progress bar
        progress_bar == true && next!(p)
    end
    return frqinst
end

"""
    itpc_s(obj; <keyword arguments>)

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Int64`
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_s::Array{Float64, 3}`: spectrogram of ITPC values
- `itpc_z_s::Array{Float64, 3}`: spectrogram ITPCz values
- `itpc_frq::Vector{Float64}`: frequencies list
"""
function itpc_s(obj::NeuroAnalyzer.NEURO; channel::Int64, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:log, w::Union{Vector{<:Real}, Nothing}=nothing)

    _check_var(frq, [:log, :lin], "frq")
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > sr(obj) ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(sr(obj) ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n must be ≥ 2."))
    if frq === :log
        frq_lim[1] == 0 && (frq_lim = (0.01, frq_lim[2]))
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    _check_channels(obj, channel)
    ep_n = epoch_n(obj)
    epoch_len = epoch_len(obj)
    ep_n < 2 && throw(ArgumentError("OBJ must contain ≥ 2 epochs."))

    itpc_s = zeros(frq_n, epoch_len)
    itpc_z_s = zeros(frq_n, epoch_len)

    # initialize progress bar
    progress_bar == true && (p = Progress(frq_n, 1))

    Threads.@threads for frq_idx in 1:frq_n
        # create Morlet wavelet
        kernel = generate_morlet(sr(obj), frq_list[frq_idx], 1, ncyc=10)
        half_kernel = floor(Int64, length(kernel) / 2) + 1
        s_conv = zeros(Float64, 1, epoch_len, ep_n)
        # convolute with Morlet wavelet
        @inbounds @simd for ep_idx in 1:ep_n
            s_conv[1, :, ep_idx] = @views DSP.conv(obj.data[channel, :, ep_idx], kernel)[(half_kernel - 1):(end - half_kernel)]
        end
        # calculate ITPC of the convoluted signals
        @inbounds @simd for t_idx in 1:epoch_len
            itpc, itpc_z, _, _ = s_itpc(s_conv, t=t_idx, w=w)
            itpc_s[frq_idx, t_idx] = itpc
            itpc_z_s[frq_idx, t_idx] = itpc_z
        end

        # update progress bar
        progress_bar == true && next!(p)
    end

    return (itpc_s=itpc_s, itpc_z_s=itpc_z_s, itpc_f=frq_list)
end

"""
    tkeo(obj; channel)

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

- `tkeo::Array{Float64, 3}`
"""
function tkeo(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    tkeo = zeros(ch_n, epoch_len(obj), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            tkeo[ch_idx, :, ep_idx] = @views s_tkeo(obj.data[channel[ch_idx], :, ep_idx])
        end
    end

    return tkeo
end

"""
    mwpsd(obj; channel, pad, norm, frq_lim, frq_n, frq, ncyc)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool`=true: normalize powers to dB
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency bounds for the spectrogram
- `frq_n::Int64=10`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin

# Returns

Named tuple containing:
- `w_pow::Array{Float64, 4}`
- `w_frq::Matrix{Float64}`
"""
function mwpsd(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), frq_n::Int64=0, frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    p_tmp, w_frq = @views s_mwpsd(obj.data[1, :, 1], fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    w_pow = zeros(ch_n, length(p_tmp), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            w_pow[ch_idx, :, ep_idx], _ = @views s_mwpsd(obj.data[channel[ch_idx], :, ep_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    return (w_pow=w_pow, w_frq=w_frq)
end

"""
    fcoherence(obj1, obj2; channel1, channel2, epoch1, epoch2, frq_lim)

Calculate coherence (mean over frequencies) and MSC (magnitude-squared coherence).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `f::Vector{Float64}`: frequencies
"""
function fcoherence(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=0, channel2::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    sr(obj1) == sr(obj2) || throw(ArgumentError("EEG1 and EEG2 must have the same sampling rate."))

    c_tmp, _, f = @views s2_fcoherence(obj1.signals[1, :, 1], obj1.signals[1, :, 1], fs=sr(obj1), frq_lim=frq_lim)
    c = zeros(length(channel1), length(c_tmp), length(epoch1))
    msc = zeros(length(channel1), length(c_tmp), length(epoch1))
    f = zeros(length(channel1), length(c_tmp), length(epoch1))
    @inbounds @simd for ep_idx in eachindex(epoch1)
        Threads.@threads for ch_idx in eachindex(channel1)
            c[ch_idx, :, ep_idx], msc[ch_idx, :, ep_idx], _ = @views s2_fcoherence(obj1.signals[channel1[ch_idx], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx], :, epoch2[ep_idx]], fs=sr(obj1), frq_lim=frq_lim)
        end
    end

    return (c=c, msc=msc, f=f)
end

"""
    vartest(obj; channel)

Calculate variance F-test.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function vartest(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    f = zeros(ch_n, ch_n, ep_n)
    p = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
       Threads.@threads for ch_idx1 in 1:ch_n
            # create half of the matrix
            for ch_idx2 in 1:ch_idx1
                ftest = @views VarianceFTest(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx])
                f[ch_idx1, ch_idx2, ep_idx] = ftest.F
                p[ch_idx1, ch_idx2, ep_idx] = pvalue(ftest)
            end
        end
        # copy to the other half
        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                f[ch_idx1, ch_idx2, ep_idx] = @views f[ch_idx2, ch_idx1, ep_idx]
                p[ch_idx1, ch_idx2, ep_idx] = @views p[ch_idx2, ch_idx1, ep_idx]
            end
        end
    end

    return (f=f, p=p)
end

"""
    vartest(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculate variance F-test.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function vartest(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    f = zeros(ch_n, ch_n, ep_n)
    p = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
       Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_n
                ftest = @views VarianceFTest(obj1.signals[channel1[ch_idx1], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx2], :, epoch2[ep_idx]])
                f[ch_idx1, ch_idx2, ep_idx] = ftest.F
                p[ch_idx1, ch_idx2, ep_idx] = pvalue(ftest)
            end
        end
    end

    return (f=f, p=p)
end

"""
    band_mpower(obj; channel, f, mt)

Calculate mean and maximum band power and its frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power [μV^2/Hz] per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power [Hz] per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency [μV^2/Hz] per channel per epoch
"""
function band_mpower(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}, mt::Bool=false)

    fs = sr(obj)
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be < $(fs / 2)."))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    mbp = zeros(ch_n, ep_n)
    maxfrq = zeros(ch_n, ep_n)
    maxbp = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            mbp[ch_idx, ep_idx], maxfrq[ch_idx, ep_idx], maxbp[ch_idx, ep_idx] = @views s_band_mpower(obj.data[channel[ch_idx], :, ep_idx], fs=fs, f=f, mt=mt)
        end
    end

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)
end

"""
    rel_psd(obj; channel, norm, mt, f)

Calculate relative power spectrum density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `f::Union{Tuple{Real, Real}, Nothing}=nothing`: calculate power relative to frequency range or total power

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`:powers
- `psd_frq::Array{Float64, 3}`: frequencies
"""
function rel_psd(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), norm::Bool=false, mt::Bool=false, f::Union{Tuple{Real, Real}, Nothing}=nothing)

    fs = sr(obj)
    if f !== nothing
        f = tuple_order(f)
        f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
        f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be < $(fs / 2)."))
    end

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    psd_tmp, psd_frq = @views s_rel_psd(obj.data[1, :, 1], fs=fs, norm=norm, mt=mt, f=f)
    psd_pow = zeros(ch_n, length(psd_tmp), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        for ch_idx in 1:ch_n
            psd_pow[ch_idx, :, ep_idx], _ = @views s_rel_psd(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=mt, f=f)
        end
    end

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    fbsplit(obj; channel, order, window)

Split EEG signal into frequency bands.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `order::Int64=8`: number of taps for FIR band-pass filter
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter; default is Hamming window, number of taps is calculated using fred harris' rule-of-thumb

# Returns

Named tuple containing:
- `band_names::Vector{Symbol}`
- `band_frq::Vector{Tuple{Real, Real}}`
- `signal_split::Array{Float64, 4}`
"""
function fbsplit(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), order::Int64=8, window::Union{Nothing, AbstractVector, Int64}=nothing)
    
    band = [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    signal_split = zeros(length(band), ch_n, epoch_len(obj), ep_n)
    band_frq = Vector{Tuple{Real, Real}}()

    @inbounds for band_idx in eachindex(band)
        band_f = band(obj, band=band[band_idx])
        push!(band_frq, band_f)
        flt = s_filter_create(fs=fs, fprototype=:fir, ftype=:bp, cutoff=band_f, order=order, window=window, n=epoch_len(obj))
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                signal_split[band_idx, ch_idx, :, ep_idx] = @views s_filter_apply(obj.data[channel[ch_idx], :, ep_idx], flt=flt)
            end
        end
    end

    return (band_names=band, band_frq=band_frq, signal_split=signal_split)
end

"""
    chdiff(obj1, obj2; channel1, channel2, epoch1, epoch2)

Subtract channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

- `ch_diff::Matrix{Float64}`
"""
function chdiff(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    ch_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ch_diff[ch_idx, :, ep_idx] = @views obj1.signals[channel1[ch_idx], :, epoch1[ep_idx]] .- obj2.signals[channel2[ch_idx], :, epoch2[ep_idx]]
        end
    end

    return ch_diff
end

"""
    cps(obj; channel, norm)

Calculate cross power spectrum.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `cps_pw::Array{Float64, 4}`: cross power spectrum power
- `cps_ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), norm::Bool=true)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    
    cps_pw_tmp, cps_ph_tmp, cps_fq = @views s2_cps(obj.data[1, :, 1], obj.data[1, :, 1], fs=fs)
    cps_pw = zeros(ch_n, ch_n, length(cps_pw_tmp), ep_n)
    cps_ph = zeros(ch_n, ch_n, length(cps_ph_tmp), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
           for ch_idx2 in 1:ch_idx1
                cps_pw[ch_idx1, ch_idx2, :, ep_idx], cps_ph[ch_idx1, ch_idx2, :, ep_idx], _ = @views s2_cps(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx], fs=fs, norm=norm)
            end

        # update progress bar
        progress_bar == true && next!(p)
        end
    end

    @inbounds @simd for time_idx in 1:size(cps_pw, 3)
        Threads.@threads for ep_idx in 1:ep_n
            for ch_idx1 in 1:(ch_n - 1)
                for ch_idx2 in (ch_idx1 + 1):ch_n
                    cps_pw[ch_idx1, ch_idx2, time_idx, ep_idx] = @views cps_pw[ch_idx2, ch_idx1, time_idx, ep_idx]
                    cps_ph[ch_idx1, ch_idx2, time_idx, ep_idx] = @views cps_ph[ch_idx2, ch_idx1, time_idx, ep_idx]
                end
            end
        end
    end

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)
end

"""
    cps(obj1, obj2; channel1, channel2, epoch1, epoch2, norm)

Calculate cross power spectrum.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `cps_pw::Vector{Float64}`: cross power spectrum power
- `cps_ph::Vector{Float64}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)), norm::Bool=true)

    sr(obj1) == sr(obj2) || throw(ArgumentError("EEG1 and EEG2 must have the same sampling rate."))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)
    fs = sr(obj1)

    cps_pw, cps_ph, cps_fq = @views s2_cps(obj1.signals[1, :, 1], obj2.signals[1, :, 1], fs=fs, norm=norm)

    cps_pw = zeros(ch_n, length(cps_pw), ep_n)
    cps_ph = zeros(ch_n, length(cps_ph), ep_n)
    cps_fq = zeros(ch_n, length(cps_fq), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            cps_pw[ch_idx, :, ep_idx], cps_ph[ch_idx, :, ep_idx], cps_fq[ch_idx, :, ep_idx] = @views s2_cps(obj1.signals[channel1[ch_idx], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx], :, epoch2[ep_idx]], fs=fs, norm=norm)
        end
    end

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)
end

"""
    phdiff(obj; channel, pad, h)

Calculate phase difference between EEG channels and mean phase of reference `channel`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of reference channels, default is all EEG/MEG channels except the analyzed one
- `avg::Symbol=:phase`: method of averaging:
    - `:phase`: phase is calculated for each reference channel separately and then averaged
    - `:signal`: signals are averaged prior to phase calculation
- `pad::Int64=0`: pad signals with 0s
- `h::Bool=false`: use FFT or Hilbert transformation

# Returns
 
- `ph_diff::Array{Float64, 3}`
"""
function phdiff(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), avg::Symbol=:phase, pad::Int64=0, h::Bool=false)

    avg in [:phase, :signal] || throw(ArgumentError("avg must be :phase or :signal."))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    ph_diff = zeros(ch_n, epoch_len(obj), ep_n)
    if avg === :phase
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                ref_channels = setdiff(channel, ch_idx)
                ph_ref = zeros(length(ref_channels), epoch_len(obj))
                for ref_idx in eachindex(ref_channels)
                    if h == true
                        _, _, _, ph = @views s_hspectrum(obj.data[ref_channels[ref_idx], :, ep_idx], pad=pad)
                    else
                        _, _, _, ph = @views s_spectrum(obj.data[ref_channels[ref_idx], :, ep_idx], pad=pad)
                    end
                    ph_ref[ref_idx, :] = ph
                end
                ph_ref = vec(mean(ph_ref, dims=1))
                if h == true
                    _, _, _, ph = @views s_hspectrum(obj.data[channel[ch_idx], :, ep_idx], pad=pad)
                else
                    _, _, _, ph = @views s_spectrum(obj.data[channel[ch_idx], :, ep_idx], pad=pad)
                end
                ph_diff[ch_idx, :, ep_idx] = ph - ph_ref
            end
        end
    else
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                ref_channels = setdiff(channel, ch_idx)
                signal_m = @views vec(mean(obj.data[ref_channels, :, ep_idx], dims=1))
                ph_diff[ch_idx, :, ep_idx] = @views s_phdiff(obj.data[channel[ch_idx], :, ep_idx], signal_m)
            end
        end
    end

    return ph_diff
end

"""
    ampdiff(obj; channel)

Calculate amplitude difference between each channel and mean amplitude.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of reference channels, default is all EEG/MEG channels except the analyzed one

# Returns
 
- `amp_diff::Array{Float64, 3}`
"""
function ampdiff(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    amp_diff = zeros(ch_n, epoch_len(obj), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ref_channels = setdiff(channel, ch_idx)
            amp_ref = @views vec(mean(obj.data[ref_channels, :, ep_idx], dims=1))
            amp_diff[ch_idx, :, ep_idx] = @views obj.data[channel[ch_idx], :, ep_idx] - amp_ref
        end
    end

    return amp_diff
end

"""
    dwt(obj; channel, wt, type, l)

Perform discrete wavelet transformation (DWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type: 
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `l::Int64=0`: number of levels, default is maximum number of levels available or total transformation

# Returns
 
- `dwt_c::Array{Float64, 4}`: DWT coefficients cAl, cD1, ..., cDl (by rows)
"""
function dwt(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), wt::T, type::Symbol, l::Int64=0) where {T <: DiscreteWavelet}

    if l == 0
        l = maxtransformlevels(obj.data[1, :, 1])
        _info("Calculating DWT using maximum level: $l.")
    end

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    dwt_c = zeros(ch_n, (l + 1), epoch_len(obj), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            dwt_c[ch_idx, :, :, ep_idx] = @views s_dwt(obj.data[channel[ch_idx], :, ep_idx], wt=wt, type=type, l=l)
        end
    end

    return dwt_c
end

"""
    cwt(obj; channel, wt)

Perform continuous wavelet transformation (CWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns
 
- `cwt_c::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function cwt(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), wt::T) where {T <: CWT}

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    l = size(ContinuousWavelets.cwt(obj.data[1, :, 1], wt), 2)
    cwt_c = zeros(ch_n, l, epoch_len(obj), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            cwt_c[ch_idx, :, :, ep_idx] = @views s_cwt(obj.data[channel[ch_idx], :, ep_idx], wt=wt)
        end
    end

    return cwt_c
end

"""
    psdslope(obj; channel, f, norm, mt)

Calculate PSD linear fit and slope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `f::Tuple{Real, Real}=(0, sr(obj)/2)`: calculate slope of the total power (default) or frequency range f[1] to f[2]
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `lf::Array{Float64, 3}`: linear fit for each channel and epoch
- `psd_slope::Array{Float64, 2}`: slopes of each linear fit
- `frq::Vector{Float64}`: range of frequencies for the linear fits
"""
function psdslope(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}=(0, sr(obj)/2), norm::Bool=false, mt::Bool=false, nt::Int64=8)

    fs = sr(obj)
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < $(fs / 2)."))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    _, frq = s_psd(obj.data[1, :, 1], fs=fs, norm=norm, mt=mt, nt=nt)
    f1_idx = vsearch(f[1], frq)
    f2_idx = vsearch(f[2], frq)
    lf = zeros(ch_n, length(frq[f1_idx:f2_idx]), ep_n)
    psd_slope = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pow, _ = s_psd(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=mt, nt=nt)
            _, _, _, _, _, _, lf[ch_idx, :, ep_idx] = @views linreg(frq[f1_idx:f2_idx], pow[f1_idx:f2_idx])
            psd_slope[ch_idx, ep_idx] = lf[ch_idx, 2, ep_idx] - lf[ch_idx, 1, ep_idx]
        end
    end

    return (lf=lf, psd_slope=psd_slope, frq=frq[f1_idx:f2_idx])
end

"""
    henv(obj; channel, d)

Calculate Hilbert spectrum amplitude envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env::Array{Float64, 3}`: Hilbert spectrum amplitude envelope
- `s_t::Vector{Float64}`: signal time
"""
function henv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), d::Int64=32)

    _check_channels(obj, channel)
    _, signal, _, _ = @views spectrum(keep_channel(obj, channel=channel), h=true)

    ch_n = size(signal, 1)
    ep_n = size(signal, 3)
    h_env = similar(signal)
    s_t = obj.epoch_time

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s = @view signal[ch_idx, :, ep_idx]
            # find peaks
            p_idx = s_findpeaks(s, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(s))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_t[p_idx], s[p_idx])
                try
                    h_env[ch_idx, :, ep_idx] = model(s_t)
                catch
                    @error "CubicSpline error, using Loess."
                    model = loess(s_t[p_idx], s[p_idx], span=0.5)
                    h_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
                end
            else
                _info("Less than 5 peaks detected, using Loess.")
                model = loess(s_t[p_idx], s[p_idx], span=0.5)
                h_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
            end
            h_env[ch_idx, 1, ep_idx] = h_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (h_env=h_env, s_t=s_t)
end

"""
    henv_mean(obj; channel, dims, d)

Calculate Hilbert spectrum amplitude envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: mean
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function henv_mean(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = henv(obj, channel=channel, d=d)
    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # mean over channels

        h_env_m = zeros(length(s_t), ep_n)
        h_env_u = zeros(length(s_t), ep_n)
        h_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            h_env_m[:, ep_idx] = mean(s_a[:, :, ep_idx], dims=1)
            s = std(h_env_m[:, ep_idx]) / sqrt(length(h_env_m[:, ep_idx]))
            h_env_u[:, ep_idx] = @. h_env_m[:, ep_idx] + 1.96 * s
            h_env_l[:, ep_idx] = @. h_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        h_env_m = zeros(length(s_t), ch_n)
        h_env_u = zeros(length(s_t), ch_n)
        h_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            h_env_m[:, ch_idx] = mean(s_a[ch_idx, :, :], dims=2)
            s = std(h_env_m[:, ch_idx]) / sqrt(length(h_env_m[:, ch_idx]))
            h_env_u[:, ch_idx] = @. h_env_m[:, ch_idx] + 1.96 * s
            h_env_l[:, ch_idx] = @. h_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        h_env_m, h_env_u, h_env_l, _ = henv_mean(obj, dims=1, d=d)
        h_env_m = mean(h_env_m, dims=2)
        h_env_u = mean(h_env_u, dims=2)
        h_env_l = mean(h_env_l, dims=2)
        h_env_m = reshape(h_env_m, size(h_env_m, 1))
        h_env_u = reshape(h_env_u, size(h_env_u, 1))
        h_env_l = reshape(h_env_l, size(h_env_l, 1))
    end

    return (h_env_m=h_env_m, h_env_u=h_env_u, h_env_l=h_env_l, s_t=s_t)
end

"""
    henv_median(obj; channel, dims, d)

Calculate Hilbert spectrum amplitude envelope of `obj`: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: median
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function henv_median(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = henv(obj, channel=channel, d=d)
    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # median over channels

        h_env_m = zeros(length(s_t), ep_n)
        h_env_u = zeros(length(s_t), ep_n)
        h_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            h_env_m[:, ep_idx] = median(s_a[:, :, ep_idx], dims=1)
            t_idx = s_findpeaks(h_env_m[:, ep_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(h_env_m[:, ep_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], h_env_m[t_idx])
                try
                    h_env_m[:, ep_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(h_env_m[:, ep_idx]) / sqrt(length(h_env_m[:, ep_idx]))
            h_env_u[:, ep_idx] = @. h_env_m[:, ep_idx] + 1.96 * s
            h_env_l[:, ep_idx] = @. h_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        h_env_m = zeros(length(s_t), ch_n)
        h_env_u = zeros(length(s_t), ch_n)
        h_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            h_env_m[:, ch_idx] = median(s_a[ch_idx, :, :], dims=2)
            t_idx = s_findpeaks(h_env_m[:, ch_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(h_env_m[:, ch_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], h_env_m[t_idx])
                try
                    h_env_m[:, ch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(h_env_m[:, ch_idx]) / sqrt(length(h_env_m[:, ch_idx]))
            h_env_u[:, ch_idx] = @. h_env_m[:, ch_idx] + 1.96 * s
            h_env_l[:, ch_idx] = @. h_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

        h_env_m, h_env_u, h_env_l, _ = henv_median(obj, dims=1, d=d)
        h_env_m = median(h_env_m, dims=2)
        h_env_u = median(h_env_u, dims=2)
        h_env_l = median(h_env_l, dims=2)
        h_env_m = reshape(h_env_m, size(h_env_m, 1))
        h_env_u = reshape(h_env_u, size(h_env_u, 1))
        h_env_l = reshape(h_env_l, size(h_env_l, 1))
    end

    return (h_env_m=h_env_m, h_env_u=h_env_u, h_env_l=h_env_l, s_t=s_t)
end

"""
    apply(obj; channel, f)

Apply custom function.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `f::String`: function to be applied, e.g. `f="mean(obj, dims=3)"; EEG signal is given using variable `obj` here.

# Returns

- `out::Array{Float64, 3}`
"""
function apply(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), f::String)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    f_tmp = replace(f, "obj" => "$(obj.data[1, :, 1])")
    out_tmp = eval(Meta.parse(f_tmp))
    out = zeros(eltype(out_tmp), ch_n, length(out_tmp), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ch_n * ep_n, 1))
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            f_tmp = replace(f, "obj" => "$(obj.data[channel[ch_idx], :, ep_idx])")
            try
                out[ch_idx, :, ep_idx] = eval(Meta.parse(f_tmp))

            catch
                @error "Formula is incorrect."
            end
            # update progress bar
            progress_bar == true && next!(p)
        end
    end
    return out
end

"""
    channels_cluster(obj, cluster)

Return channels belonging to a `cluster` of channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `cluster::Symbol`: available clusters are:
    - `:f1`: left frontal (F1, F3, F5, F7, AF3, AF7)
    - `:f2`: right frontal (F2, F4, F6, F8, AF4, AF8)
    - `:t1`: left temporal (C3, C5, T7, FC3, FC5, FT7)
    - `:t2`: right temporal (C4, C6, T8, FC4, FC6, FT8)
    - `:c1`: anterior central (Cz, C1, C2, FC1, FC2, FCz)
    - `:c2`: posterior central (Pz, P1, P2, CP1, CP2, CPz)
    - `:p1`: left parietal (P3, P5, P7, CP3, CP5, TP7)
    - `:p2`: right parietal (P4, P6, P8, CP4, CP6, TP8)
    - `:o`: occipital (Oz, O1, O2, POz, PO3, PO4)

# Returns

- `channels::Vector{Int64}`: list of channel numbers belonging to a given cluster of channels
"""
function channel_cluster(obj::NeuroAnalyzer.NEURO; cluster::Symbol)

    length(labels(obj)) == 0 && throw(ArgumentError("EEG does not contain channel labels."))

    _check_var(cluster, [:f1, :f2, :t1, :t2, :c1, :c2, :p1, :p2, :o], "cluster")
    labels = lowercase.(labels(obj))
    channels = Int64[]

    cluster === :f1 && (cluster = ["fp1", "f1", "f3", "f5", "f7", "f9", "af3", "af7"])
    cluster === :f2 && (cluster = ["fp2", "f2", "f4", "f6", "f8", "f10", "af4", "af8"])
    cluster === :t1 && (cluster = ["c3", "c5", "t7", "t9", "fc3", "fc5", "ft7", "ft9"])
    cluster === :t2 && (cluster = ["c4", "c6", "t8", "t10", "fc4", "fc6", "ft8", "ft10"])
    cluster === :c1 && (cluster = ["cz", "c1", "c2", "fc1", "fc2", "fcz"])
    cluster === :c2 && (cluster = ["pz", "p1", "p2", "cp1", "cp2", "cpz"])
    cluster === :p1 && (cluster = ["p3", "p5", "p7", "p9", "cp3", "cp5", "tp7", "tp9"])
    cluster === :p2 && (cluster = ["p4", "p6", "p8", "p10", "cp4", "cp6", "tp8", "tp10"])
    cluster === :o && (cluster = ["o1", "o2", "poz", "po3", "po4", "po7", "po8", "po9", "po10"])

    for idx in cluster
        idx in labels && push!(channels, get_channel(obj, channel=idx))
    end

    return channels
end

"""
    erp_peaks(obj)

Detect a pair of positive and negative peaks of ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`:

# Returns
 
- `p::Array{Int64, 2}`: peaks: channels × positive peak position, negative peak position
"""
function erp_peaks(obj::NeuroAnalyzer.NEURO)

    channels = signal_channels(obj)
    erp = erp(obj).signals[channels, :]

    ch_n = size(erp, 1)
    p = zeros(Int64, ch_n, 2)
    @inbounds @simd for ch_idx in 1:ch_n
        pp_pos = @views maximum(erp[ch_idx, :])
        pp_neg = @views minimum(erp[ch_idx, :])
        p[ch_idx, :] = @views [vsearch(pp_pos, erp[ch_idx, :]), vsearch(pp_neg, erp[ch_idx, :])]
    end

    return p
end

"""
    bands_dwt(obj; channel, wt, type, n)

Split EEG channel into bands using discrete wavelet transformation (DWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64}`: channel number
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type: 
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `n::Int64=0`: number of bands, default is maximum number of bands available or total transformation

# Returns
 
- `bands::Array{Float64, 4}`: bands from lowest to highest frequency (by rows)
"""
function bands_dwt(obj::NeuroAnalyzer.NEURO; channel::Int64, wt::T, type::Symbol, n::Int64=0) where {T <: DiscreteWavelet}

    n -= 1
    if n == 0
        n = maxtransformlevels(obj.data[1, :, 1])
        _info("Calculating DWT using maximum level: $n.")
    end
    n < 2 && throw(ArgumentError("n must be ≥ 2."))

    _check_channels(obj, channel)
    ep_n = epoch_n(obj)

    dwt_c = zeros((n + 1), epoch_len(obj), ep_n)
    Threads.@threads for ep_idx in 1:ep_n
        @inbounds dwt_c[:, :, ep_idx] = @views s_dwt(obj.data[channel, :, ep_idx], wt=wt, type=type, l=n)
    end
    
    bands = similar(dwt_c)
    bands[1, :, :] = dwt_c[1, :, :]
    bands[2:end, :, :] = dwt_c[end:-1:2, :, :]

    return bands
end