export mplot_spectrogram
export mplot_spectrogram_topo

"""
    mplot_spectrogram(st, sf, sp; <keyword arguments>)

Plot single-channel spectrogram.

# Arguments

- `st::Vector{Float64}`: time
- `sf::Vector{<:Real}`: frequencies
- `sp::Matrix{Float64}`: powers
- `db::Bool=true`: whether powers are normalized to dB
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `units::String=""`
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `cb::Bool=true`: plot color bar
- `cb_title::String=""`: color bar label
- `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_spectrogram(st::Vector{Float64}, sf::Vector{<:Real}, sp::Matrix{Float64}; db::Bool=true, frq::Symbol=:lin, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, units::String="", smooth::Bool=false, n::Int64=3, cb::Bool=true, cb_title::String="", threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq, kwargs...)::GLMakie.Figure

    @assert size(sp, 2) == length(st) "Size of powers ($(size(sp, 2))) and time vector ($(length(st))) do not match."
    @assert size(sp, 1) == length(sf) "Size of powers ($(size(sp, 1))) and frequencies vector ($(length(sf))) do not match."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest
    if cb_title == ""
        cb_title = db ? "[dB $units^2/Hz]" : "[$units^2/Hz]"
    end

    if smooth
        sp = imfilter(sp, Kernel.gaussian(n))
    end

    if frq === :log && frq_lim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        frq_lim = (sf[2], frq_lim[2])
    end

    # prepare plot
    plot_size = (1200, 800)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(10),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      yticks=LinearTicks(15),
                      yminorticksvisible=true,
                      yminorticks=IntervalsBetween(10),
                      yscale=frq===:lin ? identity : log,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.xlims!(ax, (st[1], st[end]))
    GLMakie.ylims!(ax, frq_lim)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    hm = GLMakie.heatmap!(ax,
                          st,
                          sf,
                          sp',
                          colormap=pal)

    if !isnothing(threshold)
        _, bm = seg_extract(sp, threshold=threshold, threshold_type=threshold_type)
        reg = ones(size(sp)) .* minimum(sp)
        reg[bm] .= maximum(sp)
        GLMakie.contour!(ax,
                         st,
                         sf,
                         reg',
                         levels=1,
                         color=:black,
                         linewidth=2)
    end

    if cb
        Colorbar(p[1, 2],
                 hm,
                 label=cb_title,
                 labelsize=18)
    end

    return p

end

"""
    mplot_spectrogram(sf, sp; <keyword arguments>)

Plot multiple-channel spectrogram.

# Arguments

- `sf::Vector{<:Real}`: frequencies
- `sp::Matrix{Float64}`: powers
- `clabels::Vector{String}`: channel labels
- `db::Bool=true`: whether powers are normalized to dB
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the Y-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `units::String=""`
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `cb::Bool=true`: plot color bar
- `cb_title::String=""`: color bar label
- `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_spectrogram(sf::Vector{<:Real}, sp::Matrix{Float64}; clabels::Vector{String}=[""], db::Bool=true, frq::Symbol=:lin, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, units::String="", smooth::Bool=false, n::Int64=3, cb::Bool=true, cb_title::String="", threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq, kwargs...)::GLMakie.Figure

    @assert size(sp, 1) == length(clabels) "Size of powers ($(size(sp, 1))) and channels vector ($(length(clabels))) do not match."
    @assert size(sp, 2) == length(sf) "Size of powers ($(size(sp, 2))) and frequencies vector ($(length(sf))) do not match."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest
    if cb_title == ""
        cb_title = db ? "[dB $units^2/Hz]" : "[$units^2/Hz]"
    end

    if smooth
        sp = imfilter(sp, Kernel.gaussian(n))
    end

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(sp, 1)))

    if frq === :log && frq_lim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        frq_lim = (sf[2], frq_lim[2])
    end

    ch = collect(eachindex(clabels)) .- 0.5
    ch_n = length(ch)
    reverse!(sp, dims=1)

    # prepare plot
    plot_size = (1200, 800)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=LinearTicks(15),
                      xminorticksvisible=true,
                      xminorticks=IntervalsBetween(10),
                      yticks=(0.5:1:ch_n, reverse(clabels)),
                      yticksvisible=false,
                      xscale=frq===:lin ? identity : log,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.xlims!(ax, frq_lim)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    hm = GLMakie.heatmap!(ax,
                          sf,
                          ch,
                          sp',
                          colormap=pal)
    if !isnothing(threshold)
        _, bm = seg_extract(sp, threshold=threshold, threshold_type=threshold_type)
        reg = ones(size(sp)) .* minimum(sp)
        reg[bm] .= maximum(sp)
        GLMakie.contour!(ax,
                         sf,
                         ch,
                         reg',
                         levels=1,
                         color=:black,
                         linewidth=2)
    end

    if cb
        Colorbar(p[1, 2],
                 hm,
                 label=cb_title,
                 labelsize=18)
    end

    return p

end

"""
    mplot_spectrogram(obj; <keyword arguments>)

Plots spectrogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `db::Bool=true`: normalize powers to dB; for CWT scaleogram: normalize to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal
- `method::Symbol=:stft`: method of calculating spectrogram:
    - `:stft`: short-time Fourier
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
    - `:hht`: Hilbert-Huang transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `gw::Real=10`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: y-axis limits
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[channel: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or gray palette
- `markers::Bool`: draw markers if available
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `cb::Bool=true`: plot color bar
- `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
- `topo::Bool=false`: plot topographical map of spectrograms
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: plot head shape
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_spectrogram(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}=(0, 10), ep::Int64=0, ch::Union{String, Vector{String}, Regex}, db::Bool=true, method::Symbol=:stft, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.90), w::Bool=true, gw::Real=10, wt::T=wavelet(Morlet(2π), β=2), frq::Symbol=:lin, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), ncyc::Union{Int64, Tuple{Int64, Int64}}=32, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, smooth::Bool=false, n::Int64=3, cb::Bool=true, threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq, topo::Bool=false, cart::Bool=false, head::Bool=true, kwargs...)::GLMakie.Figure where {T <: CWT}

    _check_var(method, [:stft, :mt, :mw, :gh, :cwt, :hht], "method")
    @assert seg[1] != seg[2] "Signal is too short for analysis."
    @assert n > 0 "n must be ≥ 1."

    ch = get_channel(obj, ch=ch)
    if method === :cwt
        if !topo
            @assert length(ch) == 1 "For :cwt method only one channel must be selected."
        end
    end
    if topo
        @assert method !== :hht "For :hht method topographical map is not available."
    end

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    if ep != 0
        _check_epochs(obj, ep)
        if nepochs(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    # get time vector
    if seg[2] <= epoch_len(obj)
        signal = obj.data[ch, seg[1]:seg[2], 1]
    else
        signal = epoch(obj, ep_n=1).data[ch, seg[1]:seg[2], 1]
    end
    length(ch) == 1 && (signal = vec(signal))

    # t = _get_t(seg[1], seg[2], sr(obj))
    t = obj.time_pts[seg[1]:seg[2]]
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    clabels = labels(obj)[ch]

    # set units
    units = _ch_units(obj, labels(obj)[ch[1]])

    # get frequency range
    fs = sr(obj)
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    # calculate spectrogram
    length(ch) > 1 && @assert length(signal) / length(ch) >= 4 * sr(obj) "For multi-channel plot, signal length must be ≥ 4 × sampling rate (4 × $(sr(obj)) samples)."

    if length(ch) == 1
        xlabel == "default" && (xlabel = "Time [s]")
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        if method === :stft
            sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=false, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Spectrogram (short-time Fourier)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=false, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            title == "default" && (title = "Spectrogram (multi-tapered)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            _, sp, _, sf, st = NeuroAnalyzer.mwspectrogram(signal, fs=fs, ncyc=ncyc, db=false, w=w)
            title == "default" && (title = "Spectrogram (Morlet wavelet)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, _, sf, st = NeuroAnalyzer.ghtspectrogram(signal, fs=fs, db=false, gw=gw, w=w)
            title == "default" && (title = "Spectrogram (Gaussian-Hilbert)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            _log_off()
            sp, sf, st = NeuroAnalyzer.cwtspectrogram(signal, fs=fs, wt=wt)
            _log_on()
            sf[1] > frq_lim[1] && (frq_lim = (sf[1], frq_lim[2]))
            sf[end] < frq_lim[2] && (frq_lim = (frq_lim[1], sf[end]))
            title == "default" && (title = "CWT Scaleogram\n[epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :hht
            imf = emd(signal, t)
            sp, _, sf, st = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs=fs, db=false)
            title == "default" && (title = "Spectrogram (Hilbert-Huang)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        end

        f1 = vsearch(frq_lim[1], sf)
        f2 = vsearch(frq_lim[2], sf)
        sf = sf[f1:f2]
        sp = sp[f1:f2, :]

        st .+= t[1]

        method !== :cwt && db && (sp = pow2db.(sp))
        p = mplot_spectrogram(st,
                              sf,
                              sp,
                              db=db,
                              frq=frq,
                              frq_lim=frq_lim,
                              xlabel=xlabel,
                              ylabel=ylabel,
                              title=title,
                              mono=mono,
                              units=units,
                              smooth=smooth,
                              n=n,
                              cb=cb,
                              cb_title=method === :cwt ? "Magnitude" : "",
                              threshold=threshold,
                              threshold_type=threshold_type;
                              kwargs...)

        # plot markers if available
        # TODO: draw markers length
        if markers && _has_markers(obj)
            markers_pos = obj.markers[!, :start]
            markers_id = obj.markers[!, :id]
            markers_desc = obj.markers[!, :value]
            for idx in eachindex(markers_pos)
                if _in(markers_pos[idx], (st[1], st[end]))
                    GLMakie.vlines!(markers_pos[idx],
                                    linestyle=:dash,
                                    linewidth=1,
                                    color=:black)
                    GLMakie.text!(markers_pos[idx] + 0.1,
                                  0.5,
                                  text="$(markers_id[idx]) / $(markers_desc[idx])",
                                  fontsize=5,
                                  align=(:left, :top),
                                  rotation=pi/2)
                end
            end
        end

    else
        ch_t = unique(obj.header.recording[:channel_type][ch])
        @assert length(ch_t) == 1 "For multi-channel spectrogram plot, all channels must be of the same type."
        if !topo
            ylabel == "default" && (ylabel = "")
            xlabel == "default" && (xlabel = "Frequency [Hz]")
            if method === :stft
                sp, sf = psd(signal, fs=fs, db=db, method=:stft, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
                title == "default" && (title = "Spectrogram (short-time Fourier)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :mt
                sp, sf = psd(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
                title == "default" && (title = "Spectrogram (multi-tapered)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :mw
                sp, sf = psd(signal, fs=fs, db=db, method=:mw, w=w, ncyc=ncyc)
                title == "default" && (title = "Spectrogram (Morlet wavelet)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :gh
                sp, sf = psd(signal, fs=fs, db=db, method=:gh, w=w, gw=gw)
                title == "default" && (title = "Spectrogram (Gaussian-Hilbert)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :cwt
                _log_off()
                sp, sf = psd(signal, fs=fs, method=:cwt, wt=wt)
                _log_on()
                sf[1] > frq_lim[1] && (frq_lim = (sf[1], frq_lim[2]))
                sf[end] < frq_lim[2] && (frq_lim = (frq_lim[1], sf[end]))
                title == "default" && (title = "CWT Scaleogram\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :hht
                imf = emd(signal[1, :], t)
                sp_tmp, _, sf, st = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs=fs, db=db)
                sp = zeros(size(signal, 1), length(sp_tmp))
                sp[1, :] = sp_tmp
                for idx in axes(signal, 1)[(begin + 1):end]
                    imf = emd(signal[idx, :], t)
                    sp[idx, :], _, _, _ = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs=fs, db=db)
                end
                title == "default" && (title = "Spectrogram (Hilbert-Huang)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            end

            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[:, f1:f2]

            p = mplot_spectrogram(sf,
                                  sp,
                                  clabels=clabels,
                                  db=db,
                                  frq=frq,
                                  frq_lim=frq_lim,
                                  xlabel=xlabel,
                                  ylabel=ylabel,
                                  title=title,
                                  mono=mono,
                                  units=units,
                                  smooth=smooth,
                                  n=n,
                                  cb=cb,
                                  threshold=threshold,
                                  threshold_type=threshold_type;
                                  kwargs...)
        else
            @assert length(ch) > 1 "For topographical plot, the number of channels must be >1."
            _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
            _has_locs(obj)
            chs = intersect(obj.locs[!, :label], labels(obj)[ch])
            locs = Base.filter(:label => in(chs), obj.locs)
            _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
            if method === :stft
                sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=false, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
                title == "default" && (title = "Spectrogram (short-time Fourier)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :mt
                sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=false, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
                title == "default" && (title = "Spectrogram (multi-tapered)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :mw
                _, sp, _, sf, st = NeuroAnalyzer.mwspectrogram(signal, fs=fs, ncyc=ncyc, db=false, w=w)
                title == "default" && (title = "Spectrogram (Morlet wavelet)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :gh
                sp, _, sf, st = NeuroAnalyzer.ghtspectrogram(signal, fs=fs, db=false, gw=gw, w=w)
                title == "default" && (title = "Spectrogram (Gaussian-Hilbert)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :cwt
                _log_off()
                sp, sf, st = NeuroAnalyzer.cwtspectrogram(signal, fs=fs, wt=wt)
                _log_on()
                sf[1] > frq_lim[1] && (frq_lim = (sf[1], frq_lim[2]))
                sf[end] < frq_lim[2] && (frq_lim = (frq_lim[1], sf[end]))
                title == "default" && (title = "CWT Scaleogram\n[epoch: $ep, time window: $t_s1:$t_s2]")
            elseif method === :hht
                imf = emd(signal, t)
                sp, _, sf, st = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs=fs, db=false)
                title == "default" && (title = "Spectrogram (Hilbert-Huang)\n[epoch: $ep, time window: $t_s1:$t_s2]")
            end

            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[f1:f2, :, :]

            st .+= t[1]
            method !== :cwt && db && (sp = pow2db.(sp))
            p = mplot_spectrogram_topo(locs,
                                       st,
                                       sf,
                                       sp,
                                       frq=frq,
                                       frq_lim=frq_lim,
                                       title=title,
                                       mono=mono,
                                       cart=cart,
                                       smooth=smooth,
                                       n=n,
                                       head=head;
                                       kwargs...)
        end
    end

    return p

end

"""
    mplot_spectrogram(obj, c; <keyword arguments>)

Plots spectrogram of embedded or external component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Union{Symbol, AbstractArray}`: component to plot
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `ep::Int64=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `db::Bool=true`: normalize powers to dB; for CWT scaleogram: normalize to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal
- `method::Symbol=:stft`: method of calculating spectrogram:
    - `:stft`: short-time Fourier
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
    - `:hht`: Hilbert-Huang transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `gw::Real=10`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: y-axis limits
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
- `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]\n[component: 1, epoch: 1, time window: 0 ms:10 s]
- `mono::Bool=false`: use color or gray palette
- `markers::Bool`: draw markers if available
- `units::String=""`
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `cb::Bool=true`: plot color bar
- `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_spectrogram(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; seg::Tuple{Real, Real}=(0, 10), ep::Union{Int64, AbstractRange}=1, c_idx::Union{Int64, Vector{Int64}, AbstractRange}, db::Bool=true, method::Symbol=:stft, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.90), w::Bool=true, frq::Symbol=:lin, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), gw::Real=10, wt::T=wavelet(Morlet(2π), β=2), ncyc::Union{Int64, Tuple{Int64, Int64}}=32, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, markers::Bool=true, smooth::Bool=false, n::Int64=3, cb::Bool=true, threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq, kwargs...)::GLMakie.Figure where {T <: CWT}

    _check_var(method, [:stft, :mt, :mw, :gh, :cwt, :hht], "method")
    @assert seg[1] != seg[2] "Signal is too short for analysis."
    @assert n > 0 "n must be ≥ 1."

    units = "A.U."

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    if ep != 0
        _check_epochs(obj, ep)
        if nepochs(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    # select component channel, default is all channels
    c isa Symbol && (c = _get_component(obj, c))
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    clabels = _gen_clabels(c)[c_idx]
    c_idx isa Int64 && (clabels = [clabels])

    # get time vector
    if seg[2] <= epoch_len(obj)
        signal = c[c_idx, seg[1]:seg[2], 1]
    else
        signal = reshape(c, size(c, 1), :, 1)[c_idx, seg[1]:seg[2], 1]
    end
    # t = _get_t(seg[1], seg[2], sr(obj))
    t = obj.time_pts[seg[1]:seg[2]]
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    # get frequency range
    fs = sr(obj)
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    # calculate spectrogram
    length(c_idx) > 1 && @assert length(signal) / length(c_idx) >= 4 * sr(obj) "For multi-channel plot, signal length must be ≥ 4 × sampling rate (4 × $(sr(obj)) samples)."

    if length(c_idx) == 1
        xlabel == "default" && (xlabel = "Time [s]")
        ylabel == "default" && (ylabel = "Frequency [Hz]")
        if method === :stft
            sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[f1:f2, :]
            title == "default" && (title = "Spectrogram (short-time Fourier)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf, st = NeuroAnalyzer.spectrogram(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[f1:f2, :]
            title == "default" && (title = "Spectrogram (multi-tapered)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            _, sp, _, sf = NeuroAnalyzer.mwspectrogram(signal, fs=fs, ncyc=ncyc, db=db, w=w)
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            title == "default" && (title = "Spectrogram (Morlet wavelet)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, _, sf = NeuroAnalyzer.ghtspectrogram(signal, fs=fs, db=db, gw=gw, w=w)
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            title == "default" && (title = "Spectrogram (Gaussian-Hilbert)\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            _log_off()
            sp, sf = NeuroAnalyzer.cwtspectrogram(signal, fs=fs, wt=wt, norm=db)
            _log_on()
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            sf[1] > frq_lim[1] && (frq_lim = (sf[1], frq_lim[2]))
            sf[end] < frq_lim[2] && (frq_lim = (frq_lim[1], sf[end]))
            title == "default" && (title = "CWT Scaleogram\n[component: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :hht
            imf = emd(signal, t)
            sp, _, sf, st = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs=fs, db=false)
            title == "default" && (title = "Spectrogram (Hilbert-Huang)\n[epoch: $ep, time window: $t_s1:$t_s2]")
        end

        st .+= t[1]

        p = mplot_spectrogram(st,
                              sf,
                              sp,
                              db=db,
                              frq=frq,
                              frq_lim=frq_lim,
                              xlabel=xlabel,
                              ylabel=ylabel,
                              title=title,
                              mono=mono,
                              units=units,
                              smooth=smooth,
                              n=n,
                              cb=cb,
                              threshold=threshold,
                              threshold_type=threshold_type;
                              kwargs...)

        # plot markers if available
        # TODO: draw markers length
        if markers && _has_markers(obj)
            markers_pos = obj.markers[!, :start]
            markers_id = obj.markers[!, :id]
            markers_desc = obj.markers[!, :value]
            for idx in eachindex(markers_pos)
                if _in(markers_pos[idx], (st[1], st[end]))
                    GLMakie.vlines!(markers_pos[idx],
                                    linestyle=:dash,
                                    linewidth=1,
                                    label=false)
                    GLMakie.text!(markers_pos[idx] + 0.1,
                                  0.5,
                                  text="$(markers_id[idx]) / $(markers_desc[idx])",
                                  fontsize=5,
                                  align=(:left, :top),
                                  rotation=pi/2)
                end
            end
        end

    else
        ylabel == "default" && (ylabel = "")
        xlabel == "default" && (xlabel = "Frequency [Hz]")
        if method === :stft
            sp, sf = psd(signal, fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[:, f1:f2]
            title == "default" && (title = "Spectrogram (short-time Fourier)\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mt
            sp, sf = psd(signal, fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            f1 = vsearch(frq_lim[1], sf)
            f2 = vsearch(frq_lim[2], sf)
            sf = sf[f1:f2]
            sp = sp[:, f1:f2]
            title == "default" && (title = "Spectrogram (multi-tapered)\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :mw
            sp, sf = mwpsd(signal, fs=fs, ncyc=ncyc, db=db, w=w)
            sf = linspace(0, frq_lim[2], size(sp, 2))
            title == "default" && (title = "Spectrogram (Morlet wavelet)\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :gh
            sp, _, sf = NeuroAnalyzer.ghtspectrogram(signal, fs=fs, db=db, gw=gw, w=w)
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            title == "default" && (title = "Spectrogram (Gaussian-Hilbert)\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :cwt
            _log_off()
            sp, sf = NeuroAnalyzer.cwtspectrogram(signal, fs=fs, wt=wt, norm=db)
            _log_on()
            st = linspace(0, (length(signal) / fs), size(sp, 2))
            sf[1] > frq_lim[1] && (frq_lim = (sf[1], frq_lim[2]))
            sf[end] < frq_lim[2] && (frq_lim = (frq_lim[1], sf[end]))
            title == "default" && (title = "CWT Scaleogram\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        elseif method === :hht
            imf = emd(signal[1, :], t)
            sp_tmp, _, sf, st = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs=fs, db=db)
            sp = zeros(size(signal, 1), length(sp_tmp))
            sp[1, :] = sp_tmp
            for idx in axes(signal, 1)[(begin + 1):end]
                imf = emd(signal[idx, :], t)
                sp[idx, :], _, _, _ = NeuroAnalyzer.hhtspectrogram(imf[1:(end - 1), :], fs=fs, db=db)
            end
            title == "default" && (title = "Spectrogram (Hilbert-Huang)\n[components: $(_channel2channel_name(c_idx)), epoch: $ep, time window: $t_s1:$t_s2]")
        end

        p = mplot_spectrogram(sf,
                             sp,
                             clabels=clabels,
                             db=db,
                             frq=frq,
                             frq_lim=frq_lim,
                             xlabel=xlabel,
                             ylabel=ylabel,
                             title=title,
                             mono=mono,
                             units=units,
                             smooth=smooth,
                             n=n,
                             cb=cb,
                             threshold=threshold,
                             threshold_type=threshold_type;
                             kwargs...)
    end

    return p

end

"""
    mplot_spectrogram_topo(locs, st, sf, sp; <keyword arguments>)

Plot topographical map of spectrograms.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `st::Vector{Float64}`: time
- `sf::Vector{Float64}`: frequencies
- `sp::Array{Float64, 3}`: powers
- `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `mono::Bool=false`: unused, for compatibility only
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `head::Bool=true`: plot head shape
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function mplot_spectrogram_topo(locs::DataFrame, st::Vector{Float64}, sf::Vector{Float64}, sp::Array{Float64, 3}; frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), title::String="", smooth::Bool=false, n::Int64=3, mono::Bool=true, frq::Symbol=:lin, cart::Bool=false, head::Bool=true, kwargs...)::GLMakie.Figure

    @assert size(sp, 3) == DataFrames.nrow(locs) "Size of powers ($(size(sp, 3))) and number of locs ($(DataFrames.nrow(locs))) do not match."
    @assert size(sp, 2) == length(st) "Size of powers ($(size(sp, 2))) and time vector ($(length(st))) do not match."
    @assert size(sp, 1) == length(sf) "Size of powers ($(size(sp, 1))) and frequencies vector ($(length(sf))) do not match."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    if frq === :log && frq_lim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        frq_lim = (sf[2], frq_lim[2])
    end

    # plot parameters
    if size(sp, 3) <= 64
        mplot_size = 1000
        marker_size = (120, 100)
        xl = 1.2
        yl = 1.2
    elseif _in(size(sp, 3), (64, 100))
        mplot_size = 1200
        marker_size = (120, 100)
        xl = 1.5
        yl = 1.5
    else
        mplot_size = 1500
        marker_size = (100, 80)
        xl = 1.5
        yl = 1.5
    end

    # get locations
    if !cart
        loc_x = zeros(size(locs, 1))
        loc_y = zeros(size(locs, 1))
        for idx in axes(locs, 1)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
    else
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
    end

    if smooth
        sp = imfilter(sp, Kernel.gaussian(n))
    end

    # prepare spectrogram plots
    pp_vec = GLMakie.Figure[]
    for idx in axes(sp, 3)
        pp = GLMakie.Figure(size=marker_size,
                            figure_padding=0)
        ax = GLMakie.Axis(pp[1, 1],
                          xlabel="",
                          ylabel="",
                          title=locs[idx, :label],
                          xautolimitmargin=(0, 0),
                          yautolimitmargin=(0, 0);
                          kwargs...)
        hidespines!(ax)
        hidedecorations!(ax)
        GLMakie.xlims!(ax, frq_lim)
        ax.titlesize = 8
        # plot powers
        GLMakie.heatmap!(sp[:, :, idx]',
                         colormap=pal)
        push!(pp_vec, pp)
    end

    # prepare plot
    plot_size = (mplot_size, mplot_size)
    p = GLMakie.Figure(size=plot_size,
                       figure_padding=0)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel="",
                      ylabel="",
                      title=title,
                      aspect=1,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.xlims!(ax, (-xl, xl))
    GLMakie.ylims!(ax, (-yl, yl))
    hidespines!(ax)
    hidedecorations!(ax)
    ax.titlesize = 20

    if head
        # nose
        GLMakie.lines!(ax, [-0.1, 0], [0.995, 1.1], linewidth=3, color=:black)
        GLMakie.lines!(ax, [0, 0.1], [1.1, 0.995], linewidth=3, color=:black)

        # ears
        # left
        GLMakie.lines!(ax, [-0.995, -1.03], [0.1, 0.15], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.03, -1.06], [0.15, 0.16], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.06, -1.1], [0.16, 0.14], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.1, -1.12], [0.14, 0.05], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.12, -1.10], [0.05, -0.1], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.10, -1.13], [-0.1, -0.3], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.13, -1.09], [-0.3, -0.37], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.09, -1.02], [-0.37, -0.39], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-1.02, -0.98], [-0.39, -0.33], linewidth=3, color=:black)
        GLMakie.lines!(ax, [-0.98, -0.975], [-0.33, -0.22], linewidth=3, color=:black)
        # right
        GLMakie.lines!(ax, [0.995, 1.03], [0.1, 0.15], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.03, 1.06], [0.15, 0.16], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.06, 1.1], [0.16, 0.14], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.1, 1.12], [0.14, 0.05], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.12, 1.10], [0.05, -0.1], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.10, 1.13], [-0.1, -0.3], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.13, 1.09], [-0.3, -0.37], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.09, 1.02], [-0.37, -0.39], linewidth=3, color=:black)
        GLMakie.lines!(ax, [1.02, 0.98], [-0.39, -0.33], linewidth=3, color=:black)
        GLMakie.lines!(ax, [0.98, 0.975], [-0.33, -0.22], linewidth=3, color=:black)

        # head
        GLMakie.arc!(ax,(0, 0), 1, 0, 2pi, linewidth=3, color=:black)
    end

    for idx in axes(sp, 3)
        io = IOBuffer()
        show(io, MIME"image/png"(), pp_vec[idx])
        pp = FileIO.load(io)
        GLMakie.scatter!(loc_x[idx],
                         loc_y[idx],
                         marker=pp,
                         markersize=marker_size,
                         markerspace=:pixel)
    end

    return p

end
