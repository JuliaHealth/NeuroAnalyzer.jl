export plot_filter_response

"""
    plot_filter_response(<keyword arguments>)

Plot filter response.

# Arguments

- `fs::Int64`: sampling rate
- `n::Int64`: signal length in samples
- `fprototype::Symbol`: filter prototype:
    - `:fir`: FIR filter
    - `:firls`: weighted least-squares FIR filter
    - `:remez`: Remez FIR filter
    - `:butterworth`: IIR filter
    - `:chebyshev1` IIR filter
    - `:chebyshev2` IIR filter
    - `:elliptic` IIR filter
    - `:iirnotch`: second-order IIR notch filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (must be a pair of frequencies for `:bp` and `:bs`)
- `fs::Int64`: signal sampling rate
- `order::Int64`: filter order
- `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
- `mono::Bool=false`: use color or gray palette
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the X-axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_filter_response(; fs::Int64, n::Int64, fprototype::Symbol, ftype::Union{Nothing, Symbol}=nothing, cutoff::Union{Real, Tuple{Real, Real}}, order::Int64, rp::Union{Nothing, Real}=nothing, rs::Union{Nothing, Real}=nothing, bw::Union{Nothing, Real}=nothing, w::Union{Nothing, AbstractVector}=nothing, mono::Bool=false, frq_lim::Tuple{Real, Real}=(0, fs / 2), kwargs...)::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :darktest
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    flt = filter_create(fprototype=fprototype, ftype=ftype, cutoff=cutoff, n=n, fs=fs, order=order, rp=rp, rs=rs, bw=bw, w=w)

    xt = round.(linspace(frq_lim[1], frq_lim[2], 10), digits=1)
    xt[end] = frq_lim[2]

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch]
        H, w = freqresp(flt)
        # convert to dB
        H = 20 * log10.(abs.(H))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi

        if fprototype !== :iirnotch
            fname = titlecase(String(fprototype))
            title = "Filter: $(fname), type: $(uppercase(String(ftype))), cutoff: $(round.(cutoff, digits=1)) Hz, order: $order\n\nFrequency response"
        else
            fname = "IIR notch"
            title = "Filter: $(fname), cutoff: $(round.(cutoff, digits=1)) Hz, transition band width: $bw Hz\n\nFrequency response"
        end

        p1 = Plots.plot(w,
                        H,
                        title=title,
                        xlims=frq_lim,
                        ylims=(-100, 20),
                        xticks=(xt, string.(xt)),
                        ylabel="Magnitude\n[dB]",
                        xlabel="Frequency [Hz]",
                        label="",
                        bottom_margin=10Plots.px,
                        left_margin=40Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            Plots.plot!((0, cutoff),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:black)
        else
            Plots.plot!((0, cutoff[1]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:red)
            Plots.plot!((0, cutoff[2]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:green)
        end

        phi, w = phaseresp(flt)
        phi = rad2deg.(phi)
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        p2 = Plots.plot(w,
                        phi,
                        title="Phase response",
                        xlims=frq_lim,
                        xticks=(xt, string.(xt)),
                        ylabel="Phase\n[deg]",
                        xlabel="Frequency [Hz]",
                        label="",
                        bottom_margin=10Plots.px,
                        left_margin=40Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            Plots.plot!((0, cutoff),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:black)
        else
            Plots.plot!((0, cutoff[1]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:red)
            Plots.plot!((0, cutoff[2]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:green)
        end

        tau = -derivative(phi)
        p3 = Plots.plot(w,
                        tau,
                        title="Group delay",
                        xlims=frq_lim,
                        xticks=(xt, string.(xt)),
                        ylabel="Group delay\n[samples]",
                        xlabel="Frequency [Hz]",
                        label="",
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            Plots.plot!((0, cutoff),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:black)
        else
            Plots.plot!((0, cutoff[1]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:red)
            Plots.plot!((0, cutoff[2]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:green)
        end

        p = Plots.plot(p1, p2, p3, size=(1200, 800), margins=20Plots.px, layout=(3, 1), palette=pal; kwargs...)
    else
        w = range(0, stop=pi, length=1024)
        H = _fir_response(flt, w)
        # convert to dB
        H = amp2db.(abs.(H))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        if fprototype === :fir
            title = "Filter: FIR, type: $(uppercase(String(ftype))), cutoff: $(round.(cutoff, digits=1)) Hz, order: $order\n\nFrequency response"
        elseif fprototype === :firls
            title = "Filter: FIR (LS), type: $(uppercase(String(ftype))), cutoff: $(round.(cutoff, digits=1)) Hz, transition band width: $bw Hz, order: $order\n\nFrequency response"
        elseif fprototype === :remez
            title = "Filter: Remez, type: $(uppercase(String(ftype))), cutoff: $(round.(cutoff, digits=1)) Hz, transition band width: $bw Hz, order: $order\n\nFrequency response"
        end
        p1 = Plots.plot(w,
                        H,
                        title=title,
                        xlims=frq_lim,
                        ylims=(-100, 20),
                        xticks=(xt, string.(xt)),
                        ylabel="Magnitude\n[dB]",
                        xlabel="Frequency [Hz]",
                        label="",
                        bottom_margin=10Plots.px,
                        left_margin=40Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            Plots.plot!((0, cutoff),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:black)
        else
            Plots.plot!((0, cutoff[1]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:red)
            Plots.plot!((0, cutoff[2]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:green)
        end

        w = range(0, stop=pi, length=1024)
        phi = _fir_response(flt, w)
        phi = rad2deg.(-atan.(imag(phi), real(phi)))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        p2 = Plots.plot(w,
                        phi,
                        title="Phase response",
                        xlims=frq_lim,
                        xticks=(xt, string.(xt)),
                        ylabel="Phase\n[deg]",
                        xlabel="Frequency [Hz]",
                        label="",
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            Plots.plot!((0, cutoff),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:black)
        else
            Plots.plot!((0, cutoff[1]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:red)
            Plots.plot!((0, cutoff[2]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:green)
        end

        tau = -derivative(phi)
        p3 = Plots.plot(w,
                        tau,
                        title="Group delay",
                        xlims=frq_lim,
                        xticks=(xt, string.(xt)),
                        ylabel="Group delay\n[samples]",
                        xlabel="Frequency [Hz]",
                        label="",
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            Plots.plot!((0, cutoff),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:black)
        else
            Plots.plot!((0, cutoff[1]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:red)
            Plots.plot!((0, cutoff[2]),
                        seriestype=:vline,
                        linestyle=:dash,
                        label="",
                        lw=0.5,
                        lc=:green)
        end

        p = Plots.plot(p1, p2, p3, size=(1200, 800), margins=20Plots.px, layout=(3, 1), palette=pal; kwargs...)
    end

    return p

end

"""
    plot_filter_response(obj, <keyword arguments>)

Plot filter response.

# Arguments

- `fs::Int64`: sampling rate
- `n::Int64`: signal length in samples
- `fprototype::Symbol`: filter prototype:
    - `:fir`: FIR filter
    - `:firls`: weighted least-squares FIR filter
    - `:remez`: Remez FIR filter
    - `:butterworth`: IIR filter
    - `:chebyshev1` IIR filter
    - `:chebyshev2` IIR filter
    - `:elliptic` IIR filter
    - `:iirnotch`: second-order IIR notch filter
- `ftype::Union{Nothing, Symbol}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (must be a pair of frequencies for `:bp` and `:bs`)
- `fs::Int64`: signal sampling rate
- `order::Int64`: filter order
- `rp::Union{Nothing, Real}=nothing`: maximum ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Union{Nothing, Real}=nothing`: minimum ripple attenuation in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz for `:firls`, `:remez` and `:iirnotch` filters
- `w::Union{Nothing, AbstractVector}=nothing`: window for `:fir` filter (default is Hamming window) or weights for `:firls` filter
- `mono::Bool=false`: use color or gray palette
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the X-axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_filter_response(obj::NeuroAnalyzer.NEURO; fprototype::Symbol, ftype::Union{Nothing, Symbol}=nothing, cutoff::Union{Real, Tuple{Real, Real}}, order::Int64, rp::Union{Nothing, Real}=nothing, rs::Union{Nothing, Real}=nothing, bw::Union{Nothing, Real}=nothing, w::Union{Nothing, AbstractVector}=nothing, mono::Bool=false, frq_lim::Tuple{Real, Real}=(0, fs / 2), kwargs...)::Plots.Plot{Plots.GRBackend}

    p = plot_filter_response(fs=sr(obj), n=epoch_len(obj), fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, w=w, mono=mono, frq_lim=frq_lim, kwargs=kwargs)

    return p

end
