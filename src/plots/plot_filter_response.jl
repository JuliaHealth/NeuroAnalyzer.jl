export plot_filter_response

"""
    plot_filter_response(<keyword arguments>)

Plot filter response.

# Arguments

- `fs::Int64`: sampling rate
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`: second-order IIR notch filter
    - `:remez`: Remez FIR filter
- `ftype::Union{Symbol, Nothing}=nothing`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
- `n::Int64=2560`: signal length in samples
- `fs::Int64`: sampling rate
- `order::Int64=8`: filter order (6 dB/octave), number of taps for `:remez`, attenuation (× 4 dB) for `:fir` filters
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
- `bw::Real=-1`: bandwidth for `:iirnotch` and :remez filters
- `w::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter (default is Hamming window, number of taps is calculated using Fred Harris' rule-of-thumb) or weights for `:firls` filter
- `mono::Bool=false`: use color or gray palette
- `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_filter_response(; fs::Int64, n::Int64=2560, fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple}, order::Int64=8, rp::Real=8, rs::Real=-1, bw::Real=-1, w::Union{Vector{<:Real}, Nothing}=nothing, mono::Bool=false, frq_lim::Tuple{Real, Real}=(0, fs / 2), kwargs...)

    pal = mono ? :grays : :darktest
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    flt = filter_create(fprototype=fprototype, ftype=ftype, cutoff=cutoff, n=n, fs=fs, order=order, rp=rp, rs=rs, bw=bw, w=w)

    if fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch]
        H, w = freqresp(flt)
        # convert to dB
        H = 20 * log10.(abs.(H))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)

        if fprototype !== :iirnotch
            fname = titlecase(String(fprototype))
            title = "Filter: $(fname), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $order\n\nFrequency response"
        else
            fname = "IIR notch"
            title = "Filter: $(fname), cutoff: $cutoff Hz, band width: $bw\n\nFrequency response"
        end

        p1 = Plots.plot(w,
                        H,
                        title=title,
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylims=(-100, 10),
                        ylabel="Magnitude\n[dB]",
                        xlabel="Frequency [Hz]",
                        label="",
                        top_margin=10Plots.px,
                        bottom_margin=10Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            p1 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p1 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p1 = Plots.plot!((0, cutoff[2]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:green)
        end

        phi, w = phaseresp(flt)
        phi = rad2deg.(angle.(phi))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        p2 = Plots.plot(w,
                        phi,
                        title="Phase response",
                        ylims=(-180, 180),
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylabel="Phase\n[°]",
                        xlabel="Frequency [Hz]",
                        label="",
                        bottom_margin=10Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            p2 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p2 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p2 = Plots.plot!((0, cutoff[2]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:green)
        end

        tau, w = grpdelay(flt)
        tau = abs.(tau)
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        p3 = Plots.plot(w,
                        tau,
                        title="Group delay",
                        # xlims=(0, x_max),
                        xlims=frq_lim,
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
            p3 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p3 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p3 = Plots.plot!((0, cutoff[2]),
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
        H = 20 * log10.(abs.(H))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        if fprototype === :fir
            title = "Filter: FIR, type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, taps: $(length(flt)), attenuation: $(order * 15) dB\nFrequency response"
        elseif fprototype === :firls
            title = "Filter: FIR (LS), type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, order: $(length(flt))\nFrequency response"
        elseif fprototype === :remez
            title = "Filter: Remez, type: $(uppercase(String(ftype))), cutoff: $cutoff Hz, taps: $(length(flt))\nFrequency response"
        end
        p1 = Plots.plot(w,
                        H,
                        title=title,
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylims=(-100, 10),
                        ylabel="Magnitude\n[dB]",
                        xlabel="Frequency [Hz]",
                        label="",
                        bottom_margin=10Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            p1 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p1 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p1 = Plots.plot!((0, cutoff[2]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:green)
        end
        w = range(0, stop=pi, length=1024)
        phi = _fir_response(flt, w)
        phi = DSP.unwrap(-atan.(imag(phi), real(phi)))
        # convert rad/sample to Hz
        w = w .* fs / 2 / pi
        x_max = w[end]
        ftype === :hp && (x_max = cutoff * 10)
        p2 = Plots.plot(w,
                        phi,
                        title="Phase response",
                        # xlims=(0, x_max),
                        xlims=frq_lim,
                        ylabel="Phase\n[rad]",
                        xlabel="Frequency [Hz]",
                        label="",
                        titlefontsize=8,
                        xlabelfontsize=6,
                        ylabelfontsize=6,
                        xtickfontsize=6,
                        ytickfontsize=6;
                        palette=pal)
        if length(cutoff) == 1
            p2 = Plots.plot!((0, cutoff),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:black)
        else
            p2 = Plots.plot!((0, cutoff[1]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:red)
            p2 = Plots.plot!((0, cutoff[2]),
                             seriestype=:vline,
                             linestyle=:dash,
                             label="",
                             lw=0.5,
                             lc=:green)
        end

        p = Plots.plot(p1, p2, size=(1200, 800), margins=20*Plots.px, layout=(2, 1), palette=pal; kwargs...)
    end

    return p

end
