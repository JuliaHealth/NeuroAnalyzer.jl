export plot_matrix
export plot_xac
export plot_histogram
export plot_bar
export plot_line
export plot_box
export plot_violin
export plot_dots
export plot_paired
export plot_polar
export plot_eros
export plot_erop
export plot_icatopo

"""
    plot_matrix(m; <keyword arguments>)

Plot matrix.

# Arguments

- `m::Array{<:Real, 2}`
- `xlabels::Vector{String}`
- `ylabels::Vector{String}`
- `xlabel::String=""`
- `ylabel::String=""`
- `title::String=""`
- `cb_title::String=""`: color bar title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_matrix(m::Array{<:Real, 2}; xlabels::Vector{String}, ylabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", cb_title::String="", mono::Bool=false, kwargs...)

    @assert size(m, 1) == size(m, 2) "Matrix is not square."
    @assert length(xlabels) == length(ylabels) "Lengths of xlabels and ylabels must be equal."
    @assert length(xlabels) == size(m, 1) "Length of xlabels and matrix size must be equal."
    @assert length(ylabels) == size(m, 2) "Length of ylabels and matrix size must be equal."

    n = size(m, 1)
    r = maximum(length.(xlabels)) > 10 ? 45 : 0
    mar = maximum(length.(xlabels)) > 10 ? 40 : 0
    pal = mono ? :grays : :darktest

    p = Plots.heatmap(m,
                      title=title,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      xaxis=(tickfontrotation=r),
                      xticks=(1:n, xlabels),
                      yticks=(1:n, xlabels),
                      seriescolor=pal,
                      colorbar_title=cb_title,
                      size=(1200, 800),
                      left_margin=mar*Plots.px,
                      bottom_margin=mar*Plots.px,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      xtickfontsize=6,
                      ytickfontsize=6;
                      kwargs...)

    return p

end

"""
    plot_xac(m, lags; <keyword arguments>)

Plot cross/auto-covariance/correlation.

# Arguments

- `m::Abstractvector`: covariance matrix
- `lags::AbstractVector`: covariance lags, lags will be displayed in s
- `xlabel::String="lag"`
- `ylabel::String=""`
- `title::String=""`
- `cb_title::String=""`: color bar title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_xac(m::AbstractVector, lags::AbstractVector; xlabel::String="lag [s]", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    pal = mono ? :grays : :darktest
    r = length(lags) > 10 ? 90 : 0
    p = Plots.plot(lags,
                   m,
                   title=title,
                   xlabel=xlabel,
                   ylabel=ylabel,
                   xticks=[lags[1], lags[end]],
                   xaxis=(tickfontrotation=r),
                   yticks=[round(minimum(m), digits=1), round(maximum(m), digits=1)],
                   palette=pal,
                   size=(600, 200),
                   lw=0.5,
                   grid=false,
                   legend=false,
                   bottom_margin=10*Plots.px,
                   titlefontsize=5,
                   xlabelfontsize=5,
                   ylabelfontsize=5,
                   xtickfontsize=4,
                   ytickfontsize=4;
                   kwargs...)

    return p

end

"""
    plot_histogram(s; <keyword arguments>)

Plot histogram.

# Arguments

- `s::AbstractVector`
- `type::Symbol`: type of histogram: regular (`:hist`) or kernel density (`:kd`)
- `bins::Union{Int64, Symbol, AbstractVector}=(length(s) ÷ 10)`: histogram bins: number of bins, range or `:sturges`, `:sqrt`, `:rice`, `:scott` or `:fd`)
- `label::String=""`: channel label
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_histogram(s::AbstractVector; type::Symbol=:hist, bins::Union{Int64, Symbol, AbstractVector}=(length(s) ÷ 10), label::String="", xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    _check_var(type, [:hist, :kd], "type")

    type === :kd && (type = :density)

    pal = mono ? :grays : :darktest

    if mean(s) < median(s)
        xticks = [floor(minimum(s), digits=1), round(mean(s), digits=1), round(median(s), digits=1), ceil(maximum(s), digits=1)]
    else
        xticks = [floor(minimum(s), digits=1), round(median(s), digits=1), round(mean(s), digits=1), ceil(maximum(s), digits=1)]
    end        

    p = Plots.plot(s,
                   seriestype=type,
                   xlabel=xlabel,
                   ylabel=ylabel,
                   label=label,
                   title=title,
                   palette=pal,
                   bins=bins,
                   grid=false,
                   linecolor=:black,
                   fillcolor=:grey,
                   fillalpha=0.5,
                   linewidth=1,
                   size=(800, 800),
                   left_margin=30*Plots.px,
                   xaxis=(tickfontrotation=90),
                   margins=10Plots.px,
                   xticks=xticks,
                   yticks=false,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=5,
                   ytickfontsize=5;
                   kwargs...)
    p = Plots.vline!([round(mean(s), digits=1)], lw=1, ls=:dot, lc=:black, label="mean")
    p = Plots.vline!([round(median(s), digits=1)], lw=0.5, ls=:dash, lc=:grey, alpha=0.5, label="median")

    Plots.plot(p)

    return p

end

"""
    plot_bar(s; <keyword arguments>)

Bar plot.

# Arguments

- `s::AbstractVector`
- `xlabels::Vector{String}`: x-ticks labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_bar(s::AbstractVector; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    @assert length(s) == length(xlabels) "signal length ($(length(s))) must be equal to xlabels length ($(length(xlabels)))."

    pal = mono ? :grays : :darktest
    color = mono == true ? :lightgrey : :lightblue

    p = Plots.plot(s,
                   seriestype=:bar,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(eachindex(xlabels), xlabels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)

    Plots.plot(p)

    return p

end

"""
    plot_line(s; <keyword arguments>)

Line plot.

# Arguments

- `s::AbstractVector`
- `xlabels::Vector{String}`: x-ticks labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_line(s::AbstractVector; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    @assert length(s) == length(xlabels) "signal length ($(length(s))) must be equal to xlabels ($(length(xlabels)))."

    pal = mono ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(s,
                   seriestype=:line,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(eachindex(xlabels), xlabels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)

    Plots.plot(p)

    return p

end

"""
    plot_line(s; <keyword arguments>)

Line plot.

# Arguments

- `s::AbstractArray`
- `rlabels::Vector{String}`: signal rows labels
- `xlabels::Vector{String}`: x-ticks labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_line(s::AbstractArray; rlabels::Vector{String}, xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    @assert ndims(s) == 2 "signal must have 2-dimensions."
    @assert size(s, 1) == length(rlabels) "Number of signal columns ($(size(s, 1))) must be equal to labels length ($(length(rlabels)))."
    @assert size(s, 2) == length(xlabels) "Number of signal columns ($(size(s, 2))) must be equal to x-ticks length ($(length(xlabels)))."

    pal = mono ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(s[1, :],
                   seriestype=:line,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=:topright,
                   label=rlabels[1],
                   xticks=(eachindex(xlabels), xlabels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)
    for idx in 2:size(s, 1)
        p = Plots.plot!(s[idx, :],
                        seriestype=:line,
                        color=idx,
                        label=rlabels[idx])
    end
    Plots.plot(p)

    return p

end

"""
    plot_box(s; <keyword arguments>)

Box plot.

# Arguments

- `s::AbstractArray`
- `glabels::Vector{String}`: group labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_box(s::AbstractArray; glabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    @assert ndims(s) == 2 "signal must have 2-dimensions."
    @assert size(s, 1) == length(glabels) "Number of signal columns ($(size(s, 1))) must be equal to x-ticks length ($(length(gxlabels)))."

    pal = mono ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(s',
                   seriestype=:box,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(eachindex(glabels), glabels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)
    Plots.plot(p)

    return p

end

"""
    plot_violin(s; <keyword arguments>)

Violin plot.

# Arguments

- `s::AbstractArray`
- `glabels::Vector{String}`: group labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_violin(s::AbstractArray; glabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    @assert ndims(s) == 2 "signal must have 2-dimensions."
    @assert size(s, 1) == length(glabels) "Number of signal columns ($(size(s, 1))) must be equal to x-ticks length ($(length(gxlabels)))."

    pal = mono ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(s',
                   seriestype=:violin,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(eachindex(glabels), glabels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   color=color,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)
    Plots.plot(p)

    return p

end

"""
    plot_dots(s; <keyword arguments>)

Dots plot.

# Arguments

- `s::Vector{Vector{Float64}}`
- `glabels::Vector{String}`: group labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_dots(signal::Vector{Vector{Float64}}; glabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    @assert size(signal, 1) == length(glabels) "Number of signal columns ($(size(signal, 1))) must be equal to x-ticks length ($(length(xlabels)))."

    pal = mono ? :grays : :darktest

    p = Plots.plot(size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(eachindex(glabels), glabels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs...)
    for idx1 in eachindex(glabels)
        for idx2 in eachindex(signal[idx1])
            if mono == false
                p = Plots.scatter!((idx1, signal[idx1][idx2]),
                                   color=idx1)
            else
                p = Plots.scatter!((idx1, signal[idx1][idx2]),
                                   color=:black)
            end
        end
    end

    Plots.plot(p)

    return p

end

"""
    plot_paired(signal; <keyword arguments>)

Plot paired data.

# Arguments

- `signal::Vector{Vector{Float64}}`
- `glabels::Vector{String}`: group labels
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_paired(signal::Vector{Vector{Float64}}; glabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    @assert size(signal, 1) == length(glabels) "Number of signal columns ($(size(signal, 1))) must be equal to x-ticks length ($(length(xlabels)))."
    ll = Vector{Int64}()
    for idx in eachindex(glabels)
        push!(ll, length(signal[idx]))
    end
    @assert length(unique(ll)) == 1 "Each group must have the same number of values."

    pal = mono ? :grays : :darktest

    p = Plots.plot(size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(eachindex(glabels), glabels),
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   palette=pal,
                   linewidth=0.5,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs...)
    for idx1 in eachindex(signal[1])
        c_tmp = zeros(length(glabels))
        for idx2 in eachindex(glabels)
            c_tmp[idx2] = signal[idx2][idx1]
        end
        p = Plots.plot!(c_tmp,
                        color=:black)
    end
    for idx1 in eachindex(glabels)
        for idx2 in eachindex(signal[idx1])
            if mono == false
                p = Plots.scatter!((idx1, signal[idx1][idx2]),
                                   color=idx1)
            else
                p = Plots.scatter!((idx1, signal[idx1][idx2]),
                                   color=:black)
            end
        end
    end

    Plots.plot(p)

    return p

end

"""
    plot_polar(s; <keyword arguments>)

Polar plot.

# Arguments

- `s::Union{AbstractVector, AbstractArray}`
- `m::Tuple{Real, Real}=(0, 0)`: major value to plot
- `title::String=""`: plot title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_polar(s::Union{AbstractVector, AbstractArray}; m::Tuple{Real, Real}=(0, 0), title::String="", mono::Bool=false, kwargs...)

    @assert length(m) == 2 "m must have exactly 2 values: phases and lengths."
    ndims(s) > 1 && @assert size(s, 2) == 2 "signal must have exactly 2 columns: phases and lengths."

    pal = mono ? :grays : :darktest

    if ndims(s) == 1
        p = Plots.plot([0, s[1]], [0, 1],
                       size=(800, 800),
                       projection=:polar,
                       left_margin=50Plots.px,
                       right_margin=50Plots.px,
                       bottom_margin=30Plots.px,
                       legend=false,
                       xticks=false,
                       yticks=false, 
                       title=title,
                       color=:black,
                       palette=pal,
                       linewidth=0.2,
                       titlefontsize=8,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       kwargs...)
        for idx in 2:length(s)
            p = Plots.plot!([0, s[idx]], [0, 1],
                            projection=:polar,
                            color=:black)
        end
    else
        p = Plots.plot([0, s[1, 1]], [0, s[1, 2]],
                       size=(800, 800),
                       projection=:polar,
                       left_margin=50Plots.px,
                       right_margin=50Plots.px,
                       bottom_margin=30Plots.px,
                       legend=false,
                       xticks=false,
                       yticks=false, 
                       title=title,
                       color=:black,
                       palette=pal,
                       linewidth=0.2,
                       titlefontsize=8,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       kwargs...)
        for idx in 2:size(s, 1)
            p = Plots.plot!([0, s[idx, 1]], [0, s[idx, 2]],
                            projection=:polar,
                            color=:black)
        end
    end
    if m != (0, 0)
        p = Plots.plot!([0, m[1]], [0, m[2]],
                        lw=2,
                        projection=:polar,
                        color=:red)
    end

    Plots.plot(p)

    return p

end

"""
    plot_eros(s, f, t; <keyword arguments>)

Plot ERO (Event-Related Oscillations) spectrogram.

# Arguments

- `s::AbstractArray`: ERO spectrogram
- `f::AbstractVector`: ERO frequencies
- `t::AbstractVector`: ERO time
- `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points
- `xlabel::String="default"`
- `ylabel::String="default"`
- `title::String="default"`
- `cb::Bool=true`: draw color bar
- `cb_title::String="Power [dB]"`: color bar title
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_eros(s::AbstractArray, f::AbstractVector, t::AbstractVector; tm::Union{Int64, Vector{Int64}}=0, xlabel::String="default", ylabel::String="default", title::String="default", cb::Bool=true, cb_title::String="Power [dB]", mono::Bool=false, kwargs...)

    @assert size(s, 1) == length(f) "f vector length does not match spectrogram."
    @assert size(s, 2) == length(t) "t vector length does not match spectrogram."
    @assert ndims(s) == 3 "s must have 3 dimensions."
    @assert size(s, 3) <= 2 "s must contain ≤ 2 epochs."

    pal = mono ? :grays : :darktest

    # set time markers
    if tm != 0
        for tm_idx in eachindex(tm)
            @assert tm[tm_idx] / 1000 >= t[1] "tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."
            @assert tm[tm_idx] / 1000 <= t[end] "tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."
            tm[tm_idx] = vsearch(tm[tm_idx] / 1000, t)
        end
    end

    t = t .* 1000

    if size(s, 3) == 1
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of epochs")
        p = Plots.heatmap(t,
                          f,
                          s[:, :, 1],
                          title=tt,
                          xlabel=xl,
                          ylabel=yl,
                          seriescolor=pal,
                          colorbar_title=cb_title,
                          size=(1200, 800),
                          left_margin=20*Plots.px,
                          bottom_margin=20*Plots.px,
                          titlefontsize=8,
                          xlabelfontsize=8,
                          ylabelfontsize=8,
                          xtickfontsize=6,
                          ytickfontsize=6;
                          kwargs...)

        # draw time markers
        if tm != 0
            for tm_idx in eachindex(tm)
                p = Plots.vline!([t[tm[tm_idx]]],
                                 linewidth=1,
                                 linecolor=:black,
                                 label=false)
            end
        end
    else
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Frequency [Hz]", "ERP spectrogram")
        p1 = Plots.heatmap(t,
                           f,                           
                           s[:, :, 1],
                           title=tt,
                           xlabel=xl,
                           ylabel=yl,
                           seriescolor=pal,
                           colorbar_title=cb_title,
                           size=(1200, 800),
                           left_margin=20*Plots.px,
                           bottom_margin=20*Plots.px,
                           titlefontsize=8,
                           xlabelfontsize=8,
                           ylabelfontsize=8,
                           xtickfontsize=6,
                           ytickfontsize=6;
                           kwargs...)
        # plot 0 v-line
        p1 = Plots.vline!([0],
                          linestyle=:dash,
                          linewidth=0.5,
                          linecolor=:black,
                          label=false)

        # draw time markers
        if tm != 0
            for tm_idx in eachindex(tm)
                p1 = Plots.vline!([t[tm[tm_idx]]],
                                  linewidth=1,
                                  linecolor=:black,
                                  label=false)
            end
        end

        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of ERP epochs")
        p2 = Plots.heatmap(t,
                           f,                           
                           s[:, :, 2],
                           title=tt,
                           xlabel=xl,
                           ylabel=yl,
                           seriescolor=pal,
                           colorbar_title=cb_title,
                           size=(1200, 800),
                           left_margin=20*Plots.px,
                           bottom_margin=20*Plots.px,
                           titlefontsize=8,
                           xlabelfontsize=8,
                           ylabelfontsize=8,
                           xtickfontsize=6,
                           ytickfontsize=6;
                           kwargs...)
        # plot 0 v-line
        p2 = Plots.vline!([0],
                          linestyle=:dash,
                          linewidth=0.5,
                          linecolor=:black,
                          label=false)

        # draw time markers
        if tm != 0
            for tm_idx in eachindex(tm)
                p2 = Plots.vline!([t[tm[tm_idx]]],
                                  linewidth=1,
                                  linecolor=:black,
                                  label=false)
            end
        end

        p = Plots.plot(p1, p2, layout=(2, 1))

    end

    return p

end

"""
    plot_erop(p, f; <keyword arguments>)

Plot ERO (Event-Related Oscillations) power-spectrum.

# Arguments

- `p::AbstractArray`: ERO powers
- `f::AbstractVector`: ERO frequencies
- `xlabel::String="default"`
- `ylabel::String="default"`
- `title::String="default"`
- `mono::Bool=false`: Use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erop(p::AbstractArray, f::AbstractVector; xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, kwargs...)

    @assert size(p, 1) == length(f) "f vector length does not match powers."
    @assert ndims(p) == 2 "p must have 2 dimensions."
    @assert size(p, 2) <= 2 "p must contain ≤ 2 epochs."

    pal = mono ? :grays : :darktest

    if size(p, 2) == 1
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", "Power [dB]", "Averaged power-spectra of epochs")
        p = Plots.plot(f,
                       p[:, 1],
                       title=tt,
                       xlabel=xl,
                       ylabel=yl,
                       seriescolor=pal,
                       size=(1200, 800),
                       left_margin=20*Plots.px,
                       bottom_margin=20*Plots.px,
                       titlefontsize=8,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6,
                       label=false;
                       kwargs...)
    else
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", "Power [dB]", "ERP power-spectrum")
        p1 = Plots.plot(f,
                        p[:, 1],
                        title=tt,
                        xlabel=xl,
                        ylabel=yl,
                        seriescolor=pal,
                        size=(1200, 800),
                        left_margin=20*Plots.px,
                        bottom_margin=20*Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=8,
                        ylabelfontsize=8,
                        xtickfontsize=6,
                        ytickfontsize=6,
                        label=false;
                        kwargs...)

        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", "Power [dB]", "Averaged power-spectra of epochs")
        p2 = Plots.plot(f,
                        p[:, 2],
                        title=tt,
                        xlabel=xl,
                        ylabel=yl,
                        seriescolor=pal,
                        size=(1200, 800),
                        left_margin=20*Plots.px,
                        bottom_margin=20*Plots.px,
                        titlefontsize=8,
                        xlabelfontsize=8,
                        ylabelfontsize=8,
                        xtickfontsize=6,
                        ytickfontsize=6,
                        label=false;
                        kwargs...)

        p = Plots.plot(p1, p2, layout=(2, 1))

    end

    return p

end

"""
    plot_icatopo(obj; <keyword arguments>)

Topographical plot of embedded ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `cb::Bool=false`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_icatopo(obj::NeuroAnalyzer.NEURO; ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), cb::Bool=false, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, kwargs...)

    @assert :ic in keys(obj.components) "OBJ does not contain :ic component. Perform ica_decompose() first."
    @assert :ic_mw in keys(obj.components) "OBJ does not contain :ic_mw component. Perform ica_decompose() first."
    
    ic = obj.components[:ic]
    ic_mw = obj.components[:ic_mw]

    # select component channels, default is all channels
    ic_idx == 0 && (ic_idx = _select_cidx(ic, ic_idx))
    _check_cidx(ic, ic_idx)

    p_topo = Vector{Plots.Plot{Plots.GRBackend}}()

    for idx in eachindex(ic_idx)
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=signal_channels(obj), ic_idx=ic_idx[idx], keep=true)
        p_tmp = plot_topo(obj_tmp, title="IC $(ic_idx[idx])", cb=cb, cb_label=cb_label, amethod=amethod, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, seg=seg, kwargs...)
        push!(p_topo, p_tmp)
    end
    
    if length(ic_idx) <= 4
        p = plot_compose(p_topo, layout=(1, 4))
    elseif length(ic_idx) <= 8
        p = plot_compose(p_topo, layout=(2, ceil(Int64, length(ic_idx) / 2)))
    elseif length(ic_idx) <= 12
        p = plot_compose(p_topo, layout=(3, ceil(Int64, length(ic_idx) / 3)))
    else
        p = plot_compose(p_topo, layout=(4, ceil(Int64, length(ic_idx) / 4)))
    end

    return p

end

"""
    plot_icatopo(obj, ic, ic_mw; <keyword arguments>)

Topographical plot of external ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `cb::Bool=false`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_icatopo(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), cb::Bool=false, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, kwargs...)

    # select component channels, default is all channels
    ic_idx == 0 && (ic_idx = _select_cidx(ic, ic_idx))
    _check_cidx(ic, ic_idx)

    p_topo = Vector{Plots.Plot{Plots.GRBackend}}()

    for idx in eachindex(ic_idx)
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=signal_channels(obj), ic_idx=ic_idx[idx], keep=true)
        p_tmp = plot_topo(obj_tmp, title="IC $(ic_idx[idx])", cb=cb, cb_label=cb_label, amethod=amethod, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, seg=seg, kwargs...)
        push!(p_topo, p_tmp)
    end

    if length(ic_idx) <= 4
        p = plot_compose(p_topo, layout=(1, 4))
    elseif length(ic_idx) <= 8
        p = plot_compose(p_topo, layout=(2, ceil(Int64, length(ic_idx) / 2)))
    elseif length(ic_idx) <= 12
        p = plot_compose(p_topo, layout=(3, ceil(Int64, length(ic_idx) / 3)))
    else
        p = plot_compose(p_topo, layout=(4, ceil(Int64, length(ic_idx) / 4)))
    end

    return p

end