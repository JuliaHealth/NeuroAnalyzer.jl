export plot_matrix
export plot_covmatrix
export plot_histogram
export plot_bar
export plot_line
export plot_box
export plot_violin
export plot_dots
export plot_paired
export plot_polar

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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_matrix(m::Array{<:Real, 2}; xlabels::Vector{String}, ylabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", cb_title::String="", mono::Bool=false, kwargs...)

    size(m, 1) == size(m, 2) || throw(ArgumentError("Matrix is not square."))
    length(xlabels) == length(ylabels) || throw(ArgumentError("Lengths of xlabels and ylabels must be equal."))
    length(xlabels) == size(m, 1) || throw(ArgumentError("Length of xlabels and matrix size must be equal."))
    length(ylabels) == size(m, 2) || throw(ArgumentError("Length of ylabels and matrix size must be equal."))

    n = size(m, 1)
    r = maximum(length.(xlabels)) > 10 ? 45 : 0
    mar = maximum(length.(xlabels)) > 10 ? 40 : 0
    pal = mono == true ? :grays : :darktest

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
                      left_margin=mar * Plots.px,
                      bottom_margin=mar * Plots.px,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      xtickfontsize=6,
                      ytickfontsize=6;
                      kwargs...)

    return p

end

"""
    plot_covmatrix(m, lags; <keyword arguments>)

Plot cross/auto-covariance matrix.

# Arguments

- `m::Abstractvector`: covariance matrix
- `lags::AbstractVector`: covariance lags, lags will be displayed in ms
- `xlabel::String="lag"`
- `ylabel::String=""`
- `title::String=""`
- `cb_title::String=""`: color bar title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_covmatrix(m::AbstractVector, lags::AbstractVector; xlabel::String="lag [ms]", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest
    r = length(lags) > 10 ? 90 : 0
    p = Plots.plot(lags,
                   m,
                   title=title,
                   xlabel=xlabel,
                   ylabel=ylabel,
                   xticks=round.(lags, digits=1),
                   xaxis=(tickfontrotation=r),
                   yticks=false,
                   palette=pal,
                   size=(600, 200),
                   lw=0.5,
                   grid=false,
                   legend=false,
                   bottom_margin=10 * Plots.px,
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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_histogram(s::AbstractVector; type::Symbol=:hist, bins::Union{Int64, Symbol, AbstractVector}=(length(s) ÷ 10), label::String="", xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    _check_var(type, [:hist, :kd], "type")

    type === :kd && (type = :density)

    pal = mono == true ? :grays : :darktest

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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_bar(s::AbstractVector; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    length(s) == length(xlabels) || throw(ArgumentError("signal length ($(length(s))) must be equal to xlabels length ($(length(xlabels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :lightblue

    p = Plots.plot(s,
                   seriestype=:bar,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(xlabels), xlabels),
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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_line(s::AbstractVector; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    length(s) == length(xlabels) || throw(ArgumentError("signal length ($(length(s))) must be equal to xlabels ($(length(xlabels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(s,
                   seriestype=:line,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(xlabels), xlabels),
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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_line(s::AbstractArray; rlabels::Vector{String}, xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    ndims(s) != 2 && throw(ArgumentError("signal must have 2-dimensions."))
    size(s, 1) == length(rlabels) || throw(ArgumentError("Number of signal columns ($(size(s, 1))) must be equal to labels length ($(length(rlabels)))."))
    size(s, 2) == length(xlabels) || throw(ArgumentError("Number of signal columns ($(size(s, 2))) must be equal to x-ticks length ($(length(xlabels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(s[1, :],
                   seriestype=:line,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=:topright,
                   label=rlabels[1],
                   xticks=(1:length(xlabels), xlabels),
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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_box(s::AbstractArray; glabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    ndims(s) != 2 && throw(ArgumentError("signal must have 2-dimensions."))
    size(s, 1) == length(glabels) || throw(ArgumentError("Number of signal columns ($(size(s, 1))) must be equal to x-ticks length ($(length(gxlabels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(s',
                   seriestype=:box,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(glabels), glabels),
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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_violin(s::AbstractArray; glabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    ndims(s) != 2 && throw(ArgumentError("signal must have 2-dimensions."))
    size(s, 1) == length(glabels) || throw(ArgumentError("Number of signal columns ($(size(s, 1))) must be equal to x-ticks length ($(length(gxlabels)))."))

    pal = mono == true ? :grays : :darktest
    color = mono == true ? :lightgrey : :auto

    p = Plots.plot(s',
                   seriestype=:violin,
                   size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(glabels), glabels),
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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_dots(signal::Vector{Vector{Float64}}; glabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    size(signal, 1) == length(glabels) || throw(ArgumentError("Number of signal columns ($(size(signal, 1))) must be equal to x-ticks length ($(length(xlabels)))."))

    pal = mono == true ? :grays : :darktest

    p = Plots.plot(size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(glabels), glabels),
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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_paired(signal::Vector{Vector{Float64}}; glabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)

    size(signal, 1) == length(glabels) || throw(ArgumentError("Number of signal columns ($(size(signal, 1))) must be equal to x-ticks length ($(length(xlabels)))."))
    ll = Vector{Int64}()
    for idx in eachindex(glabels)
        push!(ll, length(signal[idx]))
    end
    length(unique(ll)) == 1 || throw(ArgumentError("Each group must have the same number of values."))

    pal = mono == true ? :grays : :darktest

    p = Plots.plot(size=(1200, 500),
                   margins=20Plots.px,
                   legend=false,
                   xticks=(1:length(glabels), glabels),
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
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_polar(s::Union{AbstractVector, AbstractArray}; m::Tuple{Real, Real}=(0, 0), title::String="", mono::Bool=false, kwargs...)

    length(m) > 2 && throw(ArgumentError("m must have exactly 2 values: phases and lengths."))
    ndims(s) > 1 && size(s, 2) > 2 && throw(ArgumentError("signal must have exactly 2 columns: phases and lengths."))

    pal = mono == true ? :grays : :darktest

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
