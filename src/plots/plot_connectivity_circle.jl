export plot_connectivity_circle

"""
    plot_connectivity_circle(m; <keyword arguments>)

Plot connectivity circle.

# Arguments

- `m::AbstractMatrix`: matrix of connectivities (channel vs. channel)
- `clabels=Vector{String}`: channels labels
- `title::String=""`: plot title
- `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_connectivity_circle(m::AbstractMatrix; clabels=Vector{String}, title::String="", threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert size(m, 1) == length(clabels) "Number of channels in m ($(size(m, 1))) and clabels length ($(length(clabels))) differ."
    @assert size(m, 1) >= 2 "m must contain data for ≥ 2 channels."
    @assert size(m, 1) == size(m, 2) "m must be a square matrix."

    t = linspace(pi, -pi, size(m, 1) + 1)
    pos_x = zeros(size(m, 1))
    pos_y = zeros(size(m, 1))
    for idx in axes(m, 1)
        pos_x[idx] = cos(t[idx])
        pos_y[idx] = sin(t[idx])
    end

    pos_x = round.(pos_x, digits=3)
    pos_y = round.(pos_y, digits=3)

    p = Plots.plot(grid=false,
                   framestyle=:none,
                   border=:none,
                   legend=false,
                   axisratio=:equal,
                   size=(800, 800),
                   margins=-200Plots.px,
                   title=title,
                   xlims=(-1.5, 1.5),
                   ylims=(-1.5, 1.5),
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8;
                   kwargs=kwargs)

    # draw connections
    m_norm = NeuroAnalyzer.normalize_minmax(m)
    c = (0.0, 0.0)
    s = size(m, 1)
    for idx1 in 1:s
        for idx2 in (idx1 + 1):s
            if !isnothing(threshold)
                (threshold_type === :eq && m[idx1, idx2] != threshold) && break
                (threshold_type === :neq && m[idx1, idx2] == threshold) && break
                (threshold_type === :g && m[idx1, idx2] <= threshold) && break
                (threshold_type === :l && m[idx1, idx2] >= threshold) && break
                (threshold_type === :geq && m[idx1, idx2] < threshold) && break
                (threshold_type === :leq && m[idx1, idx2] > threshold) && break
            end
            mid_x = (pos_x[idx1] + pos_x[idx2]) / 2
            mid_y = (pos_y[idx1] + pos_y[idx2]) / 2
            d = distance((pos_x[idx1], pos_y[idx1]), (pos_x[idx2], pos_y[idx2]))
            c = (mid_x * (1 - d * 0.5), mid_y * (1 - d * 0.5))
            px = [pos_x[idx1], c[1], pos_x[idx2]]
            py = [pos_y[idx1], c[2], pos_y[idx2]]
            x_vals, y_vals = _bernstein_poly(px, py; steps=50)
            col = :black
            m[idx1, idx2] < 0 && (col = :blue)
            m[idx1, idx2] > 0 && (col = :red)
            Plots.plot!(x_vals,
                        y_vals,
                        color=col,
                        lw=10 * abs(m_norm[idx1, idx2]),
                        alpha=0.5,
                        label="")
        end
    end

    # draw markers
    for idx in axes(m, 1)
        p = Plots.scatter!((pos_x[idx], pos_y[idx]),
                           mc=:black,
                           ms=5)
    end

    ang = round.(rad2deg.(t[1:end-1]), digits=2)
    # draw labels
    for idx in axes(clabels, 1)
        if _bin(ang[idx], (-90, 90))
            p = Plots.plot!(annotations=(pos_x[idx] * 1.1, pos_y[idx] * 1.1, Plots.text(" " * clabels[idx], pointsize=8, rotation=ang[idx])))
        else
            p = Plots.plot!(annotations=(pos_x[idx] * 1.1, pos_y[idx] * 1.1, Plots.text(clabels[idx] * " ", pointsize=8, rotation=ang[idx] + 180)))
        end
    end

    Plots.plot(p)

    return p

end
