export plot_connectivity_circle

"""
    plot_connectivity_circle(m; <keyword arguments>)

Plot connectivity circle.

# Arguments

  - `m::AbstractMatrix`: matrix of connectivities (channel vs. channel)
  - `clabels=Vector{String}`: channels labels
  - `title::String=""`: plot title
  - `threshold::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: if set, use threshold to mark a region
  - `threshold_type::Symbol=:neq`: rule for thresholding:
      + `:eq`: draw region is values are equal to threshold
      + `:neq`: draw region is values are not equal to threshold
      + `:geq`: draw region is values are ≥ to threshold
      + `:leq`: draw region is values are ≤ to threshold
      + `:g`: draw region is values are > to threshold
      + `:l`: draw region is values are < to threshold
      + `:in`: draw region is values are in the threshold values, including threshold boundaries
      + `:bin`: draw region is values are between the threshold values, excluding threshold boundaries

# Returns

  - `p::GLMakie.Figure`
"""
function plot_connectivity_circle(
        m::AbstractMatrix;
        clabels = Vector{String},
        title::String = "",
        threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
        threshold_type::Symbol = :neq,
    )::GLMakie.Figure

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

    pos_x = round.(pos_x, digits = 3)
    pos_y = round.(pos_y, digits = 3)

    # prepare plot
    GLMakie.activate!(title = "plot_connectivity_circle()")
    plot_size = (800, 800)
    p = GLMakie.Figure(size = plot_size, figure_padding = 0)
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = "",
        ylabel = "",
        title = title,
        aspect = 1,
        xticksvisible = false,
        yticksvisible = false,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
    )
    hidedecorations!(ax)
    GLMakie.xlims!(ax, (-1.5, 1.5))
    GLMakie.ylims!(ax, (-1.5, 1.5))
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # draw connections
    m_norm = normalize_minmax(m)
    c = (0.0, 0.0)
    s = size(m, 1)
    for idx1 in 1:s
        for idx2 in (idx1 + 1):s
            if !isnothing(threshold)
                if threshold_type in [:eq, :neq, :geq, :leq, :g, :l]
                    @assert length(threshold) == 1 "threshold must contain a single value."
                else
                    @assert length(threshold) == 2 "threshold must contain two values."
                    _check_tuple(threshold, "threshold")
                end
                (threshold_type === :eq && m[idx1, idx2] != threshold) && break
                (threshold_type === :neq && m[idx1, idx2] == threshold) && break
                (threshold_type === :g && m[idx1, idx2] <= threshold) && break
                (threshold_type === :l && m[idx1, idx2] >= threshold) && break
                (threshold_type === :geq && m[idx1, idx2] < threshold) && break
                (threshold_type === :leq && m[idx1, idx2] > threshold) && break
                (threshold_type === :in && (m[idx1, idx2] >= threshold[1] && m[idx1, idx2] <= threshold[2])) && break
                (threshold_type === :bin && (m[idx1, idx2] > threshold[1] && m[idx1, idx2] < threshold[2])) && break
            end
            mid_x = (pos_x[idx1] + pos_x[idx2]) / 2
            mid_y = (pos_y[idx1] + pos_y[idx2]) / 2
            d = distance((pos_x[idx1], pos_y[idx1]), (pos_x[idx2], pos_y[idx2]))
            c = (mid_x * (1 - d * 0.5), mid_y * (1 - d * 0.5))
            px = [pos_x[idx1], c[1], pos_x[idx2]]
            py = [pos_y[idx1], c[2], pos_y[idx2]]
            x_vals, y_vals = _bernstein_poly(px, py; steps = 50)
            col = :black
            m[idx1, idx2] < 0 && (col = :blue)
            m[idx1, idx2] > 0 && (col = :red)
            GLMakie.lines!(ax, x_vals, y_vals; color = col, linewidth = 10 * abs(m_norm[idx1, idx2]), alpha = 0.5)
        end
    end

    # draw markers
    for idx in axes(m, 1)
        GLMakie.scatter!(ax, pos_x[idx], pos_y[idx]; color = :black, markersize = 15)
    end

    # draw labels
    ang = t[1:(end - 1)]
    for idx in axes(clabels, 1)
        if _bin(ang[idx], (-pi / 2, pi / 2))
            GLMakie.text!(
                pos_x[idx] * 1.1,
                pos_y[idx] * 1.1;
                text = " " * clabels[idx],
                fontsize = 12,
                align = (:left, :center),
                rotation = ang[idx],
            )
        else
            GLMakie.text!(
                pos_x[idx] * 1.1,
                pos_y[idx] * 1.1;
                text = " " * clabels[idx],
                fontsize = 12,
                align = (:right, :center),
                rotation = (ang[idx] + pi),
            )
        end
    end

    return p

end
