export plot_connectivity_circle

"""
    plot_connectivity_circle(m; <keyword arguments>)

Plot a circular connectivity diagram for a matrix of connectivities.

# Arguments

- `m::AbstractMatrix`: connectivity matrix (must be square, channel vs. channel)
- `clabels=Vector{String}`: channels labels (must match matrix dimensions)
- `title::String=""`: plot title
- `threshold::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: threshold for marking regions
    - if `Real`, use a single threshold value
    - if `Tuple{Real, Real}`, use a range for `:in` or `:bin` thresholding
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: values equal to threshold
    - `:neq`: values not equal to threshold
    - `:geq`: values ≥ threshold
    - `:leq`: values ≤ threshold
    - `:g`: values > threshold
    - `:l`: values < threshold
    - `:in`: values in the threshold range (inclusive)
    - `:bin`: values in the threshold range (exclusive)

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_connectivity_circle(
    m::AbstractMatrix;
    clabels = Vector{String},
    title::String = "",
    threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    threshold_type::Symbol = :neq,
)::GLMakie.Figure

    # validate
    size(m, 1) == length(clabels) ||
        throw(
            ArgumentError(
                "Number of channels in m ($(size(m, 1))) and clabels length ($(length(clabels))) must match.",
            ),
        )
    size(m, 1) >= 2 ||
        throw(ArgumentError("Connectivity matrix must contain data for ≥ 2 channels."))
    size(m, 1) == size(m, 2) || throw(ArgumentError("Connectivity matrix must be square."))

    # calculate polar coordinates for each channel
    t = range(π, -π; length = size(m, 1) + 1)
    pos_x = [cos(t[idx]) for idx in 1:size(m, 1)]
    pos_y = [sin(t[idx]) for idx in 1:size(m, 1)]

    # normalize connectivity matrix for visualization
    m_norm = normalize_minmax(m)

    # prepare plot
    GLMakie.activate!(; title = "plot_connectivity_circle()")
    plot_size = (800, 800)
    fig = GLMakie.Figure(; size = plot_size, figure_padding = 0)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = "",
        ylabel = "",
        title = title,
        aspect = 1,
        xticksvisible = false,
        yticksvisible = false,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    hidedecorations!(ax)
    _style_axis!(ax)

    # set axis limits
    GLMakie.xlims!(ax, (-1.5, 1.5))
    GLMakie.ylims!(ax, (-1.5, 1.5))

    # draw connections
    m_norm = normalize_minmax(m)
    c = (0.0, 0.0)
    s = size(m, 1)
    for idx1 in 1:s
        for idx2 in (idx1 + 1):s
            # apply thresholding if specified
            if !isnothing(threshold)
                if threshold_type in [:eq, :neq, :geq, :leq, :g, :l]
                    length(threshold) == 1 ||
                        throw(ArgumentError("threshold must contain a single value."))
                else
                    length(threshold) == 2 ||
                        throw(ArgumentError("threshold must contain two values."))
                    _check_tuple(threshold, extrema(m), "threshold")
                end
                (threshold_type === :eq && m[idx1, idx2] != threshold) && break
                (threshold_type === :neq && m[idx1, idx2] == threshold) && break
                (threshold_type === :g && m[idx1, idx2] <= threshold) && break
                (threshold_type === :l && m[idx1, idx2] >= threshold) && break
                (threshold_type === :geq && m[idx1, idx2] < threshold) && break
                (threshold_type === :leq && m[idx1, idx2] > threshold) && break
                (
                    threshold_type === :in &&
                    (m[idx1, idx2] >= threshold[1] && m[idx1, idx2] <= threshold[2])
                ) && break
                (
                    threshold_type === :bin &&
                    (m[idx1, idx2] > threshold[1] && m[idx1, idx2] < threshold[2])
                ) && break
            end

            # calculate midpoint and curvature
            mid_x = (pos_x[idx1] + pos_x[idx2]) / 2
            mid_y = (pos_y[idx1] + pos_y[idx2]) / 2
            d = sqrt((pos_x[idx2] - pos_x[idx1])^2 + (pos_y[idx2] - pos_y[idx1])^2)
            c = (mid_x * (1 - d * 0.5), mid_y * (1 - d * 0.5))

            # draw curved connection
            px = [pos_x[idx1], c[1], pos_x[idx2]]
            py = [pos_y[idx1], c[2], pos_y[idx2]]
            x_vals, y_vals = _bernstein_poly(px, py; steps = 50)

            col = :black
            m[idx1, idx2] < 0 && (col = :blue)
            m[idx1, idx2] > 0 && (col = :red)
            GLMakie.lines!(
                ax,
                x_vals,
                y_vals;
                color = col,
                linewidth = 10 * abs(m_norm[idx1, idx2]),
                alpha = 0.5,
            )
        end
    end

    # draw channel markers
    for idx in 1:s
        GLMakie.scatter!(
            ax,
            pos_x[idx],
            pos_y[idx];
            color = :black,
            markersize = 15,
        )
    end

    # draw channel labels
    ang = t[1:(end - 1)]
    for idx in axes(clabels, 1)
        if _bin(ang[idx], (-π / 2, π / 2))
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

    return fig
end
