export plot_locs_nirs

"""
    plot_locs_nirs(locs; <keyword arguments>)

Preview NIRS (Near-Infrared Spectroscopy) optodes and channel locations with customizable visualization.

# Arguments

- `locs::DataFrame`: channel location data
- `opt_pairs::Matrix{Int64}`: matrix of source-detector pairs (n_pairs × 2)
- `src_n::Int64`: number of sources in the NIRS cap
- `det_n::Int64`: number of detectors in the NIRS cap
- `src_labels::Bool=false`: plot source labels
- `det_labels::Bool=false`: plot detector labels
- `opt_labels::Bool=false`: plot optode type (S for source, D for detector) and number
- `head::Bool=true`: if `true`, draw head outline
- `head_labels::Bool=false`: if `true`, draw head labels
- `mono::Bool=false`: if `true`, use a monochrome palette
- `grid::Bool=false`: if `true`, draw grid (useful for locating positions)
- `ps::Symbol`: plot size:
    - `:l`: large (800×800 px)
    - `:m`: medium (300×300 px)
    - `:s`: small (100×100 px)
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `plane::Symbol=:xy`: which anatomical plane to plot:
    - `:xy`: horizontal (top-down) view
    - `:xz`: coronary (front) view
    - `:yz`: sagittal (side) view
- `ch_info::Vector{String}=string.(1:DataFrames.nrow(locs))`: channel information details

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_locs_nirs(
    locs::DataFrame,
    opt_pairs::Matrix{Int64},
    src_n::Int64,
    det_n::Int64;
    src_labels::Bool = false,
    det_labels::Bool = false,
    opt_labels::Bool = false,
    head::Bool = true,
    head_labels::Bool = true,
    mono::Bool = false,
    grid::Bool = false,
    ps::Symbol = :l,
    cart::Bool = false,
    plane::Symbol = :xy,
    ch_info::Vector{String}=string.(1:DataFrames.nrow(locs))
)::GLMakie.Figure

    # TO DO: plot channel numbers

    # channels
    ch = 1:DataFrames.nrow(locs)

    # validate
    _check_var(ps, [:l, :m, :s], "ps")
    _check_var(plane, [:xy, :yz, :xz], "plane")

    # set color palette
    pal = mono ? :grays : :darktest

    loc_x = zeros(length(ch))
    loc_y = zeros(length(ch))

    if plane === :xy
        if cart
            loc_x = locs[ch, :loc_x]
            loc_y = locs[ch, :loc_y]
        else
            for idx in eachindex(ch)
                loc_x[idx], loc_y[idx] =
                    pol2cart(locs[ch, :loc_radius][idx], locs[ch, :loc_theta][idx])
            end
        end
    elseif plane === :xz
        if cart
            loc_x = locs[ch, :loc_x]
            loc_y = locs[ch, :loc_z]
        else
            for idx in eachindex(ch)
                loc_x[idx], _, loc_y[idx] = sph2cart(
                    locs[ch, :loc_radius_sph][idx],
                    locs[ch, :loc_theta_sph][idx],
                    locs[ch, :loc_phi_sph][idx],
                )
            end
        end
    elseif plane === :yz
        if cart
            loc_x = locs[ch, :loc_y]
            loc_y = locs[ch, :loc_z]
        else
            for idx in eachindex(ch)
                _, loc_x[idx], loc_y[idx] = sph2cart(
                    locs[ch, :loc_radius_sph][idx],
                    locs[ch, :loc_theta_sph][idx],
                    locs[ch, :loc_phi_sph][idx],
                )
            end
        end
    end

    # axis limits
    xl = (-1.2, 1.2)
    yl = (-1.2, 1.2)

    if ps === :l
        plot_size   = (800, 800)
        marker_size = length(ch) > 64 ? 10 : 20
        font_size   = 14
        sw          = 2
    elseif ps === :m
        plot_size   = (300, 300)
        marker_size = length(ch) > 64 ? 5 : 10
        font_size   = 8
        sw          = 1
        src_labels  = false
        det_labels  = false
        grid        = false
    elseif ps === :s
        plot_size   = (100, 100)
        marker_size = length(ch) > 64 ? 4 : 8
        font_size   = 8
        sw          = 0.5
        head_labels = false
        src_labels  = false
        det_labels  = false
        grid        = false
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_locs_nirs()")
    fig = GLMakie.Figure(; size = plot_size, figure_padding = 0)

    shared_axis_kwargs = (
        aspect = 1,
        xlabel = "",
        ylabel = "",
        title = "",
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        backgroundcolor = :transparent,
        _AXIS_LOCK_KWARGS...,
    )

    if grid
        xt = range(xl[1], xl[2]; step = 0.1)
        yt = range(yl[1], yl[2]; step = 0.1)
        ax = GLMakie.Axis(
            fig[1, 1];
            shared_axis_kwargs...,
            xticks             = xt,
            xminorticksvisible = true,
            xminorticks        = IntervalsBetween(2),
            yticks             = yt,
            yminorticksvisible = true,
            yminorticks        = IntervalsBetween(2),
        )
    else
        ax = GLMakie.Axis(fig[1, 1]; shared_axis_kwargs...)
        hidedecorations!(ax; grid = true)
        hidespines!(ax)
    end
    GLMakie.xlims!(ax, xl)
    GLMakie.ylims!(ax, yl)
    _style_axis!(ax)

    # draw head outline
    head && _draw_head_outline!(ax; lw = 3)

    # number of channels
    ch_n = length(ch)

    cmap = GLMakie.resample_cmap(pal, ch_n)

    # draw channel connection lines
    for idx in axes(opt_pairs, 1)
        xs = loc_x[opt_pairs[idx, 1]]
        xd = loc_x[src_n + opt_pairs[idx, 2]]
        ys = loc_y[opt_pairs[idx, 1]]
        yd = loc_y[src_n + opt_pairs[idx, 2]]
        GLMakie.lines!([xs, xd], [ys, yd]; color = mono ? :gray : :blue, alpha = 0.5)
    end

    label_offset_x = 0.0
    label_offset_y = -0.08

    # draw source markers or labels
    if src_labels
        for idx in 1:src_n
            GLMakie.text!(
                loc_x[idx] + label_offset_x,
                loc_y[idx] + label_offset_y;
                text     = locs[!, :label][idx],
                align    = (:center, :bottom),
                fontsize = font_size,
            )
        end
    elseif !opt_labels
        GLMakie.scatter!(
            loc_x[1:src_n], loc_y[1:src_n];
            markersize  = marker_size,
            color       = mono ? :black : :red,
            strokewidth = sw,
            strokecolor = :black,
        )
    end

    # draw detector markers or labels
    if det_labels
        for idx in (src_n + 1):(src_n + det_n)
            GLMakie.text!(
                loc_x[idx] + label_offset_x,
                loc_y[idx] + label_offset_y;
                text     = locs[!, :label][idx],
                align    = (:center, :bottom),
                fontsize = font_size,
            )
        end
    elseif !opt_labels
        GLMakie.scatter!(
            loc_x[(src_n + 1):end], loc_y[(src_n + 1):end];
            markersize  = marker_size,
            color       = mono ? :white : :green,
            strokewidth = sw,
            strokecolor = :black,
        )
    end

    # draw S/D type labels
    if opt_labels
        for idx in 1:src_n
            GLMakie.text!(
                loc_x[idx] + label_offset_x,
                loc_y[idx] + label_offset_y;
                text     = "S" * string(idx),
                align    = (:center, :bottom),
                fontsize = font_size,
            )
        end
        for idx in 1:det_n
            GLMakie.text!(
                loc_x[src_n + idx] + label_offset_x,
                loc_y[src_n + idx] + label_offset_y;
                text     = "D" * string(idx),
                align    = (:center, :bottom),
                fontsize = font_size,
            )
        end
    end

    # draw head labels if requested
    head_labels && _draw_head_labels!(ax; font_size = font_size)

    return fig
end
