export plot_locs_nirs

"""
    plot_locs_nirs(locs, opt_pairs, src_n, det_n; <keyword arguments>)

Preview NIRS (Near-Infrared Spectroscopy) optodes and channels (source-detector pairs): sources (red), detectors (green), a line for each source-detector pair (dashed for short-separation pairs) and, optionally, the channels at the midpoints between their source and detector.

The optodes are looked up in `locs` by their labels (`src_names`, `det_names`); if the labels are not given, the first `src_n` rows of `locs` are the sources and the next `det_n` rows the detectors. Each source-detector pair is drawn once (the rows of `opt_pairs` repeat a pair for each wavelength and data type).

# Arguments

- `locs::DataFrame`: optode location data
- `opt_pairs::Matrix{Int64}`: matrix of source-detector pairs (n × 2: source index, detector index); rows with index 0 (channels without optodes) are ignored
- `src_n::Int64`: number of sources in the NIRS cap
- `det_n::Int64`: number of detectors in the NIRS cap
- `src_names::Union{Nothing, Vector{String}}=nothing`: source labels (in the order of the source indices)
- `det_names::Union{Nothing, Vector{String}}=nothing`: detector labels (in the order of the detector indices)
- `short::Union{Nothing, AbstractVector{Bool}}=nothing`: short-separation pairs (one value per row of `opt_pairs`), drawn with dashed lines
- `spairs::Vector{Int64}=Int64[]`: rows of `opt_pairs` (pairs) to highlight
- `ch_mid::Bool=false`: if `true`, plot channel markers at the midpoints of the pairs, numbered in the order of the pairs (as in `nirs_quality()`)
- `weights::Union{Nothing, AbstractVector{<:Real}}=nothing`: values shown at the channel midpoints (one value per pair, in the order of the pairs)
- `schematic::Bool=false`: if `true`, draw a schematic montage: sources as diamonds, detectors as circles, channels as numbered squares at the midpoints of the pairs (numbered in the order of the pairs, as in `nirs_quality()`), with a legend; highlighted channels (`spairs`) are dark, `weights` are shown above the channels, optode labels (`src_labels`, `det_labels`) are shown next to the optodes
- `src_labels::Bool=false`: plot source labels
- `det_labels::Bool=false`: plot detector labels
- `opt_labels::Bool=false`: plot optode type (S for source, D for detector) and number
- `head::Bool=true`: if `true`, draw head outline
- `head_labels::Bool=true`: if `true`, draw head labels
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
- `ch_info::Vector{String}=string.(1:DataFrames.nrow(locs))`: channel information details (currently not shown)

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_locs_nirs(
    locs::DataFrame,
    opt_pairs::Matrix{Int64},
    src_n::Int64,
    det_n::Int64;
    src_names::Union{Nothing, Vector{String}} = nothing,
    det_names::Union{Nothing, Vector{String}} = nothing,
    short::Union{Nothing, AbstractVector{Bool}} = nothing,
    spairs::Vector{Int64} = Int64[],
    ch_mid::Bool = false,
    weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    schematic::Bool = false,
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
    ch_info::Vector{String} = string.(1:DataFrames.nrow(locs)),
)::GLMakie.Figure
    # locations are stored in their original units: scale the montage to the unit sphere
    locs = _locs_unit(locs)

    # validate
    _check_var(ps, [:l, :m, :s], "ps")
    _check_var(plane, [:xy, :yz, :xz], "plane")
    src_n >= 0 && det_n >= 0 || throw(ArgumentError("src_n and det_n must be ≥ 0."))
    size(opt_pairs, 2) == 2 ||
        throw(ArgumentError("opt_pairs must have 2 columns (source, detector)."))
    all(0 .<= opt_pairs[:, 1] .<= src_n) ||
        throw(ArgumentError("Sources in opt_pairs must be in [1, $src_n]."))
    all(0 .<= opt_pairs[:, 2] .<= det_n) ||
        throw(ArgumentError("Detectors in opt_pairs must be in [1, $det_n]."))
    isnothing(short) || length(short) == size(opt_pairs, 1) ||
        throw(ArgumentError("short must contain one value per row of opt_pairs."))
    all(in(axes(opt_pairs, 1)), spairs) ||
        throw(ArgumentError("spairs must be rows of opt_pairs."))

    # rows of the optodes in locs
    if isnothing(src_names) && isnothing(det_names)
        src_n + det_n <= DataFrames.nrow(locs) || throw(
            ArgumentError("locs must contain ≥ $(src_n + det_n) optodes (src_n + det_n)."),
        )
        src_row = collect(1:src_n)
        det_row = collect((src_n + 1):(src_n + det_n))
        src_txt = locs[src_row, :label]
        det_txt = locs[det_row, :label]
    else
        (isnothing(src_names) || isnothing(det_names)) &&
            throw(ArgumentError("Both src_names and det_names must be given."))
        length(src_names) == src_n ||
            throw(ArgumentError("src_names must contain $src_n labels."))
        length(det_names) == det_n ||
            throw(ArgumentError("det_names must contain $det_n labels."))
        function _row(l)
            r = findfirst(==(l), locs[!, :label])
            isnothing(r) && throw(ArgumentError("Optode $l not found in locs."))
            return r
        end
        src_row = [_row(l) for l in src_names]
        det_row = [_row(l) for l in det_names]
        src_txt = src_names
        det_txt = det_names
    end

    # unique source-detector pairs (in the order of their first occurrence)
    pair_rows = Int64[]
    seen = Tuple{Int64, Int64}[]
    for j in axes(opt_pairs, 1)
        p = (opt_pairs[j, 1], opt_pairs[j, 2])
        (p[1] == 0 || p[2] == 0 || p in seen) && continue
        push!(seen, p)
        push!(pair_rows, j)
    end
    pair_n = length(pair_rows)
    # a pair is highlighted if any of its rows is
    hl = [p in Set((opt_pairs[j, 1], opt_pairs[j, 2]) for j in spairs) for p in seen]
    isnothing(weights) || length(weights) == pair_n ||
        throw(ArgumentError("weights must contain one value per pair ($pair_n)."))

    # projected optode positions
    ch = 1:DataFrames.nrow(locs)
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
    sx, sy = loc_x[src_row], loc_y[src_row]
    dx, dy = loc_x[det_row], loc_y[det_row]

    # axis limits
    xl = (-1.2, 1.2)
    yl = (-1.2, 1.2)

    opt_n = src_n + det_n
    if ps === :l
        plot_size   = (800, 800)
        marker_size = opt_n > 64 ? 10 : 20
        font_size   = 14
        sw          = 2
        lw          = 2
    elseif ps === :m
        plot_size   = (300, 300)
        marker_size = opt_n > 64 ? 5 : 10
        font_size   = 8
        sw          = 1
        lw          = 1.5
        src_labels  = false
        det_labels  = false
        grid        = false
    elseif ps === :s
        plot_size   = (100, 100)
        marker_size = opt_n > 64 ? 4 : 8
        font_size   = 8
        sw          = 0.5
        lw          = 1
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

    if schematic
        _nirs_schematic!(
            ax, sx, sy, dx, dy, opt_pairs, pair_rows;
            short   = short,
            hl      = hl,
            weights = weights,
            src_txt = src_labels ? src_txt : nothing,
            det_txt = det_labels ? det_txt : nothing,
            mono    = mono,
            ps      = ps,
        )
        head_labels && _draw_head_labels!(ax; font_size = font_size)
        return fig
    end

    # source-detector pairs (highlighted pairs in black / thicker, short pairs dashed)
    for (k, j) in enumerate(pair_rows)
        s, d = opt_pairs[j, 1], opt_pairs[j, 2]
        is_short = !isnothing(short) && short[j]
        GLMakie.lines!(
            ax,
            [sx[s], dx[d]],
            [sy[s], dy[d]];
            color     = hl[k] ? :black : (mono ? :gray : :blue),
            alpha     = hl[k] ? 1.0 : 0.5,
            linewidth = hl[k] ? 2 * lw : lw,
            linestyle = is_short ? :dash : :solid,
        )
    end

    label_offset_x = 0.0
    label_offset_y = -0.04

    # draw source markers or labels
    if src_labels
        for idx = 1:src_n
            GLMakie.text!(
                ax,
                sx[idx] + label_offset_x,
                sy[idx] + label_offset_y;
                text     = src_txt[idx],
                align    = (:center, :bottom),
                fontsize = font_size,
            )
        end
    elseif !opt_labels
        GLMakie.scatter!(
            ax,
            sx, sy;
            markersize  = marker_size,
            color       = mono ? :black : :red,
            strokewidth = sw,
            strokecolor = :black,
        )
    end

    # draw detector markers or labels
    if det_labels
        for idx = 1:det_n
            GLMakie.text!(
                ax,
                dx[idx] + label_offset_x,
                dy[idx] + label_offset_y;
                text     = det_txt[idx],
                align    = (:center, :bottom),
                fontsize = font_size,
            )
        end
    elseif !opt_labels
        GLMakie.scatter!(
            ax,
            dx, dy;
            markersize  = marker_size,
            color       = mono ? :white : :green,
            strokewidth = sw,
            strokecolor = :black,
        )
    end

    # draw S/D type labels
    if opt_labels
        for idx = 1:src_n
            GLMakie.text!(
                ax,
                sx[idx] + label_offset_x,
                sy[idx] + label_offset_y;
                text     = "S" * string(idx),
                align    = (:center, :bottom),
                fontsize = font_size,
            )
        end
        for idx = 1:det_n
            GLMakie.text!(
                ax,
                dx[idx] + label_offset_x,
                dy[idx] + label_offset_y;
                text     = "D" * string(idx),
                align    = (:center, :bottom),
                fontsize = font_size,
            )
        end
    end

    # channels: midpoints of the pairs
    if ch_mid || !isnothing(weights) || any(hl)
        mx = [(sx[opt_pairs[j, 1]] + dx[opt_pairs[j, 2]]) / 2 for j in pair_rows]
        my = [(sy[opt_pairs[j, 1]] + dy[opt_pairs[j, 2]]) / 2 for j in pair_rows]
        vis = ch_mid ? trues(pair_n) : copy(hl)
        GLMakie.scatter!(
            ax,
            mx[vis], my[vis];
            markersize  = marker_size / 2,
            color       = [hl[k] ? :black : :white for k in findall(vis)],
            strokewidth = sw / 2,
            strokecolor = :black,
        )
        if ch_mid && ps === :l
            for k = 1:pair_n
                GLMakie.text!(
                    ax,
                    mx[k],
                    my[k];
                    text     = string(k),
                    align    = (:center, :top),
                    offset   = (0, -marker_size / 4 - 2),
                    fontsize = font_size - 2,
                )
            end
        end
        if !isnothing(weights)
            for k = 1:pair_n
                GLMakie.text!(
                    ax,
                    mx[k],
                    my[k];
                    text     = _num_label(weights[k]),
                    fontsize = font_size,
                    color    = mono ? :black : (weights[k] >= 0 ? :red : :blue),
                    align    = (:center, :bottom),
                    offset   = (0, marker_size / 4 + 2),
                )
            end
        end
    end

    # draw head labels if requested
    head_labels && _draw_head_labels!(ax; font_size = font_size)

    return fig
end

# schematic NIRS montage: sources (diamonds), detectors (circles), channels (numbered squares at
# the midpoints of the pairs) connected by lines, legend
function _nirs_schematic!(
    ax::GLMakie.Axis,
    sx::AbstractVector,
    sy::AbstractVector,
    dx::AbstractVector,
    dy::AbstractVector,
    opt_pairs::Matrix{Int64},
    pair_rows::Vector{Int64};
    short::Union{Nothing, AbstractVector{Bool}},
    hl::Vector{Bool},
    weights::Union{Nothing, AbstractVector{<:Real}},
    src_txt::Union{Nothing, Vector{String}},
    det_txt::Union{Nothing, Vector{String}},
    mono::Bool,
    ps::Symbol,
)::Nothing
    pair_n = length(pair_rows)
    dense = pair_n > 30
    opt_size, ch_size, fs, show_num = if ps === :l
        dense ? (18, 20, 9, true) : (28, 30, 14, true)
    elseif ps === :m
        (10, 11, 6, !dense)
    else
        (5, 5, 4, false)
    end
    # the colors of the optodes as in the regular plot
    col_src = mono ? :black : :red
    col_det = mono ? :white : :green
    col_ch = :gray60
    col_hl = :gray15

    mx = [(sx[opt_pairs[j, 1]] + dx[opt_pairs[j, 2]]) / 2 for j in pair_rows]
    my = [(sy[opt_pairs[j, 1]] + dy[opt_pairs[j, 2]]) / 2 for j in pair_rows]

    # source - channel - detector lines
    for (k, j) in enumerate(pair_rows)
        s, d = opt_pairs[j, 1], opt_pairs[j, 2]
        GLMakie.lines!(
            ax,
            [sx[s], dx[d]],
            [sy[s], dy[d]];
            color     = :black,
            linewidth = 1,
            linestyle = (!isnothing(short) && short[j]) ? :dash : :solid,
        )
    end

    # optodes
    GLMakie.scatter!(
        ax, sx, sy;
        marker = :diamond, markersize = opt_size, color = col_src,
        strokewidth = 1, strokecolor = :black,
    )
    GLMakie.scatter!(
        ax, dx, dy;
        marker = :circle, markersize = opt_size, color = col_det,
        strokewidth = 1, strokecolor = :black,
    )
    for (x, y, txt) in ((sx, sy, src_txt), (dx, dy, det_txt))
        isnothing(txt) && continue
        for i in eachindex(txt)
            GLMakie.text!(
                ax, x[i], y[i];
                text = txt[i], fontsize = fs, font = :bold, align = (:left, :bottom),
                offset = (opt_size / 2, opt_size / 2),
            )
        end
    end

    # channels (short-separation channels, lying next to their optodes, are smaller)
    is_short = [!isnothing(short) && short[j] for j in pair_rows]
    GLMakie.scatter!(
        ax, mx, my;
        marker = :rect,
        markersize = [is_short[k] ? 0.6 * ch_size : ch_size for k = 1:pair_n],
        color = [hl[k] ? col_hl : col_ch for k = 1:pair_n],
        strokewidth = 1, strokecolor = :black,
    )
    if show_num
        for k = 1:pair_n
            GLMakie.text!(
                ax, mx[k], my[k];
                text = string(k), fontsize = is_short[k] ? 0.6 * fs : fs, font = :bold,
                color = hl[k] ? :white : :black, align = (:center, :center),
            )
        end
    end
    if !isnothing(weights)
        for k = 1:pair_n
            GLMakie.text!(
                ax, mx[k], my[k];
                text = _num_label(weights[k]), fontsize = fs,
                color = mono ? :black : (weights[k] >= 0 ? :red : :blue),
                align = (:center, :bottom), offset = (0, ch_size / 2 + 2),
            )
        end
    end

    # legend
    if ps === :l
        el = [
            GLMakie.MarkerElement(; marker = :diamond, color = col_src, markersize = 24,
                                  strokewidth = 1, strokecolor = :black),
            GLMakie.MarkerElement(; marker = :circle, color = col_det, markersize = 24,
                                  strokewidth = 1, strokecolor = :black),
            GLMakie.MarkerElement(; marker = :rect, color = col_ch, markersize = 24,
                                  strokewidth = 1, strokecolor = :black),
        ]
        lbl = ["Source", "Detector", "Channel"]
        if any(is_short)
            push!(el, GLMakie.MarkerElement(; marker = :rect, color = col_ch, markersize = 14,
                                            strokewidth = 1, strokecolor = :black))
            push!(lbl, "Short channel")
        end
        GLMakie.axislegend(
            ax, el, lbl;
            position = :rt, labelsize = 16, framevisible = true, framewidth = 2,
        )
    end

    return nothing
end
