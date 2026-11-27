export plot_locs

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Return true if `val` passes the given threshold rule.
function _passes_threshold(val::Real, threshold, threshold_type::Symbol)::Bool
    threshold_type === :g && return val > threshold
    threshold_type === :l && return val < threshold
    threshold_type === :eq && return val == threshold
    threshold_type === :neq && return val != threshold
    threshold_type === :leq && return val <= threshold
    threshold_type === :geq && return val >= threshold
    threshold_type === :in && return val >= threshold[1] && val <= threshold[2]
    threshold_type === :bin && return val > threshold[1] && val < threshold[2]
    return false
end

# Draw a single weighted connection line between two channel positions.
function _draw_connection!(
    loc_x::AbstractVector,
    loc_y::AbstractVector,
    idx1::Int,
    idx2::Int,
    val::Real,
    weight::Real,
    mono::Bool,
    use_weights::Bool,
)
    xs = [loc_x[idx1], loc_x[idx2]]
    ys = [loc_y[idx1], loc_y[idx2]]
    if use_weights
        lw = 6 * weight
        al = 0.25 * weight
        if val > 0
            GLMakie.lines!(
                xs,
                ys;
                linewidth = lw,
                alpha = al,
                color = mono ? :black : :red,
            )
        elseif val < 0
            GLMakie.lines!(
                xs,
                ys;
                linewidth = lw,
                alpha = al,
                color = mono ? :black : :blue,
                linestyle = mono ? :dot : :solid,
            )
        end
    else
        GLMakie.lines!(
            xs,
            ys;
            linewidth = 0.2,
            color = :black,
        )
    end
end

# Text of a weight / connection label (3 significant digits for non-integers).
_num_label(v::Real)::String = v isa Integer ? string(v) : string(round(v; sigdigits = 3))

# Draw a connection weight label at the midpoint between two channel positions.
function _draw_connection_label!(
    loc_x::AbstractVector,
    loc_y::AbstractVector,
    idx1::Int,
    idx2::Int,
    val::Real,
    font_size::Int,
    mono::Bool,
)
    l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
    color = mono ? :black : (val >= 0 ? :red : :blue)
    return GLMakie.text!(
        l_pos[1],
        l_pos[2];
        align    = (:center, :center),
        text     = _num_label(val),
        fontsize = font_size,
        color    = color,
    )
end

# ---------------------------------------------------------------------------

"""
    plot_locs(locs; <keyword arguments>)

Preview channel locations with customizable visualization and connection mapping.

# Arguments

- `locs::DataFrame`: channel location data
- `ch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=1:DataFrames.nrow(locs)`: list of locations to plot, default is all locations
- `sch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=0`: significant channels to highlight
- `ch_labels::Bool=true`: if `true`, draw locations labels
- `head::Bool=true`: if `true`, draw head outline
- `head_labels::Bool=false`: draw head labels
- `mono::Bool=false`: if `true`, use a monochrome palette
- `grid::Bool=true`: if `true`, draw grid for locating positions
- `ps::Symbol`: plot size:
    - `:l`: large (800×800 px)
    - `:m`: medium (300×300 px)
    - `:s`: small (100×100 px)
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `plane::Symbol=:xy`: which plane to plot:
    - `:xy`: horizontal (top)
    - `:xz`: coronary (front)
    - `:yz`: sagittal (side)
- `connections::Union{Nothing, Matrix{<:Real}}=nothing`: matrix of connections weights, shape (channels, channels)
- `threshold::Real=0`: threshold for plotting connections
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: values equal to threshold
    - `:neq`: values not equal to threshold
    - `:geq`: values ≥ threshold
    - `:leq`: values ≤ threshold
    - `:g`: values > threshold
    - `:l`: values < threshold
    - `:in`: values in the threshold range (inclusive)
    - `:bin`: values in the threshold range (exclusive)
- `weights::Union{Bool, Vector{<:Real}}=true`: if `true`, auto-scale line widths and transparency based on connection strength; if Vector, use provided weights to place at channel location
- `ch_info::Vector{String}=string.(1:DataFrames.nrow(locs))`: channel information details
- `ch_colors::Union{Nothing, AbstractVector}=nothing`: marker color of each channel in `ch` (overrides `sch` / `mono` coloring)
- `bad::Union{Nothing, AbstractVector{Bool}}=nothing`: bad channels (one value per row of `locs`), drawn with a red ×

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_locs(
    locs::DataFrame;
    ch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = 1:DataFrames.nrow(locs),
    sch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = 0,
    ch_labels::Bool = true,
    head::Bool = true,
    head_labels::Bool = false,
    mono::Bool = false,
    grid::Bool = false,
    ps::Symbol = :l,
    cart::Bool = false,
    plane::Symbol = :xy,
    connections::Union{Nothing, Matrix{<:Real}} = nothing,
    threshold::Real = 0,
    threshold_type::Symbol = :neq,
    weights::Union{Bool, Vector{<:Real}} = true,
    ch_info::Vector{String} = string.(1:DataFrames.nrow(locs)),
    ch_colors::Union{Nothing, AbstractVector} = nothing,
    bad::Union{Nothing, AbstractVector{Bool}} = nothing,
    gui::Bool = true,
)::GLMakie.Figure
    # locations are stored in their original units: scale the montage to the unit sphere
    locs = _locs_unit(locs)
    # validate
    _check_var(ps, [:l, :m, :s], "ps")
    isnothing(ch_colors) || length(ch_colors) == length(ch) ||
        throw(ArgumentError("Length of ch_colors must equal the number of channels."))
    _check_var(plane, [:xy, :yz, :xz], "plane")

    # set color palette
    pal = _palette(:darktest; mono = mono)

    # significant channel labels
    sch_labels = ch_labels

    loc_x = zeros(length(ch))
    loc_y = zeros(length(ch))

    if plane === :xy
        if cart
            loc_x = locs.loc_x[ch]
            loc_y = locs.loc_y[ch]
        else
            for idx in eachindex(ch)
                loc_x[idx], loc_y[idx] =
                    pol2cart(locs.loc_radius[ch][idx], locs.loc_theta[ch][idx])
            end
        end
    elseif plane === :xz
        if cart
            loc_x = locs.loc_x[ch]
            loc_y = locs.loc_z[ch]
        else
            for idx in eachindex(ch)
                loc_x[idx], _, loc_y[idx] = sph2cart(
                    locs.loc_radius_sph[ch][idx], locs.loc_theta_sph[ch][idx],
                    locs.loc_phi_sph[ch][idx],
                )
            end
        end
    elseif plane === :yz
        if cart
            loc_x = locs.loc_y[ch]
            loc_y = locs.loc_z[ch]
        else
            for idx in eachindex(ch)
                _, loc_x[idx], loc_y[idx] = sph2cart(
                    locs.loc_radius_sph[ch][idx], locs.loc_theta_sph[ch][idx],
                    locs.loc_phi_sph[ch][idx],
                )
            end
        end
    end

    loc_x = _n2v(loc_x)
    loc_y = _n2v(loc_y)

    head12 =
        maximum(abs.(locs.loc_x)) < 1.2 &&
        maximum(abs.(locs.loc_y)) < 1.2 &&
        maximum(abs.(locs.loc_z)) < 1.2

    # the side view of the head extends beyond the unit sphere (nose at y ≈ 1.21)
    head12 = head12 && !(head && plane === :yz)
    xl = head12 ? (-1.2, 1.2) : (-1.4, 1.4)
    yl = head12 ? (-1.2, 1.2) : (-1.4, 1.4)

    # plot parameters
    if ps === :l
        plot_size   = (800, 800)
        marker_size = length(ch) > 64 ? 6 : length(ch) > 32 ? 9 : 12
        font_size   = length(ch) > 64 ? 10 : 14
        lw          = 3
        sw          = 2
    elseif ps === :m
        plot_size   = (300, 300)
        marker_size = length(ch) > 64 ? 3 : length(ch) > 32 ? 4 : 5
        font_size   = 8
        lw          = 2
        sw          = 1
        ch_labels   = false
        sch_labels  = false
        grid        = false
    elseif ps === :s
        plot_size   = (100, 100)
        marker_size = length(ch) > 64 ? 2 : 3
        font_size   = 8
        lw          = 1
        sw          = 0.0
        head_labels = false
        ch_labels   = false
        sch_labels  = false
        grid        = false
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_locs()")
    fig = GLMakie.Figure(;
        size = plot_size,
        figure_padding = grid ? (10, 10, 10, 10) : (0, 0, 0, 0),
    ) # L R B T

    shared_ax_kwargs = (
        aspect = 1,
        xlabel = "",
        ylabel = "",
        title = "",
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        backgroundcolor = :transparent,
        _AXIS_LOCK_KWARGS...,
    )

    # create axis with customizable properties
    if grid
        ax = GLMakie.Axis(
            fig[1, 1];
            shared_ax_kwargs...,
            xminorticksvisible = true,
            xminorticks        = IntervalsBetween(5),
            yminorticksvisible = true,
            yminorticks        = IntervalsBetween(5),
        )
    else
        ax = GLMakie.Axis(fig[1, 1]; shared_ax_kwargs...)
        hidedecorations!(ax; grid = true)
        hidespines!(ax)
    end
    GLMakie.xlims!(ax, xl)
    GLMakie.ylims!(ax, yl)

    # draw head outline
    if head
        ps === :l && (lw = 3)
        ps === :m && (lw = 2)
        ps === :s && (lw = 1)
        _draw_head_outline!(ax; plane = plane, lw = lw)
    end

    # draw connections lines
    if !isnothing(connections)
        size(connections, 1) == length(ch) || throw(
            ArgumentError("Number of connections rows must equal number of channels."),
        )
        _check_var(
            threshold_type,
            [:eq, :neq, :geq, :leq, :g, :l, :in, :bin],
            "threshold_type",
        )
        if threshold_type in [:eq, :neq, :geq, :leq, :g, :l]
            length(threshold) == 1 ||
                throw(ArgumentError("threshold must contain a single value."))
        else
            length(threshold) == 2 ||
                throw(ArgumentError("threshold must contain two values."))
            _check_tuple(threshold, extrema(connections), "threshold")
        end

        m_tmp = normalize_n(abs.(connections))
        use_weights = weights === true

        for idx1 in axes(connections, 1)
            for idx2 = (idx1 + 1):size(connections, 1)
                val = connections[idx1, idx2]
                if _passes_threshold(val, threshold, threshold_type)
                    _draw_connection!(
                        loc_x,
                        loc_y,
                        idx1,
                        idx2,
                        val,
                        m_tmp[idx1, idx2],
                        mono,
                        use_weights,
                    )
                end
            end
        end
    end

    # draw channel markers
    ch_n = length(ch)
    # (resample_cmap needs ≥ 2 colors)
    cmap = _trace_colors(ch_n; mono = mono)
    sch_set = Set(sch)

    for (i, idx) in enumerate(ch)
        if !isnothing(ch_colors)
            GLMakie.scatter!(
                loc_x[i],
                loc_y[i];
                markersize  = marker_size,
                color       = ch_colors[i],
                strokewidth = sw,
                strokecolor = :black,
            )
        elseif idx in sch_set
            GLMakie.scatter!(
                loc_x[i],
                loc_y[i];
                markersize  = marker_size,
                color       = mono ? :gray : cmap[i],
                colormap    = pal,
                colorrange  = 1:ch_n,
                strokewidth = sw,
                strokecolor = :black,
            )
        else
            GLMakie.scatter!(
                loc_x[i],
                loc_y[i];
                markersize  = marker_size,
                color       = :gray,
                strokewidth = sw,
                strokecolor = :black,
            )
        end
        # bad channel: red × over the marker
        if !isnothing(bad) && bad[idx]
            GLMakie.scatter!(
                loc_x[i],
                loc_y[i];
                marker     = :xcross,
                markersize = marker_size * 1.8,
                color      = :red,
            )
        end
    end

    # in the side (:yz, seen from the right) and coronal (:xz, seen from behind) views, channels of both
    # hemispheres (or the front and the back) project onto the same point: of overlapping channels, only
    # the one nearer to the viewer is labelled
    depth = if plane === :yz
        locs.loc_x[ch]
    elseif plane === :xz
        -locs.loc_y[ch]
    else
        zeros(length(ch))
    end
    hidden = falses(length(ch))
    if plane !== :xy
        for i in eachindex(ch), j in eachindex(ch)
            i == j && continue
            hypot(loc_x[i] - loc_x[j], loc_y[i] - loc_y[j]) < 0.06 &&
                depth[j] > depth[i] + 1e-3 && (hidden[i] = true)
        end
    end

    # draw labels below the markers (offset in pixels, independent of the plot scale)
    if ch_labels
        ch_set = Set(ch)
        for idx in eachindex(locs[!, :label])
            if idx in ch_set
                local_i = findfirst(==(idx), collect(ch))
                (isnothing(local_i) || hidden[local_i]) && continue
                GLMakie.text!(
                    loc_x[local_i],
                    loc_y[local_i];
                    text     = locs[!, :label][idx],
                    align    = (:center, :top),
                    offset   = (0, -(marker_size / 2 + 2)),
                    fontsize = font_size,
                )
            end
        end
    end

    # draw head labels
    head_labels && _draw_head_labels!(ax; plane = plane, font_size = font_size)

    # draw connection weight labels
    if !isnothing(connections)
        for idx1 in axes(connections, 1)
            for idx2 = (idx1 + 1):size(connections, 1)
                val = connections[idx1, idx2]
                if _passes_threshold(val, threshold, threshold_type)
                    _draw_connection_label!(
                        loc_x,
                        loc_y,
                        idx1,
                        idx2,
                        val,
                        font_size,
                        mono,
                    )
                end
            end
        end
    end

    # draw per-channel weight values
    if weights isa Vector
        length(weights) <= length(ch) ||
            throw(
                ArgumentError(
                    "Number of weights ($(length(weights))) must be ≤ number of channels ($(length(ch))).",
                ),
            )
        length(weights) >= 1 ||
            throw(ArgumentError("weights must contain at least one value."))

        for (i, idx) in enumerate(collect(ch))
            i > length(weights) && break
            color = mono ? :black : (weights[i] >= 0 ? :red : :blue)
            # above the electrode marker (offset in pixels; channel labels are below it)
            GLMakie.text!(
                loc_x[i],
                loc_y[i];
                text     = _num_label(weights[i]),
                fontsize = font_size,
                color    = color,
                align    = (:center, :bottom),
                offset   = (0, marker_size / 2 + 2),
            )
        end
    end

    loc_x_range = [(loc_x[i] - 0.02, loc_x[i] + 0.02) for i in eachindex(loc_x)]
    loc_y_range = [(loc_y[i] - 0.02, loc_y[i] + 0.02) for i in eachindex(loc_y)]

    # mouse events
    if gui
        println()
        on(events(fig).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press
                ax_x = mouseposition(ax)[1]
                ax_y = mouseposition(ax)[2]
                for idx in eachindex(loc_x)
                    if ax_x >= loc_x_range[idx][1] && ax_x <= loc_x_range[idx][2] &&
                       ax_y >= loc_y_range[idx][1] && ax_y <= loc_y_range[idx][2]
                        println(ch_info[idx])
                        break
                    end
                end
            end
        end
        wait(display(fig))
    end

    return fig
end

"""
    plot_locs(obj; <keyword arguments>)

Preview channel locations from a NEURO object with customizable 2D/3D visualization and optional connection mapping.

All selected channels are shown; channels marked as bad are drawn with a red × (not for NIRS objects).

For NIRS objects, the optodes (sources in red, detectors in green) and the source-detector pairs of the selected channels (`ch`; channels without optodes, e.g. auxiliary channels, are skipped) are shown; short-separation pairs (source-detector distance < 15 mm) are drawn with dashed lines. `sch` highlights the pairs of the given channels, `ch_labels` shows the channels (pairs) at the midpoints between their source and detector, numbered in the order of the pairs (as in `nirs_quality()`), and a `weights` vector gives one value per pair. `connections` are not supported for NIRS objects.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `sch::Union{String, Vector{String}, Regex}`: significant channels to highlight
- `ch_labels::Bool=true`: plot channel labels
- `src_labels::Bool=false`: if `true`, plot source labels (for NIRS data)
- `det_labels::Bool=false`: if `true`, plot detector labels (for NIRS data)
- `opt_labels::Bool=false`: if `true`, plot optode type (S for source, D for detector) and number
- `head::Bool=true`: if `true`, draw head outline
- `head_labels::Bool=false`: if `true`, draw head labels
- `mono::Bool=false`: if `true`, use a monochrome palette
- `grid::Bool=true`: if `true`, draw grid for locating positions
- `ps::Symbol`: plot size:
    - `:l`: large (800×800 px)
    - `:m`: medium (300×300 px)
    - `:s`: small (100×100 px)
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `plane::Symbol=:xy`: which plane to plot:
    - `:xy`: horizontal (top)
    - `:xz`: coronary (front)
    - `:yz`: sagittal (side)
- `connections::Union{Nothing, Matrix{<:Real}}=nothing`: matrix of connections weights, shape (channels, channels)
- `threshold::Real=0`: threshold for plotting connections
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: values equal to threshold
    - `:neq`: values not equal to threshold
    - `:geq`: values ≥ threshold
    - `:leq`: values ≤ threshold
    - `:g`: values > threshold
    - `:l`: values < threshold
    - `:in`: values in the threshold range (inclusive)
    - `:bin`: values in the threshold range (exclusive)
- `weights::Union{Bool, Vector{<:Real}}=true`: if `true`, auto-scale line widths and transparency based on connection strength; if Vector, use provided weights to place at channel location
- `schematic::Bool=false`: NIRS only: if `true`, draw a schematic montage (sources as diamonds, detectors as circles, channels as numbered squares, legend), see `plot_locs_nirs()`
- `gui::Bool=true`: if `true`, keep window open and interactive

# Returns

- `Union{GLMakie.Figure, Nothing}`
"""
function plot_locs(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    sch::Union{String, Vector{String}, Regex} = "",
    ch_labels::Bool = true,
    src_labels::Bool = false,
    det_labels::Bool = false,
    opt_labels::Bool = false,
    head::Bool = true,
    head_labels::Bool = false,
    mono::Bool = false,
    grid::Bool = false,
    ps::Symbol = :l,
    cart::Bool = false,
    plane::Symbol = :xy,
    connections::Union{Nothing, Matrix{<:Real}} = nothing,
    threshold::Real = 0,
    threshold_type::Symbol = :neq,
    weights::Union{Bool, Vector{<:Real}} = true,
    schematic::Bool = false,
    gui::Bool = true,
)::Union{GLMakie.Figure, Nothing}
    # validate

    # all selected channels are shown (bad channels are marked)
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # channel types of the selected channels (objects may contain several, e.g. EEG and NIRS)
    ct = obj.header.recording[:channel_type][ch]
    all(==("ecog"), ct) && throw(ArgumentError("Use plot_locs_ecog() for ECoG channels."))
    _nirs_n_opt(obj.header.recording) > 0 && all(t -> startswith(t, "nirs"), ct) &&
        return _plot_locs_nirs_obj(
        obj, ch;
        sch          = sch,
        ch_labels    = ch_labels,
        src_labels   = src_labels,
        det_labels   = det_labels,
        opt_labels   = opt_labels,
        head         = head,
        head_labels  = head_labels,
        mono         = mono,
        grid         = grid,
        ps           = ps,
        cart         = cart,
        plane        = plane,
        connections  = connections,
        weights      = weights,
        schematic    = schematic,
        gui          = gui,
    )

    ch_info = String[]
    for idx in eachindex(ch)
        push!(ch_info, channel_info(obj; ch = labels(obj)[ch[idx]], pr = false))
    end

    chs  = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), _locs_unit(obj.locs))
    ch   = collect(1:DataFrames.nrow(locs))

    sch_resolved = if sch == ""
        Int64[]
    else
        sch_idx = get_channel(obj; ch = sch)
        sch_chs = intersect(locs[!, :label], labels(obj)[sch_idx])
        _find_bylabel(locs, sch_chs)
    end

    if any(in(("eeg", "csd", "mag", "grad")), ct) || any(t -> startswith(t, "nirs"), ct)
        return plot_locs(
            locs;
            ch             = ch,
            sch            = sch_resolved,
            ch_labels      = ch_labels,
            head           = head,
            head_labels    = head_labels,
            grid           = grid,
            ps             = ps,
            mono           = mono,
            cart           = cart,
            plane          = plane,
            connections    = connections,
            threshold      = threshold,
            threshold_type = threshold_type,
            weights        = weights,
            ch_info        = ch_info,
            bad            = in(Set(labels(obj)[obj.header.recording[:bad_channel]])).(locs[!, :label]),
            gui            = gui,
        )
    elseif any(in(("seeg", "ieeg")), ct)
        _warn("sEEG and iEEG locs are not supported yet.")
        return nothing
    end

    return throw(
        ArgumentError(
            "Locations plot is not available for channels of type $(join(unique(ct), ", ")).",
        ),
    )
end

# plot_locs() for NIRS objects: optodes and the source-detector pairs of the channels `ch`
function _plot_locs_nirs_obj(
    obj::NeuroAnalyzer.NEURO,
    ch::Vector{Int64};
    sch::Union{String, Vector{String}, Regex},
    ch_labels::Bool,
    src_labels::Bool,
    det_labels::Bool,
    opt_labels::Bool,
    head::Bool,
    head_labels::Bool,
    mono::Bool,
    grid::Bool,
    ps::Symbol,
    cart::Bool,
    plane::Symbol,
    connections::Union{Nothing, Matrix{<:Real}},
    weights::Union{Bool, Vector{<:Real}},
    schematic::Bool,
    gui::Bool,
)::GLMakie.Figure
    rec = obj.header.recording
    opt_pairs = rec[:optode_pairs]
    n_opt = size(opt_pairs, 1)
    # channels with optodes (intensity, OD, concentration)
    ch_opt = [c for c in ch if c <= n_opt && all(>(0), opt_pairs[c, :])]
    isempty(ch_opt) && throw(ArgumentError("No NIRS channels with optodes selected."))
    isnothing(connections) || _warn("connections are not supported for NIRS objects.")

    # short-separation pairs (if the distances can be calculated)
    short = try
        nirs_distance(obj).d[ch_opt] .< 15
    catch
        nothing
    end

    # highlighted pairs: rows (in ch_opt) of the channels in sch
    spairs = Int64[]
    if sch != ""
        sch_idx = get_channel(obj; ch = sch)
        sp = Set((opt_pairs[c, 1], opt_pairs[c, 2]) for c in sch_idx if c <= n_opt)
        spairs = [k for (k, c) in enumerate(ch_opt) if (opt_pairs[c, 1], opt_pairs[c, 2]) in sp]
    end

    fig = plot_locs_nirs(
        obj.locs, opt_pairs[ch_opt, :], length(rec[:src_labels]), length(rec[:det_labels]);
        src_names   = string.(rec[:src_labels]),
        det_names   = string.(rec[:det_labels]),
        short       = short,
        spairs      = spairs,
        ch_mid      = ch_labels,
        weights     = weights isa Vector ? weights : nothing,
        schematic   = schematic,
        src_labels  = src_labels,
        det_labels  = det_labels,
        opt_labels  = opt_labels,
        ps          = ps,
        head        = head,
        head_labels = head_labels,
        cart        = cart,
        grid        = grid,
        mono        = mono,
        plane       = plane,
    )
    gui && wait(display(fig))

    return fig
end
