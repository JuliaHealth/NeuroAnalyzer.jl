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
    loc_x::AbstractVector, loc_y::AbstractVector,
    idx1::Int, idx2::Int,
    val::Real, weight::Real,
    mono::Bool, use_weights::Bool,
)
    xs = [loc_x[idx1], loc_x[idx2]]
    ys = [loc_y[idx1], loc_y[idx2]]
    if use_weights
        lw = 6 * weight
        al = 0.25 * weight
        if val > 0
            GLMakie.lines!(xs, ys; linewidth = lw, alpha = al,
                color = mono ? :black : :red)
        elseif val < 0
            GLMakie.lines!(xs, ys; linewidth = lw, alpha = al,
                color = mono ? :black : :blue,
                linestyle = mono ? :dot : :solid)
        end
    else
        GLMakie.lines!(xs, ys; linewidth = 0.2, color = :black)
    end
end

# Draw a connection weight label at the midpoint between two channel positions.
function _draw_connection_label!(
    loc_x::AbstractVector, loc_y::AbstractVector,
    idx1::Int, idx2::Int,
    val::Real, font_size::Int, mono::Bool,
)
    l_pos = _midxy(loc_x[idx1], loc_y[idx1], loc_x[idx2], loc_y[idx2])
    color = mono ? :black : (val >= 0 ? :red : :blue)
    return GLMakie.text!(
        l_pos[1], l_pos[2];
        align    = (:center, :center),
        text     = string(val),
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
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=1:DataFrames.nrow(locs)`: list of locations to plot, default is all locations
- `sch::Union{Int64, Vector{Int64}, AbstractRange}=0`: significant channels to highlight
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

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_locs(
    locs::DataFrame;
    ch::Union{Int64, Vector{Int64}, AbstractRange} = 1:DataFrames.nrow(locs),
    sch::Union{Int64, Vector{Int64}, AbstractRange} = 0,
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
    gui::Bool = true,
)::GLMakie.Figure
    # validate
    _check_var(ps, [:l, :m, :s], "ps")
    _check_var(plane, [:xy, :yz, :xz], "plane")

    # set color palette
    pal = mono ? :grays : :darktest

    # significant channel labels
    sch_labels = ch_labels

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
        maximum(abs.(locs.loc_x)) <= 1.2 &&
        maximum(abs.(locs.loc_y)) <= 1.2 &&
        maximum(abs.(locs.loc_z)) <= 1.5

    xl = head12 ? (-1.2, 1.2) : (-1.6, 1.6)
    yl = head12 ? (-1.2, 1.2) : (-1.6, 1.6)

    # plot parameters
    if ps === :l
        plot_size   = (800, 800)
        marker_size = length(ch) > 64 ? 10 : 20
        font_size   = 14
        lw          = 3
        sw          = 2
    elseif ps === :m
        plot_size   = (300, 300)
        marker_size = length(ch) > 64 ? 5 : 10
        font_size   = 8
        lw          = 2
        sw          = 1
        ch_labels   = false
        sch_labels  = false
        grid        = false
    elseif ps === :s
        plot_size   = (100, 100)
        marker_size = length(ch) > 64 ? 4 : 8
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
        if plane === :xy
            _draw_head_outline!(ax; lw = lw)
        elseif plane === :yz
            # head
            GLMakie.arc!(ax, (0, 0), 1, 0, pi; linewidth = lw, color = :black)
        elseif plane === :xz
            # head
            GLMakie.arc!(ax, (0, 0), 1, 0, pi; linewidth = lw, color = :black)
        end
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
            for idx2 in (idx1 + 1):size(connections, 1)
                val = connections[idx1, idx2]
                if _passes_threshold(val, threshold, threshold_type)
                    _draw_connection!(loc_x, loc_y, idx1, idx2,
                        val, m_tmp[idx1, idx2], mono, use_weights)
                end
            end
        end
    end

    # draw channel markers
    ch_n = length(ch)
    cmap = GLMakie.resample_cmap(pal, ch_n)
    sch_set = Set(sch)

    for (i, idx) in enumerate(ch)
        if idx in sch_set
            GLMakie.scatter!(
                loc_x[i], loc_y[i];
                markersize  = marker_size,
                color       = mono ? :gray : cmap[i],
                colormap    = pal,
                colorrange  = 1:ch_n,
                strokewidth = sw,
                strokecolor = :black,
            )
        else
            GLMakie.scatter!(
                loc_x[i], loc_y[i];
                markersize  = marker_size,
                color       = :gray,
                strokewidth = sw,
                strokecolor = :black,
            )
        end
    end

    label_offset_x = 0.0
    label_offset_y = -0.08

    # draw labels
    if ch_labels
        ch_set = Set(ch)
        for idx in eachindex(locs[!, :label])
            if idx in ch_set
                local_i = findfirst(==(idx), collect(ch))
                isnothing(local_i) && continue
                GLMakie.text!(
                    loc_x[local_i] + label_offset_x,
                    loc_y[local_i] + label_offset_y;
                    text     = locs[!, :label][idx],
                    align    = (:center, :bottom),
                    fontsize = font_size,
                )
            end
        end
    end

    # draw head labels
    if head_labels
        fid_names = ["NAS", "IN", "LPA", "RPA"]
        for idx in eachindex(NeuroAnalyzer.fiducial_points)
            pt = NeuroAnalyzer.fiducial_points[idx]
            fid_loc_x, fid_loc_y = if plane === :xy
                pt[1], pt[2]
            elseif plane === :xz
                pt[1], pt[3]
            elseif plane === :yz
                pt[2], pt[3]
            end
            GLMakie.text!(
                fid_loc_x, fid_loc_y;
                text     = fid_names[idx],
                fontsize = font_size,
                align    = (:center, :center),
            )
        end
    end

    # draw connection weight labels
    if !isnothing(connections)
        for idx1 in axes(connections, 1)
            for idx2 in (idx1 + 1):size(connections, 1)
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
        label_offset_x = 0.0
        label_offset_y = 0.07
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
            GLMakie.text!(
                loc_x[i] + label_offset_x,
                loc_y[i] + label_offset_y;
                text     = string(weights[i]),
                fontsize = font_size,
                color    = color,
                align    = (:center, :top),
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
    gui::Bool = true,
)::Union{GLMakie.Figure, Nothing}

    # validate
    datatype(obj) != "ecog" || throw(ArgumentError("Use plot_locs_ecog() for ECoG data."))

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")

    ch_info = String[]
    for idx in eachindex(ch)
        push!(ch_info, channel_info(obj; ch = labels(obj)[ch[idx]], pr = false))
    end

    chs  = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    ch   = collect(1:DataFrames.nrow(locs))

    sch_resolved = if sch == ""
        Int64[]
    else
        sch_idx =
            exclude_bads ?
            get_channel(obj; ch = sch, exclude = "bad") :
            get_channel(obj; ch = sch, exclude = "")
        sch_chs = intersect(locs[!, :label], labels(obj)[sch_idx])
        _find_bylabel(locs, sch_chs)
    end

    if datatype(obj) in ["eeg", "meg", "csd", "erp", "erf"]
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
            gui            = gui,
        )
    elseif datatype(obj) == "nirs"
        opt_pairs = obj.header.recording[:optode_pairs]
        src_n = length(source_labels(obj))
        det_n = length(detector_labels(obj))
        return plot_locs_nirs(
            obj.locs, opt_pairs, src_n, det_n;
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
            ch_info     = ch_info,
        )
    elseif datatype(obj) == "ecog"
        _warn("ECOG locs are not supported yet.")
        return nothing
    elseif datatype(obj) == "seeg"
        _warn("SEEG locs are not supported yet.")
        return nothing
    elseif datatype(obj) == "ieeg"
        _warn("iEEG locs are not supported yet.")
        return nothing
    elseif datatype(obj) in ["sensors", "eda", "mep", "tpt"]
        _warn("For $(datatype(obj)) object type locs are not available.")
        return nothing
    end

    # should never be reached, but satisfies the return type
    return fig
end
