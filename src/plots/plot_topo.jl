export plot_topo

"""
    plot_topo(s; <keyword arguments>)

Plot a topographical map of signal values across channel locations.

# Arguments

- `s::AbstractVector`: signal values to plot (one value per channel)
- `locs::DataFrame`: channel location data
- `ch::Union{Int64, Vector{Int64}}=1:DataFrames.nrow(locs)`: channels to include; defaults to all rows in `locs`
- `sch::Union{Nothing, Int64, Vector{Int64}}=nothing`: significant channels to highlight
- `cb::Bool=true`: if `true`, show colorbar
- `cb_title::String="[A.U.]"`: colorbar title
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: Inverse Multiquadratic
    - `:tp`: ThinPlate
    - `:nn`: Nearest Neighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `contours::Int64=0`: plot contours (if > 0) over topo plot, number specifies how many levels to plot
- `electrodes::Bools=true`: if `true`, plot electrode locations over topography
- `ps::Symbol`: plot size:
    - `:l`: large (800×800 px)
    - `:m`: medium (300×300 px)
    - `:s`: small (100×100 px)
- `head::Bool=true`: if `true`, draw head outline
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
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
- `threshold_method::Symbol=:reg`: thresholding method:
    - `:reg`: threshold whole topomap region (default)
    - `:loc`: threshold only at channel locations

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_topo(
    s::AbstractVector;
    locs::DataFrame,
    ch::Union{Int64, Vector{Int64}} = 1:DataFrames.nrow(locs),
    sch::Union{Nothing, Int64, Vector{Int64}} = nothing,
    cb::Bool = true,
    cb_title::String = "[A.U.]",
    title::String = "",
    mono::Bool = false,
    imethod::Symbol = :sh,
    nmethod::Symbol = :minmax,
    contours::Int64 = 0,
    electrodes::Bool = true,
    ps::Symbol = :l,
    head::Bool = true,
    cart::Bool = false,
    threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    threshold_type::Symbol = :neq,
    threshold_method::Symbol = :reg,
)::GLMakie.Figure
    # validate
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(threshold_type, [:eq, :neq, :geq, :leq, :g, :l, :in, :bin], "threshold_type")
    _check_var(ps, [:l, :m, :s], "ps")
    _check_var(threshold_method, [:reg, :loc], "threshold_method")
    contours >= 0 || throw(ArgumentError("contours must be ≥ 0."))
    if !isnothing(sch)
        length(intersect(ch, sch)) == length(sch) ||
            throw(ArgumentError("Some sch channels were not found in ch."))
        !isnothing(threshold) &&
            throw(ArgumentError("Both sch and threshold cannot be specified."))
    end

    # set color palette
    pal = mono ? :grays : :bluesreds

    # number of channels
    ch_n = length(ch)

    # channel locations by labels
    local_locs = locs[ch, :]

    if cart
        # cartesian coordinates
        loc_x = local_locs.loc_x
        loc_y = local_locs.loc_y
    else
        # polar coordinates
        loc_x = zeros(ch_n)
        loc_y = zeros(ch_n)
        for idx in eachindex(ch)
            loc_x[idx], loc_y[idx] =
                pol2cart(local_locs.loc_radius[idx], local_locs.loc_theta[idx])
        end
    end

    # plot parameters
    if ps === :l
        plot_size   = (800, 800)
        marker_size = length(ch) > 64 ? 8 : 16
        iter        = 512
        font_size   = 20
        lw          = 3
        sw          = 4
        !occursin("\n", title) && title !== "" && (title *= "\n")
    elseif ps === :m
        plot_size   = (300, 300)
        marker_size = length(ch) > 64 ? 2 : 4
        iter        = 256
        font_size   = 10
        cb_title    = ""
        lw          = 2
        sw          = 2
    elseif ps === :s
        plot_size   = (100, 100)
        marker_size = length(ch) > 64 ? 1 : 2
        iter        = 128
        font_size   = 5
        title       = ""
        cb          = false
        contours    = 0
        cb_title    = ""
        lw          = 1
        sw          = 1
    end

    # interpolate signal
    s_interpolated, interpolated_x, interpolated_y =
        _interpolate2d(s, loc_x, loc_y, iter, imethod, nmethod)
    s_interpolated = s_interpolated'[:, end:-1:1]
    s_interpolated_threshold = deepcopy(s_interpolated)

    # compute thresholded region before removing peripheral values
    threshold_idx = nothing
    if !isnothing(threshold)
        if threshold_method === :loc
            if threshold_type in [:eq, :neq, :geq, :leq, :g, :l]
                length(threshold) == 1 ||
                    throw(ArgumentError("threshold must contain a single value."))
            else
                length(threshold) == 2 ||
                    throw(ArgumentError("threshold must contain two values."))
                s_norm = normalize(s; method = nmethod)
                _check_tuple(threshold, extrema(s_norm), "threshold")
            end

            s_norm = normalize(s; method = nmethod)

            if threshold_type === :eq
                threshold_idx = findall(x -> x == threshold, s_norm)
            elseif threshold_type === :neq
                threshold_idx = findall(x -> x != threshold, s_norm)
            elseif threshold_type === :geq
                threshold_idx = findall(x -> x >= threshold, s_norm)
            elseif threshold_type === :leq
                threshold_idx = findall(x -> x <= threshold, s_norm)
            elseif threshold_type === :g
                threshold_idx = findall(x -> x > threshold, s_norm)
            elseif threshold_type === :l
                threshold_idx = findall(x -> x < threshold, s_norm)
            elseif threshold_type === :in
                threshold_idx =
                    findall(x -> x >= threshold[1] && x <= threshold[2], s_norm)
            elseif threshold_type === :bin
                threshold_idx =
                    findall(x -> x > threshold[1] && x < threshold[2], s_norm)
            end
        else
            _, bm = seg_extract(
                s_interpolated;
                threshold      = threshold,
                threshold_type = threshold_type,
            )
            s_interpolated_threshold[.!bm] .= NaN
        end
    end

    head12 =
        maximum(abs.(local_locs.loc_x)) <= 1.2 &&
        maximum(abs.(local_locs.loc_y)) <= 1.2 &&
        maximum(abs.(local_locs.loc_z)) <= 1.5

    if head12
        xl = (-1.2, 1.2)
        yl = (-1.2, 1.2)
        r  = 1.2
    else
        xl = (-1.6, 1.6)
        yl = (-1.6, 1.6)
        r  = 1.6
    end

    if head12
        # get distances from (0, 0)
        d = zeros(length(interpolated_x), length(interpolated_y))
        for idx1 in eachindex(interpolated_x)
            for idx2 in eachindex(interpolated_y)
                d[idx1, idx2] =
                    distance((0, 0), (interpolated_x[idx1], interpolated_y[idx2]))
            end
        end
        s_interpolated[d .>= xl[2]] .= NaN
        !isnothing(threshold) && (s_interpolated_threshold[d .>= xl[2]] .= NaN)
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_topo()")
    fig = GLMakie.Figure(;
        size = plot_size,
        figure_padding = ps in [:l, :m] ? (10, 10, 10, 0) : (0, 0, 0, 0), # L R B T
    )

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        aspect = 1,
        xlabel = "",
        ylabel = "",
        title = title,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        backgroundcolor = :transparent,
        titlesize = font_size,
        _AXIS_LOCK_KWARGS...,
    )
    hidedecorations!(ax)
    hidespines!(ax)
    GLMakie.xlims!(ax, xl)
    GLMakie.ylims!(ax, yl)

    if !isnothing(threshold) && threshold_method === :reg
        hm = GLMakie.heatmap!(
            ax,
            interpolated_x,
            interpolated_y,
            s_interpolated_threshold;
            colorrange = extrema(s_interpolated[.!isnan.(s_interpolated)]),
            colormap   = pal,
        )
    else
        hm = GLMakie.heatmap!(
            ax,
            interpolated_x,
            interpolated_y,
            s_interpolated;
            colormap = pal,
        )
    end

    # draw contours over the unthresholded heatmap
    # skipped when threshold_method === :reg because the thresholded heatmap
    # already marks the region visually.
    if contours > 0 &&
       (isnothing(threshold) || threshold_method === :loc)
        GLMakie.contour!(
            ax,
            interpolated_x,
            interpolated_y,
            s_interpolated;
            linestyle = :dash,
            levels    = contours,
            linewidth = 0.5,
            color     = :black,
        )
    end

    # draw head outline
    head && _draw_head_outline!(ax; lw = lw)

    # draw electrodes, highlighting thresholded or significant channels
    if electrodes
        if (isnothing(threshold) && isnothing(sch)) ||
           (!isnothing(threshold) && threshold_method === :reg)
            for idx in 1:ch_n
                GLMakie.scatter!(
                    ax, loc_x[idx], loc_y[idx];
                    markersize = marker_size,
                    color      = :black,
                )
            end

        elseif threshold_method === :loc
            for idx in 1:ch_n
                if idx in threshold_idx
                    GLMakie.scatter!(
                        ax, loc_x[idx], loc_y[idx];
                        markersize  = marker_size * 2,
                        color       = :gray,
                        strokewidth = sw,
                        strokecolor = :black,
                    )
                else
                    GLMakie.scatter!(
                        ax, loc_x[idx], loc_y[idx];
                        markersize = marker_size,
                        color      = :black,
                    )
                end
            end

        elseif !isnothing(sch)
            for idx in 1:ch_n
                if idx in sch
                    GLMakie.scatter!(
                        ax, loc_x[idx], loc_y[idx];
                        markersize  = marker_size * 2,
                        color       = :gray,
                        strokewidth = sw,
                        strokecolor = :black,
                    )
                else
                    GLMakie.scatter!(
                        ax, loc_x[idx], loc_y[idx];
                        markersize = marker_size,
                        color      = :black,
                    )
                end
            end
        end
    end

    # draw mask to crop interpolation outside the head circle
    head12 && GLMakie.arc!(ax, Point2f(0, 0), r, -pi, pi; linewidth = 5, color = :white)

    # draw colorbar
    if cb
        GLMakie.Colorbar(
            fig[1, 2],
            hm;
            label = cb_title,
            labelsize = font_size - 4,
            ticklabelsize = font_size - 4,
            height = div(plot_size[2], 2),
            width = ps === :l ? 25 : 10,
            tellheight = false,
        )
        rowsize!(fig.layout, 1, ax.scene.viewport[].widths[2])
        colgap!(fig.layout, 10)
    end

    resize_to_layout!(fig)

    return fig
end

"""
    plot_topo(obj; <keyword arguments>)

Plot a topographical map of signal values from a NEURO object with customizable visualization options.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Union{Nothing, AbstractVector, AbstractMatrix}=nothing`: external data to plot:
    - `Vector`: one value per channel
    - `Matrix`: (channels × values), will be averaged by channels
    - `nothing`: use data from NEURO object at specified time point(s) (default)
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `sch::Union{Nothing, String, Vector{String}, Regex}=nothing`: significant channels to highlight
- `tpos::Union{Nothing, Real, AbstractVector}=nothing`: time point in seconds to plot, ignored if `data` is provided
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `cb::Bool=true`: if `true`, show colorbar
- `cb_title::String="[A.U.]"`: colorbar title
- `amethod::Symbol=:mean`: averaging method for matrix data:
    - `:mean`: mean averaging 
    - `:median`: median averaging
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: Inverse Multiquadratic
    - `:tp`: ThinPlate
    - `:nn`: Nearest Neighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `contours::Int64=0`: plot contours (if > 0) over topo plot, number specifies how many levels to plot
- `electrodes::Bools=true`: if `true`, plot electrode locations over topography
- `ps::Symbol`: plot size:
    - `:l`: large (800×800 px)
    - `:m`: medium (300×300 px)
    - `:s`: small (100×100 px)
- `head::Bool=true`: if `true`, draw head outline
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
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
- `threshold_method::Symbol=:reg`: thresholding method:
    - `:reg`: threshold whole topomap region (default)
    - `:loc`: threshold only at channel locations
- `nr::Int64=0`: number of rows for arranging multiple topomaps
- `nc::Int64=0`: number of columns for arranging multiple topomaps

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_topo(
    obj::NeuroAnalyzer.NEURO;
    data::Union{Nothing, AbstractArray} = nothing,
    ch::Union{String, Vector{String}, Regex},
    sch::Union{Nothing, String, Vector{String}, Regex} = nothing,
    tpos::Union{Nothing, Real, AbstractVector} = nothing,
    title::String = "default",
    mono::Bool = false,
    cb::Bool = true,
    cb_title::String = "default",
    amethod::Symbol = :mean,
    imethod::Symbol = :sh,
    nmethod::Symbol = :minmax,
    contours::Int64 = 0,
    electrodes::Bool = true,
    ps::Symbol = :l,
    head::Bool = true,
    cart::Bool = false,
    threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    threshold_type::Symbol = :neq,
    threshold_method::Symbol = :reg,
    nr::Int64 = 1,
    nc::Int64 = 0,
)::GLMakie.Figure

    # TO DO: vector of tpos:
    # generate separate plots, put them in nr × nc matrix and add one shared colorbar
    if !isnothing(tpos) && tpos isa AbstractVector && length(tpos) > 1
        if nr == 1
            nc = length(tpos)
            # FIX: was `nr > 1 & nc == 0` — bitwise & has higher precedence than >, so this
            #      parsed as `nr > (1 & nc) == 0` rather than `(nr > 1) && (nc == 0)`
        elseif nr > 1 && nc == 0
            nc = ceil(Int64, length(tpos) / nr)
        elseif nc != 0
            nr = ceil(Int64, length(tpos) / nc)
        end
        _warn("Vector of tpos is not yet supported; using first time point only.")
        # FIX: was `collect(tpos)[1]` — allocates unnecessarily; just take first element
        tpos = first(tpos)
    end

    # validate
    contours >= 0 || throw(ArgumentError("contours must be ≥ 0."))
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(amethod, [:mean, :median], "amethod")
    _check_var(
        nmethod,
        [
            :zscore, :minmax, :log, :log10, :neglog, :neglog10,
            :neg, :pos, :perc, :gauss, :invroot, :n,
            :softmax, :sigmoid, :mad, :rank, :none,
        ],
        "nmethod",
    )

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # significant channels
    if !isnothing(sch)
        if isa(sch, String)
            length(intersect(ch, get_channel(obj; ch = sch))) == 1 ||
                throw(ArgumentError("sch channel was not found in ch."))
        else
            length(intersect(ch, get_channel(obj; ch = sch))) == length(sch) ||
                throw(ArgumentError("Some sch channels were not found in ch."))
        end
    end

    length(ch) >= 2 || throw(ArgumentError("plot_topo() requires ≥ 2 channels."))

    chs = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
    !isnothing(sch) && (sch = _find_bylabel(locs, sch))

    # prepare data or time position
    if isnothing(data)
        isnothing(tpos) && throw(ArgumentError("Either tpos or data must be provided."))
        tpos >= obj.time_pts[1] ||
            throw(ArgumentError("tpos must be ≥ $(obj.time_pts[1])"))
        tpos <= obj.time_pts[end] ||
            throw(ArgumentError("tpos must be ≤ $(obj.time_pts[end])"))

        tpos_idx = vsearch(tpos, obj.time_pts)
        title == "default" && (title = "$(obj.time_pts[tpos_idx]) s")

        data = if nepochs(obj) == 1
            obj.data[ch, tpos_idx, 1]
        else
            reshape(obj.data, size(obj.data, 1), size(obj.data, 2) * size(obj.data, 3))[
                ch, tpos_idx,
            ]
        end

    else
        !isnothing(tpos) && _info("tpos is ignored when data is provided.")
        if ndims(data) == 2
            data = amethod === :mean ? mean(data; dims = 2)[:] : median(data; dims = 2)[:]
        end
        length(data) == length(ch) || throw(
            ArgumentError(
                "Number of data values ($(length(data))) must equal the number of channels ($(length(ch))).",
            ),
        )
        title == "default" && (title = "")
    end

    # colorbar title
    cb_title == "default" && (cb_title = "[A.U.]")

    return plot_topo(
        data;
        locs             = locs,
        ch               = collect(1:DataFrames.nrow(locs)),
        sch              = sch,
        cb               = cb,
        cb_title         = cb_title,
        title            = title,
        mono             = mono,
        imethod          = imethod,
        nmethod          = nmethod,
        contours         = contours,
        electrodes       = electrodes,
        ps               = ps,
        head             = head,
        cart             = cart,
        threshold        = threshold,
        threshold_type   = threshold_type,
        threshold_method = threshold_method,
    )
end
