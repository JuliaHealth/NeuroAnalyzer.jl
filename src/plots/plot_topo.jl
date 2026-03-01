export plot_topo

"""
    plot_topo(s; <keyword arguments>)

Plot topographical view.

# Arguments

  - `s::AbstractVector`: values to plot (one value per channel)
  - `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
  - `ch::Union{Int64, Vector{Int64}}=1:DataFrames.nrow(locs)`: list of channels, default is all channels
  - `sch::Union{Nothing, Int64, Vector{Int64}}=nothing`: list of significant channels
  - `cb::Bool=true`: plot colorbar
  - `cb_title::String="[A.U.]"`: colorbar title
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette
  - `imethod::Symbol=:sh`: interpolation method:
      + `:sh`: Shepard
      + `:mq`: Multiquadratic
      + `:imq`: InverseMultiquadratic
      + `:tp`: ThinPlate
      + `:nn`: NearestNeighbour
      + `:ga`: Gaussian
  - `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
  - `contours::Int64=0`: plot contours (if > 0) over topo plot, number specifies how many levels to plot
  - `electrodes::Bools=true`: plot electrodes over topo plot
  - `ps::Symbol=:l`: plot size (`:l`: large (800×800 px), `:m`: medium (300×300 px), `:s`: small (100×100 px))
  - `head::Bool=true`: draw head
  - `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
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
  - `threshold_method::Symbol=:reg`: thresholding method: threshold the whole topomap region (`:reg`) or only signal at channels locations (`:loc`)

# Returns

  - `p::GLMakie.Figure`
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

    pal = mono ? :grays : :bluesreds
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(threshold_type, [:eq, :neq, :geq, :leq, :g, :l, :in, :bin], "threshold_type")
    _check_var(ps, [:l, :m, :s], "ps")
    _check_var(threshold_method, [:reg, :loc], "threshold_method")
    @assert contours >= 0 "contours must be ≥ 0."
    if !isnothing(sch)
        @assert length(intersect(ch, sch)) == length(sch) "Some sch channels were not found in ch."
        @assert isnothing(threshold) "Both sch and threshold cannot be specified."
    end

    ch_n = length(ch)
    locs = locs[ch, :]

    if cart == false
        loc_x = zeros(length(ch))
        loc_y = zeros(length(ch))
        for idx in eachindex(ch)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
    else
        loc_x = locs[ch, :loc_x]
        loc_y = locs[ch, :loc_y]
    end

    if ps === :l
        plot_size = (800, 800)
        marker_size = length(ch) > 64 ? 8 : 16
        length(ch) > 64 && (ch_labels = false)
        iter = 512
        font_size = 20
        !occursin("\n", title) && title !== "" && (title *= "\n")
    elseif ps === :m
        plot_size = (300, 300)
        marker_size = length(ch) > 64 ? 2 : 4
        font_size = 8
        iter = 256
        cb_title = ""
        font_size = 10
    elseif ps === :s
        plot_size = (100, 100)
        marker_size = length(ch) > 64 ? 1 : 2
        iter = 128
        title = ""
        cb = false
        contours = 0
        cb_title = ""
        font_size = 5
    end

    s_interpolated, interpolated_x, interpolated_y = _interpolate2d(s, loc_x, loc_y, iter, imethod, nmethod)
    s_interpolated = s_interpolated'[:, end:-1:1]
    s_interpolated_threshold = deepcopy(s_interpolated)

    # we need to calculate thresholded region before removing peripherals
    threshold_idx = nothing
    if !isnothing(threshold)
        if threshold_method === :loc
            if threshold_type in [:eq, :neq, :geq, :leq, :g, :l]
                @assert length(threshold) == 1 "threshold must contain a single value."
            else
                @assert length(threshold) == 2 "threshold must contain two values."
                _check_tuple(threshold, "threshold")
            end
            s_norm = normalize(s, method = nmethod)
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
                threshold_idx = findall(x -> (x >= threshold[1] && x <= threshold[2]), s_norm)
            elseif threshold_type === :bin
                threshold_idx = findall(x -> (x > threshold[1] && x < threshold[2]), s_norm)
            end
        else
            _, bm = seg_extract(s_interpolated, threshold = threshold, threshold_type = threshold_type)
            s_interpolated_threshold[.!bm] .= NaN
        end
    end

    head12 = false
    maximum(abs.(locs[:, :loc_x])) <= 1.2 &&
        maximum(abs.(locs[:, :loc_y])) <= 1.2 &&
        maximum(abs.(locs[:, :loc_z])) <= 1.5 &&
        (head12 = true)

    if head12
        xl = (-1.2, 1.2)
        yl = (-1.2, 1.2)
        r = 1.2
    else
        xl = (-1.6, 1.6)
        yl = (-1.6, 1.6)
        r = 1.6
    end

    if head12
        # get distances from (0, 0)
        d = zeros(length(interpolated_x), length(interpolated_y))
        for idx1 in eachindex(interpolated_x)
            for idx2 in eachindex(interpolated_y)
                d[idx1, idx2] = distance((0, 0), (interpolated_x[idx1], interpolated_y[idx2]))
            end
        end
        # remove everything outside the radius
        s_interpolated[d .>= xl[2]] .= NaN
        !isnothing(threshold) && (s_interpolated_threshold[d .>= xl[2]] .= NaN)
    end

    # prepare plot
    GLMakie.activate!(title = "plot_topo()")
    p = GLMakie.Figure(
        size = plot_size,
        figure_padding = ps in [:l, :m] ? (10, 10, 10, 0) : (0, 0, 0, 0), # L R B T
    )
    ax = GLMakie.Axis(
        p[1, 1];
        aspect = 1,
        xlabel = "",
        ylabel = "",
        title = title,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        backgroundcolor = :transparent,
        titlesize = font_size,
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
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
            colormap = pal,
        )
    else
        hm = GLMakie.heatmap!(ax, interpolated_x, interpolated_y, s_interpolated, colormap = pal)
    end

    # draw contours
    if contours > 0 &&
            ((isnothing(threshold) && threshold_method === :reg) || (!isnothing(threshold) && threshold_method === :loc))
        GLMakie.contour!(
            ax,
            interpolated_x,
            interpolated_y,
            s_interpolated;
            linestyle = :dash,
            levels = contours,
            linewidth = 0.5,
            color = :black,
        )
    end

    # draw head
    if head
        ps === :l && (lw = 3)
        ps === :m && (lw = 2)
        ps === :s && (lw = 1)
        # nose
        GLMakie.lines!(ax, [-0.2, 0], [0.98, 1.08], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [0.2, 0], [0.98, 1.08], linewidth = lw, color = :black)

        # ears
        # left
        GLMakie.lines!(ax, [-0.995, -1.03], [0.1, 0.15], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [-1.03, -1.06], [0.15, 0.16], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [-1.06, -1.1], [0.16, 0.14], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [-1.1, -1.12], [0.14, 0.05], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [-1.12, -1.1], [0.05, -0.1], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [-1.1, -1.13], [-0.1, -0.3], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [-1.13, -1.09], [-0.3, -0.37], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [-1.09, -1.02], [-0.37, -0.39], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [-1.02, -0.98], [-0.39, -0.33], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [-0.98, -0.975], [-0.33, -0.22], linewidth = lw, color = :black)
        # right
        GLMakie.lines!(ax, [0.995, 1.03], [0.1, 0.15], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [1.03, 1.06], [0.15, 0.16], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [1.06, 1.1], [0.16, 0.14], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [1.1, 1.12], [0.14, 0.05], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [1.12, 1.1], [0.05, -0.1], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [1.1, 1.13], [-0.1, -0.3], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [1.13, 1.09], [-0.3, -0.37], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [1.09, 1.02], [-0.37, -0.39], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [1.02, 0.98], [-0.39, -0.33], linewidth = lw, color = :black)
        GLMakie.lines!(ax, [0.98, 0.975], [-0.33, -0.22], linewidth = lw, color = :black)

        # head
        GLMakie.arc!(ax, (0, 0), 1, 0, 2pi, linewidth = lw, color = :black)
    end

    # draw electrodes
    # mark thresholded or significant channels
    if electrodes
        ps === :l && (sw = 4)
        ps === :m && (sw = 2)
        ps === :s && (sw = 1)
        if (isnothing(threshold) && isnothing(sch)) || (!isnothing(threshold) && threshold_method === :reg)
            for idx in 1:ch_n
                GLMakie.scatter!(ax, loc_x[idx], loc_y[idx], markersize = marker_size, color = :black)
            end
        elseif threshold_method === :loc
            for idx in 1:ch_n
                if idx in threshold_idx
                    GLMakie.scatter!(
                        ax,
                        loc_x[idx],
                        loc_y[idx];
                        markersize = marker_size * 2,
                        color = :gray,
                        strokewidth = sw,
                        strokecolor = :black,
                    )
                else
                    GLMakie.scatter!(ax, loc_x[idx], loc_y[idx], markersize = marker_size, color = :black)
                end
            end
        elseif !isnothing(sch)
            for idx in 1:ch_n
                if idx in sch
                    GLMakie.scatter!(
                        ax,
                        loc_x[idx],
                        loc_y[idx];
                        markersize = marker_size * 2,
                        color = :gray,
                        strokewidth = sw,
                        strokecolor = :black,
                    )
                else
                    GLMakie.scatter!(ax, loc_x[idx], loc_y[idx], markersize = marker_size, color = :black)
                end
            end
        end
    end

    # draw mask
    if head12
        GLMakie.arc!(ax, Point2f(0), r, -pi, pi, linewidth = 5, color = :white)
    end

    # draw colorbar
    if cb
        GLMakie.Colorbar(
            p[1, 2],
            hm;
            label = cb_title,
            labelsize = font_size - 4,
            ticklabelsize = font_size - 4,
            height = div(plot_size[2], 2),
            width = ps === :l ? 25 : 10,
            tellheight = false,
        )
        rowsize!(p.layout, 1, ax.scene.viewport[].widths[2])
        colgap!(p.layout, 10)
    end

    resize_to_layout!(p)

    return p

end

"""
    plot_topo(obj; <keyword arguments>)

Topographical plot.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  - `data::Union{Nothing, AbstractVector, AbstractMatrix}=nothing`: external data to plot; vector: one value per channel; matrix: channels × values, will be averaged by channels
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `sch::Union{Nothing, String, Vector{String}, Regex}=nothing`: list of significant channels
  - `tpos::Union{Nothing, Real, AbstractVector}=nothing`: time point in seconds to plot, ignored if `data` is provided
  - `title::String="default"`: plot title, default is tpos value
  - `mono::Bool=false`: use color or gray palette
  - `cb::Bool=true`: plot colorbar
  - `cb_title::String="[A.U.]"`: colorbar title
  - `amethod::Symbol=:mean`: averaging method:
      + `:mean`
      + `:median`
  - `imethod::Symbol=:sh`: interpolation method:
      + `:sh`: Shepard
      + `:mq`: Multiquadratic
      + `:imq`: InverseMultiquadratic
      + `:tp`: ThinPlate
      + `:nn`: NearestNeighbour
      + `:ga`: Gaussian
  - `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
  - `contours::Int64=0`: plot contours (if > 0) over topo plot, number specifies how many levels to plot
  - `electrodes::Bools=true`: plot electrodes over topo plot
  - `ps::Symbol=:l`: plot size (`:l`: large (800×800 px), `:m`: medium (300×300 px), `:s`: small (100×100 px))
  - `head::Bool=true`: draw head
  - `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
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
  - `threshold_method::Symbol=:reg`: thresholding method: threshold the whole topomap region (`:reg`) or only signal at channels locations (`:loc`)
  - `nr::Int64=0`: number of rows to place topomaps
  - `nc::Int64=0`: number of columns to place topomaps

# Returns

  - `p::GLMakie.Figure`
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

    # TO DO: vector of tpos: generate separate plots, put them in nr × nc matrix and add one shared colorbar
    if length(tpos) > 1
        if nr == 1
            nc = length(tpos)
        elseif nr > 1 & nc == 0
            nc = ceil(Int64, length(tpos) / nr)
        elseif nc != 0
            nr = ceil(Int64, length(tpos) / nc)
        end
        _warn("Vector of tpos is not supported yet.")
        tpos = collect(tpos)[1]
    end

    @assert contours >= 0 "contours must be ≥ 0."
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(amethod, [:mean, :median], "amethod")
    _check_var(
        nmethod,
        [
            :zscore,
            :minmax,
            :log,
            :log10,
            :neglog,
            :neglog10,
            :neg,
            :pos,
            :perc,
            :gauss,
            :invroot,
            :n,
            :softmax,
            :sigmoid,
            :mad,
            :rank,
            :none,
        ],
        "nmethod",
    )

    # get channels and selected channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    if !isnothing(sch)
        if isa(sch, String)
            @assert length(intersect(ch, get_channel(obj, ch = sch))) == 1 "sch channel was not found in ch."
        else
            @assert length(intersect(ch, get_channel(obj, ch = sch))) == length(sch) "Some sch channels were not found in ch."
        end
    end
    @assert length(ch) >= 2 "plot_topo() requires ≥ 2 channels."
    chs = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
    !isnothing(sch) && (sch = _find_bylabel(locs, sch))

    # prepare data or time position
    if isnothing(data)
        @assert !isnothing(tpos) "Either tpos or data must be provided."
        @assert tpos >= obj.time_pts[1] "tpos must be ≥ $(obj.time_pts[1])"
        @assert tpos <= obj.time_pts[end] "tpos must be ≤ $(obj.time_pts[end])"
        tpos = vsearch(tpos, obj.time_pts)
        if nepochs(obj) == 1
            data = obj.data[ch, tpos, 1]
        else
            data = reshape(obj.data, size(obj.data, 1), size(obj.data, 2) * size(obj.data, 3))[ch, tpos]
        end
        title == "default" && (title = "$(obj.time_pts[tpos]) s")
    else
        !isnothing(tpos) && _info("If data is provided, tpos is ignored")
        if ndims(data) == 2
            data = amethod === :mean ? mean(data; dims = 2)[:] : median(data, dims = 2)[:]
        end
        @assert length(data) == length(ch) "Number of channels in data ($(length(data))) must equal the number of channels to plot ($(length(ch)))."
        title == "default" && (title = "")
    end

    cb_title == "default" && (cb_title = "[A.U.]")

    p = plot_topo(
        data;
        locs = locs,
        ch = collect(1:DataFrames.nrow(locs)),
        sch = sch,
        cb = cb,
        cb_title = cb_title,
        title = title,
        mono = mono,
        imethod = imethod,
        nmethod = nmethod,
        contours = contours,
        electrodes = electrodes,
        ps = ps,
        head = head,
        cart = cart,
        threshold = threshold,
        threshold_type = threshold_type,
        threshold_method = threshold_method,
    )

    return p

end
